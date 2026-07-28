/*
  Copyright (c) 2023, Triad National Security, LLC. All rights reserved.

  This is open source software; you can redistribute it and/or modify it under
  the terms of the BSD-3 License. If software is modified to produce derivative
  works, such modified software should be clearly marked, so as not to confuse
  it with the version available from LANL. Full text of the BSD-3 License can be
  found in the LICENSE.md file, and the full assertion of copyright in the
  NOTICE.md file.
*/

/*!
  \file Ume/Partition_BFS.cc

  `bfs`: recursive bisection of the dual graph, each cut grown outwards from a
  seed on the rim of the piece being split.

  Nothing here reads a coordinate.  That is the whole point of the method --
  the parts follow the mesh's connectivity rather than its shape -- and it is
  also what makes it the zero-dependency stand-in for the library methods,
  which minimize the same quantity by better means.
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include "Ume/Partition_Impl.hh"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace Ume {
namespace Part {

namespace {

//! The vertex `dist` says is furthest away, the lowest index breaking ties
int furthest(std::vector<int> const &dist) {
  return static_cast<int>(std::ranges::max_element(dist) - dist.begin());
}

//! Scratch shared by every node of one recursive graph bisection
struct Carve_Scratch {
  //! Which node of the bisection tree each vertex currently belongs to
  std::vector<int> tag;
  //! The edge weight at each vertex counting only edges inside its node
  std::vector<int> deg;
  //! The edge weight tying each vertex to the region being grown
  std::vector<int> pull;
  //! Whether the region being grown has taken each vertex
  std::vector<char> in;
  //! Distances for the seed search, kept here only to be allocated once
  std::vector<int> dist;
  //! Numbers the nodes, so that a node is one integer rather than a set
  int next_tag{0};
  //! The weight already placed in parts that are final, in part order
  std::int64_t taken{0};
};

//! A vertex out on the rim of the untaken part of node `node`
/*! One BFS sweep through the node from the first vertex left in it: the
    vertex it ends on is on the rim rather than in the middle.  Which matters,
    because a region grown from the middle leaves the rest of the node as a
    ring around itself, and a ring is both longer to bound and easily broken
    in two, where a region grown from the rim leaves it in one piece.

    A vertex the sweep never reaches stays at `FAR` and wins outright, so a
    node that falls into several pieces gets the region restarted in one of
    the stranded ones rather than leaving them to be swept up at the end. */
int rim_of(Zone_Graph const &g, int node, Carve_Scratch &s) {
  constexpr int FAR = std::numeric_limits<int>::max();
  auto open_ = [&s, node](int v) {
    auto const i = static_cast<std::size_t>(v);
    return s.tag[i] == node && !s.in[i];
  };

  //! Everything outside the node sits below every distance, so it cannot win
  std::ranges::fill(s.dist, -1);
  int src = -1;
  for (int v = 0; v < g.nverts; ++v)
    if (open_(v)) {
      s.dist[static_cast<std::size_t>(v)] = FAR;
      if (src < 0)
        src = v;
    }
  if (src < 0)
    return -1;

  s.dist[static_cast<std::size_t>(src)] = 0;
  std::vector<int> queue{src};
  for (std::size_t head = 0; head < queue.size(); ++head) {
    int const v = queue[head];
    int const nd = s.dist[static_cast<std::size_t>(v)] + 1;
    for (int const u : g.neighbors(v))
      if (open_(u) && nd < s.dist[static_cast<std::size_t>(u)]) {
        s.dist[static_cast<std::size_t>(u)] = nd;
        queue.push_back(u);
      }
  }
  return furthest(s.dist);
}

//! Grow a region inside node `node` until it carries `target` weight
/*! The region starts on the rim and repeatedly takes whichever vertex would
    shorten its own boundary the most.  `least` and `most` bound how many
    vertices it may end up with, which is how each side of the split is
    guaranteed a vertex per part it has to fill even when a few very heavy
    zones would otherwise let one side reach its target having eaten the node.

    The result is written into `s.in`, which the caller must have cleared over
    the node. */
void grow_region(Zone_Graph const &g, std::span<int const> v, int node,
    std::int64_t target, int least, int most, Carve_Scratch &s) {
  for (int const u : v) {
    auto const i = static_cast<std::size_t>(u);
    s.pull[i] = 0;
    /* The gain is measured against the boundary *of this node*: an edge
       leaving the node was cut at a level above and cannot be cut again, so
       counting it here would make a vertex on the node's own rim look
       expensive when it is the cheapest thing available. */
    s.deg[i] = 0;
    std::span<int const> const nbr = g.neighbors(u);
    std::span<int const> const w = g.weights(u);
    for (std::size_t j = 0; j < nbr.size(); ++j)
      if (s.tag[static_cast<std::size_t>(nbr[j])] == node)
        s.deg[i] += w[j];
  }

  /*! What taking `u` into the region would do to its boundary: the edge
      weight it would bury, less the edge weight it would expose.  Ranking by
      the pull alone would take a well-connected vertex however much new
      boundary came with it, and that is the difference between a cut that
      flattens itself out and one that wanders. */
  auto gain_of = [&s](int u) {
    return 2 * s.pull[static_cast<std::size_t>(u)] -
        s.deg[static_cast<std::size_t>(u)];
  };

  //! What the region might take next, best first, stale entries and all
  std::priority_queue<std::pair<int, int>> frontier;
  std::int64_t weight = 0;
  int count = 0;

  auto take = [&](int u) {
    s.in[static_cast<std::size_t>(u)] = 1;
    weight += g.vwgt[static_cast<std::size_t>(u)];
    ++count;
    std::span<int const> const nbr = g.neighbors(u);
    std::span<int const> const w = g.weights(u);
    for (std::size_t j = 0; j < nbr.size(); ++j) {
      auto const n = static_cast<std::size_t>(nbr[j]);
      if (s.tag[n] != node || s.in[n])
        continue;
      s.pull[n] += w[j];
      /* A vertex's gain towards the region only ever rises, so the entry
         pushed here is the highest for that vertex and comes out first; the
         ones underneath are recognised by their stale gain and dropped.
         Rewriting the queue instead would cost a search per neighbor per
         step. */
      frontier.push(std::pair{gain_of(nbr[j]), -nbr[j]});
    }
  };

  //! The vertex the region would take next, or -1 if its frontier is spent
  auto best = [&] {
    while (!frontier.empty()) {
      auto const [gain, negu] = frontier.top();
      int const u = -negu;
      /* Largest gain first and, at equal gain, the lowest vertex index, which
         is what the negated index in the key buys. */
      if (!s.in[static_cast<std::size_t>(u)] && gain == gain_of(u))
        return u;
      frontier.pop();
    }
    return -1;
  };

  take(rim_of(g, node, s));
  while (count < most && (weight < target || count < least)) {
    int u = best();
    if (u < 0)
      /* The node is in more than one piece and the region has filled the one
         it started in.  Carry on in another rather than stop short: a region
         in two pieces is a worse answer than one, but a part left light is a
         worse answer than either. */
      u = rim_of(g, node, s);
    take(u);
  }
}

//! Split `v` in two, and each half again, until every part has its zones
/*! The cut at each level is grown rather than sorted, which is the whole
    point of the method: nothing here reads a coordinate, so the parts follow
    the mesh's connectivity and not its shape.  Whether that is worth having
    is a question about the mesh -- on a brick the geometric methods know the
    answer already, and on a mesh that doubles back on itself they do not.

    Each split aims at a *cumulative* target -- the weight that should have
    been placed once every part below `a + k1` is full -- rather than at a
    share of the node in front of it.  Because a region always stops on the
    first vertex that reaches its target, every one of those cumulative
    boundaries lands within one zone's weight above its ideal, and a part's
    own weight is the difference of two of them: within one zone's weight of
    the mean, which is the bound `split_sorted` gives the methods that sort.
    Aiming at a share of the node would let each level's rounding accumulate
    into the next instead. */
void carve(Zone_Graph const &g, std::vector<int> const &v, int a, int k,
    int node, int np, std::int64_t total, Carve_Scratch &s,
    std::vector<int> &part) {
  if (v.empty())
    return;
  if (k == 1) {
    for (int const u : v) {
      part[static_cast<std::size_t>(u)] = a;
      s.taken += g.vwgt[static_cast<std::size_t>(u)];
    }
    return;
  }

  int const k1 = (k + 1) / 2;
  std::int64_t const target =
      total * static_cast<std::int64_t>(a + k1) / np - s.taken;

  for (int const u : v)
    s.in[static_cast<std::size_t>(u)] = 0;
  grow_region(g, v, node, target, k1, static_cast<int>(v.size()) - (k - k1), s);

  int const tag_lo = ++s.next_tag, tag_hi = ++s.next_tag;
  std::vector<int> half[2];
  for (int const u : v) {
    bool const lo = s.in[static_cast<std::size_t>(u)] != 0;
    s.tag[static_cast<std::size_t>(u)] = lo ? tag_lo : tag_hi;
    half[lo ? 0 : 1].push_back(u);
  }

  /* Depth first and low half first, so that `s.taken` is the weight of every
     part below this node by the time the next one is split. */
  carve(g, half[0], a, k1, tag_lo, np, total, s, part);
  carve(g, half[1], a + k1, k - k1, tag_hi, np, total, s, part);
}

//! Grow every part outwards along the mesh's own connectivity
class Bfs : public Partitioner {
public:
  char const *name() const override { return "bfs"; }
  char const *describe() const override {
    return "parts grown outwards along the dual graph";
  }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    begin_assign(mesh, opt, z2r);
    Zone_Graph const g = build_zone_graph(mesh, opt);
    if (g.nverts == 0)
      return;
    /* More parts than there are zones to fill them: the low parts get the
       zones and the rest stay empty, which is what `split_sorted` does with
       the same shortage. */
    int const np = std::min(opt.nparts, g.nverts);

    auto const nvs = static_cast<std::size_t>(g.nverts);
    Carve_Scratch s;
    s.tag.assign(nvs, 0);
    s.deg.assign(nvs, 0);
    s.pull.assign(nvs, 0);
    s.in.assign(nvs, 0);
    s.dist.assign(nvs, -1);

    std::int64_t total = 0;
    for (int const w : g.vwgt)
      total += w;

    std::vector<int> all(nvs);
    std::iota(all.begin(), all.end(), 0);
    std::vector<int> part(nvs, -1);
    carve(g, all, 0, np, 0, np, total, s, part);

    for (int v = 0; v < g.nverts; ++v)
      z2r[static_cast<std::size_t>(g.v2z[static_cast<std::size_t>(v)])] =
          part[static_cast<std::size_t>(v)];
  }
};

} // namespace

std::unique_ptr<Partitioner> make_bfs(Method const &m) {
  if (std::string_view{m.name} == "bfs")
    return std::make_unique<Bfs>();
  return nullptr;
}

} // namespace Part
} // namespace Ume
