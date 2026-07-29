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
  seed on the rim of the piece being split.  `bfs-kl` is the same thing with
  every cut then walked over by Kernighan--Lin.

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
  //! What moving each vertex to the other half of its cut would buy
  std::vector<int> d;
  //! Whether a refinement sweep has already moved each vertex
  std::vector<char> lock;
  //! Numbers the nodes, so that a node is one integer rather than a set
  int next_tag{0};
  //! The weight already placed in parts that are final, in part order
  std::int64_t taken{0};
  //! Whether each grown cut is handed to Kernighan--Lin afterwards
  bool refine{false};
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

/*! @name Kernighan--Lin refinement of one grown cut

    A grown region stops the moment it is heavy enough, so its boundary is
    wherever the frontier happened to be at that step: the greedy grows a shell
    before it grows a slab, and a cut can settle on the wrong one of a piece's
    two natural planes.  What it cannot do is go back, and going back is all
    this is.

    **Why swaps rather than single moves.**  The balance this method promises
    -- every part within one zone's weight of the mean -- comes from each grown
    region stopping on the *first* vertex to reach its cumulative target, which
    leaves the weight sitting right against the bottom of the window it is
    allowed.  A Fiduccia--Mattheyses pass moves one vertex at a time and would
    step straight out of that window on its first move in one of the two
    directions, so it could only ever push the cut one way.  A swap is weight
    neutral whenever the two zones weigh the same and nearly so when they do
    not, which is why the algorithm this wants is the older one. */
///@{

//! How many of the best candidates on each side a swap is chosen from
/*! Kernighan--Lin as published examines every pair, which is quadratic in the
    boundary per swap and cubic per sweep.  Ranking each side by its own gain
    and pairing only the top few costs a constant instead, and what makes that
    a fair ranking rather than an arbitrary one is that the pair term
    `2 * c(a, b)` is a correction to `d[a] + d[b]` rather than the leading
    behaviour: it is zero for the great majority of pairs, because two vertices
    on opposite sides of a cut usually do not share a face.

    Eight is where the quality stops being worth the arithmetic rather than
    where it stops improving.  On a 40x30x24 half-tet brick at sixteen parts:
    unrefined `bfs` cuts 14583 faces, a width of 1 -- pair the best of each
    side and accept it -- cuts 14550 and is thus worth nothing at all, 8 cuts
    12239, and 32 cuts 11485.  But the pair loop is quadratic in this number
    and runs once per swap, so 32 adds ten seconds to a run in which 8 adds
    none that can be measured against the cost of reading the mesh. */
inline constexpr int KL_WIDTH = 8;

//! How many sweeps over one cut to make before giving up on it
/*! A sweep that improves a cut leaves a different cut to sweep again, and the
    next one finds more.  Measured on a 12x10x8 half-tet brick summed over 2,
    3, 4, 8 and 16 parts: one sweep per cut leaves 2147 faces cut, two 1905,
    three 1784, and from there on it sits at 1786 and stops moving.  A sweep
    that gains nothing stops the loop, so this is a ceiling rather than a
    count. */
inline constexpr int KL_SWEEPS = 8;

//! The weight of the edge between `a` and `b`, or 0 if they do not share one
int edge_weight(Zone_Graph const &g, int a, int b) {
  std::span<int const> const nbr = g.neighbors(a);
  auto const it = std::ranges::lower_bound(nbr, b);
  if (it == nbr.end() || *it != b)
    return 0;
  return g.weights(a)[static_cast<std::size_t>(it - nbr.begin())];
}

//! One Kernighan--Lin sweep over the cut through node `node`; the gain it made
/*! Every vertex may be moved once.  The sweep swaps a vertex of each half at a
    time, best first, locking both and letting the running total go negative;
    what it keeps at the end is the prefix of that sequence with the highest
    cumulative gain, which is the one thing that makes this more than a greedy.

    It therefore runs on to the end of the boundary rather than stopping after
    a run of losing swaps.  Cutting it off there was tried and gives up most of
    the method: on the 12x10x8 brick a cap of 64 consecutive losing swaps left
    1958 faces cut against the full sweep's 1786.  A constant cap is the wrong
    shape whatever its value, because how far past the first hill the winning
    position lies grows with the boundary, and the boundary grows with the
    mesh.

    `wlow` is the weight of the low half, updated in place, and is held inside
    `[lo, hi]` throughout -- so a rolled-back sweep leaves both the cut and the
    balance exactly as it found them.  The return value is the edge weight no
    longer cut, and zero means there is nothing more here to find. */
std::int64_t kl_sweep(Zone_Graph const &g, std::span<int const> v, int node,
    std::int64_t lo, std::int64_t hi, std::int64_t &wlow, Carve_Scratch &s) {
  auto const low = [&s](int u) {
    return s.in[static_cast<std::size_t>(u)] != 0;
  };
  auto const dof = [&s](int u) -> int & {
    return s.d[static_cast<std::size_t>(u)];
  };
  auto const wgt = [&g](int u) {
    return static_cast<std::int64_t>(g.vwgt[static_cast<std::size_t>(u)]);
  };

  /* Indexed by side, 0 for the low half: `q[0]` is what the low half would
     give up.  Stale entries are left in and recognised on the way out by their
     key, exactly as `grow_region` does -- but here a gain can fall as well as
     rise, so an entry is pushed on *every* change rather than only on the way
     up, and the one carrying the current value is always the survivor. */
  std::priority_queue<std::pair<int, int>> q[2];
  for (int const u : v) {
    auto const i = static_cast<std::size_t>(u);
    s.lock[i] = 0;
    int ext = 0, own = 0;
    std::span<int const> const nbr = g.neighbors(u);
    std::span<int const> const w = g.weights(u);
    for (std::size_t j = 0; j < nbr.size(); ++j) {
      if (s.tag[static_cast<std::size_t>(nbr[j])] != node)
        continue;
      (low(nbr[j]) == low(u) ? own : ext) += w[j];
    }
    s.d[i] = ext - own;
    /* A vertex with nothing across the cut cannot be worth moving yet.  One
       becomes worth considering when a neighbour crosses, and is pushed then
       along with everything else whose gain changed. */
    if (ext > 0)
      q[low(u) ? 0 : 1].push(std::pair{s.d[i], -u});
  }

  //! The next `KL_WIDTH` vertices side `k` would most like to give up
  /*! A vertex can hold several entries carrying its current gain, because
      `give_back` pushes one for every candidate it returns and the entry it
      was drawn on may not have been the stale one.  They all pass the key
      test, so the duplicate has to be caught here or it would silently narrow
      the choice this function exists to widen. */
  auto const take = [&](int k, std::vector<int> &out) {
    out.clear();
    while (!q[k].empty() && static_cast<int>(out.size()) < KL_WIDTH) {
      auto const [key, negu] = q[k].top();
      int const u = -negu;
      q[k].pop();
      if (!s.lock[static_cast<std::size_t>(u)] && key == dof(u) &&
          std::ranges::find(out, u) == out.end())
        out.push_back(u);
    }
  };

  auto const give_back = [&](int k, std::vector<int> const &c) {
    for (int const u : c)
      q[k].push(std::pair{dof(u), -u});
  };

  //! Move `u` to the other half and tell its unlocked neighbours
  auto const cross = [&](int u) {
    auto const i = static_cast<std::size_t>(u);
    bool const was = low(u);
    s.in[i] = was ? 0 : 1;
    std::span<int const> const nbr = g.neighbors(u);
    std::span<int const> const w = g.weights(u);
    for (std::size_t j = 0; j < nbr.size(); ++j) {
      int const x = nbr[j];
      auto const n = static_cast<std::size_t>(x);
      if (s.tag[n] != node || s.lock[n])
        continue;
      /* `u` has left `was`.  A neighbour still on that side has traded an
         edge inside its half for one across the cut, and gains twice the
         edge's weight for crossing itself; one on the far side has traded the
         other way round. */
      s.d[n] += (low(x) == was ? 2 : -2) * w[j];
      q[low(x) ? 0 : 1].push(std::pair{s.d[n], -x});
    }
  };

  std::vector<std::pair<int, int>> made;
  std::vector<int> cand[2];
  std::int64_t running = 0, best = 0;
  std::size_t keep = 0;

  for (;;) {
    take(0, cand[0]);
    take(1, cand[1]);

    int pick_a = -1, pick_b = -1;
    std::int64_t pick = 0;
    for (int const a : cand[0])
      for (int const b : cand[1]) {
        /* The window is the whole balance argument, so it is checked before
           the gain rather than repaired afterwards. */
        std::int64_t const w = wlow - wgt(a) + wgt(b);
        if (w < lo || w > hi)
          continue;
        std::int64_t const gain = dof(a) + dof(b) - 2 * edge_weight(g, a, b);
        if (pick_a < 0 || gain > pick) {
          pick = gain;
          pick_a = a;
          pick_b = b;
        }
      }
    give_back(0, cand[0]);
    give_back(1, cand[1]);
    if (pick_a < 0)
      break;

    s.lock[static_cast<std::size_t>(pick_a)] = 1;
    s.lock[static_cast<std::size_t>(pick_b)] = 1;
    wlow += wgt(pick_b) - wgt(pick_a);
    cross(pick_a);
    cross(pick_b);

    made.emplace_back(pick_a, pick_b);
    running += pick;
    if (running > best) {
      best = running;
      keep = made.size();
    }
  }

  //! Everything after the high-water mark was speculative; undo it
  for (std::size_t i = made.size(); i > keep; --i) {
    auto const [a, b] = made[i - 1];
    s.in[static_cast<std::size_t>(a)] = 1;
    s.in[static_cast<std::size_t>(b)] = 0;
    wlow += wgt(a) - wgt(b);
  }
  return best;
}

//! Walk Kernighan--Lin over the cut `grow_region` just left in `s.in`
/*! `target` is the weight the low half was aimed at.  The window the sweeps
    have to stay inside is that aim, plus room for one zone of the node --
    which is precisely the room the growth itself was allowed, since a region
    stops on the first vertex to reach its target.  Both ends are widened to
    admit the weight the growth actually produced, so that a cut the count
    bounds in `grow_region` forced off its aim is refined rather than frozen. */
void refine_cut(Zone_Graph const &g, std::span<int const> v, int node,
    std::int64_t target, Carve_Scratch &s) {
  std::int64_t wlow = 0, heaviest = 0;
  for (int const u : v) {
    auto const i = static_cast<std::size_t>(u);
    auto const w = static_cast<std::int64_t>(g.vwgt[i]);
    heaviest = std::max(heaviest, w);
    if (s.in[i])
      wlow += w;
  }
  std::int64_t const lo = std::min(target, wlow);
  std::int64_t const hi = std::max(target + heaviest, wlow);
  for (int i = 0; i < KL_SWEEPS; ++i)
    if (kl_sweep(g, v, node, lo, hi, wlow, s) == 0)
      break;
}

///@}

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
  /* Before the halves are tagged apart, because the sweep needs both sides of
     the cut in one node -- and afterwards each half is a node of its own. */
  if (s.refine)
    refine_cut(g, v, node, target, s);

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
  //! `refine` hands each grown cut to Kernighan--Lin before recursing
  Bfs(bool refine, Method const &m) : refine_{refine}, m_{&m} {}

  char const *name() const override { return m_->name; }
  char const *describe() const override { return m_->describe; }

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
    s.refine = refine_;
    if (refine_) {
      s.d.assign(nvs, 0);
      s.lock.assign(nvs, 0);
    }

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

private:
  bool refine_;
  Method const *m_;
};

} // namespace

std::unique_ptr<Partitioner> make_bfs(Method const &m) {
  std::string_view const name{m.name};
  if (name == "bfs")
    return std::make_unique<Bfs>(false, m);
  if (name == "bfs-kl")
    return std::make_unique<Bfs>(true, m);
  return nullptr;
}

} // namespace Part
} // namespace Ume
