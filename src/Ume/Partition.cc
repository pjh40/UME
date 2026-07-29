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
  \file Ume/Partition.cc

  The method table, `make()`, the pieces every method shares, and the three
  methods that need nothing from the mesh but the list of its real zones.

  The methods that read a coordinate are in `Partition_Geometric.cc` and the
  one that reads the dual graph is in `Partition_BFS.cc`; `Partition_Impl.hh`
  says why they are apart and how to add another file of them.
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Impl.hh"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace Ume {
namespace Part {

std::vector<int> real_zones(SOA_Idx::Mesh const &mesh) {
  auto const &zmask = mesh.zones.mask;
  std::vector<int> v;
  v.reserve(zmask.size());
  for (int z = 0; z < mesh.zones.size(); ++z)
    if (zmask[z] >= 1)
      v.push_back(z);
  return v;
}

std::vector<int> zone_weights(SOA_Idx::Mesh const &mesh, Options const &opt) {
  int const zll = mesh.zones.size();
  auto const &zmask = mesh.zones.mask;
  std::vector<int> w(static_cast<std::size_t>(zll), 0);

  switch (opt.vwgt) {
  case Vertex_Weight::UNIT:
    for (int z = 0; z < zll; ++z)
      w[static_cast<std::size_t>(z)] = (zmask[z] >= 1) ? 1 : 0;
    break;
  case Vertex_Weight::CORNERS: {
    /* The corner count is the trip count of every zone-gathering kernel, so
       balancing zones rather than corners is misleading on a hex/tet mix. */
    auto const &z2c = mesh.ds->caccess_intrr("m:z>c");
    for (int z = 0; z < zll; ++z)
      if (zmask[z] >= 1)
        w[static_cast<std::size_t>(z)] = static_cast<int>(z2c[z].size());
    break;
  }
  case Vertex_Weight::SIDES: {
    /* There is no `m:z>s`, so count the real sides the long way round. */
    auto const &s2z = mesh.ds->caccess_intv("m:s>z");
    auto const &smask = mesh.sides.mask;
    for (int s = 0; s < mesh.sides.size(); ++s) {
      if (smask[s] < 1)
        continue;
      int const z = s2z[s];
      if (z >= 0 && z < zll && zmask[z] >= 1)
        w[static_cast<std::size_t>(z)] += 1;
    }
    break;
  }
  }
  return w;
}

void split_sorted(std::span<int const> v, std::span<double const> key,
    std::span<int const> wgt, int first_part, int k, Zone_Map &z2r) {
  if (k < 1)
    throw std::invalid_argument("split_sorted: k must be positive");
  if (v.empty())
    return;

  auto weight = [wgt](int z) -> std::int64_t {
    return wgt.empty() ? 1 : wgt[static_cast<std::size_t>(z)];
  };

  std::vector<int> order(v.begin(), v.end());
  std::ranges::sort(order, [key](int a, int b) {
    double const ka = key[static_cast<std::size_t>(a)];
    double const kb = key[static_cast<std::size_t>(b)];
    return ka != kb ? ka < kb : a < b;
  });

  std::int64_t total = 0;
  for (int const z : order)
    total += weight(z);

  int const n = static_cast<int>(order.size());
  std::int64_t acc = 0;
  for (int i = 0; i < n; ++i) {
    int const z = order[static_cast<std::size_t>(i)];
    std::int64_t const w = weight(z);
    int p;
    if (n < k) {
      /* Not enough zones to go round.  Fill from the first part and leave the
         tail empty, rather than letting the weight rule scatter a handful of
         zones across the middle of the range and strand part 0. */
      p = i;
    } else {
      /* Charge the zone to the part its weight *midpoint* falls in, which
         splits a heavy zone to whichever side it mostly belongs on rather
         than always to the left. */
      p = (total > 0) ? static_cast<int>(((2 * acc + w) * k) / (2 * total))
                      : (i * k) / n;
      /* Nothing in that stops a run of heavy zones from stepping over a part
         entirely.  In sorted order part `p` can only be reached once `p`
         zones have been placed, and must be reached with `k - p` zones left,
         so those two bounds make the assignment onto `[0, k)`.  Both are
         non-decreasing in `i`, so clamping to them keeps each part
         contiguous. */
      p = std::clamp(std::clamp(p, 0, k - 1), k - (n - i), i);
    }
    z2r[static_cast<std::size_t>(z)] = first_part + p;
    acc += w;
  }
}

std::vector<int> begin_assign(
    SOA_Idx::Mesh const &mesh, Options const &opt, Zone_Map &z2r) {
  if (opt.nparts < 1)
    throw std::invalid_argument("Partitioner: nparts must be positive");
  z2r.assign(static_cast<std::size_t>(mesh.zones.size()), -1);
  return real_zones(mesh);
}

namespace {

//! Cut the mesh into contiguous runs of the original zone numbering
class Linear : public Partitioner {
public:
  char const *name() const override { return "linear"; }
  char const *describe() const override {
    return "contiguous ranges of the input zone numbering";
  }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    std::vector<int> const wgt = zone_weights(mesh, opt);
    std::vector<double> key(static_cast<std::size_t>(mesh.zones.size()), 0.0);
    for (int const z : v)
      key[static_cast<std::size_t>(z)] = static_cast<double>(z);
    split_sorted(v, key, wgt, 0, opt.nparts, z2r);
  }
};

//! Deal the real zones out to the parts one at a time
class Round_Robin : public Partitioner {
public:
  char const *name() const override { return "roundrobin"; }
  char const *describe() const override {
    return "the i'th real zone to part i % nparts";
  }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    /* Dealt over the *real* zones rather than over the raw zone index, so that
       a mesh whose boundary zones are interleaved with its real ones -- or, as
       `gen_brick` does it, appended after them -- deals the same way.  A raw
       `z % nparts` would leave the residue classes of the skipped zones short.

       There is no weighting here and there is deliberately none: this method
       exists to be the worst case, and dealing by count is the thing that
       makes almost every face a cut face and the communication graph
       complete.  It is the method the layout machinery is hardest on. */
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    int p = 0;
    for (int const z : v) {
      z2r[static_cast<std::size_t>(z)] = p;
      if (++p == opt.nparts)
        p = 0;
    }
  }
};

//! One step of SplitMix64, advancing `state`
/*! The built-in methods may not reach for `std::shuffle` or
    `std::uniform_int_distribution`: neither has a specified algorithm, so the
    same seed gives different partitions under libstdc++ and libc++, and
    `--seed` would only reproduce a run on the machine that made it. */
std::uint64_t next_random(std::uint64_t &state) {
  std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

//! A uniform draw from `[0, n)`, with `n > 0`
std::uint64_t bounded_random(std::uint64_t &state, std::uint64_t n) {
  /* 2^64 % n draws have to be thrown away for the rest to be uniform: taking
     the modulus of every draw would favour the first 2^64 % n values. */
  std::uint64_t const reject = (~n + 1u) % n;
  std::uint64_t x;
  do {
    x = next_random(state);
  } while (x < reject);
  return x % n;
}

//! Scatter the zones at random, then cut the shuffle into equal-weight runs
class Random : public Partitioner {
public:
  char const *name() const override { return "random"; }
  char const *describe() const override {
    return "a shuffle of the zones determined by the seed";
  }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> v = begin_assign(mesh, opt, z2r);
    std::vector<int> const wgt = zone_weights(mesh, opt);

    /* Fisher-Yates, which is a bijection onto the permutations however the
       draws come out, so no zone can be dropped or placed twice. */
    std::uint64_t state = opt.seed;
    for (std::size_t i = v.size(); i > 1; --i)
      std::swap(v[i - 1], v[bounded_random(state, i)]);

    /* The shuffle only fixes the *order*; the weight balance still comes from
       `split_sorted`, so `random` is a partition that is as badly shaped as
       round-robin and as well balanced as the geometric methods.  Keeping the
       two effects apart is what makes it a useful control. */
    std::vector<double> key(static_cast<std::size_t>(mesh.zones.size()), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i)
      key[static_cast<std::size_t>(v[i])] = static_cast<double>(i);

    split_sorted(v, key, wgt, 0, opt.nparts, z2r);
  }
};

//! `linear`, `roundrobin`, `random`: the methods that read no geometry
std::unique_ptr<Partitioner> make_ordered(Method const &m) {
  std::string_view const name{m.name};
  if (name == "linear")
    return std::make_unique<Linear>();
  if (name == "random")
    return std::make_unique<Random>();
  if (name == "roundrobin")
    return std::make_unique<Round_Robin>();
  return nullptr;
}

//! Was this build configured with METIS?
constexpr bool HAVE_METIS =
#ifdef NO_METIS
    false;
#else
    true;
#endif

//! Was this build configured with Mt-KaHyPar?
constexpr bool HAVE_KAHYPAR =
#ifdef NO_KAHYPAR
    false;
#else
    true;
#endif

/*! The methods this source knows about.  Entries whose `requires_pkg` is set
    are compiled in or out with the package they need, but stay in the table
    either way so that an unavailable method reads differently from a
    misspelled one. */
constexpr Method TABLE[] = {
    {"bfs", "parts grown outwards along the dual graph", nullptr, true},
    {"bfs-kl", "grown parts, with each cut refined by Kernighan-Lin swaps",
        nullptr, true},
    {"block", "a grid of boxes, shaped to the extents of the mesh", nullptr,
        true},
    {"hilbert", "equal runs of a Hilbert curve through the zones", nullptr,
        true},
    {"kahypar", "Mt-KaHyPar, minimizing the points on the wire", "Mt-KaHyPar",
        HAVE_KAHYPAR},
    {"kahypar-det", "Mt-KaHyPar, reproducibly, at some cost in quality",
        "Mt-KaHyPar", HAVE_KAHYPAR},
    {"kahypar-quality", "Mt-KaHyPar with flow-based refinement; slower",
        "Mt-KaHyPar", HAVE_KAHYPAR},
    {"linear", "contiguous ranges of the input zone numbering", nullptr, true},
    {"metis-kway-cut", "METIS k-way, minimizing the weight of the cut faces",
        "METIS", HAVE_METIS},
    {"metis-kway-vol", "METIS k-way, minimizing its communication-volume proxy",
        "METIS", HAVE_METIS},
    {"metis-rb", "METIS recursive bisection of the dual graph", "METIS",
        HAVE_METIS},
    {"morton", "equal runs of a Morton (Z-order) curve through the zones",
        nullptr, true},
    {"pencil", "columns spanning the shortest extent of the mesh", nullptr,
        true},
    {"random", "a shuffle of the zones determined by the seed", nullptr, true},
    {"rcb", "recursive bisection normal to each piece's longest axis", nullptr,
        true},
    {"rib", "recursive bisection normal to each piece's principal axis",
        nullptr, true},
    {"roundrobin", "the i'th real zone to part i % nparts", nullptr, true},
    {"slab", "slabs normal to the longest extent of the mesh", nullptr, true},
    {"slab-x", "slabs normal to x", nullptr, true},
    {"slab-y", "slabs normal to y", nullptr, true},
    {"slab-z", "slabs normal to z", nullptr, true},
};

} // namespace

std::span<Method const> all_methods() { return std::span<Method const>{TABLE}; }

std::vector<std::string> methods() {
  std::vector<std::string> names;
  for (Method const &m : TABLE)
    if (m.available)
      names.emplace_back(m.name);
  return names;
}

Method const *find_method(std::string_view name) {
  for (Method const &m : TABLE)
    if (name == m.name)
      return &m;
  return nullptr;
}

std::unique_ptr<Partitioner> make(std::string_view name) {
  Method const *const m = find_method(name);
  if (!m || !m->available)
    return nullptr;
  /* Each factory recognizes its own names and passes on the rest, so a new
     family of methods is a new file plus one line here.  A name that is in
     the table and in none of the factories is a table entry somebody forgot
     to implement, and comes back as nullptr -- which `ume_partition` reports
     the same way as an unknown method rather than crashing on it. */
  if (auto p = make_ordered(*m))
    return p;
  if (auto p = make_geometric(*m))
    return p;
  if (auto p = make_bfs(*m))
    return p;
  if (auto p = make_metis(*m))
    return p;
  if (auto p = make_kahypar(*m))
    return p;
  return nullptr;
}

} // namespace Part
} // namespace Ume
