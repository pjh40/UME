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
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>

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

std::array<double, 3> centroid_extent(
    std::span<int const> v, DS_Types::VEC3V_T const &zcoord) {
  std::array<double, 3> lo, hi;
  lo.fill(std::numeric_limits<double>::max());
  hi.fill(std::numeric_limits<double>::lowest());
  for (int const z : v)
    for (std::size_t d = 0; d < 3; ++d) {
      double const x = zcoord[static_cast<std::size_t>(z)][d];
      lo[d] = std::min(lo[d], x);
      hi[d] = std::max(hi[d], x);
    }
  std::array<double, 3> ext{0.0, 0.0, 0.0};
  if (!v.empty())
    for (std::size_t d = 0; d < 3; ++d)
      ext[d] = hi[d] - lo[d];
  return ext;
}

std::array<int, 3> dims_create(
    int nparts, std::array<double, 3> const &extent, int naxes) {
  if (nparts < 1)
    throw std::invalid_argument("dims_create: nparts must be positive");
  naxes = std::clamp(naxes, 1, 3);

  /* The axes we are allowed to divide, longest extent first. */
  std::array<int, 3> axis{0, 1, 2};
  std::ranges::stable_sort(axis, [&extent](int a, int b) {
    return extent[static_cast<std::size_t>(a)] >
        extent[static_cast<std::size_t>(b)];
  });

  std::array<int, 3> dims{1, 1, 1};
  /* Prime factors, largest first: a big factor placed early lands on the axis
     that can best absorb it, and the small ones then even the shape out. */
  std::vector<int> factors;
  for (int f = 2, n = nparts; n > 1; ++f) {
    if (f * f > n) {
      factors.push_back(n);
      break;
    }
    while (n % f == 0) {
      factors.push_back(f);
      n /= f;
    }
  }
  std::ranges::sort(factors, std::greater<int>{});

  for (int const f : factors) {
    int best = axis[0];
    for (int k = 1; k < naxes; ++k) {
      int const c = axis[static_cast<std::size_t>(k)];
      /* dims[c]/extent[c] < dims[best]/extent[best], cross-multiplied so that
         a zero extent divides nothing rather than dividing by zero.  The
         comparison is strict, so a tie keeps the longer axis. */
      if (dims[static_cast<std::size_t>(c)] *
              extent[static_cast<std::size_t>(best)] <
          dims[static_cast<std::size_t>(best)] *
              extent[static_cast<std::size_t>(c)])
        best = c;
    }
    dims[static_cast<std::size_t>(best)] *= f;
  }
  return dims;
}

namespace {

//! The unit eigenvector of the largest eigenvalue of a symmetric 3x3 matrix
/*! Cyclic Jacobi, which for a 3x3 is three rotations a sweep and converges in
    a handful of sweeps.  It is exact on a matrix that is already diagonal --
    the case an axis-aligned mesh presents -- so `rib` reproduces `rcb` there
    rather than wandering off by a rounding error. */
Vec3 dominant_eigenvector(std::array<std::array<double, 3>, 3> a) {
  std::array<Vec3, 3> vec{Vec3{{1.0, 0.0, 0.0}}, Vec3{{0.0, 1.0, 0.0}},
      Vec3{{0.0, 0.0, 1.0}}};

  for (int sweep = 0; sweep < 32; ++sweep) {
    double off = 0.0;
    for (std::size_t p = 0; p + 1 < 3; ++p)
      for (std::size_t q = p + 1; q < 3; ++q)
        off += std::abs(a[p][q]);
    double const scale =
        std::abs(a[0][0]) + std::abs(a[1][1]) + std::abs(a[2][2]);
    if (off <= std::numeric_limits<double>::epsilon() * scale)
      break;

    for (std::size_t p = 0; p + 1 < 3; ++p)
      for (std::size_t q = p + 1; q < 3; ++q) {
        if (a[p][q] == 0.0)
          continue;
        double const theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        /* The smaller of the two roots, which is the rotation that keeps the
           eigenvalues in the order they are already nearly in.  A `theta` so
           large that it overflows leaves `t` at zero, which is the right
           answer: the off-diagonal is negligible beside the diagonal. */
        double const t = std::copysign(1.0, theta) /
            (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        double const c = 1.0 / std::sqrt(t * t + 1.0);
        double const s = t * c;
        double const h = t * a[p][q];
        a[p][p] -= h;
        a[q][q] += h;
        a[p][q] = a[q][p] = 0.0;

        std::size_t const r = 3 - p - q; //!< the axis this rotation fixes
        double const arp = a[r][p], arq = a[r][q];
        a[r][p] = a[p][r] = c * arp - s * arq;
        a[r][q] = a[q][r] = s * arp + c * arq;

        for (std::size_t i = 0; i < 3; ++i) {
          double const vp = vec[p][i], vq = vec[q][i];
          vec[p][i] = c * vp - s * vq;
          vec[q][i] = s * vp + c * vq;
        }
      }
  }

  std::size_t best = 0;
  for (std::size_t i = 1; i < 3; ++i)
    if (a[i][i] > a[best][best])
      best = i;
  Vec3 dir = vec[best];
  normalize(dir);

  /* An eigenvector is only defined up to sign, and the sign the rotations
     happen to leave is an accident.  Pin it down so that a mesh and the same
     mesh with its zones renumbered give the same cuts, not mirrored ones. */
  std::size_t lead = 0;
  for (std::size_t i = 1; i < 3; ++i)
    if (std::abs(dir[i]) > std::abs(dir[lead]))
      lead = i;
  if (dir[lead] < 0.0)
    dir *= -1.0;
  return dir;
}

} // namespace

Vec3 principal_axis(std::span<int const> v, DS_Types::VEC3V_T const &zcoord,
    std::span<int const> wgt) {
  Vec3 const x_axis{{1.0, 0.0, 0.0}};
  if (v.empty())
    return x_axis;

  auto raw = [wgt](int z) {
    return wgt.empty() ? 1.0
                       : std::max(
                             static_cast<double>(wgt[static_cast<std::size_t>(
                                 z)]),
                             0.0);
  };
  double total = 0.0;
  for (int const z : v)
    total += raw(z);
  /* Weights that are all zero say nothing about which axis is longest, so fall
     back to counting zones rather than dividing by zero. */
  bool const unit = !(total > 0.0);
  auto weight = [&](int z) { return unit ? 1.0 : raw(z); };
  double const W = unit ? static_cast<double>(v.size()) : total;

  std::array<double, 3> mean{0.0, 0.0, 0.0};
  for (int const z : v) {
    double const w = weight(z);
    for (std::size_t d = 0; d < 3; ++d)
      mean[d] += w * zcoord[static_cast<std::size_t>(z)][d];
  }
  for (std::size_t d = 0; d < 3; ++d)
    mean[d] /= W;

  std::array<std::array<double, 3>, 3> cov{};
  for (int const z : v) {
    double const w = weight(z);
    std::array<double, 3> r;
    for (std::size_t d = 0; d < 3; ++d)
      r[d] = zcoord[static_cast<std::size_t>(z)][d] - mean[d];
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        cov[i][j] += w * r[i] * r[j];
  }
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j)
      cov[i][j] /= W;

  return dominant_eigenvector(cov);
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

namespace {

//! Prepare `z2r` for `assign()`, and return the zones a method may place
std::vector<int> begin_assign(
    SOA_Idx::Mesh const &mesh, Options const &opt, Zone_Map &z2r) {
  if (opt.nparts < 1)
    throw std::invalid_argument("Partitioner: nparts must be positive");
  z2r.assign(static_cast<std::size_t>(mesh.zones.size()), -1);
  return real_zones(mesh);
}

//! The coordinate axis the centroids of `v` are spread furthest along
int longest_axis(std::span<int const> v, DS_Types::VEC3V_T const &zcoord) {
  std::array<double, 3> const ext = centroid_extent(v, zcoord);
  return static_cast<int>(std::ranges::max_element(ext) - ext.begin());
}

//! The unit vector along coordinate axis `a`
Vec3 axis_vector(int a) {
  Vec3 e{{0.0, 0.0, 0.0}};
  e[static_cast<std::size_t>(a)] = 1.0;
  return e;
}

//! Cut the mesh into slabs normal to one coordinate axis
class Slab : public Partitioner {
public:
  /*! `axis` is 0, 1 or 2, or -1 to cut normal to the longest extent of the
      mesh -- which is the axis that gives the smallest cut area for a single
      slab stack. */
  Slab(int axis, char const *name, char const *describe)
      : axis_{axis}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    auto const &zcoord = zone_centroids(mesh);
    std::vector<int> const wgt = zone_weights(mesh, opt);

    int const axis = (axis_ >= 0) ? axis_ : longest_axis(v, zcoord);
    std::vector<double> key(static_cast<std::size_t>(mesh.zones.size()), 0.0);
    for (int const z : v)
      key[static_cast<std::size_t>(z)] = zcoord[static_cast<std::size_t>(z)][static_cast<std::size_t>(axis)];

    split_sorted(v, key, wgt, 0, opt.nparts, z2r);
  }

private:
  int axis_;
  char const *name_;
  char const *describe_;
};

//! One axis of a grid split, and how many ways it is divided
struct Level {
  int axis;
  int nparts;
};

//! Scratch shared by every node of one grid split
struct Grid_Scratch {
  //! The centroid coordinate of every real zone, one array per axis
  std::array<std::vector<double>, 3> key;
  //! Group numbers for the node being split, overwritten as the walk descends
  Zone_Map group;
};

//! Split `v` along each level in turn, into `first_part + [0, prod(levels))`
/*! The split at each level is by centroid *weight*, not by position, so the
    cuts of one level need not line up between the groups of the level above:
    a column with heavier zones in it gets shorter boxes.  That staircase is
    the point -- it is what keeps a grid method balanced on a mesh whose zone
    cost varies -- and it costs nothing, because the decomposition never
    assumes a part is a box. */
void split_grid(std::vector<int> const &v, std::span<int const> wgt,
    std::span<Level const> levels, int first_part, Grid_Scratch &scratch,
    Zone_Map &z2r) {
  if (v.empty())
    return;
  if (levels.empty()) {
    for (int const z : v)
      z2r[static_cast<std::size_t>(z)] = first_part;
    return;
  }

  Level const &level = levels.front();
  std::span<Level const> const rest = levels.subspan(1);
  int stride = 1;
  for (Level const &l : rest)
    stride *= l.nparts;

  std::span<double const> const key =
      scratch.key[static_cast<std::size_t>(level.axis)];
  if (stride == 1) {
    split_sorted(v, key, wgt, first_part, level.nparts, z2r);
    return;
  }

  split_sorted(v, key, wgt, 0, level.nparts, scratch.group);
  std::vector<std::vector<int>> bucket(static_cast<std::size_t>(level.nparts));
  for (int const z : v)
    bucket[static_cast<std::size_t>(scratch.group[static_cast<std::size_t>(z)])]
        .push_back(z);
  for (int i = 0; i < level.nparts; ++i)
    split_grid(bucket[static_cast<std::size_t>(i)], wgt, rest,
        first_part + i * stride, scratch, z2r);
}

//! Cut the mesh into a grid of boxes
class Block : public Partitioner {
public:
  /*! `naxes` is how many of the mesh's axes may be cut: 3 for a 3-D grid of
      boxes, 2 for columns spanning the shortest axis. */
  Block(int naxes, char const *name, char const *describe)
      : naxes_{naxes}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    auto const &zcoord = zone_centroids(mesh);
    std::vector<int> const wgt = zone_weights(mesh, opt);

    std::array<int, 3> const dims =
        dims_create(opt.nparts, centroid_extent(v, zcoord), naxes_);

    /* Undivided axes drop out, so a `pencil` of a flat mesh and a `slab` are
       the same one-level split rather than two levels one of which is a
       no-op. */
    std::vector<Level> levels;
    for (int d = 0; d < 3; ++d)
      if (dims[static_cast<std::size_t>(d)] > 1)
        levels.push_back(Level{d, dims[static_cast<std::size_t>(d)]});

    Grid_Scratch scratch;
    scratch.group.assign(static_cast<std::size_t>(mesh.zones.size()), -1);
    for (Level const &l : levels) {
      auto &key = scratch.key[static_cast<std::size_t>(l.axis)];
      key.assign(static_cast<std::size_t>(mesh.zones.size()), 0.0);
      for (int const z : v)
        key[static_cast<std::size_t>(z)] = zcoord[static_cast<std::size_t>(z)]
                                                [static_cast<std::size_t>(
                                                    l.axis)];
    }

    split_grid(v, wgt, levels, 0, scratch, z2r);
  }

private:
  int naxes_;
  char const *name_;
  char const *describe_;
};

//! Scratch shared by every node of one recursive bisection
struct Bisect_Scratch {
  //! The projection of each real zone's centroid onto the current direction
  std::vector<double> key;
  //! Group numbers for the node being split, overwritten as the walk descends
  Zone_Map group;
};

//! Cut `v` in two along `direction(v)`, and each half again, `k` parts in all
/*! The direction is chosen afresh for every piece -- that is the difference
    between this and a grid method, and what lets it follow a mesh that is not
    a box.

    An odd `k` has to be split unevenly, and the ratio is taken by splitting
    `k` ways and keeping the first `k1` groups rather than by aiming at a
    weight fraction directly.  That buys two things from `split_sorted`: the
    halves are in the ratio `k1 : k - k1` by weight, and each of the `k` groups
    is non-empty whenever there are at least `k` zones, so each half arrives at
    its own recursion with at least as many zones as it has parts to fill. */
void bisect(std::vector<int> const &v, std::span<int const> wgt,
    DS_Types::VEC3V_T const &zcoord,
    std::function<Vec3(std::span<int const>)> const &direction, int first_part,
    int k, Bisect_Scratch &scratch, Zone_Map &z2r) {
  if (v.empty())
    return;
  if (k == 1) {
    for (int const z : v)
      z2r[static_cast<std::size_t>(z)] = first_part;
    return;
  }

  Vec3 const dir = direction(v);
  for (int const z : v)
    scratch.key[static_cast<std::size_t>(z)] =
        dotprod(zcoord[static_cast<std::size_t>(z)], dir);

  int const k1 = (k + 1) / 2;
  split_sorted(v, scratch.key, wgt, 0, k, scratch.group);

  /* Both halves are read out before either recurses, because the recursion
     overwrites the scratch this just wrote. */
  std::vector<int> half[2];
  for (int const z : v)
    half[scratch.group[static_cast<std::size_t>(z)] < k1 ? 0 : 1].push_back(z);

  bisect(half[0], wgt, zcoord, direction, first_part, k1, scratch, z2r);
  bisect(half[1], wgt, zcoord, direction, first_part + k1, k - k1, scratch,
      z2r);
}

//! Cut the mesh in two, and each half again, choosing a direction each time
class Recursive : public Partitioner {
public:
  /*! `inertial` picks each cut normal to the principal axis of the piece being
      cut rather than to its longest *coordinate* axis, which is what lets it
      follow a mesh whose long direction is not an axis of the coordinate
      system. */
  Recursive(bool inertial, char const *name, char const *describe)
      : inertial_{inertial}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    auto const &zcoord = zone_centroids(mesh);
    std::vector<int> const wgt = zone_weights(mesh, opt);

    std::function<Vec3(std::span<int const>)> direction;
    if (inertial_)
      direction = [&zcoord, &wgt](std::span<int const> piece) {
        return principal_axis(piece, zcoord, wgt);
      };
    else
      direction = [&zcoord](std::span<int const> piece) {
        return axis_vector(longest_axis(piece, zcoord));
      };

    Bisect_Scratch scratch;
    scratch.key.assign(static_cast<std::size_t>(mesh.zones.size()), 0.0);
    scratch.group.assign(static_cast<std::size_t>(mesh.zones.size()), -1);
    bisect(v, wgt, zcoord, direction, 0, opt.nparts, scratch, z2r);
  }

private:
  bool inertial_;
  char const *name_;
  char const *describe_;
};

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

/*! The methods this source knows about.  Entries whose `requires_pkg` is set
    are compiled in or out with the package they need, but stay in the table
    either way so that an unavailable method reads differently from a
    misspelled one. */
constexpr Method TABLE[] = {
    {"block", "a grid of boxes, shaped to the extents of the mesh", nullptr,
        true},
    {"linear", "contiguous ranges of the input zone numbering", nullptr, true},
    {"pencil", "columns spanning the shortest extent of the mesh", nullptr,
        true},
    {"rcb", "recursive bisection normal to each piece's longest axis", nullptr,
        true},
    {"rib", "recursive bisection normal to each piece's principal axis",
        nullptr, true},
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
  if (name == "block")
    return std::make_unique<Block>(3, m->name, m->describe);
  if (name == "linear")
    return std::make_unique<Linear>();
  if (name == "pencil")
    return std::make_unique<Block>(2, m->name, m->describe);
  if (name == "rcb")
    return std::make_unique<Recursive>(false, m->name, m->describe);
  if (name == "rib")
    return std::make_unique<Recursive>(true, m->name, m->describe);
  if (name == "slab")
    return std::make_unique<Slab>(-1, m->name, m->describe);
  if (name == "slab-x")
    return std::make_unique<Slab>(0, m->name, m->describe);
  if (name == "slab-y")
    return std::make_unique<Slab>(1, m->name, m->describe);
  if (name == "slab-z")
    return std::make_unique<Slab>(2, m->name, m->describe);
  return nullptr;
}

} // namespace Part
} // namespace Ume
