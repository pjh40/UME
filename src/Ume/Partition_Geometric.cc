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
  \file Ume/Partition_Geometric.cc

  The methods that partition by where the zones are: slab stacks, grids of
  boxes, recursive bisection, and space-filling curves.  Every one of them
  reduces to sorting the zone centroids by some key and handing the order to
  `split_sorted`, and they differ only in the key.

  The building blocks the keys are built from -- `centroid_extent`,
  `dims_create`, `principal_axis`, `morton_index` and `hilbert_index` -- are
  declared in `Partition.hh` rather than kept file-local, because each is a
  self-contained numerical routine that is worth testing on its own terms
  instead of only through the partitions it produces.
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include "Ume/Partition_Impl.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace Ume {
namespace Part {

namespace {

//! The corners of the smallest box holding the centroids of the zones `v`
/*! Both corners are the origin when `v` is empty. */
std::pair<std::array<double, 3>, std::array<double, 3>> centroid_bounds(
    std::span<int const> v, DS_Types::VEC3V_T const &zcoord) {
  if (v.empty())
    return {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  std::array<double, 3> lo, hi;
  lo.fill(std::numeric_limits<double>::max());
  hi.fill(std::numeric_limits<double>::lowest());
  for (int const z : v)
    for (std::size_t d = 0; d < 3; ++d) {
      double const x = zcoord[static_cast<std::size_t>(z)][d];
      lo[d] = std::min(lo[d], x);
      hi[d] = std::max(hi[d], x);
    }
  return {lo, hi};
}

} // namespace

std::array<double, 3> centroid_extent(
    std::span<int const> v, DS_Types::VEC3V_T const &zcoord) {
  auto const [lo, hi] = centroid_bounds(v, zcoord);
  std::array<double, 3> ext;
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

std::uint64_t morton_index(std::array<std::uint32_t, 3> const &q, int bits) {
  /* One bit of each axis per level, most significant level first.  The
     bit-twiddled version of this is a well known page of magic constants; at
     one call per zone of a one-off setup pass it buys nothing worth the
     reader's time. */
  std::uint64_t m = 0;
  for (int j = bits - 1; j >= 0; --j)
    for (std::size_t i = 0; i < 3; ++i)
      m = (m << 1) | ((q[i] >> j) & 1u);
  return m;
}

std::uint64_t hilbert_index(std::array<std::uint32_t, 3> const &q, int bits) {
  std::uint32_t const mask = (1u << bits) - 1u;
  std::array<std::uint32_t, 3> x{q[0] & mask, q[1] & mask, q[2] & mask};
  std::uint32_t const M = 1u << (bits - 1);

  /* Skilling's AxestoTranspose, "Programming the Hilbert curve" (2004).  The
     first loop reflects and rotates each level's octant into the orientation
     the level above left it in -- the undo of the recursive construction --
     and the second Gray-codes the result, after which the coordinates hold
     the index in transposed form: axis i carries every third bit of it. */
  for (std::uint32_t Q = M; Q > 1; Q >>= 1) {
    std::uint32_t const P = Q - 1;
    for (std::size_t i = 0; i < 3; ++i) {
      if (x[i] & Q)
        x[0] ^= P; //!< invert
      else {
        std::uint32_t const t = (x[0] ^ x[i]) & P;
        x[0] ^= t; //!< exchange
        x[i] ^= t;
      }
    }
  }
  for (std::size_t i = 1; i < 3; ++i)
    x[i] ^= x[i - 1];
  std::uint32_t t = 0;
  for (std::uint32_t Q = M; Q > 1; Q >>= 1)
    if (x[2] & Q)
      t ^= Q - 1;
  for (std::size_t i = 0; i < 3; ++i)
    x[i] ^= t;

  /* Undoing the transpose is the same interleave Morton wants, which is why
     `morton_index` puts x in the high bit of each triple. */
  return morton_index(x, bits);
}

namespace {

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

//! Thread a space-filling curve through the mesh and cut it into runs
class Curve : public Partitioner {
public:
  /*! `hilbert` picks the curve that never jumps, as against the Z-order curve
      that jumps a whole block whenever it crosses a power-of-two boundary. */
  Curve(bool hilbert, char const *name, char const *describe)
      : hilbert_{hilbert}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const v = begin_assign(mesh, opt, z2r);
    auto const &zcoord = zone_centroids(mesh);
    std::vector<int> const wgt = zone_weights(mesh, opt);

    /* Put the centroids on the curve's lattice.  The box is the mesh's own,
       so the resolution goes where the zones are rather than where the origin
       is; an axis with no extent at all quantizes to zero and contributes
       nothing to the ordering, which is what a flat mesh wants. */
    auto const [lo, hi] = centroid_bounds(v, zcoord);
    double const top = static_cast<double>((1u << CURVE_BITS) - 1u);
    std::array<double, 3> scale{0.0, 0.0, 0.0};
    for (std::size_t d = 0; d < 3; ++d)
      if (hi[d] > lo[d])
        scale[d] = top / (hi[d] - lo[d]);

    std::vector<std::pair<std::uint64_t, int>> along;
    along.reserve(v.size());
    for (int const z : v) {
      std::array<std::uint32_t, 3> q{};
      for (std::size_t d = 0; d < 3; ++d) {
        double const t =
            (zcoord[static_cast<std::size_t>(z)][d] - lo[d]) * scale[d];
        q[d] = static_cast<std::uint32_t>(
            std::clamp(std::llround(t), 0LL, static_cast<long long>(top)));
      }
      along.emplace_back(hilbert_ ? hilbert_index(q) : morton_index(q), z);
    }
    /* Lexicographic, so zones sharing a lattice cell are ordered by zone index
       and the cuts do not depend on the sort algorithm. */
    std::ranges::sort(along);

    /* Hand `split_sorted` each zone's *position* along the curve rather than
       its index: a 63-bit index does not survive the trip through a `double`
       key, and only the order matters. */
    std::vector<double> key(static_cast<std::size_t>(mesh.zones.size()), 0.0);
    for (std::size_t i = 0; i < along.size(); ++i)
      key[static_cast<std::size_t>(along[i].second)] = static_cast<double>(i);

    split_sorted(v, key, wgt, 0, opt.nparts, z2r);
  }

private:
  bool hilbert_;
  char const *name_;
  char const *describe_;
};

} // namespace

std::unique_ptr<Partitioner> make_geometric(Method const &m) {
  std::string_view const name{m.name};
  if (name == "block")
    return std::make_unique<Block>(3, m.name, m.describe);
  if (name == "hilbert")
    return std::make_unique<Curve>(true, m.name, m.describe);
  if (name == "morton")
    return std::make_unique<Curve>(false, m.name, m.describe);
  if (name == "pencil")
    return std::make_unique<Block>(2, m.name, m.describe);
  if (name == "rcb")
    return std::make_unique<Recursive>(false, m.name, m.describe);
  if (name == "rib")
    return std::make_unique<Recursive>(true, m.name, m.describe);
  if (name == "slab")
    return std::make_unique<Slab>(-1, m.name, m.describe);
  if (name == "slab-x")
    return std::make_unique<Slab>(0, m.name, m.describe);
  if (name == "slab-y")
    return std::make_unique<Slab>(1, m.name, m.describe);
  if (name == "slab-z")
    return std::make_unique<Slab>(2, m.name, m.describe);
  return nullptr;
}

} // namespace Part
} // namespace Ume
