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
  \file Ume/Brick_Mesh.cc
*/

#include "Ume/Brick_Mesh.hh"
#include "Ume/Comm_Transport.hh"
#include "Ume/mesh_audit.hh"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <numeric>
#include <ostream>
#include <random>
#include <stdexcept>

namespace Ume {
namespace Brick {

using SOA_Idx::Entity;
using SOA_Idx::Mesh;

namespace {

//! An ordered ring of point indices bounding one face of one zone
using Loop = std::vector<int>;

/*! The six Kuhn tetrahedra of a cell are indexed by a permutation of the
    axes: tetrahedron `p` has the vertices
    `{(0,0,0), e_p0, e_p0 + e_p1, (1,1,1)}`.  Every tetrahedron contains the
    body diagonal (0,0,0)-(1,1,1), which makes the induced triangulation of
    each cell face translation invariant: the diagonal of the face always runs
    between local coordinates (0,0) and (1,1) of the two in-face axes, taken in
    ascending order.  Adjacent cells therefore always agree on how a shared
    face is split, and a hex cell can use the same rule for a face it shares
    with a tet cell. */
constexpr int KUHN_PERM[6][3] = {
    {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};

//! Vertices (as indices into a tet's four vertices) of each tet face
constexpr int TET_FACES[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

//! The two axes of a cell face normal to `axis`, in ascending order
constexpr std::array<int, 2> face_axes(int const axis) {
  return axis == 0 ? std::array<int, 2>{1, 2}
      : axis == 1  ? std::array<int, 2>{0, 2}
                   : std::array<int, 2>{0, 1};
}

//! The four corners of a cell face, counterclockwise about the face normal
constexpr int QUAD_UV[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

//! Sorted, -1 padded point list; identifies a face independent of orientation
std::array<int, 4> face_key(Loop const &loop) {
  std::array<int, 4> key{-1, -1, -1, -1};
  std::copy(loop.begin(), loop.end(), key.begin());
  std::sort(key.begin(), key.begin() + static_cast<long>(loop.size()));
  return key;
}

class Builder {
public:
  Builder(Spec const &spec, std::vector<Cell_Type> const &cells)
      : spec_{spec}, cells_{cells} {}

  Stats run(Mesh &mesh);

private:
  /* --- lattice helpers --- */
  int point_id(int const i, int const j, int const k) const {
    return i + (spec_.nx + 1) * (j + (spec_.ny + 1) * k);
  }
  int cell_id(int const i, int const j, int const k) const {
    return i + spec_.nx * (j + spec_.ny * k);
  }
  bool is_tet(int const i, int const j, int const k) const {
    return cells_[static_cast<size_t>(cell_id(i, j, k))] == Cell_Type::TET;
  }
  /*! Is the cell face normal to `axis` at lattice index (i,j,k) triangulated?
      A face is split if either of the cells it bounds is subdivided. */
  bool is_split(int const axis, int const i, int const j, int const k) const;

  /* --- geometry --- */
  Vec3 point(int const p) const { return pcoord_[static_cast<size_t>(p)]; }
  Vec3 loop_centroid(Loop const &loop) const;
  Vec3 loop_normal(Loop const &loop, Vec3 const &center) const;
  //! Reverse `loop` if needed so its normal (right-hand rule) points at `zc`
  void orient(Loop &loop, Vec3 const &zc) const;

  /* --- entity construction --- */
  void make_points();
  std::vector<Loop> hex_loops(int i, int j, int k) const;
  std::vector<Loop> tet_loops(int i, int j, int k, int tet) const;
  int edge_id(int p1, int p2);
  int face_id(Loop const &loop, int owner);
  int add_zone(short mask);
  //! Add corners for `z` at every point of `loops`, ascending by point index
  void add_corners(int z, std::vector<Loop> const &loops, short mask);
  //! Index of the corner of zone `z` at point `p`
  int corner_of(int z, int p) const;
  //! Mean of the points of `z`; this is what `VAR_zcoord` computes
  Vec3 zone_centroid(int z) const;
  int new_side(int z, int p1, int p2, int e, int f, int c1, int c2, short mask);
  void link_s2(int f, int e, int s);
  /*! Emit the sides of zone `z`.  `face_owner` is recorded in `m:f>z1`/`z2`;
      pass -1 to leave the face's zone links alone. */
  void emit_zone(
      int z, std::vector<Loop> const &loops, int face_owner, short side_mask);
  void make_boundary();
  void fill(Mesh &mesh) const;

  Spec spec_;
  std::vector<Cell_Type> const &cells_;

  DS_Types::VEC3V_T pcoord_;
  std::vector<short> pmask_;

  std::vector<int> ep1_, ep2_;
  std::vector<short> emask_;
  std::map<std::array<int, 2>, int> edge_ids_;

  std::vector<int> fz1_, fz2_;
  std::vector<short> fmask_;
  std::vector<Loop> floop_; //!< the face ring, as oriented by `m:f>z1`
  std::map<std::array<int, 4>, int> face_ids_;

  std::vector<short> zmask_;

  std::vector<int> cp_, cz_;
  std::vector<short> cmask_;
  std::vector<int> zbegin_; //!< first corner of each zone (corners are dense)

  std::vector<int> sz_, sp1_, sp2_, se_, sf_, sc1_, sc2_, ss2_, ss3_, ss4_,
      ss5_;
  std::vector<short> smask_;
  //! (face, edge) pairs waiting for their second side
  std::map<std::pair<int, int>, int> open_fe_;

  Stats stats_;
};

bool Builder::is_split(
    int const axis, int const i, int const j, int const k) const {
  int const lo[3] = {i - (axis == 0), j - (axis == 1), k - (axis == 2)};
  int const n[3] = {spec_.nx, spec_.ny, spec_.nz};
  int const hi[3] = {i, j, k};
  bool split{false};
  if (lo[axis] >= 0)
    split = split || is_tet(lo[0], lo[1], lo[2]);
  if (hi[axis] < n[axis])
    split = split || is_tet(hi[0], hi[1], hi[2]);
  return split;
}

Vec3 Builder::loop_centroid(Loop const &loop) const {
  Vec3 c(0.0);
  for (int const p : loop)
    c += point(p);
  return c / static_cast<double>(loop.size());
}

Vec3 Builder::loop_normal(Loop const &loop, Vec3 const &center) const {
  Vec3 n(0.0);
  size_t const len = loop.size();
  for (size_t i = 0; i < len; ++i) {
    Vec3 const a = point(loop[i]) - center;
    Vec3 const b = point(loop[(i + 1) % len]) - center;
    n += crossprod(a, b);
  }
  return n / 2.0;
}

void Builder::orient(Loop &loop, Vec3 const &zc) const {
  Vec3 const c = loop_centroid(loop);
  /* `side_vol` is positive when the ring is traversed counterclockwise as seen
     from inside the zone, i.e. when the right-hand-rule normal of the ring
     points *at* the zone center. */
  if (dotprod(loop_normal(loop, c), zc - c) < 0.0)
    std::reverse(loop.begin(), loop.end());
}

void Builder::make_points() {
  int const np = (spec_.nx + 1) * (spec_.ny + 1) * (spec_.nz + 1);
  pcoord_.resize(static_cast<size_t>(np));
  pmask_.assign(static_cast<size_t>(np), 1);
  for (int k = 0; k <= spec_.nz; ++k)
    for (int j = 0; j <= spec_.ny; ++j)
      for (int i = 0; i <= spec_.nx; ++i) {
        pcoord_[static_cast<size_t>(point_id(i, j, k))] =
            Vec3{{spec_.x0 + i * spec_.dx, spec_.y0 + j * spec_.dy,
                spec_.z0 + k * spec_.dz}};
      }
}

std::vector<Loop> Builder::hex_loops(
    int const i, int const j, int const k) const {
  std::vector<Loop> loops;
  loops.reserve(8);
  int const base[3] = {i, j, k};
  for (int axis = 0; axis < 3; ++axis) {
    auto const [b, c] = face_axes(axis);
    for (int t = 0; t < 2; ++t) {
      int quad[4];
      for (int v = 0; v < 4; ++v) {
        int off[3];
        off[axis] = t;
        off[b] = QUAD_UV[v][0];
        off[c] = QUAD_UV[v][1];
        quad[v] =
            point_id(base[0] + off[0], base[1] + off[1], base[2] + off[2]);
      }
      if (is_split(axis, base[0] + (axis == 0 ? t : 0),
              base[1] + (axis == 1 ? t : 0), base[2] + (axis == 2 ? t : 0))) {
        /* Split along the (0,0)-(1,1) diagonal, matching the Kuhn split. */
        loops.push_back(Loop{quad[0], quad[1], quad[2]});
        loops.push_back(Loop{quad[0], quad[2], quad[3]});
      } else {
        loops.push_back(Loop{quad[0], quad[1], quad[2], quad[3]});
      }
    }
  }
  return loops;
}

std::vector<Loop> Builder::tet_loops(
    int const i, int const j, int const k, int const tet) const {
  int const *const perm = KUHN_PERM[tet];
  int off[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {1, 1, 1}};
  off[1][perm[0]] = 1;
  off[2][perm[0]] = 1;
  off[2][perm[1]] = 1;

  int v[4];
  for (int n = 0; n < 4; ++n)
    v[n] = point_id(i + off[n][0], j + off[n][1], k + off[n][2]);

  std::vector<Loop> loops;
  loops.reserve(4);
  for (auto const &f : TET_FACES)
    loops.push_back(Loop{v[f[0]], v[f[1]], v[f[2]]});
  return loops;
}

int Builder::edge_id(int const p1, int const p2) {
  std::array<int, 2> const key{std::min(p1, p2), std::max(p1, p2)};
  auto const [it, added] =
      edge_ids_.try_emplace(key, static_cast<int>(ep1_.size()));
  if (added) {
    ep1_.push_back(key[0]);
    ep2_.push_back(key[1]);
    emask_.push_back(1);
  }
  return it->second;
}

int Builder::face_id(Loop const &loop, int const owner) {
  auto const [it, added] =
      face_ids_.try_emplace(face_key(loop), static_cast<int>(fz1_.size()));
  int const f = it->second;
  if (added) {
    floop_.push_back(loop);
    fz1_.push_back(owner);
    fz2_.push_back(-1);
    fmask_.push_back(1);
  } else if (owner >= 0) {
    if (fz2_[static_cast<size_t>(f)] != -1)
      throw std::logic_error("brick: a face is shared by more than two zones");
    fz2_[static_cast<size_t>(f)] = owner;
  }
  return f;
}

int Builder::add_zone(short const mask) {
  zmask_.push_back(mask);
  zbegin_.push_back(static_cast<int>(cp_.size()));
  return static_cast<int>(zmask_.size()) - 1;
}

void Builder::add_corners(
    int const z, std::vector<Loop> const &loops, short const mask) {
  std::vector<int> pts;
  for (auto const &loop : loops)
    pts.insert(pts.end(), loop.begin(), loop.end());
  std::sort(pts.begin(), pts.end());
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
  for (int const p : pts) {
    cp_.push_back(p);
    cz_.push_back(z);
    cmask_.push_back(mask);
  }
}

int Builder::corner_of(int const z, int const p) const {
  size_t const first = static_cast<size_t>(zbegin_[static_cast<size_t>(z)]);
  size_t const last = (static_cast<size_t>(z) + 1 < zbegin_.size())
      ? static_cast<size_t>(zbegin_[static_cast<size_t>(z) + 1])
      : cp_.size();
  auto const b = cp_.begin() + static_cast<long>(first);
  auto const e = cp_.begin() + static_cast<long>(last);
  auto const it = std::lower_bound(b, e, p);
  assert(it != e && *it == p);
  return static_cast<int>(it - cp_.begin());
}

Vec3 Builder::zone_centroid(int const z) const {
  size_t const first = static_cast<size_t>(zbegin_[static_cast<size_t>(z)]);
  size_t const last = (static_cast<size_t>(z) + 1 < zbegin_.size())
      ? static_cast<size_t>(zbegin_[static_cast<size_t>(z) + 1])
      : cp_.size();
  Vec3 zc(0.0);
  for (size_t c = first; c < last; ++c)
    zc += point(cp_[c]);
  return zc / static_cast<double>(last - first);
}

int Builder::new_side(int const z, int const p1, int const p2, int const e,
    int const f, int const c1, int const c2, short const mask) {
  sz_.push_back(z);
  sp1_.push_back(p1);
  sp2_.push_back(p2);
  se_.push_back(e);
  sf_.push_back(f);
  sc1_.push_back(c1);
  sc2_.push_back(c2);
  ss2_.push_back(-1);
  ss3_.push_back(-1);
  ss4_.push_back(-1);
  ss5_.push_back(-1);
  smask_.push_back(mask);
  return static_cast<int>(smask_.size()) - 1;
}

void Builder::link_s2(int const f, int const e, int const s) {
  auto const key = std::pair{f, e};
  auto const it = open_fe_.find(key);
  if (it == open_fe_.end()) {
    open_fe_.emplace(key, s);
  } else {
    int const other = it->second;
    ss2_[static_cast<size_t>(other)] = s;
    ss2_[static_cast<size_t>(s)] = other;
    open_fe_.erase(it);
  }
}

void Builder::emit_zone(int const z, std::vector<Loop> const &loops,
    int const face_owner, short const side_mask) {
  //! sides of this zone that contain a given edge; a closed polyhedron has two
  std::map<int, std::vector<int>> edge_sides;

  for (auto const &loop : loops) {
    int const f = face_id(loop, face_owner);
    int const n = static_cast<int>(loop.size());
    int const first = static_cast<int>(smask_.size());
    for (int i = 0; i < n; ++i) {
      int const p1 = loop[static_cast<size_t>(i)];
      int const p2 = loop[static_cast<size_t>((i + 1) % n)];
      int const e = edge_id(p1, p2);
      int const s = new_side(
          z, p1, p2, e, f, corner_of(z, p1), corner_of(z, p2), side_mask);
      /* s3/s4: the next and previous sides around this face of this zone. */
      ss3_[static_cast<size_t>(s)] = first + (i + 1) % n;
      ss4_[static_cast<size_t>(s)] = first + (i + n - 1) % n;
      link_s2(f, e, s);
      edge_sides[e].push_back(s);
    }
  }

  /* s5: the other side of this zone on the same edge, i.e. the one on the
     zone's other face containing that edge.  A degenerate boundary zone has a
     single face, so its sides refer to themselves. */
  for (auto const &[e, sides] : edge_sides) {
    if (sides.size() == 2) {
      ss5_[static_cast<size_t>(sides[0])] = sides[1];
      ss5_[static_cast<size_t>(sides[1])] = sides[0];
    } else if (sides.size() == 1) {
      ss5_[static_cast<size_t>(sides[0])] = sides[0];
    } else {
      throw std::logic_error("brick: an edge appears in more than two faces "
                             "of the same zone");
    }
  }
}

void Builder::make_boundary() {
  /* Whatever is left unpaired in `open_fe_` is a (face, edge) with only one
     zone behind it: the physical mesh boundary.  std::map keeps these ordered
     by face, then edge. */
  std::vector<int> bfaces;
  for (auto const &[fe, s] : open_fe_) {
    (void)s;
    if (bfaces.empty() || bfaces.back() != fe.first)
      bfaces.push_back(fe.first);
  }
  stats_.boundary_faces = static_cast<int>(bfaces.size());

  for (int const f : bfaces) {
    Loop const &ring = floop_[static_cast<size_t>(f)];
    fmask_[static_cast<size_t>(f)] = -1;
    size_t const len = ring.size();
    for (size_t i = 0; i < len; ++i) {
      pmask_[static_cast<size_t>(ring[i])] = -1;
      emask_[static_cast<size_t>(edge_id(ring[i], ring[(i + 1) % len]))] = -1;
    }

    /* The phantom sides traverse the face in the opposite direction from the
       real zone, exactly as the second zone of an interior face would. */
    std::vector<Loop> loops{Loop(ring.rbegin(), ring.rend())};

    if (spec_.boundary == Boundary_Style::BOUNDARY_ZONES) {
      int const z = add_zone(-1);
      add_corners(z, loops, -1);
      emit_zone(z, loops, z, -1);
      stats_.boundary_zones += 1;
    } else {
      emit_zone(fz1_[static_cast<size_t>(f)], loops, -1, -1);
    }
  }

  if (!open_fe_.empty())
    throw std::logic_error("brick: unpaired sides remain after boundary "
                           "construction");
}

void Builder::fill(Mesh &mesh) const {
  mesh.ivtag = UME_VERSION_2;
  mesh.version_header = true;
  mesh.mype = 0;
  mesh.numpe = 1;
  mesh.geo = Mesh::CARTESIAN;
  mesh.dump_iotas = false;

  auto const n = [](auto const &v) { return static_cast<int>(v.size()); };
  auto const internal = [](Entity &e) {
    std::fill(e.comm_type.begin(), e.comm_type.end(), Entity::INTERNAL);
  };

  mesh.points.resize(n(pmask_), n(pmask_), 0);
  mesh.points.mask = pmask_;
  internal(mesh.points);
  mesh.ds->access_vec3v("pcoord") = pcoord_;

  mesh.edges.resize(n(emask_), n(emask_), 0);
  mesh.edges.mask = emask_;
  internal(mesh.edges);
  mesh.ds->access_intv("m:e>p1") = ep1_;
  mesh.ds->access_intv("m:e>p2") = ep2_;

  mesh.faces.resize(n(fmask_), n(fmask_), 0);
  mesh.faces.mask = fmask_;
  internal(mesh.faces);
  mesh.ds->access_intv("m:f>z1") = fz1_;
  mesh.ds->access_intv("m:f>z2") = fz2_;

  mesh.sides.resize(n(smask_), n(smask_), 0);
  mesh.sides.mask = smask_;
  internal(mesh.sides);
  mesh.ds->access_intv("m:s>z") = sz_;
  mesh.ds->access_intv("m:s>p1") = sp1_;
  mesh.ds->access_intv("m:s>p2") = sp2_;
  mesh.ds->access_intv("m:s>e") = se_;
  mesh.ds->access_intv("m:s>f") = sf_;
  mesh.ds->access_intv("m:s>c1") = sc1_;
  mesh.ds->access_intv("m:s>c2") = sc2_;
  mesh.ds->access_intv("m:s>s2") = ss2_;
  mesh.ds->access_intv("m:s>s3") = ss3_;
  mesh.ds->access_intv("m:s>s4") = ss4_;
  mesh.ds->access_intv("m:s>s5") = ss5_;

  mesh.corners.resize(n(cmask_), n(cmask_), 0);
  mesh.corners.mask = cmask_;
  internal(mesh.corners);
  mesh.ds->access_intv("m:c>p") = cp_;
  mesh.ds->access_intv("m:c>z") = cz_;

  mesh.zones.resize(n(zmask_), n(zmask_), 0);
  mesh.zones.mask = zmask_;
  internal(mesh.zones);
}

Stats Builder::run(Mesh &mesh) {
  make_points();

  stats_.cells = spec_.num_cells();
  for (int k = 0; k < spec_.nz; ++k)
    for (int j = 0; j < spec_.ny; ++j)
      for (int i = 0; i < spec_.nx; ++i) {
        std::vector<std::vector<Loop>> cell_zones;
        if (is_tet(i, j, k)) {
          stats_.tet_cells += 1;
          for (int t = 0; t < 6; ++t)
            cell_zones.push_back(tet_loops(i, j, k, t));
          stats_.tet_zones += 6;
        } else {
          cell_zones.push_back(hex_loops(i, j, k));
          stats_.split_faces += static_cast<int>(cell_zones[0].size()) - 6;
          stats_.hex_zones += 1;
        }
        for (auto &loops : cell_zones) {
          int const z = add_zone(1);
          add_corners(z, loops, 1);
          Vec3 const zc = zone_centroid(z);
          for (auto &loop : loops)
            orient(loop, zc);
          emit_zone(z, loops, z, 1);
        }
      }

  stats_.real_sides = static_cast<int>(smask_.size());
  make_boundary();
  stats_.phantom_sides = static_cast<int>(smask_.size()) - stats_.real_sides;

  fill(mesh);

  stats_.points = static_cast<int>(pmask_.size());
  stats_.edges = static_cast<int>(emask_.size());
  stats_.faces = static_cast<int>(fmask_.size());
  stats_.corners = static_cast<int>(cmask_.size());
  return stats_;
}

//! Number of face-connected components of the tet-cell set
int count_clusters(Spec const &spec, std::vector<Cell_Type> const &cells) {
  int const n = spec.num_cells();
  std::vector<int> parent(static_cast<size_t>(n));
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&parent](int a) {
    while (parent[static_cast<size_t>(a)] != a) {
      parent[static_cast<size_t>(a)] =
          parent[static_cast<size_t>(parent[static_cast<size_t>(a)])];
      a = parent[static_cast<size_t>(a)];
    }
    return a;
  };
  auto join = [&](int a, int b) {
    parent[static_cast<size_t>(find(a))] = find(b);
  };

  auto tet = [&](int i, int j, int k) {
    return cells[static_cast<size_t>(i + spec.nx * (j + spec.ny * k))] ==
        Cell_Type::TET;
  };
  int count{0};
  for (int k = 0; k < spec.nz; ++k)
    for (int j = 0; j < spec.ny; ++j)
      for (int i = 0; i < spec.nx; ++i) {
        if (!tet(i, j, k))
          continue;
        ++count;
        int const c = i + spec.nx * (j + spec.ny * k);
        if (i > 0 && tet(i - 1, j, k))
          join(c, c - 1);
        if (j > 0 && tet(i, j - 1, k))
          join(c, c - spec.nx);
        if (k > 0 && tet(i, j, k - 1))
          join(c, c - spec.nx * spec.ny);
      }
  if (count == 0)
    return 0;
  std::vector<int> roots;
  for (int c = 0; c < n; ++c)
    if (cells[static_cast<size_t>(c)] == Cell_Type::TET)
      roots.push_back(find(c));
  std::sort(roots.begin(), roots.end());
  return static_cast<int>(
      std::unique(roots.begin(), roots.end()) - roots.begin());
}

} // namespace

double Stats::tet_cell_fraction() const {
  return cells ? static_cast<double>(tet_cells) / cells : 0.0;
}

double Stats::tet_zone_fraction() const {
  int const z = hex_zones + tet_zones;
  return z ? static_cast<double>(tet_zones) / z : 0.0;
}

void Stats::print(std::ostream &os) const {
  os << "\tCells:           " << cells << " (" << tet_cells << " tet in "
     << tet_clusters << " cluster(s), " << cells - tet_cells << " hex)\n";
  os << "\tTet fraction:    " << tet_cell_fraction() << " of cells, "
     << tet_zone_fraction() << " of zones\n";
  os << "\tZones:           " << hex_zones + tet_zones << " real (" << hex_zones
     << " hex, " << tet_zones << " tet), " << boundary_zones << " boundary\n";
  os << "\tSplit hex faces: " << split_faces << '\n';
  os << "\tPoints:          " << points << '\n';
  os << "\tEdges:           " << edges << '\n';
  os << "\tFaces:           " << faces << " (" << boundary_faces
     << " on the boundary)\n";
  os << "\tCorners:         " << corners << '\n';
  os << "\tSides:           " << real_sides << " real, " << phantom_sides
     << " phantom\n";
}

std::vector<Cell_Type> select_cells(Spec const &spec) {
  int const ncells = spec.num_cells();
  if (ncells < 1)
    throw std::invalid_argument("brick: nx, ny and nz must all be positive");
  std::vector<Cell_Type> cells(static_cast<size_t>(ncells), Cell_Type::HEX);

  double f = std::clamp(spec.tet_fraction, 0.0, 1.0);
  if (spec.basis == Fraction_Basis::ZONES) {
    /* A subdivided cell contributes 6 zones instead of 1, so a zone fraction
       t corresponds to a cell fraction t / (6 - 5t). */
    f = (f >= 1.0) ? 1.0 : f / (6.0 - 5.0 * f);
  }
  int const target =
      std::clamp(static_cast<int>(std::llround(f * ncells)), 0, ncells);
  if (target == 0)
    return cells;

  double const lambda = std::clamp(spec.clustering, 0.0, 1.0);
  int const nseeds = std::clamp(
      static_cast<int>(
          std::llround(std::pow(static_cast<double>(target), 1.0 - lambda))),
      1, target);

  std::mt19937_64 rng{spec.seed};
  std::vector<int> order(static_cast<size_t>(ncells));
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), rng);

  std::vector<char> in_frontier(static_cast<size_t>(ncells), 0);
  std::vector<int> frontier;
  int marked{0};

  auto neighbors = [&](int const c, int out[6]) {
    int const i = c % spec.nx;
    int const j = (c / spec.nx) % spec.ny;
    int const k = c / (spec.nx * spec.ny);
    int n{0};
    if (i > 0)
      out[n++] = c - 1;
    if (i + 1 < spec.nx)
      out[n++] = c + 1;
    if (j > 0)
      out[n++] = c - spec.nx;
    if (j + 1 < spec.ny)
      out[n++] = c + spec.nx;
    if (k > 0)
      out[n++] = c - spec.nx * spec.ny;
    if (k + 1 < spec.nz)
      out[n++] = c + spec.nx * spec.ny;
    return n;
  };

  auto mark = [&](int const c) {
    cells[static_cast<size_t>(c)] = Cell_Type::TET;
    ++marked;
    int nbr[6];
    int const cnt = neighbors(c, nbr);
    for (int m = 0; m < cnt; ++m) {
      size_t const u = static_cast<size_t>(nbr[m]);
      if (cells[u] == Cell_Type::HEX && !in_frontier[u]) {
        in_frontier[u] = 1;
        frontier.push_back(nbr[m]);
      }
    }
  };

  /* Plant the seeds, then grow them simultaneously (an Eden-style cluster
     growth) until the requested number of cells is subdivided. */
  for (int s = 0; s < nseeds && marked < target; ++s) {
    int const c = order[static_cast<size_t>(s)];
    if (cells[static_cast<size_t>(c)] == Cell_Type::HEX)
      mark(c);
  }
  while (marked < target && !frontier.empty()) {
    size_t const pick = static_cast<size_t>(rng() % frontier.size());
    int const c = frontier[pick];
    frontier[pick] = frontier.back();
    frontier.pop_back();
    in_frontier[static_cast<size_t>(c)] = 0;
    if (cells[static_cast<size_t>(c)] == Cell_Type::HEX)
      mark(c);
  }
  return cells;
}

Stats build(Spec const &spec, std::vector<Cell_Type> const &cells, Mesh &mesh) {
  if (static_cast<int>(cells.size()) != spec.num_cells())
    throw std::invalid_argument("brick: cell list has the wrong length");
  Builder builder{spec, cells};
  Stats stats = builder.run(mesh);
  stats.tet_clusters = count_clusters(spec, cells);
  return stats;
}

bool validate(Mesh &mesh, Spec const &spec, std::ostream &os) {
  /* Everything that is true of any mesh lives in the shared audit; what is
     left here needs to know that this mesh is a brick of known cell size. */
  SOA_Idx::Audit_Options aopts;
  aopts.census = true;
  bool ok = SOA_Idx::audit(mesh, os, aopts).ok();

  auto fail = [&os, &ok](char const *const what) -> std::ostream & {
    ok = false;
    return os << "  [" << what << "] ";
  };

  auto const &z2c = mesh.ds->caccess_intrr("m:z>c");
  auto const &corner_vol = mesh.ds->caccess_dblv("corner_vol");
  auto const &point_norm = mesh.ds->caccess_vec3v("point_norm");
  auto const &zmask = mesh.zones.mask;
  auto const &pmask = mesh.points.mask;

  int const pll = mesh.points.size();
  int const zll = mesh.zones.size();

  double const cell_vol = spec.dx * spec.dy * spec.dz;
  double const eps = 1e-11 * cell_vol;

  /* Each zone is a whole cell (8 corners) or one of six Kuhn tetrahedra (4
     corners), whatever its faces were split into. */
  double total_vol{0.0};
  bool named_a_zone{false};
  for (int z = 0; z < zll; ++z) {
    if (zmask[z] < 1)
      continue;
    double vol{0.0};
    int ncorners{0};
    for (int const c : z2c[z]) {
      vol += corner_vol[c];
      ncorners += 1;
    }
    total_vol += vol;
    double const want = (ncorners == 4) ? cell_vol / 6.0 : cell_vol;
    if ((ncorners != 4 && ncorners != 8) || std::abs(vol - want) > eps) {
      /* Only the first bad zone is named -- a systematic error would name
         every zone -- but the sum runs to the end regardless.  The
         whole-mesh check below is a different question, and answering it
         from a partial sum reported a second failure that said nothing:
         a mesh whose volume is exactly right would be called wrong. */
      if (!named_a_zone) {
        fail("zone volume") << "zone " << z << " has " << ncorners
                            << " corners and volume " << vol << ", expected "
                            << want << '\n';
        named_a_zone = true;
      }
    }
  }
  double const want_vol = spec.num_cells() * cell_vol;
  if (std::abs(total_vol - want_vol) > 1e-10 * want_vol)
    fail("mesh volume") << total_vol << ", expected " << want_vol << '\n';

  /* Point masks follow the lattice, and a point in the interior of a face of
     the brick has exactly the outward axis normal. */
  int const np[3] = {spec.nx, spec.ny, spec.nz};
  for (int p = 0; p < pll; ++p) {
    int const ijk[3] = {p % (spec.nx + 1), (p / (spec.nx + 1)) % (spec.ny + 1),
        p / ((spec.nx + 1) * (spec.ny + 1))};
    int nfaces{0}, axis{-1}, side{0};
    for (int d = 0; d < 3; ++d) {
      if (ijk[d] == 0) {
        ++nfaces;
        axis = d;
        side = -1;
      } else if (ijk[d] == np[d]) {
        ++nfaces;
        axis = d;
        side = 1;
      }
    }
    if ((nfaces > 0) != (pmask[p] < 0)) {
      fail("point mask") << "point " << p << " has mask " << pmask[p]
                         << " but lies on " << nfaces << " brick face(s)\n";
      break;
    }
    if (nfaces == 1) {
      Vec3 want(0.0);
      want[static_cast<size_t>(axis)] = side;
      if (vectormag(point_norm[p] - want) > 1e-12) {
        fail("point_norm") << "point " << p << " normal " << point_norm[p]
                           << ", expected " << want << '\n';
        break;
      }
    }
  }

  return ok;
}

} // namespace Brick
} // namespace Ume
