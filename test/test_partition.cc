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
  \file test_partition.cc

  Tests for the partitioning strategies.  What is checked here is the contract
  every method has to keep -- that it covers exactly the real zones, stays
  inside `[0, nparts)`, and balances the weight it was given -- rather than the
  geometry of any one method, so the sweep picks its methods up from
  `Part::methods()` and grows as methods land.
*/

#include "Ume/Brick_Mesh.hh"
#include "Ume/Decompose.hh"
#include "Ume/Global_Ids.hh"
#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include "Ume/Partition_Metrics.hh"
#include "Ume/mesh_audit.hh"
#include "Ume/partition_check.hh"
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace Ume;
using Ume::SOA_Idx::Mesh;

namespace {

//! Build a brick, ready to partition
Brick::Stats brick(Brick::Spec const &spec, Mesh &mesh) {
  Ume::attach_serial_transport(mesh);
  return Brick::build(spec, mesh);
}

Brick::Spec mixed_spec() {
  Brick::Spec spec;
  spec.nx = 4;
  spec.ny = 3;
  spec.nz = 2;
  spec.tet_fraction = 0.5;
  return spec;
}

//! The total weight each part was given
std::vector<long> part_weights(Part::Zone_Map const &z2r,
    std::vector<int> const &wgt, int nparts) {
  std::vector<long> w(static_cast<std::size_t>(nparts), 0);
  for (std::size_t z = 0; z < z2r.size(); ++z)
    if (z2r[z] >= 0)
      w[static_cast<std::size_t>(z2r[z])] += wgt[z];
  return w;
}

} // namespace

TEST_CASE("partition: the method table distinguishes unknown from unavailable",
    "[partition]") {
  auto const names = Part::methods();
  CHECK(!names.empty());

  for (std::string const &n : names) {
    INFO(n);
    Part::Method const *const m = Part::find_method(n);
    REQUIRE(m != nullptr);
    CHECK(m->available);
    CHECK(m->describe != nullptr);
    auto const p = Part::make(n);
    REQUIRE(p != nullptr);
    CHECK(n == p->name());
  }

  /* An unknown name and an unavailable one both fail to construct, but the
     table tells them apart -- which is what lets the driver say "requires
     METIS" rather than "unknown method". */
  CHECK(Part::find_method("no-such-method") == nullptr);
  CHECK(Part::make("no-such-method") == nullptr);

  for (Part::Method const &m : Part::all_methods())
    if (!m.available) {
      INFO(m.name);
      CHECK(m.requires_pkg != nullptr);
      CHECK(Part::make(m.name) == nullptr);
      CHECK(Part::find_method(m.name) != nullptr);
    }
}

TEST_CASE("partition: every method covers the real zones exactly",
    "[partition]") {
  Mesh mesh;
  brick(mixed_spec(), mesh);
  auto const &zmask = mesh.zones.mask;

  Part::Options opt;
  Part::Zone_Map z2r;
  for (std::string const &name : Part::methods())
    for (int nparts : {1, 2, 3, 4, 7, 8}) {
      INFO(name << " nparts=" << nparts);
      opt.nparts = nparts;
      Part::make(name)->assign(mesh, opt, z2r);

      REQUIRE(z2r.size() == static_cast<std::size_t>(mesh.zones.size()));
      std::vector<int> used(static_cast<std::size_t>(nparts), 0);
      for (int z = 0; z < mesh.zones.size(); ++z) {
        if (zmask[z] >= 1) {
          REQUIRE(z2r[static_cast<std::size_t>(z)] >= 0);
          REQUIRE(z2r[static_cast<std::size_t>(z)] < nparts);
          used[static_cast<std::size_t>(z2r[static_cast<std::size_t>(z)])] += 1;
        } else {
          /* Boundary zones follow their real zone, which is the
             decomposition's business rather than the strategy's. */
          CHECK(z2r[static_cast<std::size_t>(z)] == -1);
        }
      }
      /* There are far more zones than parts here, so nothing may be empty. */
      CHECK(std::ranges::none_of(used, [](int n) { return n == 0; }));
    }
}

TEST_CASE("partition: the split balances the weight it was given",
    "[partition]") {
  Mesh mesh;
  brick(mixed_spec(), mesh);

  Part::Options opt;
  Part::Zone_Map z2r;
  for (auto vwgt : {Part::Vertex_Weight::UNIT, Part::Vertex_Weight::CORNERS,
           Part::Vertex_Weight::SIDES}) {
    opt.vwgt = vwgt;
    std::vector<int> const wgt = Part::zone_weights(mesh, opt);
    long const total = std::accumulate(wgt.begin(), wgt.end(), 0L);
    REQUIRE(total > 0);

    for (std::string const &name : Part::methods())
      for (int nparts : {2, 3, 4, 8}) {
        INFO(name << " nparts=" << nparts << " vwgt="
                  << static_cast<int>(vwgt));
        opt.nparts = nparts;
        Part::make(name)->assign(mesh, opt, z2r);

        auto const w = part_weights(z2r, wgt, nparts);
        CHECK(std::accumulate(w.begin(), w.end(), 0L) == total);

        /* One zone's worth of slop per part: the cut can only fall between
           zones, and the heaviest zone here is a hex with eight corners. */
        int const heaviest = *std::ranges::max_element(wgt);
        double const mean = static_cast<double>(total) / nparts;
        CHECK(static_cast<double>(*std::ranges::max_element(w)) <=
            mean + heaviest);
        CHECK(static_cast<double>(*std::ranges::min_element(w)) >=
            mean - heaviest);
      }
  }
}

TEST_CASE("partition: a slab is contiguous along its axis", "[partition]") {
  Mesh mesh;
  brick(mixed_spec(), mesh);
  auto const &zcoord = Part::zone_centroids(mesh);
  auto const &zmask = mesh.zones.mask;

  Part::Options opt;
  opt.nparts = 3;
  Part::Zone_Map z2r;

  for (std::size_t axis = 0; axis < 3; ++axis) {
    char const *const name = (axis == 0) ? "slab-x" : (axis == 1) ? "slab-y"
                                                                 : "slab-z";
    INFO(name);
    Part::make(name)->assign(mesh, opt, z2r);

    /* Part numbers rise monotonically along the axis, so each part is a
       contiguous stack of zones.  The test has to order ties the same way
       `split_sorted` does: several zones of one Kuhn-subdivided cell share a
       centroid coordinate exactly, and such a tie group can and does straddle
       a cut, so comparing bare coordinates would be comparing overlapping
       closed intervals. */
    std::vector<int> order = Part::real_zones(mesh);
    std::ranges::sort(order, [&zcoord, axis](int a, int b) {
      double const ka = zcoord[static_cast<std::size_t>(a)][axis];
      double const kb = zcoord[static_cast<std::size_t>(b)][axis];
      return ka != kb ? ka < kb : a < b;
    });
    CHECK(std::ranges::is_sorted(order, {}, [&z2r](int z) {
      return z2r[static_cast<std::size_t>(z)];
    }));
  }
}

TEST_CASE("partition: the grid shape follows the mesh extents", "[partition]") {
  using Dims = std::array<int, 3>;
  using Extent = std::array<double, 3>;

  /* Whatever the shape and however many axes may be cut, the grid has to
     account for every part: a `dims` whose product is less than `nparts`
     silently produces empty ranks. */
  for (int n : {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 17, 24, 36, 64, 1000})
    for (Extent const &ext : {Extent{1.0, 1.0, 1.0}, Extent{4.0, 2.0, 1.0},
             Extent{0.0, 3.0, 0.0}, Extent{0.0, 0.0, 0.0}})
      for (int naxes : {1, 2, 3}) {
        INFO("nparts=" << n << " naxes=" << naxes);
        Dims const d = Part::dims_create(n, ext, naxes);
        CHECK(d[0] * d[1] * d[2] == n);
        CHECK(std::ranges::all_of(d, [](int k) { return k >= 1; }));
      }

  /* A cube divides evenly; a long mesh divides along its length; a prime
     count has nowhere to go but one axis. */
  CHECK(Part::dims_create(8, Extent{1.0, 1.0, 1.0}) == Dims{2, 2, 2});
  CHECK(Part::dims_create(4, Extent{4.0, 1.0, 1.0}) == Dims{4, 1, 1});
  CHECK(Part::dims_create(7, Extent{1.0, 1.0, 1.0}) == Dims{7, 1, 1});
  CHECK(Part::dims_create(12, Extent{4.0, 2.0, 1.0}) == Dims{6, 2, 1});

  /* A single-cell-thick mesh has no extent to cut along its third axis, and
     cutting it there would buy communication for nothing. */
  CHECK(Part::dims_create(8, Extent{3.0, 3.0, 0.0}) == Dims{4, 2, 1});

  /* `pencil` is `block` less its shortest axis, and the shortest axis is the
     one with the least extent rather than the last one. */
  CHECK(Part::dims_create(4, Extent{1.0, 2.0, 3.0}, 2) == Dims{1, 2, 2});
  CHECK(Part::dims_create(4, Extent{1.0, 2.0, 3.0}, 1) == Dims{1, 1, 4});
}

TEST_CASE("partition: a pencil spans the shortest axis", "[partition]") {
  /* 4x4x2 divided four ways puts both cuts on the long axes, and both of them
     fall between coordinate planes rather than inside one, so the columns are
     exact: every zone of a column belongs to one part. */
  Brick::Spec spec;
  spec.nx = 4;
  spec.ny = 4;
  spec.nz = 2;
  Mesh mesh;
  brick(spec, mesh);
  auto const &zcoord = Part::zone_centroids(mesh);

  Part::Options opt;
  opt.nparts = 4;
  Part::Zone_Map z2r;
  Part::make("pencil")->assign(mesh, opt, z2r);

  std::map<std::pair<double, double>, int> column;
  for (int const z : Part::real_zones(mesh)) {
    auto const u = static_cast<std::size_t>(z);
    std::pair<double, double> const xy{zcoord[u][0], zcoord[u][1]};
    auto const it = column.emplace(xy, z2r[u]).first;
    INFO("zone " << z << " at (" << xy.first << ", " << xy.second << ')');
    CHECK(it->second == z2r[u]);
  }
  CHECK(column.size() == 16);
}

TEST_CASE("partition: boxes cost less to communicate than slabs",
    "[partition]") {
  /* This is the reason `block` exists.  Four slabs of a 4x4x4 brick cut three
     whole 4x4 planes; four boxes cut one plane and then half of another, so
     the surface they share is smaller -- and it is the surface, not the
     method, that ends up on the wire. */
  constexpr int NPARTS = 4;
  Brick::Spec spec;
  spec.nx = 4;
  spec.ny = 4;
  spec.nz = 4;
  Mesh mesh;
  brick(spec, mesh);

  Part::Options opt;
  opt.nparts = NPARTS;
  auto measure = [&](char const *method) {
    Part::Zone_Map z2r;
    Part::make(method)->assign(mesh, opt, z2r);
    return Part::measure(
        mesh, Part::build_decomposition(mesh, std::move(z2r), NPARTS));
  };
  Part::Metrics const box = measure("block");
  Part::Metrics const slab = measure("slab-x");

  CHECK(box.cut_faces < slab.cut_faces);
  CHECK(box.kinds[Part::POINT].duplicated < slab.kinds[Part::POINT].duplicated);
  CHECK(box.total_comm_elems() < slab.total_comm_elems());
  CHECK(box.max_surface_to_volume < slab.max_surface_to_volume);

  /* Cheaper communication is only worth having if the parts are still even:
     both methods cut at weight boundaries, so a brick of identical hexes
     divides exactly. */
  CHECK(box.kinds[Part::ZONE].imbalance == 1.0);
}

TEST_CASE("partition: nparts must be positive", "[partition]") {
  Mesh mesh;
  brick(mixed_spec(), mesh);
  Part::Options opt;
  Part::Zone_Map z2r;
  for (int bad : {0, -1}) {
    opt.nparts = bad;
    CHECK_THROWS_AS(Part::make("linear")->assign(mesh, opt, z2r),
        std::invalid_argument);
  }
}

/* --------------------------- the decomposition ---------------------------- */

namespace {

//! The decompositions the layout sweep runs over
struct Case {
  Brick::Spec spec;
  std::string method;
  int nparts;
  std::string describe() const {
    std::ostringstream os;
    os << spec.nx << 'x' << spec.ny << 'x' << spec.nz
       << " tets=" << spec.tet_fraction << ' ' << method << " nparts="
       << nparts;
    return os.str();
  }
};

std::vector<Case> const &cases() {
  static std::vector<Case> const v = [] {
    std::vector<Case> c;
    for (auto const &spec : {[] {
           Brick::Spec s;
           s.nx = 4;
           s.ny = 3;
           s.nz = 2;
           return s;
         }(),
             [] {
               Brick::Spec s;
               s.nx = 3;
               s.ny = 3;
               s.nz = 2;
               s.tet_fraction = 0.5;
               return s;
             }(),
             [] {
               Brick::Spec s;
               s.nx = 3;
               s.ny = 2;
               s.nz = 2;
               s.tet_fraction = 1.0;
               s.boundary = Brick::Boundary_Style::PHANTOM_SIDES_ONLY;
               return s;
             }()})
      /* Every method this build can run goes through the audit, the cross-rank
         check and the metrics, so a method joins the sweep by being registered
         rather than by being remembered. */
      for (std::string const &method : Part::methods())
        for (int n : {1, 2, 3, 4, 7, 8})
          c.push_back({spec, method, n});
    return c;
  }();
  return v;
}

Part::Decomposition decompose(Mesh &mesh, Case const &c,
    Part::Decompose_Options const &dopt = {}) {
  brick(c.spec, mesh);
  Part::Options opt;
  opt.nparts = c.nparts;
  Part::Zone_Map z2r;
  Part::make(c.method)->assign(mesh, opt, z2r);
  return Part::build_decomposition(mesh, std::move(z2r), c.nparts, dopt);
}

//! Check that one rank's list is strictly ascending
bool ascending(std::vector<int> const &v, int lo, int hi) {
  return std::is_sorted(v.begin() + lo, v.begin() + hi) &&
      std::adjacent_find(v.begin() + lo, v.begin() + hi) == v.begin() + hi;
}

} // namespace

TEST_CASE("decompose: the layouts cover the mesh exactly once", "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);

    struct Kind {
      char const *name;
      Part::Entity_Layout const *layout;
      SOA_Idx::Entity const *entity;
      bool exclusive; //!< strictly partitioned rather than duplicated on cuts
    };
    Kind const kinds[] = {{"point", &d.points, &mesh.points, false},
        {"edge", &d.edges, &mesh.edges, false},
        {"face", &d.faces, &mesh.faces, false},
        {"zone", &d.zones, &mesh.zones, true},
        {"corner", &d.corners, &mesh.corners, true},
        {"side", &d.sides, &mesh.sides, true},
        {"iota", &d.iotas, &mesh.iotas, true}};

    for (Kind const &k : kinds) {
      INFO(k.name);
      REQUIRE(k.layout->nparts() == c.nparts);
      REQUIRE(k.layout->nglobal() == k.entity->size());

      std::vector<int> copies(static_cast<std::size_t>(k.entity->size()), 0);
      for (int r = 0; r < c.nparts; ++r) {
        int const lo = k.layout->start[static_cast<std::size_t>(r)];
        int const hi = k.layout->start[static_cast<std::size_t>(r) + 1];
        CHECK(ascending(k.layout->gid, lo, hi));
        for (int i = lo; i < hi; ++i)
          copies[static_cast<std::size_t>(
              k.layout->gid[static_cast<std::size_t>(i)])] += 1;
      }

      for (int g = 0; g < k.entity->size(); ++g) {
        /* Every entity of the serial mesh is somewhere, exactly once if it is
           strictly partitioned, and the SOURCE is one of the ranks it is
           on. */
        int const n = copies[static_cast<std::size_t>(g)];
        CHECK(n >= 1);
        if (k.exclusive)
          CHECK(n == 1);
        int const owner = k.layout->owner[static_cast<std::size_t>(g)];
        REQUIRE(owner >= 0);
        REQUIRE(owner < c.nparts);
        CHECK(k.layout->local_index(owner, g) < k.layout->local_size(owner));
        CHECK(k.layout->local_index(owner, g) >= 0);
      }
    }
  }
}

TEST_CASE("decompose: one part reproduces the input numbering", "[partition]") {
  for (Case const &c : cases()) {
    if (c.nparts != 1)
      continue;
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    for (Part::Entity_Layout const *l : d.layouts()) {
      REQUIRE(l->ghost_local_size(0) == 0);
      REQUIRE(l->local_size(0) == l->nglobal());
      for (int g = 0; g < l->nglobal(); ++g) {
        REQUIRE(l->local_index(0, g) == g);
        REQUIRE(l->owner[static_cast<std::size_t>(g)] == 0);
      }
    }
  }
}

TEST_CASE("decompose: local and ghost indices round trip", "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    for (Part::Entity_Layout const *l : d.layouts())
      for (int r = 0; r < c.nparts; ++r) {
        int const glo = l->gstart[static_cast<std::size_t>(r)];
        int const ghi = l->gstart[static_cast<std::size_t>(r) + 1];
        CHECK(ascending(l->ggid, glo, ghi));
        for (int i = 0; i < l->size(r); ++i)
          REQUIRE(l->local_index(r, l->global_index(r, i)) == i);

        /* Nothing is both local and a ghost here.  A shared point that this
           rank holds but does not source is *not* a ghost -- it is a COPY in
           the real index range -- so the two lists have to be compared
           directly rather than through `owner`. */
        std::vector<int> both;
        std::ranges::set_intersection(
            std::span{l->gid}.subspan(
                static_cast<std::size_t>(l->start[static_cast<std::size_t>(r)]),
                static_cast<std::size_t>(l->local_size(r))),
            std::span{l->ggid}.subspan(static_cast<std::size_t>(glo),
                static_cast<std::size_t>(ghi - glo)),
            std::back_inserter(both));
        CHECK(both.empty());
      }
  }
}

TEST_CASE("decompose: ghosts name only entities the rank already holds",
    "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    auto const &s2zm = mesh.ds->caccess_intv("m:s>z");
    auto const &s2f = mesh.ds->caccess_intv("m:s>f");
    auto const &s2e = mesh.ds->caccess_intv("m:s>e");
    auto const &s2p1 = mesh.ds->caccess_intv("m:s>p1");
    auto const &s2p2 = mesh.ds->caccess_intv("m:s>p2");
    auto const &s2c1 = mesh.ds->caccess_intv("m:s>c1");
    auto const &s2c2 = mesh.ds->caccess_intv("m:s>c2");
    auto const &s2s3 = mesh.ds->caccess_intv("m:s>s3");
    auto const &s2s4 = mesh.ds->caccess_intv("m:s>s4");
    auto const &s2s2 = mesh.ds->caccess_intv("m:s>s2");
    auto const &c2z = mesh.ds->caccess_intv("m:c>z");
    auto const &c2p = mesh.ds->caccess_intv("m:c>p");

    for (int r = 0; r < c.nparts; ++r) {
      auto here = [r](Part::Entity_Layout const &l, int g) {
        return l.local_index(r, g) >= 0;
      };

      for (int i = d.zones.gstart[static_cast<std::size_t>(r)];
          i < d.zones.gstart[static_cast<std::size_t>(r) + 1]; ++i) {
        int const z = d.zones.ggid[static_cast<std::size_t>(i)];
        /* A ghost is a real zone that lives somewhere else. */
        CHECK(mesh.zones.mask[z] >= 1);
        CHECK(d.zones.owner[static_cast<std::size_t>(z)] != r);
      }

      for (int i = d.sides.gstart[static_cast<std::size_t>(r)];
          i < d.sides.gstart[static_cast<std::size_t>(r) + 1]; ++i) {
        int const s = d.sides.ggid[static_cast<std::size_t>(i)];
        INFO("rank " << r << " ghost side " << s);
        CHECK(here(d.zones, s2zm[static_cast<std::size_t>(s)]));
        /* Everything an imported side names has to be here already: there are
           no ghost points, edges or faces to fall back on. */
        CHECK(here(d.faces, s2f[static_cast<std::size_t>(s)]));
        CHECK(here(d.edges, s2e[static_cast<std::size_t>(s)]));
        CHECK(here(d.points, s2p1[static_cast<std::size_t>(s)]));
        CHECK(here(d.points, s2p2[static_cast<std::size_t>(s)]));
        CHECK(here(d.corners, s2c1[static_cast<std::size_t>(s)]));
        CHECK(here(d.corners, s2c2[static_cast<std::size_t>(s)]));
        /* `m:s>s2` and the `s3`/`s4` ring have to resolve; only `m:s>s5`,
           which leaves the face, is allowed to dangle. */
        CHECK(here(d.sides, s2s2[static_cast<std::size_t>(s)]));
        CHECK(here(d.sides, s2s3[static_cast<std::size_t>(s)]));
        CHECK(here(d.sides, s2s4[static_cast<std::size_t>(s)]));
      }

      for (int i = d.corners.gstart[static_cast<std::size_t>(r)];
          i < d.corners.gstart[static_cast<std::size_t>(r) + 1]; ++i) {
        int const cc = d.corners.ggid[static_cast<std::size_t>(i)];
        CHECK(here(d.zones, c2z[static_cast<std::size_t>(cc)]));
        CHECK(d.zones.owner[static_cast<std::size_t>(
                  c2z[static_cast<std::size_t>(cc)])] != r);
        CHECK(here(d.points, c2p[static_cast<std::size_t>(cc)]));
      }
    }
  }
}

namespace {

//! The number of (face, edge) facets that do not have both their sides
/*! This is the audit invariant a cut face is most likely to break -- it is
    what the R3 fix in `mesh_audit.cc` exists to state correctly -- and it is
    decided entirely by which sides get imported, so it can be checked before
    any mesh is built out of the layouts. */
long lonely_facets(Mesh const &mesh, Part::Decomposition const &d) {
  auto const &s2f = mesh.ds->caccess_intv("m:s>f");
  auto const &s2e = mesh.ds->caccess_intv("m:s>e");
  long lonely = 0;
  for (int r = 0; r < d.nparts; ++r) {
    std::vector<std::pair<int, int>> fe;
    for (int i = 0; i < d.sides.size(r); ++i) {
      int const s = d.sides.global_index(r, i);
      fe.emplace_back(s2f[static_cast<std::size_t>(s)],
          s2e[static_cast<std::size_t>(s)]);
    }
    std::ranges::sort(fe);
    for (std::size_t i = 0; i < fe.size();) {
      std::size_t j = i;
      while (j < fe.size() && fe[j] == fe[i])
        ++j;
      if (j - i != 2)
        lonely += 1;
      i = j;
    }
  }
  return lonely;
}

} // namespace

TEST_CASE("decompose: every (face, edge) keeps both its sides", "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    CHECK(lonely_facets(mesh, d) == 0);
  }
}

TEST_CASE("decompose: without ghosts a cut face loses a side", "[partition]") {
  /* Ghost_Depth::NONE is offered as a diagnostic, and this is what it costs:
     the invariant above stops holding the moment a face is cut, which is the
     whole reason a decomposition carries ghost zones at all. */
  Case c{mixed_spec(), "slab-x", 3};
  Mesh mesh;
  Part::Decompose_Options dopt;
  dopt.ghosts = Part::Ghost_Depth::NONE;
  Part::Decomposition const d = decompose(mesh, c, dopt);
  for (Part::Entity_Layout const *l : d.layouts())
    for (int r = 0; r < c.nparts; ++r)
      CHECK(l->ghost_local_size(r) == 0);
  CHECK(lonely_facets(mesh, d) > 0);

  Mesh m2;
  Part::Decomposition const full = decompose(m2, c);
  long ghost_zones = 0, ghost_sides = 0, ghost_corners = 0;
  for (int r = 0; r < c.nparts; ++r) {
    ghost_zones += full.zones.ghost_local_size(r);
    ghost_sides += full.sides.ghost_local_size(r);
    ghost_corners += full.corners.ghost_local_size(r);
  }
  CHECK(ghost_zones > 0);
  CHECK(ghost_sides > 0);
  CHECK(ghost_corners > 0);
  CHECK(lonely_facets(m2, full) == 0);
}

TEST_CASE("decompose: boundary zones follow their real zone", "[partition]") {
  Case c{mixed_spec(), "slab-y", 3};
  Mesh mesh;
  Part::Decomposition const d = decompose(mesh, c);
  auto const &s2z = mesh.ds->caccess_intv("m:s>z");
  auto const &s2s2 = mesh.ds->caccess_intv("m:s>s2");
  auto const &smask = mesh.sides.mask;
  auto const &zmask = mesh.zones.mask;

  int checked = 0;
  for (int s = 0; s < mesh.sides.size(); ++s) {
    if (smask[s] >= 0)
      continue;
    int const z = s2z[static_cast<std::size_t>(s)];
    if (zmask[z] >= 0)
      continue; // phantom-sides-only: the side is in the real zone already
    int const real = s2z[static_cast<std::size_t>(
        s2s2[static_cast<std::size_t>(s)])];
    CHECK(d.zones.owner[static_cast<std::size_t>(z)] ==
        d.zones.owner[static_cast<std::size_t>(real)]);
    checked += 1;
  }
  CHECK(checked > 0);
}

TEST_CASE("decompose: the owner rules all pick a rank that holds the entity",
    "[partition]") {
  Case c{mixed_spec(), "slab-x", 4};
  for (auto rule : {Part::Owner_Rule::LOWEST_RANK,
           Part::Owner_Rule::LEAST_LOADED, Part::Owner_Rule::HASH}) {
    INFO("rule " << static_cast<int>(rule));
    Mesh mesh;
    Part::Decompose_Options dopt;
    dopt.owner_rule = rule;
    Part::Decomposition const d = decompose(mesh, c, dopt);
    for (Part::Entity_Layout const *l : d.layouts())
      for (int g = 0; g < l->nglobal(); ++g) {
        int const owner = l->owner[static_cast<std::size_t>(g)];
        REQUIRE(owner >= 0);
        CHECK(l->local_index(owner, g) >= 0);
        CHECK(l->local_index(owner, g) < l->local_size(owner));
      }
  }

  /* Shared points are what the rule is about, and the default piles them onto
     the low ranks -- which is the imbalance the alternatives exist to fix. */
  Mesh mesh;
  Part::Decompose_Options low, balanced;
  balanced.owner_rule = Part::Owner_Rule::LEAST_LOADED;
  Part::Decomposition const a = decompose(mesh, c, low);
  Mesh mesh2;
  Part::Decomposition const b = decompose(mesh2, c, balanced);

  auto sourced = [](Part::Entity_Layout const &l, int nparts) {
    std::vector<long> n(static_cast<std::size_t>(nparts), 0);
    for (int owner : l.owner)
      if (owner >= 0)
        n[static_cast<std::size_t>(owner)] += 1;
    return n;
  };
  auto const na = sourced(a.points, c.nparts);
  auto const nb = sourced(b.points, c.nparts);
  auto spread = [](std::vector<long> const &v) {
    return *std::ranges::max_element(v) - *std::ranges::min_element(v);
  };
  CHECK(spread(nb) <= spread(na));
}

TEST_CASE("decompose: rejects a zone map it cannot use", "[partition]") {
  Mesh mesh;
  brick(mixed_spec(), mesh);
  Part::Options opt;
  opt.nparts = 2;
  Part::Zone_Map z2r;
  Part::make("linear")->assign(mesh, opt, z2r);

  CHECK_THROWS_AS(Part::build_decomposition(mesh, z2r, 0), std::invalid_argument);
  CHECK_THROWS_AS(
      Part::build_decomposition(mesh, Part::Zone_Map{0, 1}, 2),
      std::invalid_argument);

  Part::Zone_Map bad = z2r;
  for (int z = 0; z < mesh.zones.size(); ++z)
    if (mesh.zones.mask[z] >= 1) {
      bad[static_cast<std::size_t>(z)] = 2;
      break;
    }
  CHECK_THROWS_AS(
      Part::build_decomposition(mesh, bad, 2), std::invalid_argument);
}

TEST_CASE("partition: scarce zones fill the low parts", "[partition]") {
  /* Asking for more parts than there are zones is a legitimate thing to do by
     accident, and the answer has to be "the first n parts get one zone each".
     The weight rule alone would put a lone zone at its weight midpoint --
     part k/2 -- leaving part 0 empty, and an empty part 0 is the one that
     tends to be looked at first. */
  Brick::Spec spec;
  spec.nx = spec.ny = spec.nz = 1;
  Mesh mesh;
  brick(spec, mesh);
  int const zones = static_cast<int>(Part::real_zones(mesh).size());
  REQUIRE(zones == 1);

  Part::Options opt;
  Part::Zone_Map z2r;
  for (std::string const &name : Part::methods())
    for (int nparts : {1, 2, 3, 8}) {
      INFO(name << " nparts=" << nparts);
      opt.nparts = nparts;
      Part::make(name)->assign(mesh, opt, z2r);
      for (int z = 0; z < mesh.zones.size(); ++z)
        if (mesh.zones.mask[z] >= 1)
          CHECK(z2r[static_cast<std::size_t>(z)] == 0);
    }
}
/* ------------------------------- extraction ------------------------------- */

TEST_CASE("extract: one part reproduces the input mesh", "[partition]") {
  /* With the `(is_ghost, global index)` ordering rule the one-rank map is the
     identity, so this compares every stored array of a real mesh against
     itself through the whole renumbering path -- the cheapest test available
     that a map was not quietly permuted, and the reason `Mesh::operator==`
     had to learn about faces first. */
  for (Case const &c : cases()) {
    if (c.nparts != 1)
      continue;
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);

    Mesh part;
    Ume::attach_serial_transport(part);
    Part::extract(mesh, d, 0, part);

    CHECK(part.numpe == 1);
    CHECK(part.mype == 0);
    CHECK(part == mesh);
  }
}

TEST_CASE("extract: every partition passes the mesh audit", "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    for (int r = 0; r < c.nparts; ++r) {
      INFO("rank " << r);
      Mesh part;
      Ume::attach_serial_transport(part);
      Part::extract(mesh, d, r, part);
      std::ostringstream report;
      bool const ok = SOA_Idx::audit(part, report, SOA_Idx::Audit_Options{}).ok();
      INFO(report.str());
      CHECK(ok);
    }
  }
}

TEST_CASE("extract: degenerate part counts still produce valid meshes",
    "[partition]") {
  /* One cell, and more parts than there are zones to fill them.  An empty
     part is a real outcome of asking for more ranks than the mesh can feed,
     and it has to come out as a well-formed empty mesh rather than a crash. */
  Brick::Spec spec;
  spec.nx = spec.ny = spec.nz = 1;

  for (int nparts : {1, 2, 8}) {
    INFO("nparts " << nparts);
    Mesh mesh;
    brick(spec, mesh);
    Part::Options opt;
    opt.nparts = nparts;
    Part::Zone_Map z2r;
    Part::make("slab-x")->assign(mesh, opt, z2r);
    Part::Decomposition const d =
        Part::build_decomposition(mesh, std::move(z2r), nparts);

    int empty = 0;
    for (int r = 0; r < nparts; ++r) {
      INFO("rank " << r);
      Mesh part;
      Ume::attach_serial_transport(part);
      Part::extract(mesh, d, r, part);
      if (part.zones.size() == 0)
        empty += 1;
      std::ostringstream report;
      bool const ok =
          SOA_Idx::audit(part, report, SOA_Idx::Audit_Options{}).ok();
      INFO(report.str());
      CHECK(ok);
    }
    CHECK(empty == nparts - 1);
  }
}

/* --------------------------- the global-id side-car ----------------------- */

TEST_CASE("global ids: the side-car round trips", "[partition]") {
  /* The side-car is the only record of what a partition's entities used to be
     called -- a `.ume` file has nowhere to put it -- so a run comparing itself
     against a serial one is trusting this file completely.  What is checked
     here is that it survives the trip through text and still agrees with the
     layouts it came from, kind by kind and rank by rank. */
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    auto const layouts = d.layouts();

    for (int r = 0; r < c.nparts; ++r) {
      INFO("rank " << r);
      std::ostringstream out;
      Part::write_global_ids(d, r, out);
      std::istringstream in(out.str());
      Part::Global_Ids const g = Part::read_global_ids(in);

      CHECK(g.mype == r);
      CHECK(g.numpe == c.nparts);
      for (int k = 0; k < Part::NKINDS; ++k) {
        INFO(Part::kind_name(k));
        Part::Entity_Layout const &l = *layouts[static_cast<std::size_t>(k)];
        Part::Kind_Ids const &ids = g[k];
        REQUIRE(ids.local == l.local_size(r));
        REQUIRE(ids.ghost() == l.ghost_local_size(r));
        for (int i = 0; i < static_cast<int>(ids.gid.size()); ++i)
          REQUIRE(ids.gid[static_cast<std::size_t>(i)] == l.global_index(r, i));
      }
    }
  }
}

TEST_CASE("global ids: one part is the identity", "[partition]") {
  /* The same property `extract: one part reproduces the input mesh` rests on,
     seen from the side-car: at one rank nothing is renumbered, so a reference
     run reads its own zone indices back unchanged and the comparison it
     anchors is against the serial mesh rather than against a relabeling. */
  Mesh mesh;
  Case const c{mixed_spec(), "linear", 1};
  Part::Decomposition const d = decompose(mesh, c);

  std::ostringstream out;
  Part::write_global_ids(d, 0, out);
  std::istringstream in(out.str());
  Part::Global_Ids const g = Part::read_global_ids(in);

  Part::Kind_Ids const &zones = g[Part::ZONE];
  REQUIRE(zones.local == mesh.zones.local_size());
  REQUIRE(zones.ghost() == 0);
  for (int z = 0; z < zones.local; ++z)
    REQUIRE(zones.gid[static_cast<std::size_t>(z)] == z);
}

TEST_CASE("global ids: a truncated side-car is rejected", "[partition]") {
  /* Reading the wrong file, or a partly written one, has to be an error
     rather than a short list of plausible indices: the comparison downstream
     would otherwise line the wrong zones up and disagree about physics. */
  Mesh mesh;
  Case const c{mixed_spec(), "slab-x", 3};
  Part::Decomposition const d = decompose(mesh, c);

  std::ostringstream out;
  Part::write_global_ids(d, 1, out);
  std::string const full = out.str();

  std::istringstream truncated(full.substr(0, full.size() / 2));
  CHECK_THROWS_AS(Part::read_global_ids(truncated), std::runtime_error);

  std::istringstream empty("");
  CHECK_THROWS_AS(Part::read_global_ids(empty), std::runtime_error);

  std::istringstream wrong("ume_gid 2\n");
  CHECK_THROWS_AS(Part::read_global_ids(wrong), std::runtime_error);
}

/* ---------------------------- cross-rank checks --------------------------- */

namespace {

//! Extract every rank of `d`.  A Mesh cannot be moved, so they are held by
//! pointer -- which is also the shape `check_partitions` wants them in.
std::vector<std::unique_ptr<Mesh>> extract_all(
    Mesh const &mesh, Part::Decomposition const &d) {
  std::vector<std::unique_ptr<Mesh>> parts;
  for (int r = 0; r < d.nparts; ++r) {
    parts.push_back(std::make_unique<Mesh>());
    Ume::attach_serial_transport(*parts.back());
    Part::extract(mesh, d, r, *parts.back());
  }
  return parts;
}

std::vector<Mesh const *> view(
    std::vector<std::unique_ptr<Mesh>> const &parts) {
  std::vector<Mesh const *> v;
  v.reserve(parts.size());
  for (auto const &p : parts)
    v.push_back(p.get());
  return v;
}

//! A three-way split of the mixed brick, for the tests that break one on purpose
struct Split {
  Mesh mesh;
  std::vector<std::unique_ptr<Mesh>> parts;
  Split() {
    Case const c{mixed_spec(), "slab-x", 3};
    Part::Decomposition const d = decompose(mesh, c);
    parts = extract_all(mesh, d);
  }
  //! Run the checks, leaving what they said in `last_report` for the caller
  Part::Check_Result check() const {
    std::ostringstream report;
    Part::Check_Result const res =
        Part::check_partitions(mesh, view(parts), report);
    last_report = report.str();
    return res;
  }
  mutable std::string last_report;
};

} // namespace

TEST_CASE("check: the ranks of a decomposition agree with each other",
    "[partition]") {
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    auto const parts = extract_all(mesh, d);

    std::ostringstream report;
    Part::Check_Result const res =
        Part::check_partitions(mesh, view(parts), report);
    INFO(report.str());
    CHECK(res.ok());
    CHECK(res.problems == 0);
    CHECK(res.checks > 0);
  }
}

TEST_CASE("check: the ranks alone are checkable without the serial mesh",
    "[partition]") {
  /* This is the form `ume_serial` gets: N files and nothing to compare them
     against but each other, so the reassembly checks do not run. */
  Case const c{mixed_spec(), "slab-y", 4};
  Mesh mesh;
  Part::Decomposition const d = decompose(mesh, c);
  auto const parts = extract_all(mesh, d);

  std::ostringstream with, without;
  Part::Check_Result const a = Part::check_partitions(mesh, view(parts), with);
  Part::Check_Result const b = Part::check_partitions(view(parts), without);
  INFO(with.str() << without.str());
  CHECK(a.ok());
  CHECK(b.ok());
  CHECK(b.checks < a.checks);
}

TEST_CASE("check: a permuted copy list is caught", "[partition]") {
  /* The failure this file exists for.  Swapping two entries of one copy list
     leaves every length, comm type and mask right, and leaves each rank
     passing `mesh_audit` on its own -- but the two ranks now disagree about
     which entity position 0 of the buffer holds, so every gather would put
     the values on the wrong points.  Nothing local can see it. */
  Split s;
  REQUIRE(s.check().ok());

  Ume::Comm::Neighbor *victim{nullptr};
  for (auto const &p : s.parts)
    for (auto &nbr : p->points.myCpys)
      if (victim == nullptr && nbr.elements.size() >= 2)
        victim = &nbr;
  REQUIRE(victim != nullptr);
  std::swap(victim->elements[0], victim->elements[1]);

  Part::Check_Result const res = s.check();
  INFO(s.last_report);
  CHECK(!res.ok());
  CHECK(res.failed_checks == 1);
  CHECK(res.problems >= 2);
}

TEST_CASE("check: a copy list that outgrows its source list is caught",
    "[partition]") {
  Split s;
  REQUIRE(s.check().ok());

  Ume::Comm::Neighbor *victim{nullptr};
  for (auto const &p : s.parts)
    for (auto &nbr : p->points.mySrcs)
      if (victim == nullptr && !nbr.elements.empty())
        victim = &nbr;
  REQUIRE(victim != nullptr);
  victim->elements.pop_back();

  Part::Check_Result const res = s.check();
  INFO(s.last_report);
  CHECK(!res.ok());
  CHECK(res.failed_checks == 1);
}

TEST_CASE("check: a ghost pointed at the wrong source is caught",
    "[partition]") {
  Split s;
  REQUIRE(s.check().ok());

  bool moved{false};
  for (auto const &p : s.parts) {
    if (p->zones.ghost_local_size() == 0)
      continue;
    int const pe = p->zones.src_pe[0];
    int const n = s.parts[static_cast<std::size_t>(pe)]->zones.local_size();
    REQUIRE(n > 1);
    p->zones.src_idx[0] = (p->zones.src_idx[0] + 1) % n;
    moved = true;
    break;
  }
  REQUIRE(moved);

  Part::Check_Result const res = s.check();
  INFO(s.last_report);
  CHECK(!res.ok());
  CHECK(res.failed_checks == 1);
}

TEST_CASE("check: a zone no rank claims is caught", "[partition]") {
  /* Counting is not enough on its own, so both reassembly checks fire: one
     rank fewer zone in the total, and one zone of the serial mesh that
     nothing covers. */
  Split s;
  REQUIRE(s.check().ok());

  bool disowned{false};
  auto &zones = s.parts.front()->zones;
  for (int i = 0; i < zones.local_size(); ++i)
    if (zones.comm_type[static_cast<std::size_t>(i)] ==
        SOA_Idx::Entity::INTERNAL) {
      zones.comm_type[static_cast<std::size_t>(i)] = SOA_Idx::Entity::COPY;
      disowned = true;
      break;
    }
  REQUIRE(disowned);

  Part::Check_Result const res = s.check();
  INFO(s.last_report);
  CHECK(!res.ok());
  CHECK(res.failed_checks == 2);
}

TEST_CASE("check: two sides on one edge are told apart", "[partition]") {
  /* A side's two endpoints do not name it: the two sides of a zone that meet
     along an edge, one on each of the two faces sharing that edge, have the
     same pair.  Swapping such a pair in a copy list is the permutation that a
     coordinate comparison is most likely to wave through, so it is the one
     worth pinning down. */
  Split s;
  REQUIRE(s.check().ok());

  Ume::Comm::Neighbor *victim{nullptr};
  std::size_t a{0}, b{0};
  for (auto const &p : s.parts) {
    auto const &s2p1 = p->ds->caccess_intv("m:s>p1");
    auto const &s2p2 = p->ds->caccess_intv("m:s>p2");
    auto ends = [&](int i) {
      auto const u = static_cast<std::size_t>(i);
      return std::minmax(s2p1[u], s2p2[u]);
    };
    for (auto &nbr : p->sides.myCpys)
      for (std::size_t i = 0; i < nbr.elements.size() && victim == nullptr; ++i)
        for (std::size_t j = i + 1; j < nbr.elements.size(); ++j)
          if (ends(nbr.elements[i]) == ends(nbr.elements[j])) {
            victim = &nbr;
            a = i;
            b = j;
            break;
          }
    if (victim != nullptr)
      break;
  }
  REQUIRE(victim != nullptr);
  std::swap(victim->elements[a], victim->elements[b]);

  Part::Check_Result const res = s.check();
  INFO(s.last_report);
  CHECK(!res.ok());
  CHECK(res.failed_checks == 1);
}

/* -------------------------------- metrics --------------------------------- */

namespace {

SOA_Idx::Entity const &entity_of(Mesh const &m, int const kind) {
  switch (kind) {
  case Part::POINT:
    return m.points;
  case Part::EDGE:
    return m.edges;
  case Part::FACE:
    return m.faces;
  case Part::ZONE:
    return m.zones;
  case Part::CORNER:
    return m.corners;
  case Part::SIDE:
    return m.sides;
  default:
    return m.iotas;
  }
}

long elements_in(Ume::Comm::Neighbors const &n) {
  long total = 0;
  for (auto const &nbr : n)
    total += static_cast<long>(nbr.elements.size());
  return total;
}

//! z % nparts over the real zones: the worst zone map there is
Part::Zone_Map scattered(Mesh const &mesh, int nparts) {
  Part::Zone_Map z2r(static_cast<std::size_t>(mesh.zones.size()), -1);
  int next = 0;
  for (int z = 0; z < mesh.zones.size(); ++z)
    if (mesh.zones.mask[z] >= 1)
      z2r[static_cast<std::size_t>(z)] = next++ % nparts;
  return z2r;
}

} // namespace

TEST_CASE("metrics: one rank costs nothing", "[partition]") {
  for (Case const &c : cases()) {
    if (c.nparts != 1)
      continue;
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    Part::Metrics const m = Part::measure(mesh, d);

    REQUIRE(m.nparts == 1);
    for (int k = 0; k < Part::NKINDS; ++k) {
      INFO(Part::kind_name(k));
      auto const &km = m.kinds[static_cast<std::size_t>(k)];
      CHECK(km.local.front() == km.serial);
      CHECK(km.total_ghost == 0);
      CHECK(km.duplicated == 0);
      CHECK(km.src_elems.front() == 0);
      CHECK(km.cpy_elems.front() == 0);
      CHECK(km.max_neighbors == 0);
      CHECK(km.imbalance == (km.total_load > 0 ? 1.0 : 0.0));
    }
    CHECK(m.cut_faces == 0);
    CHECK(m.cut_face_points == 0);
    CHECK(m.interior_faces > 0);
    CHECK(m.max_neighbors == 0);
    CHECK(m.total_comm_elems() == 0);
    CHECK(m.max_surface_to_volume == 0.0);
  }
}

TEST_CASE("metrics: the numbers are what extract produces", "[partition]") {
  /* The whole point of measuring the `Decomposition` rather than the meshes
     is that a method sweep never has to build them.  That is only worth
     anything if the two agree exactly, so every counter is pinned against the
     mesh it predicts -- including the send matrix, which is the one a reader
     is most likely to take on trust. */
  for (Case const &c : cases()) {
    INFO(c.describe());
    Mesh mesh;
    Part::Decomposition const d = decompose(mesh, c);
    Part::Metrics const m = Part::measure(mesh, d);
    auto const parts = extract_all(mesh, d);

    REQUIRE(m.nparts == d.nparts);
    std::vector<long> sends(
        static_cast<std::size_t>(d.nparts) * static_cast<std::size_t>(d.nparts),
        0);

    for (int pe = 0; pe < d.nparts; ++pe) {
      auto const u = static_cast<std::size_t>(pe);
      INFO("rank " << pe);
      std::set<int> all_nbrs;
      for (int k = 0; k < Part::NKINDS; ++k) {
        INFO(Part::kind_name(k));
        auto const &km = m.kinds[static_cast<std::size_t>(k)];
        SOA_Idx::Entity const &e = entity_of(*parts[u], k);

        CHECK(km.local[u] == e.local_size());
        CHECK(km.ghost[u] == e.ghost_local_size());
        CHECK(km.src_elems[u] == elements_in(e.mySrcs));
        CHECK(km.cpy_elems[u] == elements_in(e.myCpys));

        std::set<int> nbrs;
        for (auto const *n : {&e.mySrcs, &e.myCpys})
          for (auto const &nbr : *n) {
            nbrs.insert(nbr.pe);
            all_nbrs.insert(nbr.pe);
          }
        CHECK(km.neighbors[u] == static_cast<int>(nbrs.size()));

        for (auto const &nbr : e.mySrcs)
          sends[u * static_cast<std::size_t>(d.nparts) +
              static_cast<std::size_t>(nbr.pe)] +=
              static_cast<long>(nbr.elements.size());
      }
      CHECK(m.neighbors[u] == static_cast<int>(all_nbrs.size()));
    }
    CHECK(m.comm_matrix == sends);
  }
}

TEST_CASE(
    "metrics: a scattered zone map costs more than a slab", "[partition]") {
  /* Reporting these numbers is only useful if they can tell a good
     decomposition from a bad one, so here is the worst one available.  Every
     face of the mesh is cut, so every rank talks to every other -- which is
     the latency half of the story -- and every point it holds is shared. */
  constexpr int NPARTS = 4;
  Mesh mesh;
  brick(mixed_spec(), mesh);

  Part::Options opt;
  opt.nparts = NPARTS;
  Part::Zone_Map slab;
  Part::make("slab-x")->assign(mesh, opt, slab);

  Part::Metrics const good = Part::measure(
      mesh, Part::build_decomposition(mesh, std::move(slab), NPARTS));
  Part::Metrics const bad = Part::measure(
      mesh, Part::build_decomposition(mesh, scattered(mesh, NPARTS), NPARTS));

  CHECK(bad.kinds[Part::POINT].duplicated > good.kinds[Part::POINT].duplicated);
  CHECK(bad.cut_faces > good.cut_faces);
  CHECK(bad.total_comm_elems() > good.total_comm_elems());
  CHECK(bad.max_neighbors == NPARTS - 1);
  CHECK(bad.max_neighbors >= good.max_neighbors);
  CHECK(bad.mean_surface_to_volume >= good.mean_surface_to_volume);

  /* Load balance is the one thing round robin is good at, and the metric has
     to say so rather than condemning it across the board. */
  CHECK(bad.kinds[Part::ZONE].imbalance <= good.kinds[Part::ZONE].imbalance);
}

TEST_CASE("metrics: the CSV row matches its header", "[partition]") {
  /* A sweep concatenates rows from separately configured builds under one
     header, so a row that has drifted from the header is silently wrong data
     rather than a broken file. */
  Mesh mesh;
  Part::Decomposition const d =
      decompose(mesh, Case{mixed_spec(), "slab-y", 3});
  Part::Metrics const m = Part::measure(mesh, d);

  auto fields = [](std::string const &line) {
    return 1 + std::count(line.begin(), line.end(), ',');
  };

  std::ostringstream head, row;
  Part::report_csv_header(head);
  Part::report_csv(m, "slab-y", row);
  INFO(head.str() << row.str());

  std::string const h = head.str(), r = row.str();
  REQUIRE(h.back() == '\n');
  REQUIRE(r.back() == '\n');
  CHECK(std::count(h.begin(), h.end(), '\n') == 1);
  CHECK(std::count(r.begin(), r.end(), '\n') == 1);
  CHECK(fields(h) == fields(r));
  CHECK(r.starts_with("slab-y,3,"));

  /* And the matrix is nparts rows of nparts columns, whatever is in it. */
  std::ostringstream matrix;
  Part::write_comm_matrix(m, matrix);
  std::istringstream in(matrix.str());
  int rows = 0;
  for (std::string line; std::getline(in, line); ++rows)
    CHECK(fields(line) == d.nparts);
  CHECK(rows == d.nparts);
}
