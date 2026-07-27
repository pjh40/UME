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
#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <iterator>
#include <numeric>
#include <span>
#include <sstream>

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
  char const *method;
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
      for (char const *method : {"linear", "slab-x", "slab-z"})
        for (int n : {1, 2, 3, 5})
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
