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
#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <numeric>
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
