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
  \file test_mesh.cc

  Tests for the invariants of SOA_Idx::Mesh itself, as distinct from any
  particular mesh it might hold.
*/

#include "Ume/Mesh_Base.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include "Ume/face_area.hh"
#include "Ume/gradient.hh"
#include "Ume/mem_exec_spaces.hh"
#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <ranges>
#include <vector>

using Ume::SOA_Idx::Mesh;

/* A partition can legitimately come out empty -- ask for more ranks than the
   mesh has zones and some of them get nothing -- so the derived variables have
   to survive a mesh with no entities in it.  They take the address of element
   zero of the arrays they read, which is a hard error on an empty one. */
TEST_CASE("an empty mesh can still derive its variables", "[mesh]") {
  Mesh mesh;
  Ume::Comm::Dummy_Transport transport;
  mesh.comm = &transport;
  mesh.ivtag = UME_VERSION_2;
  mesh.version_header = true;
  mesh.mype = 0;
  mesh.numpe = 1;
  mesh.geo = Mesh::CARTESIAN;
  mesh.dump_iotas = false;
  Ume::SOA_Idx::Entity *const entities[] = {&mesh.points, &mesh.edges,
      &mesh.faces, &mesh.zones, &mesh.corners, &mesh.sides, &mesh.iotas};
  for (Ume::SOA_Idx::Entity *e : entities)
    e->resize(0, 0, 0);

  for (char const *name : {"zcoord", "fcoord", "ecoord", "side_surf",
           "side_surz", "corner_csurf", "point_norm"})
    CHECK(mesh.ds->caccess_vec3v(name).empty());
  for (char const *name : {"side_vol", "corner_vol"})
    CHECK(mesh.ds->caccess_dblv(name).empty());
  /* The ragged derivations index their accumulators with `at()`, so all they
     have to do here is not throw -- there are no rows to look at. */
  for (char const *name :
      {"m:z>c", "m:z>p", "m:z>pz", "m:c>s", "m:p>zs", "m:p>rc"})
    CHECK_NOTHROW(mesh.ds->caccess_intrr(name));

  /* The kernels that consume those variables build the same kind of View over
     their own arrays, so they have to survive the empty mesh too.  The gradient
     outputs start non-empty so that the size checks see the kernels resize
     them. */
  Ume::DS_Types::DBLV_T const zone_field;
  Ume::DS_Types::VEC3V_T zone_gradient{Ume::DS_Types::VEC3_T(1.0)};
  Ume::DS_Types::VEC3V_T point_gradient{Ume::DS_Types::VEC3_T(1.0)};
  Ume::gradzatz(mesh, zone_field, zone_gradient, point_gradient);
  CHECK(zone_gradient.empty());
  CHECK(point_gradient.empty());

  zone_gradient.assign(1, Ume::DS_Types::VEC3_T(1.0));
  point_gradient.assign(1, Ume::DS_Types::VEC3_T(1.0));
  Ume::gradzatz_invert(mesh, zone_field, zone_gradient, point_gradient);
  CHECK(zone_gradient.empty());
  CHECK(point_gradient.empty());

  /* calc_face_area fills a caller-sized output, so there is nothing to check
     beyond surviving the call. */
  Ume::DS_Types::DBLV_T face_area;
  Ume::calc_face_area(mesh, face_area);
}

namespace {

/* The number of faces whose area is not `reference` to rounding.  Concurrent
   atomic adds may sum a face's sides in any order, so the comparison is not
   exact. */
std::size_t count_wrong_areas(Ume::DS_Types::DBLV_T const &face_area,
    Ume::DS_Types::DBLV_T const &reference) {
  std::size_t wrong = 0;
  for (std::size_t f = 0; f < reference.size(); ++f)
    if (!(std::abs(face_area[f] - reference[f]) <= 1e-12 * reference[f]))
      ++wrong;
  return wrong;
}

} // namespace

/* calc_face_area walks the sides, and each face is covered twice, once by the
   sides of each zone it separates; a side and its partner `m:s>s2` are the
   same triangle seen from either zone, so the face area takes one of the two.
   On a threaded backend the two members of a pair can run concurrently, and
   the choice must not depend on which got there first.  The mesh here is
   disjoint triangular faces, three sides a zone, with the side numbering
   shuffled so that the members of a pair land at arbitrary distances on
   arbitrary threads, and the call is repeated because a race need not show
   on any one run.  Every fourth face is a non-source copy of a face owned
   elsewhere, which does not count its sides here.  Every eighth is on the
   problem boundary, its partner sides not real (`mask < 1`), and every eighth
   has sides that are their own partners, which the side tag the kernel used
   to keep counted once. */
TEST_CASE("face areas take one side of each pair, whichever thread runs it",
    "[mesh][face_area]") {
  constexpr int nface = 200'000;
  constexpr int sides_per_zone = 3;
  constexpr int npoint = sides_per_zone * nface;
  constexpr int nside = 2 * sides_per_zone * nface;
  constexpr int repetitions = 20;

  Mesh mesh;
  Ume::Comm::Dummy_Transport transport;
  mesh.comm = &transport;
  mesh.ivtag = UME_VERSION_2;
  mesh.version_header = true;
  mesh.mype = 0;
  mesh.numpe = 1;
  mesh.geo = Mesh::CARTESIAN;
  mesh.dump_iotas = false;
  Ume::SOA_Idx::Entity *const entities[] = {&mesh.points, &mesh.edges,
      &mesh.faces, &mesh.zones, &mesh.corners, &mesh.sides, &mesh.iotas};
  for (Ume::SOA_Idx::Entity *e : entities)
    e->resize(0, 0, 0);
  mesh.points.resize(npoint, npoint, 0);
  mesh.faces.resize(nface, nface, 0);
  mesh.sides.resize(nside, nside, 0);

  auto const is_copy = [](int const f) { return f % 4 == 3; };
  auto const is_boundary = [](int const f) { return f % 8 == 1; };
  auto const is_self_paired = [](int const f) { return f % 8 == 5; };

  /* Face f is the triangle on points 3f, 3f+1, 3f+2, of area (1 + f % 3)/2. */
  auto &pcoord = mesh.ds->access_vec3v("pcoord");
  for (int f = 0; f < nface; ++f) {
    double const z = f;
    double const base = 1.0 + f % 3;
    pcoord[3 * f + 0] = Ume::Vec3(std::array{0.0, 0.0, z});
    pcoord[3 * f + 1] = Ume::Vec3(std::array{base, 0.0, z});
    pcoord[3 * f + 2] = Ume::Vec3(std::array{0.0, 1.0, z});
  }

  std::fill(mesh.faces.mask.begin(), mesh.faces.mask.end(), short{1});
  for (int f = 0; f < nface; ++f)
    mesh.faces.comm_type[f] = is_copy(f) ? Ume::SOA_Idx::Entity::COPY
                                         : Ume::SOA_Idx::Entity::INTERNAL;

  /* Slot k < 3 of face f is the side on edge {k, k+1} in the first zone; slot
     k + 3 is the same edge, reversed, in the second. */
  std::vector<int> side_of_slot(nside);
  std::iota(side_of_slot.begin(), side_of_slot.end(), 0);
  std::shuffle(side_of_slot.begin(), side_of_slot.end(), std::mt19937{8});
  auto &s_to_f = mesh.ds->access_intv("m:s>f");
  auto &s_to_p1 = mesh.ds->access_intv("m:s>p1");
  auto &s_to_p2 = mesh.ds->access_intv("m:s>p2");
  auto &s_to_s2 = mesh.ds->access_intv("m:s>s2");
  for (int f = 0; f < nface; ++f) {
    for (int k = 0; k < sides_per_zone; ++k) {
      int const s = side_of_slot[2 * sides_per_zone * f + k];
      int const s2 = side_of_slot[2 * sides_per_zone * f + sides_per_zone + k];
      int const a = sides_per_zone * f + k;
      int const b = sides_per_zone * f + (k + 1) % sides_per_zone;
      s_to_f[s] = s_to_f[s2] = f;
      s_to_p1[s] = s_to_p2[s2] = a;
      s_to_p2[s] = s_to_p1[s2] = b;
      s_to_s2[s] = is_self_paired(f) ? s : s2;
      s_to_s2[s2] = is_self_paired(f) ? s2 : s;
      mesh.sides.mask[s] = 1;
      mesh.sides.mask[s2] = is_self_paired(f) ? short{0}
          : is_boundary(f)                    ? short{-1}
                                              : short{1};
    }
  }

  /* The single-threaded reference: the first zone's sides of each face that
     this rank owns.  It reads the same side_surz the kernel does, so it holds
     whatever the face centers came out as. */
  auto const &surz = mesh.ds->caccess_vec3v("side_surz");
  Ume::DS_Types::DBLV_T reference(nface, 0.0);
  int degenerate = 0;
  for (int f = 0; f < nface; ++f) {
    if (is_copy(f))
      continue;
    for (int k = 0; k < sides_per_zone; ++k)
      reference[f] +=
          Ume::vectormag(surz[side_of_slot[2 * sides_per_zone * f + k]]);
    if (!(reference[f] > 0.0))
      ++degenerate;
  }
  REQUIRE(degenerate == 0);

  /* Positive control: one face counted twice is one wrong face. */
  Ume::DS_Types::DBLV_T doubled = reference;
  doubled[0] *= 2.0;
  REQUIRE(count_wrong_areas(doubled, reference) == 1);

  Ume::DS_Types::DBLV_T face_area(nface);
  for (int rep = 0; rep < repetitions; ++rep) {
    std::fill(face_area.begin(), face_area.end(), -1.0);
    Ume::calc_face_area(mesh, face_area);
    INFO("repetition " << rep);
    CHECK(count_wrong_areas(face_area, reference) == 0);
  }
}

namespace {

/* The number of entries of `centroid` that differ from `reference`.  The
   coordinates are small integers, so every partial sum is exact and the
   centroids are equal bit for bit whatever order the points were added in. */
std::size_t count_wrong_centroids(Ume::DS_Types::VEC3V_T const &centroid,
    Ume::DS_Types::VEC3V_T const &reference) {
  std::size_t wrong = 0;
  for (std::size_t i = 0; i < reference.size(); ++i)
    if (!(centroid[i] == reference[i]))
      ++wrong;
  return wrong;
}

/* The centroid of the points `elem_to_pt` hangs on each of `nelem` elements,
   over the members whose `mask` is set, summed in one thread. */
Ume::DS_Types::VEC3V_T serial_centroids(int const nelem,
    std::vector<int> const &member_to_elem,
    std::vector<int> const &member_to_pt, std::vector<short> const &mask,
    Ume::DS_Types::VEC3V_T const &pcoord) {
  Ume::DS_Types::VEC3V_T centroid(nelem, Ume::Vec3(0.0));
  std::vector<int> count(nelem, 0);
  for (std::size_t m = 0; m < member_to_elem.size(); ++m) {
    if (mask[m] == 0)
      continue;
    centroid[member_to_elem[m]] += pcoord[member_to_pt[m]];
    ++count[member_to_elem[m]];
  }
  for (int e = 0; e < nelem; ++e)
    centroid[e] /= static_cast<double>(count[e]);
  return centroid;
}

} // namespace

/* zcoord and fcoord sum the points of each zone (face) over its corners
   (sides) in a parallel loop, and divide by how many there were.  With a
   threaded host backend, two corners of one zone counted at once must both be
   counted.  The mesh here has few zones and faces and many corners and sides
   on each, numbered round-robin so that every thread's block of the loop
   touches every zone, and every sixteenth corner and side is masked out.  The
   derived variables are cached once computed, so each repetition builds a
   fresh mesh. */
TEST_CASE("zone and face centroids count every point, whichever thread "
          "counts it",
    "[mesh][centroid]") {
  constexpr int nelem = 64;
  constexpr int members_per_elem = 8192;
  constexpr int nmember = nelem * members_per_elem;
  constexpr int repetitions = 10;
  INFO("host concurrency " << HostExecSpace().concurrency());

  auto const elem_of = [](int const m) { return m % nelem; };
  auto const is_masked = [](int const m) { return m / nelem % 16 == 15; };

  for (int rep = 0; rep < repetitions; ++rep) {
    Ume::Comm::Dummy_Transport transport;
    Mesh mesh;
    mesh.comm = &transport;
    mesh.ivtag = UME_VERSION_2;
    mesh.version_header = true;
    mesh.mype = 0;
    mesh.numpe = 1;
    mesh.geo = Mesh::CARTESIAN;
    mesh.dump_iotas = false;
    Ume::SOA_Idx::Entity *const entities[] = {&mesh.points, &mesh.edges,
        &mesh.faces, &mesh.zones, &mesh.corners, &mesh.sides, &mesh.iotas};
    for (Ume::SOA_Idx::Entity *e : entities)
      e->resize(0, 0, 0);
    mesh.points.resize(nmember, nmember, 0);
    mesh.zones.resize(nelem, nelem, 0);
    mesh.faces.resize(nelem, nelem, 0);
    mesh.corners.resize(nmember, nmember, 0);
    mesh.sides.resize(nmember, nmember, 0);

    /* Point p hangs on corner p and side p. */
    auto &pcoord = mesh.ds->access_vec3v("pcoord");
    for (int p = 0; p < nmember; ++p)
      pcoord[p] = Ume::Vec3(std::array{static_cast<double>(p % 5),
          static_cast<double>(p % 7), static_cast<double>(p % 11)});
    std::fill(mesh.zones.mask.begin(), mesh.zones.mask.end(), short{1});
    std::fill(mesh.faces.mask.begin(), mesh.faces.mask.end(), short{1});

    auto &c_to_z = mesh.ds->access_intv("m:c>z");
    auto &c_to_p = mesh.ds->access_intv("m:c>p");
    auto &s_to_f = mesh.ds->access_intv("m:s>f");
    auto &s_to_p1 = mesh.ds->access_intv("m:s>p1");
    for (int m = 0; m < nmember; ++m) {
      c_to_z[m] = s_to_f[m] = elem_of(m);
      c_to_p[m] = s_to_p1[m] = m;
      mesh.corners.mask[m] = mesh.sides.mask[m] =
          is_masked(m) ? short{0} : short{1};
    }

    Ume::DS_Types::VEC3V_T const zone_reference =
        serial_centroids(nelem, c_to_z, c_to_p, mesh.corners.mask, pcoord);
    Ume::DS_Types::VEC3V_T const face_reference =
        serial_centroids(nelem, s_to_f, s_to_p1, mesh.sides.mask, pcoord);

    if (rep == 0) {
      /* Positive control: point 0 is the origin, so dropping corner 0 leaves
         zone 0's sum as it was and its count one short, which is what a
         lost increment does.  That is one wrong centroid. */
      std::vector<short> one_lost = mesh.corners.mask;
      one_lost[0] = 0;
      REQUIRE(count_wrong_centroids(
                  serial_centroids(nelem, c_to_z, c_to_p, one_lost, pcoord),
                  zone_reference) == 1);
    }

    INFO("repetition " << rep);
    auto const &zcoord = mesh.ds->caccess_vec3v("zcoord");
    REQUIRE(zcoord.size() == zone_reference.size());
    CHECK(count_wrong_centroids(zcoord, zone_reference) == 0);
    auto const &fcoord = mesh.ds->caccess_vec3v("fcoord");
    REQUIRE(fcoord.size() == face_reference.size());
    CHECK(count_wrong_centroids(fcoord, face_reference) == 0);
  }
}

namespace {

/* The indices `range` yields, in order. */
std::vector<int> indices_of(std::ranges::common_range auto const &range) {
  return {std::ranges::begin(range), std::ranges::end(range)};
}

/* A mesh with no entities in it, for the index range tests to resize one
   entity of. */
struct Bare_Mesh {
  Ume::Comm::Dummy_Transport transport;
  Mesh mesh;
  Bare_Mesh() {
    mesh.comm = &transport;
    mesh.ivtag = UME_VERSION_2;
    mesh.version_header = true;
    mesh.mype = 0;
    mesh.numpe = 1;
    mesh.geo = Mesh::CARTESIAN;
    mesh.dump_iotas = false;
    Ume::SOA_Idx::Entity *const entities[] = {&mesh.points, &mesh.edges,
        &mesh.faces, &mesh.zones, &mesh.corners, &mesh.sides, &mesh.iotas};
    for (Ume::SOA_Idx::Entity *e : entities)
      e->resize(0, 0, 0);
  }
};

} // namespace

/* The index ranges are half-open: every element of an Entity is in exactly
   one of local_indices() and ghost_indices(), and all_indices() is their
   concatenation.  The sizes are checked before the contents, so that a range
   whose bound lies below its start fails rather than iterates. */
TEST_CASE("entity index ranges cover every element", "[mesh][entity]") {
  Bare_Mesh bare;
  auto &points = bare.mesh.points;
  points.resize(5, 8, 3);

  std::vector<int> const every_point{0, 1, 2, 3, 4, 5, 6, 7};
  REQUIRE(std::ranges::size(points.all_indices()) == 8);
  CHECK(indices_of(points.all_indices()) == every_point);
  REQUIRE(std::ranges::size(points.local_indices()) == 5);
  CHECK(indices_of(points.local_indices()) == std::vector<int>{0, 1, 2, 3, 4});
  REQUIRE(std::ranges::size(points.ghost_indices()) == 3);
  CHECK(indices_of(points.ghost_indices()) == std::vector<int>{5, 6, 7});
  CHECK(points.ghost_size() == 3);
  REQUIRE(std::ranges::size(points.ghost_indices_offset()) == 3);
  CHECK(indices_of(points.ghost_indices_offset()) == std::vector<int>{0, 1, 2});

  /* Positive control: the closed-interval bound the ranges used to have,
     `size() - 1`, drops the last element, and the comparison sees it. */
  CHECK(
      indices_of(std::ranges::iota_view{0, points.size() - 1}) != every_point);
}

/* An empty Entity -- an empty rank's, say -- yields no indices at all. */
TEST_CASE(
    "entity index ranges are empty on an empty entity", "[mesh][entity]") {
  Bare_Mesh bare;
  auto &points = bare.mesh.points;
  REQUIRE(points.size() == 0);

  CHECK(std::ranges::size(points.all_indices()) == 0);
  CHECK(points.all_indices().empty());
  CHECK(std::ranges::size(points.local_indices()) == 0);
  CHECK(points.local_indices().empty());
  CHECK(std::ranges::size(points.ghost_indices()) == 0);
  CHECK(points.ghost_indices().empty());
  CHECK(points.ghost_size() == 0);
  CHECK(std::ranges::size(points.ghost_indices_offset()) == 0);
  CHECK(points.ghost_indices_offset().empty());

  /* Positive control: one local and one ghost element make every range
     non-empty, so the checks above see the entity's size. */
  points.resize(1, 2, 1);
  CHECK(indices_of(points.all_indices()) == std::vector<int>{0, 1});
  CHECK(indices_of(points.local_indices()) == std::vector<int>{0});
  CHECK(indices_of(points.ghost_indices()) == std::vector<int>{1});
  CHECK(points.ghost_size() == 1);
  CHECK(indices_of(points.ghost_indices_offset()) == std::vector<int>{0});
}

/* An Entity with no ghosts has no ghost indices, offset or not; a loop over
   ghost_indices_offset() indexing the ghost arrays (sized by resize()'s ghost
   argument, here zero) must not run. */
TEST_CASE(
    "entity ghost ranges are empty on a ghost-free entity", "[mesh][entity]") {
  Bare_Mesh bare;
  auto &points = bare.mesh.points;
  points.resize(4, 4, 0);
  REQUIRE(points.ghost_mask.empty());

  REQUIRE(std::ranges::size(points.all_indices()) == 4);
  CHECK(indices_of(points.all_indices()) == std::vector<int>{0, 1, 2, 3});
  CHECK(indices_of(points.local_indices()) == std::vector<int>{0, 1, 2, 3});
  CHECK(points.ghost_indices().empty());
  CHECK(points.ghost_size() == 0);
  CHECK(std::ranges::size(points.ghost_indices_offset()) == 0);
  CHECK(points.ghost_indices_offset().empty());

  /* Positive control: a single ghost makes both ghost ranges non-empty. */
  points.resize(4, 5, 1);
  CHECK(indices_of(points.ghost_indices()) == std::vector<int>{4});
  CHECK(points.ghost_size() == 1);
  CHECK(indices_of(points.ghost_indices_offset()) == std::vector<int>{0});
}
