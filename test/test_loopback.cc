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
  \file test_loopback.cc

  What a decomposition cannot say about itself: does it still compute the same
  physics?

  `test_partition.cc` checks that a partition is *well formed* -- the audit, the
  cross-rank check, the layouts, the metrics -- and a decomposition can be all
  of those and still be wrong.  Until now the only thing that would notice was
  `partition_mpi_compare`, which needs `USE_MPI=ON`, an `mpiexec`, and files on
  disk.  A `Comm::Loopback::Group` runs every rank of the decomposition in this
  process, so the same statement is made in every build.

  The reference is the serial mesh, in memory, and the two answers are lined up
  through `Entity_Layout::global_index` rather than through a `.gid` side-car --
  which is the one thing this can do that the MPI test cannot.

  The tolerance is relative and is not zero.  A partition genuinely reorders the
  accumulation in `gradzatz`, so the answers agree to roundoff and
  `test_brick.cc`'s exact comparison would be the wrong test here.
*/

#include "Ume/Brick_Mesh.hh"
#include "Ume/Comm_Loopback.hh"
#include "Ume/Decompose.hh"
#include "Ume/Partition.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include "Ume/gradient.hh"
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using namespace Ume;
using Ume::SOA_Idx::Mesh;
using Group = Ume::Comm::Loopback::Group;
using DBLV_T = Ume::DS_Types::DBLV_T;
using INTV_T = Ume::DS_Types::INTV_T;
using VEC3V_T = Ume::DS_Types::VEC3V_T;
using VEC3_T = Ume::DS_Types::VEC3_T;

namespace {

Brick::Spec mixed_spec() {
  Brick::Spec spec;
  spec.nx = 4;
  spec.ny = 3;
  spec.nz = 2;
  spec.tet_fraction = 0.5;
  return spec;
}

double norm(VEC3_T const &v) { return std::sqrt(dotprod(v, v)); }

//! The real zones of `mesh`, in local index order
std::vector<int> real_zones(Mesh const &mesh) {
  std::vector<int> v;
  for (int z = 0; z < mesh.zones.local_size(); ++z)
    if (mesh.zones.mask[static_cast<std::size_t>(z)] >= 1)
      v.push_back(z);
  return v;
}

//! One decomposition of one brick, with every rank extracted and ready to run
/*! The ranks are held by pointer because a `Mesh` is pinned to its
    construction address, and every one of them has a `Loopback` Transport
    attached rather than the serial no-op one -- so the moment any of their
    derived variables is touched, it has to be inside `group.run`. */
struct Ranks {
  Mesh serial;
  Part::Decomposition d;
  Group group;
  std::vector<std::unique_ptr<Mesh>> parts;

  Ranks(Brick::Spec const &spec, std::string_view method, int const nparts)
      : d{split(spec, method, nparts)}, group{nparts} {
    for (int r = 0; r < nparts; ++r) {
      parts.push_back(std::make_unique<Mesh>());
      parts.back()->comm = &group[r];
      Part::extract(serial, d, r, *parts.back());
    }
  }

  int nparts() const { return d.nparts; }
  Mesh &operator[](int const r) { return *parts[static_cast<std::size_t>(r)]; }

private:
  //! Build the brick into `serial` and decompose it
  Part::Decomposition split(
      Brick::Spec const &spec, std::string_view method, int const nparts) {
    Ume::attach_serial_transport(serial);
    Brick::build(spec, serial);
    Part::Options opt;
    opt.nparts = nparts;
    Part::Zone_Map z2r;
    Part::make(method)->assign(serial, opt, z2r);
    return Part::build_decomposition(serial, std::move(z2r), nparts);
  }
};

//! The zone-centered gradient of a field that is `spike` in one zone alone
/*! `spiked` is a local zone index, or -1 on a rank that does not hold the
    zone being perturbed -- which still has to take the gradient, because its
    corners contribute to points that the spiked zone's rank owns. */
VEC3V_T gradient_of_a_spike(Mesh &mesh, int const spiked) {
  DBLV_T zfield(static_cast<std::size_t>(mesh.zones.size()), 0.0);
  if (spiked >= 0)
    zfield[static_cast<std::size_t>(spiked)] = 100000.0;
  VEC3V_T zgrad, pgrad;
  Ume::gradzatz(mesh, zfield, zgrad, pgrad);
  return zgrad;
}

} // namespace

/* --------------------------- the transport alone -------------------------- */

TEST_CASE("loopback: a rank cannot communicate outside a run", "[loopback]") {
  Group group(2);
  Ume::Comm::Neighbors none;
  Ume::Comm::Buffers<INTV_T> sends(none), recvs(none);

  /* A Loopback rank has no idea what the other ranks are doing outside `run`,
     and there is nothing it can usefully wait for.  Saying so beats
     deadlocking -- and this is the same diagnostic a mesh gets if one of its
     derived variables is touched before the run starts. */
  CHECK_THROWS_AS(group[0].exchange(sends, recvs), std::logic_error);
}

TEST_CASE("loopback: a group of one has nobody to talk to", "[loopback]") {
  Group group(1);
  CHECK(group.nranks() == 1);
  CHECK(group[0].id() == 0);
  CHECK_THROWS_AS(group[1], std::out_of_range);
  CHECK_THROWS_AS(Group(0), std::invalid_argument);

  int visited = -1;
  group.run([&](int const r) { visited = r; });
  CHECK(visited == 0);
}

TEST_CASE(
    "loopback: a gathered sum reaches every copy of a point", "[loopback]") {
  /* The pattern `ume_mpi`'s `test_point_gathscat` uses, and the cheapest
     statement that the pairing between one rank's copies and another's sources
     is positionally right: a copy holds -1, a source holds the number of copies
     of it there are, so a SUM leaves every shared point at zero and every
     private one at the background value.  `linear` is the method here because
     contiguous ranges of zone index cut far more faces than a slab does. */
  Ranks ranks(mixed_spec(), "linear", 4);
  int const background = 100;

  std::vector<INTV_T> field(static_cast<std::size_t>(ranks.nparts()));
  ranks.group.run([&](int const r) {
    Mesh &mesh = ranks[r];
    INTV_T &f = field[static_cast<std::size_t>(r)];
    f.assign(static_cast<std::size_t>(mesh.points.size()), background);
    for (auto const &n : mesh.points.myCpys)
      for (int const e : n.elements)
        f[static_cast<std::size_t>(e)] = -1;
    std::map<int, int> copies;
    for (auto const &n : mesh.points.mySrcs)
      for (int const e : n.elements)
        copies[e] += 1;
    for (auto const &[e, count] : copies)
      f[static_cast<std::size_t>(e)] = count;

    mesh.points.gathscat(Ume::Comm::Op::SUM, f);
  });

  int shared = 0;
  for (int r = 0; r < ranks.nparts(); ++r) {
    INFO("rank " << r);
    Mesh const &mesh = ranks[r];
    INTV_T const &f = field[static_cast<std::size_t>(r)];
    for (int p = 0; p < mesh.points.local_size(); ++p) {
      auto const type = mesh.points.comm_type[static_cast<std::size_t>(p)];
      int const got = f[static_cast<std::size_t>(p)];
      if (type == Ume::SOA_Idx::Entity::INTERNAL) {
        CHECK(got == background);
      } else {
        CHECK(got == 0);
        shared += 1;
      }
    }
  }
  /* If the decomposition had no shared points at all the checks above would
     pass without having looked at the transport. */
  CHECK(shared > 0);
}

/* ------------------------- the numerical equivalence ---------------------- */

TEST_CASE(
    "loopback: a partitioned gradient is the serial gradient", "[loopback]") {
  /* Every method this build can run, at part counts that divide the brick and
     part counts that do not.  A method joins this sweep by being registered,
     exactly as in `test_partition.cc`. */
  for (std::string const &method : Part::methods()) {
    for (int const nparts : {1, 2, 3, 5}) {
      INFO(method << " nparts=" << nparts);
      Ranks ranks(mixed_spec(), method, nparts);

      /* The reference, and the zone both runs perturb.  No rank could have
         picked this zone locally -- `interior[interior.size() / 2]` names a
         different zone in every decomposition -- so it is chosen here, in
         serial index space, and each rank looks it up in the layout. */
      std::vector<int> const interior = real_zones(ranks.serial);
      REQUIRE(!interior.empty());
      int const spiked = interior[interior.size() / 2];
      VEC3V_T const want = gradient_of_a_spike(ranks.serial, spiked);

      double scale = 0.0;
      for (int const z : interior)
        scale = std::max(scale, norm(want[static_cast<std::size_t>(z)]));
      REQUIRE(scale > 0.0);

      std::vector<VEC3V_T> got(static_cast<std::size_t>(nparts));
      ranks.group.run([&](int const r) {
        Mesh &mesh = ranks[r];
        int local = -1;
        for (int z = 0; z < mesh.zones.local_size(); ++z)
          if (ranks.d.zones.global_index(r, z) == spiked &&
              mesh.zones.mask[static_cast<std::size_t>(z)] >= 1) {
            local = z;
            break;
          }
        got[static_cast<std::size_t>(r)] = gradient_of_a_spike(mesh, local);
      });

      /* Every real zone of every rank, against the serial answer for the zone
         it came from.  That the ranks between them cover every zone of the
         serial mesh is `check_partitions`'s job and is not re-established
         here. */
      double worst = 0.0;
      int compared = 0;
      for (int r = 0; r < nparts; ++r) {
        Mesh const &mesh = ranks[r];
        VEC3V_T const &mine = got[static_cast<std::size_t>(r)];
        for (int const z : real_zones(mesh)) {
          int const g = ranks.d.zones.global_index(r, z);
          worst = std::max(worst,
              norm(mine[static_cast<std::size_t>(z)] -
                  want[static_cast<std::size_t>(g)]));
          compared += 1;
        }
      }
      CHECK(compared == static_cast<int>(interior.size()));
      /* Not an exact comparison even at one rank, where the partition is the
         input mesh: `gradzatz` accumulates its point sums with
         `Kokkos::atomic_add` in every backend but Serial, so the order is not
         fixed and the *reference* is only reproducible to roundoff either. */
      CHECK(worst <= 1e-12 * scale);
    }
  }
}

/* ---------------------------- misuse is reported -------------------------- */

TEST_CASE("loopback: a rank that leaves an exchange is reported, not hung",
    "[loopback]") {
  /* The failure mode a rendezvous transport has and MPI does not: `exchange`
     is collective, and a rank that returns early leaves the others waiting on
     a message that will never be posted.  It has to come back as a thrown
     error, because the alternative is a test suite that hangs. */
  Ranks ranks(mixed_spec(), "slab-x", 2);
  CHECK_THROWS_AS(ranks.group.run([&](int const r) {
    if (r == 0)
      return;
    DBLV_T f(static_cast<std::size_t>(ranks[r].points.size()), 1.0);
    ranks[r].points.gathscat(Ume::Comm::Op::SUM, f);
  }),
      std::runtime_error);
}

TEST_CASE("loopback: what a body throws comes out of run", "[loopback]") {
  Group group(3);
  /* The rank that throws is not the first one to run, so the ones after it
     have to be released by its unwinding rather than by their own progress. */
  CHECK_THROWS_AS(group.run([](int const r) {
    if (r == 1)
      throw std::domain_error("no");
  }),
      std::domain_error);

  /* And the Group is reusable afterwards: `run` is not a one-shot. */
  int ran = 0;
  group.run([&](int const) { ran += 1; });
  CHECK(ran == 3);
}

TEST_CASE("loopback: an aborting rank throws instead of killing the process",
    "[loopback]") {
  Group group(2);
  CHECK_THROWS_AS(group.run([&](int const r) {
    if (r == 0)
      group[r].abort("something a kernel could not cope with");
  }),
      std::runtime_error);
}
