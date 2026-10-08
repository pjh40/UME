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
#include <catch2/catch_test_macros.hpp>

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
