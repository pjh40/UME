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
  \file test_brick.cc

  Tests for the synthetic hex/tet brick generator, and the gradient operators
  applied to the meshes it produces.
*/

#include "Ume/Brick_Mesh.hh"
#include "Ume/gradient.hh"
#include "Ume/mesh_audit.hh"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <sstream>

using namespace Ume;
using Ume::SOA_Idx::Mesh;

namespace {

/*! Tolerance for comparing gradients.  On a `UME_SERIAL` build every
    accumulation runs in a fixed order, so the push and pull formulations agree
    bitwise (this is the property the unsorted `m:z>c` ordering exists to
    preserve).  Every other backend accumulates with `Kokkos::atomic_add`,
    whose ordering is not reproducible, so only a tolerance is meaningful. */
#if defined(UME_SERIAL)
constexpr double GRAD_TOL = 0.0;
#else
constexpr double GRAD_TOL = 1e-12;
#endif

//! The specs used by the topology/geometry sweeps
std::vector<Brick::Spec> const &sweep() {
  static std::vector<Brick::Spec> const specs = [] {
    std::vector<Brick::Spec> v;
    auto add = [&v](int nx, int ny, int nz, double frac, double cluster) {
      Brick::Spec s;
      s.nx = nx;
      s.ny = ny;
      s.nz = nz;
      s.tet_fraction = frac;
      s.clustering = cluster;
      v.push_back(s);
    };
    add(1, 1, 1, 0.0, 0.0); // one hexahedron
    add(1, 1, 1, 1.0, 0.0); // one Kuhn-subdivided cell
    add(3, 2, 2, 0.0, 0.0); // pure hex
    add(3, 2, 2, 1.0, 0.0); // pure tet
    add(4, 4, 4, 0.25, 0.0); // scattered tets
    add(4, 4, 4, 0.5, 1.0); // one connected tet blob
    add(4, 3, 5, 0.75, 0.5); // a mix of blobs
    return v;
  }();
  return specs;
}

//! The mixed brick the audit tests corrupt
Brick::Spec audit_spec() {
  Brick::Spec spec;
  spec.nx = 3;
  spec.ny = 3;
  spec.nz = 2;
  spec.tet_fraction = 0.4;
  return spec;
}

/*! Give `mesh` the iota entity the generator does not build: one tetrahedron
    `(zone center, face center, edge center, point)` per side endpoint. */
void add_iotas(Mesh &mesh) {
  int const sll = mesh.sides.size();
  mesh.iotas.resize(2 * sll, 2 * sll, 0);
  std::fill(mesh.iotas.mask.begin(), mesh.iotas.mask.end(), short{1});
  std::fill(mesh.iotas.comm_type.begin(), mesh.iotas.comm_type.end(),
      static_cast<int>(SOA_Idx::Entity::INTERNAL));

  auto const &s2z = mesh.ds->caccess_intv("m:s>z");
  auto const &s2f = mesh.ds->caccess_intv("m:s>f");
  auto const &s2e = mesh.ds->caccess_intv("m:s>e");
  auto const &s2p1 = mesh.ds->caccess_intv("m:s>p1");
  auto const &s2p2 = mesh.ds->caccess_intv("m:s>p2");
  auto &a2z = mesh.ds->access_intv("m:a>z");
  auto &a2f = mesh.ds->access_intv("m:a>f");
  auto &a2e = mesh.ds->access_intv("m:a>e");
  auto &a2p = mesh.ds->access_intv("m:a>p");
  auto &a2s = mesh.ds->access_intv("m:a>s");

  for (int s = 0; s < sll; ++s)
    for (int h = 0; h < 2; ++h) {
      int const a = 2 * s + h;
      a2z[a] = s2z[s];
      a2f[a] = s2f[s];
      a2e[a] = s2e[s];
      a2p[a] = h ? s2p2[s] : s2p1[s];
      a2s[a] = s;
    }
}

std::string describe(Brick::Spec const &s) {
  std::ostringstream os;
  os << s.nx << 'x' << s.ny << 'x' << s.nz << " tet_fraction=" << s.tet_fraction
     << " clustering=" << s.clustering;
  return os.str();
}

//! Face-connected components of the tetrahedral cells
int clusters(Brick::Spec const &spec, std::vector<Brick::Cell_Type> const &c) {
  Mesh mesh;
  Brick::attach_serial_transport(mesh);
  return Brick::build(spec, c, mesh).tet_clusters;
}

} // namespace

TEST_CASE("brick: cell selection honors the requested fraction", "[brick]") {
  Brick::Spec spec;
  spec.nx = spec.ny = spec.nz = 8;
  spec.basis = Brick::Fraction_Basis::CELLS;

  for (double f : {0.0, 0.1, 0.25, 0.5, 1.0}) {
    spec.tet_fraction = f;
    auto const cells = Brick::select_cells(spec);
    int const tets = static_cast<int>(
        std::count(cells.begin(), cells.end(), Brick::Cell_Type::TET));
    INFO("cell fraction " << f);
    CHECK(tets == std::llround(f * spec.num_cells()));
  }

  /* A zone fraction t needs a cell fraction t / (6 - 5t), because each
     subdivided cell contributes six zones instead of one. */
  spec.basis = Brick::Fraction_Basis::ZONES;
  spec.tet_fraction = 0.5;
  auto const cells = Brick::select_cells(spec);
  int const tets = static_cast<int>(
      std::count(cells.begin(), cells.end(), Brick::Cell_Type::TET));
  int const zones = 6 * tets + (spec.num_cells() - tets);
  CHECK_THAT(static_cast<double>(6 * tets) / zones,
      Catch::Matchers::WithinAbs(0.5, 0.01));
}

TEST_CASE("brick: clustering controls the number of components", "[brick]") {
  Brick::Spec spec;
  spec.nx = spec.ny = spec.nz = 8;
  spec.basis = Brick::Fraction_Basis::CELLS;
  spec.tet_fraction = 0.2; // 102 cells

  spec.clustering = 1.0;
  CHECK(clusters(spec, Brick::select_cells(spec)) == 1);

  /* A pure scatter puts every cell down independently; some land next to each
     other, but nothing like a single blob. */
  spec.clustering = 0.0;
  int const scattered = clusters(spec, Brick::select_cells(spec));
  CHECK(scattered > 20);

  spec.clustering = 0.5;
  int const mixed = clusters(spec, Brick::select_cells(spec));
  CHECK(mixed > 1);
  CHECK(mixed < scattered);
}

TEST_CASE("brick: generated meshes satisfy the mesh invariants", "[brick]") {
  for (auto const &spec : sweep()) {
    INFO(describe(spec));
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::Stats const stats = Brick::build(spec, mesh);
    CHECK(stats.hex_zones + stats.tet_zones ==
        stats.cells - stats.tet_cells + 6 * stats.tet_cells);
    std::ostringstream report;
    bool const valid = Brick::validate(mesh, spec, report);
    INFO(report.str());
    CHECK(valid);
  }
}

TEST_CASE("brick: the phantom-side boundary style also validates", "[brick]") {
  Brick::Spec spec;
  spec.nx = 3;
  spec.ny = 3;
  spec.nz = 2;
  spec.tet_fraction = 0.5;
  spec.boundary = Brick::Boundary_Style::PHANTOM_SIDES_ONLY;
  Mesh mesh;
  Brick::attach_serial_transport(mesh);
  Brick::Stats const stats = Brick::build(spec, mesh);
  CHECK(stats.boundary_zones == 0);
  std::ostringstream report;
  bool const valid = Brick::validate(mesh, spec, report);
  INFO(report.str());
  CHECK(valid);
}

TEST_CASE("brick: survives a write/read round trip", "[brick]") {
  Brick::Spec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.tet_fraction = 0.5;
  spec.clustering = 0.5;

  Mesh mesh;
  Brick::attach_serial_transport(mesh);
  Brick::build(spec, mesh);

  std::stringstream buf;
  mesh.write(buf);
  Mesh copy;
  copy.read(buf);
  CHECK(mesh == copy);
}

TEST_CASE("brick: gradient of a constant field vanishes", "[brick][gradient]") {
  for (auto const &spec : sweep()) {
    INFO(describe(spec));
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(spec, mesh);

    DS_Types::DBLV_T zone_field(static_cast<size_t>(mesh.zones.size()), 3.25);
    DS_Types::VEC3V_T zgrad, pgrad;
    gradzatz(mesh, zone_field, zgrad, pgrad);

    double worst{0.0};
    for (int p = 0; p < mesh.points.local_size(); ++p)
      worst = std::max(worst, vectormag(pgrad[p]));
    for (int z = 0; z < mesh.zones.local_size(); ++z)
      worst = std::max(worst, vectormag(zgrad[z]));
    CHECK(worst < 1e-13);
  }
}

/*! The corner gradient reproduces a linear field exactly only where the
    stencil around a point is symmetric enough for the per-zone facet moments
    to cancel: on a uniform lattice of one cell type they do, but across a
    hex/tet interface they do not.  So exactness is required only at points
    whose surrounding cells are all of the same type, which is still every
    point of a pure-hex or pure-tet brick. */
TEST_CASE("brick: gradient of a linear field is exact on a uniform stencil",
    "[brick][gradient]") {
  Vec3 const slope{{0.5, -1.25, 2.0}};

  for (auto const &spec : sweep()) {
    INFO(describe(spec));
    auto const cells = Brick::select_cells(spec);
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(spec, cells, mesh);

    auto const &zcoord = mesh.ds->caccess_vec3v("zcoord");
    auto const &zmask = mesh.zones.mask;
    auto const &pmask = mesh.points.mask;
    auto const &z2c = mesh.ds->caccess_intrr("m:z>c");
    auto const &c2p = mesh.ds->caccess_intv("m:c>p");

    //! Is every cell touching this point of the same type?
    auto uniform_stencil = [&](int const p) {
      int const i = p % (spec.nx + 1);
      int const j = (p / (spec.nx + 1)) % (spec.ny + 1);
      int const k = p / ((spec.nx + 1) * (spec.ny + 1));
      Brick::Cell_Type first = cells[0];
      bool seen{false}, same{true};
      for (int dk = -1; dk < 1; ++dk)
        for (int dj = -1; dj < 1; ++dj)
          for (int di = -1; di < 1; ++di) {
            int const ci = i + di, cj = j + dj, ck = k + dk;
            if (ci < 0 || cj < 0 || ck < 0 || ci >= spec.nx || cj >= spec.ny ||
                ck >= spec.nz)
              continue;
            auto const t =
                cells[static_cast<size_t>(ci + spec.nx * (cj + spec.ny * ck))];
            if (!seen) {
              first = t;
              seen = true;
            }
            same = same && t == first;
          }
      return same;
    };

    DS_Types::DBLV_T zone_field(static_cast<size_t>(mesh.zones.size()), 0.0);
    for (int z = 0; z < mesh.zones.size(); ++z)
      zone_field[z] = dotprod(slope, zcoord[z]);

    DS_Types::VEC3V_T zgrad, pgrad;
    gradzatz(mesh, zone_field, zgrad, pgrad);

    int checked{0};
    for (int p = 0; p < mesh.points.local_size(); ++p) {
      if (pmask[p] < 1)
        continue; // the boundary treatment is a projection, not a gradient
      if (!uniform_stencil(p))
        continue;
      INFO("point " << p);
      REQUIRE(vectormag(pgrad[p] - slope) < 1e-11);
      ++checked;
    }
    for (int z = 0; z < mesh.zones.local_size(); ++z) {
      if (zmask[z] < 1)
        continue;
      bool ok{true};
      for (int const c : z2c[z])
        ok = ok && pmask[c2p[c]] >= 1 && uniform_stencil(c2p[c]);
      if (!ok)
        continue;
      INFO("zone " << z);
      REQUIRE(vectormag(zgrad[z] - slope) < 1e-11);
      ++checked;
    }
    if (spec.nx > 1 && spec.ny > 1 && spec.nz > 1)
      CHECK(checked > 0); // a 1x1x1 brick has no interior points at all
  }
}

TEST_CASE("brick: gradzatz and gradzatz_invert agree", "[brick][gradient]") {
  for (auto const &spec : sweep()) {
    INFO(describe(spec));
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(spec, mesh);

    /* Put a delta on a zone that touches the physical boundary, so that the
       comparison actually exercises the boundary stencil. */
    auto const &pmask = mesh.points.mask;
    auto const &zmask = mesh.zones.mask;
    auto const &z2c = mesh.ds->caccess_intrr("m:z>c");
    auto const &c2p = mesh.ds->caccess_intv("m:c>p");
    int delta_zone{-1};
    for (int z = 0; z < mesh.zones.local_size() && delta_zone < 0; ++z) {
      if (zmask[z] < 1)
        continue;
      for (int const c : z2c[z])
        if (pmask[c2p[c]] < 0) {
          delta_zone = z;
          break;
        }
    }
    REQUIRE(delta_zone >= 0);

    DS_Types::DBLV_T zone_field(static_cast<size_t>(mesh.zones.size()), 0.0);
    zone_field[delta_zone] = 1.0;

    DS_Types::VEC3V_T zgrad, pgrad, zgrad_inv, pgrad_inv;
    gradzatz(mesh, zone_field, zgrad, pgrad);
    gradzatz_invert(mesh, zone_field, zgrad_inv, pgrad_inv);

    for (int p = 0; p < mesh.points.local_size(); ++p) {
      INFO("point " << p);
      REQUIRE(vectormag(pgrad[p] - pgrad_inv[p]) <= GRAD_TOL);
    }
    for (int z = 0; z < mesh.zones.local_size(); ++z) {
      INFO("zone " << z);
      REQUIRE(vectormag(zgrad[z] - zgrad_inv[z]) <= GRAD_TOL);
    }
  }
}

TEST_CASE("audit: a generated brick is clean", "[brick][audit]") {
  for (auto const &spec : sweep()) {
    INFO(describe(spec));
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(spec, mesh);

    std::ostringstream report;
    SOA_Idx::Audit_Options opts;
    opts.census = true;
    SOA_Idx::Audit_Result const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK(res.ok());
    CHECK(res.problems == 0);
    CHECK(res.checks > 30);
  }
}

TEST_CASE("audit: a facet whose second side is a ghost still has two sides",
    "[brick][audit]") {
  /* A partition boundary cuts a face: the zone on the far side is imported as
     a ghost, and so are its sides, which carry mask 0.  The facet still has
     exactly two sides -- one followable, one not -- and the side-count check
     has to see both, or every cut face in a decomposed mesh reports a
     spurious "has 1 sides". */
  Mesh mesh;
  Brick::attach_serial_transport(mesh);
  Brick::build(audit_spec(), mesh);

  /* Ghost a whole zone, the way a decomposition does: the zone, every side
     that names it, and every corner that names it all go to mask 0
     together. */
  auto const &s2z = mesh.ds->caccess_intv("m:s>z");
  auto const &c2z = mesh.ds->caccess_intv("m:c>z");
  int ghost_zone = -1;
  for (int z = 0; z < mesh.zones.size(); ++z)
    if (mesh.zones.mask[static_cast<size_t>(z)] > 0) {
      ghost_zone = z;
      break;
    }
  REQUIRE(ghost_zone >= 0);
  mesh.zones.mask[static_cast<size_t>(ghost_zone)] = 0;
  for (int s = 0; s < mesh.sides.size(); ++s)
    if (s2z[static_cast<size_t>(s)] == ghost_zone)
      mesh.sides.mask[static_cast<size_t>(s)] = 0;
  for (int c = 0; c < mesh.corners.size(); ++c)
    if (c2z[static_cast<size_t>(c)] == ghost_zone)
      mesh.corners.mask[static_cast<size_t>(c)] = 0;

  SOA_Idx::Audit_Options opts;
  opts.geometry = false;
  opts.masks = false;
  opts.derived = false;
  opts.parallel = false;

  std::ostringstream report;
  [[maybe_unused]] auto const res = SOA_Idx::audit(mesh, report, opts);
  INFO(report.str());
  CHECK_THAT(report.str(),
      !Catch::Matchers::ContainsSubstring("(face, edge) side counts"));
}

TEST_CASE("audit: the iota maps are checked", "[brick][audit]") {
  SECTION("a consistent set of iotas passes") {
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(audit_spec(), mesh);
    add_iotas(mesh);

    std::ostringstream report;
    SOA_Idx::Audit_Options opts;
    opts.census = true;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK(res.ok());
    CHECK_THAT(report.str(), Catch::Matchers::ContainsSubstring("iotas:"));
  }

  SECTION("an iota at a point its side does not touch") {
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(audit_spec(), mesh);
    add_iotas(mesh);
    auto &a2p = mesh.ds->access_intv("m:a>p");
    a2p[0] = a2p[1]; // both endpoints of side 0 are now the same point

    std::ostringstream report;
    auto const res = SOA_Idx::audit(mesh, report, {});
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("iotas "
                                           "per side"));
  }

  SECTION("an iota that disagrees with its side") {
    Mesh mesh;
    Brick::attach_serial_transport(mesh);
    Brick::build(audit_spec(), mesh);
    add_iotas(mesh);
    auto &a2f = mesh.ds->access_intv("m:a>f");
    a2f[0] = a2f[0] == 0 ? 1 : 0;

    std::ostringstream report;
    auto const res = SOA_Idx::audit(mesh, report, {});
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("m:a> maps agree with m:a>s"));
  }
}

TEST_CASE("audit: corruption is detected", "[brick][audit]") {
  Brick::Spec const spec = audit_spec();

  Mesh mesh;
  Brick::attach_serial_transport(mesh);
  Brick::build(spec, mesh);

  //! Topology only: most of the corruptions below are in the stored maps
  SOA_Idx::Audit_Options opts;
  opts.geometry = false;
  opts.masks = false;
  opts.derived = false;
  opts.parallel = false;

  std::ostringstream report;

  SECTION("broken m:s>s2 involution") {
    auto &s2s2 = mesh.ds->access_intv("m:s>s2");
    s2s2[0] = static_cast<int>(s2s2.size()) - 1;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("m:s>s2 pairing"));
  }

  SECTION("a facet that does not have two sides") {
    /* Moving a side to another face leaves its old (face, edge) with one
       side and its new one with three. */
    auto &s2f = mesh.ds->access_intv("m:s>f");
    s2f[0] = s2f[0] == 1 ? 2 : 1;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("(face, edge) side counts"));
  }

  SECTION("m:s>c1 pointing at the wrong corner") {
    auto &s2c1 = mesh.ds->access_intv("m:s>c1");
    s2c1[0] = s2c1[0] == 0 ? 1 : 0;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("m:s>c1/c2 agree"));
  }

  SECTION("out-of-range map entry") {
    mesh.ds->access_intv("m:s>z")[3] = mesh.zones.size() + 100;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("side map ranges"));
  }

  SECTION("a side whose points are not the points of its edge") {
    auto &s2e = mesh.ds->access_intv("m:s>e");
    s2e[0] = s2e[0] == 0 ? 1 : 0;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("m:s>e agrees with m:s>p1/p2"));
  }

  SECTION("a face ring that does not close") {
    auto &s2s3 = mesh.ds->access_intv("m:s>s3");
    s2s3[0] = s2s3[s2s3[0]]; // skip a side of the ring
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(), Catch::Matchers::ContainsSubstring("m:s>s3"));
  }

  SECTION("a side that is its own m:s>s5") {
    mesh.ds->access_intv("m:s>s5")[0] = 0;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("m:s>s5 pairing"));
  }

  SECTION("a face that does not name the zone of its sides") {
    auto const &s2f = mesh.ds->caccess_intv("m:s>f");
    auto &f2z1 = mesh.ds->access_intv("m:f>z1");
    auto const &f2z2 = mesh.ds->caccess_intv("m:f>z2");
    int const f = s2f[0];
    int const other = mesh.zones.size() - 1;
    f2z1[f] = (other == f2z2[f]) ? other - 1 : other;
    auto const res = SOA_Idx::audit(mesh, report, opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("sides agree with m:f>z1/z2"));
  }

  SECTION("a stale derived map") {
    /* Force the zone-to-corner map, then move a corner to another zone.  The
       derived map is cached, so the operators would keep using the old one. */
    [[maybe_unused]] auto const &z2c = mesh.ds->caccess_intrr("m:z>c");
    auto &c2z = mesh.ds->access_intv("m:c>z");
    c2z[0] = mesh.zones.size() - 1;

    SOA_Idx::Audit_Options derived_opts;
    derived_opts.topology = false;
    derived_opts.masks = false;
    derived_opts.geometry = false;
    derived_opts.parallel = false;
    auto const res = SOA_Idx::audit(mesh, report, derived_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("m:z>c inverts m:c>z"));
  }

  /* A map that would index out of bounds has to disable the groups that would
     do the indexing, rather than let them run into undefined behavior. */
  SECTION("an unusable map disables the derived and geometry groups") {
    mesh.ds->access_intv("m:c>p")[0] = -1;
    auto const res = SOA_Idx::audit(mesh, report, {});
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("derived maps are index-safe"));
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("[derived] skipped"));
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("[geometry] skipped"));
  }

  SECTION("a point with no real corner starves gradzatp") {
    /* Demote every corner at one point to a boundary corner.  This is what a
       "a corner is a boundary corner whenever its point is on the boundary"
       mask convention would produce, and it makes gradzatp's point_volume
       zero at that point. */
    auto const &c2p = mesh.ds->caccess_intv("m:c>p");
    int const victim = c2p[0];
    for (int c = 0; c < mesh.corners.local_size(); ++c)
      if (c2p[c] == victim)
        mesh.corners.mask[c] = -1;

    SOA_Idx::Audit_Options mask_opts;
    mask_opts.topology = false;
    mask_opts.geometry = false;
    mask_opts.derived = false;
    auto const res = SOA_Idx::audit(mesh, report, mask_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("divide by zero"));
  }

  SECTION("a face with no live side starves fcoord") {
    auto const &s2f = mesh.ds->caccess_intv("m:s>f");
    int const victim = s2f[0];
    for (int s = 0; s < mesh.sides.local_size(); ++s)
      if (s2f[s] == victim)
        mesh.sides.mask[s] = 0;

    SOA_Idx::Audit_Options mask_opts;
    mask_opts.topology = false;
    mask_opts.geometry = false;
    mask_opts.derived = false;
    auto const res = SOA_Idx::audit(mesh, report, mask_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(),
        Catch::Matchers::ContainsSubstring("fcoord will divide by zero"));
  }

  SECTION("a live entity in the ghost index range") {
    mesh.zones.mask[0] = 0;
    SOA_Idx::Audit_Options mask_opts;
    mask_opts.topology = false;
    mask_opts.geometry = false;
    mask_opts.derived = false;
    auto const res = SOA_Idx::audit(mesh, report, mask_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("below local_size"));
  }

  SECTION("a copy record that names this rank") {
    mesh.points.cpy_idx.push_back(0);
    mesh.points.src_pe.push_back(mesh.mype);
    mesh.points.src_idx.push_back(0);
    mesh.points.ghost_mask.push_back(0);

    SOA_Idx::Audit_Options comm_opts;
    comm_opts.topology = false;
    comm_opts.masks = false;
    comm_opts.geometry = false;
    comm_opts.derived = false;
    auto const res = SOA_Idx::audit(mesh, report, comm_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(
        report.str(), Catch::Matchers::ContainsSubstring("copy records"));
  }

  SECTION("an invalid communication type") {
    mesh.sides.comm_type[0] = 42;

    SOA_Idx::Audit_Options comm_opts;
    comm_opts.topology = false;
    comm_opts.masks = false;
    comm_opts.geometry = false;
    comm_opts.derived = false;
    auto const res = SOA_Idx::audit(mesh, report, comm_opts);
    INFO(report.str());
    CHECK_FALSE(res.ok());
    CHECK_THAT(report.str(), Catch::Matchers::ContainsSubstring("comm types"));
  }
}
