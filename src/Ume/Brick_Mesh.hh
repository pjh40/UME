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
  \file Ume/Brick_Mesh.hh

  Synthetic brick-mesh generator for testing.

  Generates a structured `nx * ny * nz` lattice of cells, where each *cell* is
  either left as a single hexahedral zone or subdivided into six tetrahedral
  zones (a Kuhn subdivision about the body diagonal).  Two knobs control the
  mixture:

  - `Spec::tet_fraction`: what portion of the mesh is tetrahedral.
  - `Spec::clustering`: how the tetrahedral cells are arranged, from a uniform
    random scatter (0.0) through to a single face-connected component (1.0).

  Mixing the two cell types is only possible because Ume's side/corner
  decomposition supports general polyhedra: a hex cell that neighbors a tet
  cell has the shared quadrilateral face split into two triangles, so that hex
  zone becomes a 7-, 8-, ... faced polyhedron rather than a cube.  The split
  diagonal is the one induced by the Kuhn subdivision, which is translation
  invariant, so the two sides of every shared face always agree.

  The generator only supplies topology, `pcoord`, masks and communication
  structure; all geometry (`side_vol`, `corner_csurf`, `zcoord`, ...) is
  derived lazily by the `VAR_*::init_()` methods in the Datastore.
*/

#ifndef UME_BRICK_MESH_HH
#define UME_BRICK_MESH_HH 1

#include "Ume/SOA_Idx_Mesh.hh"
#include <cstdint>
#include <iosfwd>
#include <vector>

namespace Ume {
namespace Brick {

//! How a single lattice cell is discretized
enum class Cell_Type : std::uint8_t {
  HEX = 0, //!< one hexahedral zone
  TET = 1 //!< six tetrahedral zones (Kuhn subdivision)
};

//! What `Spec::tet_fraction` is a fraction *of*
enum class Fraction_Basis : std::uint8_t {
  CELLS, //!< fraction of lattice cells that are subdivided
  ZONES //!< fraction of resulting zones that are tetrahedra
};

//! How the physical mesh boundary is represented
/*! Either style produces a mesh that every consumer in the library accepts;
    see `doc/brick_mesh.md` for the reasoning behind the default. */
enum class Boundary_Style : std::uint8_t {
  /*! Each boundary face gets a degenerate "boundary zone" outside the mesh
      (`zmask == -1`) owning one zero-volume corner per face point
      (`cmask == -1`) and the phantom sides (`smask == -1`).  This is what the
      library's own consumers imply: `VAR_side_surf` treats a phantom side as
      having no usable `zcoord`, which only holds if `m:s>z` names a degenerate
      zone, and `VAR_zcoord` averages every `cmask != 0` corner, which would
      shift the centroid if phantom corners were added to the *real* zone. */
  BOUNDARY_ZONES,
  /*! No extra zones or corners: phantom sides point back at the adjacent real
      zone and its corners, so the mesh has no `mask < 0` zones or corners at
      all. */
  PHANTOM_SIDES_ONLY
};

//! Description of a brick mesh to generate
struct Spec {
  //! Number of lattice cells in each direction
  int nx{1}, ny{1}, nz{1};
  //! Cell size in each direction
  double dx{1.0}, dy{1.0}, dz{1.0};
  //! Coordinates of the lattice origin
  double x0{0.0}, y0{0.0}, z0{0.0};
  //! Portion of the mesh that is tetrahedral, in [0, 1]
  double tet_fraction{0.0};
  //! Interpretation of `tet_fraction`
  Fraction_Basis basis{Fraction_Basis::ZONES};
  /*! Spatial arrangement of the tetrahedral cells, in [0, 1].  The tet cells
      are grown from `K = round(T^(1 - clustering))` random seeds, where `T` is
      the number of tet cells requested, so 0.0 gives `T` independent cells (a
      uniform random scatter) and 1.0 gives a single face-connected component.
      Intermediate values give a mix of blobs of random size. */
  double clustering{0.0};
  //! Seed for the (deterministic) cell selection
  std::uint64_t seed{0x9E3779B97F4A7C15ULL};
  //! Physical-boundary representation
  Boundary_Style boundary{Boundary_Style::BOUNDARY_ZONES};

  constexpr int num_cells() const { return nx * ny * nz; }
};

//! Summary of a generated mesh
struct Stats {
  int cells{0};
  int tet_cells{0};
  //! Number of face-connected components of the tet-cell set
  int tet_clusters{0};
  int hex_zones{0};
  int tet_zones{0};
  int boundary_zones{0};
  int points{0}, edges{0}, faces{0}, corners{0};
  int boundary_faces{0};
  int real_sides{0}, phantom_sides{0};
  //! Faces of hex zones that had to be split because of a tet neighbor
  int split_faces{0};

  double tet_cell_fraction() const;
  double tet_zone_fraction() const;
  void print(std::ostream &os) const;
};

//! Choose which lattice cells are tetrahedral (deterministic given `spec`)
std::vector<Cell_Type> select_cells(Spec const &spec);

//! Build a brick mesh into `mesh` using an explicit cell-type list
/*! `cells` must have `spec.num_cells()` entries, indexed
    `i + nx * (j + ny * k)`. */
Stats build(
    Spec const &spec, std::vector<Cell_Type> const &cells, SOA_Idx::Mesh &mesh);

//! Build a brick mesh into `mesh`, selecting cells from `spec`
inline Stats build(Spec const &spec, SOA_Idx::Mesh &mesh) {
  return build(spec, select_cells(spec), mesh);
}

//! Check the topological and geometric invariants of a generated brick
/*! Writes a line per failed check to `os`.  Requires a Transport (see
    `Ume::attach_serial_transport`). */
bool validate(SOA_Idx::Mesh &mesh, Spec const &spec, std::ostream &os);

} // namespace Brick
} // namespace Ume

#endif
