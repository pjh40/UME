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
  \file Ume/Decompose.hh

  Turning a zone-to-part map into N per-rank meshes.

  ### What lives where

  Zones, corners and sides are strictly partitioned: each belongs to exactly
  one rank, and a rank imports a layer of its neighbors' zones -- with their
  corners and the sides on the cut -- as *ghosts*, in the trailing part of the
  index range with `mask == 0`.

  Points, edges and faces on a cut are not ghosts.  They are duplicated into
  the real index range of every rank that touches them, one instance is
  designated the SOURCE and the rest are COPYs, and the kernels reconcile them
  with `gathscat`.  This is the distinction the file format was built around
  (`cpy_idx`/`src_pe`/`src_idx` are indexed by *ghost offset*; `myCpys`/`mySrcs`
  carry the shared entities), and a decomposition needs both mechanisms.

  Ghost zones exist even though the shared-point gathscat is what makes the
  kernels correct, because `calc_face_area` writes through `m:s>s2` for every
  real side including one on a cut face, and the audit requires every corner's
  `m:c>p` to be in range.

  ### Local index order

  For every entity kind the local order is `(is_ghost, global index)`
  lexicographic: locals first in ascending global order, then ghosts in
  ascending global order.  It has to be lexicographic and not plain ascending,
  because a ghost zone is an ordinary real zone somewhere else and its global
  index interleaves with the local ones.  Two things follow: `mask == 0` and
  `index >= local_size()` mean the same thing, which the library assumes
  throughout; and at `nparts == 1` the map is the identity, so a one-rank
  decomposition reproduces its input bit for bit.

  ### Why `Decomposition` holds every rank at once

  Building rank R's `mySrcs` needs to know what local index a shared entity has
  on rank P.  `Entity_Layout::local_index(P, g)` answers that without P's mesh
  existing, so `extract` can produce one mesh at a time and the driver never
  holds more than one `SOA_Idx::Mesh`.
*/

#ifndef UME_DECOMPOSE_HH
#define UME_DECOMPOSE_HH 1

#include "Ume/Partition.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include <vector>

namespace Ume {
namespace Part {

//! How far a rank reaches into its neighbors for ghost zones
enum class Ghost_Depth {
  /*! No ghost zones at all.  A cut face then has only one zone behind it and
      its sides have no `m:s>s2` partner, which trips the audit and breaks
      `calc_face_area`; diagnostic use only. */
  NONE,
  //! Every zone sharing a face with a local zone
  FACE_ADJACENT,
  //! Every zone sharing a point with a local zone
  POINT_ADJACENT
};

//! Which instance of a shared point, edge or face is designated the SOURCE
enum class Owner_Rule {
  /*! The lowest-numbered rank holding it.  Deterministic and symmetric, so no
      rank has to negotiate; it does concentrate SOURCE duty on low ranks,
      which is why the alternatives are here to measure against. */
  LOWEST_RANK,
  //! Whichever candidate rank is so far sourcing the fewest of that kind
  LEAST_LOADED,
  //! A hash of the global index, spreading SOURCE duty without regard to load
  HASH
};

struct Decompose_Options {
  Ghost_Depth ghosts{Ghost_Depth::FACE_ADJACENT};
  Owner_Rule owner_rule{Owner_Rule::LOWEST_RANK};
  //! Give ghost zones, corners and sides `myCpys`/`mySrcs` entries too
  bool ghost_entity_comm{true};
  //! Carry the input's subsets through to the partitions
  bool split_subsets{true};
};

//! Where every instance of one entity kind ended up
/*! A CSR of locals and a CSR of ghosts, both by rank and both listing global
    indices in ascending order, plus the SOURCE rank of each global entity.
    Together those are enough to build any rank's mesh, and any rank's view of
    any other rank's mesh. */
struct Entity_Layout {
  //! Locals of rank `r` are `gid[start[r] .. start[r+1])`, ascending
  std::vector<int> start, gid;
  //! Ghosts of rank `r` are `ggid[gstart[r] .. gstart[r+1])`, ascending
  std::vector<int> gstart, ggid;
  //! Global index -> the rank that sources it, or -1 if it is on no rank
  /*! For zones, corners and sides -- which are strictly partitioned -- this is
      simply the rank that holds the entity, so it is also what a ghost's
      `src_pe` should be. */
  std::vector<int> owner;

  int nparts() const { return static_cast<int>(start.size()) - 1; }
  //! The number of global entities of this kind in the serial mesh
  int nglobal() const { return static_cast<int>(owner.size()); }

  int local_size(int r) const;
  int ghost_local_size(int r) const;
  int size(int r) const { return local_size(r) + ghost_local_size(r); }

  //! The local index of global entity `g` on rank `r`, or -1 if it is absent
  int local_index(int r, int g) const;
  //! The global index of local index `l` on rank `r`
  int global_index(int r, int l) const;
};

//! The entity kinds, in the order `Mesh` declares them
/*! This is the order `Decomposition::layouts()` returns and the order
    everything that walks the kinds generically -- metrics, the cross-rank
    check, the global-id side-car -- relies on, so it is declared here beside
    the thing that fixes it. */
enum Kind { POINT, EDGE, FACE, ZONE, CORNER, SIDE, IOTA, NKINDS };

//! The singular name of `kind`, for reports and for the side-car
char const *kind_name(int kind);

//! True if `kind` is duplicated across a cut rather than ghosted over it
constexpr bool is_shared_kind(int kind) {
  return kind == POINT || kind == EDGE || kind == FACE;
}

//! Where everything ended up, for every rank at once
struct Decomposition {
  int nparts{0};
  //! The zone map this was built from, unchanged
  Zone_Map z2r;
  Entity_Layout points, edges, faces, zones, corners, sides, iotas;
  Decompose_Options opt;

  //! The layouts in the order `Mesh` declares the entities, for generic code
  std::vector<Entity_Layout const *> layouts() const {
    return {&points, &edges, &faces, &zones, &corners, &sides, &iotas};
  }
};

//! Work out what every rank holds, given a zone-to-part map
/*! `src` must be a serial mesh -- one with no ghost entities of its own.
    Zones that are not real do not appear in `z2r`: a boundary zone follows the
    real zone behind its single face, which is resolved here rather than by a
    `Partitioner`.

    Throws `std::invalid_argument` if `nparts < 1`, if `z2r` is the wrong
    length or names a part outside `[0, nparts)`, if `src` already has ghosts,
    or if a boundary zone cannot be traced back to a real one. */
Decomposition build_decomposition(SOA_Idx::Mesh const &src, Zone_Map z2r,
    int nparts, Decompose_Options const &opt = {});

//! Build rank `pe`'s mesh out of `src` and `d`
/*! `dst` must be a freshly constructed `Mesh` with a Transport attached: the
    Datastore's derived variables are initialized lazily and never invalidated,
    so a mesh that has already been read from is not safe to overwrite.

    Only stored maps are written.  Everything derived -- `m:z>c`, `zcoord`,
    `point_norm` and the rest -- is left for the Datastore to work out from
    them, which is also why `Mesh::write` needs nothing else.

    A cut is not a mesh boundary: points on a cut keep their interior mask and
    no phantom sides are invented for cut faces.  Marking them would make
    `VAR_point_norm` and the boundary branch of `gradzatp` subtract an outward
    normal that is not there, and the gradient would be quietly wrong along
    every partition boundary.

    Throws `std::out_of_range` if `pe` is not a part of `d`, and
    `std::logic_error` if a map cannot be renumbered -- which means the
    decomposition and the mesh disagree, not that the input was bad. */
void extract(SOA_Idx::Mesh const &src, Decomposition const &d, int pe,
    SOA_Idx::Mesh &dst);

} // namespace Part
} // namespace Ume

#endif
