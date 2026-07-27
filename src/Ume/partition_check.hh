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
  \file Ume/partition_check.hh

  Consistency checks that need more than one partition to see.

  `mesh_audit` reads one mesh and can only ask whether its parallel bookkeeping
  is well formed -- indices in range, comm types in the enumeration, no rank
  naming itself.  It cannot ask the question that actually matters, which is
  whether the ranks agree with each other, because the entity a `src_idx`
  names lives in a mesh the audit never sees.  That is what these checks are
  for, and they are the reason a decomposition can be believed.

  Four things are checked, in increasing order of how badly they hurt when
  wrong:

  - **The neighbor lists pair up.**  Every `myCpys` entry naming rank P has
    exactly one `mySrcs` counterpart on P naming this rank, of the same
    length, and the entities it names are COPY or GHOST here and SOURCE there.
  - **The pairing is positional.**  `myCpys[P].elements[j]` here and
    `mySrcs[R].elements[j]` there are the *same entity*.  `Comm::Buffers`
    packs and unpacks by position, so an agreement on lengths with a
    disagreement on order is a silent permutation of the data on the wire --
    the one failure MPI would turn into quietly wrong physics rather than a
    crash, and the only one worth building this file for.
  - **Ghost sources are real.**  A ghost's `{src_pe, src_idx}` names an entity
    in the real index range of that rank, which is that same entity.
  - **The ranks reassemble.**  Counting each entity once -- on its SOURCE rank
    for a shared entity, on its owner for a ghosted one -- reproduces the
    serial mesh exactly, with nothing lost and nothing counted twice.

  ### Identifying an entity across two meshes

  Nothing in the file format carries a global index, so "the same entity" has
  to be answered geometrically.  Every check here identifies an entity by the
  set of point coordinates it is built on, taken from the stored maps and
  `pcoord` alone:

  | kind | coordinates |
  |---|---|
  | point | its own |
  | edge, side | its two endpoints |
  | corner, iota | its point |
  | face | the endpoints of the real sides on it |
  | zone | the points of its corners |

  Those are bit-identical between ranks, because `extract` copies `pcoord`
  rather than recomputing it, so the comparisons are exact and no tolerance
  is involved.  Derived geometry (`fcoord`, `zcoord`) is deliberately not used:
  it is filled on ghosts by a scatter that has not happened here, so it would
  compare garbage against a real value.

  A ghost is compared by containment rather than equality.  A rank imports
  only the part of a ghost zone that touches points it already holds, so the
  ghost's point set is a subset of the source zone's -- a proper one whenever
  the zone straddles the cut.
*/

#ifndef UME_PARTITION_CHECK_HH
#define UME_PARTITION_CHECK_HH 1

#include "Ume/SOA_Idx_Mesh.hh"
#include <iosfwd>
#include <span>

namespace Ume {
namespace Part {

//! Which groups of cross-rank checks `check_partitions` should run
struct Check_Options {
  //! `myCpys`/`mySrcs` pair up one-to-one, with matching lengths and comm types
  bool neighbors{true};
  /*! Paired positions name the same entity.  Requires `neighbors`, which
      establishes that there is a pairing to walk. */
  bool pairing{true};
  //! Every ghost's `{src_pe, src_idx}` names that same entity, real, elsewhere
  bool ghosts{true};
  //! The ranks put the serial mesh back together; ignored without one
  bool reassembly{true};
  //! Maximum number of offending entities reported per failing check
  int max_reports{5};
};

//! What `check_partitions` found
struct Check_Result {
  //! Number of checks run
  int checks{0};
  //! Number of checks that found at least one problem
  int failed_checks{0};
  //! Total number of offending entities across all checks
  int problems{0};

  constexpr bool ok() const { return failed_checks == 0; }
  explicit constexpr operator bool() const { return ok(); }
};

//! Check that `ranks` are consistent with each other, writing problems to `os`
/*! `ranks` must be indexed by rank: `ranks[r]->mype == r`, and every mesh must
    agree that there are `ranks.size()` of them.  If they do not, that is
    reported and nothing else is run -- every other check indexes `ranks` by a
    rank number read out of a mesh. */
Check_Result check_partitions(std::span<SOA_Idx::Mesh const *const> ranks,
    std::ostream &os, Check_Options const &opts = {});

//! The same, and additionally that `ranks` reassemble into `serial`
Check_Result check_partitions(SOA_Idx::Mesh const &serial,
    std::span<SOA_Idx::Mesh const *const> ranks, std::ostream &os,
    Check_Options const &opts = {});

} // namespace Part
} // namespace Ume

#endif
