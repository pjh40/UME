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
  \file Ume/mesh_audit.hh

  Mesh-independent consistency checks for a `SOA_Idx::Mesh`.

  These hold for any mesh the library is willing to compute on, whatever
  produced it, so the same routine can be pointed at a generated brick or at a
  heritage DELFI dump.  Anything that depends on knowing the intended geometry
  (analytic volumes, expected normals) belongs to the generator instead.

  The audit covers the whole connectivity of the mesh: every stored map
  (`m:e>p*`, `m:f>z*`, `m:c>p`, `m:c>z`, `m:s>*`, `m:a>*`), every map the
  library derives from them (`m:c>s`, `m:z>c`, `m:z>p`, `m:z>pz`, `m:p>zs`,
  `m:p>rc`), the mask conventions those maps are read with, the derived
  geometry, and the parallel bookkeeping in each `Entity`.

  Two caveats apply to a mesh that came out of a decomposition:

  - The audit sees one partition.  Ghost entities (mask 0) carry no local
    connectivity, so checks that would have to follow a map into a neighboring
    partition skip them.
  - Where the file format does not fix a convention -- the direction in which
    `m:s>s2` partners traverse their shared edge, the direction of the `m:s>s3`
    ring around a face -- the audit infers the convention this mesh uses from
    the majority of its sides and reports the ones that disagree.
*/

#ifndef UME_MESH_AUDIT_HH
#define UME_MESH_AUDIT_HH 1

#include "Ume/SOA_Idx_Mesh.hh"
#include <iosfwd>

namespace Ume {
namespace SOA_Idx {

//! Which groups of checks `audit` should run
struct Audit_Options {
  /*! Connectivity: the index ranges of every stored map, and the relations
      between them (`m:s>s2`/`s3`/`s4`/`s5`, side-corner-zone-face agreement,
      the iota maps).  Index-range failures disable the `derived` and
      `geometry` groups, which would otherwise index out of bounds. */
  bool topology{true};
  //! Mask invariants, including the ones `gradzatp`/`gradzatz` depend on
  bool masks{true};
  /*! The maps the library derives from the stored ones (`m:c>s`, `m:z>c`,
      `m:z>p`, `m:z>pz`, `m:p>zs`, `m:p>rc`), checked against their
      definitions.  These are built on the host and need no Transport. */
  bool derived{true};
  /*! Derived geometry: side volumes, zone volumes, corner surface closure,
      point normals.  This forces the lazy `VAR_*::init_()` derivations, so the
      mesh needs a Transport attached. */
  bool geometry{true};
  //! The per-`Entity` parallel bookkeeping: comm types, ghost and copy records
  bool parallel{true};
  //! Write a census of entity counts by mask category
  bool census{false};
  //! Maximum number of offending entities reported per failing check
  int max_reports{5};
};

//! What `audit` found
struct Audit_Result {
  //! Number of checks run
  int checks{0};
  //! Number of checks that found at least one problem
  int failed_checks{0};
  //! Total number of offending entities across all checks
  int problems{0};

  constexpr bool ok() const { return failed_checks == 0; }
  explicit constexpr operator bool() const { return ok(); }
};

//! Check the internal consistency of `mesh`, writing problems to `os`
Audit_Result audit(
    Mesh const &mesh, std::ostream &os, Audit_Options const &opts = {});

} // namespace SOA_Idx
} // namespace Ume

#endif
