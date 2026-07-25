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
*/

#ifndef UME_MESH_AUDIT_HH
#define UME_MESH_AUDIT_HH 1

#include "Ume/SOA_Idx_Mesh.hh"
#include <iosfwd>

namespace Ume {
namespace SOA_Idx {

//! Which groups of checks `audit` should run
struct Audit_Options {
  //! Connectivity: index ranges, `m:s>s2` pairing, side/corner agreement
  bool topology{true};
  //! Mask invariants, including the ones `gradzatp`/`gradzatz` depend on
  bool masks{true};
  /*! Derived geometry: side volumes, zone volumes, corner surface closure,
      point normals.  This forces the lazy `VAR_*::init_()` derivations, so the
      mesh needs a Transport attached. */
  bool geometry{true};
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
