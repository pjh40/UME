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
  \file Ume/Partition_Graph.hh

  The mesh views that the partitioning methods share: zone centroids now, and
  the dual graph and hypergraph when the methods that need them land.
*/

#ifndef UME_PARTITION_GRAPH_HH
#define UME_PARTITION_GRAPH_HH 1

#include "Ume/DS_Types.hh"
#include "Ume/SOA_Idx_Mesh.hh"

namespace Ume {
namespace Part {

//! The centroid of every zone, indexed by zone index
/*! This is just `zcoord`, which is derived lazily and whose `init_()` calls
    `zones().scatter()` -- so `mesh` must have a Transport attached (see
    `Ume::attach_serial_transport`) before this is called.  Only real zones
    have a meaningful centroid; the entry for a boundary zone is whatever the
    degenerate zone's points average to.

    Geometric methods work off these centroids rather than off lattice indices
    both because a `.ume` file has no lattice to speak of and because sorting
    centroids reproduces the lattice cuts exactly on a mesh that does. */
DS_Types::VEC3V_T const &zone_centroids(SOA_Idx::Mesh const &mesh);

} // namespace Part
} // namespace Ume

#endif
