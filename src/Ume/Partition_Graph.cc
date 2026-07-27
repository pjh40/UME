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
  \file Ume/Partition_Graph.cc
*/

#include "Ume/Partition_Graph.hh"

namespace Ume {
namespace Part {

DS_Types::VEC3V_T const &zone_centroids(SOA_Idx::Mesh const &mesh) {
  /* `caccess_*` is const but forces the lazy `init_()`, which is why a
     Partitioner can take the mesh by const reference and still ask for a
     derived variable. */
  return mesh.ds->caccess_vec3v("zcoord");
}

} // namespace Part
} // namespace Ume
