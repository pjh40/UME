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

  The mesh views that the partitioning methods share: zone centroids, the dual
  graph and the hypergraph.
*/

#ifndef UME_PARTITION_GRAPH_HH
#define UME_PARTITION_GRAPH_HH 1

#include "Ume/DS_Types.hh"
#include "Ume/Partition.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include <cstddef>
#include <span>
#include <vector>

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

//! The dual graph of the mesh: a vertex per real zone, an edge per shared face
/*! Held in the compressed-row form every graph partitioner in sight wants, so
    that a library method can hand `xadj`/`adjncy`/`vwgt`/`adjwgt` straight
    over rather than translating a friendlier structure at the boundary.

    Only real (`mask >= 1`) zones are vertices, so the vertex numbering is a
    compaction of the zone numbering rather than the zone numbering itself;
    `z2v` and `v2z` carry the correspondence both ways.  Two zones that share
    more than one face -- which a general polyhedral mesh allows even though a
    brick does not -- give one edge carrying the sum of the faces' weights,
    because it is the pair that has to be separated, not each face. */
struct Zone_Graph {
  //! The number of real zones
  int nverts{0};
  //! Row starts into `adjncy` and `adjwgt`, with length `nverts + 1`
  std::vector<int> xadj;
  //! The neighbors of each vertex, ascending within a row
  std::vector<int> adjncy;
  //! What cutting each edge of `adjncy` costs, per `Options::ewgt`
  std::vector<int> adjwgt;
  //! What owning each vertex costs, per `Options::vwgt`
  std::vector<int> vwgt;
  //! Zone index -> vertex, -1 for a zone that is not real
  std::vector<int> z2v;
  //! Vertex -> zone index, ascending
  std::vector<int> v2z;

  //! The neighbors of vertex `v`
  std::span<int const> neighbors(int v) const;
  //! The edge weights of vertex `v`, positionally matching `neighbors(v)`
  std::span<int const> weights(int v) const;
};

//! Build the dual graph of `mesh`, weighted as `opt` asks
/*! `Edge_Weight::FACE_AREA` needs `side_surz`, which is derived, so `mesh`
    must have a Transport attached as for `zone_centroids`.  Areas are scaled
    onto integers with the largest face weighing `AREA_SCALE`: the weights are
    only ever compared with one another, so the scale is arbitrary, and fixing
    the top keeps three significant figures of every ratio while leaving the
    sums nowhere near overflowing. */
Zone_Graph build_zone_graph(SOA_Idx::Mesh const &mesh, Options const &opt);

//! What the largest face weighs under `Edge_Weight::FACE_AREA`
inline constexpr int AREA_SCALE = 1000;

//! The mesh as a hypergraph: a vertex per real zone, a net per shared point
/*! Where the dual graph asks which *faces* a decomposition cuts, this asks
    which *points* it duplicates -- and that is the question UME's
    communication actually poses.  A point held by \f$\lambda\f$ ranks needs
    \f$\lambda-1\f$ copies, each of which is one element in one direction of
    the `points.gathscat()` that every gradient kernel ends with, so
    \f$\sum(\lambda-1)\f$ over the nets is not a proxy for the traffic: it is
    the traffic, and `Metrics::kinds[POINT].duplicated` is the same number
    counted from the other end.  A graph edge cut cannot express it, because a
    point on three cut faces of one rank pair is charged three times by the
    faces and once by the point.

    The vertices are the real zones, compacted exactly as `Zone_Graph`
    compacts them, so `z2v`/`v2z` mean the same thing in both and a zone map
    built from one can be compared against the other.  The nets are the points
    of `m:p>zs` that touch two or more real zones; a point touching one is a
    net that can never be cut, and leaving it out changes no objective. */
struct Zone_Hypergraph {
  //! The number of real zones
  int nverts{0};
  //! The number of points that touch more than one real zone
  int nnets{0};
  /*! Row starts into `net_pins`, with length `nnets + 1`.  `std::size_t`
      rather than `int` because that is the width Mt-KaHyPar's
      `hyperedge_indices` is, and it is the one array that can be handed over
      without a copy. */
  std::vector<std::size_t> net_index;
  //! The vertices on each net, ascending within a net
  std::vector<int> net_pins;
  /*! What each net costs to cut, with length `nnets`.  Every point carries
      the same payload, so this is 1 everywhere -- which is precisely what
      makes the connectivity objective over this hypergraph equal to the
      element count on the wire rather than proportional to it.
      `Options::ewgt` weighs *faces* and has no bearing here. */
  std::vector<int> netwgt;
  //! What owning each vertex costs, per `Options::vwgt`, at least 1
  std::vector<int> vwgt;
  //! Zone index -> vertex, -1 for a zone that is not real
  std::vector<int> z2v;
  //! Vertex -> zone index, ascending
  std::vector<int> v2z;
  //! Net -> the point index it came from, ascending
  std::vector<int> net2p;

  //! The vertices on net `e`
  std::span<int const> pins(int e) const;
};

//! Build the hypergraph of `mesh`, with vertices weighted as `opt` asks
/*! Reads the derived `m:p>zs`, so `mesh` must have a Transport attached as
    for `zone_centroids`. */
Zone_Hypergraph build_zone_hypergraph(
    SOA_Idx::Mesh const &mesh, Options const &opt);

} // namespace Part
} // namespace Ume

#endif
