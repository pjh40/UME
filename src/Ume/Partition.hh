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
  \file Ume/Partition.hh

  The partitioning strategy interface.

  A `Partitioner` answers exactly one question: which part does each real zone
  belong to?  Everything downstream of that -- ownership of the points, edges
  and faces on a cut, the ghost layer, the local renumbering, the neighbor
  lists -- is a pure function of the resulting `Zone_Map` and lives in
  `Ume/Decompose.hh`.  Keeping the contract this narrow is what makes the
  methods interchangeable, and lets a zone map be dumped to a file and read
  back so that the strategy and extraction layers can be tested apart.

  Not every method is available in every build: the library-backed ones
  disappear when the library was not found at configure time.  `all_methods()`
  therefore describes every method the source knows about, whether or not this
  build can run it, so that a misspelled method name and an unavailable one can
  be reported differently.
*/

#ifndef UME_PARTITION_HH
#define UME_PARTITION_HH 1

#include "Ume/DS_Types.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace Ume {

//! Mesh partitioning strategies
namespace Part {

//! Zone index -> part number, with length `mesh.zones.size()`
/*! Zones that are not real (`mask < 1`) hold -1: a boundary zone follows the
    real zone it hangs off, which is a decision for the decomposition rather
    than for a strategy. */
using Zone_Map = std::vector<int>;

//! What a partitioner should charge for owning a zone
enum class Vertex_Weight {
  UNIT, //!< every zone costs the same
  CORNERS, //!< the length of the zone's `m:z>c` row
  SIDES //!< the number of real sides of the zone
};

//! What a partitioner should charge for cutting a face
enum class Edge_Weight {
  UNIT, //!< every face costs the same
  FACE_POINTS, //!< the number of distinct points on the face
  FACE_AREA //!< the area of the face
};

//! Everything a `Partitioner` is allowed to look at besides the mesh
struct Options {
  //! Number of parts to produce; parts are numbered `[0, nparts)`
  int nparts{1};
  //! Seed for any method that makes random choices
  std::uint64_t seed{0x9E3779B97F4A7C15ULL};
  Vertex_Weight vwgt{Vertex_Weight::CORNERS};
  Edge_Weight ewgt{Edge_Weight::FACE_POINTS};
  //! Allowed ratio of the heaviest part to the average part
  double imbalance{1.03};
  //! Method-specific settings, passed through from `--set key=value`
  std::map<std::string, std::string> extra;
};

//! A strategy for assigning zones to parts
class Partitioner {
public:
  virtual ~Partitioner() = default;
  //! The name this method is selected by
  virtual char const *name() const = 0;
  //! One line describing what the method does
  virtual char const *describe() const = 0;

  //! Assign every zone with `mask >= 1` to a part in `[0, opt.nparts)`
  /*! `z2r` is resized to `mesh.zones.size()` and every entry written: -1 for
      zones that are not real, a part number for the rest.  Throws
      `std::invalid_argument` if `opt.nparts < 1`. */
  virtual void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const = 0;
};

//! One entry in the compile-time table of partitioning methods
struct Method {
  //! The name `make()` selects on
  char const *name;
  //! One line describing what the method does
  char const *describe;
  //! The external package the method needs, or nullptr if it needs none
  char const *requires_pkg;
  //! False when `requires_pkg` was not found when this build was configured
  bool available;
};

//! Every method this source knows about, available in this build or not
std::span<Method const> all_methods();

//! The names of the methods this build can actually run
std::vector<std::string> methods();

//! Look `name` up in the method table; nullptr if there is no such method
Method const *find_method(std::string_view name);

//! Construct a partitioner; nullptr if `name` is unknown *or* unavailable
/*! Use `find_method()` to tell those two cases apart. */
std::unique_ptr<Partitioner> make(std::string_view name);

/*! @name Building blocks shared by the built-in methods */
///@{

//! The indices of the real (`mask >= 1`) zones of `mesh`, ascending
std::vector<int> real_zones(SOA_Idx::Mesh const &mesh);

//! The weight `opt.vwgt` assigns to each zone, 0 for zones that are not real
/*! Indexed by zone index, with length `mesh.zones.size()`. */
std::vector<int> zone_weights(SOA_Idx::Mesh const &mesh, Options const &opt);

//! The extent of the centroids of the zones `v` along each axis
/*! An axis along which every zone of `v` has the same centroid coordinate --
    the third axis of a single-cell-thick mesh, say -- has extent 0, and an
    empty `v` has extent 0 everywhere. */
std::array<double, 3> centroid_extent(
    std::span<int const> v, DS_Types::VEC3V_T const &zcoord);

//! Shape `nparts` into a grid of per-axis division counts
/*! The product of the result is exactly `nparts`.  Each prime factor of
    `nparts`, largest first, is given to the axis that is currently the least
    divided *relative to its extent*, so the boxes come out as cubical as the
    factorization allows -- in physical rather than in index space, which is
    what makes this work on a mesh with no lattice behind it and on one whose
    cells are not cubes.

    Only the `naxes` longest axes are divided; the rest keep a single
    division, which is the whole difference between `block` (3) and `pencil`
    (2).  An axis of zero extent is never divided unless every axis has zero
    extent.  Throws `std::invalid_argument` if `nparts < 1`. */
std::array<int, 3> dims_create(
    int nparts, std::array<double, 3> const &extent, int naxes = 3);

//! Split the zones `v` into `k` groups of equal weight, ordered by `key`
/*! `key` and `wgt` are indexed by *zone index*, not by position in `v`, so a
    recursive method can split a subset without recompacting anything; `wgt`
    may be empty, which means unit weights.  Group `i` is written to `z2r` as
    part `first_part + i`.

    The split is made at cumulative-*weight* boundaries rather than at uniform
    spatial or count boundaries, which is what gives exact load balance on a
    mesh whose zones are not all the same cost.  When `v` holds at least `k`
    zones every group gets at least one, even if the weights are wildly uneven.
    Ties in `key` are broken by zone index, so the result does not depend on
    the sort algorithm. */
void split_sorted(std::span<int const> v, std::span<double const> key,
    std::span<int const> wgt, int first_part, int k, Zone_Map &z2r);

///@}

} // namespace Part
} // namespace Ume

#endif
