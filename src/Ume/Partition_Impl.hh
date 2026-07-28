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
  \file Ume/Partition_Impl.hh

  The seam between `Partition.hh` and the files that implement the methods.

  `Partition.hh` is what the rest of UME sees: a name goes in and a
  `Partitioner` comes out.  The methods themselves are spread across several
  translation units -- one per family -- and this header is how `make()`
  reaches classes that are otherwise file-local.  Splitting them up is not
  tidiness: a family backed by an external library has to be compilable out
  along with the library, and a `#ifndef NO_METIS` around a whole file is a
  great deal easier to be sure of than one around a region of a large one.

  Adding a family is three things: a `make_*` declaration here, its definition
  in the new file, and a line in `make()`.  Each factory returns nullptr for a
  name it does not implement, so `make()` is a walk down the list rather than
  a table that has to agree with itself.

  This header is deliberately not installed: nothing outside `libume` should
  be constructing a `Partitioner` except through `make()`.
*/

#ifndef UME_PARTITION_IMPL_HH
#define UME_PARTITION_IMPL_HH 1

#include "Ume/Partition.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include <memory>
#include <vector>

namespace Ume {
namespace Part {

//! Prepare `z2r` for `assign()`, and return the zones a method may place
/*! Every `Partitioner::assign()` opens with this: it rejects a nonsensical
    part count, fills `z2r` with -1 so that a method only has to write the
    zones it actually places, and hands back the real zones in ascending index
    order.  Throws `std::invalid_argument` if `opt.nparts < 1`. */
std::vector<int> begin_assign(
    SOA_Idx::Mesh const &mesh, Options const &opt, Zone_Map &z2r);

/*! @name The method factories, one per implementation file

    `m` is the table entry `find_method()` matched, already checked for
    availability; the geometric and grid methods take their `name()` and
    `describe()` straight from it, so one class serves several table entries.
    A factory handed a name it does not implement returns nullptr. */
///@{

//! `slab`, `slab-x/y/z`, `block`, `pencil`, `rcb`, `rib`, `morton`, `hilbert`
std::unique_ptr<Partitioner> make_geometric(Method const &m);

//! `bfs`
std::unique_ptr<Partitioner> make_bfs(Method const &m);

//! `metis-kway-cut`, `metis-kway-vol`, `metis-rb`; nothing without METIS
std::unique_ptr<Partitioner> make_metis(Method const &m);

///@}

} // namespace Part
} // namespace Ume

#endif
