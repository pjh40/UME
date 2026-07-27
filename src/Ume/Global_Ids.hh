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
  \file Ume/Global_Ids.hh

  The `<base>.NNNNN.gid` side-car: what each local entity of a partition was
  called in the mesh it came from.

  A `.ume` file says nothing about where its entities sat in the serial mesh --
  the format has no place to put it, and adding one would change what every
  existing file means.  But comparing a partitioned run against a serial one
  needs exactly that, on both ends: the two runs have to perturb the *same*
  zone, and the answers have to be lined up zone for zone.  So the mapping goes
  beside the mesh in its own file, written only when asked for
  (`ume_partition --global-ids`), and nothing in the library depends on it.

  The format is text, one section per entity kind in the order
  `Decomposition::layouts()` lists them:

  ```
  ume_gid 1
  mype 2
  numpe 3
  point 24 0
  0 1 2 3 4 5 6 7 8 9
  10 11 ...
  edge 44 0
  ...
  ```

  Each section header is the kind name, the number of real entities, and the
  number of ghosts; the indices that follow are in local index order, so the
  first `local` of them are the real ones.  Being text rather than binary is
  deliberate: a side-car exists to be looked at when a comparison disagrees.
*/

#ifndef UME_GLOBAL_IDS_HH
#define UME_GLOBAL_IDS_HH 1

#include "Ume/Decompose.hh"
#include <array>
#include <iosfwd>
#include <string>
#include <string_view>

namespace Ume {
namespace Part {

//! One kind's local-index-to-serial-index map
struct Kind_Ids {
  //! `gid[l]` for `l < local` is a real entity; the rest are ghosts
  int local{0};
  std::vector<int> gid;

  int ghost() const { return static_cast<int>(gid.size()) - local; }
};

//! One rank's side-car
struct Global_Ids {
  int mype{-1}, numpe{0};
  //! Indexed by `Kind`, which is also the order they appear in the file
  std::array<Kind_Ids, NKINDS> kinds;

  Kind_Ids const &operator[](int k) const {
    return kinds[static_cast<std::size_t>(k)];
  }
  Kind_Ids &operator[](int k) { return kinds[static_cast<std::size_t>(k)]; }
};

//! Write rank `pe`'s side-car
/*! Throws `std::out_of_range` if `pe` is not a part of `d`. */
void write_global_ids(Decomposition const &d, int pe, std::ostream &os);

//! Read a side-car back
/*! Throws `std::runtime_error` if the stream does not hold one: a partition
    read together with the wrong side-car would silently compare the wrong
    entities, which is the one outcome worth being loud about. */
Global_Ids read_global_ids(std::istream &is);

//! `<base>.NNNNN.<ext>`, the per-rank naming `txt2bin` established
/*! A partition and its side-car have to be found by the same rule, so the rule
    lives in one place. */
std::string rank_file(std::string_view base, int pe, std::string_view ext);

} // namespace Part
} // namespace Ume

#endif
