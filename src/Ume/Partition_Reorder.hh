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
  \file Ume/Partition_Reorder.hh

  Which rank each part is given, once the parts are decided.

  A partitioning method answers *what goes together*; it says nothing about
  where the resulting parts sit on a machine, and the launcher numbers ranks
  `0..N-1` in whatever order it likes.  On a machine where a hop between two
  nodes costs more than a hop inside one -- which is every machine UME runs on
  -- those two questions have different answers, and this file holds the
  second one: a permutation `part -> rank` that puts heavy neighbours close
  together.

  ### It is a relabelling, but not quite only a relabelling

  `apply_rank_order` renames the parts in a `Zone_Map`, so the same zones stay
  together, the same points are duplicated, and the same number of elements
  goes on the wire.  Everything the methods compete over is therefore
  invariant under it, and this cannot make a decomposition worse -- or
  better -- in any of those terms.

  What is *not* invariant is which rank sends and which receives.
  `Owner_Rule::LOWEST_RANK` names the SOURCE of a shared entity by rank
  number, so an entity held by three ranks is a star centred on the
  lowest-numbered of them, and renaming moves the centre: the comm matrix
  changes shape, and with it the neighbour counts, while its total does not.

  That is why `place_ranks` rebuilds and re-measures rather than trusting
  `hop_cost` of the order it chose.  The prediction is close -- within a
  fraction of a percent on everything measured here -- but "never worse than
  the numbering you had" is worth having exactly, and it costs one more
  `measure()` to have it as a fact rather than as a model.

  ### The machine model

  A `Machine` is a tree of equal-sized groups over consecutively numbered
  ranks: `4x2` is four ranks per node and two nodes per switch.  That is
  coarse, deliberately.  A real network is not a balanced tree, its distances
  are not integers, and a site's rank-to-node mapping is a scheduler policy
  rather than a property of the hardware -- so a more detailed model would be
  more precisely wrong.  What the tree does capture is the one distinction
  that dominates: traffic that stays inside a node against traffic that does
  not.
*/

#ifndef UME_PARTITION_REORDER_HH
#define UME_PARTITION_REORDER_HH 1

#include "Ume/Decompose.hh"
#include "Ume/Partition.hh"
#include "Ume/Partition_Metrics.hh"
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace Ume {
namespace Part {

//! A machine as a tree of equal-sized groups over consecutive rank numbers
struct Machine {
  /*! Group sizes from the leaves up: `{ranks per node, nodes per switch}`.
      An empty list is a flat machine, on which every pair of ranks is equally
      far apart -- so every order costs the same and there is nothing to
      choose. */
  std::vector<int> levels;
  //! How many ranks the machine has room for; at least the part count
  int nranks{1};

  /*! Hops between two ranks: 0 when they share a leaf group, 1 when they
      share the group above it, and so on up to `depth()` for two ranks with
      nothing in common.  Symmetric, and 0 on the diagonal. */
  int distance(int a, int b) const;

  //! The most any pair of ranks can be apart
  int depth() const { return static_cast<int>(levels.size()); }
};

//! Read a machine from an `x`-separated list of group sizes, leaves first
/*! `"4"` is four ranks per node, `"4x2"` four per node and two nodes per
    switch.  A trailing level covering everything is implied, so a spec need
    not multiply out to `nranks` and `"64"` on a 16-rank job is one node.
    Throws `std::invalid_argument` on an empty spec, a malformed one, or a
    group size below 1. */
Machine parse_machine(std::string_view spec, int nranks);

//! What one exchange sweep costs on `mach`, in element-hops
/*! \f$\sum_{a \neq b} \mathrm{sent}(a,b) \cdot d(\mathrm{order}[a],
    \mathrm{order}[b])\f$ -- the elements on the wire weighted by how far each
    one travels.  `order` maps part to rank and must be a permutation of
    `m.nparts` distinct ranks in `[0, mach.nranks)`.

    Zero on a flat machine, whatever the order, because nothing there is
    farther away than anything else. */
std::int64_t hop_cost(
    Metrics const &m, Machine const &mach, std::span<int const> order);

//! Choose a rank order, and measure what it really achieved
/*! `reorder_ranks` chooses from the comm matrix of the decomposition as it is
    numbered now, and applying its answer moves the SOURCE of every shared
    entity -- see the note at the top of this file.  So the chosen order is
    applied, the decomposition rebuilt and measured, and the identity kept
    instead unless the rebuilt cost really is lower.  `before` and `after` are
    both measured numbers, and `after <= before` always.

    `z2r` must be a map `build_decomposition` accepts. */
struct Rank_Order {
  //! Part -> rank, and the identity when nothing was gained
  std::vector<int> order;
  //! Element-hops per exchange sweep, as measured before and after
  std::int64_t before{0}, after{0};
};
Rank_Order place_ranks(SOA_Idx::Mesh const &src, Zone_Map const &z2r,
    int nparts, Decompose_Options const &dopt, Machine const &mach);

//! Choose which rank each part should be given
/*! Greedy placement -- heaviest part first, then whichever unplaced part is
    most strongly tied to what is already down, into the free rank that costs
    least -- followed by pairwise swaps until no swap improves.  The identity
    order is refined the same way, and the cheaper of the two is returned, so
    the result is never worse than the order the method produced; ties go to
    the identity, which is what makes a flat machine leave the numbering
    alone.

    The refinement terminates because the cost is a non-negative integer that
    strictly decreases at every swap, so no iteration cap appears here.  It
    costs \f$O(\mathrm{nparts}^3)\f$, over the rank count rather than over the
    mesh -- a thousand ranks is well under a second, and a thousand ranks is
    a large job.

    This is a quadratic assignment problem, so nothing here is optimal; what
    is promised is only `hop_cost(result) <= hop_cost(identity)`. */
std::vector<int> reorder_ranks(Metrics const &m, Machine const &mach);

//! Rename the parts of `z2r` through `order`, leaving the -1 entries alone
/*! Throws `std::invalid_argument` if `z2r` names a part `order` has no entry
    for. */
void apply_rank_order(std::span<int const> order, Zone_Map &z2r);

} // namespace Part
} // namespace Ume

#endif
