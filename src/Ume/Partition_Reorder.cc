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
  \file Ume/Partition_Reorder.cc
*/

#include "Ume/Partition_Reorder.hh"
#include <algorithm>
#include <charconv>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Ume {
namespace Part {

int Machine::distance(int a, int b) const {
  long span = 1;
  for (std::size_t i = 0; i < levels.size(); ++i) {
    span *= levels[i];
    if (a / span == b / span)
      return static_cast<int>(i);
  }
  return depth();
}

Machine parse_machine(std::string_view spec, int nranks) {
  if (nranks < 1)
    throw std::invalid_argument("a machine needs at least one rank");

  Machine mach;
  mach.nranks = nranks;
  if (spec.empty())
    throw std::invalid_argument("empty machine topology");

  for (std::size_t at = 0; at <= spec.size();) {
    std::size_t const sep = std::min(spec.find('x', at), spec.size());
    std::string_view const field = spec.substr(at, sep - at);
    int size = 0;
    auto const [ptr, ec] =
        std::from_chars(field.data(), field.data() + field.size(), size);
    if (ec != std::errc{} || ptr != field.data() + field.size() || size < 1) {
      std::ostringstream os;
      os << "bad machine topology \"" << spec
         << "\": expected x-separated group sizes, each at least 1";
      throw std::invalid_argument(os.str());
    }
    mach.levels.push_back(size);
    at = sep + 1;
  }
  return mach;
}

namespace {

//! Elements crossing between `a` and `b` in either direction
std::int64_t traffic(Metrics const &m, int a, int b) {
  return static_cast<std::int64_t>(m.sent(a, b)) +
      static_cast<std::int64_t>(m.sent(b, a));
}

/*! What swapping the ranks of parts `a` and `b` would add to the cost.  The
    `a`-`b` term itself cancels: the two are the same distance apart either
    way round. */
std::int64_t swap_delta(Metrics const &m, Machine const &mach,
    std::span<int const> order, int a, int b) {
  int const sa = order[static_cast<std::size_t>(a)];
  int const sb = order[static_cast<std::size_t>(b)];
  std::int64_t delta = 0;
  for (int c = 0; c < m.nparts; ++c) {
    if (c == a || c == b)
      continue;
    int const sc = order[static_cast<std::size_t>(c)];
    std::int64_t const shift = mach.distance(sb, sc) - mach.distance(sa, sc);
    delta += (traffic(m, a, c) - traffic(m, b, c)) * shift;
  }
  return delta;
}

//! Swap pairs of parts until no single swap improves the cost
/*! The cost is a non-negative integer and every accepted swap strictly
    lowers it, so this terminates on its own and needs no sweep count.  The
    lesson `Partition_BFS.cc` learned about bounded patience applies here
    too: a cap would be a quality knob set by the size of the last machine
    anybody tested on. */
void refine(Metrics const &m, Machine const &mach, std::vector<int> &order) {
  for (bool moved = true; moved;) {
    moved = false;
    for (int a = 0; a < m.nparts; ++a)
      for (int b = a + 1; b < m.nparts; ++b)
        if (swap_delta(m, mach, order, a, b) < 0) {
          std::swap(order[static_cast<std::size_t>(a)],
              order[static_cast<std::size_t>(b)]);
          moved = true;
        }
  }
}

//! Place the heaviest part first, then whatever is most tied to what is down
std::vector<int> greedy(Metrics const &m, Machine const &mach) {
  auto const nps = static_cast<std::size_t>(m.nparts);
  std::vector<int> order(nps, -1);
  std::vector<char> taken(static_cast<std::size_t>(mach.nranks), 0);

  /* Affinity to the parts already placed, which is what decides who goes
     next; before anything is placed it is the part's whole traffic, so the
     busiest part is seeded first and the rest gather around it. */
  std::vector<std::int64_t> affinity(nps, 0);
  for (int a = 0; a < m.nparts; ++a)
    for (int b = 0; b < m.nparts; ++b)
      if (a != b)
        affinity[static_cast<std::size_t>(a)] += traffic(m, a, b);

  for (int placed = 0; placed < m.nparts; ++placed) {
    int next = -1;
    for (int p = 0; p < m.nparts; ++p)
      if (order[static_cast<std::size_t>(p)] < 0 &&
          (next < 0 ||
              affinity[static_cast<std::size_t>(p)] >
                  affinity[static_cast<std::size_t>(next)]))
        next = p;

    int best_rank = -1;
    std::int64_t best_cost = 0;
    for (int r = 0; r < mach.nranks; ++r) {
      if (taken[static_cast<std::size_t>(r)])
        continue;
      std::int64_t cost = 0;
      for (int q = 0; q < m.nparts; ++q) {
        int const sq = order[static_cast<std::size_t>(q)];
        if (sq >= 0)
          cost += traffic(m, next, q) * mach.distance(r, sq);
      }
      if (best_rank < 0 || cost < best_cost) {
        best_rank = r;
        best_cost = cost;
      }
    }

    order[static_cast<std::size_t>(next)] = best_rank;
    taken[static_cast<std::size_t>(best_rank)] = 1;
    /* `next` is placed, so it now counts towards everyone else's affinity
       and stops counting towards its own -- which is why the loop above
       tests `order` rather than trusting the affinity to have gone stale. */
    affinity[static_cast<std::size_t>(next)] = -1;
    for (int p = 0; p < m.nparts; ++p)
      if (order[static_cast<std::size_t>(p)] < 0)
        affinity[static_cast<std::size_t>(p)] += traffic(m, p, next);
  }
  return order;
}

} // namespace

std::int64_t hop_cost(
    Metrics const &m, Machine const &mach, std::span<int const> order) {
  if (static_cast<int>(order.size()) != m.nparts)
    throw std::invalid_argument("a rank order covers every part exactly once");

  std::int64_t cost = 0;
  for (int a = 0; a < m.nparts; ++a)
    for (int b = 0; b < m.nparts; ++b)
      if (a != b)
        cost += static_cast<std::int64_t>(m.sent(a, b)) *
            mach.distance(order[static_cast<std::size_t>(a)],
                order[static_cast<std::size_t>(b)]);
  return cost;
}

namespace {

std::vector<int> identity_order(int nparts) {
  std::vector<int> v(static_cast<std::size_t>(nparts));
  std::iota(v.begin(), v.end(), 0);
  return v;
}

} // namespace

std::vector<int> reorder_ranks(Metrics const &m, Machine const &mach) {
  if (mach.nranks < m.nparts)
    throw std::invalid_argument("the machine has fewer ranks than parts");

  std::vector<int> identity = identity_order(m.nparts);
  if (m.nparts < 2 || mach.depth() == 0)
    return identity;

  std::vector<int> built = greedy(m, mach);
  refine(m, mach, built);
  refine(m, mach, identity);

  /* Ties go to the identity: relabelling the ranks rewrites every file for
     nothing unless it buys something, and a machine that cannot tell its
     ranks apart should leave the numbering the method chose. */
  if (hop_cost(m, mach, built) < hop_cost(m, mach, identity))
    return built;
  return identity;
}

Rank_Order place_ranks(SOA_Idx::Mesh const &src, Zone_Map const &z2r,
    int nparts, Decompose_Options const &dopt, Machine const &mach) {
  Rank_Order out;
  out.order = identity_order(nparts);

  Metrics const before =
      measure(src, build_decomposition(src, z2r, nparts, dopt));
  out.before = hop_cost(before, mach, out.order);
  out.after = out.before;

  std::vector<int> const chosen = reorder_ranks(before, mach);
  if (chosen == out.order)
    return out;

  /* What the chosen order predicts and what it produces are two numbers, and
     this is the one the machine will actually pay: renaming the ranks moves
     the owner of every shared entity, so the matrix it was chosen from is
     not quite the matrix it leaves behind. */
  Zone_Map moved = z2r;
  apply_rank_order(chosen, moved);
  Metrics const after =
      measure(src, build_decomposition(src, std::move(moved), nparts, dopt));
  std::int64_t const got = hop_cost(after, mach, out.order);
  if (got < out.before) {
    out.order = chosen;
    out.after = got;
  }
  return out;
}

void apply_rank_order(std::span<int const> order, Zone_Map &z2r) {
  for (int &part : z2r) {
    if (part < 0)
      continue;
    if (static_cast<std::size_t>(part) >= order.size())
      throw std::invalid_argument("a zone names a part the rank order lacks");
    part = order[static_cast<std::size_t>(part)];
  }
}

} // namespace Part
} // namespace Ume
