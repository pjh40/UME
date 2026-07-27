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
  \file Ume/Decompose.cc
*/

#include "Ume/Decompose.hh"
#include <algorithm>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>

namespace Ume {
namespace Part {

int Entity_Layout::local_size(int r) const {
  return start[static_cast<std::size_t>(r) + 1] -
      start[static_cast<std::size_t>(r)];
}

int Entity_Layout::ghost_local_size(int r) const {
  return gstart[static_cast<std::size_t>(r) + 1] -
      gstart[static_cast<std::size_t>(r)];
}

int Entity_Layout::local_index(int r, int g) const {
  auto const lo = gid.begin() + start[static_cast<std::size_t>(r)];
  auto const hi = gid.begin() + start[static_cast<std::size_t>(r) + 1];
  auto const it = std::lower_bound(lo, hi, g);
  if (it != hi && *it == g)
    return static_cast<int>(it - lo);

  auto const glo = ggid.begin() + gstart[static_cast<std::size_t>(r)];
  auto const ghi = ggid.begin() + gstart[static_cast<std::size_t>(r) + 1];
  auto const git = std::lower_bound(glo, ghi, g);
  if (git != ghi && *git == g)
    return local_size(r) + static_cast<int>(git - glo);
  return -1;
}

int Entity_Layout::global_index(int r, int l) const {
  int const nl = local_size(r);
  if (l < nl)
    return gid[static_cast<std::size_t>(start[static_cast<std::size_t>(r)] + l)];
  return ggid[static_cast<std::size_t>(
      gstart[static_cast<std::size_t>(r)] + l - nl)];
}

namespace {

//! (global index, rank), so that a sort groups the ranks holding one entity
using Pairs = std::vector<std::pair<int, int>>;

void sort_unique(Pairs &pairs) {
  std::ranges::sort(pairs);
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
}

//! Lay sorted, deduplicated `pairs` out as a rank-major CSR
/*! Because the input is keyed by global index, a stable pass over it leaves
    each rank's list in ascending global order -- which is the local index
    order the whole design rests on. */
void to_csr(Pairs const &pairs, int nparts, std::vector<int> &start,
    std::vector<int> &gid) {
  start.assign(static_cast<std::size_t>(nparts) + 1, 0);
  for (auto const &[g, r] : pairs)
    start[static_cast<std::size_t>(r) + 1] += 1;
  for (int r = 0; r < nparts; ++r)
    start[static_cast<std::size_t>(r) + 1] += start[static_cast<std::size_t>(r)];

  gid.resize(pairs.size());
  std::vector<int> fill(start.begin(), start.end() - 1);
  for (auto const &[g, r] : pairs)
    gid[static_cast<std::size_t>(fill[static_cast<std::size_t>(r)]++)] = g;
}

std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

//! Pick the SOURCE rank of every global entity from the ranks holding it
/*! `pairs` must be sorted and deduplicated, so the ranks holding one entity
    form an ascending run.  `LEAST_LOADED` counts across the whole entity kind,
    which is the only scale at which balancing SOURCE duty means anything. */
void assign_owners(Pairs const &pairs, int nglobal, int nparts, Owner_Rule rule,
    std::vector<int> &owner) {
  owner.assign(static_cast<std::size_t>(nglobal), -1);
  std::vector<long> load(static_cast<std::size_t>(nparts), 0);

  std::size_t i = 0;
  while (i < pairs.size()) {
    std::size_t j = i;
    while (j < pairs.size() && pairs[j].first == pairs[i].first)
      ++j;
    int const g = pairs[i].first;
    std::span<std::pair<int, int> const> const candidates{
        pairs.data() + i, j - i};

    int chosen = candidates.front().second;
    switch (rule) {
    case Owner_Rule::LOWEST_RANK:
      break;
    case Owner_Rule::HASH:
      chosen = candidates[splitmix64(static_cast<std::uint64_t>(g)) %
          candidates.size()]
                   .second;
      break;
    case Owner_Rule::LEAST_LOADED:
      for (auto const &[gg, r] : candidates)
        if (load[static_cast<std::size_t>(r)] <
            load[static_cast<std::size_t>(chosen)])
          chosen = r;
      break;
    }
    owner[static_cast<std::size_t>(g)] = chosen;
    load[static_cast<std::size_t>(chosen)] += 1;
    i = j;
  }
}

//! The real elements of each zone, from an element-to-zone map
/*! The derived `m:z>c` would do for corners, but there is no `m:z>s` at all,
    and a counting sort over the stored map is both cheaper than a derived
    variable and says exactly what is meant here. */
class Zone_Elements {
public:
  Zone_Elements(char const *what, DS_Types::INTV_T const &e2z,
      std::vector<short> const &emask, int nzones) {
    int const n = static_cast<int>(emask.size());
    auto zone_of = [&](int e) {
      int const z = e2z[static_cast<std::size_t>(e)];
      if (z < 0 || z >= nzones)
        throw std::invalid_argument(std::string("build_decomposition: ") +
            what + ' ' + std::to_string(e) + " names zone " +
            std::to_string(z) + ", outside [0, " + std::to_string(nzones) +
            ")");
      return z;
    };

    start_.assign(static_cast<std::size_t>(nzones) + 1, 0);
    for (int e = 0; e < n; ++e)
      if (emask[static_cast<std::size_t>(e)] >= 1)
        start_[static_cast<std::size_t>(zone_of(e)) + 1] += 1;
    for (int z = 0; z < nzones; ++z)
      start_[static_cast<std::size_t>(z) + 1] +=
          start_[static_cast<std::size_t>(z)];

    elem_.resize(static_cast<std::size_t>(start_.back()));
    std::vector<int> fill(start_.begin(), start_.end() - 1);
    for (int e = 0; e < n; ++e)
      if (emask[static_cast<std::size_t>(e)] >= 1)
        elem_[static_cast<std::size_t>(
            fill[static_cast<std::size_t>(zone_of(e))]++)] = e;
  }

  std::span<int const> operator[](int z) const {
    auto const b = static_cast<std::size_t>(start_[static_cast<std::size_t>(z)]);
    auto const e =
        static_cast<std::size_t>(start_[static_cast<std::size_t>(z) + 1]);
    return std::span<int const>{elem_.data() + b, e - b};
  }

private:
  std::vector<int> start_, elem_;
};

//! True if any entity of `src` is a ghost
bool has_ghosts(SOA_Idx::Mesh const &src) {
  SOA_Idx::Entity const *const entities[] = {&src.points, &src.edges,
      &src.faces, &src.zones, &src.corners, &src.sides, &src.iotas};
  return std::ranges::any_of(entities, [](SOA_Idx::Entity const *e) {
    return e->local_size() != e->size();
  });
}

//! Give every zone that is not a ghost a rank
/*! Real zones take theirs from `z2r`.  A boundary zone has one face, and the
    real zone behind it is what its phantom sides pair with -- which works for
    either `Brick::Boundary_Style`, since with phantom sides only there are no
    boundary zones to place. */
std::vector<int> zone_ranks(SOA_Idx::Mesh const &src, Zone_Map const &z2r) {
  int const zll = src.zones.size(), sll = src.sides.size();
  auto const &zmask = src.zones.mask;
  auto const &smask = src.sides.mask;
  auto const &s2z = src.ds->caccess_intv("m:s>z");
  auto const &s2s2 = src.ds->caccess_intv("m:s>s2");

  std::vector<int> zrank(static_cast<std::size_t>(zll), -1);
  for (int z = 0; z < zll; ++z)
    if (zmask[z] >= 1)
      zrank[static_cast<std::size_t>(z)] = z2r[static_cast<std::size_t>(z)];

  for (int s = 0; s < sll; ++s) {
    if (smask[s] >= 0)
      continue;
    int const z = s2z[static_cast<std::size_t>(s)];
    if (z < 0 || z >= zll || zmask[z] >= 0 ||
        zrank[static_cast<std::size_t>(z)] >= 0)
      continue;
    int const partner = s2s2[static_cast<std::size_t>(s)];
    if (partner < 0 || partner >= sll)
      continue;
    int const real = s2z[static_cast<std::size_t>(partner)];
    if (real >= 0 && real < zll && zmask[real] >= 1)
      zrank[static_cast<std::size_t>(z)] = z2r[static_cast<std::size_t>(real)];
  }

  for (int z = 0; z < zll; ++z)
    if (zmask[z] != 0 && zrank[static_cast<std::size_t>(z)] < 0)
      throw std::invalid_argument("build_decomposition: zone " +
          std::to_string(z) + " (mask " + std::to_string(zmask[z]) +
          ") could not be traced back to a real zone");
  return zrank;
}

} // namespace

Decomposition build_decomposition(SOA_Idx::Mesh const &src, Zone_Map z2r,
    int nparts, Decompose_Options const &opt) {
  int const pll = src.points.size(), ell = src.edges.size();
  int const fll = src.faces.size(), zll = src.zones.size();
  int const cll = src.corners.size(), sll = src.sides.size();
  int const all = src.iotas.size();

  if (nparts < 1)
    throw std::invalid_argument("build_decomposition: nparts must be positive");
  if (z2r.size() != static_cast<std::size_t>(zll))
    throw std::invalid_argument(
        "build_decomposition: the zone map is the wrong length");
  if (has_ghosts(src))
    throw std::invalid_argument(
        "build_decomposition: the input mesh already has ghosts");

  auto const &zmask = src.zones.mask;
  auto const &cmask = src.corners.mask;
  auto const &smask = src.sides.mask;
  auto const &fmask = src.faces.mask;
  auto const &amask = src.iotas.mask;
  for (int z = 0; z < zll; ++z)
    if (zmask[z] >= 1 &&
        (z2r[static_cast<std::size_t>(z)] < 0 ||
            z2r[static_cast<std::size_t>(z)] >= nparts))
      throw std::invalid_argument("build_decomposition: zone " +
          std::to_string(z) + " is assigned to part " +
          std::to_string(z2r[static_cast<std::size_t>(z)]) + ", outside [0, " +
          std::to_string(nparts) + ")");

  auto const &c2z = src.ds->caccess_intv("m:c>z");
  auto const &c2p = src.ds->caccess_intv("m:c>p");
  auto const &s2z = src.ds->caccess_intv("m:s>z");
  auto const &s2f = src.ds->caccess_intv("m:s>f");
  auto const &s2e = src.ds->caccess_intv("m:s>e");
  auto const &s2p1 = src.ds->caccess_intv("m:s>p1");
  auto const &s2p2 = src.ds->caccess_intv("m:s>p2");
  auto const &a2z = src.ds->caccess_intv("m:a>z");
  auto const &f2z1 = src.ds->caccess_intv("m:f>z1");
  auto const &f2z2 = src.ds->caccess_intv("m:f>z2");

  std::vector<int> const zrank = zone_ranks(src, z2r);
  auto rank_of_zone = [&zrank](int z) {
    return zrank[static_cast<std::size_t>(z)];
  };

  /* Locals.  Zones, corners, sides and iotas go where their zone goes; a
     point, edge or face is present wherever a side names it, which is every
     rank holding a zone that touches it. */
  Pairs zp, cp, sp, ap, pp, ep, fp;
  for (int z = 0; z < zll; ++z)
    if (zmask[z] != 0)
      zp.emplace_back(z, rank_of_zone(z));
  for (int c = 0; c < cll; ++c)
    if (cmask[c] != 0)
      cp.emplace_back(c, rank_of_zone(c2z[static_cast<std::size_t>(c)]));
  for (int a = 0; a < all; ++a)
    if (amask[a] != 0)
      ap.emplace_back(a, rank_of_zone(a2z[static_cast<std::size_t>(a)]));
  for (int s = 0; s < sll; ++s) {
    if (smask[s] == 0)
      continue;
    int const r = rank_of_zone(s2z[static_cast<std::size_t>(s)]);
    sp.emplace_back(s, r);
    pp.emplace_back(s2p1[static_cast<std::size_t>(s)], r);
    pp.emplace_back(s2p2[static_cast<std::size_t>(s)], r);
    ep.emplace_back(s2e[static_cast<std::size_t>(s)], r);
    fp.emplace_back(s2f[static_cast<std::size_t>(s)], r);
  }

  Decomposition d;
  d.nparts = nparts;
  d.z2r = std::move(z2r);
  d.opt = opt;

  struct Kind {
    Pairs *pairs;
    Entity_Layout *layout;
    int nglobal;
  };
  Kind const kinds[] = {{&pp, &d.points, pll}, {&ep, &d.edges, ell},
      {&fp, &d.faces, fll}, {&zp, &d.zones, zll}, {&cp, &d.corners, cll},
      {&sp, &d.sides, sll}, {&ap, &d.iotas, all}};
  for (Kind const &k : kinds) {
    sort_unique(*k.pairs);
    assign_owners(*k.pairs, k.nglobal, nparts, opt.owner_rule, k.layout->owner);
    to_csr(*k.pairs, nparts, k.layout->start, k.layout->gid);
    /* Nothing has ghosts until the next stage, and the empty CSR has to be
       well formed before then: the ghost search asks `local_index` whether a
       point or face is already here. */
    k.layout->gstart.assign(static_cast<std::size_t>(nparts) + 1, 0);
  }

  /* Ghost zones.  A cut face has a real zone on each side on different ranks,
     and each rank imports the other's; the point-adjacent depth widens that to
     every zone sharing a point, which is a strict superset. */
  Pairs gz;
  switch (opt.ghosts) {
  case Ghost_Depth::NONE:
    break;
  case Ghost_Depth::FACE_ADJACENT:
    for (int f = 0; f < fll; ++f) {
      if (fmask[f] == 0)
        continue;
      int const z1 = f2z1[static_cast<std::size_t>(f)];
      int const z2 = f2z2[static_cast<std::size_t>(f)];
      if (z1 < 0 || z2 < 0 || z1 >= zll || z2 >= zll)
        continue;
      if (zmask[z1] < 1 || zmask[z2] < 1 ||
          rank_of_zone(z1) == rank_of_zone(z2))
        continue;
      gz.emplace_back(z2, rank_of_zone(z1));
      gz.emplace_back(z1, rank_of_zone(z2));
    }
    break;
  case Ghost_Depth::POINT_ADJACENT: {
    auto const &p2zs = src.ds->caccess_intrr("m:p>zs");
    for (int r = 0; r < nparts; ++r)
      for (int i = d.points.start[static_cast<std::size_t>(r)];
          i < d.points.start[static_cast<std::size_t>(r) + 1]; ++i)
        for (int const z : p2zs[d.points.gid[static_cast<std::size_t>(i)]])
          if (z >= 0 && z < zll && zmask[z] >= 1 && rank_of_zone(z) != r)
            gz.emplace_back(z, r);
    break;
  }
  }
  sort_unique(gz);

  /* Ghost sides and corners hang off the ghost zones: the whole side ring on
     every face this rank already holds, so that each (face, edge) still has
     its two sides, and one corner per point of the ghost zone that is local
     here, so that every imported side's `m:s>c1/c2` lands on something.

     Both conditions are "is it already here", so nothing is imported that
     names a point, edge or face this rank does not hold -- which is what lets
     there be no ghost points, edges or faces at all.  For a face-adjacent
     ghost that reduces to the plan's "the ring on the cut face", since a cut
     face is exactly a face of the ghost zone that this rank holds. */
  Zone_Elements const zone_sides{"side", s2z, smask, zll};
  Zone_Elements const zone_corners{"corner", c2z, cmask, zll};
  Pairs gs, gc;
  for (auto const &[z, r] : gz) {
    for (int const s : zone_sides[z])
      if (d.faces.local_index(r, s2f[static_cast<std::size_t>(s)]) >= 0)
        gs.emplace_back(s, r);
    for (int const c : zone_corners[z])
      if (d.points.local_index(r, c2p[static_cast<std::size_t>(c)]) >= 0)
        gc.emplace_back(c, r);
  }
  sort_unique(gs);
  sort_unique(gc);

  to_csr(gz, nparts, d.zones.gstart, d.zones.ggid);
  to_csr(gs, nparts, d.sides.gstart, d.sides.ggid);
  to_csr(gc, nparts, d.corners.gstart, d.corners.ggid);

  return d;
}

} // namespace Part
} // namespace Ume
