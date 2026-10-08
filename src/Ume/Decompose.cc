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
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

namespace Ume {
namespace Part {

char const *kind_name(int const kind) {
  constexpr char const *NAMES[NKINDS]{
      "point", "edge", "face", "zone", "corner", "side", "iota"};
  return (kind >= 0 && kind < NKINDS) ? NAMES[kind] : "?";
}

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

  /* A point, edge or face that no side names has no zone to follow, and
     without one it would be enrolled nowhere: the entity would vanish from
     every rank and every index above it would shift down, so even a one-rank
     decomposition would stop reproducing its input.  The audit counts these
     as non-fatal, so a mesh that carries them is one this has to carry
     through.  Rank 0 is as good a home as any -- nothing references them, so
     the choice costs no communication -- and at `nparts == 1` it is the only
     rank, which is what keeps the identity. */
  auto keep_unreferenced = [](Pairs &p, auto const &mask, int const n) {
    std::vector<char> named(static_cast<std::size_t>(n), 0);
    for (auto const &[g, r] : p)
      if (g >= 0 && g < n)
        named[static_cast<std::size_t>(g)] = 1;
    for (int g = 0; g < n; ++g)
      if (mask[g] != 0 && !named[static_cast<std::size_t>(g)])
        p.emplace_back(g, 0);
  };
  keep_unreferenced(pp, src.points.mask, pll);
  keep_unreferenced(ep, src.edges.mask, ell);
  keep_unreferenced(fp, fmask, fll);

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

/* ------------------------------- extract --------------------------------- */

namespace {

//! For each global entity, the ranks holding it, ascending
/*! The layouts are keyed the other way round -- rank to entities -- which is
    what a rank's own mesh needs.  Turning one round is what the *other* half
    of `extract` needs: which ranks to tell about an entity this rank sources,
    and how many ranks share it at all. */
class Presence {
public:
  Presence(std::vector<int> const &start, std::vector<int> const &gid,
      int nparts, int nglobal) {
    start_.assign(static_cast<std::size_t>(nglobal) + 1, 0);
    for (int const g : gid)
      start_[static_cast<std::size_t>(g) + 1] += 1;
    for (int g = 0; g < nglobal; ++g)
      start_[static_cast<std::size_t>(g) + 1] +=
          start_[static_cast<std::size_t>(g)];

    rank_.resize(gid.size());
    std::vector<int> fill(start_.begin(), start_.end() - 1);
    for (int r = 0; r < nparts; ++r)
      for (int i = start[static_cast<std::size_t>(r)];
          i < start[static_cast<std::size_t>(r) + 1]; ++i)
        rank_[static_cast<std::size_t>(
            fill[static_cast<std::size_t>(gid[static_cast<std::size_t>(i)])]++)] =
            r;
  }

  std::span<int const> operator[](int g) const {
    auto const b = static_cast<std::size_t>(start_[static_cast<std::size_t>(g)]);
    auto const e =
        static_cast<std::size_t>(start_[static_cast<std::size_t>(g) + 1]);
    return std::span<int const>{rank_.data() + b, e - b};
  }

private:
  std::vector<int> start_, rank_;
};

//! Global index to this rank's local index, or -1
class Renumber {
public:
  Renumber(Entity_Layout const &l, int pe, int nglobal)
      : g2l_(static_cast<std::size_t>(nglobal), -1) {
    for (int i = 0; i < l.size(pe); ++i)
      g2l_[static_cast<std::size_t>(l.global_index(pe, i))] = i;
  }
  int operator[](int g) const {
    return (g >= 0 && g < static_cast<int>(g2l_.size()))
        ? g2l_[static_cast<std::size_t>(g)]
        : -1;
  }
  //! The same, but a map that fails to resolve is a bug rather than a value
  int require(int g, char const *what, int at) const {
    int const l = (*this)[g];
    if (l < 0)
      throw std::logic_error(std::string("extract: ") + what + " of " +
          std::to_string(at) + " names " + std::to_string(g) +
          ", which this rank does not hold");
    return l;
  }

private:
  std::vector<int> g2l_;
};

//! Everything about one entity kind that does not involve its own maps
/*! Masks, comm types, the ghost records and the neighbor lists are the same
    shape of problem for all seven kinds; only whether the kind is duplicated
    on a cut or ghosted across one changes the answer. */
void fill_common(SOA_Idx::Entity const &se, SOA_Idx::Entity &de,
    Entity_Layout const &l, int pe, bool shared_kind, bool comm) {
  int const nl = l.local_size(pe), ng = l.ghost_local_size(pe);
  Presence const locals{l.start, l.gid, l.nparts(), l.nglobal()};
  Presence const ghosts{l.gstart, l.ggid, l.nparts(), l.nglobal()};

  std::map<int, std::vector<int>> cpys, srcs;

  for (int i = 0; i < nl; ++i) {
    int const g = l.global_index(pe, i);
    int const owner = l.owner[static_cast<std::size_t>(g)];
    de.mask[static_cast<std::size_t>(i)] = se.mask[static_cast<std::size_t>(g)];

    if (shared_kind) {
      /* A point, edge or face on a cut is a real entity on every rank that
         holds it; one of them is the SOURCE and the rest are COPYs. */
      auto const on = locals[g];
      if (on.size() == 1)
        de.comm_type[static_cast<std::size_t>(i)] = SOA_Idx::Entity::INTERNAL;
      else if (owner == pe) {
        de.comm_type[static_cast<std::size_t>(i)] = SOA_Idx::Entity::SOURCE;
        for (int const r : on)
          if (r != pe)
            srcs[r].push_back(i);
      } else {
        de.comm_type[static_cast<std::size_t>(i)] = SOA_Idx::Entity::COPY;
        cpys[owner].push_back(i);
      }
    } else {
      /* A zone, corner or side is real on exactly one rank; it is a SOURCE
         only when somebody else keeps a ghost of it. */
      auto const shadows = comm ? ghosts[g] : std::span<int const>{};
      de.comm_type[static_cast<std::size_t>(i)] = shadows.empty()
          ? SOA_Idx::Entity::INTERNAL
          : SOA_Idx::Entity::SOURCE;
      for (int const r : shadows)
        srcs[r].push_back(i);
    }
  }

  for (int j = 0; j < ng; ++j) {
    int const i = nl + j;
    int const g = l.global_index(pe, i);
    int const owner = l.owner[static_cast<std::size_t>(g)];
    de.mask[static_cast<std::size_t>(i)] = 0;
    de.comm_type[static_cast<std::size_t>(i)] = SOA_Idx::Entity::GHOST;
    de.cpy_idx[static_cast<std::size_t>(j)] = i;
    de.src_pe[static_cast<std::size_t>(j)] = owner;
    de.src_idx[static_cast<std::size_t>(j)] = l.local_index(owner, g);
    /* `renumbering.cc` reads ghost_mask == 0 as "null ghost", so the owner's
       mask is the right thing to propagate. */
    de.ghost_mask[static_cast<std::size_t>(j)] =
        se.mask[static_cast<std::size_t>(g)];
    if (comm)
      cpys[owner].push_back(i);
  }

  /* One Neighbor per rank, never two: `Comm_MPI` matches a message by (pe,
     tag), so a second entry for the same rank would pair with the wrong
     receive.  Both sides built their list in ascending global order, which
     makes the correspondence positional without anything being exchanged. */
  auto to_neighbors = [](std::map<int, std::vector<int>> &m) {
    Comm::Neighbors n;
    n.reserve(m.size());
    for (auto &[rank, elements] : m)
      n.push_back(Comm::Neighbor{rank, std::move(elements)});
    return n;
  };
  de.myCpys = to_neighbors(cpys);
  de.mySrcs = to_neighbors(srcs);
}

//! Carry `src`'s subsets through, dropping what this rank does not hold
void filter_subsets(
    SOA_Idx::Entity const &se, SOA_Idx::Entity &de, Renumber const &g2l) {
  de.subsets.clear();
  for (auto const &s : se.subsets) {
    SOA_Idx::Entity::Subset out;
    out.name = s.name;
    /* Sorting by the new index is what keeps the subset's own locals ahead of
       its ghosts, which is what `lsize` means.  On a serial input, where no
       subset has a ghost in it, this leaves the order alone. */
    std::vector<std::pair<int, short>> kept;
    for (std::size_t i = 0; i < s.elements.size(); ++i) {
      int const l = g2l[s.elements[i]];
      if (l >= 0)
        kept.emplace_back(l, i < s.mask.size() ? s.mask[i] : short{1});
    }
    std::ranges::sort(kept);
    for (auto const &[l, m] : kept) {
      out.elements.push_back(l);
      out.mask.push_back(m);
    }
    out.lsize = static_cast<int>(std::ranges::count_if(
        out.elements, [&de](int l) { return l < de.local_size(); }));
    de.subsets.push_back(std::move(out));
  }
}

} // namespace

void extract(SOA_Idx::Mesh const &src, Decomposition const &d, int pe,
    SOA_Idx::Mesh &dst) {
  if (pe < 0 || pe >= d.nparts)
    throw std::out_of_range("extract: rank " + std::to_string(pe) +
        " is not one of the " + std::to_string(d.nparts) + " parts");

  dst.ivtag = src.ivtag;
  dst.version_header = src.version_header;
  dst.mype = pe;
  dst.numpe = d.nparts;
  dst.geo = src.geo;
  dst.dump_iotas = src.dump_iotas;

  /* `resize` is virtual and sizes the stored maps as well as the masks, so
     nothing below may run before all seven have been sized. */
  struct Kind {
    SOA_Idx::Entity const *se;
    SOA_Idx::Entity *de;
    Entity_Layout const *layout;
    bool shared;
  };
  Kind const kinds[] = {{&src.points, &dst.points, &d.points, true},
      {&src.edges, &dst.edges, &d.edges, true},
      {&src.faces, &dst.faces, &d.faces, true},
      {&src.zones, &dst.zones, &d.zones, false},
      {&src.corners, &dst.corners, &d.corners, false},
      {&src.sides, &dst.sides, &d.sides, false},
      {&src.iotas, &dst.iotas, &d.iotas, false}};
  for (Kind const &k : kinds) {
    int const nl = k.layout->local_size(pe), ng = k.layout->ghost_local_size(pe);
    k.de->resize(nl, nl + ng, ng);
  }
  for (Kind const &k : kinds)
    fill_common(
        *k.se, *k.de, *k.layout, pe, k.shared, k.shared || d.opt.ghost_entity_comm);

  Renumber const p{d.points, pe, src.points.size()};
  Renumber const e{d.edges, pe, src.edges.size()};
  Renumber const f{d.faces, pe, src.faces.size()};
  Renumber const z{d.zones, pe, src.zones.size()};
  Renumber const c{d.corners, pe, src.corners.size()};
  Renumber const s{d.sides, pe, src.sides.size()};
  Renumber const a{d.iotas, pe, src.iotas.size()};

  //! Copy `name` from `src`, renumbering each entry through `to`
  auto remap = [&](char const *name, Renumber const &from, Renumber const &to,
                   char const *what, int n) {
    auto const &in = src.ds->caccess_intv(name);
    auto &out = dst.ds->access_intv(name);
    for (int g = 0; g < n; ++g) {
      int const l = from[g];
      if (l >= 0)
        out[static_cast<std::size_t>(l)] =
            to.require(in[static_cast<std::size_t>(g)], what, g);
    }
  };

  int const pll = src.points.size(), ell = src.edges.size();
  int const fll = src.faces.size(), cll = src.corners.size();
  int const sll = src.sides.size(), all = src.iotas.size();

  {
    auto const &in = src.ds->caccess_vec3v("pcoord");
    auto &out = dst.ds->access_vec3v("pcoord");
    for (int g = 0; g < pll; ++g)
      if (p[g] >= 0)
        out[static_cast<std::size_t>(p[g])] = in[static_cast<std::size_t>(g)];
  }

  remap("m:e>p1", e, p, "m:e>p1", ell);
  remap("m:e>p2", e, p, "m:e>p2", ell);

  /* `m:f>z1/z2` is the one map allowed to dangle: -1 means "no zone on that
     side of this face", which is what a cut face looks like when there is no
     ghost layer to resolve it to.  Which of the two was z1 is preserved. */
  for (char const *name : {"m:f>z1", "m:f>z2"}) {
    auto const &in = src.ds->caccess_intv(name);
    auto &out = dst.ds->access_intv(name);
    for (int g = 0; g < fll; ++g)
      if (f[g] >= 0)
        out[static_cast<std::size_t>(f[g])] = z[in[static_cast<std::size_t>(g)]];
  }

  remap("m:c>p", c, p, "m:c>p", cll);
  remap("m:c>z", c, z, "m:c>z", cll);

  remap("m:s>z", s, z, "m:s>z", sll);
  remap("m:s>f", s, f, "m:s>f", sll);
  remap("m:s>e", s, e, "m:s>e", sll);
  remap("m:s>p1", s, p, "m:s>p1", sll);
  remap("m:s>p2", s, p, "m:s>p2", sll);
  remap("m:s>c1", s, c, "m:s>c1", sll);
  remap("m:s>c2", s, c, "m:s>c2", sll);

  /*! The side-to-side maps are the only ones that can leave this rank.  A
      ghost side's `m:s>s5` is on a face of the ghost zone we did not import;
      it becomes -1, the one spelling of "no partner" the mesh format has
      (`m:f>z2` on a face with nothing behind it, `m:s>s5` on a phantom side).
      Nothing dereferences s5, and every `Auditor::connectivity` check on it
      gates on `usable()`, which a ghost side fails.

      With `Ghost_Depth::NONE` a cut side's `m:s>s2` has nothing to pair with
      at all.  That one stays a self-reference, not -1: `calc_face_area` and
      `VAR_point_norm` index through s2 without a bounds check, and NONE is a
      diagnostic that still has to be computable.  The audit's `(face, edge)`
      pairing check is what fails for it.  A phantom side's s5 source value is
      already -1 and `s[-1]` is never read. */
  for (char const *name : {"m:s>s2", "m:s>s3", "m:s>s4", "m:s>s5"}) {
    bool const is_s5 = std::string_view{name} == "m:s>s5";
    auto const &in = src.ds->caccess_intv(name);
    auto &out = dst.ds->access_intv(name);
    for (int g = 0; g < sll; ++g) {
      int const l = s[g];
      if (l < 0)
        continue;
      int const gp = in[static_cast<std::size_t>(g)];
      int const partner = (gp >= 0) ? s[gp] : -1;
      out[static_cast<std::size_t>(l)] =
          (partner >= 0) ? partner : (is_s5 ? -1 : l);
    }
  }

  if (all > 0) {
    remap("m:a>z", a, z, "m:a>z", all);
    remap("m:a>f", a, f, "m:a>f", all);
    remap("m:a>p", a, p, "m:a>p", all);
    remap("m:a>e", a, e, "m:a>e", all);
    remap("m:a>s", a, s, "m:a>s", all);
  }

  if (d.opt.split_subsets) {
    filter_subsets(src.points, dst.points, p);
    filter_subsets(src.edges, dst.edges, e);
    filter_subsets(src.faces, dst.faces, f);
    filter_subsets(src.zones, dst.zones, z);
    filter_subsets(src.corners, dst.corners, c);
    filter_subsets(src.sides, dst.sides, s);
    filter_subsets(src.iotas, dst.iotas, a);
  }
}

} // namespace Part
} // namespace Ume
