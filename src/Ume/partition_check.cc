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
  \file Ume/partition_check.cc
*/

#include "Ume/partition_check.hh"
#include <algorithm>
#include <array>
#include <ostream>
#include <set>
#include <vector>

namespace Ume {
namespace Part {

using SOA_Idx::Entity;
using SOA_Idx::Mesh;

namespace {

using INTV = DS_Types::INTV_T;
using VEC3V = DS_Types::VEC3V_T;
using Coord = DS_Types::VEC3_T;

//! The coordinates an entity is built on, sorted and deduplicated
using Point_Set = std::vector<Coord>;

enum Kind { POINT, EDGE, FACE, ZONE, CORNER, SIDE, IOTA, NKINDS };

constexpr std::array<char const *, NKINDS> KIND_NAME{
    "point", "edge", "face", "zone", "corner", "side", "iota"};

Entity const &entity(Mesh const &m, int const k) {
  switch (k) {
  case POINT:
    return m.points;
  case EDGE:
    return m.edges;
  case FACE:
    return m.faces;
  case ZONE:
    return m.zones;
  case CORNER:
    return m.corners;
  case SIDE:
    return m.sides;
  default:
    return m.iotas;
  }
}

Comm::Neighbors const &neighbors_of(Entity const &e, bool const copies) {
  return copies ? e.myCpys : e.mySrcs;
}

//! A total order on coordinates, so a point set can be sorted and compared
/*! Exact, and meant to be: two ranks holding the same point hold the same
    bits, because `extract` copies `pcoord` rather than recomputing it. */
bool coord_less(Coord const &a, Coord const &b) {
  for (std::size_t k = 0; k < 3; ++k) {
    if (a[k] < b[k])
      return true;
    if (b[k] < a[k])
      return false;
  }
  return false;
}

bool set_less(Point_Set const &a, Point_Set const &b) {
  return std::lexicographical_compare(
      a.begin(), a.end(), b.begin(), b.end(), coord_less);
}

//! Is every coordinate of `part` also in `whole`?  Both must be sorted.
bool contains(Point_Set const &whole, Point_Set const &part) {
  return std::includes(
      whole.begin(), whole.end(), part.begin(), part.end(), coord_less);
}

std::ostream &operator<<(std::ostream &os, Point_Set const &s) {
  os << '{';
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (i == 3 && s.size() > 4) {
      os << " ... (" << s.size() << " points)";
      break;
    }
    os << (i ? " " : "") << s[i];
  }
  return os << '}';
}

//! What identifies one entity: its own points, then its container's
/*! One point set is not enough for every kind.  Two sides of a zone that
    share an edge across two different faces have the same endpoints, and two
    ghost corners of two ghost zones meeting at a local point have the same
    point -- so an entity that is named by where it sits has to carry where it
    sits.  Components are compared in order and each stands on its own, which
    is what lets the containing zone of a *ghost* be compared by containment
    while the entity's own points still have to match exactly. */
using Signature = std::vector<Point_Set>;

bool sig_less(Signature const &a, Signature const &b) {
  return std::lexicographical_compare(
      a.begin(), a.end(), b.begin(), b.end(), set_less);
}

std::ostream &operator<<(std::ostream &os, Signature const &s) {
  for (std::size_t i = 0; i < s.size(); ++i)
    os << (i ? " in " : "") << s[i];
  return os;
}

/*! What identifies an entity of one mesh, in terms that another mesh holding
    the same entity will agree with.  See the header for why this is
    coordinates and not derived geometry. */
class Identity {
public:
  explicit Identity(Mesh const &m);

  //! The signature of entity `i` of `kind`, or empty if a map is out of range
  /*! An empty answer means "cannot say", not "no points": every kind here has
      at least one point.  Callers skip the comparison rather than fail it,
      because an out-of-range map is `mesh_audit`'s finding to report, not
      this file's, and reporting it twice buries the cross-rank problem. */
  Signature operator()(int kind, int i) const;

private:
  //! The distinct points of the real sides on face `f`
  Point_Set face_points(int f) const;
  //! The distinct points of the corners of zone `z`
  Point_Set zone_points(int z) const;

  VEC3V const *pc_;
  INTV const *e2p1_, *e2p2_, *c2p_, *c2z_, *s2p1_, *s2p2_, *s2z_, *s2f_, *a2p_,
      *a2z_;
  std::vector<short> const *smask_;
  std::array<int, NKINDS> n_;
  //! face -> the real sides on it, and zone -> its corners, as CSR
  std::vector<int> fstart_, fsides_, zstart_, zcorners_;
};

//! Group `key[i]` -> `i` into a CSR over `nkeys` keys, for `i` in `[0, n)`
template <class Keep>
void group_by(INTV const &key, int const n, int const nkeys, Keep const &keep,
    std::vector<int> &start, std::vector<int> &value) {
  start.assign(static_cast<std::size_t>(nkeys) + 1, 0);
  auto usable = [&](int const i) {
    return keep(i) && key[static_cast<std::size_t>(i)] >= 0 &&
        key[static_cast<std::size_t>(i)] < nkeys;
  };
  for (int i = 0; i < n; ++i)
    if (usable(i))
      start[static_cast<std::size_t>(key[static_cast<std::size_t>(i)]) + 1] += 1;
  for (int k = 0; k < nkeys; ++k)
    start[static_cast<std::size_t>(k) + 1] += start[static_cast<std::size_t>(k)];
  value.resize(static_cast<std::size_t>(start.back()));
  std::vector<int> fill(start.begin(), start.end() - 1);
  for (int i = 0; i < n; ++i)
    if (usable(i))
      value[static_cast<std::size_t>(
          fill[static_cast<std::size_t>(key[static_cast<std::size_t>(i)])]++)] =
          i;
}

Identity::Identity(Mesh const &m)
    : pc_{&m.ds->caccess_vec3v("pcoord")},
      e2p1_{&m.ds->caccess_intv("m:e>p1")},
      e2p2_{&m.ds->caccess_intv("m:e>p2")}, c2p_{&m.ds->caccess_intv("m:c>p")},
      c2z_{&m.ds->caccess_intv("m:c>z")}, s2p1_{&m.ds->caccess_intv("m:s>p1")},
      s2p2_{&m.ds->caccess_intv("m:s>p2")}, s2z_{&m.ds->caccess_intv("m:s>z")},
      s2f_{&m.ds->caccess_intv("m:s>f")},
      a2p_{m.iotas.size() > 0 ? &m.ds->caccess_intv("m:a>p") : nullptr},
      a2z_{m.iotas.size() > 0 ? &m.ds->caccess_intv("m:a>z") : nullptr},
      smask_{&m.sides.mask} {
  for (int k = 0; k < NKINDS; ++k)
    n_[static_cast<std::size_t>(k)] = entity(m, k).size();

  /* A face is identified by the ring of real sides on it.  One zone's ring
     spans every edge of the face, so a rank that holds either zone gets the
     same answer -- and so does the serial mesh, where both rings are present
     and their union is the same set of points.  Phantom sides are left out
     for the same reason on every mesh, so they cancel. */
  group_by(*s2f_, n_[SIDE], n_[FACE],
      [this](int const s) { return (*smask_)[static_cast<std::size_t>(s)] >= 1; },
      fstart_, fsides_);

  /* Ghost corners are kept: they are how a ghost zone is identified at all. */
  group_by(
      *c2z_, n_[CORNER], n_[ZONE], [](int) { return true; }, zstart_, zcorners_);
}

//! Put a point set in the canonical form the comparisons assume
void canonicalize(Point_Set &s) {
  std::sort(s.begin(), s.end(), coord_less);
  s.erase(std::unique(s.begin(), s.end()), s.end());
}

Point_Set Identity::face_points(int const f) const {
  Point_Set out;
  if (f < 0 || f >= n_[FACE])
    return out;
  for (int j = fstart_[static_cast<std::size_t>(f)];
      j < fstart_[static_cast<std::size_t>(f) + 1]; ++j) {
    auto const s = static_cast<std::size_t>(fsides_[static_cast<std::size_t>(j)]);
    int const p[2] = {(*s2p1_)[s], (*s2p2_)[s]};
    for (int const q : p) {
      if (q < 0 || q >= n_[POINT])
        return Point_Set{};
      out.push_back((*pc_)[static_cast<std::size_t>(q)]);
    }
  }
  canonicalize(out);
  return out;
}

Point_Set Identity::zone_points(int const z) const {
  Point_Set out;
  if (z < 0 || z >= n_[ZONE])
    return out;
  for (int j = zstart_[static_cast<std::size_t>(z)];
      j < zstart_[static_cast<std::size_t>(z) + 1]; ++j) {
    int const q =
        (*c2p_)[static_cast<std::size_t>(zcorners_[static_cast<std::size_t>(j)])];
    if (q < 0 || q >= n_[POINT])
      return Point_Set{};
    out.push_back((*pc_)[static_cast<std::size_t>(q)]);
  }
  canonicalize(out);
  return out;
}

Signature Identity::operator()(int const kind, int const i) const {
  if (i < 0 || i >= n_[static_cast<std::size_t>(kind)])
    return Signature{};

  auto at = [i](INTV const *v) { return (*v)[static_cast<std::size_t>(i)]; };
  Point_Set own;
  bool ok{true};
  auto add = [&](int const p) {
    if (p < 0 || p >= n_[POINT])
      ok = false;
    else
      own.push_back((*pc_)[static_cast<std::size_t>(p)]);
  };

  Signature sig;
  switch (kind) {
  case POINT:
    add(i);
    break;
  case EDGE:
    add(at(e2p1_));
    add(at(e2p2_));
    break;
  case FACE:
    own = face_points(i);
    break;
  case ZONE:
    own = zone_points(i);
    break;
  case CORNER:
    /* Two ghost zones meeting at a local point give that point two ghost
       corners, so the zone is part of the answer to "which corner". */
    add(at(c2p_));
    sig.push_back(zone_points(at(c2z_)));
    break;
  case SIDE:
    /* Likewise two sides of one zone on either side of an edge: same two
       endpoints, different faces. */
    add(at(s2p1_));
    add(at(s2p2_));
    sig.push_back(face_points(at(s2f_)));
    sig.push_back(zone_points(at(s2z_)));
    break;
  default: // IOTA
    if (a2p_ == nullptr || a2z_ == nullptr)
      ok = false;
    else {
      add(at(a2p_));
      sig.push_back(zone_points(at(a2z_)));
    }
    break;
  }

  canonicalize(own);
  sig.insert(sig.begin(), std::move(own));
  if (!ok || std::ranges::any_of(sig, [](Point_Set const &s) {
        return s.empty();
      }))
    return Signature{};
  return sig;
}

/*! Accumulates problems for one named check, printing at most
    `Check_Options::max_reports` of them.  Deliberately the same shape as
    `mesh_audit`'s `Auditor`, so the two read alike in a driver's output. */
class Checker {
public:
  Checker(Mesh const *serial, std::span<Mesh const *const> ranks,
      std::ostream &os, Check_Options const &opts);

  Check_Result run();

private:
  void start(char const *const name) {
    finish();
    name_ = name;
    hits_ = 0;
    res_.checks += 1;
  }
  void finish() {
    if (name_ == nullptr)
      return;
    if (hits_ > 0) {
      res_.failed_checks += 1;
      if (hits_ > opt_.max_reports)
        os_ << "  [" << name_ << "] ... and " << hits_ - opt_.max_reports
            << " more\n";
    }
    name_ = nullptr;
  }
  //! Report one problem; the caller streams the detail and a newline
  std::ostream &bad() {
    hits_ += 1;
    res_.problems += 1;
    if (hits_ > opt_.max_reports)
      return null_;
    return os_ << "  [" << name_ << "] ";
  }
  bool enough() const { return hits_ > opt_.max_reports; }

  //! `rank 2 face 17`, the way every message here names an entity
  static std::ostream &where(
      std::ostream &os, int const r, int const k, int const i) {
    return os << "rank " << r << ' ' << KIND_NAME[static_cast<std::size_t>(k)]
              << ' ' << i;
  }
  static char const *kname(int const k) {
    return KIND_NAME[static_cast<std::size_t>(k)];
  }

  static constexpr bool in(int const i, int const n) { return i >= 0 && i < n; }
  //! Is `i` of `kind` a ghost on rank `r`?
  bool is_ghost(int const r, int const k, int const i) const {
    return i >= entity(*ranks_[static_cast<std::size_t>(r)], k).local_size();
  }

  //! The `mySrcs`/`myCpys` entry on `pe` that answers one on `r`, or nullptr
  Comm::Neighbor const *counterpart(
      int const r, int const pe, int const k, bool const copies) const;

  bool ranks_line_up();
  void check_neighbors();
  void check_pairing();
  void check_ghosts();
  void check_reassembly();

  //! Compare two instances of what should be one entity
  /*! `ghost` says the first is an imported fragment rather than the whole
      thing: a rank takes only the part of a ghost zone that touches points it
      already holds, so its point set is a subset of the source's -- a proper
      one whenever the zone straddles the cut.  An unidentifiable entity
      passes, quietly; `mesh_audit` is where that gets reported. */
  static bool agree(
      Signature const &a, Signature const &b, bool const ghost) {
    if (a.empty() || b.empty() || a.size() != b.size())
      return true;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (ghost ? !contains(b[i], a[i]) : a[i] != b[i])
        return false;
    return true;
  }

  Mesh const *serial_;
  std::span<Mesh const *const> ranks_;
  std::ostream &os_;
  Check_Options opt_;
  Check_Result res_;
  char const *name_{nullptr};
  int hits_{0};
  std::ostream null_{nullptr};
  int const n_;
  std::vector<Identity> id_;
};

Checker::Checker(Mesh const *const serial,
    std::span<Mesh const *const> const ranks, std::ostream &os,
    Check_Options const &opts)
    : serial_{serial}, ranks_{ranks}, os_{os}, opt_{opts},
      n_{static_cast<int>(ranks.size())} {}

Comm::Neighbor const *Checker::counterpart(
    int const r, int const pe, int const k, bool const copies) const {
  /* A copy list on r is answered by a source list on pe, and vice versa. */
  Comm::Neighbors const &other =
      neighbors_of(entity(*ranks_[static_cast<std::size_t>(pe)], k), !copies);
  Comm::Neighbor const *found{nullptr};
  for (Comm::Neighbor const &nbr : other) {
    if (nbr.pe != r)
      continue;
    if (found != nullptr)
      return nullptr; // ambiguous; check_neighbors() reports it
    found = &nbr;
  }
  return found;
}

bool Checker::ranks_line_up() {
  start("the ranks are a complete set");
  if (ranks_.empty()) {
    bad() << "no partitions were given\n";
    finish();
    return false;
  }
  for (int r = 0; r < n_; ++r) {
    Mesh const &m = *ranks_[static_cast<std::size_t>(r)];
    if (m.mype != r)
      bad() << "the mesh at position " << r << " says it is rank " << m.mype
            << '\n';
    if (m.numpe != n_)
      bad() << "rank " << r << " says there are " << m.numpe << " ranks, but "
            << n_ << " were given\n";
  }
  bool const ok = hits_ == 0;
  finish();
  return ok;
}

void Checker::check_neighbors() {
  start("neighbor lists pair up");
  for (int r = 0; r < n_ && !enough(); ++r)
    for (int k = 0; k < NKINDS && !enough(); ++k) {
      Entity const &e = entity(*ranks_[static_cast<std::size_t>(r)], k);
      for (bool const copies : {true, false}) {
        char const *const which = copies ? "myCpys" : "mySrcs";
        std::set<int> seen;
        for (Comm::Neighbor const &nbr : neighbors_of(e, copies)) {
          if (enough())
            break;
          if (!in(nbr.pe, n_) || nbr.pe == r) {
            bad() << "rank " << r << ' ' << kname(k) << ' ' << which
                  << " names rank " << nbr.pe << '\n';
            continue;
          }
          /* Two entries for one rank would post two sends with the same tag,
             and `Comm_MPI` matches a receive on (pe, tag) alone. */
          if (!seen.insert(nbr.pe).second) {
            bad() << "rank " << r << ' ' << kname(k) << ' ' << which
                  << " has two entries for rank " << nbr.pe << '\n';
            continue;
          }
          Comm::Neighbor const *const other = counterpart(r, nbr.pe, k, copies);
          if (other == nullptr) {
            bad() << "rank " << r << ' ' << kname(k) << ' ' << which
                  << " with rank " << nbr.pe << " has no unique counterpart in "
                  << (copies ? "mySrcs" : "myCpys") << " there\n";
          } else if (other->elements.size() != nbr.elements.size()) {
            bad() << "rank " << r << ' ' << kname(k) << ' ' << which
                  << " with rank " << nbr.pe << " has " << nbr.elements.size()
                  << " elements, its counterpart there has "
                  << other->elements.size() << '\n';
          }
        }
      }
    }
  finish();
}

void Checker::check_pairing() {
  /* Walking from the copy side visits each pair once.  `check_neighbors` has
     already established the pairing, so anything it rejected is skipped here
     rather than reported twice. */
  for (int pass = 0; pass < 2; ++pass) {
    /* The second pass is the one that matters.  `Comm::Buffers` packs and
       unpacks by position, so a pairing that agrees on length and disagrees
       on order moves every value to the wrong entity, quietly. */
    start(pass == 0 ? "copies and sources agree on the comm type"
                    : "paired copies and sources are the same entity");
    for (int r = 0; r < n_ && !enough(); ++r)
      for (int k = 0; k < NKINDS && !enough(); ++k) {
        Entity const &e = entity(*ranks_[static_cast<std::size_t>(r)], k);
        for (Comm::Neighbor const &nbr : e.myCpys) {
          if (enough())
            break;
          if (!in(nbr.pe, n_) || nbr.pe == r)
            continue;
          Comm::Neighbor const *const other = counterpart(r, nbr.pe, k, true);
          if (other == nullptr || other->elements.size() != nbr.elements.size())
            continue;
          Entity const &oe =
              entity(*ranks_[static_cast<std::size_t>(nbr.pe)], k);
          for (std::size_t j = 0; j < nbr.elements.size() && !enough(); ++j) {
            int const li = nbr.elements[j], lo = other->elements[j];
            if (!in(li, e.size()) || !in(lo, oe.size()))
              continue;
            if (pass == 0) {
              int const t = e.comm_type[static_cast<std::size_t>(li)];
              if (t != Entity::COPY && t != Entity::GHOST)
                where(bad(), r, k, li) << " is a copy of rank " << nbr.pe
                                       << " but its comm type is " << t << '\n';
              int const ot = oe.comm_type[static_cast<std::size_t>(lo)];
              if (ot != Entity::SOURCE)
                where(bad(), nbr.pe, k, lo) << " sources a copy on rank " << r
                                            << " but its comm type is " << ot
                                            << '\n';
              continue;
            }
            Signature const a = id_[static_cast<std::size_t>(r)](k, li);
            Signature const b = id_[static_cast<std::size_t>(nbr.pe)](k, lo);
            if (!agree(a, b, is_ghost(r, k, li)))
              where(bad(), r, k, li)
                  << ' ' << a << " is element " << j << " of its copy list for "
                  << "rank " << nbr.pe << ", which offers " << kname(k) << ' '
                  << lo << ' ' << b << '\n';
          }
        }
      }
  }
  finish();
}

void Checker::check_ghosts() {
  start("ghost sources are the same entity, real, on the source rank");
  for (int r = 0; r < n_ && !enough(); ++r)
    for (int k = 0; k < NKINDS && !enough(); ++k) {
      Entity const &e = entity(*ranks_[static_cast<std::size_t>(r)], k);
      int const nl = e.local_size(), ng = e.ghost_local_size();
      if (static_cast<int>(e.cpy_idx.size()) != ng) {
        bad() << "rank " << r << " has " << e.cpy_idx.size()
              << " copy records for " << ng << ' ' << kname(k) << " ghosts\n";
        continue;
      }
      for (int g = 0; g < ng && !enough(); ++g) {
        auto const ug = static_cast<std::size_t>(g);
        int const l = nl + g;
        if (e.cpy_idx[ug] != l) {
          where(bad(), r, k, l)
              << " is ghost " << g << ", but cpy_idx says " << e.cpy_idx[ug]
              << "; ghosts are the tail of the index range, in order\n";
          continue;
        }
        int const pe = e.src_pe[ug], si = e.src_idx[ug];
        if (!in(pe, n_) || pe == r) {
          where(bad(), r, k, l) << " has source rank " << pe << '\n';
          continue;
        }
        Entity const &oe = entity(*ranks_[static_cast<std::size_t>(pe)], k);
        if (!in(si, oe.local_size())) {
          where(bad(), r, k, l)
              << " has source index " << si << " on rank " << pe
              << ", outside its real range [0, " << oe.local_size() << ")\n";
          continue;
        }
        /* Without `ghost_entity_comm` the source of a ghost gets no neighbor
           entry and so stays INTERNAL; what is never right is for it to be
           somebody else's copy. */
        int const t = oe.comm_type[static_cast<std::size_t>(si)];
        if (t == Entity::COPY || t == Entity::GHOST)
          where(bad(), r, k, l)
              << " is sourced by " << kname(k) << ' ' << si << " on rank " << pe
              << ", whose comm type is " << t << '\n';
        if (e.ghost_mask[ug] != oe.mask[static_cast<std::size_t>(si)])
          where(bad(), r, k, l)
              << " has ghost_mask " << e.ghost_mask[ug]
              << ", but its source on rank " << pe << " has mask "
              << oe.mask[static_cast<std::size_t>(si)] << '\n';
        Signature const a = id_[static_cast<std::size_t>(r)](k, l);
        Signature const b = id_[static_cast<std::size_t>(pe)](k, si);
        if (!agree(a, b, true))
          where(bad(), r, k, l)
              << ' ' << a << " is sourced by " << kname(k) << ' ' << si
              << " on rank " << pe << ", which is " << b << '\n';
      }
    }
  finish();
}

void Checker::check_reassembly() {
  Identity const whole{*serial_};

  /* An entity is counted on the one rank that is authoritative for it: the
     SOURCE of a shared point, the owner of a ghosted zone.  Every other
     instance is a COPY or a GHOST, so the rule is the same for both
     mechanisms and needs no per-kind special casing. */
  auto authoritative = [](Entity const &e, int const i) {
    int const t = e.comm_type[static_cast<std::size_t>(i)];
    return t != Entity::COPY && t != Entity::GHOST;
  };

  start("every entity is accounted for exactly once");
  std::array<std::vector<Signature>, NKINDS> found;
  for (int k = 0; k < NKINDS; ++k) {
    long total{0};
    for (int r = 0; r < n_; ++r) {
      Entity const &e = entity(*ranks_[static_cast<std::size_t>(r)], k);
      for (int i = 0; i < e.size(); ++i)
        if (authoritative(e, i)) {
          total += 1;
          found[static_cast<std::size_t>(k)].push_back(
              id_[static_cast<std::size_t>(r)](k, i));
        }
    }
    long const want = entity(*serial_, k).size();
    if (total != want)
      bad() << "the ranks hold " << total << ' '
            << kname(k) << "s between them, the serial mesh has " << want
            << '\n';
  }

  /* Counting right is not the same as covering right: a partition that drops
     one zone and duplicates another has the same total. */
  start("the ranks cover the serial mesh");
  for (int k = 0; k < NKINDS && !enough(); ++k) {
    std::vector<Signature> want;
    Entity const &se = entity(*serial_, k);
    want.reserve(static_cast<std::size_t>(se.size()));
    for (int i = 0; i < se.size(); ++i)
      want.push_back(whole(k, i));

    std::vector<Signature> &have = found[static_cast<std::size_t>(k)];
    std::sort(have.begin(), have.end(), sig_less);
    std::sort(want.begin(), want.end(), sig_less);

    /* Walk the two sorted lists together: what is in one and not the other is
       an entity that got lost, or one that got counted twice. */
    std::size_t i{0}, j{0};
    while ((i < have.size() || j < want.size()) && !enough()) {
      if (j == want.size() || (i < have.size() && sig_less(have[i], want[j]))) {
        bad() << "the ranks hold a " << kname(k) << ' ' << have[i]
              << " that the serial mesh does not\n";
        ++i;
      } else if (i == have.size() || sig_less(want[j], have[i])) {
        bad() << "no rank holds the " << kname(k) << ' ' << want[j]
              << " of the serial mesh\n";
        ++j;
      } else {
        ++i;
        ++j;
      }
    }
  }
  finish();
}

Check_Result Checker::run() {
  if (!ranks_line_up())
    return res_;

  id_.reserve(ranks_.size());
  for (Mesh const *const m : ranks_)
    id_.emplace_back(*m);

  if (opt_.neighbors)
    check_neighbors();
  if (opt_.pairing)
    check_pairing();
  if (opt_.ghosts)
    check_ghosts();
  if (opt_.reassembly && serial_ != nullptr)
    check_reassembly();
  finish();
  return res_;
}

} // namespace

Check_Result check_partitions(std::span<Mesh const *const> const ranks,
    std::ostream &os, Check_Options const &opts) {
  return Checker{nullptr, ranks, os, opts}.run();
}

Check_Result check_partitions(Mesh const &serial,
    std::span<Mesh const *const> const ranks, std::ostream &os,
    Check_Options const &opts) {
  return Checker{&serial, ranks, os, opts}.run();
}

} // namespace Part
} // namespace Ume
