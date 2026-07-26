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
  \file Ume/mesh_audit.cc
*/

#include "Ume/mesh_audit.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <span>
#include <vector>

namespace Ume {
namespace SOA_Idx {

namespace {

using INTV = DS_Types::INTV_T;
using INTRR = DS_Types::INTRR_T;
using SHORTV = std::vector<short>;

//! A row of `N` integers: `nkey` sort keys followed by a payload
template <size_t N> using Row = std::array<int, N>;

/*! Sort `rows` on their first `nkey` entries and hand each run of equal keys to
    `fn`, stopping early if it returns false.  This is how the audit groups
    sides by (face, edge), corners by (zone, point) and so on: a sort and a
    linear scan, rather than a node-based container. */
template <size_t N, class F>
void for_each_group(std::vector<Row<N>> &rows, int const nkey, F const &fn) {
  std::sort(rows.begin(), rows.end());
  size_t i = 0;
  while (i < rows.size()) {
    size_t j = i + 1;
    while (j < rows.size() &&
        std::equal(rows[i].cbegin(), rows[i].cbegin() + nkey, rows[j].cbegin()))
      ++j;
    if (!fn(std::span<Row<N> const>{rows.data() + i, j - i}))
      return;
    i = j;
  }
}

//! Categories used by the census
struct Census {
  int real{0}; //!< mask >= 1
  int boundary{0}; //!< mask < 0
  int ghost{0}; //!< mask == 0
};

Census tally(SHORTV const &mask, int const upto) {
  Census c;
  for (int i = 0; i < upto; ++i) {
    if (mask[i] >= 1)
      c.real += 1;
    else if (mask[i] < 0)
      c.boundary += 1;
    else
      c.ghost += 1;
  }
  return c;
}

std::ostream &operator<<(std::ostream &os, Census const &c) {
  return os << c.real << " real, " << c.boundary << " boundary, " << c.ghost
            << " ghost";
}

//! One mesh entity, for the checks that treat them all alike
struct Kind {
  char const *name;
  Entity const *entity;
};

/*! Accumulates problems for one named check, printing at most
    `Audit_Options::max_reports` of them. */
class Auditor {
public:
  Auditor(Mesh const &mesh, std::ostream &os, Audit_Options const &opts);

  Audit_Result run();

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
  //! Report one offending entity; the caller streams the detail and a newline
  std::ostream &bad() {
    hits_ += 1;
    res_.problems += 1;
    if (hits_ > opt_.max_reports)
      return null_;
    return os_ << "  [" << name_ << "] ";
  }
  //! Stop scanning once a check has produced everything it is going to print
  bool enough() const { return hits_ > opt_.max_reports; }
  //! Say that a whole group could not be run
  void skip(char const *const group, char const *const why) {
    finish();
    os_ << "  [" << group << "] skipped: " << why << '\n';
  }

  //! Is `i` a usable index into an array of `n` entries?
  static constexpr bool in(int const i, int const n) { return i >= 0 && i < n; }
  //! Can the maps of side `s` be followed?  Ghost sides carry no connectivity.
  bool usable(int const s) const {
    return smask_[s] != 0 && in(s2z_[s], zll_) && in(s2p1_[s], pll_) &&
        in(s2p2_[s], pll_) && in(s2e_[s], ell_) && in(s2f_[s], fll_) &&
        in(s2c1_[s], cll_) && in(s2c2_[s], cll_);
  }

  //! Print at most eight entries of a list
  void list(std::ostream &os, std::span<int const> v) const;
  /*! Check a derived ragged map against the (key, value) rows that define it.
      `want` is sorted in place, and is left that way for later use. */
  void compare_map(char const *const name, char const *const what,
      INTRR const &have, std::vector<Row<2>> &want, int const n);

  void ranges(bool const report);
  void connectivity();
  void iota_maps();
  void masks();
  void derived();
  void geometry();
  void parallel();
  void census();

  Mesh const &m_;
  std::ostream &os_;
  Audit_Options opt_;
  Audit_Result res_;
  char const *name_{nullptr};
  int hits_{0};
  std::ostream null_{nullptr};

  //! Entity sizes: `X_ll` counts the ghosts, `X_l` does not
  int const pll_, pl_, ell_, fll_, fl_, zll_, zl_, cll_, cl_, sll_, sl_, all_;

  INTV const &s2z_, &s2p1_, &s2p2_, &s2e_, &s2f_, &s2c1_, &s2c2_, &s2s2_,
      &s2s3_, &s2s4_, &s2s5_;
  INTV const &c2p_, &c2z_, &e2p1_, &e2p2_, &f2z1_, &f2z2_;
  INTV const &a2z_, &a2f_, &a2p_, &a2e_, &a2s_;
  SHORTV const &pmask_, &emask_, &fmask_, &zmask_, &cmask_, &smask_, &amask_;

  //! Every index a lazy derivation will use is in range
  bool ranges_ok_{true};
  //! ... including the ones the ragged-map derivations pass to `std::at()`
  bool derived_safe_{true};
};

Auditor::Auditor(Mesh const &mesh, std::ostream &os, Audit_Options const &opts)
    : m_{mesh}, os_{os}, opt_{opts}, pll_{mesh.points.size()},
      pl_{mesh.points.local_size()}, ell_{mesh.edges.size()},
      fll_{mesh.faces.size()}, fl_{mesh.faces.local_size()},
      zll_{mesh.zones.size()}, zl_{mesh.zones.local_size()},
      cll_{mesh.corners.size()}, cl_{mesh.corners.local_size()},
      sll_{mesh.sides.size()}, sl_{mesh.sides.local_size()},
      all_{mesh.iotas.size()}, s2z_{mesh.ds->caccess_intv("m:s>z")},
      s2p1_{mesh.ds->caccess_intv("m:s>p1")},
      s2p2_{mesh.ds->caccess_intv("m:s>p2")},
      s2e_{mesh.ds->caccess_intv("m:s>e")},
      s2f_{mesh.ds->caccess_intv("m:s>f")},
      s2c1_{mesh.ds->caccess_intv("m:s>c1")},
      s2c2_{mesh.ds->caccess_intv("m:s>c2")},
      s2s2_{mesh.ds->caccess_intv("m:s>s2")},
      s2s3_{mesh.ds->caccess_intv("m:s>s3")},
      s2s4_{mesh.ds->caccess_intv("m:s>s4")},
      s2s5_{mesh.ds->caccess_intv("m:s>s5")},
      c2p_{mesh.ds->caccess_intv("m:c>p")},
      c2z_{mesh.ds->caccess_intv("m:c>z")},
      e2p1_{mesh.ds->caccess_intv("m:e>p1")},
      e2p2_{mesh.ds->caccess_intv("m:e>p2")},
      f2z1_{mesh.ds->caccess_intv("m:f>z1")},
      f2z2_{mesh.ds->caccess_intv("m:f>z2")},
      a2z_{mesh.ds->caccess_intv("m:a>z")},
      a2f_{mesh.ds->caccess_intv("m:a>f")},
      a2p_{mesh.ds->caccess_intv("m:a>p")},
      a2e_{mesh.ds->caccess_intv("m:a>e")},
      a2s_{mesh.ds->caccess_intv("m:a>s")}, pmask_{mesh.points.mask},
      emask_{mesh.edges.mask}, fmask_{mesh.faces.mask}, zmask_{mesh.zones.mask},
      cmask_{mesh.corners.mask}, smask_{mesh.sides.mask},
      amask_{mesh.iotas.mask} {}

void Auditor::list(std::ostream &os, std::span<int const> v) const {
  os << '{';
  size_t const shown = std::min<size_t>(v.size(), 8);
  for (size_t i = 0; i < shown; ++i)
    os << (i ? " " : "") << v[i];
  if (shown < v.size())
    os << " ...";
  os << '}';
}

/* -------------------------------- ranges --------------------------------- */

/*! Every index in every stored map, plus the two stronger conditions the lazy
    derivations need: `VAR_zcoord` accumulates into an array of *local* zones
    indexed by `m:c>z`, and `VAR_fcoord` into an array of *local* faces indexed
    by `m:s>f`, both without bounds checking.

    Called with `report == false` when the caller turned the topology group
    off, because the `derived` and `geometry` groups have to know the answer
    either way. */
void Auditor::ranges(bool const report) {
  bool ok{true}, safe{true};

  auto open = [&](char const *const name) {
    if (report)
      start(name);
  };
  auto flag = [&]() -> std::ostream & {
    ok = false;
    return report ? bad() : null_;
  };
  auto unsafe = [&]() -> std::ostream & {
    safe = false;
    return report ? bad() : null_;
  };
  auto stop = [&]() { return report && enough(); };

  open("side map ranges");
  for (int s = 0; s < sll_ && !stop(); ++s) {
    if (smask_[s] == 0)
      continue; // ghost sides carry no local connectivity
    if (!in(s2z_[s], zll_) || !in(s2p1_[s], pll_) || !in(s2p2_[s], pll_) ||
        !in(s2e_[s], ell_) || !in(s2f_[s], fll_) || !in(s2c1_[s], cll_) ||
        !in(s2c2_[s], cll_) || !in(s2s2_[s], sll_) || !in(s2s3_[s], sll_) ||
        !in(s2s4_[s], sll_) || !in(s2s5_[s], sll_))
      flag() << "side " << s << " has an out-of-range map entry\n";
    else if (s2p1_[s] == s2p2_[s])
      flag() << "side " << s << " has p1 == p2 == " << s2p1_[s] << '\n';
  }

  open("edge map ranges");
  for (int e = 0; e < ell_ && !stop(); ++e) {
    if (emask_[e] == 0)
      continue;
    if (!in(e2p1_[e], pll_) || !in(e2p2_[e], pll_) || e2p1_[e] == e2p2_[e])
      flag() << "edge " << e << " has bad endpoints (" << e2p1_[e] << ", "
             << e2p2_[e] << ")\n";
  }

  open("corner map ranges");
  for (int c = 0; c < cll_ && !stop(); ++c) {
    if (!in(c2p_[c], pll_))
      flag() << "corner " << c << " points at point " << c2p_[c] << '\n';
    else if (cmask_[c] != 0 && !in(c2z_[c], zll_))
      flag() << "corner " << c << " (mask " << cmask_[c] << ") points at zone "
             << c2z_[c] << '\n';
  }

  /* `m:f>z2` is -1 on a face that only has a zone behind one of its sides. */
  open("face map ranges");
  for (int f = 0; f < fll_ && !stop(); ++f) {
    if (fmask_[f] == 0)
      continue;
    int const z1 = f2z1_[f], z2 = f2z2_[f];
    bool const has1 = in(z1, zll_), has2 = in(z2, zll_);
    if ((z1 >= 0 && !has1) || (z2 >= 0 && !has2))
      flag() << "face " << f << " names zones (" << z1 << ", " << z2
             << "), out of [0, " << zll_ << ")\n";
    else if (!has1 && !has2)
      flag() << "face " << f << " (mask " << fmask_[f] << ") has no zone\n";
    else if (has1 && has2 && z1 == z2)
      flag() << "face " << f << " has zone " << z1 << " on both sides\n";
  }

  if (all_ > 0) {
    open("iota map ranges");
    for (int a = 0; a < all_ && !stop(); ++a) {
      if (amask_[a] == 0)
        continue;
      if (!in(a2z_[a], zll_) || !in(a2f_[a], fll_) || !in(a2p_[a], pll_) ||
          !in(a2e_[a], ell_) || !in(a2s_[a], sll_))
        flag() << "iota " << a << " has an out-of-range map entry\n";
    }
  }

  /*! `VAR_zcoord` counts the points of each *local* zone in an array of length
      `zones.local_size()`, indexed by `m:c>z` over the local corners. */
  open("local corners name local zones");
  for (int c = 0; c < cl_ && !stop(); ++c)
    if (cmask_[c] != 0 && in(c2z_[c], zll_) && c2z_[c] >= zl_)
      flag() << "local corner " << c << " is in ghost zone " << c2z_[c]
             << " (VAR_zcoord writes past its local zone array)\n";

  //! `VAR_fcoord` does the same over the local faces, indexed by `m:s>f`.
  open("local sides name local faces");
  for (int s = 0; s < sl_ && !stop(); ++s)
    if (smask_[s] != 0 && in(s2f_[s], fll_) && s2f_[s] >= fl_)
      flag() << "local side " << s << " is on ghost face " << s2f_[s]
             << " (VAR_fcoord writes past its local face array)\n";

  /*! The ragged-map derivations index accumulators with `std::vector::at`, and
      they do it for every side and every corner, ghosts included: a negative
      index there throws instead of being ignored. */
  open("derived maps are index-safe");
  for (int s = 0; s < sll_ && !stop(); ++s)
    if (!in(s2c1_[s], cll_) || !in(s2c2_[s], cll_))
      unsafe() << "side " << s << " (mask " << smask_[s] << ") has corners ("
               << s2c1_[s] << ", " << s2c2_[s] << "), out of [0, " << cll_
               << "): m:c>s cannot be built\n";
  for (int c = 0; c < cll_ && !stop(); ++c)
    if (c2p_[c] < 0 || c2z_[c] < 0)
      unsafe() << "corner " << c << " (mask " << cmask_[c] << ") names ("
               << "point " << c2p_[c] << ", zone " << c2z_[c]
               << "): m:z>c, m:z>p, m:p>zs and m:z>pz cannot be built\n";

  if (report)
    finish();
  ranges_ok_ = ok;
  derived_safe_ = ok && safe;
}

/* ----------------------------- connectivity ------------------------------ */

void Auditor::connectivity() {
  std::vector<char> ok(static_cast<size_t>(sll_), 0);
  for (int s = 0; s < sll_; ++s)
    ok[s] = static_cast<char>(usable(s));

  /* A ghost corner's `m:c>z` is allowed to be garbage, but only in the
     recognizable way the library already guards against: out of range, or
     naming a zone that is itself a ghost. */
  start("ghost corners do not name real zones");
  for (int c = 0; c < cll_ && !enough(); ++c)
    if (cmask_[c] == 0 && in(c2z_[c], zll_) && zmask_[c2z_[c]] != 0)
      bad() << "ghost corner " << c << " names real zone " << c2z_[c] << '\n';

  /* `m:s>p1`/`p2` are a denormalization of the side's edge, and every operator
     that uses one assumes it agrees with the other. */
  start("m:s>e agrees with m:s>p1/p2");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s])
      continue;
    int const e = s2e_[s], p1 = s2p1_[s], p2 = s2p2_[s];
    if (!in(e2p1_[e], pll_) || !in(e2p2_[e], pll_))
      continue; // already reported
    if (!((e2p1_[e] == p1 && e2p2_[e] == p2) ||
            (e2p1_[e] == p2 && e2p2_[e] == p1)))
      bad() << "side " << s << " has points (" << p1 << ", " << p2
            << ") but edge " << e << " has (" << e2p1_[e] << ", " << e2p2_[e]
            << ")\n";
  }

  start("m:s>s2 pairing");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s] || !in(s2s2_[s], sll_))
      continue;
    int const s2 = s2s2_[s];
    if (s2 == s)
      bad() << "side " << s << " is its own s2\n";
    else if (s2s2_[s2] != s)
      bad() << "side " << s << " -> " << s2 << " -> " << s2s2_[s2]
            << " is not an involution\n";
    else if (s2f_[s2] != s2f_[s] || s2e_[s2] != s2e_[s])
      bad() << "side " << s << " and " << s2 << " share s2 but not their "
            << "(face, edge): (" << s2f_[s] << ", " << s2e_[s] << ") vs ("
            << s2f_[s2] << ", " << s2e_[s2] << ")\n";
  }

  /*! The two sides on a facet traverse their shared edge in opposite
      directions, so that their `side_surz` cancel across the facet.  Which of
      the two orders that is is not fixed by the file format, so take the
      convention from the mesh and report the sides that break it. */
  {
    int opposed{0}, aligned{0};
    auto paired = [&](int const s) {
      return ok[s] && in(s2s2_[s], sll_) && ok[s2s2_[s]] && s2s2_[s] != s;
    };
    for (int s = 0; s < sll_; ++s)
      if (paired(s)) {
        if (s2p1_[s2s2_[s]] == s2p2_[s])
          opposed += 1;
        else if (s2p1_[s2s2_[s]] == s2p1_[s])
          aligned += 1;
      }
    bool const expect_opposed = opposed >= aligned;
    start("m:s>s2 edge direction");
    for (int s = 0; s < sll_ && !enough(); ++s) {
      if (!paired(s))
        continue;
      int const s2 = s2s2_[s];
      bool const is_opposed = (s2p1_[s2] == s2p2_[s] && s2p2_[s2] == s2p1_[s]);
      bool const is_aligned = (s2p1_[s2] == s2p1_[s] && s2p2_[s2] == s2p2_[s]);
      if (expect_opposed ? !is_opposed : !is_aligned)
        bad() << "sides " << s << " (" << s2p1_[s] << ", " << s2p2_[s]
              << ") and " << s2 << " (" << s2p1_[s2] << ", " << s2p2_[s2]
              << ") do not traverse their shared edge "
              << (expect_opposed ? "in opposite directions" : "alike")
              << ", as the rest of the mesh does\n";
    }
  }

  /* Every (face, edge) is shared by exactly two sides: one per zone across
     the face, or the real side and its phantom on the mesh boundary. */
  start("(face, edge) side counts");
  {
    std::vector<Row<3>> fe;
    fe.reserve(static_cast<size_t>(sll_));
    for (int s = 0; s < sll_; ++s)
      if (ok[s])
        fe.push_back({s2f_[s], s2e_[s], s});
    for_each_group(fe, 2, [&](std::span<Row<3> const> g) {
      if (enough())
        return false;
      if (g.size() != 2)
        bad() << "face " << g[0][0] << ", edge " << g[0][1] << " has "
              << g.size() << " sides\n";
      return true;
    });
  }

  start("m:s>c1/c2 agree with (zone, p1/p2)");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s])
      continue;
    int const c1 = s2c1_[s], c2 = s2c2_[s];
    if (c2z_[c1] != s2z_[s] || c2p_[c1] != s2p1_[s])
      bad() << "side " << s << " c1 " << c1 << " is (zone " << c2z_[c1]
            << ", point " << c2p_[c1] << "), expected (zone " << s2z_[s]
            << ", point " << s2p1_[s] << ")\n";
    else if (c2z_[c2] != s2z_[s] || c2p_[c2] != s2p2_[s])
      bad() << "side " << s << " c2 " << c2 << " is (zone " << c2z_[c2]
            << ", point " << c2p_[c2] << "), expected (zone " << s2z_[s]
            << ", point " << s2p2_[s] << ")\n";
  }

  start("no duplicate corners in a zone");
  {
    std::vector<Row<3>> zp;
    zp.reserve(static_cast<size_t>(cl_));
    for (int c = 0; c < cl_; ++c)
      if (cmask_[c] != 0 && in(c2z_[c], zll_) && in(c2p_[c], pll_))
        zp.push_back({c2z_[c], c2p_[c], c});
    for_each_group(zp, 2, [&](std::span<Row<3> const> g) {
      if (enough())
        return false;
      if (g.size() > 1)
        bad() << "corners " << g[0][2] << " and " << g[1][2]
              << " are both (zone " << g[0][0] << ", point " << g[0][1]
              << ")\n";
      return true;
    });
  }

  /*! `m:s>s3` and `m:s>s4` walk the ring of sides around one face of one zone.
      A boundary face carries two rings when the mesh has no boundary zones:
      the real zone's, and the phantom one that shadows it, which is why the
      rings are grouped by mask class as well as by (zone, face). */
  auto ring_class = [this](int const s) { return smask_[s] > 0 ? 1 : 0; };

  start("m:s>s3/s4 inverse");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s] || !in(s2s3_[s], sll_) || !in(s2s4_[s], sll_))
      continue;
    int const nxt = s2s3_[s], prv = s2s4_[s];
    if (!in(s2s4_[nxt], sll_) || s2s4_[nxt] != s || s2s3_[prv] != s)
      bad() << "side " << s << " has s3 " << nxt << " and s4 " << prv
            << ", but s4(s3) is " << s2s4_[nxt] << " and s3(s4) is "
            << s2s3_[prv] << '\n';
    else if (s2z_[nxt] != s2z_[s] || s2f_[nxt] != s2f_[s])
      bad() << "side " << s << " and its s3 " << nxt
            << " are not on the same (zone, face): (" << s2z_[s] << ", "
            << s2f_[s] << ") vs (" << s2z_[nxt] << ", " << s2f_[nxt] << ")\n";
    else if (ring_class(nxt) != ring_class(s))
      bad() << "side " << s << " (mask " << smask_[s] << ") has s3 " << nxt
            << " (mask " << smask_[nxt] << "): a real side and a phantom side "
            << "are not on the same ring\n";
  }

  start("m:s>s3 face rings");
  {
    std::vector<Row<4>> zf;
    zf.reserve(static_cast<size_t>(sll_));
    for (int s = 0; s < sll_; ++s)
      if (ok[s])
        zf.push_back({s2z_[s], s2f_[s], ring_class(s), s});
    std::vector<int> orbit, want;
    for_each_group(zf, 3, [&](std::span<Row<4> const> g) {
      if (enough())
        return false;
      int const first = g[0][3];
      orbit.clear();
      int s = first;
      bool closed{false};
      while (orbit.size() <= g.size()) {
        orbit.push_back(s);
        if (!in(s2s3_[s], sll_))
          break;
        s = s2s3_[s];
        if (s == first) {
          closed = true;
          break;
        }
      }
      want.clear();
      for (auto const &r : g)
        want.push_back(r[3]);
      std::sort(orbit.begin(), orbit.end());
      if (!closed || orbit != want) {
        bad() << "zone " << g[0][0] << ", face " << g[0][1] << " has "
              << g.size() << " sides ";
        list(os_, want);
        os_ << ", but the m:s>s3 ring through side " << first << " is ";
        list(os_, orbit);
        os_ << (closed ? "\n" : " and does not close\n");
      } else if (g.size() < 3)
        bad() << "zone " << g[0][0] << ", face " << g[0][1] << " has only "
              << g.size() << " sides\n";
      return true;
    });
  }

  /*! Which way the ring turns is again a convention, but a mesh has to pick
      one: the ring is a directed cycle of the face's points. */
  {
    int forward{0}, backward{0};
    auto ringed = [&](int const s) {
      return ok[s] && in(s2s3_[s], sll_) && ok[s2s3_[s]] && s2s3_[s] != s;
    };
    for (int s = 0; s < sll_; ++s)
      if (ringed(s)) {
        forward += (s2p1_[s2s3_[s]] == s2p2_[s]) ? 1 : 0;
        backward += (s2p2_[s2s3_[s]] == s2p1_[s]) ? 1 : 0;
      }
    bool const expect_forward = forward >= backward;
    start("m:s>s3 ring direction");
    for (int s = 0; s < sll_ && !enough(); ++s) {
      if (!ringed(s))
        continue;
      int const nxt = s2s3_[s];
      bool const match =
          expect_forward ? (s2p1_[nxt] == s2p2_[s]) : (s2p2_[nxt] == s2p1_[s]);
      if (!match)
        bad() << "side " << s << " (" << s2p1_[s] << ", " << s2p2_[s]
              << ") does not meet its s3 " << nxt << " (" << s2p1_[nxt] << ", "
              << s2p2_[nxt] << ") at the point the rest of the mesh would\n";
    }
  }

  /*! `m:s>s5` is the other side of the same body on the same edge: the one on
      the second of the two faces that meet there.  A body with a single face
      -- a boundary zone, or the phantom ring shadowing one boundary face --
      has no second face, and its sides refer to themselves. */
  start("m:s>s5 pairing");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s] || !in(s2s5_[s], sll_))
      continue;
    int const s5 = s2s5_[s];
    if (s5 == s)
      continue; // a one-faced body; the grouping below is what catches these
    if (s2s5_[s5] != s)
      bad() << "side " << s << " -> " << s5 << " -> " << s2s5_[s5]
            << " is not an involution\n";
    else if (s2z_[s5] != s2z_[s] || s2e_[s5] != s2e_[s])
      bad() << "side " << s << " and " << s5 << " share s5 but not their "
            << "(zone, edge): (" << s2z_[s] << ", " << s2e_[s] << ") vs ("
            << s2z_[s5] << ", " << s2e_[s5] << ")\n";
    else if (s2f_[s5] == s2f_[s])
      bad() << "sides " << s << " and " << s5 << " are s5 partners on face "
            << s2f_[s] << ": the two must be on the faces that meet at edge "
            << s2e_[s] << '\n';
  }

  /* Within a real zone the statement is complete: a closed polyhedron meets
     each of its edges on exactly two of its faces.  (Phantom sides are left
     out: two boundary faces of one zone can meet at an edge, and their
     phantom rings are separate one-faced bodies that share the zone index.) */
  start("real zones meet each edge twice");
  {
    std::vector<Row<3>> ze;
    ze.reserve(static_cast<size_t>(sll_));
    for (int s = 0; s < sll_; ++s)
      if (ok[s] && smask_[s] > 0)
        ze.push_back({s2z_[s], s2e_[s], s});
    for_each_group(ze, 2, [&](std::span<Row<3> const> g) {
      if (enough())
        return false;
      if (g.size() != 2)
        bad() << "zone " << g[0][0] << " has " << g.size() << " real sides on "
              << "edge " << g[0][1] << '\n';
      else if (s2s5_[g[0][2]] != g[1][2] || s2s5_[g[1][2]] != g[0][2])
        bad() << "sides " << g[0][2] << " and " << g[1][2]
              << " are the sides of zone " << g[0][0] << " on edge " << g[0][1]
              << ", but their s5 are (" << s2s5_[g[0][2]] << ", "
              << s2s5_[g[1][2]] << ")\n";
      return true;
    });
  }

  start("sides agree with m:f>z1/z2");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!ok[s])
      continue;
    int const f = s2f_[s], z = s2z_[s];
    if (fmask_[f] == 0)
      continue; // a ghost face keeps no zone links
    if (f2z1_[f] != z && f2z2_[f] != z)
      bad() << "side " << s << " is in zone " << z << " on face " << f
            << ", which names zones (" << f2z1_[f] << ", " << f2z2_[f] << ")\n";
  }

  /* The smallest closed polyhedron is a tetrahedron. */
  start("real zones are closed");
  {
    std::vector<int> ncorners(static_cast<size_t>(zll_), 0);
    for (int c = 0; c < cl_; ++c)
      if (cmask_[c] != 0 && in(c2z_[c], zll_))
        ncorners[c2z_[c]] += 1;
    std::vector<Row<2>> zf;
    zf.reserve(static_cast<size_t>(sll_));
    for (int s = 0; s < sll_; ++s)
      if (ok[s] && smask_[s] > 0)
        zf.push_back({s2z_[s], s2f_[s]});
    std::sort(zf.begin(), zf.end());
    zf.erase(std::unique(zf.begin(), zf.end()), zf.end());
    std::vector<int> nfaces(static_cast<size_t>(zll_), 0);
    for (auto const &r : zf)
      nfaces[r[0]] += 1;
    for (int z = 0; z < zll_ && !enough(); ++z) {
      if (zmask_[z] < 1)
        continue;
      if (ncorners[z] < 4 || nfaces[z] < 4)
        bad() << "zone " << z << " has " << ncorners[z] << " corners on "
              << nfaces[z] << " faces\n";
    }
  }
}

/* -------------------------------- iotas ---------------------------------- */

/*! An iota is one endpoint of one side: the tetrahedron `(zone center, face
    center, edge center, point)`.  Every one of its maps is a copy of the
    corresponding map of its side. */
void Auditor::iota_maps() {
  start("m:a> maps agree with m:a>s");
  for (int a = 0; a < all_ && !enough(); ++a) {
    if (amask_[a] == 0 || !in(a2s_[a], sll_) || !in(a2p_[a], pll_))
      continue;
    int const s = a2s_[a];
    if (!usable(s))
      continue;
    if (a2z_[a] != s2z_[s] || a2f_[a] != s2f_[s] || a2e_[a] != s2e_[s])
      bad() << "iota " << a << " names (zone " << a2z_[a] << ", face "
            << a2f_[a] << ", edge " << a2e_[a] << "), but its side " << s
            << " names (zone " << s2z_[s] << ", face " << s2f_[s] << ", edge "
            << s2e_[s] << ")\n";
    else if (a2p_[a] != s2p1_[s] && a2p_[a] != s2p2_[s])
      bad() << "iota " << a << " is at point " << a2p_[a] << ", but its side "
            << s << " runs (" << s2p1_[s] << ", " << s2p2_[s] << ")\n";
  }

  start("iotas per side");
  {
    std::vector<Row<3>> sa;
    sa.reserve(static_cast<size_t>(all_));
    for (int a = 0; a < all_; ++a)
      if (amask_[a] != 0 && in(a2s_[a], sll_))
        sa.push_back({a2s_[a], a2p_[a], a});
    for_each_group(sa, 1, [&](std::span<Row<3> const> g) {
      if (enough())
        return false;
      if (g.size() > 2)
        bad() << "side " << g[0][0] << " carries " << g.size()
              << " iotas; a side has one per endpoint\n";
      else if (g.size() == 2 && g[0][1] == g[1][1])
        bad() << "iotas " << g[0][2] << " and " << g[1][2] << " are both at "
              << "point " << g[0][1] << " of side " << g[0][0] << '\n';
      return true;
    });
  }
}

/* -------------------------------- masks ---------------------------------- */

void Auditor::masks() {
  start("real corners belong to real zones");
  for (int c = 0; c < cl_ && !enough(); ++c)
    if (cmask_[c] >= 1 && in(c2z_[c], zll_) && zmask_[c2z_[c]] < 1)
      bad() << "corner " << c << " (mask " << cmask_[c] << ") is in zone "
            << c2z_[c] << " (mask " << zmask_[c2z_[c]] << ")\n";

  start("real sides belong to real zones");
  for (int s = 0; s < sll_ && !enough(); ++s)
    if (smask_[s] > 0 && in(s2z_[s], zll_) && zmask_[s2z_[s]] < 1)
      bad() << "side " << s << " (mask " << smask_[s] << ") is in zone "
            << s2z_[s] << " (mask " << zmask_[s2z_[s]] << ")\n";

  /* A side that is not itself a ghost is a piece of a local zone, so nothing
     it names can be a ghost either. */
  start("sides do not reference ghosts");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (!usable(s))
      continue;
    char const *what{nullptr};
    if (zmask_[s2z_[s]] == 0)
      what = "zone";
    else if (fmask_[s2f_[s]] == 0)
      what = "face";
    else if (emask_[s2e_[s]] == 0)
      what = "edge";
    else if (pmask_[s2p1_[s]] == 0 || pmask_[s2p2_[s]] == 0)
      what = "point";
    else if (cmask_[s2c1_[s]] == 0 || cmask_[s2c2_[s]] == 0)
      what = "corner";
    if (what)
      bad() << "side " << s << " (mask " << smask_[s] << ") references a ghost "
            << what << '\n';
  }

  /*! `gradzatp` accumulates `point_volume` only over `cmask >= 1` corners and
      then divides every `pmask != 0` point by it, so a point without such a
      corner produces 0/0 and poisons the gradient. */
  start("every real point owns a cmask >= 1 corner");
  {
    std::vector<char> owned(static_cast<size_t>(pll_), 0);
    for (int c = 0; c < cl_; ++c)
      if (cmask_[c] >= 1 && in(c2p_[c], pll_))
        owned[c2p_[c]] = 1;
    for (int p = 0; p < pl_ && !enough(); ++p)
      if (pmask_[p] != 0 && !owned[p])
        bad() << "point " << p << " (mask " << pmask_[p]
              << ") has no cmask >= 1 corner: gradzatp will divide by zero\n";
  }

  //! `VAR_zcoord` divides each zone by the number of corners that reached it.
  start("every non-ghost zone owns a non-ghost corner");
  {
    std::vector<char> owned(static_cast<size_t>(zll_), 0);
    for (int c = 0; c < cl_; ++c)
      if (cmask_[c] != 0 && in(c2z_[c], zll_))
        owned[c2z_[c]] = 1;
    for (int z = 0; z < zl_ && !enough(); ++z)
      if (zmask_[z] != 0 && !owned[z])
        bad() << "zone " << z << " (mask " << zmask_[z]
              << ") owns no corner: zcoord will divide by zero\n";
  }

  //! ... and `VAR_fcoord` by the number of sides that reached each face.
  start("every non-ghost face owns a non-ghost side");
  {
    std::vector<char> owned(static_cast<size_t>(fll_), 0);
    for (int s = 0; s < sl_; ++s)
      if (smask_[s] != 0 && in(s2f_[s], fll_))
        owned[s2f_[s]] = 1;
    for (int f = 0; f < fl_ && !enough(); ++f)
      if (fmask_[f] != 0 && !owned[f])
        bad() << "face " << f << " (mask " << fmask_[f]
              << ") owns no side: fcoord will divide by zero\n";
  }

  /*! The library takes "ghost" and "in the upper part of the index range" to
      mean the same thing: iteration bounds are `local_size()`, but the mask is
      what the kernels branch on. */
  start("ghosts occupy the trailing indices");
  for (auto const &k : {Kind{"point", &m_.points}, Kind{"edge", &m_.edges},
           Kind{"face", &m_.faces}, Kind{"zone", &m_.zones},
           Kind{"corner", &m_.corners}, Kind{"side", &m_.sides},
           Kind{"iota", &m_.iotas}}) {
    int const n = k.entity->size(), l = k.entity->local_size();
    if (l < 0 || l > n) {
      bad() << k.name << "s: local_size " << l << " is not within [0, " << n
            << "]\n";
      continue;
    }
    for (int i = 0; i < n && !enough(); ++i) {
      if (i >= l && k.entity->mask[i] != 0)
        bad() << k.name << ' ' << i << " is a ghost (>= local_size " << l
              << ") but has mask " << k.entity->mask[i] << '\n';
      else if (i < l && k.entity->mask[i] == 0)
        bad() << k.name << ' ' << i << " has mask 0 but is below local_size "
              << l << '\n';
    }
  }

  /*! A side outside the real mesh is on the problem boundary, and so is
      everything it touches: this is the marking `point_norm` and the boundary
      term of `gradzatp` read. */
  start("the mesh boundary is consistently marked");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (smask_[s] >= 0 || !usable(s))
      continue;
    char const *what{nullptr};
    if (fmask_[s2f_[s]] >= 1)
      what = "face";
    else if (emask_[s2e_[s]] >= 1)
      what = "edge";
    else if (pmask_[s2p1_[s]] >= 1 || pmask_[s2p2_[s]] >= 1)
      what = "point";
    if (what)
      bad() << "side " << s << " (mask " << smask_[s]
            << ") lies outside the real mesh, but its " << what
            << " is marked interior\n";
  }

  /*! `VAR_point_norm` reaches through a boundary side to the real side behind
      it, and takes that side's `side_surz` as the outward facet. */
  start("phantom sides pair with real sides");
  for (int s = 0; s < sll_ && !enough(); ++s) {
    if (smask_[s] >= 0 || !in(s2s2_[s], sll_))
      continue;
    if (smask_[s2s2_[s]] <= 0)
      bad() << "side " << s << " (mask " << smask_[s] << ") has s2 " << s2s2_[s]
            << " (mask " << smask_[s2s2_[s]]
            << "): point_norm needs a real side there\n";
  }
}

/* ------------------------------- derived --------------------------------- */

void Auditor::compare_map(char const *const name, char const *const what,
    INTRR const &have, std::vector<Row<2>> &want, int const n) {
  start(name);
  std::sort(want.begin(), want.end());
  size_t i = 0;
  std::vector<int> got, expect;
  for (int k = 0; k < n; ++k) {
    size_t j = i;
    while (j < want.size() && want[j][0] == k)
      ++j;
    if (enough())
      return;
    auto const entry = have[k];
    got.assign(entry.begin(), entry.end());
    std::sort(got.begin(), got.end());
    bool same = got.size() == j - i;
    for (size_t t = 0; same && t < got.size(); ++t)
      same = (got[t] == want[i + t][1]);
    if (!same) {
      expect.clear();
      for (size_t t = i; t < j; ++t)
        expect.push_back(want[t][1]);
      bad() << what << ' ' << k << " has ";
      list(os_, got);
      os_ << ", expected ";
      list(os_, expect);
      os_ << '\n';
    }
    i = j;
  }
}

/*! The maps the library derives from the stored ones, each against the
    definition its `VAR_*::init_()` uses.  A derived map is cached, and the
    operators use it in place of the stored maps it came from, so it is worth
    confirming that the two still say the same thing. */
void Auditor::derived() {
  std::vector<Row<2>> want;

  /* Every side contributes its two corners, ghost sides included. */
  want.clear();
  want.reserve(2 * static_cast<size_t>(sll_));
  for (int s = 0; s < sll_; ++s) {
    want.push_back({s2c1_[s], s});
    want.push_back({s2c2_[s], s});
  }
  compare_map("m:c>s inverts m:s>c1/c2", "corner",
      m_.ds->caccess_intrr("m:c>s"), want, cll_);

  want.clear();
  want.reserve(static_cast<size_t>(cll_));
  for (int c = 0; c < cll_; ++c)
    if (c2z_[c] < zll_)
      want.push_back({c2z_[c], c});
  compare_map(
      "m:z>c inverts m:c>z", "zone", m_.ds->caccess_intrr("m:z>c"), want, zll_);
  std::vector<Row<2>> const zone_corners{want};

  want.clear();
  want.reserve(static_cast<size_t>(cll_));
  for (int c = 0; c < cll_; ++c)
    if (c2p_[c] < pll_ && c2z_[c] < zll_)
      want.push_back({c2z_[c], c2p_[c]});
  compare_map("m:z>p lists the zone's points", "zone",
      m_.ds->caccess_intrr("m:z>p"), want, zll_);

  want.clear();
  want.reserve(static_cast<size_t>(cll_));
  for (int c = 0; c < cll_; ++c)
    if (c2p_[c] < pll_ && c2z_[c] < zll_)
      want.push_back({c2p_[c], c2z_[c]});
  compare_map("m:p>zs inverts (m:c>p, m:c>z)", "point",
      m_.ds->caccess_intrr("m:p>zs"), want, pll_);
  std::vector<Row<2>> const point_zones{want}; // sorted by (point, zone)

  want.clear();
  want.reserve(static_cast<size_t>(cl_));
  for (int c = 0; c < cl_; ++c)
    if (cmask_[c] >= 1 && c2p_[c] < pll_)
      want.push_back({c2p_[c], c});
  compare_map("m:p>rc lists the real corners at a point", "point",
      m_.ds->caccess_intrr("m:p>rc"), want, pll_);

  /*! `m:z>pz` is the set of zones that share a point with a zone, which is one
      hop further than the rest: walk it through the (zone, corner) and (point,
      zone) rows built above, marking each neighbor once. */
  start("m:z>pz lists the point-sharing zones");
  {
    auto offsets = [](std::vector<Row<2>> const &rows, int const n) {
      std::vector<size_t> off(static_cast<size_t>(n) + 1, 0);
      for (auto const &r : rows)
        off[static_cast<size_t>(r[0]) + 1] += 1;
      for (int k = 0; k < n; ++k)
        off[static_cast<size_t>(k) + 1] += off[static_cast<size_t>(k)];
      return off;
    };
    std::vector<size_t> const zoff = offsets(zone_corners, zll_);
    std::vector<size_t> const poff = offsets(point_zones, pll_);
    auto const &z2pz = m_.ds->caccess_intrr("m:z>pz");
    std::vector<int> stamp(static_cast<size_t>(zll_), -1);
    for (int z = 0; z < zll_ && !enough(); ++z) {
      int n{0};
      for (size_t i = zoff[z]; i < zoff[z + 1]; ++i) {
        int const p = c2p_[zone_corners[i][1]];
        if (p >= pll_)
          continue;
        for (size_t t = poff[p]; t < poff[p + 1]; ++t) {
          int const other = point_zones[t][1];
          if (other != z && stamp[other] != z) {
            stamp[other] = z;
            n += 1;
          }
        }
      }
      auto const have = z2pz[z];
      bool same = static_cast<int>(have.size()) == n;
      for (size_t t = 0; same && t < have.size(); ++t)
        same = in(have[t], zll_) && stamp[have[t]] == z &&
            (t == 0 || have[t - 1] < have[t]);
      if (!same) {
        bad() << "zone " << z << " has ";
        list(os_, have);
        os_ << ", expected " << n << " distinct point-sharing zones\n";
      }
    }
  }
}

/* ------------------------------- geometry -------------------------------- */

void Auditor::geometry() {
  auto const &pcoord = m_.ds->caccess_vec3v("pcoord");
  auto const &side_vol = m_.ds->caccess_dblv("side_vol");
  auto const &side_surz = m_.ds->caccess_vec3v("side_surz");
  auto const &corner_vol = m_.ds->caccess_dblv("corner_vol");
  auto const &csurf = m_.ds->caccess_vec3v("corner_csurf");
  auto const &point_norm = m_.ds->caccess_vec3v("point_norm");

  start("pcoord is finite");
  for (int p = 0; p < pll_ && !enough(); ++p) {
    if (pmask_[p] == 0)
      continue;
    if (!(std::isfinite(pcoord[p][0]) && std::isfinite(pcoord[p][1]) &&
            std::isfinite(pcoord[p][2])))
      bad() << "point " << p << " is at " << pcoord[p] << '\n';
  }

  start("edges have positive length");
  for (int e = 0; e < ell_ && !enough(); ++e) {
    if (emask_[e] == 0 || !in(e2p1_[e], pll_) || !in(e2p2_[e], pll_))
      continue;
    if (!(vectormag(pcoord[e2p1_[e]] - pcoord[e2p2_[e]]) > 0.0))
      bad() << "edge " << e << " has zero length at " << pcoord[e2p1_[e]]
            << '\n';
  }

  start("side_vol > 0 on real sides");
  for (int s = 0; s < sl_ && !enough(); ++s)
    if (smask_[s] > 0 && !(side_vol[s] > 0.0))
      bad() << "side " << s << " has volume " << side_vol[s]
            << " (check the p1/p2 orientation convention)\n";

  start("corner_vol > 0 on real corners");
  for (int c = 0; c < cl_ && !enough(); ++c)
    if (cmask_[c] >= 1 && !(corner_vol[c] > 0.0))
      bad() << "corner " << c << " has volume " << corner_vol[c] << '\n';

  /* Accumulate per zone straight from `m:c>z` rather than through `m:z>c`, so
     that this stays valid whatever the zone-to-corner map is filtered on. */
  std::vector<double> cvol(static_cast<size_t>(zll_), 0.0);
  std::vector<Vec3> zsurf(static_cast<size_t>(zll_), Vec3(0.0));
  std::vector<double> zscale(static_cast<size_t>(zll_), 0.0);
  for (int c = 0; c < cl_; ++c) {
    if (cmask_[c] == 0 || !in(c2z_[c], zll_))
      continue;
    int const z = c2z_[c];
    cvol[z] += corner_vol[c];
    zsurf[z] += csurf[c];
    zscale[z] += vectormag(csurf[c]);
  }

  /* The same quantities gathered over the real sides of each zone. */
  std::vector<double> svol(static_cast<size_t>(zll_), 0.0);
  std::vector<Vec3> zsurz(static_cast<size_t>(zll_), Vec3(0.0));
  std::vector<double> zscalz(static_cast<size_t>(zll_), 0.0);
  for (int s = 0; s < sl_; ++s) {
    if (smask_[s] <= 0 || !in(s2z_[s], zll_))
      continue;
    int const z = s2z_[s];
    svol[z] += side_vol[s];
    zsurz[z] += side_surz[s];
    zscalz[z] += vectormag(side_surz[s]);
  }

  start("real zones have positive volume");
  for (int z = 0; z < zll_ && !enough(); ++z)
    if (zmask_[z] >= 1 && !(cvol[z] > 0.0))
      bad() << "zone " << z << " has volume " << cvol[z] << '\n';

  /*! `corner_vol` is the side volumes of the zone, split between the two
      corners each side meets, so the two ways of adding them up agree.  They
      stop agreeing if `m:s>c1/c2` and `m:s>z` disagree about which zone a side
      belongs to. */
  start("zone volume agrees between sides and corners");
  for (int z = 0; z < zll_ && !enough(); ++z) {
    if (zmask_[z] < 1)
      continue;
    if (std::abs(cvol[z] - svol[z]) > 1e-12 * std::abs(svol[z]))
      bad() << "zone " << z << " has volume " << cvol[z]
            << " through its corners and " << svol[z] << " through its sides\n";
  }

  /*! The corner surfaces of a zone are its internal facets, each counted once
      with either sign, so they have to cancel.  If they do not, `gradzatp`
      does not reproduce a constant field. */
  start("corner_csurf closes over each zone");
  for (int z = 0; z < zll_ && !enough(); ++z) {
    if (zmask_[z] == 0)
      continue;
    if (vectormag(zsurf[z]) > 1e-10 * zscale[z])
      bad() << "zone " << z << " corner surfaces sum to " << zsurf[z]
            << " (facet area scale " << zscale[z] << ")\n";
  }

  /*! `side_surz` is the outward facet of the zone's own boundary, so those
      have to close too -- over the surface of the zone rather than its
      interior.  A face ring traversed the wrong way about one of the zones
      shows up here. */
  start("side_surz closes over each zone");
  for (int z = 0; z < zll_ && !enough(); ++z) {
    if (zmask_[z] < 1)
      continue;
    if (vectormag(zsurz[z]) > 1e-10 * zscalz[z])
      bad() << "zone " << z << " side surfaces sum to " << zsurz[z]
            << " (facet area scale " << zscalz[z] << ")\n";
  }

  start("point_norm is a unit vector on boundary points");
  for (int p = 0; p < pl_ && !enough(); ++p) {
    double const mag = vectormag(point_norm[p]);
    if (pmask_[p] < 0 && std::abs(mag - 1.0) > 1e-12)
      bad() << "point " << p << " (mask " << pmask_[p] << ") has |point_norm| "
            << mag << '\n';
    else if (pmask_[p] >= 1 && mag > 1e-12)
      bad() << "interior point " << p << " has a nonzero point_norm "
            << point_norm[p] << '\n';
  }
}

/* ------------------------------- parallel -------------------------------- */

/*! The minimal parallel connectivity each `Entity` carries.  Only one
    partition is present, so what can be checked is that the records are
    self-consistent and point at entities and ranks that exist. */
void Auditor::parallel() {
  std::array<Kind, 7> const kinds{{{"point", &m_.points}, {"edge", &m_.edges},
      {"face", &m_.faces}, {"zone", &m_.zones}, {"corner", &m_.corners},
      {"side", &m_.sides}, {"iota", &m_.iotas}}};

  start("communication array lengths");
  for (auto const &k : kinds) {
    Entity const &e = *k.entity;
    size_t const n = static_cast<size_t>(e.size());
    if (e.comm_type.size() != n)
      bad() << k.name << "s: " << e.comm_type.size() << " comm types for " << n
            << " entities\n";
    if (e.src_pe.size() != e.cpy_idx.size() ||
        e.src_idx.size() != e.cpy_idx.size() ||
        e.ghost_mask.size() != e.cpy_idx.size())
      bad() << k.name << "s: copy records disagree in length (cpy_idx "
            << e.cpy_idx.size() << ", src_pe " << e.src_pe.size()
            << ", src_idx " << e.src_idx.size() << ", ghost_mask "
            << e.ghost_mask.size() << ")\n";
  }

  start("comm types");
  for (auto const &k : kinds) {
    Entity const &e = *k.entity;
    int const n = std::min(e.size(), static_cast<int>(e.comm_type.size()));
    for (int i = 0; i < n && !enough(); ++i) {
      int const t = e.comm_type[i];
      if (t < Entity::INTERNAL || t > Entity::GHOST)
        bad() << k.name << ' ' << i << " has comm type " << t << '\n';
      else if (i >= e.local_size() && t != Entity::GHOST)
        bad() << k.name << ' ' << i << " is beyond local_size "
              << e.local_size() << " but its comm type is " << t << '\n';
    }
  }

  start("copy records");
  for (auto const &k : kinds) {
    Entity const &e = *k.entity;
    size_t const n = std::min({e.cpy_idx.size(), e.src_pe.size(),
        e.src_idx.size(), e.ghost_mask.size()});
    for (size_t i = 0; i < n && !enough(); ++i) {
      if (!in(e.cpy_idx[i], e.size()))
        bad() << k.name << " copy " << i << " is local index " << e.cpy_idx[i]
              << ", out of [0, " << e.size() << ")\n";
      else if (!in(e.src_pe[i], m_.numpe) || e.src_pe[i] == m_.mype)
        bad() << k.name << " copy " << i << " has source rank " << e.src_pe[i]
              << " (this is rank " << m_.mype << " of " << m_.numpe << ")\n";
      else if (e.src_idx[i] < 0)
        bad() << k.name << " copy " << i << " has source index " << e.src_idx[i]
              << '\n';
    }
  }

  start("neighbor lists");
  for (auto const &k : kinds) {
    Entity const &e = *k.entity;
    for (auto const *nbrs : {&e.myCpys, &e.mySrcs}) {
      char const *const which = (nbrs == &e.myCpys) ? "myCpys" : "mySrcs";
      for (auto const &nbr : *nbrs) {
        if (enough())
          break;
        if (!in(nbr.pe, m_.numpe) || nbr.pe == m_.mype) {
          bad() << k.name << ' ' << which << ": rank " << nbr.pe
                << " (this is rank " << m_.mype << " of " << m_.numpe << ")\n";
          continue;
        }
        for (int const i : nbr.elements)
          if (!in(i, e.size())) {
            bad() << k.name << ' ' << which << " with rank " << nbr.pe
                  << " names " << k.name << ' ' << i << ", out of [0, "
                  << e.size() << ")\n";
            break;
          }
      }
    }
  }
}

/* -------------------------------- census --------------------------------- */

void Auditor::census() {
  os_ << "  points:  " << tally(pmask_, pll_) << '\n';
  os_ << "  edges:   " << tally(emask_, ell_) << '\n';
  os_ << "  faces:   " << tally(fmask_, fll_) << '\n';
  os_ << "  zones:   " << tally(zmask_, zll_) << '\n';
  os_ << "  corners: " << tally(cmask_, cll_) << '\n';
  os_ << "  sides:   " << tally(smask_, sll_) << '\n';
  if (all_ > 0)
    os_ << "  iotas:   " << tally(amask_, all_) << '\n';

  /* Entities that nothing points at.  These are not errors -- a partition can
     carry an entity only its ghost neighbors use -- but a mesh generator that
     leaves them behind is usually leaving something else behind too. */
  std::vector<char> upoint(static_cast<size_t>(pll_), 0),
      uedge(static_cast<size_t>(ell_), 0), uface(static_cast<size_t>(fll_), 0),
      uzone(static_cast<size_t>(zll_), 0),
      ucorner(static_cast<size_t>(cll_), 0);
  for (int c = 0; c < cll_; ++c) {
    if (cmask_[c] == 0)
      continue;
    if (in(c2p_[c], pll_))
      upoint[c2p_[c]] = 1;
    if (in(c2z_[c], zll_))
      uzone[c2z_[c]] = 1;
  }
  for (int s = 0; s < sll_; ++s) {
    if (smask_[s] == 0)
      continue;
    if (in(s2e_[s], ell_))
      uedge[s2e_[s]] = 1;
    if (in(s2f_[s], fll_))
      uface[s2f_[s]] = 1;
    if (in(s2c1_[s], cll_))
      ucorner[s2c1_[s]] = 1;
    if (in(s2c2_[s], cll_))
      ucorner[s2c2_[s]] = 1;
  }
  auto orphans = [](std::vector<char> const &used, SHORTV const &mask) {
    int n{0};
    for (size_t i = 0; i < used.size(); ++i)
      if (mask[i] != 0 && !used[i])
        n += 1;
    return n;
  };
  os_ << "  unreferenced non-ghost entities: " << orphans(upoint, pmask_)
      << " points, " << orphans(uedge, emask_) << " edges, "
      << orphans(uface, fmask_) << " faces, " << orphans(ucorner, cmask_)
      << " corners, " << orphans(uzone, zmask_) << " zones\n";

  /* The corner count of a zone is the shape of the zone: 4 for a tetrahedron,
     8 for a hexahedron, and whatever a polyhedral mesh brings. */
  std::vector<int> ncorners(static_cast<size_t>(zll_), 0);
  for (int c = 0; c < cl_; ++c)
    if (cmask_[c] != 0 && in(c2z_[c], zll_))
      ncorners[c2z_[c]] += 1;
  int lo{0}, hi{0};
  for (int z = 0; z < zll_; ++z) {
    if (zmask_[z] < 1)
      continue;
    lo = (lo == 0) ? ncorners[z] : std::min(lo, ncorners[z]);
    hi = std::max(hi, ncorners[z]);
  }
  os_ << "  corners per real zone: " << lo << " to " << hi << '\n';

  if (all_ > 0) {
    std::vector<int> per_side(static_cast<size_t>(sll_), 0);
    for (int a = 0; a < all_; ++a)
      if (amask_[a] != 0 && in(a2s_[a], sll_))
        per_side[a2s_[a]] += 1;
    int none{0}, one{0}, two{0}, more{0};
    for (int s = 0; s < sll_; ++s) {
      if (smask_[s] == 0)
        continue;
      switch (per_side[s]) {
      case 0:
        none += 1;
        break;
      case 1:
        one += 1;
        break;
      case 2:
        two += 1;
        break;
      default:
        more += 1;
      }
    }
    os_ << "  non-ghost sides by iota count: " << none << " with none, " << one
        << " with one, " << two << " with two, " << more << " with more\n";
  }

  /* Does any local zone have zmask < 1, and if so does it own real corners?
     This is what decides whether the zone filter in gradzatz_invert and the
     absence of one in gradzatz can disagree (finding F2). */
  std::vector<char> owns_real(static_cast<size_t>(zll_), 0);
  for (int c = 0; c < cl_; ++c)
    if (in(c2z_[c], zll_) && cmask_[c] >= 1)
      owns_real[c2z_[c]] = 1;
  int unreal_zones{0}, unreal_with_real_corners{0};
  for (int z = 0; z < zl_; ++z)
    if (zmask_[z] < 1) {
      unreal_zones += 1;
      unreal_with_real_corners += owns_real[z];
    }
  os_ << "  local zones with zmask < 1: " << unreal_zones << " ("
      << unreal_with_real_corners << " owning cmask >= 1 corners)\n";

  /* A cmask < 0 corner only changes gradzatz_invert if it carries volume:
     that is the whole of finding F1.  Corners fed exclusively by phantom
     (smask < 0) sides have corner_vol == 0 and are numerically inert. */
  if (opt_.geometry && ranges_ok_) {
    auto const &corner_vol = m_.ds->caccess_dblv("corner_vol");
    int neg_vol_corners{0};
    for (int c = 0; c < cl_; ++c)
      if (cmask_[c] < 0 && corner_vol[c] != 0.0)
        neg_vol_corners += 1;
    os_ << "  cmask < 0 corners carrying volume: " << neg_vol_corners
        << " (nonzero => m:z>c filtering changes gradzatz_invert)\n";
  }
}

Audit_Result Auditor::run() {
  ranges(opt_.topology);
  if (opt_.topology) {
    connectivity();
    if (all_ > 0)
      iota_maps();
  } else if (!ranges_ok_) {
    /* Nothing below can run, and the caller turned off the group that would
       have said why. */
    start("index ranges");
    bad() << "one or more maps have out-of-range entries; enable the topology "
             "group for the details\n";
  }
  if (opt_.masks)
    masks();
  if (opt_.derived) {
    if (derived_safe_)
      derived();
    else {
      finish();
      skip("derived", "the maps it is built from are not index-safe");
    }
  }
  if (opt_.geometry) {
    if (ranges_ok_)
      geometry();
    else {
      finish();
      skip("geometry", "the maps it is derived from are out of range");
    }
  }
  if (opt_.parallel)
    parallel();
  finish();
  if (opt_.census)
    census();
  return res_;
}

} // namespace

Audit_Result audit(
    Mesh const &mesh, std::ostream &os, Audit_Options const &opts) {
  return Auditor{mesh, os, opts}.run();
}

} // namespace SOA_Idx
} // namespace Ume
