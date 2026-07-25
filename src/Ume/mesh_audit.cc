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
#include <cmath>
#include <map>
#include <ostream>
#include <vector>

namespace Ume {
namespace SOA_Idx {

namespace {

//! Categories used by the census
struct Census {
  int real{0}; //!< mask >= 1
  int boundary{0}; //!< mask < 0
  int ghost{0}; //!< mask == 0
};

Census tally(std::vector<short> const &mask, int const upto) {
  Census c;
  for (int i = 0; i < upto; ++i) {
    if (mask[static_cast<size_t>(i)] >= 1)
      c.real += 1;
    else if (mask[static_cast<size_t>(i)] < 0)
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

/*! Accumulates problems for one named check, printing at most
    `Audit_Options::max_reports` of them. */
class Auditor {
public:
  Auditor(Mesh const &mesh, std::ostream &os, Audit_Options const &opts)
      : m_{mesh}, os_{os}, opt_{opts} {}

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

  void topology();
  void masks();
  void geometry();
  void census();

  Mesh const &m_;
  std::ostream &os_;
  Audit_Options opt_;
  Audit_Result res_;
  char const *name_{nullptr};
  int hits_{0};
  std::ostream null_{nullptr};
};

void Auditor::topology() {
  int const pll = m_.points.size();
  int const zll = m_.zones.size();
  int const cll = m_.corners.size();
  int const sll = m_.sides.size();
  int const fll = m_.faces.size();
  int const ell = m_.edges.size();
  int const cl = m_.corners.local_size();

  auto const &s2z = m_.ds->caccess_intv("m:s>z");
  auto const &s2p1 = m_.ds->caccess_intv("m:s>p1");
  auto const &s2p2 = m_.ds->caccess_intv("m:s>p2");
  auto const &s2e = m_.ds->caccess_intv("m:s>e");
  auto const &s2f = m_.ds->caccess_intv("m:s>f");
  auto const &s2c1 = m_.ds->caccess_intv("m:s>c1");
  auto const &s2c2 = m_.ds->caccess_intv("m:s>c2");
  auto const &s2s2 = m_.ds->caccess_intv("m:s>s2");
  auto const &s2s3 = m_.ds->caccess_intv("m:s>s3");
  auto const &s2s4 = m_.ds->caccess_intv("m:s>s4");
  auto const &s2s5 = m_.ds->caccess_intv("m:s>s5");
  auto const &c2p = m_.ds->caccess_intv("m:c>p");
  auto const &c2z = m_.ds->caccess_intv("m:c>z");
  auto const &e2p1 = m_.ds->caccess_intv("m:e>p1");
  auto const &e2p2 = m_.ds->caccess_intv("m:e>p2");
  auto const &smask = m_.sides.mask;
  auto const &cmask = m_.corners.mask;
  auto const &zmask = m_.zones.mask;

  auto in = [](int const i, int const n) { return i >= 0 && i < n; };

  start("side map ranges");
  for (int s = 0; s < sll && !enough(); ++s) {
    if (smask[s] == 0)
      continue; // ghost sides carry no local connectivity
    if (!in(s2z[s], zll) || !in(s2p1[s], pll) || !in(s2p2[s], pll) ||
        !in(s2e[s], ell) || !in(s2f[s], fll) || !in(s2c1[s], cll) ||
        !in(s2c2[s], cll) || !in(s2s2[s], sll) || !in(s2s3[s], sll) ||
        !in(s2s4[s], sll) || !in(s2s5[s], sll))
      bad() << "side " << s << " has an out-of-range map entry\n";
    else if (s2p1[s] == s2p2[s])
      bad() << "side " << s << " has p1 == p2 == " << s2p1[s] << '\n';
  }

  start("edge map ranges");
  for (int e = 0; e < ell && !enough(); ++e)
    if (!in(e2p1[e], pll) || !in(e2p2[e], pll) || e2p1[e] == e2p2[e])
      bad() << "edge " << e << " has bad endpoints\n";

  start("corner map ranges");
  for (int c = 0; c < cll && !enough(); ++c) {
    if (!in(c2p[c], pll))
      bad() << "corner " << c << " points at point " << c2p[c] << '\n';
    else if (cmask[c] != 0 && !in(c2z[c], zll))
      bad() << "corner " << c << " (mask " << cmask[c] << ") points at zone "
            << c2z[c] << '\n';
  }

  /* A ghost corner's `m:c>z` is allowed to be garbage, but only in the
     recognizable way the library already guards against: out of range, or
     naming a zone that is itself a ghost. */
  start("ghost corners do not name real zones");
  for (int c = 0; c < cll && !enough(); ++c)
    if (cmask[c] == 0 && in(c2z[c], zll) && zmask[c2z[c]] != 0)
      bad() << "ghost corner " << c << " names real zone " << c2z[c] << '\n';

  start("m:s>s2 pairing");
  for (int s = 0; s < sll && !enough(); ++s) {
    if (smask[s] == 0)
      continue;
    int const s2 = s2s2[s];
    if (!in(s2, sll))
      continue; // already reported by the range check
    if (s2 == s)
      bad() << "side " << s << " is its own s2\n";
    else if (s2s2[s2] != s)
      bad() << "side " << s << " -> " << s2 << " -> " << s2s2[s2]
            << " is not an involution\n";
    else if (s2f[s2] != s2f[s] || s2e[s2] != s2e[s])
      bad() << "side " << s << " and " << s2 << " share s2 but not their "
            << "(face, edge): (" << s2f[s] << ", " << s2e[s] << ") vs ("
            << s2f[s2] << ", " << s2e[s2] << ")\n";
  }

  /* Every (face, edge) is shared by exactly two sides: one per zone across
     the face, or the real side and its phantom on the mesh boundary. */
  start("(face, edge) side counts");
  {
    std::map<std::pair<int, int>, int> counts;
    for (int s = 0; s < sll; ++s)
      if (smask[s] != 0 && in(s2f[s], fll) && in(s2e[s], ell))
        counts[{s2f[s], s2e[s]}] += 1;
    for (auto const &[fe, n] : counts) {
      if (enough())
        break;
      if (n != 2)
        bad() << "face " << fe.first << ", edge " << fe.second << " has " << n
              << " sides\n";
    }
  }

  start("m:s>c1/c2 agree with (zone, p1/p2)");
  for (int s = 0; s < sll && !enough(); ++s) {
    if (smask[s] == 0)
      continue;
    int const c1 = s2c1[s], c2 = s2c2[s];
    if (!in(c1, cll) || !in(c2, cll))
      continue; // already reported
    if (c2z[c1] != s2z[s] || c2p[c1] != s2p1[s])
      bad() << "side " << s << " c1 " << c1 << " is (zone " << c2z[c1]
            << ", point " << c2p[c1] << "), expected (zone " << s2z[s]
            << ", point " << s2p1[s] << ")\n";
    else if (c2z[c2] != s2z[s] || c2p[c2] != s2p2[s])
      bad() << "side " << s << " c2 " << c2 << " is (zone " << c2z[c2]
            << ", point " << c2p[c2] << "), expected (zone " << s2z[s]
            << ", point " << s2p2[s] << ")\n";
  }

  start("no duplicate corners in a zone");
  {
    std::map<std::pair<int, int>, int> seen;
    for (int c = 0; c < cl; ++c) {
      if (cmask[c] == 0 || !in(c2z[c], zll) || !in(c2p[c], pll))
        continue;
      if (enough())
        break;
      auto const [it, added] = seen.try_emplace({c2z[c], c2p[c]}, c);
      if (!added)
        bad() << "corners " << it->second << " and " << c << " are both (zone "
              << c2z[c] << ", point " << c2p[c] << ")\n";
    }
  }
}

void Auditor::masks() {
  int const cl = m_.corners.local_size();
  int const pll = m_.points.size();
  int const zll = m_.zones.size();
  auto const &c2p = m_.ds->caccess_intv("m:c>p");
  auto const &c2z = m_.ds->caccess_intv("m:c>z");
  auto const &cmask = m_.corners.mask;
  auto const &zmask = m_.zones.mask;
  auto const &pmask = m_.points.mask;

  start("real corners belong to real zones");
  for (int c = 0; c < cl && !enough(); ++c)
    if (cmask[c] >= 1 && c2z[c] >= 0 && c2z[c] < zll && zmask[c2z[c]] < 1)
      bad() << "corner " << c << " (mask " << cmask[c] << ") is in zone "
            << c2z[c] << " (mask " << zmask[c2z[c]] << ")\n";

  /*! `gradzatp` accumulates `point_volume` only over `cmask >= 1` corners and
      then divides every `pmask != 0` point by it, so a point without such a
      corner produces 0/0 and poisons the gradient. */
  start("every real point owns a cmask >= 1 corner");
  {
    std::vector<char> owned(static_cast<size_t>(pll), 0);
    for (int c = 0; c < cl; ++c)
      if (cmask[c] >= 1 && c2p[c] >= 0 && c2p[c] < pll)
        owned[static_cast<size_t>(c2p[c])] = 1;
    for (int p = 0; p < m_.points.local_size() && !enough(); ++p)
      if (pmask[p] != 0 && !owned[static_cast<size_t>(p)])
        bad() << "point " << p << " (mask " << pmask[p]
              << ") has no cmask >= 1 corner: gradzatp will divide by zero\n";
  }
}

void Auditor::geometry() {
  int const sl = m_.sides.local_size();
  int const cl = m_.corners.local_size();
  int const zll = m_.zones.size();
  int const pl = m_.points.local_size();

  auto const &side_vol = m_.ds->caccess_dblv("side_vol");
  auto const &corner_vol = m_.ds->caccess_dblv("corner_vol");
  auto const &csurf = m_.ds->caccess_vec3v("corner_csurf");
  auto const &point_norm = m_.ds->caccess_vec3v("point_norm");
  auto const &c2z = m_.ds->caccess_intv("m:c>z");
  auto const &smask = m_.sides.mask;
  auto const &cmask = m_.corners.mask;
  auto const &zmask = m_.zones.mask;
  auto const &pmask = m_.points.mask;

  start("side_vol > 0 on real sides");
  for (int s = 0; s < sl && !enough(); ++s)
    if (smask[s] > 0 && !(side_vol[s] > 0.0))
      bad() << "side " << s << " has volume " << side_vol[s]
            << " (check the p1/p2 orientation convention)\n";

  /* Accumulate per zone straight from `m:c>z` rather than through `m:z>c`, so
     that this stays valid whatever the zone-to-corner map is filtered on. */
  std::vector<double> zvol(static_cast<size_t>(zll), 0.0);
  std::vector<Vec3> zsurf(static_cast<size_t>(zll), Vec3(0.0));
  std::vector<double> zscale(static_cast<size_t>(zll), 0.0);
  for (int c = 0; c < cl; ++c) {
    if (cmask[c] == 0 || c2z[c] < 0 || c2z[c] >= zll)
      continue;
    size_t const z = static_cast<size_t>(c2z[c]);
    zvol[z] += corner_vol[c];
    zsurf[z] += csurf[c];
    zscale[z] += vectormag(csurf[c]);
  }

  start("real zones have positive volume");
  for (int z = 0; z < zll && !enough(); ++z)
    if (zmask[z] >= 1 && !(zvol[static_cast<size_t>(z)] > 0.0))
      bad() << "zone " << z << " has volume " << zvol[static_cast<size_t>(z)]
            << '\n';

  /*! The corner surfaces of a zone are its internal facets, each counted once
      with either sign, so they have to cancel.  If they do not, `gradzatp`
      does not reproduce a constant field. */
  start("corner_csurf closes over each zone");
  for (int z = 0; z < zll && !enough(); ++z) {
    if (zmask[z] == 0)
      continue;
    size_t const zi = static_cast<size_t>(z);
    if (vectormag(zsurf[zi]) > 1e-10 * zscale[zi])
      bad() << "zone " << z << " corner surfaces sum to " << zsurf[zi]
            << " (facet area scale " << zscale[zi] << ")\n";
  }

  start("point_norm is a unit vector on boundary points");
  for (int p = 0; p < pl && !enough(); ++p) {
    double const mag = vectormag(point_norm[p]);
    if (pmask[p] < 0 && std::abs(mag - 1.0) > 1e-12)
      bad() << "point " << p << " (mask " << pmask[p] << ") has |point_norm| "
            << mag << '\n';
    else if (pmask[p] >= 1 && mag > 1e-12)
      bad() << "interior point " << p << " has a nonzero point_norm "
            << point_norm[p] << '\n';
  }
}

void Auditor::census() {
  int const cl = m_.corners.local_size();
  int const zl = m_.zones.local_size();
  auto const &c2z = m_.ds->caccess_intv("m:c>z");
  auto const &cmask = m_.corners.mask;
  auto const &zmask = m_.zones.mask;

  os_ << "  points:  " << tally(m_.points.mask, m_.points.size()) << '\n';
  os_ << "  edges:   " << tally(m_.edges.mask, m_.edges.size()) << '\n';
  os_ << "  faces:   " << tally(m_.faces.mask, m_.faces.size()) << '\n';
  os_ << "  zones:   " << tally(m_.zones.mask, m_.zones.size()) << '\n';
  os_ << "  corners: " << tally(m_.corners.mask, m_.corners.size()) << '\n';
  os_ << "  sides:   " << tally(m_.sides.mask, m_.sides.size()) << '\n';

  /* Does any local zone have zmask < 1, and if so does it own real corners?
     This is what decides whether the zone filter in gradzatz_invert and the
     absence of one in gradzatz can disagree (finding F2). */
  std::vector<char> owns_real(static_cast<size_t>(m_.zones.size()), 0);
  int neg_vol_corners{0};
  for (int c = 0; c < cl; ++c) {
    if (c2z[c] < 0 || c2z[c] >= m_.zones.size())
      continue;
    if (cmask[c] >= 1)
      owns_real[static_cast<size_t>(c2z[c])] = 1;
  }
  int unreal_zones{0}, unreal_with_real_corners{0};
  for (int z = 0; z < zl; ++z)
    if (zmask[z] < 1) {
      unreal_zones += 1;
      unreal_with_real_corners += owns_real[static_cast<size_t>(z)];
    }
  os_ << "  local zones with zmask < 1: " << unreal_zones << " ("
      << unreal_with_real_corners << " owning cmask >= 1 corners)\n";

  /* A cmask < 0 corner only changes gradzatz_invert if it carries volume:
     that is the whole of finding F1.  Corners fed exclusively by phantom
     (smask < 0) sides have corner_vol == 0 and are numerically inert. */
  if (opt_.geometry) {
    auto const &corner_vol = m_.ds->caccess_dblv("corner_vol");
    for (int c = 0; c < cl; ++c)
      if (cmask[c] < 0 && corner_vol[c] != 0.0)
        neg_vol_corners += 1;
    os_ << "  cmask < 0 corners carrying volume: " << neg_vol_corners
        << " (nonzero => m:z>c filtering changes gradzatz_invert)\n";
  }
}

Audit_Result Auditor::run() {
  if (opt_.topology)
    topology();
  if (opt_.masks)
    masks();
  if (opt_.geometry)
    geometry();
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
