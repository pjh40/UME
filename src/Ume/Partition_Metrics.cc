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
  \file Ume/Partition_Metrics.cc
*/

#include "Ume/Partition_Metrics.hh"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <set>
#include <span>
#include <utility>

namespace Ume {
namespace Part {

using SOA_Idx::Entity;
using SOA_Idx::Mesh;

namespace {

//! For each global entity, the ranks holding it, ascending
/*! The layouts are keyed rank-to-entities, which is what building a rank's
    mesh needs; every question here is the other way round -- how many ranks
    hold this entity, and which. */
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
        rank_[static_cast<std::size_t>(fill[static_cast<std::size_t>(
            gid[static_cast<std::size_t>(i)])]++)] = r;
  }

  std::span<int const> operator[](int g) const {
    auto const b =
        static_cast<std::size_t>(start_[static_cast<std::size_t>(g)]);
    auto const e =
        static_cast<std::size_t>(start_[static_cast<std::size_t>(g) + 1]);
    return std::span<int const>{rank_.data() + b, e - b};
  }

private:
  std::vector<int> start_, rank_;
};

//! Heaviest entry over the mean, or 0 when there is nothing to balance
double imbalance_of(std::vector<long> const &load, long total) {
  if (total <= 0 || load.empty())
    return 0.0;
  long const heaviest = *std::ranges::max_element(load);
  double const mean =
      static_cast<double>(total) / static_cast<double>(load.size());
  return static_cast<double>(heaviest) / mean;
}

//! Which of the two zones behind each face, if both are real, and their ranks
struct Face_Cut {
  bool interior{false}; //!< a real zone on each side
  bool cut{false}; //!< ...and they are on different ranks
};

std::vector<Face_Cut> face_cuts(Mesh const &src, Decomposition const &d) {
  int const fll = src.faces.size(), zll = src.zones.size();
  auto const &fmask = src.faces.mask;
  auto const &zmask = src.zones.mask;
  auto const &f2z1 = src.ds->caccess_intv("m:f>z1");
  auto const &f2z2 = src.ds->caccess_intv("m:f>z2");

  std::vector<Face_Cut> cuts(static_cast<std::size_t>(fll));
  for (int f = 0; f < fll; ++f) {
    if (fmask[f] == 0)
      continue;
    int const z1 = f2z1[static_cast<std::size_t>(f)];
    int const z2 = f2z2[static_cast<std::size_t>(f)];
    if (z1 < 0 || z2 < 0 || z1 >= zll || z2 >= zll)
      continue;
    if (zmask[z1] < 1 || zmask[z2] < 1)
      continue;
    Face_Cut &c = cuts[static_cast<std::size_t>(f)];
    c.interior = true;
    c.cut = d.z2r[static_cast<std::size_t>(z1)] !=
        d.z2r[static_cast<std::size_t>(z2)];
  }
  return cuts;
}

//! The distinct points on the cut faces, summed
/*! This is the `Edge_Weight::FACE_POINTS` edge cut, and the reason to report
    it next to `duplicated` for points is that the two differ: a point on
    several cut faces is charged once by the second and once per face by the
    first, which is exactly the gap between what a graph partitioner minimizes
    and what goes on the wire. */
long cut_face_points(Mesh const &src, std::vector<Face_Cut> const &cuts) {
  int const sll = src.sides.size();
  auto const &smask = src.sides.mask;
  auto const &s2f = src.ds->caccess_intv("m:s>f");
  auto const &s2p1 = src.ds->caccess_intv("m:s>p1");
  auto const &s2p2 = src.ds->caccess_intv("m:s>p2");

  std::vector<std::pair<int, int>> fp;
  for (int s = 0; s < sll; ++s) {
    if (smask[s] < 1)
      continue;
    int const f = s2f[static_cast<std::size_t>(s)];
    if (f < 0 || f >= static_cast<int>(cuts.size()) ||
        !cuts[static_cast<std::size_t>(f)].cut)
      continue;
    fp.emplace_back(f, s2p1[static_cast<std::size_t>(s)]);
    fp.emplace_back(f, s2p2[static_cast<std::size_t>(s)]);
  }
  std::ranges::sort(fp);
  return static_cast<long>(std::unique(fp.begin(), fp.end()) - fp.begin());
}

} // namespace

char const *kind_name(int const kind) {
  constexpr char const *NAMES[NKINDS]{
      "point", "edge", "face", "zone", "corner", "side", "iota"};
  return (kind >= 0 && kind < NKINDS) ? NAMES[kind] : "?";
}

long Metrics::total_comm_elems() const {
  long total = 0;
  for (long const n : comm_matrix)
    total += n;
  return total;
}

Metrics measure(Mesh const &src, Decomposition const &d) {
  int const np = d.nparts;
  auto const nps = static_cast<std::size_t>(np);

  Metrics m;
  m.nparts = np;
  m.comm_matrix.assign(nps * nps, 0);
  m.neighbors.assign(nps, 0);

  auto const layouts = d.layouts();
  Entity const *const ents[NKINDS]{&src.points, &src.edges, &src.faces,
      &src.zones, &src.corners, &src.sides, &src.iotas};

  std::vector<std::set<int>> all_nbr(nps);

  for (int k = 0; k < NKINDS; ++k) {
    Entity_Layout const &l = *layouts[static_cast<std::size_t>(k)];
    Kind_Metrics &km = m.kinds[static_cast<std::size_t>(k)];
    auto const &mask = ents[k]->mask;
    int const nglobal = l.nglobal();

    /* Ghost zones, corners and sides are only exchanged when the
       decomposition was asked to; points, edges and faces on a cut always
       are, because a COPY that is never filled is simply wrong. */
    bool const shared = is_shared_kind(k);
    bool const comm = shared || d.opt.ghost_entity_comm;

    km.local.assign(nps, 0);
    km.real.assign(nps, 0);
    km.ghost.assign(nps, 0);
    km.load.assign(nps, 0);
    km.shared.assign(nps, 0);
    km.src_elems.assign(nps, 0);
    km.cpy_elems.assign(nps, 0);
    km.neighbors.assign(nps, 0);
    km.serial = nglobal;

    for (int g = 0; g < nglobal; ++g)
      if (mask[g] >= 1)
        km.serial_real += 1;

    for (int r = 0; r < np; ++r) {
      km.local[static_cast<std::size_t>(r)] = l.local_size(r);
      km.ghost[static_cast<std::size_t>(r)] = l.ghost_local_size(r);
      km.total_local += l.local_size(r);
      km.total_ghost += l.ghost_local_size(r);
      for (int i = 0; i < l.local_size(r); ++i)
        if (mask[l.global_index(r, i)] >= 1)
          km.real[static_cast<std::size_t>(r)] += 1;
      km.total_real += km.real[static_cast<std::size_t>(r)];
      km.load[static_cast<std::size_t>(r)] = shared
          ? km.local[static_cast<std::size_t>(r)]
          : km.real[static_cast<std::size_t>(r)];
      km.total_load += km.load[static_cast<std::size_t>(r)];
    }

    std::vector<std::set<int>> nbr(nps);
    auto note = [&](int const from, int const to) {
      m.comm_matrix[static_cast<std::size_t>(from) * nps +
          static_cast<std::size_t>(to)] += 1;
      nbr[static_cast<std::size_t>(from)].insert(to);
      nbr[static_cast<std::size_t>(to)].insert(from);
      all_nbr[static_cast<std::size_t>(from)].insert(to);
      all_nbr[static_cast<std::size_t>(to)].insert(from);
    };

    /* This mirrors `fill_common` in `Decompose.cc` exactly, entity for
       entity, which is what makes the numbers below predictions of what
       `extract` will produce rather than an independent estimate of it. */
    Presence const locals{l.start, l.gid, np, nglobal};
    Presence const ghosts{l.gstart, l.ggid, np, nglobal};
    for (int g = 0; g < nglobal; ++g) {
      std::span<int const> const on = shared ? locals[g] : ghosts[g];
      if (on.empty())
        continue;
      int const owner = l.owner[static_cast<std::size_t>(g)];
      if (owner < 0)
        continue;

      if (shared) {
        if (on.size() == 1)
          continue;
        km.duplicated += static_cast<long>(on.size()) - 1;
        for (int const r : on)
          km.shared[static_cast<std::size_t>(r)] += 1;
      } else {
        km.duplicated += static_cast<long>(on.size());
        km.shared[static_cast<std::size_t>(owner)] += 1;
      }
      if (!comm)
        continue;

      for (int const r : on)
        if (r != owner) {
          km.src_elems[static_cast<std::size_t>(owner)] += 1;
          km.cpy_elems[static_cast<std::size_t>(r)] += 1;
          note(owner, r);
        }
    }

    for (int r = 0; r < np; ++r)
      km.neighbors[static_cast<std::size_t>(r)] =
          static_cast<int>(nbr[static_cast<std::size_t>(r)].size());
    km.imbalance = imbalance_of(km.load, km.total_load);
    km.max_neighbors =
        km.neighbors.empty() ? 0 : *std::ranges::max_element(km.neighbors);
    km.mean_neighbors = np == 0
        ? 0.0
        : static_cast<double>(
              std::reduce(km.neighbors.begin(), km.neighbors.end(), long{0})) /
            static_cast<double>(np);
  }

  for (int r = 0; r < np; ++r)
    m.neighbors[static_cast<std::size_t>(r)] =
        static_cast<int>(all_nbr[static_cast<std::size_t>(r)].size());
  m.max_neighbors =
      m.neighbors.empty() ? 0 : *std::ranges::max_element(m.neighbors);
  m.mean_neighbors = np == 0
      ? 0.0
      : static_cast<double>(
            std::reduce(m.neighbors.begin(), m.neighbors.end(), long{0})) /
          static_cast<double>(np);

  /* Ranks with no points at all -- which happens as soon as there are more
     parts than zones -- have no surface-to-volume ratio rather than a ratio
     of zero, so they are left out of both the maximum and the mean. */
  Kind_Metrics const &pts = m.kinds[POINT];
  double sum = 0.0;
  int counted = 0;
  for (int r = 0; r < np; ++r) {
    long const volume = pts.local[static_cast<std::size_t>(r)];
    if (volume <= 0)
      continue;
    double const ratio =
        static_cast<double>(pts.shared[static_cast<std::size_t>(r)]) /
        static_cast<double>(volume);
    m.max_surface_to_volume = std::max(m.max_surface_to_volume, ratio);
    sum += ratio;
    counted += 1;
  }
  m.mean_surface_to_volume =
      counted == 0 ? 0.0 : sum / static_cast<double>(counted);

  std::vector<Face_Cut> const cuts = face_cuts(src, d);
  for (Face_Cut const &c : cuts) {
    m.interior_faces += c.interior ? 1 : 0;
    m.cut_faces += c.cut ? 1 : 0;
  }
  m.cut_face_points = cut_face_points(src, cuts);

  return m;
}

/* -------------------------------- reports -------------------------------- */

namespace {

long minimum(std::vector<long> const &v) {
  return v.empty() ? 0 : *std::ranges::min_element(v);
}
long maximum(std::vector<long> const &v) {
  return v.empty() ? 0 : *std::ranges::max_element(v);
}

} // namespace

void report(Metrics const &m, std::ostream &os) {
  auto const flags = os.flags();
  auto const precision = os.precision();

  os << "Decomposition metrics for " << m.nparts << " parts\n";
  os << "  " << std::setw(7) << std::left << "kind" << std::right
     << std::setw(9) << "serial" << std::setw(8) << "real" << std::setw(10)
     << "min load" << std::setw(10) << "max load" << std::setw(8) << "imbal"
     << std::setw(9) << "ghost" << std::setw(10) << "lambda-1" << std::setw(6)
     << "nbrs" << '\n';
  for (int k = 0; k < NKINDS; ++k) {
    Kind_Metrics const &km = m.kinds[static_cast<std::size_t>(k)];
    if (km.serial == 0)
      continue;
    os << "  " << std::setw(7) << std::left << kind_name(k) << std::right
       << std::setw(9) << km.serial << std::setw(8) << km.serial_real
       << std::setw(10) << minimum(km.load) << std::setw(10) << maximum(km.load)
       << std::setw(8) << std::fixed << std::setprecision(3) << km.imbalance
       << std::defaultfloat << std::setw(9) << km.total_ghost << std::setw(10)
       << km.duplicated << std::setw(6) << km.max_neighbors << '\n';
  }

  long const lambda1 = m.kinds[POINT].duplicated;
  /* `gradzatz` moves a DBLV and a VEC3V through `points.gathscat`, each of
     which sends and receives, plus a VEC3V `scatter`: (8 + 24) * 2 + 24
     bytes per duplicated point. */
  os << "  cut faces:           " << m.cut_faces << " of " << m.interior_faces
     << " interior, " << m.cut_face_points << " points on them\n"
     << "  neighbors:           max " << m.max_neighbors << ", mean "
     << std::fixed << std::setprecision(2) << m.mean_neighbors << '\n'
     << "  surface/volume:      max " << m.max_surface_to_volume << ", mean "
     << m.mean_surface_to_volume << std::defaultfloat << '\n'
     << "  point lambda-1:      " << lambda1 << ", about " << 88 * lambda1
     << " bytes on the wire per gradzatz\n"
     << "  elements exchanged:  " << m.total_comm_elems()
     << " in one direction of one sweep over every kind\n";

  os.flags(flags);
  os.precision(precision);
}

void report_csv_header(std::ostream &os) {
  os << "method,nparts,zones,zone_imbalance,corner_imbalance,cut_faces,"
     << "interior_faces,cut_face_points,point_lambda1,edge_lambda1,"
     << "face_lambda1,ghost_zones,ghost_corners,ghost_sides,max_neighbors,"
     << "mean_neighbors,max_surface_to_volume,mean_surface_to_volume,"
     << "comm_elems\n";
}

void report_csv(Metrics const &m, std::string_view method, std::ostream &os) {
  auto const flags = os.flags();
  auto const precision = os.precision();
  os << std::defaultfloat << std::setprecision(6);

  os << method << ',' << m.nparts << ',' << m.kinds[ZONE].serial_real << ','
     << m.kinds[ZONE].imbalance << ',' << m.kinds[CORNER].imbalance << ','
     << m.cut_faces << ',' << m.interior_faces << ',' << m.cut_face_points
     << ',' << m.kinds[POINT].duplicated << ',' << m.kinds[EDGE].duplicated
     << ',' << m.kinds[FACE].duplicated << ',' << m.kinds[ZONE].total_ghost
     << ',' << m.kinds[CORNER].total_ghost << ',' << m.kinds[SIDE].total_ghost
     << ',' << m.max_neighbors << ',' << m.mean_neighbors << ','
     << m.max_surface_to_volume << ',' << m.mean_surface_to_volume << ','
     << m.total_comm_elems() << '\n';

  os.flags(flags);
  os.precision(precision);
}

void write_comm_matrix(Metrics const &m, std::ostream &os) {
  for (int r = 0; r < m.nparts; ++r) {
    for (int c = 0; c < m.nparts; ++c)
      os << (c ? "," : "") << m.sent(r, c);
    os << '\n';
  }
}

} // namespace Part
} // namespace Ume
