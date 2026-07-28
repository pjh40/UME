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
  \file Ume/Partition_Graph.cc
*/

#include "Ume/Partition_Graph.hh"
#include "Ume/VecN.hh"
#include <algorithm>
#include <cmath>
#include <utility>

namespace Ume {
namespace Part {

DS_Types::VEC3V_T const &zone_centroids(SOA_Idx::Mesh const &mesh) {
  /* `caccess_*` is const but forces the lazy `init_()`, which is why a
     Partitioner can take the mesh by const reference and still ask for a
     derived variable. */
  return mesh.ds->caccess_vec3v("zcoord");
}

std::span<int const> Zone_Graph::neighbors(int v) const {
  auto const b = static_cast<std::size_t>(xadj[static_cast<std::size_t>(v)]);
  auto const e =
      static_cast<std::size_t>(xadj[static_cast<std::size_t>(v) + 1]);
  return std::span<int const>{adjncy.data() + b, e - b};
}

std::span<int const> Zone_Graph::weights(int v) const {
  auto const b = static_cast<std::size_t>(xadj[static_cast<std::size_t>(v)]);
  auto const e =
      static_cast<std::size_t>(xadj[static_cast<std::size_t>(v) + 1]);
  return std::span<int const>{adjwgt.data() + b, e - b};
}

namespace {

//! The number of distinct points on each face, indexed by face
/*! A face's points are the endpoints of its sides.  Both of the zones behind
    an interior face contribute a ring of sides naming the same points, which
    the deduplication absorbs. */
std::vector<int> face_point_counts(SOA_Idx::Mesh const &mesh) {
  int const sll = mesh.sides.size(), fll = mesh.faces.size();
  auto const &smask = mesh.sides.mask;
  auto const &s2f = mesh.ds->caccess_intv("m:s>f");
  auto const &s2p1 = mesh.ds->caccess_intv("m:s>p1");
  auto const &s2p2 = mesh.ds->caccess_intv("m:s>p2");

  std::vector<std::pair<int, int>> fp;
  fp.reserve(static_cast<std::size_t>(2 * sll));
  for (int s = 0; s < sll; ++s) {
    if (smask[s] < 1)
      continue;
    int const f = s2f[static_cast<std::size_t>(s)];
    if (f < 0 || f >= fll)
      continue;
    fp.emplace_back(f, s2p1[static_cast<std::size_t>(s)]);
    fp.emplace_back(f, s2p2[static_cast<std::size_t>(s)]);
  }
  std::ranges::sort(fp);
  fp.erase(std::unique(fp.begin(), fp.end()), fp.end());

  std::vector<int> n(static_cast<std::size_t>(fll), 0);
  for (auto const &[f, p] : fp)
    n[static_cast<std::size_t>(f)] += 1;
  return n;
}

//! The flat area of each face, indexed by face
/*! The side surface vectors of one face come in pairs, one side from each of
    the two zones behind it, so the ring sums to twice the area.  Halving is
    the same answer `calc_face_area` reaches by tagging each side's `m:s>s2`
    partner, and here every face of interest has both of its rings real. */
std::vector<double> face_areas(SOA_Idx::Mesh const &mesh) {
  int const sll = mesh.sides.size(), fll = mesh.faces.size();
  auto const &smask = mesh.sides.mask;
  auto const &s2f = mesh.ds->caccess_intv("m:s>f");
  auto const &surz = mesh.ds->caccess_vec3v("side_surz");

  std::vector<double> area(static_cast<std::size_t>(fll), 0.0);
  for (int s = 0; s < sll; ++s) {
    if (smask[s] < 1)
      continue;
    int const f = s2f[static_cast<std::size_t>(s)];
    if (f < 0 || f >= fll)
      continue;
    area[static_cast<std::size_t>(f)] +=
        0.5 * vectormag(surz[static_cast<std::size_t>(s)]);
  }
  return area;
}

//! What cutting each face costs, indexed by face
std::vector<int> face_weights(SOA_Idx::Mesh const &mesh, Options const &opt) {
  int const fll = mesh.faces.size();
  switch (opt.ewgt) {
  case Edge_Weight::UNIT:
    return std::vector<int>(static_cast<std::size_t>(fll), 1);
  case Edge_Weight::FACE_POINTS:
    /* The points on a cut face are exactly the points that may have to be
       duplicated across that cut, which is what makes this the default. */
    return face_point_counts(mesh);
  case Edge_Weight::FACE_AREA:
    break;
  }

  std::vector<double> const area = face_areas(mesh);
  double const top = area.empty() ? 0.0 : *std::ranges::max_element(area);
  std::vector<int> w(static_cast<std::size_t>(fll), 1);
  if (!(top > 0.0))
    return w; //!< a mesh with no area to speak of weighs every face the same
  for (std::size_t f = 0; f < area.size(); ++f)
    w[f] =
        std::max(1, static_cast<int>(std::lround(area[f] * AREA_SCALE / top)));
  return w;
}

} // namespace

Zone_Graph build_zone_graph(SOA_Idx::Mesh const &mesh, Options const &opt) {
  int const zll = mesh.zones.size(), fll = mesh.faces.size();
  auto const &zmask = mesh.zones.mask;
  auto const &fmask = mesh.faces.mask;
  auto const &f2z1 = mesh.ds->caccess_intv("m:f>z1");
  auto const &f2z2 = mesh.ds->caccess_intv("m:f>z2");

  Zone_Graph g;
  g.v2z = real_zones(mesh);
  g.nverts = static_cast<int>(g.v2z.size());
  g.z2v.assign(static_cast<std::size_t>(zll), -1);
  for (int v = 0; v < g.nverts; ++v)
    g.z2v[static_cast<std::size_t>(g.v2z[static_cast<std::size_t>(v)])] = v;

  std::vector<int> const zwgt = zone_weights(mesh, opt);
  g.vwgt.reserve(g.v2z.size());
  for (int const z : g.v2z)
    g.vwgt.push_back(zwgt[static_cast<std::size_t>(z)]);

  std::vector<int> const fwgt = face_weights(mesh, opt);

  /* Both directions of every edge, then sorted and merged: one pass, no
     per-vertex containers, and the merge is what folds a repeated zone pair
     into a single weighted edge. */
  std::vector<std::pair<std::pair<int, int>, int>> arcs;
  for (int f = 0; f < fll; ++f) {
    if (fmask[f] == 0)
      continue;
    int const z1 = f2z1[static_cast<std::size_t>(f)];
    int const z2 = f2z2[static_cast<std::size_t>(f)];
    if (z1 < 0 || z2 < 0 || z1 >= zll || z2 >= zll || z1 == z2)
      continue;
    if (zmask[z1] < 1 || zmask[z2] < 1)
      continue;
    int const v1 = g.z2v[static_cast<std::size_t>(z1)];
    int const v2 = g.z2v[static_cast<std::size_t>(z2)];
    int const w = fwgt[static_cast<std::size_t>(f)];
    arcs.push_back({{v1, v2}, w});
    arcs.push_back({{v2, v1}, w});
  }
  std::ranges::sort(arcs);

  g.xadj.assign(static_cast<std::size_t>(g.nverts) + 1, 0);
  g.adjncy.reserve(arcs.size());
  g.adjwgt.reserve(arcs.size());
  for (std::size_t i = 0; i < arcs.size();) {
    std::size_t j = i;
    int w = 0;
    while (j < arcs.size() && arcs[j].first == arcs[i].first)
      w += arcs[j++].second;
    g.adjncy.push_back(arcs[i].first.second);
    g.adjwgt.push_back(w);
    g.xadj[static_cast<std::size_t>(arcs[i].first.first) + 1] += 1;
    i = j;
  }
  for (int v = 0; v < g.nverts; ++v)
    g.xadj[static_cast<std::size_t>(v) + 1] +=
        g.xadj[static_cast<std::size_t>(v)];
  return g;
}

} // namespace Part
} // namespace Ume
