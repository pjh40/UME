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
  \file Ume/Partition_Metrics.hh

  What a decomposition costs, in the terms the kernels actually pay.

  Everything here is computed from a `Decomposition` and the serial mesh
  alone.  No partition mesh is built, so sweeping a dozen methods across a
  dozen rank counts with `--no-write` is cheap, and the numbers are available
  before deciding whether the split is worth writing out.

  ### The headline number is `duplicated` for points

  A point held by \f$\lambda(p)\f$ ranks needs \f$\lambda(p)-1\f$ copies, and
  each copy is one element in one direction of `points.gathscat`.  So
  \f$\sum_p (\lambda(p)-1)\f$ is not a surrogate for UME's communication -- it
  is exactly the number of point elements on the wire.  `gradzatz` does a
  `DBLV` and a `VEC3V` gathscat plus a `VEC3V` scatter, so bytes per call are
  about \f$88\sum_p(\lambda(p)-1)\f$; `report` prints that, because it is a
  number to put next to a stopwatch.

  `max_neighbors` deserves equal billing.  It is the latency-bound half of the
  story: a scattered zone map drives it to \f$N-1\f$ while a blocked one holds
  it near 26, at very similar volume.

  ### One measure, two mechanisms

  Points, edges and faces on a cut are duplicated into the real index range of
  every rank that touches them, so their \f$\lambda\f$ is how many ranks hold
  them.  Zones, corners and sides are strictly partitioned and shadowed by
  ghosts, so their \f$\lambda-1\f$ is how many ghost copies exist.  Both are
  reported in the same column because both are elements moved per exchange;
  the distinction is `Decompose.hh`'s subject, not this file's.
*/

#ifndef UME_PARTITION_METRICS_HH
#define UME_PARTITION_METRICS_HH 1

#include "Ume/Decompose.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include <array>
#include <iosfwd>
#include <string_view>
#include <vector>

namespace Ume {
namespace Part {

/* `Kind`, `kind_name` and `is_shared_kind` come from Decompose.hh, which is
   where the order they name is decided. */

//! How one entity kind spread itself over the ranks
/*! Every per-rank vector has length `Metrics::nparts`. */
struct Kind_Metrics {
  //! Entities in each rank's real index range
  std::vector<long> local;
  //! Of those, the ones the serial mesh calls real (`mask >= 1`)
  std::vector<long> real;
  /*! What each rank is charged with, which is what `imbalance` measures.

      For a shared kind this is `local`, because every point a rank holds is
      one the kernels visit -- `gradzatp` runs over all of them and takes the
      `mask == -1` branch on the boundary ones.  For a ghosted kind it is
      `real`, because a `mask < 1` zone is a phantom that the kernels skip.

      The distinction is not pedantry: on a brick only the strictly interior
      points have `mask >= 1` -- six of sixty on a 4x3x2 -- so a point
      imbalance taken over `real` would be measuring the mesh's surface
      rather than the decomposition. */
  std::vector<long> load;
  //! Entities in each rank's ghost index range
  std::vector<long> ghost;
  //! Local entities that some other rank also holds, as a copy or as a ghost
  std::vector<long> shared;
  //! Elements each rank sends per exchange, summed over its `mySrcs`
  std::vector<long> src_elems;
  //! Elements each rank receives per exchange, summed over its `myCpys`
  std::vector<long> cpy_elems;
  //! Ranks each rank exchanges this kind with, in either direction
  std::vector<int> neighbors;

  //! How many of this kind the serial mesh had, in total and real
  long serial{0}, serial_real{0};
  long total_local{0}, total_real{0}, total_load{0}, total_ghost{0};
  /*! \f$\sum(\lambda-1)\f$: elements moved in one direction of one exchange
      of this kind.  Equal to `total_ghost` for a ghosted kind, and to the
      number of surplus instances for a shared one. */
  long duplicated{0};
  //! Heaviest rank's `load` over the mean, or 0 if there is no load at all
  double imbalance{0.0};
  int max_neighbors{0};
  double mean_neighbors{0.0};
};

//! What a decomposition costs
struct Metrics {
  int nparts{0};
  std::array<Kind_Metrics, NKINDS> kinds;

  //! Faces with a real zone on each side, whose two zones are on one rank
  long interior_faces{0};
  //! Faces with a real zone on each side, on different ranks
  long cut_faces{0};
  //! Distinct points on the cut faces, summed -- the `FACE_POINTS` edge cut
  long cut_face_points{0};

  //! Ranks each rank exchanges anything at all with
  std::vector<int> neighbors;
  int max_neighbors{0};
  double mean_neighbors{0.0};

  /*! Shared points over held points, per rank: the fraction of a rank's
      points that sit on a cut.  Dimensionless, so it compares across rank
      counts, and it is what falls when a method makes rounder subdomains. */
  double max_surface_to_volume{0.0}, mean_surface_to_volume{0.0};

  //! Elements sent from each rank to each other, all kinds, row-major
  std::vector<long> comm_matrix;
  //! Elements `from` sends to `to` per exchange sweep
  long sent(int from, int to) const {
    return comm_matrix[static_cast<std::size_t>(from) *
            static_cast<std::size_t>(nparts) +
        static_cast<std::size_t>(to)];
  }
  //! Every element on the wire in one direction of one sweep over all kinds
  long total_comm_elems() const;
};

//! Measure `d`, which must have been built from `src`
Metrics measure(SOA_Idx::Mesh const &src, Decomposition const &d);

//! Write a human-readable report
void report(Metrics const &m, std::ostream &os);

//! Write the column names `report_csv` writes values for
/*! Kept separate so a sweep writes one header and then appends a row per
    `(method, nparts)`, which is directly plottable. */
void report_csv_header(std::ostream &os);

//! Append one comma-separated row, labelled with the method that produced it
void report_csv(Metrics const &m, std::string_view method, std::ostream &os);

//! Write `comm_matrix` as `nparts` comma-separated rows
void write_comm_matrix(Metrics const &m, std::ostream &os);

} // namespace Part
} // namespace Ume

#endif
