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
  \file ume_serial.cc

  This is an example of a serial driver for Ume.  It simply reads one or more
  binary Ume file(s) and accesses some computed variables in the first one.
*/

#include "Ume/SOA_Idx_Mesh.hh"
#include "Ume/mesh_audit.hh"
#include "Ume/partition_check.hh"
#include "Ume/process_mgmt.hh"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

using namespace Ume::SOA_Idx;

using Mesh_Ptr = std::unique_ptr<Mesh>;

std::vector<Mesh_Ptr> read_meshes(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  /* The derived mesh variables are computed with Kokkos kernels, so Kokkos has
     to be up before anything touches the Datastore.  It also has to be up
     before the command line is read at all: Ume::initialize() takes argc by
     reference and strips the --kokkos-* arguments, so until it has run,
     argv[1] may be one of those rather than the caller's first argument --
     and `ume_serial --kokkos-num-threads=4 --audit part.00000.ume` would
     silently skip the audit and then try to open "--audit" as a mesh. */
  Ume::initialize(argc, argv);

  bool audit_only{false};
  if (argc > 1 && std::strcmp(argv[1], "--audit") == 0) {
    audit_only = true;
    /* Drop the flag so read_meshes() sees only file names. */
    for (int i = 1; i < argc - 1; ++i)
      argv[i] = argv[i + 1];
    --argc;
  }
  if (argc == 1) {
    std::cerr << "Usage: ume_serial [--audit] <ume file>+" << std::endl;
    Ume::finalize();
    return 1;
  }

  int status = 0;
  {
    std::vector<Mesh_Ptr> ranks{read_meshes(argc, argv)};
    if (ranks.empty()) {
      status = 1;
    } else {
      for (Mesh_Ptr const &m : ranks)
        Ume::attach_serial_transport(*m);

      if (audit_only) {
        /* Every rank, not just the first: an audit of one partition is an
           audit of one partition, and a decomposition can put its mistake in
           any of them. */
        Audit_Options opts;
        opts.census = true;
        for (Mesh_Ptr const &m : ranks) {
          std::cout << "\nAudit of rank " << m->mype << '\n';
          Audit_Result const res = audit(*m, std::cout, opts);
          std::cout << "  " << res.checks - res.failed_checks << '/'
                    << res.checks << " checks passed";
          if (!res)
            std::cout << ", " << res.problems << " problems";
          std::cout << std::endl;
          if (!res)
            status = 1;
        }

        /* What no per-rank audit can reach: whether the ranks agree with each
           other about the entities they share.  It needs the whole set, so a
           partial one is skipped rather than reported as a disagreement. */
        int const numpe = ranks[0]->numpe;
        if (numpe > 1 && static_cast<int>(ranks.size()) == numpe) {
          std::vector<Mesh const *> all;
          all.reserve(ranks.size());
          for (Mesh_Ptr const &m : ranks)
            all.push_back(m.get());

          std::cout << "\nCross-rank check\n";
          Ume::Part::Check_Result const res =
              Ume::Part::check_partitions(all, std::cout);
          std::cout << "  " << res.checks - res.failed_checks << '/'
                    << res.checks << " checks passed";
          if (!res)
            std::cout << ", " << res.problems << " problems";
          std::cout << std::endl;
          if (!res)
            status = 1;
        } else if (numpe > 1) {
          std::cout << "\nCross-rank check skipped: " << ranks.size() << " of "
                    << numpe << " ranks were given" << std::endl;
        }
      } else {
        [[maybe_unused]] auto const &test =
            ranks[0]->ds->caccess_vec3v("corner_csurf");
        [[maybe_unused]] auto const &test2 =
            ranks[0]->ds->caccess_vec3v("side_surz");
        [[maybe_unused]] auto const &test3 =
            ranks[0]->ds->caccess_vec3v("point_norm");
      }
    }
  }
  Ume::finalize();
  return status;
}

/*! A Mesh cannot be moved (see the note on Mesh's deleted move constructor),
    so the meshes are held by pointer: sorting the ranks below then permutes
    pointers and leaves every Mesh at the address it was constructed at. */
std::vector<Mesh_Ptr> read_meshes(int const argc, char *argv[]) {
  std::vector<Mesh_Ptr> ranks;
  bool need_sort{false};
  /* No argc check here: main() now counts the arguments after
     Ume::initialize() has stripped the --kokkos-* ones, so what arrives is
     what it saw. */
  for (int i = 1; i < argc; ++i) {
    std::cout << "Reading: " << argv[i] << '\n';
    std::ifstream is(argv[i]);
    if (!is) {
      std::cerr << "Unable to open file \"" << argv[i] << "\" for reading."
                << std::endl;
      return std::vector<Mesh_Ptr>{};
    }
    ranks.push_back(std::make_unique<Mesh>());
    ranks.back()->read(is);
    if (!is) {
      std::cerr << "Unable to read a mesh from file \"" << argv[i] << "\"."
                << std::endl;
      return std::vector<Mesh_Ptr>{};
    }
    if (ranks.back()->mype != i - 1)
      need_sort = true;
  }
  size_t const numpe = static_cast<size_t>(ranks[0]->numpe);
  if (numpe != ranks.size()) {
    std::cerr << "Warning: the initial mesh had " << numpe
              << " ranks, but only " << ranks.size() << " were read"
              << std::endl;
  }
  if (need_sort) {
    std::cerr << "Warning: sorting input ranks" << std::endl;
    std::sort(ranks.begin(), ranks.end(),
        [](Mesh_Ptr const &a, Mesh_Ptr const &b) { return a->mype < b->mype; });
  }

  return ranks;
}
