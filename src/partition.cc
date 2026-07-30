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
  \file partition.cc

  Split a serial `.ume` mesh into N per-rank `.ume` files.  See
  `Ume/Partition.hh` for the strategies and `Ume/Decompose.hh` for what the
  split actually means.
*/

#include "Ume/Decompose.hh"
#include "Ume/Global_Ids.hh"
#include "Ume/Partition.hh"
#include "Ume/Partition_Metrics.hh"
#include "Ume/Partition_Reorder.hh"
#include "Ume/mesh_audit.hh"
#include "Ume/partition_check.hh"
#include "Ume/process_mgmt.hh"
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using namespace Ume::SOA_Idx;
namespace Part = Ume::Part;

namespace {

void usage(std::ostream &os) {
  os << R"(Usage: ume_partition [options] --ranks N <input.ume> <output-base>

Reads a serial .ume mesh and writes one <output-base>.NNNNN.ume per rank.

Options:
  --ranks N            Number of partitions (required)
  --method NAME        Partitioning method [slab]
  --list-methods       Print the methods this build knows about, and exit
  --vertex-weight W    unit | corners | sides [corners]
  --edge-weight W      unit | face-points | face-area [face-points].  Weighs
                       faces, so the kahypar* methods, which weigh points, do
                       not read it
  --imbalance F        Allowed heaviest-part to average-part ratio [1.03]
  --seed N             Seed for methods that make random choices.  kahypar and
                       kahypar-quality refine in parallel and are not
                       reproducible whatever it is; kahypar-det is
  --ghosts D           none | face | point [face]
  --owner R            low | balanced | hash [low]
  --topology SPEC      Renumber the parts onto a machine described as
                       x-separated group sizes, leaves first: 4 is four ranks
                       per node, 4x2 four per node and two nodes per switch.
                       Heavy neighbors are moved close together.  The parts
                       themselves are untouched; only which rank each is
                       given
  --zone-map FILE      Read a zone->part map instead of running a method
  --dump-zone-map FILE Write the zone->part map the method produced
  --metrics            Report what the decomposition costs
  --metrics-csv FILE   Append one row of metrics to FILE, with a header if it
                       is new, so that a sweep builds a plottable table
  --comm-matrix FILE   Write the rank-to-rank element counts to FILE, as
                       --ranks comma-separated rows
  --global-ids         Also write <output-base>.NNNNN.gid, saying what each
                       local entity was called in the input mesh.  A .ume file
                       has nowhere to record that, and comparing a partitioned
                       run against a serial one needs it.
  --audit              Run mesh_audit on each partition before writing
  --check              Check the partitions against each other afterwards.
                       Unlike --audit, this holds every partition in memory.
  --no-write           Do everything except write the files
  -h, --help           This message
)";
}

template <typename T>
void parse(std::string_view text, T &out, char const *what) {
  auto const last = text.data() + text.size();
  auto const [ptr, ec] = std::from_chars(text.data(), last, out);
  if (ec != std::errc{} || ptr != last) {
    std::cerr << "ume_partition: bad value for " << what << ": \"" << text
              << "\"\n";
    std::exit(2);
  }
}

//! Print every method, saying which of them this build can run
void list_methods(std::ostream &os) {
  for (Part::Method const &m : Part::all_methods()) {
    os << "  " << m.name;
    for (std::size_t i = std::char_traits<char>::length(m.name); i < 20; ++i)
      os << ' ';
    os << m.describe;
    if (!m.available)
      os << " [unavailable: needs " << m.requires_pkg << ']';
    os << '\n';
  }
}

bool read_zone_map(std::string const &path, Part::Zone_Map &z2r) {
  std::ifstream is(path);
  if (!is) {
    std::cerr << "ume_partition: cannot open \"" << path << "\"\n";
    return false;
  }
  z2r.clear();
  for (int part{}; is >> part;)
    z2r.push_back(part);
  return true;
}

bool write_zone_map(std::string const &path, Part::Zone_Map const &z2r) {
  std::ofstream os(path);
  if (!os) {
    std::cerr << "ume_partition: cannot open \"" << path << "\" for writing\n";
    return false;
  }
  for (int const part : z2r)
    os << part << '\n';
  return static_cast<bool>(os);
}

} // namespace

int main(int argc, char *argv[]) {
  Part::Options opt;
  Part::Decompose_Options dopt;
  std::string method{"slab"}, topology, zone_map_in, zone_map_out;
  std::string metrics_csv, comm_matrix_file;
  bool do_audit{false}, do_check{false}, do_write{true}, have_ranks{false};
  bool do_metrics{false}, do_global_ids{false}, have_topology{false};
  std::vector<std::string_view> positional;

  for (int i = 1; i < argc; ++i) {
    std::string_view const arg{argv[i]};
    auto next = [&](char const *what) -> std::string_view {
      if (i + 1 >= argc) {
        std::cerr << "ume_partition: " << what << " requires an argument\n";
        std::exit(2);
      }
      return std::string_view{argv[++i]};
    };
    auto choice = [&](std::string_view value, char const *what,
                      std::initializer_list<std::pair<char const *, int>> opts) {
      for (auto const &[name, v] : opts)
        if (value == name)
          return v;
      std::cerr << "ume_partition: " << what << " must be one of";
      for (auto const &[name, v] : opts)
        std::cerr << ' ' << name;
      std::cerr << '\n';
      std::exit(2);
    };

    if (arg == "-h" || arg == "--help") {
      usage(std::cout);
      return 0;
    } else if (arg == "--list-methods") {
      list_methods(std::cout);
      return 0;
    } else if (arg == "--ranks") {
      parse(next("--ranks"), opt.nparts, "--ranks");
      have_ranks = true;
    } else if (arg == "--method") {
      method = next("--method");
    } else if (arg == "--seed") {
      parse(next("--seed"), opt.seed, "--seed");
    } else if (arg == "--imbalance") {
      parse(next("--imbalance"), opt.imbalance, "--imbalance");
    } else if (arg == "--vertex-weight") {
      opt.vwgt = static_cast<Part::Vertex_Weight>(
          choice(next("--vertex-weight"), "--vertex-weight",
              {{"unit", static_cast<int>(Part::Vertex_Weight::UNIT)},
                  {"corners", static_cast<int>(Part::Vertex_Weight::CORNERS)},
                  {"sides", static_cast<int>(Part::Vertex_Weight::SIDES)}}));
    } else if (arg == "--edge-weight") {
      opt.ewgt = static_cast<Part::Edge_Weight>(choice(next("--edge-weight"),
          "--edge-weight",
          {{"unit", static_cast<int>(Part::Edge_Weight::UNIT)},
              {"face-points", static_cast<int>(Part::Edge_Weight::FACE_POINTS)},
              {"face-area", static_cast<int>(Part::Edge_Weight::FACE_AREA)}}));
    } else if (arg == "--ghosts") {
      dopt.ghosts = static_cast<Part::Ghost_Depth>(choice(next("--ghosts"),
          "--ghosts",
          {{"none", static_cast<int>(Part::Ghost_Depth::NONE)},
              {"face", static_cast<int>(Part::Ghost_Depth::FACE_ADJACENT)},
              {"point", static_cast<int>(Part::Ghost_Depth::POINT_ADJACENT)}}));
    } else if (arg == "--owner") {
      dopt.owner_rule = static_cast<Part::Owner_Rule>(choice(next("--owner"),
          "--owner",
          {{"low", static_cast<int>(Part::Owner_Rule::LOWEST_RANK)},
              {"balanced", static_cast<int>(Part::Owner_Rule::LEAST_LOADED)},
              {"hash", static_cast<int>(Part::Owner_Rule::HASH)}}));
    } else if (arg == "--topology") {
      topology = next("--topology");
      have_topology = true;
    } else if (arg == "--zone-map") {
      zone_map_in = next("--zone-map");
    } else if (arg == "--dump-zone-map") {
      zone_map_out = next("--dump-zone-map");
    } else if (arg == "--metrics") {
      do_metrics = true;
    } else if (arg == "--metrics-csv") {
      metrics_csv = next("--metrics-csv");
    } else if (arg == "--comm-matrix") {
      comm_matrix_file = next("--comm-matrix");
    } else if (arg == "--global-ids") {
      do_global_ids = true;
    } else if (arg == "--audit") {
      do_audit = true;
    } else if (arg == "--check") {
      do_check = true;
    } else if (arg == "--no-write") {
      do_write = false;
    } else if (arg.starts_with("-") && arg.size() > 1) {
      std::cerr << "ume_partition: unknown option \"" << arg << "\"\n";
      return 2;
    } else {
      positional.push_back(arg);
    }
  }

  if (!have_ranks || opt.nparts < 1) {
    std::cerr << "ume_partition: --ranks N is required, with N positive\n";
    return 2;
  }
  if (positional.size() != (do_write ? 2u : 1u)) {
    usage(std::cerr);
    return 2;
  }
  /* A side-car names the partition it belongs to, so there is nothing for it
     to sit beside if the partitions are not written. */
  if (do_global_ids && !do_write) {
    std::cerr << "ume_partition: --global-ids and --no-write disagree\n";
    return 2;
  }
  std::string const infile{positional[0]};
  std::string const outbase{do_write ? std::string{positional[1]} : ""};

  /* An unavailable method and a misspelled one are different mistakes, and
     saying which is the whole reason the table carries both. */
  if (zone_map_in.empty()) {
    Part::Method const *const m = Part::find_method(method);
    if (m == nullptr) {
      std::cerr << "ume_partition: unknown method \"" << method
                << "\".  Available methods:\n";
      list_methods(std::cerr);
      return 2;
    }
    if (!m->available) {
      std::cerr << "ume_partition: method \"" << method << "\" requires "
                << m->requires_pkg << ", which this build was configured "
                << "without\n";
      return 2;
    }
  }

  /* A malformed topology is a typing mistake like a misspelled method, so it
     is caught here rather than after a mesh has been read. */
  Part::Machine machine;
  if (have_topology) {
    try {
      machine = Part::parse_machine(topology, opt.nparts);
    } catch (std::invalid_argument const &e) {
      std::cerr << "ume_partition: --topology: " << e.what() << '\n';
      return 2;
    }
  }

  /* The centroids the geometric methods run on are a Kokkos kernel, so Kokkos
     has to be up before anything touches the Datastore. */
  Ume::initialize(argc, argv);
  int status = 0;
  {
    Mesh mesh;
    Ume::attach_serial_transport(mesh);
    {
      std::ifstream is(infile);
      if (!is) {
        std::cerr << "ume_partition: cannot open \"" << infile << "\"\n";
        Ume::finalize();
        return 3;
      }
      mesh.read(is);
    }
    std::cout << "Read \"" << infile << "\"\n";
    mesh.print_stats(std::cout);

    Part::Zone_Map z2r;
    if (!zone_map_in.empty()) {
      if (!read_zone_map(zone_map_in, z2r)) {
        Ume::finalize();
        return 3;
      }
      std::cout << "Zone map from \"" << zone_map_in << "\"\n";
    } else {
      auto const partitioner = Part::make(method);
      if (!partitioner) {
        /* The name is in the table and marked available, but no factory
           builds it -- a table entry somebody forgot to implement.
           Part::make() documents that it comes back as nullptr for this
           driver to report rather than be dereferenced. */
        std::cerr << "ume_partition: method \"" << method
                  << "\" is in the method table, but this build has no "
                     "implementation for it\n";
        Ume::finalize();
        return 2;
      }
      partitioner->assign(mesh, opt, z2r);
      std::cout << "Method " << partitioner->name() << ": "
                << partitioner->describe() << '\n';
    }
    /* Ahead of everything downstream, so that the map which gets dumped,
       decomposed, measured and written is the renamed one.  It costs two
       decompositions of its own: what the ranks send each other is what
       decides where they should sit, and the placement is measured after
       the fact rather than predicted.  See Partition_Reorder.hh. */
    if (have_topology) {
      Part::Rank_Order const placed =
          Part::place_ranks(mesh, z2r, opt.nparts, dopt, machine);
      Part::apply_rank_order(placed.order, z2r);

      std::cout << "Topology " << topology << ": " << placed.before
                << " element-hops per exchange sweep becomes " << placed.after;
      if (placed.after == placed.before)
        std::cout << " -- the numbering the method chose is kept";
      std::cout << '\n';
    }

    if (!zone_map_out.empty() && !write_zone_map(zone_map_out, z2r)) {
      Ume::finalize();
      return 3;
    }

    Part::Decomposition const d =
        Part::build_decomposition(mesh, std::move(z2r), opt.nparts, dopt);

    int empty = 0;
    for (int pe = 0; pe < d.nparts; ++pe)
      if (d.zones.local_size(pe) == 0)
        empty += 1;
    if (empty > 0)
      std::cerr << "ume_partition: warning: " << empty << " of " << d.nparts
                << " parts have no zones; the mesh has fewer zones than that\n";

    /* Metrics come out of the decomposition alone, so they are available
       before -- and with --no-write, instead of -- building any partition. */
    if (do_metrics || !metrics_csv.empty() || !comm_matrix_file.empty()) {
      Part::Metrics const met = Part::measure(mesh, d);
      if (do_metrics)
        Part::report(met, std::cout);

      if (!metrics_csv.empty()) {
        /* Appending a row per (method, ranks) is the point, so the header is
           written only when the file is not already a table. */
        std::error_code ec;
        auto const size = std::filesystem::file_size(metrics_csv, ec);
        bool const fresh = static_cast<bool>(ec) || size == 0;
        std::ofstream os(metrics_csv, std::ios::app);
        if (!os) {
          std::cerr << "ume_partition: cannot open \"" << metrics_csv
                    << "\" for writing\n";
          Ume::finalize();
          return 3;
        }
        if (fresh)
          Part::report_csv_header(os);
        Part::report_csv(met, zone_map_in.empty() ? method : "zone-map", os);
      }

      if (!comm_matrix_file.empty()) {
        std::ofstream os(comm_matrix_file);
        if (!os) {
          std::cerr << "ume_partition: cannot open \"" << comm_matrix_file
                    << "\" for writing\n";
          Ume::finalize();
          return 3;
        }
        Part::write_comm_matrix(met, os);
      }
    }

    /* One Mesh alive at a time, unless --check asks to compare them: the
       layouts already know everything about the other ranks that building
       this one needs.  A Mesh cannot be moved, so they are held by pointer. */
    std::vector<std::unique_ptr<Mesh>> kept;
    for (int pe = 0; pe < d.nparts; ++pe) {
      auto held = std::make_unique<Mesh>();
      Mesh &part = *held;
      Ume::attach_serial_transport(part);
      Part::extract(mesh, d, pe, part);

      std::cout << "  rank " << pe << ": " << part.zones.local_size() << '+'
                << part.zones.ghost_local_size() << " zones, "
                << part.points.local_size() << " points, "
                << part.sides.local_size() << '+'
                << part.sides.ghost_local_size() << " sides\n";

      if (do_audit) {
        Audit_Result const res = audit(part, std::cout, Audit_Options{});
        if (!res.ok()) {
          std::cerr << "ume_partition: rank " << pe << " failed the audit\n";
          status = 1;
        }
      }
      if (do_check)
        kept.push_back(std::move(held));

      if (do_write && status == 0) {
        std::string const path = Part::rank_file(outbase, pe, "ume");
        std::ofstream os(path);
        if (!os) {
          std::cerr << "ume_partition: cannot open \"" << path
                    << "\" for writing\n";
          status = 3;
          break;
        }
        part.write(os);
        os.close();

        /* A write/read round trip is nearly free, and Mesh::operator== covers
           everything that lands in the file. */
        Mesh check;
        std::ifstream is(path);
        check.read(is);
        if (!(part == check)) {
          std::cerr << "ume_partition: write/read round trip failed for rank "
                    << pe << '\n';
          status = 4;
          break;
        }

        /* The side-car comes out of the decomposition rather than the mesh:
           the mesh no longer knows what its entities used to be called, which
           is the whole reason this file exists. */
        if (do_global_ids) {
          std::string const gpath = Part::rank_file(outbase, pe, "gid");
          std::ofstream gos(gpath);
          if (!gos) {
            std::cerr << "ume_partition: cannot open \"" << gpath
                      << "\" for writing\n";
            status = 3;
            break;
          }
          Part::write_global_ids(d, pe, gos);
        }
      }
    }

    /* The audit above sees one partition and can only ask whether its
       bookkeeping is well formed.  This is the check that the ranks agree
       with each other, and that between them they are still the input. */
    if (do_check && status == 0) {
      std::vector<Mesh const *> ranks;
      ranks.reserve(kept.size());
      for (auto const &m : kept)
        ranks.push_back(m.get());

      std::cout << "Cross-rank check\n";
      Part::Check_Result const res =
          Part::check_partitions(mesh, ranks, std::cout);
      std::cout << "  " << res.checks - res.failed_checks << '/' << res.checks
                << " checks passed";
      if (!res)
        std::cout << ", " << res.problems << " problems";
      std::cout << std::endl;
      if (!res)
        status = 5;
    }
  }
  Ume::finalize();
  return status;
}
