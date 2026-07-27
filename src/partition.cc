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
#include "Ume/Partition.hh"
#include "Ume/mesh_audit.hh"
#include "Ume/process_mgmt.hh"
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
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
  --edge-weight W      unit | face-points | face-area [face-points]
  --imbalance F        Allowed heaviest-part to average-part ratio [1.03]
  --seed N             Seed for methods that make random choices
  --ghosts D           none | face | point [face]
  --owner R            low | balanced | hash [low]
  --zone-map FILE      Read a zone->part map instead of running a method
  --dump-zone-map FILE Write the zone->part map the method produced
  --audit              Run mesh_audit on each partition before writing
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

//! `<base>.NNNNN.ume`, the naming `txt2bin` already produces
std::string rank_file(std::string const &base, int pe) {
  char suffix[16];
  std::snprintf(suffix, sizeof(suffix), ".%05d.ume", pe);
  return base + suffix;
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
  std::string method{"slab"}, zone_map_in, zone_map_out;
  bool do_audit{false}, do_write{true}, have_ranks{false};
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
    } else if (arg == "--zone-map") {
      zone_map_in = next("--zone-map");
    } else if (arg == "--dump-zone-map") {
      zone_map_out = next("--dump-zone-map");
    } else if (arg == "--audit") {
      do_audit = true;
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
      partitioner->assign(mesh, opt, z2r);
      std::cout << "Method " << partitioner->name() << ": "
                << partitioner->describe() << '\n';
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

    /* One Mesh alive at a time: the layouts already know everything about the
       other ranks that building this one needs. */
    for (int pe = 0; pe < d.nparts; ++pe) {
      Mesh part;
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

      if (do_write && status == 0) {
        std::string const path = rank_file(outbase, pe);
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
      }
    }
  }
  Ume::finalize();
  return status;
}
