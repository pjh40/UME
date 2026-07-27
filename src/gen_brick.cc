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
  \file gen_brick.cc

  Generate a synthetic hex/tet brick mesh and write it as a binary `.ume`
  file.  See `Ume/Brick_Mesh.hh` for the generator itself.
*/

#include "Ume/Brick_Mesh.hh"
#include "Ume/process_mgmt.hh"
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

using namespace Ume::SOA_Idx;
namespace Brick = Ume::Brick;

namespace {

void usage(std::ostream &os) {
  os << R"(Usage: gen_brick [options] <nx> <ny> <nz> <output.ume>

Generates a brick of nx*ny*nz cells.  Each cell is either a single hexahedral
zone or six tetrahedra (a Kuhn subdivision); faces shared with a subdivided
cell are triangulated on both sides, so hex zones next to tet cells become
general polyhedra.

Options:
  --tet-fraction F   Portion of the mesh that is tetrahedral [0]
  --basis zones      F is the fraction of resulting zones (default)
  --basis cells      F is the fraction of lattice cells
  --clustering L     Arrangement of the tet cells, 0 = uniform random scatter,
                     1 = a single face-connected component [0]
  --seed N           Seed for the cell selection [deterministic default]
  --spacing DX DY DZ Cell size [1 1 1]
  --origin X Y Z     Coordinates of the lattice origin [0 0 0]
  --phantom-sides    Represent the boundary with phantom sides only, instead
                     of degenerate boundary zones (see Boundary_Style)
  --validate         Check the mesh invariants before writing
  --no-write         Skip the output file (with --validate, a pure check run)
  -h, --help         This message
)";
}

//! Parse `text` into `out`, exiting with a message on failure
template <typename T>
void parse(std::string_view text, T &out, char const *what) {
  auto const last = text.data() + text.size();
  auto const [ptr, ec] = std::from_chars(text.data(), last, out);
  if (ec != std::errc{} || ptr != last) {
    std::cerr << "gen_brick: bad value for " << what << ": \"" << text
              << "\"\n";
    std::exit(2);
  }
}

} // namespace

int main(int argc, char *argv[]) {
  Brick::Spec spec;
  bool do_validate{false}, do_write{true};
  std::string outfile;
  std::vector<std::string_view> positional;

  for (int i = 1; i < argc; ++i) {
    std::string_view const arg{argv[i]};
    auto next = [&](char const *what) -> std::string_view {
      if (i + 1 >= argc) {
        std::cerr << "gen_brick: " << what << " requires an argument\n";
        std::exit(2);
      }
      return std::string_view{argv[++i]};
    };
    if (arg == "-h" || arg == "--help") {
      usage(std::cout);
      return 0;
    } else if (arg == "--tet-fraction") {
      parse(next("--tet-fraction"), spec.tet_fraction, "--tet-fraction");
    } else if (arg == "--clustering") {
      parse(next("--clustering"), spec.clustering, "--clustering");
    } else if (arg == "--seed") {
      parse(next("--seed"), spec.seed, "--seed");
    } else if (arg == "--basis") {
      std::string_view const b = next("--basis");
      if (b == "cells")
        spec.basis = Brick::Fraction_Basis::CELLS;
      else if (b == "zones")
        spec.basis = Brick::Fraction_Basis::ZONES;
      else {
        std::cerr << "gen_brick: --basis must be \"cells\" or \"zones\"\n";
        return 2;
      }
    } else if (arg == "--spacing") {
      parse(next("--spacing"), spec.dx, "--spacing");
      parse(next("--spacing"), spec.dy, "--spacing");
      parse(next("--spacing"), spec.dz, "--spacing");
    } else if (arg == "--origin") {
      parse(next("--origin"), spec.x0, "--origin");
      parse(next("--origin"), spec.y0, "--origin");
      parse(next("--origin"), spec.z0, "--origin");
    } else if (arg == "--phantom-sides") {
      spec.boundary = Brick::Boundary_Style::PHANTOM_SIDES_ONLY;
    } else if (arg == "--validate") {
      do_validate = true;
    } else if (arg == "--no-write") {
      do_write = false;
    } else if (arg.starts_with("-") && arg.size() > 1) {
      std::cerr << "gen_brick: unknown option \"" << arg << "\"\n";
      return 2;
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.size() != (do_write ? 4u : 3u)) {
    usage(std::cerr);
    return 2;
  }
  parse(positional[0], spec.nx, "nx");
  parse(positional[1], spec.ny, "ny");
  parse(positional[2], spec.nz, "nz");
  if (do_write)
    outfile = positional[3];

  if (spec.nx < 1 || spec.ny < 1 || spec.nz < 1) {
    std::cerr << "gen_brick: nx, ny and nz must all be positive\n";
    return 2;
  }

  /* The derived mesh variables run through Kokkos, so it has to be up before
     anything touches the Datastore. */
  Ume::initialize(argc, argv);
  int status = 0;
  {
    Mesh mesh;
    Ume::attach_serial_transport(mesh);
    Brick::Stats const stats = Brick::build(spec, mesh);

    std::cout << "Brick " << spec.nx << 'x' << spec.ny << 'x' << spec.nz
              << '\n';
    stats.print(std::cout);

    if (do_validate) {
      std::cout << "Validating..." << std::endl;
      if (Brick::validate(mesh, spec, std::cout)) {
        std::cout << "  all checks passed\n";
      } else {
        status = 1;
      }
    }

    if (do_write && status == 0) {
      std::cout << "Writing \"" << outfile << "\"" << std::endl;
      std::ofstream os(outfile);
      if (!os) {
        std::cerr << "gen_brick: cannot open \"" << outfile
                  << "\" for writing\n";
        status = 3;
      } else {
        mesh.write(os);
        os.close();
        /* A write/read round trip is nearly free, and Mesh::operator== is a
           complete check of everything that lands in the file. */
        Mesh check;
        std::ifstream is(outfile);
        check.read(is);
        if (!(mesh == check)) {
          std::cerr << "gen_brick: write/read round trip failed\n";
          status = 4;
        }
      }
    }
  }
  Ume::finalize();
  return status;
}
