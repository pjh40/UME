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
  \file write_zone_mask_mesh.cc

  Writes a single-rank Ume binary mesh that has local zones with the masks
  given on the command line and no other entities:

      write_zone_mask_mesh <output.ume> [<zone mask>...]

  With no masks the mesh is empty, as a rank of an over-decomposed mesh is.
  The ume_mpi tests use these to place the interior zones (mask >= 1) where the
  driver's choice of a zone to seed its gradient has to find them.
*/

#include "Ume/SOA_Idx_Mesh.hh"
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <span>
#include <string_view>
#include <system_error>
#include <vector>

namespace {

/* The zone mask in `text`, or nothing if it is not entirely a `short`. */
std::optional<short> parse_mask(std::string_view const text) {
  short value{};
  auto const [end, ec] =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (ec != std::errc{} || end != text.data() + text.size())
    return std::nullopt;
  return value;
}

} // namespace

int main(int argc, char *argv[]) {
  std::span<char *const> const args(argv, static_cast<std::size_t>(argc));
  if (args.size() < 2) {
    std::cerr << "Usage: write_zone_mask_mesh <output.ume> [<zone mask>...]\n";
    return EXIT_FAILURE;
  }

  std::vector<short> masks;
  for (std::string_view const text : args.subspan(2)) {
    auto const mask = parse_mask(text);
    if (!mask) {
      std::cerr << "write_zone_mask_mesh: \"" << text
                << "\" is not a zone mask\n";
      return EXIT_FAILURE;
    }
    masks.push_back(*mask);
  }

  Ume::SOA_Idx::Mesh mesh;
  mesh.ivtag = UME_VERSION_2;
  mesh.mype = 0;
  mesh.numpe = 1;
  mesh.geo = Ume::SOA_Idx::Mesh::CARTESIAN;
  mesh.dump_iotas = false;
  Ume::SOA_Idx::Entity *const empty[] = {&mesh.points, &mesh.edges, &mesh.faces,
      &mesh.sides, &mesh.corners, &mesh.iotas};
  for (Ume::SOA_Idx::Entity *const e : empty)
    e->resize(0, 0, 0);
  int const nzones = static_cast<int>(masks.size());
  mesh.zones.resize(nzones, nzones, 0);
  mesh.zones.mask = masks;

  std::ofstream os(args[1], std::ios::binary);
  if (!os) {
    std::cerr << "write_zone_mask_mesh: cannot open \"" << args[1]
              << "\" for writing\n";
    return EXIT_FAILURE;
  }
  mesh.write(os);
  os.close();
  if (!os) {
    std::cerr << "write_zone_mask_mesh: failed writing \"" << args[1] << "\"\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
