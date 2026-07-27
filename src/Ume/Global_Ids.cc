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
  \file Ume/Global_Ids.cc
*/

#include "Ume/Global_Ids.hh"
#include <cstdio>
#include <istream>
#include <ostream>
#include <stdexcept>

namespace Ume {
namespace Part {

namespace {

//! How many indices go on one line
constexpr int PER_LINE = 10;

[[noreturn]] void bad(std::string const &what) {
  throw std::runtime_error("global-id side-car: " + what);
}

//! Read `name` and the value that must follow it
int keyword(std::istream &is, char const *name) {
  std::string tok;
  int value{};
  if (!(is >> tok >> value) || tok != name)
    bad(std::string{"expected \""} + name + '"');
  return value;
}

} // namespace

void write_global_ids(Decomposition const &d, int pe, std::ostream &os) {
  if (pe < 0 || pe >= d.nparts)
    throw std::out_of_range("write_global_ids: no such part");

  os << "ume_gid 1\n";
  os << "mype " << pe << '\n';
  os << "numpe " << d.nparts << '\n';

  auto const layouts = d.layouts();
  for (std::size_t k = 0; k < NKINDS; ++k) {
    Entity_Layout const &l = *layouts[k];
    int const nlocal = l.local_size(pe), nghost = l.ghost_local_size(pe);
    os << kind_name(static_cast<int>(k)) << ' ' << nlocal << ' ' << nghost
       << '\n';
    for (int i = 0; i < nlocal + nghost; ++i)
      os << l.global_index(pe, i) << ((i + 1) % PER_LINE == 0 ? '\n' : ' ');
    if ((nlocal + nghost) % PER_LINE != 0)
      os << '\n';
  }
}

Global_Ids read_global_ids(std::istream &is) {
  Global_Ids g;
  if (keyword(is, "ume_gid") != 1)
    bad("only version 1 is understood");
  g.mype = keyword(is, "mype");
  g.numpe = keyword(is, "numpe");

  for (std::size_t k = 0; k < NKINDS; ++k) {
    std::string tok;
    int nlocal{}, nghost{};
    char const *const name = kind_name(static_cast<int>(k));
    if (!(is >> tok >> nlocal >> nghost) || tok != name)
      bad(std::string{"expected the "} + name + " section");
    if (nlocal < 0 || nghost < 0)
      bad(std::string{"negative "} + name + " count");

    Kind_Ids &ids = g.kinds[k];
    ids.local = nlocal;
    ids.gid.resize(static_cast<std::size_t>(nlocal + nghost));
    for (int &v : ids.gid)
      if (!(is >> v))
        bad(std::string{"the "} + name + " section is short");
  }
  return g;
}

std::string rank_file(std::string_view base, int pe, std::string_view ext) {
  char suffix[16];
  std::snprintf(suffix, sizeof(suffix), ".%05d.", pe);
  return std::string{base} + suffix + std::string{ext};
}

} // namespace Part
} // namespace Ume
