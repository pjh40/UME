/*
  Copyright (c) 2023, Triad National Security, LLC. All rights reserved.

  This is open source software; you can redistribute it and/or modify it under
  the terms of the BSD-3 License. If software is modified to produce derivative
  works, such modified software should be clearly marked, so as not to confuse
  it with the version available from LANL. Full text of the BSD-3 License can be
  found in the LICENSE.md file, and the full assertion of copyright in the
  NOTICE.md file.
*/

#include "Ume/VecN.hh"
#include "Ume/utils.hh"
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <string>
#include <vector>

template <class T> void wr(T const &in_val, T &out_val) {
  std::stringstream filebuf;
  Ume::write_bin(filebuf, in_val);
  Ume::read_bin(filebuf, out_val);
}

TEST_CASE("size_t binary r/w", "[utils]") {
  const size_t in_size_t{1234567890123};
  size_t out_size_t;
  wr(in_size_t, out_size_t);
  REQUIRE(in_size_t == out_size_t);
}

TEST_CASE("int binary r/w", "[utils]") {
  const int in_int{-1234};
  int out_int;
  wr(in_int, out_int);
  REQUIRE(in_int == out_int);
}

TEST_CASE("std::string binary r/w", "[utils]") {
  const std::string in_str{"  This is a test  "};
  std::string out_str;
  wr(in_str, out_str);
  REQUIRE(in_str == out_str);
}

TEST_CASE("vector<int> binary r/w", "[utils]") {
  const std::vector<int> in_vec{-1, 1, -2, 2, -4, 4};
  std::vector<int> out_vec;
  wr(in_vec, out_vec);
  REQUIRE(in_vec == out_vec);
}

TEST_CASE("vector<short> binary r/w", "[utils]") {
  const std::vector<short> in_vec{-1, 1, -2, 2, -4, 4};
  std::vector<short> out_vec;
  wr(in_vec, out_vec);
  REQUIRE(in_vec == out_vec);
}

TEST_CASE("vector<Vec3> binary r/w", "[utils]") {
  using Ume::Vec3;
  std::vector<Vec3> in_vec{2}, out_vec;
  in_vec[0] = Vec3({1.0, 1.0e-12, 10});
  in_vec[1] = Vec3({-5.0, -5.0e-12, -50});
  wr(in_vec, out_vec);
  REQUIRE(in_vec == out_vec);
}

namespace {

/* Leaves a recognizable non-zero pattern in the stack below the caller, where
   the next call's locals land: an uninitialized local in that call then holds
   the pattern rather than whatever zeros happened to be there.  At -O0, which
   is how the test trees build, this is reliable; optimized, it may not be. */
[[gnu::noinline]] void scribble_stack() {
  volatile unsigned char junk[4096];
  for (auto &byte : junk)
    byte = 0xA5;
}

} // namespace

/* read_bin's length had no initializer, so a stream that had already failed
   (a file that could not be opened) left it holding stack garbage, which then
   sized the allocation and the read. */
TEST_CASE("std::string binary read from an empty stream", "[utils]") {
  std::string out{"stale"};
  std::istringstream empty;
  scribble_stack();
  Ume::read_bin(empty, out);
  CHECK(empty.fail());
  CHECK(out.empty());

  // Positive control: a stream that holds a string reads it into `out`.
  std::stringstream full;
  Ume::write_bin(full, std::string{"fresh"});
  Ume::read_bin(full, out);
  CHECK(out == "fresh");
}

TEST_CASE("vector<int> binary read from an empty stream", "[utils]") {
  std::vector<int> out{1, 2, 3};
  std::istringstream empty;
  scribble_stack();
  Ume::read_bin(empty, out);
  CHECK(empty.fail());
  CHECK(out.empty());

  // Positive control: a stream that holds a vector reads it into `out`.
  std::stringstream full;
  Ume::write_bin(full, std::vector<int>{4, 5});
  Ume::read_bin(full, out);
  CHECK(out == std::vector<int>{4, 5});
}

TEST_CASE("ltrim", "[utils]") {
  std::string a{"  \tThis is a test\t    "};
  std::string b;
  b = Ume::ltrim(a);
  REQUIRE(b == "This is a test\t    ");
}

TEST_CASE("rtrim", "[utils]") {
  std::string a{"  \tThis is a test\t    "};
  std::string b;
  b = Ume::rtrim(a);
  REQUIRE(b == "  \tThis is a test");
}

TEST_CASE("trim", "[utils]") {
  std::string a{"  \tThis is a test\t    "};
  std::string b;
  b = Ume::trim(a);
  REQUIRE(b == "This is a test");
}

TEST_CASE("mesh_filename with a 200-character basename", "[utils]") {
  // The drivers formatted these names into a char[80]; a basename this long
  // overflowed it by more than a hundred bytes.
  std::string base;
  for (int i = 0; i < 49; ++i)
    base += "dir/";
  base += "mesh";
  REQUIRE(base.size() == 200);

  REQUIRE(Ume::mesh_filename(base, 7) == base + ".00007.ume");
  REQUIRE(Ume::mesh_filename(base, 7, 4) == base + ".00004.00007.ume");
  // Positive controls: the rank, and the order of rank and scale, show up in
  // the name.
  REQUIRE(Ume::mesh_filename(base, 8) != Ume::mesh_filename(base, 7));
  REQUIRE(Ume::mesh_filename(base, 4, 7) == base + ".00007.00004.ume");
  // A field wider than five digits widens rather than truncates, as %05d did.
  REQUIRE(Ume::mesh_filename(base, 123456) == base + ".123456.ume");
  REQUIRE(Ume::mesh_filename(base, 0, 131072) == base + ".131072.00000.ume");
}
