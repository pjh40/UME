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
  \file test_comm_buffers.cc

  Tests for the pack/unpack layer that every gather and scatter runs through.
*/

#include "Ume/Comm_Buffers.hh"
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <numeric>

using namespace Ume;
using namespace Ume::Comm;

namespace {

//! Two remotes, so that the second one's offset into the buffer matters
Neighbors two_remotes() { return Neighbors{{1, {0, 2, 4}}, {3, {1, 5}}}; }

} // namespace

TEMPLATE_TEST_CASE("comm buffers: pack then unpack is the identity", "[comm]",
    DS_Types::INTV_T, DS_Types::DBLV_T, DS_Types::VEC3V_T) {
  using Field = TestType;
  using Elem = typename Field::value_type;

  Field field(8);
  for (std::size_t i = 0; i < field.size(); ++i)
    field[i] =
        Elem(static_cast<typename DS_Type_Info<Field>::base_type>(i + 1));
  Field const original = field;

  Buffers<Field> bufs{two_remotes()};
  REQUIRE(bufs.num_entries() == 5);
  bufs.pack(field);

  /* Scribble over the field so the unpack has to put the values back rather
     than merely leave them alone. */
  for (auto &e : field)
    e = Elem(0);
  // Positive control: the scribble must flip the comparison checked below.
  for (std::size_t i : {0u, 2u, 4u, 1u, 5u})
    REQUIRE_FALSE(field[i] == original[i]);
  bufs.unpack(field, Op::OVERWRITE);

  for (std::size_t i : {0u, 2u, 4u, 1u, 5u})
    CHECK(field[i] == original[i]);
}

TEST_CASE(
    "comm buffers: every component of a Vec3 gets its own slot", "[comm]") {
  /* `pack` has to advance the buffer cursor by the whole element.  When it did
     not, every Vec3 landed on top of the first one, so a buffer holding N
     points held one point's coordinates and N-1 zeros -- which is invisible
     until two entities are exchanged with the same remote. */
  DS_Types::VEC3V_T field(8, Vec3(0.0));
  field[0] = Vec3({1.0, 2.0, 3.0});
  field[2] = Vec3({4.0, 5.0, 6.0});
  field[4] = Vec3({7.0, 8.0, 9.0});

  Buffers<DS_Types::VEC3V_T> bufs{Neighbors{{1, {0, 2, 4}}}};
  bufs.pack(field);

  double const *const buf = bufs.get_buf();
  std::vector<double> const packed(buf, buf + 9);
  CHECK(packed == std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9});
}
