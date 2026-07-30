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
  \file Ume/random.hh

  A seeded generator whose sequence is fixed here rather than by the standard
  library.

  Anything a `--seed` is supposed to reproduce must draw from these and not
  from `std::shuffle` or `std::uniform_int_distribution`: neither has a
  specified algorithm, so the same seed gives different results under
  libstdc++ and libc++, and a seed would only reproduce a run on the machine
  that made it.  Where a shuffle is wanted, write Fisher-Yates out -- it is a
  bijection onto the permutations however the draws come out, so nothing can
  be dropped or placed twice.
*/

#ifndef UME_RANDOM_HH
#define UME_RANDOM_HH 1

#include <cstdint>

namespace Ume {

//! One step of SplitMix64, advancing `state`
inline std::uint64_t next_random(std::uint64_t &state) {
  std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

//! A uniform draw from `[0, n)`, with `n > 0`
inline std::uint64_t bounded_random(std::uint64_t &state, std::uint64_t n) {
  /* 2^64 % n draws have to be thrown away for the rest to be uniform: taking
     the modulus of every draw would favour the first 2^64 % n values. */
  std::uint64_t const reject = (~n + 1u) % n;
  std::uint64_t x;
  do {
    x = next_random(state);
  } while (x < reject);
  return x % n;
}

} // namespace Ume

#endif
