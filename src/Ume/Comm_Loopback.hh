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
  \file Ume/Comm_Loopback.hh

  Every rank of a decomposition, communicating inside one process.

  ### Why this exists

  The strongest thing that can be said about a decomposition is that it still
  computes the same answer, and until now saying it needed `USE_MPI=ON`, an
  `mpiexec`, and files on disk (see `partition_mpi_compare` in
  `test/CMakeLists.txt`).  A `Loopback::Group` says it in-process, so the test
  runs in every build.

  ### Why threads, and why only one of them runs at a time

  `Transport::exchange` is called synchronously from the middle of
  `Entity::gather`, which is itself called from the middle of a kernel: a rank
  cannot hand over its send buffer and come back for the answer later.  So each
  rank needs its own stack, and that means a thread.

  But only *one* of those threads may execute at a time.  Kokkos parallel
  dispatch is not thread-safe -- the default execution space here is OpenMP,
  whose backend instance carries shared scratch -- and `gradzatz` dispatches on
  either side of the `gathscat` in the middle of it.  So `Group` is a
  cooperative scheduler rather than a parallel one: exactly one rank holds the
  token, and `exchange` is where it changes hands.  A rank posts its buffers,
  passes the token to the next rank that has not yet reached this exchange, and
  blocks; the rank that completes the set performs every copy and hands the
  token back.  That is a correctness harness, not a performance model, and it
  is all this is for.

  Two consequences worth knowing:

  - **`exchange` is collective**, exactly as it is under MPI: every rank must
    make the same sequence of calls.  A rank that returns early, or that throws
    while the others are waiting, is reported rather than allowed to hang the
    run -- `run()` throws instead.
  - **Assertions belong on the main thread.**  Bodies run one at a time, but
    they are still not the thread `run()` was called on; record per-rank
    results and check them after `run()` returns.
*/

#ifndef UME_COMM_LOOPBACK_HH
#define UME_COMM_LOOPBACK_HH 1

#include "Ume/Comm_Transport.hh"
#include <functional>
#include <memory>
#include <vector>

namespace Ume {
namespace Comm {

//! One rank's view of a `Loopback::Group`
class Loopback : public Transport {
public:
  class Group;

  void exchange(Buffers<DS_Types::INTV_T> const &sends,
      Buffers<DS_Types::INTV_T> &recvs) override;
  void exchange(Buffers<DS_Types::DBLV_T> const &sends,
      Buffers<DS_Types::DBLV_T> &recvs) override;
  void exchange(Buffers<DS_Types::VEC3V_T> const &sends,
      Buffers<DS_Types::VEC3V_T> &recvs) override;

  int id() const override { return rank_; }
  //! The exchanges above really move data between the group's ranks
  bool does_exchanges() const override { return true; }
  //! Nothing to shut down; the Group outlives every rank that uses it
  int stop() override { return 0; }

  //! Throw rather than `std::abort()`, so one rank cannot take down the suite
  /*! The other ranks are waiting in an `exchange` this rank will now never
      reach.  Unwinding out of the body is what tells `Group` to release them,
      and `run()` rethrows this on the calling thread. */
  [[noreturn]] void abort(char const *const message) override;

private:
  friend class Group;
  Loopback(Group &group, int rank) : group_{&group}, rank_{rank} {}

  //! Post this rank's buffers, yield the token, and return once they are filled
  void rendezvous(int tag, void const *sends, void *recvs);

  Group *group_;
  int rank_;
};

//! The whole decomposition: one `Transport` per rank, and the scheduler
class Loopback::Group {
public:
  //! Throws `std::invalid_argument` unless `nranks >= 1`
  explicit Group(int nranks);
  ~Group();
  Group(Group const &) = delete;
  Group &operator=(Group const &) = delete;

  int nranks() const { return static_cast<int>(ranks_.size()); }

  //! Rank `r`'s Transport, to be attached as that rank's `Mesh::comm`
  /*! The reference is stable for the lifetime of the Group.  Throws
      `std::out_of_range` if `r` is not a rank of this Group. */
  Transport &operator[](int r);

  //! Run `body(r)` for every rank, and return once all of them have finished
  /*! Bodies run one at a time and interleave only where they communicate.

      Throws whatever a body threw -- the lowest-numbered rank's, if more than
      one did -- after every thread has been joined.  Throws
      `std::runtime_error` if the ranks disagreed about what to exchange, and
      `std::logic_error` if called from inside a body. */
  void run(std::function<void(int)> const &body);

private:
  friend class Loopback;
  //! The scheduler.  Defined in the .cc: nothing outside needs its shape.
  struct State;

  std::vector<std::unique_ptr<Loopback>> ranks_;
  std::unique_ptr<State> state_;
};

} // namespace Comm
} // namespace Ume

#endif
