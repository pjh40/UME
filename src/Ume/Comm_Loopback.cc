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
  \file Ume/Comm_Loopback.cc
*/

#include "Ume/Comm_Loopback.hh"
#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

namespace Ume {
namespace Comm {

namespace {

//! Which `exchange` overload a rank posted
/*! The ranks are all inside one call of one kernel, so they must agree.  The
    tag is what turns a disagreement into a complaint instead of a
    `Buffers<double>` read as a `Buffers<int>`. */
enum Buffer_Tag { TAG_INT, TAG_DBL, TAG_VEC3 };

char const *tag_name(int const tag) {
  switch (tag) {
  case TAG_INT:
    return "int";
  case TAG_DBL:
    return "double";
  case TAG_VEC3:
    return "Vec3";
  default:
    return "unknown";
  }
}

} // namespace

//! The cooperative scheduler shared by every rank of one Group
struct Loopback::Group::State {
  //! What one rank handed to `exchange`, valid only while that rank is posted
  struct Post {
    int tag{-1};
    void const *sends{nullptr};
    void *recvs{nullptr};
  };

  std::mutex m;
  std::condition_variable cv;

  int nranks{0};
  int active{-1}; //!< the one rank allowed to run, or -1 if none is
  int nlive{0}; //!< ranks whose body has not returned
  int nposted{0}; //!< ranks waiting in the exchange now being assembled
  std::vector<char> live, posted;
  std::vector<Post> post;
  std::vector<std::exception_ptr> err;
  //! Every way the ranks failed to agree about what they were exchanging
  std::vector<std::string> complaints;

  /* Every one of these is called with `m` held. */

  //! Give the token to the next rank that has not reached the current exchange
  void hand_off();
  //! Fill every waiting rank's receive buffers, then end the round
  void finish_round();
  //! Copy into rank `r`'s receive buffers from the ranks that are sending to it
  template <class T> void deliver(int r);
};

void Loopback::Group::State::hand_off() {
  int const from = active;
  active = -1;
  /* Round-robin from the rank that just yielded, so a rank that keeps
     communicating cannot starve the others.  `from` is -1 before the first
     hand-off, which starts the search at rank 0. */
  for (int i = 1; i <= nranks; ++i) {
    int const r = (from + i) % nranks;
    if (live[static_cast<std::size_t>(r)] &&
        !posted[static_cast<std::size_t>(r)]) {
      active = r;
      break;
    }
  }
  cv.notify_all();
}

template <class T> void Loopback::Group::State::deliver(int const r) {
  auto &recvs =
      *static_cast<Buffers<T> *>(post[static_cast<std::size_t>(r)].recvs);
  for (auto const &in : recvs.remotes) {
    int const p = in.pe;
    auto complain = [&](std::string const &why) {
      std::ostringstream os;
      os << "rank " << r << " expects " << in.buf_len
         << " element(s) from rank " << p << ", but " << why;
      complaints.push_back(os.str());
    };

    if (p < 0 || p >= nranks || !posted[static_cast<std::size_t>(p)]) {
      complain("that rank is not in this exchange");
      continue;
    }
    auto const &sends = *static_cast<Buffers<T> const *>(
        post[static_cast<std::size_t>(p)].sends);
    auto const out = std::find_if(sends.remotes.begin(), sends.remotes.end(),
        [r](auto const &rem) { return rem.pe == r; });
    if (out == sends.remotes.end()) {
      complain("it is not sending to this one");
      continue;
    }
    if (out->buf_len != in.buf_len) {
      complain("it is sending " + std::to_string(out->buf_len));
      continue;
    }
    /* The send buffer is another rank's, and it is filled: that rank packed it
       before posting, and is blocked until this round is over. */
    std::copy_n(
        sends.buf.begin() + static_cast<std::ptrdiff_t>(out->buf_offset),
        in.buf_len,
        recvs.buf.begin() + static_cast<std::ptrdiff_t>(in.buf_offset));
  }
}

void Loopback::Group::State::finish_round() {
  int tag = -1;
  bool agreed = true;
  for (int r = 0; r < nranks; ++r) {
    if (!posted[static_cast<std::size_t>(r)])
      continue;
    int const rtag = post[static_cast<std::size_t>(r)].tag;
    if (tag < 0) {
      tag = rtag;
    } else if (rtag != tag) {
      std::ostringstream os;
      os << "rank " << r << " is exchanging " << tag_name(rtag)
         << " elements while another rank is exchanging " << tag_name(tag);
      complaints.push_back(os.str());
      agreed = false;
    }
  }

  /* Without agreement there is no safe way to read the other ranks' buffers,
     so nothing is delivered and every rank sees the values it arrived with. */
  if (agreed) {
    for (int r = 0; r < nranks; ++r) {
      if (!posted[static_cast<std::size_t>(r)])
        continue;
      switch (tag) {
      case TAG_INT:
        deliver<DS_Types::INTV_T>(r);
        break;
      case TAG_DBL:
        deliver<DS_Types::DBLV_T>(r);
        break;
      case TAG_VEC3:
        deliver<DS_Types::VEC3V_T>(r);
        break;
      default:
        break;
      }
    }
  }

  std::fill(posted.begin(), posted.end(), char{0});
  nposted = 0;
}

/* ------------------------------- one rank -------------------------------- */

void Loopback::rendezvous(int const tag, void const *sends, void *recvs) {
  Group::State &s = *group_->state_;
  std::unique_lock lk(s.m);
  if (s.active != rank_)
    throw std::logic_error(
        "Comm::Loopback::exchange: this rank is not running; "
        "a Loopback Transport is only usable from inside "
        "Loopback::Group::run");

  s.post[static_cast<std::size_t>(rank_)] =
      Group::State::Post{tag, sends, recvs};
  s.posted[static_cast<std::size_t>(rank_)] = 1;
  ++s.nposted;

  /* The rank that completes the set does every copy: it is the only one that
     knows the set is complete, and it can reach every other rank's buffers. */
  if (s.nposted == s.nlive)
    s.finish_round();

  s.hand_off();
  s.cv.wait(lk, [&s, this] { return s.active == rank_; });
}

void Loopback::exchange(
    Buffers<DS_Types::INTV_T> const &sends, Buffers<DS_Types::INTV_T> &recvs) {
  rendezvous(TAG_INT, &sends, &recvs);
}

void Loopback::exchange(
    Buffers<DS_Types::DBLV_T> const &sends, Buffers<DS_Types::DBLV_T> &recvs) {
  rendezvous(TAG_DBL, &sends, &recvs);
}

void Loopback::exchange(Buffers<DS_Types::VEC3V_T> const &sends,
    Buffers<DS_Types::VEC3V_T> &recvs) {
  rendezvous(TAG_VEC3, &sends, &recvs);
}

void Loopback::abort(char const *const message) {
  std::ostringstream os;
  os << "Comm::Loopback: rank " << rank_ << " aborted: " << message;
  throw std::runtime_error(os.str());
}

/* -------------------------------- the group ------------------------------- */

Loopback::Group::Group(int const nranks) : state_{std::make_unique<State>()} {
  if (nranks < 1)
    throw std::invalid_argument(
        "Comm::Loopback::Group: a group needs at least one rank");

  ranks_.reserve(static_cast<std::size_t>(nranks));
  for (int r = 0; r < nranks; ++r)
    ranks_.emplace_back(new Loopback(*this, r));

  State &s = *state_;
  s.nranks = nranks;
  s.live.assign(static_cast<std::size_t>(nranks), 0);
  s.posted.assign(static_cast<std::size_t>(nranks), 0);
  s.post.assign(static_cast<std::size_t>(nranks), State::Post{});
  s.err.assign(static_cast<std::size_t>(nranks), std::exception_ptr{});
}

Loopback::Group::~Group() = default;

Transport &Loopback::Group::operator[](int const r) {
  if (r < 0 || r >= nranks())
    throw std::out_of_range("Comm::Loopback::Group: no such rank");
  return *ranks_[static_cast<std::size_t>(r)];
}

void Loopback::Group::run(std::function<void(int)> const &body) {
  State &s = *state_;
  int const n = nranks();

  {
    std::lock_guard lk(s.m);
    if (s.nlive != 0)
      throw std::logic_error("Comm::Loopback::Group::run is not reentrant");
    std::fill(s.live.begin(), s.live.end(), char{1});
    std::fill(s.posted.begin(), s.posted.end(), char{0});
    std::fill(s.err.begin(), s.err.end(), std::exception_ptr{});
    s.complaints.clear();
    s.nlive = n;
    s.nposted = 0;
    s.active = 0;
  }

  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(n));
  for (int r = 0; r < n; ++r) {
    threads.emplace_back([&s, &body, r] {
      {
        std::unique_lock lk(s.m);
        s.cv.wait(lk, [&s, r] { return s.active == r; });
      }
      /* This rank holds the token from here until it hands it back below, so
         nothing else is running and `body` may touch process-wide state --
         Kokkos dispatch included -- as freely as a serial program would. */
      try {
        body(r);
      } catch (...) {
        s.err[static_cast<std::size_t>(r)] = std::current_exception();
      }

      std::lock_guard lk(s.m);
      s.live[static_cast<std::size_t>(r)] = 0;
      --s.nlive;
      /* Every rank still running is waiting in an exchange this one was
         supposed to join, and will wait forever.  Saying so and releasing them
         is the difference between a reported failure and a hung test. */
      if (s.nlive > 0 && s.nposted == s.nlive) {
        std::ostringstream os;
        os << "rank " << r << " returned while " << s.nposted
           << " other rank(s) were waiting in an exchange";
        s.complaints.push_back(os.str());
        s.finish_round();
      }
      s.hand_off();
    });
  }

  for (auto &t : threads)
    t.join();

  /* A body that threw is the more informative failure, so it goes first --
     a rank that aborts mid-exchange also produces a complaint about the round
     it abandoned, and that complaint is a consequence rather than a cause. */
  for (int r = 0; r < n; ++r)
    if (s.err[static_cast<std::size_t>(r)])
      std::rethrow_exception(s.err[static_cast<std::size_t>(r)]);

  if (!s.complaints.empty()) {
    std::ostringstream os;
    os << "Comm::Loopback: the ranks did not agree on what to exchange:";
    for (std::string const &c : s.complaints)
      os << "\n  " << c;
    throw std::runtime_error(os.str());
  }
}

} // namespace Comm
} // namespace Ume
