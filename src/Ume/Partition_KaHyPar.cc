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
  \file Ume/Partition_KaHyPar.cc

  `kahypar`, `kahypar-quality`, `kahypar-det`: the hypergraph handed to
  Mt-KaHyPar, minimizing the connectivity metric.

  This is the method the whole exercise exists to measure against.  Every
  other partitioner here minimizes something *correlated* with UME's
  communication -- cut faces, a volume surrogate, surface area -- whereas
  \f$\sum(\lambda-1)\f$ over `Zone_Hypergraph`'s nets is the count of doubles
  on the wire itself, so this method optimizes the number the report prints
  rather than a stand-in for it.  Whether that is worth its dependency is the
  experiment, and `test_partition.cc` runs it.

  Three methods rather than one because Mt-KaHyPar's presets differ in a way
  that matters to this project beyond quality: `DEFAULT` and `QUALITY` refine
  in parallel and race, so the same seed on the same mesh gives a *different*
  partition run to run, and only `DETERMINISTIC` keeps the promise
  `Options::seed` makes.  What determinism costs in quality is exactly the
  sort of thing worth being able to measure.

  Two cautions if you touch this file.  Mt-KaHyPar's C header exports
  unscoped ALL-CAPS enumerators -- `CUT`, `KM1`, `DEFAULT`, `QUALITY`,
  `EPSILON`, `SUCCESS` and, most awkwardly, `METIS` -- so this file must not
  include `metis.h` as well, and must not name anything of its own briefly.
  And `mt_kahypar_hypernode_id_t` is 64 bits wide against UME's `int`, which
  makes the marshalling discipline mandatory rather than defensive.

  The whole file compiles away when the build was configured without
  Mt-KaHyPar, leaving a factory that recognizes nothing; `Partition_METIS.cc`
  explains why that is the shape.
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include "Ume/Partition_Impl.hh"
#include <memory>

#ifdef NO_KAHYPAR

namespace Ume {
namespace Part {

std::unique_ptr<Partitioner> make_kahypar(Method const &) { return nullptr; }

} // namespace Part
} // namespace Ume

#else

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <mtkahypar.h>

namespace Ume {
namespace Part {

namespace {

//! Copy an `int` array into a wider library one, raising anything below `low`
/*! `low` is 1 for the weights, which a partitioner divides by, and 0 for the
    pins, where vertex 0 is an ordinary vertex and clamping it would quietly
    move every net that touches it onto vertex 1. */
template <class T> std::vector<T> marshal(std::vector<int> const &v, int low) {
  std::vector<T> out(v.size());
  std::ranges::transform(v, out.begin(),
      [low](int x) { return static_cast<T>(std::max(x, low)); });
  return out;
}

//! Bring up Mt-KaHyPar's thread pool, once per process
/*! `mt_kahypar_initialize` must be called before anything else in the library
    and must not be called twice, which makes it the one piece of global state
    this file has.  The pool is sized to the machine; a `DEFAULT` or `QUALITY`
    partition is not reproducible anyway, and `DETERMINISTIC` is specified to
    give the same answer whatever the pool size. */
void start_library() {
  static std::once_flag once;
  std::call_once(once, [] {
    unsigned const hw = std::thread::hardware_concurrency();
    mt_kahypar_initialize(hw > 0 ? hw : 1u, false);
  });
}

//! Throw with whatever the library said, and free what it allocated to say it
[[noreturn]] void fail(
    char const *name, char const *what, mt_kahypar_error_t &err) {
  std::string msg{name};
  msg += ": ";
  msg += what;
  if (err.msg != nullptr && err.msg_len > 0)
    msg += ": " + std::string(err.msg, err.msg_len);
  mt_kahypar_free_error_content(&err);
  throw std::runtime_error(msg);
}

/*! @name Owners for the three library handles

    Each is freed by its own entry point rather than by `delete`, and there is
    a throw between every acquisition and the next. */
///@{
struct Context_Handle {
  mt_kahypar_context_t *p;
  explicit Context_Handle(mt_kahypar_preset_type_t preset)
      : p{mt_kahypar_context_from_preset(preset)} {}
  ~Context_Handle() {
    if (p != nullptr)
      mt_kahypar_free_context(p);
  }
  Context_Handle(Context_Handle const &) = delete;
  Context_Handle &operator=(Context_Handle const &) = delete;
};

struct Hypergraph_Handle {
  mt_kahypar_hypergraph_t h;
  explicit Hypergraph_Handle(mt_kahypar_hypergraph_t hg) : h{hg} {}
  ~Hypergraph_Handle() { mt_kahypar_free_hypergraph(h); }
  Hypergraph_Handle(Hypergraph_Handle const &) = delete;
  Hypergraph_Handle &operator=(Hypergraph_Handle const &) = delete;
};

struct Partitioned_Handle {
  mt_kahypar_partitioned_hypergraph_t h;
  explicit Partitioned_Handle(mt_kahypar_partitioned_hypergraph_t phg)
      : h{phg} {}
  ~Partitioned_Handle() { mt_kahypar_free_partitioned_hypergraph(h); }
  Partitioned_Handle(Partitioned_Handle const &) = delete;
  Partitioned_Handle &operator=(Partitioned_Handle const &) = delete;
};
///@}

//! Hand the hypergraph to Mt-KaHyPar and take its answer
class KaHyPar : public Partitioner {
public:
  KaHyPar(
      mt_kahypar_preset_type_t preset, char const *name, char const *describe)
      : preset_{preset}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    std::vector<int> const real = begin_assign(mesh, opt, z2r);
    Zone_Hypergraph const h = build_zone_hypergraph(mesh, opt);
    if (h.nverts == 0)
      return;

    /* The degenerate part counts, answered the way every other method here
       answers them rather than by asking the library what it makes of them:
       one part, and more parts than there are zones to fill them. */
    if (opt.nparts == 1 || h.nverts <= opt.nparts) {
      for (int v = 0; v < h.nverts; ++v)
        z2r[static_cast<std::size_t>(h.v2z[static_cast<std::size_t>(v)])] =
            std::min(v, opt.nparts - 1);
      return;
    }

    /* A hypergraph with no nets has nothing to minimize -- no two real zones
       of this mesh share a point -- so every balanced answer is optimal and
       the cheapest one to reach is `linear`'s. */
    if (h.nnets == 0) {
      std::vector<double> key(z2r.size(), 0.0);
      for (std::size_t z = 0; z < key.size(); ++z)
        key[z] = static_cast<double>(z);
      split_sorted(real, key, zone_weights(mesh, opt), 0, opt.nparts, z2r);
      return;
    }

    start_library();
    /* Not thread safe, and global to the library rather than to the context,
       which is why nothing here may be run from two threads at once. */
    mt_kahypar_set_seed(static_cast<std::size_t>(opt.seed));

    Context_Handle ctx{preset_};
    if (ctx.p == nullptr)
      throw std::runtime_error(
          std::string{name_} + ": Mt-KaHyPar rejected its own preset");
    /* Mt-KaHyPar reads `epsilon` as the fraction a block may exceed the
       average by, which is `imbalance - 1`, and takes the balance as a hard
       constraint rather than as a hint. */
    mt_kahypar_set_partitioning_parameters(ctx.p,
        static_cast<mt_kahypar_partition_id_t>(opt.nparts),
        std::max(0.0, opt.imbalance - 1.0), KM1);

    auto const pins = marshal<mt_kahypar_hypernode_id_t>(h.net_pins, 0);
    auto const netwgt = marshal<mt_kahypar_hyperedge_weight_t>(h.netwgt, 1);
    auto const vwgt = marshal<mt_kahypar_hypernode_weight_t>(h.vwgt, 1);

    /* Each handle is taken only once the call that produced it has been
       checked: a failed call has nothing to free, and freeing the null
       handle it hands back would be a second guess at the library. */
    mt_kahypar_error_t err{};
    /* The hypergraph is built through the context because the library keeps
       different data structures for different presets, so the two have to
       agree; `net_index` is the one array already the width it wants. */
    mt_kahypar_hypergraph_t const raw_hg = mt_kahypar_create_hypergraph(ctx.p,
        static_cast<mt_kahypar_hypernode_id_t>(h.nverts),
        static_cast<mt_kahypar_hyperedge_id_t>(h.nnets), h.net_index.data(),
        pins.data(), netwgt.data(), vwgt.data(), &err);
    if (err.status != SUCCESS)
      fail(name_, "could not build the hypergraph", err);
    Hypergraph_Handle hg{raw_hg};

    mt_kahypar_partitioned_hypergraph_t const raw_part =
        mt_kahypar_partition(hg.h, ctx.p, &err);
    if (err.status != SUCCESS)
      fail(name_, "could not partition the hypergraph", err);
    Partitioned_Handle part{raw_part};

    std::vector<mt_kahypar_partition_id_t> block(
        static_cast<std::size_t>(h.nverts), 0);
    mt_kahypar_get_partition(part.h, block.data());
    for (int v = 0; v < h.nverts; ++v)
      z2r[static_cast<std::size_t>(h.v2z[static_cast<std::size_t>(v)])] =
          static_cast<int>(block[static_cast<std::size_t>(v)]);
  }

private:
  mt_kahypar_preset_type_t preset_;
  char const *name_;
  char const *describe_;
};

} // namespace

std::unique_ptr<Partitioner> make_kahypar(Method const &m) {
  std::string_view const name{m.name};
  if (name == "kahypar")
    return std::make_unique<KaHyPar>(DEFAULT, m.name, m.describe);
  if (name == "kahypar-quality")
    return std::make_unique<KaHyPar>(QUALITY, m.name, m.describe);
  if (name == "kahypar-det")
    return std::make_unique<KaHyPar>(DETERMINISTIC, m.name, m.describe);
  return nullptr;
}

} // namespace Part
} // namespace Ume

#endif
