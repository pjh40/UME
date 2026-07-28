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
  \file Ume/Partition_METIS.cc

  `metis-kway-cut`, `metis-kway-vol`, `metis-rb`: the dual graph handed to
  METIS.

  `Zone_Graph` is already in the compressed-row form METIS wants, so there is
  no algorithm here at all -- what there is is a marshalling layer and a
  handful of degenerate cases METIS declines to take.  Three methods rather
  than one because the *objective* is the experiment: `metis-kway-cut`
  minimizes the weight of the cut faces, `metis-kway-vol` minimizes METIS's
  own communication-volume surrogate, and putting the two next to
  `Metrics::point_lambda1` -- the number of doubles that actually go on the
  wire -- is how we find out how good a proxy either of them is.

  The whole file compiles away when the build was configured without METIS,
  leaving only a factory that recognizes nothing.  `make()` calls it
  unconditionally and never reaches it for a `metis-*` name in such a build,
  because the table entry says the method is unavailable and `make()` stops
  there.
*/

#include "Ume/Partition.hh"
#include "Ume/Partition_Graph.hh"
#include "Ume/Partition_Impl.hh"
#include <memory>

#ifdef NO_METIS

namespace Ume {
namespace Part {

std::unique_ptr<Partitioner> make_metis(Method const &) { return nullptr; }

} // namespace Part
} // namespace Ume

#else

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <metis.h>

namespace Ume {
namespace Part {

namespace {

/* `idx_t` is 32 or 64 bits wide depending on how METIS was configured, and
   nothing in the header makes a mismatch against our `int` arrays a compile
   error -- it would be silent memory corruption instead.  So every array
   crosses the boundary through a `std::vector<idx_t>` and comes back the same
   way, which is correct at either width; this assertion documents the width
   this build found rather than gating on it. */
static_assert(sizeof(idx_t) == 4 || sizeof(idx_t) == 8,
    "METIS idx_t is neither 32 nor 64 bits wide");

//! Copy an `int` array into an `idx_t` one, raising anything below `floor`
std::vector<idx_t> marshal(std::vector<int> const &v, int floor) {
  std::vector<idx_t> out(v.size());
  std::ranges::transform(v, out.begin(),
      [floor](int x) { return static_cast<idx_t>(std::max(x, floor)); });
  return out;
}

//! Which METIS entry point a method calls, and what it asks it to minimize
struct Metis_Kind {
  //! `METIS_PartGraphKway` if true, `METIS_PartGraphRecursive` if false
  bool kway;
  //! `METIS_OBJTYPE_CUT` or `METIS_OBJTYPE_VOL`; meaningless unless `kway`
  int objtype;
};

//! Hand the dual graph to METIS and take its answer
class Metis : public Partitioner {
public:
  Metis(Metis_Kind kind, char const *name, char const *describe)
      : kind_{kind}, name_{name}, describe_{describe} {}

  char const *name() const override { return name_; }
  char const *describe() const override { return describe_; }

  void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
      Zone_Map &z2r) const override {
    begin_assign(mesh, opt, z2r);
    Zone_Graph const g = build_zone_graph(mesh, opt);
    if (g.nverts == 0)
      return;

    auto place = [&g, &z2r](int v, int part) {
      z2r[static_cast<std::size_t>(g.v2z[static_cast<std::size_t>(v)])] = part;
    };

    /* METIS declines both degenerate part counts -- one part, and more parts
       than there are vertices to fill them -- and each has exactly one
       sensible answer anyway: everything in part 0, and one vertex per part
       from part 0 up, which is what `split_sorted` and `bfs` do with the same
       shortage. */
    if (opt.nparts == 1 || g.nverts <= opt.nparts) {
      for (int v = 0; v < g.nverts; ++v)
        place(v, std::min(v, opt.nparts - 1));
      return;
    }

    /* METIS divides by the total weight of each constraint, so a zero weight
       anywhere is a division by zero rather than a free vertex: hence the
       floor of 1 on both. `Edge_Weight::FACE_AREA` really can round a sliver
       face down to zero, and `Vertex_Weight::SIDES` gives 0 to a zone whose
       sides are all phantom. */
    std::vector<idx_t> xadj = marshal(g.xadj, 0);
    std::vector<idx_t> adjncy = marshal(g.adjncy, 0);
    std::vector<idx_t> vwgt = marshal(g.vwgt, 1);
    std::vector<idx_t> adjwgt = marshal(g.adjwgt, 1);
    std::vector<idx_t> part(static_cast<std::size_t>(g.nverts), 0);

    idx_t nvtxs = static_cast<idx_t>(g.nverts);
    idx_t ncon = 1;
    idx_t nparts = static_cast<idx_t>(opt.nparts);
    idx_t objval = 0;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    /* METIS reads `ufactor` as an imbalance of `1 + ufactor/1000`, and takes
       it as an error below 1 -- so `--imbalance 1.0`, which is a perfectly
       reasonable thing to ask for, becomes the tightest balance METIS has
       rather than a rejected run. */
    options[METIS_OPTION_UFACTOR] = static_cast<idx_t>(
        std::max(1L, std::lround((opt.imbalance - 1.0) * 1000.0)));
    options[METIS_OPTION_SEED] = static_cast<idx_t>(opt.seed & 0x7FFFFFFFULL);
    if (kind_.kway)
      options[METIS_OPTION_OBJTYPE] = kind_.objtype;

    int const rc = kind_.kway
        ? METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(),
              vwgt.data(), nullptr, adjwgt.data(), &nparts, nullptr, nullptr,
              options, &objval, part.data())
        : METIS_PartGraphRecursive(&nvtxs, &ncon, xadj.data(), adjncy.data(),
              vwgt.data(), nullptr, adjwgt.data(), &nparts, nullptr, nullptr,
              options, &objval, part.data());
    if (rc != METIS_OK)
      throw std::runtime_error(std::string{name_} +
          ": METIS failed with status " + std::to_string(rc));

    for (int v = 0; v < g.nverts; ++v)
      place(v, static_cast<int>(part[static_cast<std::size_t>(v)]));
  }

private:
  Metis_Kind kind_;
  char const *name_;
  char const *describe_;
};

} // namespace

std::unique_ptr<Partitioner> make_metis(Method const &m) {
  std::string_view const name{m.name};
  /* Each method calls the entry point it is named after, at every part count.
     Karypis recommends recursive bisection below eight parts and k-way above,
     but silently switching would make `metis-kway-vol` ignore the objective
     it exists for: `METIS_OPTION_OBJTYPE` is a k-way option, and
     `METIS_PartGraphRecursive` minimizes the edge cut whatever it is set to.
     `metis-rb` is how you ask for the recursive one. */
  if (name == "metis-kway-cut")
    return std::make_unique<Metis>(
        Metis_Kind{true, METIS_OBJTYPE_CUT}, m.name, m.describe);
  if (name == "metis-kway-vol")
    return std::make_unique<Metis>(
        Metis_Kind{true, METIS_OBJTYPE_VOL}, m.name, m.describe);
  if (name == "metis-rb")
    return std::make_unique<Metis>(
        Metis_Kind{false, METIS_OBJTYPE_CUT}, m.name, m.describe);
  return nullptr;
}

} // namespace Part
} // namespace Ume

#endif
