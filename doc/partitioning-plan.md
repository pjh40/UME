# Partitioning a serial `SOA_Idx::Mesh` across MPI ranks

## Status

**Phase 0 complete.** Four pre-existing defects fixed, one commit each, all
three build trees green at 63/63:

| Commit | Fix |
|---|---|
| `18f02df` | R0 — `MPISupport.cmake` appends to `COMMON_LINK_LIBRARIES` instead of clobbering it |
| `3b98a77` | R1 — `Mesh` copy/move deleted; it is pinned to its construction address. `ume_serial` holds ranks by `unique_ptr` |
| `02685d1` | R4 — `faces` compared in `Mesh::operator==` |
| `cbc5beb` | R3 — the facet audit counts ghost sides |

Two notes for whoever picks this up:

- R1 landed *stronger* than this plan specifies. Rather than sorting an index
  permutation, `Mesh`'s copy and move operations are deleted outright, so the
  bug class is a compile error. `ume_serial.cc` was the only place a `Mesh` had
  value semantics.
- R3's scope was verified, not assumed. Ghosting a single side produces three
  extra audit failures, but that is not what a decomposition does; ghosting a
  whole zone together with its sides and corners leaves `(face, edge) side
  counts` as the only firing check, as claimed below.

**Phase 1, commit 5 complete.** `Partition.hh/.cc` (interface, method table,
`linear`, `slab`, `slab-x/y/z`, `real_zones`, `zone_weights`, `split_sorted`),
`Partition_Graph.hh/.cc` (centroids only), `Ume::attach_serial_transport`
promoted out of `Brick` into `Mesh_Base.hh`/`Comm_Transport.cc`, and
`test/test_partition.cc`. All three trees green at 68/68.

Three notes for whoever picks this up:

- `Partitioner::assign` takes `Mesh const &` and the centroid path still works,
  because `Datastore::caccess_*` is const but calls `init_()` — the lazy
  derivation runs through the `mutable` state. So no partitioning code needs a
  mutable mesh.
- `split_sorted`'s `key`/`wgt` spans are indexed by *zone index*, not by
  position in `v`, so `rcb`/`rib` can split a subset without recompacting.
- Several zones of one Kuhn-subdivided cell share a centroid coordinate
  exactly, so a tie group straddles a slab cut. Contiguity is an invariant of
  the `(key, zone index)` sort order, not of the bare coordinate — the interval
  test that looked obvious is wrong, and `test_partition.cc` says why.

**Phase 1, commit 6 complete.** `Decompose.hh/.cc` — `Ghost_Depth`,
`Owner_Rule`, `Decompose_Options`, `Entity_Layout`, `build_decomposition`, and
nine layout tests in `test_partition.cc`. All three trees green at 77/77.

Three notes:

- The ghost side/corner rule generalized cleanly. Rather than "the ring on each
  cut face", the implementation imports every side of a ghost zone whose *face
  this rank already holds*, and every corner of a ghost zone whose *point this
  rank already holds*. For `FACE_ADJACENT` that is exactly the plan's rule (a
  cut face is precisely a face of the ghost zone that this rank holds), it
  extends to `POINT_ADJACENT` without a second rule, and it makes "no ghost
  points, edges or faces" true by construction rather than by argument.
- `lonely_facets()` in `test_partition.cc` checks the `(face, edge)` invariant
  straight off the layouts, before any mesh exists. `Ghost_Depth::NONE` is
  tested to *break* it, so the ghost machinery is demonstrably load-bearing
  rather than merely present.
- A shared point that a rank holds but does not source is a COPY in the real
  index range, **not** a ghost. `owner[g] != r` therefore does not mean "ghost
  on r" — the two lists have to be compared directly. Getting this wrong is
  easy and the layouts will not complain.

**Phase 1, commit 7 complete.** `extract()` in `Decompose.cc`, plus three
pre-existing defects it uncovered, one commit each. All three trees green at
86/86.

| Commit | Fix |
|---|---|
| `4fac608` | R8 — `Comm::pack` discarded `std::copy`'s return, so every element of a `VEC3V` field packed over the first one. Every Vec3 gather/scatter was wrong; invisible until a mesh is actually decomposed |
| `4df2fd2` | R9 — the derived variables built Kokkos Views from `&array[0]`, which aborts on an empty array. An empty partition is a legitimate answer |
| `f84dd7a` | `split_sorted` put a lone zone in part k/2, stranding part 0, when there were fewer zones than parts |

R8 is the one worth remembering: `point_norm` came out nonzero on interior
points of a partition, which reads exactly like a partitioning bug and is not
one. What made it findable was that clearing `myCpys`/`mySrcs` before touching
`point_norm` made it go away — the accumulation loop was right and the
gathscat was destroying the answer.

Notes on `extract` itself:

- Every partition of the sweep passes the *full* `mesh_audit` unchanged, and
  N=1 is bit-identical to the input. The plan's claim that R3 is the only audit
  change needed held up.
- The neighbor lists get their positional correspondence for free: both sides
  build in ascending serial-global order, so `mySrcs[P]` on R and `myCpys[R]`
  on P are the same set enumerated the same way, with nothing exchanged.
- Only `m:s>s5` on a ghost side and `m:s>s2` under `Ghost_Depth::NONE` are
  allowed to dangle; both become self-references. Every other map is required
  to resolve, and `Renumber::require` throws if it does not.

**Phase 1, commit 8 complete.** `src/partition.cc` → `ume_partition`, with
`--ranks --method --list-methods --vertex-weight --edge-weight --imbalance
--seed --ghosts --owner --zone-map --dump-zone-map --audit --no-write`, plus
end-to-end CTests. All three trees green at 92/92.

`gen_brick → ume_partition --audit → ume_serial --audit` on the written files
is the first thing to prove the format works from both ends. The CTests wire
exactly that up with CTest fixtures.

`--metrics`/`--comm-matrix`/`--check` are deliberately absent; they land with
the code behind them in phase 2.

**Next: Phase 2, commit 10** — `partition_check.hh/.cc` and
`ume_partition --check`. Phase 1's commit 9 (the audit sweep over
`test_brick.cc` specs) landed early, inside commit 7: `test_partition.cc`
already sweeps three specs × four methods × N ∈ {1,2,3,4,7,8} through the full
`mesh_audit`, plus a 1×1×1 brick split eight ways.

Build trees: `build/` (`UME_SERIAL=OFF`), `build-serial/` (`UME_SERIAL=ON`),
`build-mpi/` (`USE_MPI=ON`, created in commit 1). All three are configured with
`-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16` against the spack `ume`
environment.

---

## Context

`gen_brick` produces a single-rank mesh: `mype=0`, `numpe=1`, every entity
`INTERNAL`, empty neighbor lists (`Brick_Mesh.cc:440-447`, and `doc/brick_mesh.md`
§9 says so outright). The only way to get a multi-rank `.ume` set today is a DELFI
`UmeDump` through `txt2bin` — so the parallel gather/scatter paths, which are the
whole point of UME as a proxy app, cannot be exercised on a synthetic mesh at all.

We want to fix that *and* make the decomposition itself a variable: slabs, blocks,
space-filling curves, round-robin, and a true minimum-communication-volume
partition, so we can measure how the communication pattern affects performance.

The deliverable is a `ume_partition` executable that reads any serial `.ume` and
writes N per-rank `<base>.%05d.ume` files, with the partitioning *strategy* behind
a one-function interface so methods are interchangeable.

---

## The three facts that determine the design

**(a) `txt2bin.cc` fixes what the parallel arrays mean.** Verified in
`read_points` (`txt2bin.cc:253-265`) and `Entity::resize` (`SOA_Entity.cc:102-110`):
`cpy_idx`, `src_pe`, `src_idx`, `ghost_mask` all have length *ghost count* and are
indexed by ghost offset, not entity index. Shared points/edges/faces are **not**
ghosts — they live in `[0, lsize)` with a normal mask and `comm_type ∈ {SOURCE,
COPY}`, recorded only in `myCpys`/`mySrcs`. Two distinct mechanisms; a partition
uses both. `read_mpi` (`txt2bin.cc:503`) confirms `myCpys` = "Recv From",
`mySrcs` = "Send To".

**(b) The kernels don't need ghost zones for correctness — the shared-point
gathscat does the work.** `gradzatz` accumulates over corners with `cmask >= 1`
then does `points.gathscat(SUM, ...)` (`gradient.cc:104-123`). If ranks A and B
each own some zones around point *p*, each contributes its own corners and the SUM
reproduces the serial value exactly.

**(c) But two things force ghost *sides* to exist anyway.** `calc_face_area`
(`face_area.cc:87`) writes `d_side_tag(s2)` unconditionally, so `m:s>s2` must be
in range for every real side including one on a cut face. And `Auditor::ranges`
(`mesh_audit.cc:268-274, 322-331`) requires `m:c>p` in range and `m:c>z >= 0` for
*every* corner including ghosts.

⇒ **Zones/corners/sides are strictly partitioned and get one face-adjacent ghost
layer; points/edges/faces on the cut are duplicated in the real index range as
SOURCE/COPY.** That is exactly the DELFI model the format was built for.

---

## Architecture

```
src/Ume/Partition.hh / .cc          strategy interface + registry
src/Ume/Partition_Graph.hh / .cc    shared inputs: dual graph, hypergraph, centroids
src/Ume/Partition_Geometric.cc      slab/pencil/block/rcb/rib/morton/hilbert/rr/random
src/Ume/Partition_BFS.cc            zero-dependency graph-growing partitioner
src/Ume/Partition_METIS.cc          #ifdef'd out when NO_METIS
src/Ume/Partition_KaHyPar.cc        #ifdef'd out when NO_KAHYPAR
src/Ume/Decompose.hh / .cc          Decomposition + extract   <-- bulk of the work
src/Ume/Partition_Metrics.hh / .cc  reporting
src/Ume/partition_check.hh / .cc    cross-rank verification
src/partition.cc                    the ume_partition driver
```

### Strategy interface (`Partition.hh`)

```cpp
namespace Ume { namespace Part {

using Zone_Map = std::vector<int>;   //!< zone index -> rank, length zones.size()

enum class Vertex_Weight { UNIT, CORNERS, SIDES };
enum class Edge_Weight   { UNIT, FACE_POINTS, FACE_AREA };

struct Options {
  int nparts{1};
  std::uint64_t seed{0x9E3779B97F4A7C15ULL};
  Vertex_Weight vwgt{Vertex_Weight::CORNERS};
  Edge_Weight   ewgt{Edge_Weight::FACE_POINTS};
  double imbalance{1.03};
  std::map<std::string, std::string> extra;
};

class Partitioner {
public:
  virtual ~Partitioner() = default;
  virtual char const *name() const = 0;
  virtual char const *describe() const = 0;
  //! Assign every zone with mask >= 1 to a part in [0, opt.nparts)
  virtual void assign(SOA_Idx::Mesh const &mesh, Options const &opt,
                      Zone_Map &z2r) const = 0;
};

std::vector<std::string> methods();                  //!< what this build supports
std::unique_ptr<Partitioner> make(std::string_view); //!< nullptr if unknown
}}
```

`methods()`/`make()` are a plain table in `Partition.cc` bracketed by
`#ifndef NO_METIS` / `#ifndef NO_KAHYPAR`. No self-registration, no static-init
ordering games — the set is known at compile time.

**Method availability is a first-class concept.** The four builds (neither
library, METIS only, Mt-KaHyPar only, both) each expose a different method set,
and every consumer has to cope:

- `methods()` returns only what this build supports; `make(name)` returns
  `nullptr` otherwise.
- `ume_partition --list-methods` prints the available set, one per line with
  `describe()`, and exits 0. Keep a compile-time-complete list of *all* known
  names so an unavailable method is distinguishable from a misspelled one:
  `--method metis-kway-cut` on a no-METIS build must say "requires METIS, which
  this build was configured without", not "unknown method".
- The Catch2 sweep (`test_partition.cc`) iterates `Part::methods()` rather than a
  hardcoded list, so one test binary covers whatever is linked.
- `report_csv` writes the method name it actually ran, so sweeps from differently
  configured builds concatenate into one table cleanly.

**`zone_to_rank` is the entire contract.** Every library reduces to it (METIS
`part[]`, KaHyPar `get_partition`), it serializes trivially (`--dump-zone-map` /
`--zone-map` make the strategy and extraction layers independently testable), and
everything downstream — ownership, ghosts, renumbering, comm patterns — is a pure
function of it. Boundary zones (`mask == -1`) are assigned by the framework, not
by the `Partitioner`.

### Shared inputs (`Partition_Graph.hh`)

```cpp
struct Zone_Graph {                 //!< CSR dual graph over real zones
  int nverts{0};
  std::vector<int> xadj, adjncy, vwgt, adjwgt;
  std::vector<int> z2v, v2z;
};
Zone_Graph build_zone_graph(SOA_Idx::Mesh const &, Options const &);

struct Zone_Hypergraph {            //!< vertices = zones, nets = points
  int nverts{0}, nnets{0};
  std::vector<std::size_t> net_index;
  std::vector<int> net_pins, vwgt, netwgt;
  std::vector<int> z2v, v2z;
};
Zone_Hypergraph build_zone_hypergraph(SOA_Idx::Mesh const &, Options const &);

DS_Types::VEC3V_T const &zone_centroids(SOA_Idx::Mesh &);  //!< = zcoord
```

- Dual graph from `m:f>z1`/`m:f>z2`, edge when both zones are real and distinct.
  `FACE_POINTS` edge weight = distinct points on the face = points that may have to
  be duplicated across that cut.
- `CORNERS` vertex weight = `m:z>c` row length, the actual trip count of every
  `gradzatz` kernel. On a hex/tet mix, balancing *zones* is misleading (hex 8
  corners, tet 4).
- Hypergraph is literally `m:p>zs` filtered to real zones — reuse the existing
  derived `Points::VAR_point_to_zones`, don't rebuild it.
- **Geometric methods work off `zcoord` centroids, not lattice indices.** Sorting
  by centroid reproduces the lattice cuts exactly on a brick and generalizes to any
  `.ume`. `zcoord` is lazily derived and its `init_()` calls `zones().scatter()`,
  so a `Dummy_Transport` must be attached first — reuse `Brick::attach_serial_transport`
  (`Brick_Mesh.cc:703`), promoting it to `Ume::attach_serial_transport`.

---

## Built-in methods (no dependencies)

One helper does most of the work:

```cpp
//! Split compact zone ids into k weighted-equal groups by key, writing
//! parts [first_part, first_part+k) into z2r
void split_sorted(std::span<int const> v, std::span<double const> key,
                  std::span<int const> wgt, int first_part, int k, Zone_Map &);
```

Cutting at cumulative-*weight* boundaries rather than uniform spatial bins gives
exact load balance on a non-uniform mesh.

| method | notes |
|---|---|
| `slab-x/-y/-z/slab` | one `split_sorted`; bare `slab` picks the longest extent |
| `pencil` | factor N over the two longest axes, nested split |
| `block` | own `dims_create` (prime-factorize N, assign each factor to the axis with smallest `dims[i]/extent[i]`) — works without MPI *and* gives cube-ish subdomains in physical rather than index space |
| `rcb` | recursive coordinate bisection along each subset's own longest extent; handles prime N |
| `rib` | recursive inertial bisection on the principal eigenvector of the weighted centroid covariance; beats `rcb` on rotated/diagonal meshes |
| `morton` | 21-bit-per-axis interleave into `uint64_t`, sort, split |
| `hilbert` | Skilling's transpose algorithm. Worth having next to Morton — Morton's long jumps at power-of-two boundaries cost measurable surface-to-volume, and showing that gap is the kind of result we're after |
| `linear` | contiguous ranges of original zone index; the natural baseline for a `renumber_mesh`-ordered mesh |
| `roundrobin` | `z % N`. Deliberate worst case: nearly every face cut, comm graph complete |
| `random` | deterministic shuffle from `Options::seed` |
| `bfs` | greedy growing on the dual graph from farthest-point-sampled seeds, weight-targeted round-robin expansion. Gives a *non-geometric, topology-respecting* partition with zero dependencies, and is the natural seed for a later Kernighan–Lin pass |

These cover everything Zoltan's RCB/RIB/HSFC/Block would give us, without
Zoltan's `MPI_Init` requirement.

## Library methods

Both libraries get a `cmake/pkgs/<X>Support.cmake` following
`libunwindSupport.cmake` verbatim: `SimplePkg(<PKG> <header> <lib>)`, then define
`NO_<PKG>` on failure rather than erroring, else append to
`COMMON_LINK_LIBRARIES` and `<PKG>_INCLUDE_DIRS`. `include()` each in the
top-level `CMakeLists.txt` before `add_subdirectory(src)`;
`src/Ume/CMakeLists.txt` picks up the include dirs alongside
`${LIBUNWIND_INCLUDE_DIRS}`. Both are genuinely optional and independent — all
four combinations must build and pass the sweep.

> The `MPISupport.cmake` `set`-instead-of-`list(APPEND)` bug (R0 below) must be
> fixed before either module lands, or the second one to be `include()`d will
> silently drop the first.

### METIS — `metis-kway-cut`, `metis-kway-vol`, `metis-rb`

`SimplePkg(METIS metis.h metis)`. Input is `Zone_Graph` directly.
`METIS_PartGraphKway` for N ≥ 8, `METIS_PartGraphRecursive` below (Karypis's own
guidance — recursive bisection is better for small k). `METIS_SetDefaultOptions`,
then `options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT` or `METIS_OBJTYPE_VOL`,
`options[METIS_OPTION_UFACTOR] = lround((imbalance - 1.0) * 1000)`,
`options[METIS_OPTION_SEED] = seed`, `options[METIS_OPTION_NUMBERING] = 0`.

Exposing CUT and VOL as two separate methods is the point: it gives the
edgecut-vs-volume-surrogate-vs-λ−1 comparison directly, on the same graph, in one
sweep.

**The footgun:** `idx_t` is 32- or 64-bit depending on how METIS was configured
(`IDXTYPEWIDTH` in `metis.h`), and `real_t` likewise. Never pass
`std::vector<int>::data()` — marshal into `std::vector<idx_t>` at the boundary and
copy the result back. A 64-bit-idx METIS against our `int` arrays is silent memory
corruption, not a compile error. Add a `static_assert` on `sizeof(idx_t)` only to
document, not to gate; the marshalling makes either width correct.

METIS 5.1+ (KarypisLab) is Apache-2.0 — fine to depend on optionally. **ParMETIS
is not freely redistributable for commercial use; do not add it.**

### Mt-KaHyPar — `kahypar`, `kahypar-quality`

`SimplePkg(KAHYPAR libmtkahypar.h mtkahypar)`. C API: `mt_kahypar_initialize`,
`mt_kahypar_context_from_preset(DEFAULT|QUALITY)`,
`mt_kahypar_set_partitioning_parameters(ctx, k, epsilon, KM1)`,
`mt_kahypar_create_hypergraph(...)`, `mt_kahypar_partition`,
`mt_kahypar_get_partition`. Net weight = payload bytes per point; vertex weight =
corners.

**Why this is the one to beat.** The λ−1 (connectivity-minus-one) objective on
{vertices = zones, nets = points} is *not a surrogate* for UME's communication —
it is exactly UME's communication. A point touched by λ(p) distinct ranks needs
λ(p)−1 copies, and each copy is one element in one direction of
`points.gathscat`, the dominant traffic in `gradzatz`. Minimizing Σ(λ(p)−1)
minimizes the doubles on the wire, full stop. Graph edge-cut minimizes faces
instead; `METIS_OBJTYPE_VOL` minimizes a vertex-based volume surrogate. Running
all three against the same mesh and putting `point_lambda1` next to a stopwatch is
the headline experiment.

MIT-licensed (original KaHyPar is GPL-3 — do not add that path; we are BSD-3).

---

## The mesh-splitting machinery (`Decompose.hh` / `.cc`)

### Ownership and ghost rules

- **Zones** with `mask >= 1` go where `z2r` says. `mask == -1` boundary zones
  follow the real zone across their single face (resolve via the phantom side's
  `m:s>s2` partner's `m:s>z`, which works for both `Boundary_Style` values).
- **Corners, sides** follow their zone (`m:c>z`, `m:s>z`).
- **Points, edges, faces** are present on rank R if any zone of R touches them.
  **SOURCE rank = lowest-numbered rank containing them** — deterministic,
  symmetric, needs no negotiation. It does concentrate SOURCE duty on low ranks,
  which is itself measurable: expose `Owner_Rule ∈ {LOWEST_RANK, LEAST_LOADED, HASH}`,
  default `LOWEST_RANK`.
- **Ghost zones: one face-adjacent layer.** For each cut face (both `m:f>z1/z2`
  real, different ranks), the remote zone becomes a ghost. `Ghost_Depth ∈ {NONE,
  FACE_ADJACENT, POINT_ADJACENT}`, default `FACE_ADJACENT`. `NONE` forces cut-face
  `m:s>s2` to self and trips the audit — offer it, document as diagnostic-only.
- **Ghost corners:** one per *local* point each ghost zone touches. Never import a
  point only a remote zone touches — the audit demands every corner's `m:c>p` be
  in range.
- **Ghost sides:** the remote zone's whole side ring on each cut face, imported
  with their global identity rather than synthesized.
- **No ghost points/edges/faces.** Everything a real side names is real.

### Local index layout — state the rule once

For every entity kind, local index order is **`(is_ghost, global index)`
lexicographic**: locals first in ascending global order, then ghosts in ascending
global order.

Two consequences: `mask == 0 ⟺ index >= local_size()` holds automatically; and at
N=1 the mapping is the identity, so the partition is bit-identical to the input —
a free and very strong regression test.

It must be *lexicographic*, not plain ascending global: ghost zones are ordinary
`mask == 1` zones elsewhere, so their global indices interleave with the local
ones. Getting this wrong puts a ghost below `local_size()`, which makes
`VAR_zcoord` write past its `local_size()`-sized accumulator.

### The `Decomposition` object

```cpp
struct Entity_Layout {
  std::vector<int> start, gid;    //!< locals: gid[start[r]..start[r+1]) ascending
  std::vector<int> gstart, ggid;  //!< ghosts
  std::vector<int> owner;         //!< global index -> SOURCE rank, -1 if absent
  int local_size(int r) const;
  int size(int r) const;
  int local_index(int r, int g) const;   //!< binary search; -1 if absent
};

struct Decompose_Options {
  Ghost_Depth ghosts{Ghost_Depth::FACE_ADJACENT};
  Owner_Rule owner_rule{Owner_Rule::LOWEST_RANK};
  bool ghost_entity_comm{true};
  bool split_subsets{true};
};

struct Decomposition {
  int nparts{0};
  Part::Zone_Map z2r;
  Entity_Layout points, edges, faces, zones, corners, sides, iotas;
  Decompose_Options opt;
};

Decomposition build_decomposition(SOA_Idx::Mesh const &src, Part::Zone_Map z2r,
                                  int nparts, Decompose_Options const &);
void extract(SOA_Idx::Mesh const &src, Decomposition const &d, int pe,
             SOA_Idx::Mesh &dst);
```

`Decomposition` holds *all* ranks' layouts while `extract` produces *one* mesh.
That is the key structural decision: building rank R's `mySrcs` needs the local
index on rank P of a shared entity, and `Entity_Layout::local_index(P, g)` gives
it without P's mesh existing. So we never hold N `Mesh` objects at once (see R1).

### `extract()` step by step

1. **Header.** `ivtag = UME_VERSION_2`, `version_header = true`, `mype = pe`,
   `numpe = d.nparts`, `geo`/`dump_iotas` copied.
2. **Size, then build `g2l`.** `dst.X.resize(l, l+g, g)` first — `resize` is
   virtual and also sizes the datastore arrays (`Points::resize` does `pcoord`,
   `Sides::resize` all eleven `m:s>*`). Then a dense `std::vector<int> g2l(src.X.size(), -1)`.
3. **Masks and comm types.** `mask` copied verbatim for real entities, `0` for
   ghosts. `comm_type = GHOST` for ghosts; for points/edges/faces `INTERNAL` if on
   exactly one rank, else `SOURCE` if `owner == pe` else `COPY`; for real
   zones/corners/sides `SOURCE` if a ghost of it exists elsewhere, else `INTERNAL`.

   > **A cut is not a mesh boundary.** Do not set `pmask = -1` on cut points and do
   > not emit `mask == -1` phantom sides on cut faces. `VAR_point_norm` and the
   > `point_type == -1` branch of `gradzatp` would subtract a spurious outward-normal
   > component and the gradient would be quietly wrong near every partition boundary.
   > This is the most plausible-looking wrong answer available.

4. **Renumber the stored maps** — straight `g2l`, plus:
   - `m:f>z1/z2`: `g2l` of the global zone (real *or* ghost), else `-1`. `-1` is
     legal here per the audit, but with `FACE_ADJACENT` it always resolves. Preserve
     which was `z1`.
   - `m:s>s3/s4` (ring within one face of one zone) and `m:s>s5` (same zone) always
     resolve.
   - **`m:s>s2` across a cut** resolves to the imported ghost side. Because ghost
     sides are the genuine remote sides from the global mesh, `s2` stays an
     involution and `m:s>f`/`m:s>e` agree between the pair — every audit `s2` check
     holds.
   - Ghost sides' own maps: `m:s>z` = ghost zone, `m:s>f/e/p1/p2` = local reals,
     `m:s>c1/c2` = ghost corners, `m:s>s2` = the real partner, `m:s>s3/s4` = the
     imported ring. `m:s>s5` points into faces we didn't import — set it to self.
     Safe: every `Auditor::connectivity` check gates on `usable(s)`, which requires
     `smask != 0`.
5. **Ghost records.** For ghost offset `g` at local index `l = lsize + g`:
   `cpy_idx[g] = l`, `src_pe[g] = owner`, `src_idx[g] = layout.local_index(owner, gid)`,
   `ghost_mask[g] = src.X.mask[gid]` (`renumbering.cc` reads `ghost_mask == 0` as
   "null ghost", so propagating the owner's mask is the right semantics).
6. **`myCpys` / `mySrcs`.** Per entity kind, per remote rank P, **at most one
   `Neighbor` per pe** (two `Isend`s to the same rank with the same tag would pair
   with the wrong `Irecv` — `Comm_MPI.cc:69-102` matches by `(pe, tag)`, not by
   position in `remotes`). `mySrcs[P]` = local indices that are SOURCE here with a
   copy/ghost on P; `myCpys[P]` = local indices that are COPY/GHOST here sourced on
   P. **Both sorted by global index in the serial mesh** — element order *within* a
   neighbor is load-bearing, and sorting both sides by the same global key makes the
   correspondence positional by construction. Sort the neighbor list by pe too, for
   deterministic `write`/`operator==`.
7. **Subsets.** Keep elements whose zone is local, renumber, filter `mask` in
   lockstep, recompute `lsize`. `gen_brick` emits none; `txt2bin` meshes carry
   materials.

### Where the existing code will fight us

Each of R0–R4 is a pre-existing defect independent of this feature, and each gets
its own commit (see Phasing).

- **R0 — `MPISupport.cmake` clobbers `COMMON_LINK_LIBRARIES`.** It does
  `set(COMMON_LINK_LIBRARIES MPI::MPI_CXX)` where every other module appends,
  dropping the `Kokkos::kokkos` entry `KokkosSupport.cmake` just added. Currently
  harmless only because `src/Ume/CMakeLists.txt` re-links Kokkos explicitly. Adding
  two more optional modules makes it a real bug.
- **R1 — `Mesh` is move-unsafe.** `Entity` holds a raw `Mesh *mesh_` set at
  construction, so moving a `Mesh` leaves every `Entity::mesh_` dangling. **This is
  already a live bug at `src/ume_serial.cc:112-116`**, which `std::sort`s a
  `std::vector<Mesh>`. Fix before building on it (sort an index permutation).
- **R2 — `Entity::lsize_` is private with no setter**; only `resize(local, total,
  ghost)` sets it. Call it first, then fill.
- **R3 — the audit's `(face, edge) side counts` check fires on every cut face.**
  It forms groups from `usable(s)` sides only, so a cut face reports "has 1 sides".
  Fix in `mesh_audit.cc`: include mask-0 sides whose `(f, e)` are in range. This
  preserves the invariant's meaning and is the *only* audit change required — every
  other check in `ranges`/`connectivity`/`masks`/`parallel` passes under this layout.
- **R4 — `Mesh::operator==` omits `faces`** (`SOA_Idx_Mesh.cc:101-106` compares
  points, edges, sides, corners, zones, iotas). Verified. The N=1 identity test
  won't catch an `m:f>z1/z2` regression until this one-line fix lands.
- **R5 — Datastore derived entries are lazily initialized and never invalidated.**
  Fill every raw `m:*` array on `dst` before reading any derived key from `dst`.
- **R6 — `VAR_zcoord`/`VAR_fcoord` index `local_size()`-sized accumulators** by
  `m:c>z`/`m:s>f` with no bounds check. The ordering rule above is not stylistic.
- **R7 — never write derived maps.** `Mesh::write` only emits stored ones, so this
  is automatic — just don't get clever and copy `m:z>c` across.

---

## Metrics (`Partition_Metrics.hh`)

Computed from the `Decomposition` alone — no partition meshes needed, so a method
sweep with `--no-write` is cheap.

```cpp
struct Metrics {
  int nparts;
  std::vector<long> zones, real_zones, corners, sides, points, faces, edges;
  double zone_imbalance, corner_imbalance;
  long cut_faces, cut_faces_weighted, edge_cut, edge_cut_weighted;
  long point_lambda1, edge_lambda1, face_lambda1;   //!< sum over e of (lambda(e)-1)
  long ghost_zones, ghost_sides, ghost_corners;
  std::vector<long> cpy_elems, src_elems;
  std::vector<int>  neighbors;
  int max_neighbors; double mean_neighbors;
  double max_surface_to_volume, mean_surface_to_volume;
  std::vector<long> comm_matrix;                    //!< nparts x nparts
};
Metrics measure(SOA_Idx::Mesh const &, Decomposition const &);
void report(Metrics const &, std::ostream &);
void report_csv(Metrics const &, std::string_view method, std::ostream &);
```

**The headline number is `point_lambda1`** = Σ_p (λ(p) − 1). It is exactly the
number of point elements moved in one direction of `points.gathscat`, and
`gradzatz` does two `gathscat`s (a `DBLV` and a `VEC3V`) plus a `scatter`, so wire
bytes per `gradzatz` ≈ `point_lambda1 * (8 + 24) * 2 + point_lambda1 * 24` — a
number to put next to a stopwatch. Edge cut (`metis-kway-cut`) and METIS's volume
surrogate (`metis-kway-vol`) both approximate it; reporting all three next to
measured time *is* the experiment.

`max_neighbors` deserves equal billing — it is the latency-bound half of the
story, where `roundrobin`/`random` blow up to N−1 while `block` stays at 26.

`--metrics-csv` appends one row per `(method, nparts)` so a shell loop produces a
plottable table.

---

## Driver: `src/partition.cc` → `ume_partition`

```
ume_partition [options] <input.ume> <output-base>
  --ranks N                 number of partitions (required)
  --method NAME             see --list-methods
  --list-methods            print methods this build supports, exit
  --vertex-weight u|corners|sides
  --edge-weight   u|face-points|face-area
  --imbalance F             [1.03]
  --seed N
  --ghosts none|face|point  [face]
  --owner low|balanced|hash [low]
  --zone-map FILE           read zone->rank instead of running a method
  --dump-zone-map FILE      write the zone->rank a method produced
  --metrics / --metrics-csv FILE / --comm-matrix FILE
  --audit                   run mesh_audit on each partition before writing
  --check                   cross-rank consistency check
  --no-write                metrics only
```

`main()` mirrors `gen_brick.cc`: `Ume::initialize(argc, argv)` before anything
touches the Datastore (the centroid path forces `zcoord`, a Kokkos kernel), attach
`Dummy_Transport`, read, partition, then loop `for pe in [0,N)` { `Mesh dst;`
attach transport; `extract`; optional `audit`; write `<base>.%05d.ume`; optional
read-back `operator==` round trip } — one `Mesh` alive at a time except under
`--check`.

Then ~15 lines in `gen_brick.cc`: `--ranks N --method M` calls `Part::make`,
`build_decomposition`, and the same write loop. No duplicated logic.

---

## Verification

1. **N == 1 is the identity.** With the ordering rule, `extract(src, d, 0, dst)`
   gives `dst == src` except `numpe`. Catch2, near-free, catches renumbering bugs
   immediately. Fix R4 first or it won't check faces.
2. **`mesh_audit` on every partition**, reusing `Ume::SOA_Idx::audit` unchanged
   apart from R3. Sweep the existing `test_brick.cc` specs × N ∈ {1,2,3,4,7,8} ×
   {slab-x, slab-y, slab-z, block, roundrobin}. `roundrobin` is the important one —
   it maximizes cut faces and shared points and finds the layout bugs the nice
   methods hide. Also N > zone count (warn, don't crash) and a 1×1×1 brick.
3. **Cross-rank consistency — `partition_check.hh`.** Doesn't exist today and is
   worth the most:
   ```cpp
   bool check_partitions(std::vector<SOA_Idx::Mesh *> const &ranks, std::ostream &);
   ```
   Per entity kind: every `myCpys[i]` with `pe == P` has exactly one matching
   `mySrcs` entry on P **of the same length**, with `comm_type` COPY/GHOST here and
   SOURCE there; **each paired position is geometrically the same entity** (`pcoord`
   bit-equal for points, equal centroids through local maps otherwise) — this is the
   only thing that catches an off-by-one in the positional pairing, which is the
   most insidious bug available and the one MPI would turn into silently wrong
   physics; every ghost's `src_pe`/`src_idx` names a SOURCE entity that is the same
   entity geometrically; and **global reassembly** — summing non-COPY, non-GHOST
   entities over all ranks reproduces the serial counts exactly, and the real zones
   partition the serial real zones with no gaps or overlaps.

   Wire as `ume_partition --check`, and extend `ume_serial --audit` to audit *every*
   rank rather than only `ranks[0]` and then run `check_partitions` (fix R1 first).
4. **Numerical equivalence under MPI.** CTest guarded by `USE_MPI` +
   `MPIEXEC_EXECUTABLE`: `gen_brick`, `ume_partition --ranks 2`, `mpirun -n 2
   ume_mpi`, compare `zgrad` against serial. Needs a global-id side-car
   (`<base>.%05d.gid` under a flag) rather than a file-format change. Relative
   tolerance ~1e-12 — do **not** reuse `test_brick.cc`'s `GRAD_TOL == 0.0`, because
   partitioned accumulation order genuinely differs from serial.
5. **Threaded loopback transport (optional).** `Comm_Loopback.hh/.cc`: a
   `Loopback::Group(nranks)` handing out one `Transport*` per rank, `exchange()` as
   a mutex+condvar rendezvous, one thread per partition. ~80 lines. `Transport::exchange`
   is called synchronously from inside `Entity::gather`, so threads are unavoidable.
   Converts the most valuable test from "needs an MPI build and mpirun" to "runs in
   CI".

---

## Phasing and commit discipline

**Small, single-concern commits.** Every pre-existing-defect fix is its own commit
with no feature code mixed in, so each can be reviewed, bisected, or reverted on
its own. New code lands in commits that each leave the tree building and the test
suite green.

### Build trees

Two trees exist, `build/` (`UME_SERIAL=OFF`) and `build-serial/` (`UME_SERIAL=ON`),
both `Debug` and — importantly — **both `USE_MPI=OFF`**. So neither currently
exercises `MPISupport.cmake` at all, and the R0 fix would land unverified.

Add a third tree, `build-mpi/` (`USE_MPI=ON`), as the first step of commit 1. It is
needed three times over: to show R0 actually matters (configure must succeed and
`ume_mpi` must link Kokkos *and* MPI), to build `ume_mpi` for the phase-2 numerical
test, and to catch `#ifdef HAVE_MPI` bit-rot in `Comm_MPI.cc`, which nothing
compiles today.

Per-commit verification runs `ninja && ctest` in `build/` and `build-serial/`
throughout; `build-mpi/` from commit 1 onward. Phase 4 additionally needs the four
library combinations — configure each as a throwaway tree rather than reconfiguring
in place, since `SimplePkg` caches its result in `CHECK_<PKG>_DONE` and will not
re-probe.

### Phase 0 — prep. Four independent commits, no new features.

| # | Commit | Files |
|---|---|---|
| 1 | `Append rather than replace the common link libraries` (R0) | `cmake/pkgs/MPISupport.cmake` |
| 2 | `Stop sorting a vector of Mesh` (R1) | `src/ume_serial.cc` — sort an index permutation; `Entity::mesh_` back-pointers dangle on move |
| 3 | `Compare faces in Mesh::operator==` (R4) | `src/Ume/SOA_Idx_Mesh.cc` |
| 4 | `Count ghost sides in the facet audit` (R3) | `src/Ume/mesh_audit.cc`, `test/test_brick.cc` |

Commits 2–4 each need a test that fails before and passes after. Commit 4 is the
only one whose motivation is this feature, but it is a correctness fix to the audit
in its own right and reads as one.

### Phase 1 — minimum working 2-rank slab.

The extractor is the entire risk, so it lands with exactly one trivial method
rather than waiting on method variety.

| # | Commit | Contents |
|---|---|---|
| 5 | `Add the partitioner interface and centroid inputs` | `Partition.hh/.cc` (interface, registry, `linear`, `slab-x/y/z`), `Partition_Graph.hh/.cc` centroids only, `Ume::attach_serial_transport` promoted out of `Brick` |
| 6 | `Add the mesh decomposition layout` | `Decompose.hh/.cc` — `Entity_Layout`, `build_decomposition`; no `extract` yet, tested directly on layouts |
| 7 | `Extract a partition from a decomposition` | `extract()`; test: N=1 identity round trip |
| 8 | `Add the ume_partition driver` | `src/partition.cc`, `src/CMakeLists.txt`; `--ranks --method --list-methods --audit --no-write` |
| 9 | `Test partitioned bricks against the mesh audit` | `test/test_partition.cc` — sweep `test_brick.cc` specs × N ∈ {1,2,3,4,7,8} |

Commit 7 is the big one and cannot usefully be split further — the renumbering,
ghost records, and neighbor lists are mutually dependent and nothing passes until
all three are right.

### Phase 2 — proof. A decomposition cannot be trusted until this lands.

| # | Commit | Contents |
|---|---|---|
| 10 | `Add cross-rank partition consistency checks` | `partition_check.hh/.cc`, `ume_partition --check` |
| 11 | `Audit every rank in ume_serial` | `src/ume_serial.cc` — all ranks, then `check_partitions` |
| 12 | `Add partition metrics and reporting` | `Partition_Metrics.hh/.cc`, `--metrics --metrics-csv --comm-matrix` |
| 13 | `Compare partitioned and serial gradients under MPI` | CTest guarded by `USE_MPI`, global-id side-car, `ume_mpi` reference-compare option |

### Phase 3 — the zero-dependency methods.

One commit per coherent group, each dropping into the phase 1–2 harness unchanged:
`block`/`pencil`/`dims_create`; `rcb`/`rib`; `morton`/`hilbert`;
`roundrobin`/`random`; `bfs`. After this the full experiment runs on any machine
with no external libraries.

### Phase 4 — the two libraries, independently.

| # | Commit | Contents |
|---|---|---|
| — | `Add optional METIS support` | `cmake/pkgs/METISSupport.cmake`, `Partition_METIS.cc`, registry entries for `metis-kway-cut`/`metis-kway-vol`/`metis-rb` |
| — | `Add optional Mt-KaHyPar support` | `cmake/pkgs/KaHyParSupport.cmake`, `Partition_KaHyPar.cc`, hypergraph construction in `Partition_Graph.cc`, `kahypar`/`kahypar-quality` |

METIS first — it exercises the `Zone_Graph` path and the availability-gating
machinery, so the Mt-KaHyPar commit only has to add the hypergraph. Verify all four
configure combinations build and that `--list-methods` reports correctly in each
before moving on.

### Phase 5 — polish.

Threaded loopback transport + in-CI N-rank numerical test; `doc/partitioning.md`
(model, conventions, per-library licences, how to read the metrics); optional rank
reordering onto machine topology; Kernighan–Lin refinement on `bfs`. Separate
commits.
