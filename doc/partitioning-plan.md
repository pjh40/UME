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

Phase 1's commit 9 (the audit sweep over `test_brick.cc` specs) landed early,
inside commit 7: `test_partition.cc` already sweeps three specs × four methods
× N ∈ {1,2,3,4,7,8} through the full `mesh_audit`, plus a 1×1×1 brick split
eight ways.

**Phase 2, commit 10 complete.** `partition_check.hh/.cc` — seven cross-rank
checks, `ume_partition --check`, six Catch2 tests and one CTest. All three
trees green at 100/100.

Three notes:

- **An entity is identified by a *tuple* of point sets, not one.** A side's
  two endpoints do not name it: the two sides of a zone that meet along an
  edge, one on each of the two faces sharing it, have the same pair. Ghost
  corners have the same problem — two ghost zones meeting at a local point
  give it two. So the signature is `(own points, face points, zone points)`,
  compared component by component. The mixed brick really does contain such a
  side pair, and `two sides on one edge are told apart` swaps exactly one, so
  the weaker key is a failing test rather than a latent hole.
- Derived geometry (`fcoord`, `zcoord`) is exactly the wrong thing to compare
  with. It is filled on ghosts by a scatter that has not happened in a
  `Dummy_Transport` process, so it would compare garbage against a real value.
  Stored maps and `pcoord` are bit-identical across ranks because `extract`
  copies them, so every comparison here is exact and no tolerance appears.
- A ghost is compared by *containment*, not equality: a rank imports only the
  part of a ghost zone touching points it already holds. Ghost corners and
  ghost sides still match exactly on their own points; it is the zone
  component that is a proper subset.

**Phase 2, commit 11 complete.** `ume_serial --audit` gives every rank a
Transport, audits every rank rather than only `ranks[0]`, and then runs
`check_partitions` over the whole set. All three trees green at 101/101.

One note: a *partial* set is skipped, not failed. `ume_serial --audit
part.00001.ume` is a legitimate thing to type, and `numpe` says 3, so
`check_partitions` would have reported "1 of 3 ranks" as a disagreement. The
driver compares `ranks.size()` to `numpe` first and says it is skipping.
`partition_audit_all` is now the strongest end-to-end test in the suite: an
unrelated program reads all three files back, audits each, and cross-checks
them.

**Phase 2, commit 12 complete.** `Partition_Metrics.hh/.cc` — `Kind_Metrics`,
`Metrics`, `measure`, `report`, `report_csv[_header]`, `write_comm_matrix`;
`ume_partition --metrics --metrics-csv --comm-matrix`; four Catch2 tests and
four CTests. All three trees green at 109/109.

Three notes:

- **`mask >= 1` is not "carries work" for a point, edge or face.** On a 4×3×2
  brick only six of sixty points are real; the rest are the mesh's own outer
  surface. An imbalance taken over those measures the brick, not the
  decomposition. So `Kind_Metrics::load` is `local` for the shared kinds —
  `gradzatp` visits every point a rank holds and takes the `mask == -1` branch
  on the boundary ones — and `real` for the ghosted kinds, where a phantom
  zone genuinely is skipped. The table prints `load`, not `real`, and says so.
- `measure` walks the layouts the same way `fill_common` walks them, entity for
  entity, which is what makes it a *prediction* of `extract` rather than an
  estimate. `metrics: the numbers are what extract produces` pins every
  counter — including the send matrix — against the meshes across the whole
  case sweep, so the two cannot drift apart silently.
- `cut_face_points` (102 on the sample) and `point_lambda1` (41) are
  deliberately both reported: a point on several cut faces is charged once by
  the second and once per face by the first, and that gap is exactly the gap
  between what a graph edge-cut minimizes and what goes on the wire. It is the
  reason `metis-kway-cut` is expected to lose to `kahypar` in phase 4.

**Phase 2, commit 13's CMake prerequisite complete.** `add_subdirectory(test)`
moved out of `cmake/pkgs/Catch2Support.cmake` and into the top-level
`CMakeLists.txt` after `add_subdirectory(src)`, plus the first scope-statement
`message(STATUS)` in `test/CMakeLists.txt`. All three trees green at 109/109.

The defect was reproduced rather than reasoned about: with the message in place
but the move not yet made, a fresh `-DUSE_MPI=ON` configure prints `Tests: no
MPI` and then `Found MPI: TRUE` three lines later. After the move the same fresh
configure prints `Tests: MPI is available (.../mpiexec)`, and a fresh configure
with MPI off prints the `no MPI` line. Both branches are therefore known-correct
on a *fresh* tree, which is the case that was broken.

**Phase 2, commit 13's second prerequisite complete.** R10 — `ume_mpi` picked
its spike zone by scanning *upward from the middle of the whole zone range*, so
it found nothing on a mesh whose interior zones are all below the midpoint. It
now takes the middle of the interior zones. `partition_mpi_run` — the first
CTest under the `TARGET ume_mpi` guard, and the first time any test hands a
partition to the code it was written for — fails before and passes after. All
three trees green, 109/109 and 110/110 with MPI.

`gen_brick` is exactly such a mesh: `Builder::make_boundary` appends every
`mask == -1` zone after every real one. In a Debug build the existing assert
fired; in a Release build `czi == local_size()` would have silently spiked a
*ghost* zone, or written past `zfield` on a rank with no ghosts. Found by
running `mpiexec -n 2 ume_mpi` on a partitioned brick, which is the first thing
commit 13 has to do.

**Phase 2, commit 13 complete.** `Global_Ids.hh/.cc` (the `<base>.%05d.gid`
side-car), `ume_partition --global-ids`, `ume_mpi --dump-grad --compare --tol`
with a real option parser, three Catch2 tests and three CTests. All three trees
green, 112/112 and 116/116 with MPI.

The shape is the one design review settled on, and it worked unchanged: a
`--ranks 1` partition is bit-identical to its input, so `mpiexec -n 1 ume_mpi
--dump-grad` *is* the serial answer computed by the same binary along the same
code path — no second reference implementation, and no `ume_serial` change.
`--compare` reads the spike zone's *serial* index back out of the reference, so
both runs perturb the same zone even though no rank could have picked it
locally, and the `.gid` side-car lines the answers up zone for zone afterwards.

Four notes:

- **The comparison is the first test that could catch a decomposition that is
  well formed and still wrong.** The audit, the cross-rank check and the round
  trip all check bookkeeping. On the sample the worst disagreement is 8e-12
  against gradients of order 5e4 — a relative 1.6e-16, one ulp — but the
  tolerance is *relative and nonzero* on purpose: partitioning genuinely
  reorders the accumulation in `gradzatz`, and `test_brick.cc`'s `GRAD_TOL ==
  0.0` would be the wrong test here. It was checked to fail, not assumed to:
  `--tol 0`, a reference whose `spike` line names a different zone, and a
  missing side-car each exit nonzero.
- `Kind`/`kind_name`/`is_shared_kind` moved from `Partition_Metrics.hh` down to
  `Decompose.hh`. The side-car is the third thing to walk the kinds generically
  and the order is `Decomposition::layouts()`'s to fix, not the metrics'.
- No MPI reduction is involved and none is needed. Each rank compares the real
  zones it holds and returns its own verdict; `mpiexec` fails the job if any
  rank does. That the ranks between them cover every zone of the serial mesh is
  `check_partitions`'s job, already done, and is not re-established here — but
  it is visible in the run: 12 + 13 + 9 = 34 real zones on the sample.
- `--dump-grad` refuses to run on more than one rank rather than gathering.
  Gathering would be a second, unchecked path to something one rank already
  produces exactly.

**Phase 3, commit 14 complete.** `centroid_extent`, `dims_create`, and the
`block`/`pencil` methods in `Partition.cc`, plus three Catch2 tests. All three
trees green at 115/115, and 119/119 with MPI.

Three notes:

- The plan's claim that the sweep picks methods up from `Part::methods()` was
  only half true: the *strategy* tests did, but `cases()` — which feeds
  `extract`, the audit, the cross-rank check and the metrics — carried a
  hardcoded `{"linear", "slab-x", "slab-y", "slab-z"}`. It now iterates
  `Part::methods()`, so the rest of phase 3 gets the whole harness for free.
  That is the change to make when adding a method, not a new test.
- **The grid cuts do not have to line up between columns, and it matters that
  they don't.** Each level splits its own group at a cumulative-*weight*
  boundary, so a column of heavier zones gets shorter boxes. The parts are
  therefore not boxes in general — which costs nothing, because nothing
  downstream assumes they are, and it is what keeps the existing "no part is
  more than one zone off the mean" contract test passing at three nesting
  levels.
- `dims_create` works off the *centroid* extent, so an axis one cell thick has
  extent 0 and is never divided while any other axis has extent left. That is
  the case a `dims[i]/extent[i]` rule would divide by zero on, and the reason
  the comparison is cross-multiplied rather than divided.

On an 8×6×4 brick split eight ways, `point_lambda1` is 275 for `slab-x` and 198
for both `block` and `pencil` — the two agree here because the shape rule gives
(4,2,1) either way, z being too thin to be worth cutting.

Both new methods live in `Partition.cc` alongside `slab` and `linear` rather
than in the `Partition_Geometric.cc` the architecture section names. Splitting
the file needs a private header to reach the classes from `make()`; do it when
`morton`/`hilbert` land and the file has earned it.

**Phase 3, commit 15 complete.** `principal_axis` (weighted centroid
covariance, cyclic Jacobi) and the `rcb`/`rib` methods in `Partition.cc`, plus
three Catch2 tests. All three trees green at 118/118, and 122/122 with MPI.

Three notes:

- **The uneven half of an odd `k` is taken by splitting `k` ways and keeping
  the first `k1` groups**, rather than by aiming `split_sorted` at a weight
  fraction. That buys two things at once: the halves come out in the ratio
  `k1 : k - k1` by weight, and because `split_sorted` gives every one of its
  `k` groups a zone whenever there are `k` to give, each half reaches its own
  recursion holding at least as many zones as it has parts to fill. Aiming at
  a fraction directly would need that second property re-established by hand at
  every level.
- **`rib` reproduces `rcb` exactly only while the pieces are still boxes.** The
  obvious test — a 4×3×2 brick, every part count — fails, and is right to: a
  cut that falls *inside* a column of zones leaves a piece with a step in it,
  and a stepped piece has a genuinely off-axis principal direction. The test
  that survives is a cube split by powers of two, and what it pins is the eigen
  solver: an already-diagonal covariance has to come back as a clean coordinate
  direction, ties included, or the projections would differ in the last bit and
  the tie-breaking would wander.
- No `Brick::Spec` knob rotates a mesh, so the test that shows `rib` earning
  its place rotates `pcoord` itself, after `Brick::build` and before anything
  derived is read — the Datastore never invalidates. A 12×2×2 brick turned 45°
  about z has equal x and y extents, so every axis `rcb` can choose cuts it
  diagonally, and `rib` finds the direction the mesh is actually long in.

On a pure-hex 8×6×4 brick split eight ways, `rcb` and `rib` both give
`point_lambda1` 165 — identical to `block`, because recursive bisection of a
box by powers of two *is* the grid — against 245 for `slab-x`.

**Phase 3, commit 16 complete.** `CURVE_BITS`, `morton_index`,
`hilbert_index` (Skilling's transpose) and the `morton`/`hilbert` methods in
`Partition.cc`, plus three Catch2 tests. All three trees green at 121/121, and
125/125 with MPI.

Four notes:

- **`bits` is a parameter of both index functions, and that is what makes them
  testable.** At the production 21 bits the lattice has 2^63 cells, so the
  defining properties — bijection onto `[0, 8^bits)`, and consecutive indices
  always one lattice step apart — can only be checked by enumeration on a small
  lattice. The test does `bits` ∈ {1,2,3} and pins Morton's worst step at
  `2*n-1` in the same loop, so the one property the two curves differ in is a
  measured number rather than a claim.
- `morton_index` interleaves with **x in the high bit of each triple**, not the
  usual low one. That is deliberate: it is exactly the order Skilling's
  transpose has to be read out in, so `hilbert_index` ends with
  `return morton_index(x, bits)` and there is one interleave rather than two.
- **The key handed to `split_sorted` is the zone's *position* along the curve,
  not its index.** A 63-bit index does not survive a `double`, and only the
  order matters. Sorting `(index, zone)` lexicographically also gives the
  tie-break the rest of the file uses.
- The headline claim needed weakening against measurement, twice. Hilbert beats
  Morton **over a sweep of part counts, not pointwise** — on an 8-cube it wins
  at 3, 5, 7, 9, 10, 11, ties at 12, and *loses* at 6 (345 vs 337); summed over
  {3,5,6,7,9,10,11,12} it is 2993 against 3226. And a curve only beats a slab
  stack from **four parts up**: three parts of a cube is two cuts, and the two
  cheapest cuts available are the slab's own planes (182, against Hilbert's
  200). Both are in the test as written, with the reason.

At a part count that divides the curve into whole blocks the two curves agree
*exactly* — on a 4-cube split eight ways both produce the same eight zone sets
`block` does, which is the sharpest available check that the quantization,
ordering and run-splitting line up.

**Phase 3, commit 17 complete.** The `roundrobin` and `random` methods in
`Partition.cc`, plus four Catch2 tests. All three trees green at 125/125, and
129/129 with MPI.

Four notes:

- **"Round-robin makes the communication graph complete" is false, and the
  exception is structural rather than lucky.** The deal has period `nparts`;
  two zones of a brick share a face only if their positions differ by 1, `nx`
  or `nx*ny`; so whenever `nparts` divides one of those strides the deal
  *aliases* with the mesh. On a 4-wide brick at four parts `roundrobin`
  produces `slab-x`'s zone map exactly — not approximately, the identical map —
  and beats `block` on neighbor count (2 against 3). Both the aliased and the
  unaliased case are in the tests, with the arithmetic. A worst case that is
  only a worst case on average has to be swept, not sampled.
- `random` has no period to alias with, and it is the one that reliably reaches
  the complete graph: `mean_neighbors == nparts - 1` at every count measured.
  `roundrobin` is the worse of the two on *volume* when it does not alias — on
  5×4×4 at eight parts it cuts 184 of 184 interior faces, against `random`'s
  158 and `block`'s 72. The two are worst cases of different halves of the cost
  model, which is why both are worth having.
- **The shuffle is hand-rolled — SplitMix64 and Fisher–Yates — on purpose.**
  Neither `std::shuffle` nor `std::uniform_int_distribution` has a specified
  algorithm, so the same seed gives different partitions under libstdc++ and
  libc++, and `--seed` would only reproduce a run on the machine that made it.
  No test here can check that, which is exactly why it is written down.
- `roundrobin` is the first method that does not go through `split_sorted`: it
  balances the *count* exactly and the weight only incidentally. It still meets
  the existing `the split balances the weight it was given` bound on the mixed
  brick, so nothing needed relaxing — but that is a measurement, not a
  guarantee, and a mesh with a strong periodic weight pattern could break it.

One pre-existing defect found, fixed in its own commit afterwards: the doc
comment on `Metrics::interior_faces` said "whose two zones are on one rank",
but `face_cuts` sets `interior` on every face with a real zone on each side and
`cut` on the subset of those that are cut. The number, `report()` and the CSV
were right; only the comment was wrong. It cost a failing test to find, and the
test it cost is the one that now pins the property the comment got backwards —
`interior_faces` is the same for every method, because it is the mesh's
denominator rather than the decomposition's numerator.

**Phase 3, commit 18 complete.** `Zone_Graph`/`build_zone_graph` in
`Partition_Graph.hh/.cc` — the CSR dual graph, with all three `Edge_Weight`
choices behind it — and the `bfs` method in `Partition.cc`, plus five Catch2
tests. All three trees green at 130/130, and 134/134 with MPI.

Five notes:

- **The growth rule went through three shapes, and the two that failed failed
  structurally rather than by a bad constant.** Growing every part at once
  from farthest-point-sampled seeds — which is what this plan specifies — is
  the worst of them: regions grown together race each other and meet along the
  perpendicular bisector of their seeds, and for a box seeded at its two far
  corners that is the diagonal plane, the most expensive cut it has. On an
  8×6×4 brick at two parts it cut 59 faces against `block`'s 24. Carving one
  region at a time fixed two parts; doing it by recursive bisection rather
  than one part after another fixed the rest. `bfs` now ties `block`
  *exactly* — cut faces and duplicated points both — at 2, 3 and 8 parts, and
  is within a third of it at 4 and 16.
- **The gain is `2 * pull - degree`, and the degree is counted inside the
  current node only.** Ranking a candidate by how much edge weight ties it to
  the region, without subtracting what it would newly expose, takes
  well-connected vertices however much boundary comes with them; that alone
  was 96 cut faces against 77 at three parts. Counting the degree over the
  whole graph is the other half: an edge leaving the node was cut at a level
  above and cannot be cut again, so charging for it makes a vertex on the
  node's own rim look expensive when it is the cheapest thing there is.
- **The balance contract caught the first version, and the fix is where each
  split aims.** A part walled in by its neighbours simply stops growing and
  ends light — 16 against a mean of 30 on the mixed brick at eight parts, well
  outside the one-zone slop every other method keeps. Aiming each split at the
  *cumulative* weight that should have been placed by the end of its low half,
  rather than at a share of the node in front of it, puts every one of those
  boundaries within one zone's weight above its ideal; a part's weight is the
  difference of two of them, so it lands within one zone's weight of the mean.
  Rounding accumulates into the last part rather than compounding down the
  tree.
- `--edge-weight` was parsed and advertised from commit 8 and did nothing:
  `bfs` is the first method to use an edge weight at all. It earns the option
  end to end — on a brick of 1×2×4 cells at four parts, `face-area` cuts 66
  faces where `unit` and `face-points` both cut 70.
- **A grown part is one piece in the common case and not by theorem.** Each
  bisection carves a region out of the piece it was given, and the *remainder*
  of that piece can fall in two, after which the next region down that branch
  restarts. Measured over 2..32 parts of an 8×6×4 brick: every part is one
  piece except at 13, 14, 22, 29 and 31, where one part is two. The test pins
  both the property and its limit, because a method that quietly stopped being
  contiguous would pass everything else here.

`bfs` lives in `Partition.cc` alongside every other built-in method, not in the
`Partition_BFS.cc` the architecture section names. The file is now around 1100
lines and the case for splitting it is a private header to reach the classes
from `make()`; that is a refactor of its own, not part of a method commit.

`build_zone_graph`'s arc count is checked against `Metrics::interior_faces`,
which walks the faces and shares no code with it. That they agree also says
that no zone pair in any of the sweep's meshes shares two faces — the case the
edge merge exists for, and the one that would make the two numbers differ.

Next: phase 4, `metis-kway-cut`/`metis-kway-vol`/`metis-rb` behind
`METISSupport.cmake`. `Zone_Graph` is already in the shape `METIS_PartGraphKway`
wants, so that commit is the CMake module, the `idx_t` marshalling, and the
registry entries.

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
| `bfs` | recursive bisection of the dual graph, each cut grown from a seed on the rim of the piece being split and taking whichever zone shortens the boundary most. Gives a *non-geometric, topology-respecting* partition with zero dependencies, and is the natural seed for a later Kernighan–Lin pass. (This started as "grow every part at once from farthest-point-sampled seeds"; see commit 18 for why that does not work) |

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
4. **Numerical equivalence under MPI.** CTest guarded by `TARGET ume_mpi` +
   `MPIEXEC_EXECUTABLE` — *not* by `USE_MPI`, and not from where
   `test/CMakeLists.txt` sits today; see "Commit 13 needs a CMake fix first"
   under Phasing: `gen_brick`, `ume_partition --ranks 2`, `mpirun -n 2
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

These three trees are what *this* development environment can cover, not what a
user is expected to have. Options are enabled where the hardware and the site
support them, so a build with none of MPI, OpenACC or the partitioning libraries
is an ordinary configuration and must be green on its own terms. Nothing in this
plan may make a test's *presence* depend on an option without also making its
absence visible — see "Commit 13 needs a CMake fix first" under Phase 2.

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

#### Commit 13 needs a CMake fix first, in its own commit

**Which options a tree is built with is a property of the machine, not of the
project.** MPI, `UME_SERIAL`, OpenACC, and later METIS and Mt-KaHyPar are each
built by the users whose hardware and site support them, so *every* combination
is a configuration somebody runs, and every one of them has to be green and has
to say what it covered. Commit 13 is the first test that is not universal, so
it is where that gets settled.

The obvious guard does not work, and fails in the worst available way. Three
facts, all verified rather than reasoned about:

- **An unguarded `add_test` fails; it does not skip.** `add_test(NAME x COMMAND
  ume_mpi ...)` in a `USE_MPI=OFF` tree configures cleanly — CMake resolves
  target names in `COMMAND` at generate time and passes an unknown name through
  as a literal — and then ctest says `Unable to find executable: ume_mpi` and
  reports `***Not Run`. A red suite on every build that does not have MPI.
- **`test/` is processed before `find_package(MPI)` runs.** The top-level
  `CMakeLists.txt` has `include(pkgs/Catch2Support)` at line 31 and
  `include(pkgs/MPISupport)` at line 34, and it is `Catch2Support.cmake` that
  calls `add_subdirectory(test)`. So inside `test/CMakeLists.txt`:

  | | `USE_MPI` | `MPIEXEC_EXECUTABLE` |
  |---|---|---|
  | fresh `-DUSE_MPI=ON` | `ON` | **empty** |
  | reconfigure of that same tree | `ON` | populated |

  `USE_MPI` survives because `-D` writes it to the cache before any
  `CMakeLists.txt` runs; `MPIEXEC_EXECUTABLE` does not exist yet on the first
  pass and is only read back from the cache on the second. So
  `if(USE_MPI AND MPIEXEC_EXECUTABLE)` — what this plan said above — silently
  drops the MPI tests on a fresh CI configure and quietly enables them on a
  developer's rebuild. Absent-and-green is worse than either failing or
  skipping, because nobody goes looking.
- **`if(TARGET ume_mpi)` is always false there** for the same reason:
  `add_subdirectory(src)` is line 38.

Nothing about this is specific to MPI. `Catch2Support` is line 31 and
`add_subdirectory(src)` is line 38, so *every* variable and target that
`KokkosSupport`, `MPISupport`, `OpenACCSupport` and `libunwindSupport` produce
is invisible to `test/CMakeLists.txt` on a fresh configure and visible on the
next one. (`UME_SERIAL` is the exception, and only by luck: it is set at line 25,
before Catch2.) Phase 4's `METISSupport` and `KaHyParSupport` will want to gate
tests the same way, so fixing it here fixes it for them.

**The fix.** Move `add_subdirectory(test)` out of `cmake/pkgs/Catch2Support.cmake`
and into the top-level `CMakeLists.txt` after `add_subdirectory(src)`. A `pkgs/`
module adding a subdirectory is the actual layering violation — every other one
only finds a package — and moving it makes every optional feature's variables
right on a fresh configure *and* makes target checks usable. Then, in
`test/CMakeLists.txt`, gate on the target rather than on the option, because the
target existing is what the test actually needs:

```cmake
if(TARGET ume_mpi AND MPIEXEC_EXECUTABLE)
  # ... the numerical-equivalence tests ...
else()
  message(STATUS "No MPI: skipping the partitioned-gradient equivalence tests")
endif()
```

Register conditionally rather than registering always and skipping with
`SKIP_RETURN_CODE`: a test that cannot run on this machine is not a pending
test, and not existing is already what the suite does for optional things. What
the `message(STATUS ...)` buys is that a green run states its own scope — the
one thing a user comparing their results against someone else's needs, and the
thing a silent `if()` alone would take away. Every optional-feature guard added
from here on gets one.

This is a pre-existing defect on the R0–R7 pattern, so it lands as its own
commit ahead of the feature, with the `message` lines as its visible effect.

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
