# Partitioning a mesh across ranks

`ume_partition` reads one serial `.ume` file and writes N per-rank `.ume` files
that together hold the same mesh.  Which zone goes to which rank is chosen by
an interchangeable *method* — eighteen families at the time of writing, or
twenty-one names counting the per-axis `slab` variants, from `slab` through
space-filling curves to METIS and Mt-KaHyPar — and what each choice costs is
measured in the terms the kernels actually pay.

It exists because `gen_brick` produces a single-rank mesh: `mype = 0` (the
number of this rank), `numpe = 1` (how many there are), every entity
`INTERNAL`, empty neighbor lists.  Before this the only way to get a multi-rank
`.ume` set was a DELFI `UmeDump` through `txt2bin`, so the parallel
gather/scatter paths — the whole point of UME (Unstructured Mesh Explorations)
as a proxy app — could not be exercised on a synthetic mesh at all.

Source: `src/Ume/Partition.{hh,cc}` and `Partition_{Geometric,BFS,METIS,
KaHyPar}.cc` (the methods), `src/Ume/Partition_Graph.{hh,cc}` (the dual graph
and the hypergraph), `src/Ume/Decompose.{hh,cc}` (the mesh splitting),
`src/Ume/Partition_Metrics.{hh,cc}` (what it cost),
`src/Ume/partition_check.{hh,cc}` (cross-rank verification),
`src/Ume/Global_Ids.{hh,cc}` (the side-car), `src/Ume/Comm_Loopback.{hh,cc}`
(N ranks in one process), `src/partition.cc` (the command-line interface, or
CLI), `test/test_partition.cc` and `test/test_loopback.cc` (tests).

`doc/partitioning-plan.md` is the working record of how this was built, what
each commit measured, and which claims had to be weakened against
measurement.  This file describes what is there now.

## 1. Vocabulary and notation

The entity kinds and the two naming conventions come from `doc/brick_mesh.md`
§1 and are used here without restating them: **point**, **edge**, **face** and
**zone**, plus the subzonal **corner** and **side**; connectivity maps named
`m:a>b`, read "the map from `a` to `b`"; and masks named `pmask`, `zmask` and
so on, which classify an entity as interior, on the physical boundary, or a
*ghost*.  *Kernel* means one of the library's mesh operators, in practice the
gradient operators `gradzatp` and `gradzatz`.

UME runs in parallel under the Message Passing Interface (MPI) [MPI Forum
2023], which numbers the processes of a run `0..N-1`.  Each is a **rank**, and
each holds one piece of the mesh.  Everything below is about deciding what
those pieces are and what the split costs.

| term | what it means here |
|---|---|
| part | one of the N pieces the mesh is cut into.  A part becomes a rank; which part becomes which rank is `--topology`'s question, and only that |
| zone map | the array a method produces: one part number per zone of the serial mesh.  It is the entire output of a method (§3) |
| real | an entity the mesh counts as genuine: `mask >= 1`.  This is the sense the `real` column of the metrics uses |
| boundary | `mask == -1`.  A point, edge or face on the physical boundary is a genuine entity the kernels still visit, on a branch of their own; a zone, corner or side masked `-1` is a phantom standing in for the outside of the mesh (`doc/brick_mesh.md` §6).  Neither is counted `real` |
| ghost | a read-only local copy of a zone, corner or side owned by another rank, carried so that local operators have complete stencils.  Ghosts sit at the end of the local index range with `mask == 0` |
| phantom | anything the kernels skip: `mask < 1`, which covers ghosts and the boundary zones, corners and sides |
| shared | a point, edge or face that lies on a cut and so exists, as a real entity, on more than one rank |
| SOURCE, COPY | of the several instances of a shared entity, one is the SOURCE and holds the authoritative value; the rest are COPYs |
| `gathscat` | gather–scatter: the exchange that reconciles the instances of a shared entity — every COPY sends its contribution to the SOURCE, which combines them (here with `SUM`) and sends the result back |
| element | one entity's worth of data in one message: the unit the communication metrics count, independent of what is being exchanged.  Its size depends on the field — 8 bytes for a `DBLV` (scalar `double`) field, 24 for a `VEC3V` (3-vector) one |
| cut | the surface between two parts.  A *cut face* is an interior face whose two zones went to different parts |
| dual graph | the graph a mesh induces for graph partitioners: a vertex per real zone, an edge per interior face |
| hypergraph, net | the hypergraph a mesh induces for hypergraph partitioners: a vertex per real zone, and a *net* (hyperedge) per point, holding every zone that touches that point |
| λ(*p*) | the *connectivity* of entity *p*: how many parts hold it.  Σ(λ−1) summed over an entity kind is the number of elements that kind puts on the wire per exchange [Çatalyürek and Aykanat 1999] |

## 2. Using it

```
ume_partition [options] --ranks N <input.ume> <output-base>

  --ranks N            Number of partitions (required)
  --method NAME        Partitioning method [slab]
  --list-methods       Print the methods this build knows about, and exit
  --vertex-weight W    unit | corners | sides [corners]
  --edge-weight W      unit | face-points | face-area [face-points]
  --imbalance F        Allowed heaviest-part to average-part ratio [1.03]
  --seed N             Seed for methods that make random choices
  --ghosts D           none | face | point [face]
  --owner R            low | balanced | hash [low]
  --topology SPEC      Renumber the parts onto a machine, given as x-separated
                       group sizes, leaves first: 4x2 is four ranks per node
                       and two nodes per switch
  --zone-map FILE      Read a zone->part map instead of running a method
  --dump-zone-map FILE Write the zone->part map the method produced
  --metrics            Report what the decomposition costs
  --metrics-csv FILE   Append one row of metrics to FILE, with a header if new
  --comm-matrix FILE   Write the rank-to-rank element counts to FILE
  --global-ids         Also write <output-base>.NNNNN.gid
  --audit              Run mesh_audit on each partition before writing
  --check              Check the partitions against each other afterwards
  --no-write           Do everything except write the files
```

The whole loop, from nothing to a checked four-rank mesh set:

```sh
gen_brick --tet-fraction 0.4 --seed 7 8 6 4 brick.ume
ume_partition --ranks 4 --method block --metrics --check --global-ids \
    brick.ume part
mpiexec -n 4 ume_mpi part
```

which leaves `part.0000{0,1,2,3}.ume` and their `.gid` side-cars.  With
`--no-write` the output base is omitted and nothing is written, which is the
form a method sweep uses.

From C++, skipping files entirely — the two layers are separable and the tests
use them apart:

```c++
#include "Ume/Partition.hh"
#include "Ume/Decompose.hh"

Ume::Part::Options opt;
opt.nparts = 4;

Ume::Part::Zone_Map z2r;
Ume::Part::make("block")->assign(src, opt, z2r);

auto const d = Ume::Part::build_decomposition(src, z2r, opt.nparts);

Ume::SOA_Idx::Mesh rank2;
Ume::attach_serial_transport(rank2);
Ume::Part::extract(src, d, 2, rank2);
```

`build_decomposition` holds every rank at once — rank R's `mySrcs` needs to
know what local index a shared entity has on rank P — but `extract` produces
one mesh at a time, so the driver never holds more than one `Mesh`.  Only
`--check` breaks that, and says so in its help text.

## 3. The parallel model

**A `Partitioner` answers exactly one question: which part does each real zone
belong to?**  Everything downstream — ownership of the points, edges and faces
on a cut, the ghost layer, the local renumbering, the neighbor lists — is a
pure function of the resulting `Zone_Map`.  Keeping the contract that narrow is
what makes the methods interchangeable, and it is why a zone map can be dumped
to a file and read back so the two layers can be tested apart.

### Two mechanisms, not one

Zones, corners and sides are **strictly partitioned**: each belongs to exactly
one rank, and a rank imports a layer of its neighbors' zones — with their
corners, and the sides on the cut — as *ghosts*, in the trailing part of the
index range with `mask == 0`.

Points, edges and faces on a cut are **not ghosts**.  They are duplicated into
the real index range of every rank that touches them; one instance is
designated the SOURCE and the rest are COPYs, and the kernels reconcile them
with `gathscat`.  That distinction is the one the file format was built around
— `cpy_idx`/`src_pe`/`src_idx`, which say where each imported entity comes
from, are indexed by *ghost offset*, while the neighbor lists `myCpys`/`mySrcs`
carry the shared entities — and a decomposition needs both.

A consequence that is easy to get backwards: a shared point a rank holds but
does not source is a COPY *in the real index range*, not a ghost.  `owner[g] !=
r` therefore does not mean "ghost on r".

### Why ghost zones exist at all

The kernels do not need them for correctness.  `gradzatz` accumulates over
corners with `cmask >= 1` and then does `points.gathscat(SUM, ...)`, so if
ranks A and B each own some zones around point *p*, each contributes its own
corners and the SUM reproduces the serial value exactly.

Two other things force ghost sides to exist anyway.  `calc_face_area` writes
`d_side_tag(s2)` unconditionally, so `m:s>s2` must be in range for every real
side including one on a cut face; and `Auditor::ranges` requires `m:c>p` in
range and `m:c>z >= 0` for *every* corner, ghosts included.

`--ghosts none` turns the layer off.  It is a diagnostic: it breaks
`calc_face_area` and trips the audit, and it is in the tests precisely to show
that the ghost machinery is load-bearing rather than merely present.

### What a rank imports

`--ghosts face` (the default) imports every zone sharing a face with a local
zone; `--ghosts point` every zone sharing a point.  In both cases the rule
applied is the same one: **import every side of a ghost zone whose face this
rank already holds, and every corner of a ghost zone whose point this rank
already holds.**  For face adjacency that is exactly "the ring on each cut
face", it extends to point adjacency without a second rule, and it makes "no
ghost points, edges or faces" true by construction.

### Who sources a shared entity

`--owner low` gives it to the lowest-numbered rank holding it: deterministic
and symmetric, so no rank has to negotiate.  It does concentrate SOURCE duty on
low ranks, which is why `--owner balanced` (whichever candidate is so far
sourcing the fewest of that kind) and `--owner hash` are there to measure
against.

### Which rank a part is given

The method decides what goes together; it says nothing about where the parts
sit on a machine, and the launcher numbers ranks `0..N-1` however it likes.
`--topology` renames the parts so that heavy neighbors land close together,
modelling the machine as a tree of equal-sized groups over consecutive rank
numbers — `4x2` is four ranks per node and two nodes per switch, and a
trailing level covering everything is implied.  The cost it minimizes is
**element-hops**: the elements one exchange sweep puts on the wire, each
weighted by how many levels of the tree it has to climb.

This is the classical topology-mapping problem, and the classical heuristic for
it [Hoefler and Snir 2011]: seed with the busiest part, place each remaining
part next to whatever it is most tied to, then swap pairs of parts as long as a
swap lowers the cost.  Stated exactly, it is a quadratic assignment problem
[Koopmans and Beckmann 1957], which is NP-hard; what is implemented is a
heuristic for it and nothing more (§9).

The machine model is coarse deliberately.  A real network is not a balanced
tree, its distances are not integers, and the rank-to-node mapping is a
scheduler policy rather than a property of the hardware, so a more detailed
model would be more precisely wrong.  What the tree does capture is the
distinction that dominates: traffic that stays inside a node against traffic
that does not.

**Renaming the ranks is not quite only a renaming**, and it is worth knowing
why.  The zones that were together are still together, so `point lambda-1`, the
cut faces, the load balance and the total on the wire are all identical.  But
`--owner low` names the SOURCE of a shared entity by rank *number*, so an entity
held by three ranks is a star centred on the lowest-numbered of them, and
renaming moves the centre: the comm matrix changes shape while its total does
not.  So the placement is chosen from one matrix and then rebuilt, re-measured,
and kept only if the measured cost really fell.  The two figures printed are
both measurements.

**Most methods here have nothing to gain from it**, which is itself the
result.  On an 8×6×4 `--tet-fraction 0.4` brick at sixteen ranks and `4x2`,
`slab-x`, `block`, `hilbert`, `bfs-kl`, `metis-*` and `kahypar-det` all keep the
numbering they produced: a method that splits recursively numbers its parts so
that a contiguous block of ranks is a contiguous piece of mesh, which is the
machine tree's own shape.  Where it pays is where the numbering means something
else — `linear`, which numbers by zone index, goes from 5150 element-hops to
4063, and `random` from 19099 to 17452.

### A cut is not a boundary

Points on a cut keep their interior mask, and no phantom sides are invented for
cut faces.  Marking them would make `VAR_point_norm` and the boundary branch of
`gradzatp` subtract an outward normal that is not there, and the gradient would
be quietly wrong along every partition boundary.

## 4. Conventions

**File naming.**  `<output-base>.%05d.ume`, and `<output-base>.%05d.gid` for
the side-car.  `ume_mpi` and `ume_serial` both take the base and append the
suffix themselves.

**Local index order** is `(is_ghost, global index)` lexicographic for every
entity kind: locals first in ascending serial order, then ghosts in ascending
serial order.  It has to be lexicographic rather than plain ascending, because
a ghost zone is an ordinary real zone somewhere else and its global index
interleaves with the local ones.  Two things follow.  `mask == 0` and
`index >= local_size()` mean the same thing, which the library assumes
throughout — and at `--ranks 1` the map is the identity, so a one-rank
partition is **bit-identical to its input**.  That is what makes
`mpiexec -n 1 ume_mpi --dump-grad` a serial reference computed by the same
binary along the same code path, rather than a second implementation to keep in
step.

**Only stored maps are written.**  Everything derived — `m:z>c`, `zcoord`,
`point_norm` and the rest — is left for the Datastore (`doc/datastore.md`) to
work out, so a partition is built by exactly the code paths a mesh read from a
file uses.

**The zone map** (`--dump-zone-map`, `--zone-map`) is one integer per line, one
line per zone of the serial mesh in index order: the part number for a real
zone, `-1` for one that is not real.  A boundary zone follows the real zone it
hangs off, which is the decomposition's decision rather than a method's.

**The `.gid` side-car** records the global identifier (GID) of each local
entity — where it sat in the serial mesh.  It is text, one section per entity
kind, and is written only when asked for, because the `.ume` format has no
place to record such a thing and adding one would change what every existing
file means:

```
ume_gid 1
mype 1
numpe 4
point 107 0
18 27 28 29 30 31 36 37 38 39
40 45 46 47 48 49 54 55 56 57
...
```

Each header is the kind name, the number of real entities and the number of
ghosts; the indices follow in local index order, so the first `local` of them
are the real ones.  Nothing in the library depends on the side-car — it exists
so a partitioned run and a serial run can be lined up entity for entity, and
being text is deliberate, because a side-car is looked at when a comparison
disagrees.

**The comm matrix** (`--comm-matrix`) is `--ranks` comma-separated rows: row
`from`, column `to`, elements sent per exchange sweep over every kind.

```
0,302,180,7
191,0,13,172
111,9,0,248
0,118,166,0
```

It is not symmetric, and should not be: an element is charged to the direction
it moves in.

## 5. The methods

`--list-methods` reports what *this* build can run; the table below is what the
source knows about.  A misspelled name and one whose library is missing are
reported differently, so a sweep script can tell "you typed it wrong" from
"this machine does not have it".

The `based on` column cites the algorithm each method implements; full entries
are in the References.  A dash means the method is a naming or an ordering
convention rather than a published algorithm.

| method | what it does | based on | needs |
|---|---|---|---|
| `linear` | contiguous ranges of the input zone numbering | — | — |
| `slab`, `slab-x/y/z` | slabs normal to the longest extent, or to a named axis | — | — |
| `block` | a grid of boxes, shaped to the extents of the mesh | — | — |
| `pencil` | columns spanning the shortest extent | — | — |
| `rcb` | recursive coordinate bisection (RCB): each cut normal to the piece's longest axis | [Berger and Bokhari 1987], [Williams 1991] | — |
| `rib` | recursive inertial bisection (RIB): each cut normal to the piece's principal axis, the dominant eigenvector of its coordinate covariance | [Simon 1991] | — |
| `morton` | equal runs of a Morton (Z-order) curve | [Morton 1966], [Warren and Salmon 1993] | — |
| `hilbert` | equal runs of a Hilbert curve | [Hilbert 1891], [Skilling 2004] | — |
| `roundrobin` | the i'th real zone to part `i % nparts` | — | — |
| `random` | a shuffle determined by the seed | [Durstenfeld 1964], [Steele et al. 2014] | — |
| `bfs` | parts grown outwards along the dual graph, each from a seed a breadth-first search (BFS) puts out on the rim of what is left | [Farhat 1988] | — |
| `bfs-kl` | the same, with each cut refined by Kernighan–Lin (KL) swaps | [Kernighan and Lin 1970] | — |
| `metis-kway-cut` | METIS k-way, minimizing the weight of the cut faces | [Karypis and Kumar 1998b] | METIS |
| `metis-kway-vol` | METIS k-way, minimizing its communication-volume proxy | [Karypis 2013] | METIS |
| `metis-rb` | METIS recursive bisection of the dual graph | [Karypis and Kumar 1998a] | METIS |
| `kahypar` | Mt-KaHyPar, minimizing the points on the wire | [Gottesbüren et al. 2021] | Mt-KaHyPar |
| `kahypar-quality` | the same with flow-based refinement; slower | [Gottesbüren et al. 2024] | Mt-KaHyPar |
| `kahypar-det` | the same, reproducibly, at some cost in quality | [Gottesbüren and Hamann 2022] | Mt-KaHyPar |

The library methods are *multilevel*: they coarsen the graph or hypergraph by
merging vertices, partition the small version, and then uncoarsen it, refining
the cut at each level.  That scheme is why they beat the direct methods on
quality, and [Buluç et al. 2016] is a readable survey of it.

`bfs` and `bfs-kl` are the only methods that read neither a coordinate nor a
library.  `bfs` grows each half of a bisection outwards and cannot go back;
`bfs-kl` then walks the boundary it left, swapping a zone of each half at a
time.  Pairs rather than single moves, because the balance the built-in methods
promise leaves no room for a move that is not weight-neutral.

What the refinement is worth depends on how much the growth had to get wrong.
Summed over 2, 3, 4, 8 and 16 parts of an 8×6×4 brick at `--tet-fraction 0.4`,
where the zones differ eightfold in weight: 783 duplicated points for `bfs`,
707 for `bfs-kl`, and 718 for `block` — so the tenth it wins back is the
difference between losing to a method that reads the coordinates and beating
one.  On the same brick with no tets it is worth a hundredth (700 against 694),
because a box cut into boxes leaves a greedy little room to go astray.  Neither
figure survives being taken at a single part count; the sweep below is one such
point, and at it `block` is still ahead.

`roundrobin` and `random` are worst cases rather than proposals, and they are
worst cases of *different halves* of the cost model: `random` reliably reaches
the complete communication graph, while `roundrobin` moves more volume when it
does not alias with the mesh.  Both are useful precisely because a method that
is only bad on average has to be swept rather than sampled.

### Determinism

Every method is reproducible from `--seed` with one advertised exception:
`kahypar` and `kahypar-quality` refine in parallel and race, so the same seed on
the same mesh gives a different partition run to run.  `kahypar-det` is the
preset that keeps the promise, by an algorithm designed for it
[Gottesbüren and Hamann 2022].  Determinism did not cost quality in the
measurements taken here — `kahypar-det` ties `kahypar-quality`'s midpoint and
beats the fast preset outright — but that is a measurement, not a guarantee.

`random`'s shuffle is hand-rolled — SplitMix64 [Steele et al. 2014] and
Fisher–Yates [Durstenfeld 1964] — on purpose.  Neither `std::shuffle` nor
`std::uniform_int_distribution` has a specified algorithm, so the same seed
would give different partitions under libstdc++ and libc++, and `--seed` would
only reproduce a run on the machine that made it.

### Balance

Every built-in method splits a sorted order at cumulative-*weight* boundaries,
which makes "no part is more than one zone's weight off the mean" a theorem
about them rather than a hope.  METIS balances to a ratio and reaches it
heuristically, so it can and does miss that bound while honouring
`--imbalance`.  The test file states each family's own promise; do not carry
the stronger one over to a library method.

## 6. How to read the metrics

`--metrics` computes everything from the decomposition and the serial mesh
alone.  No partition mesh is built, so sweeping a dozen methods across a dozen
rank counts with `--no-write` is cheap, and the numbers are available before
deciding whether a split is worth writing out.

```
Decomposition metrics for 4 parts
  kind      serial    real  min load  max load   imbal    ghost  lambda-1  nbrs
  point        315     105        98       107   1.052        0        92     3
  edge         930     498       270       281   1.026        0       165     3
  face         903     679       241       250   1.021        0        76     3
  zone         511     287        70        74   1.031      128       128     3
  corner      2704    1840       456       464   1.009      506       506     3
  side        6560    5696      1420      1432   1.006      550       550     3
  cut faces:           76 of 679 interior, 275 points on them
  neighbors:           max 3, mean 3.00
  surface/volume:      max 0.45, mean 0.43
  point lambda-1:      92, about 8096 bytes on the wire per gradzatz
  elements exchanged:  1517 in one direction of one sweep over every kind
```

**`point lambda-1` is the headline number.**  A point held by λ(*p*) ranks needs
λ(*p*)−1 copies, and each copy is one element in one direction of
`points.gathscat`.  So Σ(λ−1) is not a surrogate for UME's communication — it
*is* the number of point elements on the wire.  `gradzatz` does a `DBLV` and a
`VEC3V` gathscat plus a `VEC3V` scatter, so the byte figure beside it is about
88 Σ(λ−1) — 2×8 for the scalar gathscat, 2×24 for the vector one and 24 for the
one-way scatter — and it is printed because it is a number to put next to a
stopwatch.

**`nbrs` is the latency-bound half**, and deserves equal billing.  A scattered
zone map drives `max_neighbors` to N−1 while a blocked one holds it near 26, at
very similar volume.  In the sweep below `slab` moves 65% more elements than
`block` and yet talks to fewer ranks; which of those costs more is the machine's
answer, not the metric's.

**The `load` columns are not `real`.**  For points, edges and faces, load is
everything a rank holds, because `gradzatp` visits every point a rank holds and
takes the `mask == -1` branch on the boundary ones.  For zones, corners and
sides it is the real count, because a `mask < 1` zone is a phantom the kernels
skip.  The distinction is not pedantry: on a 4×3×2 brick only six of sixty
points are real, so a point imbalance taken over `real` would be measuring the
mesh's surface rather than the decomposition.

**`cut_face_points` and `point lambda-1` are deliberately both reported.**  A
point on several cut faces is charged once per face by the first and once
altogether by the second, and that gap is exactly the gap between what a graph
edge cut minimizes and what actually goes on the wire — the observation the
hypergraph model was introduced to fix [Çatalyürek and Aykanat 1999],
[Hendrickson and Kolda 2000].  It is the reason the hypergraph methods beat the
graph ones.

**`interior_faces` is the mesh's denominator, not the decomposition's
numerator** — faces with a real zone on each side, cut or not.  It is the same
for every method, which is what makes `cut faces: 76 of 679` a comparable
fraction.

### A sweep

`--metrics-csv` appends one row per run, with a header if the file is new, so a
sweep is a shell loop and the result is directly plottable:

```sh
for m in $(ume_partition --list-methods | awk '{print $1}'); do
  ume_partition --ranks 8 --method "$m" --metrics-csv sweep.csv \
      --no-write brick.ume
done
```

On the 8×6×4 brick at `--tet-fraction 0.4` above, eight parts, seven of the
nineteen columns:

| method | zone_imbalance | cut_faces | point_lambda1 | max_nbrs | mean_nbrs | comm_elems |
|---|---|---|---|---|---|---|
| `kahypar-quality` | 1.087 | 134 | 169 | 6 | 3.75 | 2720 |
| `metis-kway-cut` | 1.087 | 129 | 169 | 6 | 4.50 | 2675 |
| `kahypar-det` | 1.087 | 134 | 170 | 5 | 3.50 | 2728 |
| `metis-rb` | 1.115 | 132 | 171 | 5 | 3.50 | 2704 |
| `block`, `pencil` | 1.059 | 145 | 179 | 5 | 3.50 | 2889 |
| `rcb` | 1.059 | 150 | 185 | 5 | 4.00 | 2971 |
| `morton` | 1.087 | 148 | 194 | 7 | 6.00 | 3004 |
| `hilbert` | 1.143 | 155 | 195 | 7 | 5.50 | 3081 |
| `bfs-kl` | 1.115 | 158 | 200 | 7 | 4.75 | 3168 |
| `rib` | 1.087 | 161 | 202 | 6 | 4.75 | 3158 |
| `metis-kway-vol` | 1.087 | 150 | 203 | 6 | 4.00 | 3037 |
| `bfs` | 1.059 | 174 | 225 | 7 | 5.25 | 3492 |
| `slab` | 1.059 | 228 | 260 | 4 | 3.00 | 4514 |
| `linear` | 1.143 | 227 | 293 | 6 | 4.25 | 4692 |
| `roundrobin` | 1.003 | 588 | 927 | 7 | 7.00 | 12597 |
| `random` | 1.059 | 602 | 977 | 7 | 7.00 | 13164 |

Four things in it are worth reading twice.  `bfs-kl` moves three places up the
table from `bfs` and pays for it in balance — 1.115 against 1.059, because the
refinement is allowed to spend the slack the growth left unspent, and a swap of
two zones of different weight spends some.  `metis-kway-vol` optimizes a
volume proxy and loses to `metis-kway-cut` on the volume that matters, on this
mesh and on every other one measured here — its objective is vertex-based, and
the gap between it and Σ(λ−1) is the gap the hypergraph closes.  `roundrobin`
has the *best* balance in the table and costs five times what anything sensible
does.  And the top of the table is close: one mesh at one part count separates
`kahypar-quality` from `block` by 6%, which is why the margin claims in the
tests are made on a larger mesh summed over several part counts, and why the
racing presets are never pinned to a number.

## 7. Verification

Five layers, weakest first.  The first three check bookkeeping; only the last
two could catch a decomposition that is well formed and still wrong.

1. **`--audit`** runs `mesh_audit` on each partition before writing.  Every
   partition of every mesh in the test sweep passes the *full* audit
   unchanged — the only audit change this work needed was teaching the facet
   check to count ghost sides.
2. **`--check`** (`check_partitions`) compares the ranks against each other,
   which the audit structurally cannot: the entity a `src_idx` names lives in a
   mesh the audit never sees.  Four groups, hardest-hitting last — the neighbor
   lists pair up; **the pairing is positional**, because `Comm::Buffers` packs
   and unpacks by position and an agreement on lengths with a disagreement on
   order is a silent permutation of the data on the wire; ghost sources are
   real; and the ranks reassemble into the serial mesh with nothing lost or
   double-counted.
3. **`ume_serial --audit part.00000.ume ...`** does both from an unrelated
   program that reads the files back.  Handed a *partial* set it says so and
   skips the cross-rank check rather than failing it.
4. **The loopback transport** runs every rank of a decomposition in one
   process and compares the partitioned gradient against the serial one.  It
   needs no MPI, no `mpiexec` and no files, so it runs in every build, over
   every method at N ∈ {1,2,3,5}.
5. **`mpiexec -n N ume_mpi --compare`** does the same thing for real.
   `mpiexec -n 1 ume_mpi --dump-grad ref` writes the reference; `--compare ref`
   reads back which zone that run spiked so both runs perturb the same one, and
   the `.gid` side-cars line the answers up zone for zone.

All five audit meshes this project generated, so all five agree with this
project's conventions by construction.  What calibrates the checkers themselves
is the reference data: meshes the reference application wrote and `txt2bin`
converted, where a disagreement means the checker is wrong.  They are too large
to carry in the tree, so the tests that read them are registered only when
asked for:

```sh
cmake -DUME_DATA_DIR=/path/to/UME-data ...
ctest -R reference_audit
```

Every directory under it holding `.ume` files becomes one test, given every
file of that directory at once, so a multi-rank dataset gets the cross-rank
comparison as well as a per-rank audit.  Expect minutes rather than seconds:
these are meshes of tens of millions of sides.

**Both numerical comparisons use a relative, nonzero tolerance, on purpose.**
Partitioning genuinely reorders the accumulation in `gradzatz`, and addition in
floating point is not associative, so a reordered sum is a different sum
[Goldberg 1991].  `gradzatz` also accumulates its point sums with
`Kokkos::atomic_add` [Trott et al. 2022] in every backend but Serial, so even
the *reference* is only reproducible to roundoff.  An exact comparison is wrong
even at one rank, where the partition is the input mesh.  In practice the
disagreement is about one unit in the last place (ulp): the worst seen across
the four build trees is 7.3e-12 against gradients of order 2.0e4.

## 8. Optional libraries

Neither library is required.  A build with neither is an ordinary
configuration, must be green on its own terms, and says what it covered:
`--list-methods` drops the methods it cannot run, and the CTests that need one
are not registered at all, with a `message(STATUS ...)` at configure time
saying so.  The Catch2 tests are a single binary covering both configurations,
so those gate at *runtime* on `Method::available` and ctest reports them as
`Skipped`.  All four configure combinations are built and run before a library
commit lands.

| | METIS | Mt-KaHyPar |
|---|---|---|
| what it is | multilevel graph partitioner [Karypis and Kumar 1998a] | Multi-Threaded Karlsruhe Hypergraph Partitioner [Gottesbüren et al. 2024] |
| version used here | 5.1.0 | 1.6.2 |
| licence | Apache-2.0 | MIT |
| found by | `cmake/pkgs/METISSupport.cmake` (`SimplePkg`) | `cmake/pkgs/KaHyParSupport.cmake` (`find_package(MtKaHyPar CONFIG)`) |
| disabled with | `-DCHECK_METIS_DONE=TRUE -DMETIS_NOT_FOUND=TRUE` | `-DCMAKE_DISABLE_FIND_PACKAGE_MtKaHyPar=ON` |

Both licences are permissive and compatible with UME's BSD-3, and both
libraries are linked rather than vendored — no third-party source is in this
tree.

`KaHyParSupport.cmake` is the one `pkgs/` module that does not use `SimplePkg`:
upstream ships a config package, and `libmtkahypar` links oneTBB (Intel's
oneAPI Threading Building Blocks) transitively, which a bare
header-and-library probe would not carry.  METIS's `SimplePkg` caches its
result in `CHECK_METIS_DONE` and will not re-probe, so a tree configured
without it has to be a throwaway; the Mt-KaHyPar switch needs no such thing.

Mt-KaHyPar is in neither Spack's builtin repo nor Homebrew, so
`scripts/spack-repo/` carries a package for it as a stopgap:

```sh
spack repo add --scope=user scripts/spack-repo/spack_repo/ume
spack install mt-kahypar
```

Its README says why each line of the recipe is there, including why oneTBB
forces the whole stack to clang on macOS.  That is a macOS-only problem —
gcc + oneTBB is the supported combination on Linux — and a clang-built
`libmtkahypar` linked from g++ works, because the C application programming
interface (API) is `extern "C"` over plain-old-data (POD) types.

### Why the hypergraph

METIS partitions the **dual graph**: a vertex per real zone, an edge per
interior face.  Minimizing the weight of the cut edges minimizes cut *faces*,
which is a proxy for what UME moves.

It is a proxy and not the thing itself, and the reason is well known
[Çatalyürek and Aykanat 1999], [Hendrickson and Kolda 2000]: an edge cut counts
each cut face separately, but a point shared by four cut faces is still one
point, sent once.  METIS's own `metis-kway-vol` objective addresses this from
the vertex side, which is why it appears in §6's table optimizing something
adjacent to, but not equal to, the number of elements on the wire.

Mt-KaHyPar partitions a **hypergraph**: a vertex per real zone, and a net per
point holding the zones that touch it.  Its connectivity objective is
Σ(λ−1) over the nets — which, net for net, *is* `point lambda-1`.  That
equality is checked in every build, over every method and part count, by
computing it from the hypergraph and from the decomposition by two paths that
share no code.

`--edge-weight` weighs faces and has no meaning for a hypergraph whose nets are
points, so `build_zone_hypergraph` ignores it and every net weighs 1.  That is
not a simplification: it is precisely what makes the objective *equal* the
element count on the wire rather than proportional to it.

## 9. Limitations

The input must be a serial mesh — one with no ghosts of its own.
Repartitioning an already-decomposed set, and partitioning a mesh too large for
one process, are both out of scope; `build_decomposition` throws rather than
guessing.

`--topology` places the parts on a tree of equal-sized groups and nothing
finer.  A machine whose nodes are not all the same size, or whose network is
not a tree, cannot be described to it, and the placement is a heuristic for a
quadratic assignment problem rather than a solution to one.  At eight parts,
where brute force over all 40320 permutations is affordable, it found the
optimum in fourteen of the sixteen cases checked and came within one percent
in the other two — but nothing promises that.  What is promised is only that
the measured cost never rises.

Nothing tells the launcher about any of this.  The renaming assumes rank *r*
goes where the spec says rank *r* goes, which is `--map-by` on the `mpiexec`
command line and is the user's to get right.

`bfs-kl` refines each bisection as it is made and never revisits it.  A
refinement pass over the *finished* k-way partition — one that could move a
zone between any two parts rather than only across the cut that separated
them — is not implemented, and it is where the rest of the gap to the library
methods is.

`--imbalance` is passed through to the library methods and honoured by them.
The built-in methods balance exactly by construction and ignore it.

## References

**Graph and hypergraph partitioning**

- [Buluç et al. 2016] A. Buluç, H. Meyerhenke, I. Safro, P. Sanders and
  C. Schulz, "Recent advances in graph partitioning," in *Algorithm
  Engineering*, Lecture Notes in Computer Science 9220, 117–158, Springer,
  2016.  A survey of the multilevel scheme the library methods use; the
  readable entry point for this section.
- [Çatalyürek and Aykanat 1999] Ü. V. Çatalyürek and C. Aykanat,
  "Hypergraph-partitioning-based decomposition for parallel sparse-matrix
  vector multiplication," *IEEE Transactions on Parallel and Distributed
  Systems* 10(7), 673–693, 1999.  The Σ(λ−1) connectivity metric, and why an
  edge cut is not a communication volume.
- [Hendrickson and Kolda 2000] B. Hendrickson and T. G. Kolda, "Graph
  partitioning models for parallel computing," *Parallel Computing* 26(12),
  1519–1534, 2000.  What the standard graph model does and does not measure.
- [Kernighan and Lin 1970] B. W. Kernighan and S. Lin, "An efficient heuristic
  procedure for partitioning graphs," *Bell System Technical Journal* 49(2),
  291–307, 1970.  The pairwise-swap refinement `bfs-kl` runs over each cut.
- [Farhat 1988] C. Farhat, "A simple and efficient automatic FEM domain
  decomposer," *Computers & Structures* 28(5), 579–602, 1988.  The
  greedy/breadth-first region growth `bfs` implements.

**Geometric methods**

- [Berger and Bokhari 1987] M. J. Berger and S. H. Bokhari, "A partitioning
  strategy for nonuniform problems on multiprocessors," *IEEE Transactions on
  Computers* C-36(5), 570–580, 1987.  Recursive coordinate bisection.
- [Simon 1991] H. D. Simon, "Partitioning of unstructured problems for parallel
  processing," *Computing Systems in Engineering* 2(2–3), 135–148, 1991.
  Inertial and spectral bisection compared on unstructured meshes.
- [Williams 1991] R. D. Williams, "Performance of dynamic load balancing
  algorithms for unstructured mesh calculations," *Concurrency: Practice and
  Experience* 3(5), 457–481, 1991.  Orthogonal recursive bisection, which
  `rcb` is, measured against spectral bisection and simulated annealing on
  unstructured meshes.
- [Hilbert 1891] D. Hilbert, "Ueber die stetige Abbildung einer Linie auf ein
  Flächenstück," *Mathematische Annalen* 38(3), 459–460, 1891.
- [Morton 1966] G. M. Morton, *A Computer Oriented Geodetic Data Base and a New
  Technique in File Sequencing*, technical report, IBM Ltd., Ottawa, 1966.  The
  Z-order interleaving `morton_index` computes.
- [Skilling 2004] J. Skilling, "Programming the Hilbert curve," *AIP Conference
  Proceedings* 707, 381–387, 2004.  The `AxestoTranspose` algorithm
  `hilbert_index` implements.
- [Warren and Salmon 1993] M. S. Warren and J. K. Salmon, "A parallel hashed
  oct-tree N-body algorithm," in *Proceedings of Supercomputing '93*, 12–21,
  ACM/IEEE, 1993.  Space-filling-curve orderings used as a domain
  decomposition.

**Libraries**

- [Karypis and Kumar 1998a] G. Karypis and V. Kumar, "A fast and high quality
  multilevel scheme for partitioning irregular graphs," *SIAM Journal on
  Scientific Computing* 20(1), 359–392, 1998.  METIS and its multilevel
  recursive bisection.
- [Karypis and Kumar 1998b] G. Karypis and V. Kumar, "Multilevel k-way
  partitioning scheme for irregular graphs," *Journal of Parallel and
  Distributed Computing* 48(1), 96–129, 1998.  The direct k-way formulation
  `metis-kway-*` uses.
- [Karypis 2013] G. Karypis, *METIS: A Software Package for Partitioning
  Unstructured Graphs, Partitioning Meshes, and Computing Fill-Reducing
  Orderings of Sparse Matrices, Version 5.1.0*, University of Minnesota, 2013.
  Where `METIS_OBJTYPE_VOL`, the communication-volume proxy, is defined.
- [Gottesbüren et al. 2021] L. Gottesbüren, T. Heuer, P. Sanders and S. Schlag,
  "Scalable shared-memory hypergraph partitioning," in *Proceedings of ALENEX
  2021*, 16–30, SIAM, 2021.  Mt-KaHyPar's default preset.
- [Gottesbüren and Hamann 2022] L. Gottesbüren and M. Hamann, "Deterministic
  parallel hypergraph partitioning," in *Euro-Par 2022*, Lecture Notes in
  Computer Science 13440, 301–316, Springer, 2022.  The `DETERMINISTIC` preset
  `kahypar-det` selects.
- [Gottesbüren et al. 2024] L. Gottesbüren, T. Heuer, N. Maas, P. Sanders and
  S. Schlag, "Scalable high-quality hypergraph partitioning," *ACM Transactions
  on Algorithms* 20(1), article 9, 2024.  The framework as a whole, including
  the flow-based refinement `kahypar-quality` turns on.

**Task placement**

- [Hoefler and Snir 2011] T. Hoefler and M. Snir, "Generic topology mapping
  strategies for large-scale parallel architectures," in *Proceedings of ICS
  '11*, 75–84, ACM, 2011.  The greedy placement and pairwise-swap refinement
  `--topology` performs.
- [Koopmans and Beckmann 1957] T. C. Koopmans and M. Beckmann, "Assignment
  problems and the location of economic activities," *Econometrica* 25(1),
  53–76, 1957.  The quadratic assignment problem `--topology` approximates.

**Implementation**

- [MPI Forum 2023] Message Passing Interface Forum, *MPI: A Message-Passing
  Interface Standard, Version 4.1*, 2023.
- [Steele et al. 2014] G. L. Steele Jr., D. Lea and C. H. Flood, "Fast
  splittable pseudorandom number generators," in *Proceedings of OOPSLA 2014*,
  453–472, ACM, 2014.  SplitMix64, the generator in `src/Ume/random.hh`.
- [Durstenfeld 1964] R. Durstenfeld, "Algorithm 235: Random permutation,"
  *Communications of the ACM* 7(7), 420, 1964.  The in-place Fisher–Yates
  shuffle `random` uses.
- [Trott et al. 2022] C. R. Trott et al., "Kokkos 3: Programming model
  extensions for the exascale era," *IEEE Transactions on Parallel and
  Distributed Systems* 33(4), 805–817, 2022.
- [Goldberg 1991] D. Goldberg, "What every computer scientist should know
  about floating-point arithmetic," *ACM Computing Surveys* 23(1), 5–48, 1991.
  Why a change in summation order changes the answer.
