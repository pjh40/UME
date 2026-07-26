# Brick mesh generation

`Ume::Brick` builds a `SOA_Idx::Mesh` in memory from a handful of parameters:
a structured `nx * ny * nz` lattice of cells, where each cell is either a
single hexahedral zone or six tetrahedra.  The mixture and the spatial
arrangement of the two cell types are both controllable, so the generator
produces anything from a pure hex brick, through a hex-dominated mesh with
tetrahedral inclusions, to a pure tet mesh.

It exists because the library's operators need meshes to run on: meshes with a
physical boundary, with a topology that is not a uniform hex grid, and with
geometry that is known analytically so results can be checked rather than
merely compared against themselves.

Source: `src/Ume/Brick_Mesh.{hh,cc}` (library), `src/gen_brick.cc` (CLI),
`src/Ume/mesh_audit.{hh,cc}` (mesh-independent consistency checks),
`test/test_brick.cc` (tests).

The generator supplies only topology, `pcoord`, masks and communication
structure.  All geometry — `side_vol`, `side_surf`, `corner_vol`,
`corner_csurf`, `zcoord`, `fcoord`, `ecoord`, `point_norm`, and every inverted
map — is derived lazily by the `VAR_*::init_()` methods in the Datastore, so
the mesh a test gets is built by exactly the code paths a mesh read from a file
would use.

## 1. Using it

From the command line:

```
gen_brick [options] <nx> <ny> <nz> <output.ume>

  --tet-fraction F   Portion of the mesh that is tetrahedral [0]
  --basis zones      F is the fraction of resulting zones (default)
  --basis cells      F is the fraction of lattice cells
  --clustering L     Arrangement of the tet cells, 0 = uniform random
                     scatter, 1 = a single face-connected component [0]
  --seed N           Seed for the cell selection [fixed default]
  --spacing DX DY DZ Cell size [1 1 1]
  --origin X Y Z     Coordinates of the lattice origin [0 0 0]
  --phantom-sides    Use the alternative boundary representation (see §5)
  --validate         Check the mesh invariants before writing
  --no-write         Skip the output file; with --validate, a pure check run
```

The driver writes a binary `.ume` file and then reads it back and compares with
`Mesh::operator==`, so a successful exit also certifies the round trip.

From C++, skipping files entirely:

```c++
#include "Ume/Brick_Mesh.hh"

Ume::SOA_Idx::Mesh mesh;
Ume::Brick::attach_serial_transport(mesh);   // derived vars call scatter()

Ume::Brick::Spec spec;
spec.nx = spec.ny = spec.nz = 4;
spec.tet_fraction = 0.25;
spec.clustering = 0.5;

Ume::Brick::Stats const stats = Ume::Brick::build(spec, mesh);
```

`Spec` also carries `dx/dy/dz`, `x0/y0/z0`, `basis`, `seed` and `boundary`.
`Stats` reports the cell, zone, point, edge, face, corner and side counts, how
many hex faces had to be split, and how many face-connected components the
tetrahedral cells actually formed.

Cell selection is a separate, deterministic step, so a test that needs an exact
pattern can bypass the parameters entirely:

```c++
std::vector<Ume::Brick::Cell_Type> cells(spec.num_cells(), Cell_Type::HEX);
cells[13] = Ume::Brick::Cell_Type::TET;
Ume::Brick::build(spec, cells, mesh);
```

`attach_serial_transport` installs a process-local no-op `Transport`.  A mesh
without one crashes the moment a derived variable calls `scatter()`.

## 2. The two mixture controls

### `tet_fraction` and `basis`

The underlying decision is per lattice *cell*: subdivided or not.  Because a
subdivided cell contributes six zones where a hex cell contributes one, "20%
tets" is ambiguous, and `basis` selects the denominator:

| `basis` | `tet_fraction` means | conversion to a cell fraction |
|---|---|---|
| `ZONES` (default) | fraction of resulting zones | `f_cell = t / (6 - 5t)` |
| `CELLS` | fraction of lattice cells | — |

The cell target is met exactly.  The zone fraction is therefore quantized in
steps of one cell — five zones — so a requested zone fraction of 0.5 on a
4×4×4 brick lands at 0.4954.  Both numbers are reported.

### `clustering`

`T` tetrahedral cells are grown from `K = round(T^(1 - clustering))` random
seeds by Eden-style growth: repeatedly take a uniformly random cell from the
frontier of the marked set and mark it, using 6-connectivity.

| `clustering` | `K` | result |
|---|---|---|
| 0.0 | `T` | every tet cell is its own seed: a uniform random scatter |
| 0.5 | `sqrt(T)` | a handful of blobs of random size |
| 1.0 | 1 | a single face-connected component |

Seeds that collide during growth merge, so the realized component count is at
most `K`.  `Stats::tet_clusters` reports the true count, computed with a
union-find over the selected cells, rather than the number requested.

Selection is fully determined by `seed`, so a given `Spec` always produces the
same mesh.

## 3. Why hex and tet cells can be mixed

Zones in this representation are general polyhedra: a zone is defined by its
corners and sides, not by a fixed vertex count.  A side is the tetrahedron
`(zone center, face center, p1, p2)`, and a zone's volume is exactly the sum of
its side volumes for any set of planar faces.

The constraint is conformity.  A subdivided cell presents two triangles on each
of its six cell faces, so a hex neighbor has to present the same two triangles
— its quadrilateral face must be split.  That is legal (the hex zone becomes a
7- to 12-faced polyhedron of unchanged volume), and it is consistent because of
a property of the Kuhn subdivision:

> Every Kuhn tetrahedron of a cell contains the body diagonal
> `(0,0,0)-(1,1,1)`.  The triangulation this induces on a cell face is
> therefore always the diagonal between local coordinates `(0,0)` and `(1,1)`
> of the two in-face axes taken in ascending order — the same for every cell,
> every face and every orientation.

So the rule is purely local: **a cell face is split if and only if either cell
it bounds is subdivided**, and the two sides of a shared face always choose the
same diagonal.  Faces are registered globally by their sorted point set, which
matches the two sides up regardless of how each zone orients its ring.

Per-cell inventory:

| cell | zones | corners | faces contributed | sides |
|---|---|---|---|---|
| hex, no split faces | 1 | 8 | 6 quadrilaterals | 24 |
| hex, all faces split | 1 | 8 | 12 triangles | 36 |
| subdivided | 6 | 24 | 12 boundary triangles + 6 internal | 72 |

## 4. Conventions

### Orientation

`side_vol` is

```
((p2 - zc) x (p1 - zc)) . (fc - zc) / 6
```

which is positive when a face ring is traversed counterclockwise *as seen from
inside the zone* — equivalently, when the ring's right-hand-rule normal points
at the zone center.  Rather than deriving this per cell type, the builder
orients every ring at runtime against the zone centroid, and validation asserts
`side_vol > 0` on every real side.

`m:s>c1` and `m:s>c2` are the corners at `p1` and `p2` respectively.
`corner_csurf` adds `side_surf` to `c1` and subtracts it from `c2`, and
`side_surf` points from `p1` toward `p2`, so the two agree.

### Side adjacency

| map | meaning |
|---|---|
| `m:s>s2` | the side with the same `(face, edge)` in the zone on the other side of that face.  The only one the library consumes; built by pairing off `(face, edge)` keys. |
| `m:s>s3`, `m:s>s4` | the next and previous sides around the same face of the same zone. |
| `m:s>s5` | the other side of the same zone on the same edge.  An edge of a closed polyhedron lies in exactly two of its faces; a degenerate boundary zone has one face, so its sides refer to themselves. |

`s3`, `s4` and `s5` are written to the file and renumbered, but are not read by
any calculation; the conventions above are this generator's choice.

### Masks

| entity | interior | physical boundary | note |
|---|---|---|---|
| points | 1 | -1 | boundary = belongs to a boundary face |
| edges | 1 | -1 | boundary = belongs to a boundary face |
| faces | 1 | -1 | boundary = bounded by only one real zone |
| zones | 1 | -1 | boundary zones only (see §5) |
| corners | 1 | -1 | boundary zones' corners only (see §5) |
| sides | 1 | -1 | phantom exterior sides only; real sides of a boundary zone stay 1 |

No entity is masked 0: a serial mesh has no ghosts.  `m:f>z1` is the real zone
behind a face; `m:f>z2` is the other zone, or -1 for a boundary face under the
alternative boundary style.

## 5. Boundary representation

The default, `Boundary_Style::BOUNDARY_ZONES`, gives each boundary face a
degenerate zone outside the mesh:

| entity | mask | notes |
|---|---|---|
| boundary zone | -1 | one per boundary face; `zcoord` degenerates to the face centroid |
| its corners | -1 | one per face point; carry no volume |
| its sides | -1 | one per face edge; `m:s>s2` is the paired real side |

Two consumers in the library make this the natural reading:

1. `VAR_side_surf` handles `smask < 0` with a separate formula because "there
   isn't really a zx here".  If a phantom side's `m:s>z` named the adjacent
   real zone there would be a perfectly good `zcoord` to use.  The special case
   only makes sense if the side's zone is degenerate — which is exactly what a
   zone whose corners are the coplanar points of one face gives you, and it is
   also why the formula cannot use it.
2. `VAR_zcoord` averages the points of every corner with `cmask != 0`.  Adding
   phantom corners to the *real* zone would count the points on its boundary
   face twice and shift its centroid.  Putting them in a separate zone leaves
   real centroids exact.

Phantom sides are required, not optional: `VAR_point_norm` accumulates
boundary point normals by iterating sides with `smask == -1` and following
`m:s>s2` to the corresponding real side.  Without them every boundary point
normalizes a zero vector.

`Boundary_Style::PHANTOM_SIDES_ONLY` is the alternative: no extra zones or
corners, with phantom sides pointing back at the adjacent real zone and its
corners.  It validates too, and produces a mesh with no `mask < 0` zones or
corners at all.  Keeping both makes it possible to reproduce whichever
convention a given input file uses.

One constraint rules out a third, tempting convention — masking a corner as a
boundary corner whenever its point lies on the boundary.  `gradzatp`
accumulates `point_volume` only over `cmask >= 1` corners and then divides
every `pmask != 0` point by it.  Since every corner at a boundary point would
be excluded, that convention makes `point_volume` zero and the gradient `0/0`
at every boundary point.  The audit checks for this directly.

## 6. Validation

Checks are split by what they have to know.

**`Ume::SOA_Idx::audit`** (`src/Ume/mesh_audit.hh`) holds for any mesh the
library will compute on, whatever produced it, so it can be pointed at a mesh
read from a file as easily as at a generated brick:

```c++
Ume::SOA_Idx::Audit_Options opts;
opts.census = true;
Ume::SOA_Idx::Audit_Result const res = Ume::SOA_Idx::audit(mesh, std::cout, opts);
if (!res) { /* res.failed_checks, res.problems */ }
```

`Audit_Options` switches the groups (`topology`, `masks`, `derived`,
`geometry`, `parallel`, `census`) and caps how many offending entities each
check prints.  Together they cover every map the mesh carries: the stored ones
(`m:e>p*`, `m:f>z*`, `m:c>p`, `m:c>z`, `m:s>*`, `m:a>*`), the ones the library
derives from them (`m:c>s`, `m:z>c`, `m:z>p`, `m:z>pz`, `m:p>zs`, `m:p>rc`),
the masks they are read with, and the parallel bookkeeping in each `Entity`.

*Topology* — every map entry in range, including the two stronger conditions
the lazy derivations need: `VAR_zcoord` and `VAR_fcoord` accumulate into
arrays of *local* zones and faces without bounds checking, and the ragged-map
derivations index with `std::vector::at`, which throws rather than ignores a
negative index at a ghost.  Then the relations between the maps:

| what | invariant |
|---|---|
| `m:s>p1`/`p2` | distinct, and the endpoints of the side's own `m:s>e` |
| `m:s>s2` | a self-inverse pairing of sides agreeing on `(face, edge)`, traversing that edge in opposite directions; exactly two sides per `(face, edge)` |
| `m:s>s3`/`s4` | mutually inverse; a single closed directed ring per `(zone, face)`, covering every side of that face of that zone |
| `m:s>s5` | a self-inverse pairing within a zone, on one edge and two faces; a real zone meets each of its edges on exactly two faces |
| `m:s>c1`/`c2` | corners of the side's own zone, at `p1` and `p2` |
| `m:s>f`, `m:f>z*` | a side's zone is one of the two its face names; a face does not name the same zone twice |
| `m:c>p`, `m:c>z` | no two corners of a zone at the same point; a real zone has at least four corners on at least four faces; a ghost corner's zone is out of range or itself a ghost, never real |
| `m:a>*` | each iota copies the zone, face and edge of its `m:a>s`, sits at one of that side's two endpoints, and shares the side with at most one other iota, at the other endpoint |

*Masks* — every `cmask >= 1` corner belongs to a `zmask >= 1` zone, and every
`smask > 0` side likewise; a non-ghost side references nothing that is a ghost;
ghosts are exactly the entities at or above `local_size()`; a side outside the
real mesh has its face, edge and points marked as boundary too, and its
`m:s>s2` is a real side, which is what `VAR_point_norm` follows.  Three checks
name the division that fails without them: every real point owns a
`cmask >= 1` corner (`gradzatp`'s `point_volume`), every non-ghost zone owns a
non-ghost corner (`zcoord`), and every non-ghost face owns a non-ghost side
(`fcoord`).

*Derived* — each derived map against the definition its `VAR_*::init_()` uses:
`m:c>s` inverts `m:s>c1`/`c2`, `m:z>c` inverts `m:c>z`, `m:z>p` and `m:p>zs`
are the two readings of `(m:c>p, m:c>z)`, `m:p>rc` keeps the `cmask >= 1`
corners at each point, and `m:z>pz` is the set of zones sharing a point with a
zone.  These are cached on first use, so this is also what catches a derived
map that was built before a map it came from changed.

*Geometry* — coordinates are finite and no edge has zero length; every real
`side_vol > 0` and every real `corner_vol > 0`; every real zone has positive
volume, and the same volume whether it is added up over its sides or over its
corners; `corner_csurf` closes (sums to zero) over each zone, which is what
makes `gradzatp` reproduce a constant field, and `side_surz` closes over the
surface of each real zone, which is what a face ring traversed the wrong way
about one of its two zones breaks; `point_norm` is a unit vector at
`pmask < 0` points and zero elsewhere.  Per-zone quantities are accumulated
through `m:c>z` rather than `m:z>c`, so the geometry checks do not depend on
how the zone-to-corner map is filtered.

*Parallel* — the `Entity` communication records: matching array lengths, comm
types in range and `GHOST` above `local_size()`, copy records naming a local
index and a rank other than this one, and neighbor lists whose elements are
local indices.  Only one partition is present, so this is what can be checked
without the others.

*Census* (opt-in) — entity counts by mask category, the number of entities
nothing references, the corner counts of the real zones, and two summary
lines:

```
  local zones with zmask < 1: 99 (0 owning cmask >= 1 corners)
  cmask < 0 corners carrying volume: 0
```

Both describe situations where two formulations of the same operator can
disagree.  The first counts zones that a zone-mask filter would skip but that
still own real corners, so an unfiltered loop would write to them.  The second
counts corners that a corner-mask filter would skip but that carry volume, so
including or excluding them changes any volume-weighted zone quantity.  On a
generated brick both are zero: boundary corners are fed only by phantom sides,
and `corner_vol` accumulates only from sides with `smask > 0`.

**`Ume::Brick::validate`** runs the audit and adds what only the generator
knows: every real zone's `corner_vol` sums to `dx*dy*dz` (8 corners) or
`dx*dy*dz/6` (4 corners); the mesh total is the brick volume; `point_norm` is
exactly the outward axis normal in the interior of a brick face; `pmask` agrees
with the lattice.

Two things are worth knowing when the audit is pointed at a partition of a
decomposed mesh rather than at a whole one.  First, it sees only that
partition: ghost entities carry no local connectivity, so a check that would
have to follow a map into a neighboring partition skips them.  Second, where
the file format does not fix a convention — the direction in which `m:s>s2`
partners traverse their shared edge, the direction of the `m:s>s3` ring around
a face — the audit takes the convention from the majority of the mesh's own
sides and reports the ones that disagree, rather than imposing the
generator's.

The audit is tested in both directions.  `[audit]` in `test/test_brick.cc`
checks that generated meshes come out clean, with and without iotas, and that
a broken `s2` link, a mis-pointed `c1`, an out-of-range `m:s>z`, a side whose
points are not those of its edge, a face ring that does not close, a
self-referential `s5`, a face that disowns the zone of its sides, an iota that
disagrees with its side, a stale `m:z>c`, a point stripped of its real
corners, a face left with no live side, a live entity in the ghost index range
and a copy record naming its own rank are each caught with the expected
message.  A map that would index out of bounds is expected to *disable* the
`derived` and `geometry` groups, which would otherwise run into undefined
behavior; the report says which groups it skipped and why.

`ume_serial --audit <file>+` runs the audit on a stored mesh and exits nonzero
if anything fails.  It works under `UME_SERIAL`, which is the configuration in
which the audit's own accumulations are deterministic.

## 7. Gradient accuracy on a mixed mesh

The corner gradient reproduces a linear field `phi = a.x` exactly on a uniform
lattice of one cell type, but *not* at a hex/tet interface — the error is about
4% of `|a|` at the interface points of a 4×4×4 mix, and it does not shrink
under refinement.

This is a property of the discretization rather than of the mesh.  The scheme
replaces `phi(x)` on each median facet with `phi` at the zone center.  Per zone
those errors cancel exactly, because each side contributes `+surf` to `c1` and
`-surf` to `c2` with the same `phi_z`.  Per *point* they cancel only when the
stencil around the point is symmetric enough, which a uniform lattice is and an
interface is not.

Accordingly the linear-field test requires exactness only at points whose
surrounding cells are all of the same type — which is every interior point of a
pure hex or pure tet brick, plus the uniform regions of a mixed one.  A mixed
brick is a good test of topological generality; it is not a good convergence
study.

## 8. Build coverage

The generator itself uses no Kokkos; it is plain C++ and backend agnostic.
The tests reach Kokkos only through the library's derived variables and the
gradient operators, so the one build-dependent thing in `test_brick.cc` is a
tolerance:

```c++
#if defined(UME_SERIAL)
constexpr double GRAD_TOL = 0.0;   // fixed accumulation order => bitwise
#else
constexpr double GRAD_TOL = 1e-12; // Kokkos::atomic_add ordering is not stable
#endif
```

`gradzatz` accumulates with `atomic_add` on every non-`UME_SERIAL` backend —
host OpenMP included, not only GPUs — while `gradzatz_invert` is host-only and
ordered, so an exact comparison between the two is only justified under
`UME_SERIAL`.  Both configurations are exercised:

| build | brick tests | push vs. pull gradient |
|---|---|---|
| default (OpenMP host, multithreaded) | pass | equal to 1e-12 |
| `-DUME_SERIAL=ON` (`Kokkos::Serial`) | pass | bitwise |

`-DUME_SERIAL=ON` needs a Kokkos built with the Serial backend; without one the
library itself fails to compile in `mem_exec_spaces.hh`.

## 9. Limitations

No MPI decomposition.  Every entity is emitted as `INTERNAL` with no ghosts,
neighbor lists are empty, and `Spec` has no partitioning fields, so a generated
mesh is a single-rank mesh.  Adding slab partitioning — ghost layers,
`cpy_idx`/`src_pe`/`src_idx`, and the `myCpys`/`mySrcs` neighbor lists — is the
obvious next step, and would make the generator usable for the parallel
gather/scatter paths as well.

Faces are planar and cells are rectangular, so the generator cannot produce
warped or badly conditioned zones.  Non-uniform `dx/dy/dz` is supported, but
grading and distortion are not.
