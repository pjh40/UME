# A Spack repo for UME's optional dependencies

At the moment this holds one package, `mt-kahypar`, because Mt-KaHyPar is not
in Spack's builtin repo and is not in Homebrew either — and UME's
`kahypar`/`kahypar-quality` partitioning methods need it. METIS, by contrast,
is builtin (`spack install metis`) and needs nothing from here.

```sh
spack repo add --scope=user scripts/spack-repo/spack_repo/ume
spack install mt-kahypar
```

The path ends in `spack_repo/ume` and not at this directory: Spack's repo API
v2 requires `spack_repo/<namespace>/` in the path, and rejects anything else.
`--scope=user` puts it where environments inherit it, rather than inside one
environment's `spack.yaml`.

`cmake/pkgs/KaHyParSupport.cmake` then finds it through the environment view
with no `-D` flags, much as `METISSupport.cmake` finds METIS — though through
`find_package(MtKaHyPar CONFIG)` rather than `SimplePkg`, because upstream
ships a config package and a bare header/lib probe would not carry the
transitive TBB link. A build without it is an ordinary configuration: the
`kahypar*` methods report themselves as unavailable and nothing else changes.

### On macOS you will need two more things

`spack compiler find` first, if only GCC is registered. The recipe refuses to
build with g++ on darwin (see below), so with no clang registered the
concretizer has nothing to choose and fails with that message.

And if you want it *in* an environment that pins its compiler — as UME's own
`ume` environment does, with `all: require: ...c=gcc-16` — mt-kahypar needs an
exception, since it cannot satisfy that on darwin:

```yaml
packages:
  mt-kahypar:
    require: '%[when=%c]c=apple-clang %[when=%cxx]cxx=apple-clang'
```

UME itself stays on GCC. The boundary is a C API and that combination is
tested: a binary linking UME, Kokkos, METIS and a clang-built libmtkahypar,
with libstdc++ and libc++ both loaded, runs correctly.

## Why the recipe looks the way it does

Each of these is a measured obstacle rather than a precaution; the reasoning
is in `doc/partitioning-plan.md` under "Phase 4, Mt-KaHyPar".

- **`KAHYPAR_DOWNLOAD_TBB` is off.** Upstream can fetch and build oneTBB
  itself, but that copy is `EXCLUDE_FROM_ALL`, so `install-mtkahypar` does not
  install it and the installed library is left with a dangling
  `@rpath/libtbb.12.dylib`. The dependency has to be a real package.
- **`%gcc` conflicts on macOS.** oneTBB's macOS export list names symbols as
  libc++ mangles them, so a libstdc++ build fails to link with an undefined
  vtable for `tbb::detail::r1::unsafe_wait`. Build it `%apple-clang` there.
  UME itself can stay on GCC — the boundary is a C API, and that combination
  is tested. On Linux the question does not arise.
- **`-march=native` is off.** Spack's compiler wrapper supplies the target
  flags; letting the package add its own would make the installed library
  depend on the machine that built it rather than on the spec.
- **The x86 instruction flags are gated on the target family.** They are
  `-mcx16 -msse4.2 -mcrc32`, which are fatal to name on aarch64.
- **`install-mtkahypar`, not `install`.** The plain target also wants the CLI
  and the Python module; the library interface is all UME links against.

## This is a stopgap

The right home for this recipe is Spack's builtin repo. If `mt-kahypar`
appears there, delete this directory and drop the `spack repo add` step —
nothing in UME's CMake refers to it.
