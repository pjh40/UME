# Copyright (c) 2023, Triad National Security, LLC. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class MtKahypar(CMakePackage):
    """Mt-KaHyPar is a shared-memory multilevel graph and hypergraph
    partitioner.

    UME uses its connectivity (lambda-1) objective on the hypergraph whose
    vertices are zones and whose nets are points: a point held by lambda ranks
    needs lambda-1 copies, so minimizing the sum of (lambda-1) minimizes the
    doubles moved by `points.gathscat` exactly, rather than approximating it
    the way a graph edge cut does.
    """

    homepage = "https://github.com/kahypar/mt-kahypar"
    git = "https://github.com/kahypar/mt-kahypar.git"

    license("MIT")

    # Pinned by commit rather than by a release tarball's sha256: upstream's
    # CMake expects a git checkout, and a commit is checkable without
    # fetching an archive first.
    version("1.6.2", tag="v1.6.2", commit="5495f2a4768ef070f2cac1fea298be53a8de85d6")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.26:", type="build")

    depends_on("tbb")
    # Upstream requires oneTBB 2021.5 or newer.
    conflicts("^intel-tbb@:2021.4")

    variant("shared", default=True, description="Build a shared library")
    variant("hwloc", default=True, description="Use hwloc for NUMA awareness")
    variant("graph", default=False, description="Graph partitioning features")
    variant(
        "highest_quality",
        default=True,
        description="n-level coarsening, i.e. the HIGHEST_QUALITY preset",
    )
    variant("large_k", default=False, description="Features for k >= 1024")
    variant("soed", default=False, description="Sum-of-external-degree metric")

    depends_on("hwloc", when="+hwloc")

    # oneTBB's macOS export list (src/tbb/def/mac64-tbb.def) names symbols as
    # libc++ mangles them, so a libstdc++ build fails to link with an
    # undefined vtable for tbb::detail::r1::unsafe_wait.  Its CMake also reads
    # a "GNU assembler version" out of `-Wa,-v`, which Apple's assembler does
    # not print.  Neither is Mt-KaHyPar's doing and neither applies on Linux.
    #
    # `%[virtuals=cxx] gcc` and not a bare `%gcc`: the latter matches gcc
    # anywhere in the DAG, so a dependency's gcc-runtime would trip it and no
    # spec at all would concretize on darwin.  What is wrong here is building
    # *this* with g++.
    conflicts(
        "%[virtuals=cxx] gcc",
        when="platform=darwin",
        msg="oneTBB does not support GCC on macOS; build this with apple-clang",
    )


class CMakeBuilder(cmake.CMakeBuilder):
    # `install` would also want the CLI and the Python module; the library
    # interface has an install target of its own, and it is the only thing
    # UME links against.
    build_targets = ["mtkahypar"]
    install_targets = ["install-mtkahypar"]

    def cmake_args(self):
        spec = self.spec
        return [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES", "graph"),
            self.define_from_variant(
                "KAHYPAR_ENABLE_HIGHEST_QUALITY_FEATURES", "highest_quality"
            ),
            self.define_from_variant("KAHYPAR_ENABLE_LARGE_K_PARTITIONING_FEATURES", "large_k"),
            self.define_from_variant("KAHYPAR_ENABLE_SOED_METRIC", "soed"),
            self.define("KAHYPAR_ENABLE_TESTING", False),
            # Take TBB from the dependency rather than fetching it: the
            # downloaded copy is EXCLUDE_FROM_ALL and so is never installed,
            # which leaves the installed library with a dangling
            # @rpath/libtbb.12.dylib.
            self.define("KAHYPAR_DOWNLOAD_TBB", False),
            # Spack's compiler wrapper already supplies the target flags, and
            # -march=native would make the build depend on the machine that
            # ran it rather than on the spec.
            self.define("KAHYPAR_ENABLE_ARCH_COMPILE_OPTIMIZATIONS", False),
            self.define("KAHYPAR_DISABLE_HWLOC", spec.satisfies("~hwloc")),
            # -mcx16 -msse4.2 -mcrc32: x86 only, and fatal to name elsewhere.
            self.define(
                "KAHYPAR_ENABLE_EXTENDED_INSTRUCTIONS", spec.target.family == "x86_64"
            ),
            # Upstream makes the Steiner tree metric without graph features a
            # configure error rather than ignoring it.
            self.define("KAHYPAR_ENABLE_STEINER_TREE_METRIC", spec.satisfies("+graph")),
        ]
