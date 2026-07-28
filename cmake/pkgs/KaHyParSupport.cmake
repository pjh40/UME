
# Optional support for the Mt-KaHyPar hypergraph partitioner.
#
# Mt-KaHyPar is not required, and a build without it is an ordinary
# configuration: the `kahypar*` partitioning methods simply cannot run.  They
# stay in the method table either way -- see `Partition.hh` -- so
# `--method kahypar` on a build configured without it reports that it needs
# Mt-KaHyPar rather than that it is a misspelling.
#
# Mt-KaHyPar is MIT-licensed, so depending on it optionally is fine.  The
# original KaHyPar is GPL-3 and must not be added; UME is BSD-3.
#
# This one uses `find_package(... CONFIG)` rather than `SimplePkg`, which the
# other optional libraries here use, for two reasons: upstream ships a config
# package, and `libmtkahypar` links oneTBB, which a bare header-and-library
# probe would not carry over.  It is not built by a `FetchContent` for the
# reason the rest of this project is arranged the way it is -- which options a
# tree is built with is a property of the machine, not of the project -- and
# `scripts/spack-repo/` carries a Spack recipe for installing it.
#
# The config package does not set INTERFACE_INCLUDE_DIRECTORIES on its
# imported target, so the header directory is found separately; a prefix with
# the library but not the header is a broken install and is treated as an
# absent one.

find_package(MtKaHyPar CONFIG QUIET)

if (MtKaHyPar_FOUND)
  get_filename_component(KAHYPAR_PREFIX "${MtKaHyPar_DIR}/../../.." ABSOLUTE)
  find_path(KAHYPAR_INCLUDE_DIRS mtkahypar.h
    HINTS "${KAHYPAR_PREFIX}/include")
endif()

if (MtKaHyPar_FOUND AND KAHYPAR_INCLUDE_DIRS)
  set(KAHYPAR_NOT_FOUND FALSE)
  list(APPEND COMMON_LINK_LIBRARIES MtKaHyPar::mtkahypar)
  message(STATUS
    "Partitioning: Mt-KaHyPar ${MtKaHyPar_VERSION} is available (${MtKaHyPar_DIR})")
else()
  set(KAHYPAR_NOT_FOUND TRUE)
  set(COMMON_COMPILE_DEFINITIONS
    ${COMMON_COMPILE_DEFINITIONS}
    NO_KAHYPAR
    )
  set(KAHYPAR_INCLUDE_DIRS "")
  message(STATUS "Partitioning: no Mt-KaHyPar, so the kahypar* methods are absent")
endif()
