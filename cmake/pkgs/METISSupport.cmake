
# Optional support for the METIS graph partitioner.
#
# METIS is not required, and a build without it is an ordinary configuration:
# the `metis-*` partitioning methods simply cannot run.  They stay in the
# method table either way -- see `Partition.hh` -- so `--method metis-rb` on a
# build configured without METIS reports that it needs METIS rather than that
# it is a misspelling.
#
# METIS 5.1+ (KarypisLab) is Apache-2.0, so depending on it optionally is
# fine.  ParMETIS is *not* freely redistributable for commercial use and must
# not be added.

include(pkgs/helpers)
SimplePkg(METIS metis.h metis)
if (METIS_NOT_FOUND)
  set(COMMON_COMPILE_DEFINITIONS
    ${COMMON_COMPILE_DEFINITIONS}
    NO_METIS
    )
  set(METIS_INCLUDE_DIRS "")
  message(STATUS "Partitioning: no METIS, so the metis-* methods are absent")
else()
  list(APPEND COMMON_LINK_LIBRARIES ${METIS_LIBRARIES})
  message(STATUS "Partitioning: METIS is available (${METIS_LIBRARIES})")
endif()
