# Ume/test/CMakeLists.txt

# test/ links the targets src/ defines and may gate on the variables the
# pkgs/ modules set.  Added any earlier, it sees neither on a fresh configure
# and reads the variables back from the cache on a reconfigure.
if(NOT TARGET Ume)
  message(FATAL_ERROR
    "test/ must be added after src/, from the top-level CMakeLists.txt")
endif()

# UME standard tests

add_executable(ume_tests
  test_comm_buffers.cc
  test_datastore.cc
  test_entity_comm.cc
  test_io.cc
  test_utils.cc
  test_vecn.cc
)

target_link_libraries(ume_tests PRIVATE
  Ume Catch2::Catch2WithMain)

# HAVE_MPI, among others: test_comm_buffers.cc exchanges through Comm::MPI when
# the library was built with it.
target_compile_definitions(ume_tests PRIVATE ${COMMON_COMPILE_DEFINITIONS})

catch_discover_tests(ume_tests)

# UME GPU-enabled tests

# Note: UME_SERIAL is defined in the top-level CMakeLists.txt.

add_executable(ume_gpu_tests
  test_mesh.cc
  test_scratch_arrays.cc
  custom_main.cc
)

target_link_libraries(ume_gpu_tests PRIVATE
  Ume Catch2::Catch2)

catch_discover_tests(ume_gpu_tests)

# UME build-configuration tests

# Every translation unit has to agree on UME_SERIAL; check_ume_serial.cmake
# reads the compilation database, which only the Ninja and Makefile generators
# write.
if(CMAKE_EXPORT_COMPILE_COMMANDS AND CMAKE_GENERATOR MATCHES "Ninja|Makefiles")
  add_test(NAME cmake_ume_serial_consistent
    COMMAND ${CMAKE_COMMAND}
      -DCOMPILE_COMMANDS=${CMAKE_BINARY_DIR}/compile_commands.json
      -DUME_SERIAL=${UME_SERIAL}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/check_ume_serial.cmake
    )
endif()

# Driver test meshes

# Writes a one-partition mesh holding only zones with the given masks; the
# scale_mesh and ume_mpi driver tests below read it.
add_executable(write_zone_mask_mesh write_zone_mask_mesh.cc)
target_link_libraries(write_zone_mask_mesh PRIVATE Ume)
target_compile_definitions(write_zone_mask_mesh PRIVATE
  ${COMMON_COMPILE_DEFINITIONS})

# scale_mesh driver tests

# A basename with no partition file behind it has to stop scale_mesh before it
# reads the unopened stream as a mesh.  "Aborting." is printed only on the way
# to main's `return 1`; an unchecked open went on to Mesh::read, which in a
# Debug build stops at an assertion instead.
add_test(NAME scale_mesh_unreadable_mesh
  COMMAND scale_mesh ${CMAKE_CURRENT_BINARY_DIR}/no-such-mesh 2)
set_tests_properties(scale_mesh_unreadable_mesh PROPERTIES
  PASS_REGULAR_EXPRESSION
    "Unable to open file \"[^\"]*/no-such-mesh\\.00000\\.ume\" for reading\\.[\r\n]+Aborting\\."
  TIMEOUT 90)

# An output file that cannot be opened has to fail the run too: an unchecked
# open wrote nothing, printed "Done." and exited 0.  A directory sits where
# scale 1 writes its partition file.
set(_scale_mesh_dir ${CMAKE_CURRENT_BINARY_DIR}/scale_mesh_meshes)
file(MAKE_DIRECTORY ${_scale_mesh_dir}/unwritable.00001.00000.ume)
add_test(NAME scale_mesh_write_input_mesh
  COMMAND write_zone_mask_mesh ${_scale_mesh_dir}/unwritable.00000.ume 1)
set_tests_properties(scale_mesh_write_input_mesh PROPERTIES
  FIXTURES_SETUP scale_mesh_input_mesh)
add_test(NAME scale_mesh_unwritable_output
  COMMAND scale_mesh ${_scale_mesh_dir}/unwritable 1)
set_tests_properties(scale_mesh_unwritable_output PROPERTIES
  FIXTURES_REQUIRED scale_mesh_input_mesh
  PASS_REGULAR_EXPRESSION
    "Unable to open file \"[^\"]*/unwritable\\.00001\\.00000\\.ume\" for writing\\.[\r\n]+Aborting\\."
  TIMEOUT 90)

# UME driver startup tests

# ume_mpi run as an MPI singleton: a missing basename, one that the
# --kokkos-* argument stripping leaves missing, and a partition file that
# cannot be read all have to end in Comm::MPI::abort.  One rank cannot show
# the others being taken down rather than left in a gathscat; the abort
# message is what tells abort() from a return out of main().
if(TARGET ume_mpi)
  set(_ume_mpi_check
    -DUME_MPI=$<TARGET_FILE:ume_mpi>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_ume_mpi_startup.cmake)
  add_test(NAME ume_mpi_no_basename
    COMMAND ${CMAKE_COMMAND}
      "-DEXPECT=Usage: ume_mpi <basename>.*Transport::abort: no mesh basename"
      ${_ume_mpi_check})
  add_test(NAME ume_mpi_only_kokkos_args
    COMMAND ${CMAKE_COMMAND}
      -DARGS=--kokkos-num-threads=1
      "-DEXPECT=Usage: ume_mpi <basename>.*Transport::abort: no mesh basename"
      ${_ume_mpi_check})
  add_test(NAME ume_mpi_unreadable_mesh
    COMMAND ${CMAKE_COMMAND}
      -DARGS=${CMAKE_CURRENT_BINARY_DIR}/no-such-mesh
      "-DEXPECT=Transport::abort: this rank could not read its mesh file"
      ${_ume_mpi_check})
  set_tests_properties(
    ume_mpi_no_basename ume_mpi_only_kokkos_args ume_mpi_unreadable_mesh
    PROPERTIES TIMEOUT 90)

  # ume_mpi seeds its gradient in one interior zone (mask >= 1) of each rank.
  # Searching upward from the middle of the zone range ran off the end when the
  # upper half held none; a rank with no interior zone at all, which an
  # over-decomposed mesh has, has to take the job down.
  set(_zone_mask_dir ${CMAKE_CURRENT_BINARY_DIR}/zone_mask_meshes)
  file(MAKE_DIRECTORY ${_zone_mask_dir})
  # ume_mpi_long_basename names its mesh by an absolute basename more than
  # 112 characters long, whatever the length of the build directory's path:
  # read_mesh() once formatted the file name into a char[80].
  string(REPEAT "/long-basename" 8 _long_subdir)
  set(_long_basename_dir ${_zone_mask_dir}${_long_subdir})
  file(MAKE_DIRECTORY ${_long_basename_dir})
  add_test(NAME ume_mpi_write_no_zone_mesh
    COMMAND write_zone_mask_mesh no-zones.00000.ume
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_write_lower_interior_mesh
    COMMAND write_zone_mask_mesh lower-interior.00000.ume 1 1 -1 -1
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_write_long_basename_mesh
    COMMAND write_zone_mask_mesh
      ${_long_basename_dir}/lower-interior.00000.ume 1 1 -1 -1)
  set_tests_properties(
    ume_mpi_write_no_zone_mesh ume_mpi_write_lower_interior_mesh
    ume_mpi_write_long_basename_mesh PROPERTIES
    FIXTURES_SETUP ume_mpi_zone_mask_meshes)

  add_test(NAME ume_mpi_no_interior_zone
    COMMAND ${CMAKE_COMMAND}
      -DARGS=no-zones
      "-DEXPECT=Transport::abort: .*ume_mpi cannot run an over-decomposed mesh"
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_interior_zones_in_lower_half
    COMMAND ${CMAKE_COMMAND}
      -DARGS=lower-interior
      -DEXPECT_SUCCESS=ON
      "-DEXPECT=Checking gradient result.*Done\\."
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_long_basename
    COMMAND ${CMAKE_COMMAND}
      -DARGS=${_long_basename_dir}/lower-interior
      -DEXPECT_SUCCESS=ON
      "-DEXPECT=Checking gradient result.*Done\\."
      ${_ume_mpi_check})
  set_tests_properties(
    ume_mpi_no_interior_zone ume_mpi_interior_zones_in_lower_half
    ume_mpi_long_basename PROPERTIES
    FIXTURES_REQUIRED ume_mpi_zone_mask_meshes
    TIMEOUT 90)
endif()
