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

# error_stop tests

# Pool exhaustion has to print its message.  error_stop once finalized first,
# and the pool's Finalize asserts that no claims are outstanding, which the
# exhaustion leaves behind: a Debug build died on the assertion and the
# message that explained it never printed.  With a 2 MB pool the same probe's
# claims fit, and the message must not appear.
add_executable(error_stop_probe error_stop_probe.cc)
target_link_libraries(error_stop_probe PRIVATE Ume)
target_compile_definitions(error_stop_probe PRIVATE
  ${COMMON_COMPILE_DEFINITIONS})
set(_error_stop_check
  -DPROBE=$<TARGET_FILE:error_stop_probe>
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_error_stop.cmake)
add_test(NAME error_stop_pool_exhausted
  COMMAND ${CMAKE_COMMAND}
    -DPOOL_MB=1
    "-DEXPECT=Memory pool cannot claim memory\\."
    "-DREJECT=both claims fit"
    ${_error_stop_check})
add_test(NAME error_stop_pool_fits
  COMMAND ${CMAKE_COMMAND}
    -DPOOL_MB=2
    -DEXPECT_SUCCESS=ON
    "-DEXPECT=error_stop_probe: both claims fit"
    "-DREJECT=Memory pool cannot claim memory"
    ${_error_stop_check})
set_tests_properties(error_stop_pool_exhausted error_stop_pool_fits
  PROPERTIES TIMEOUT 90)

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

# A partition file that ends early has to fail the run too: Mesh::read goes
# quiet once its stream fails, and scale_mesh scaled and wrote what it had and
# printed "Done.".  Ten bytes off the end take the zones' subset count and the
# newlines after it; the cut has to fall after the last entity and neighbor
# tag, which a Debug build asserts on.  The same mesh written whole reaches
# "Done.", the control.
add_test(NAME scale_mesh_write_truncated_mesh
  COMMAND write_zone_mask_mesh --drop-last=10
    ${_scale_mesh_dir}/truncated.00000.ume 1 1 -1 -1)
add_test(NAME scale_mesh_write_intact_mesh
  COMMAND write_zone_mask_mesh ${_scale_mesh_dir}/intact.00000.ume 1 1 -1 -1)
set_tests_properties(scale_mesh_write_truncated_mesh PROPERTIES
  FIXTURES_SETUP scale_mesh_truncated_input)
set_tests_properties(scale_mesh_write_intact_mesh PROPERTIES
  FIXTURES_SETUP scale_mesh_intact_input)
add_test(NAME scale_mesh_truncated_mesh
  COMMAND scale_mesh ${_scale_mesh_dir}/truncated 1)
set_tests_properties(scale_mesh_truncated_mesh PROPERTIES
  FIXTURES_REQUIRED scale_mesh_truncated_input
  PASS_REGULAR_EXPRESSION
    "Unable to read a mesh from file \"[^\"]*/truncated\\.00000\\.ume\"\\.[\r\n]+Aborting\\."
  TIMEOUT 90)
add_test(NAME scale_mesh_intact_mesh
  COMMAND scale_mesh ${_scale_mesh_dir}/intact 1)
set_tests_properties(scale_mesh_intact_mesh PROPERTIES
  FIXTURES_REQUIRED scale_mesh_intact_input
  PASS_REGULAR_EXPRESSION "Done\\."
  FAIL_REGULAR_EXPRESSION "Aborting\\."
  TIMEOUT 90)

# An output file that opens but cannot be written, as on a full disk, has to
# fail the run as well: an unchecked write left an empty file, printed "Done."
# and exited 0.  A file-size limit of zero stands in for the full disk, with
# SIGXFSZ ignored so that the write fails with EFBIG instead of killing the
# process; the limit does not apply to the read.  scale_mesh_intact_mesh, the
# same mesh with no limit, is the control.
if(UNIX)
  add_test(NAME scale_mesh_write_size_limited_mesh
    COMMAND write_zone_mask_mesh
      ${_scale_mesh_dir}/size_limited.00000.ume 1 1 -1 -1)
  set_tests_properties(scale_mesh_write_size_limited_mesh PROPERTIES
    FIXTURES_SETUP scale_mesh_size_limited_input)
  add_test(NAME scale_mesh_output_over_size_limit
    COMMAND sh -c "trap '' XFSZ; ulimit -f 0; exec \"$0\" \"$@\""
      $<TARGET_FILE:scale_mesh> ${_scale_mesh_dir}/size_limited 1)
  set_tests_properties(scale_mesh_output_over_size_limit PROPERTIES
    FIXTURES_REQUIRED scale_mesh_size_limited_input
    PASS_REGULAR_EXPRESSION
      "Unable to write a mesh to file \"[^\"]*/size_limited\\.00001\\.00000\\.ume\"\\.[\r\n]+Aborting\\."
    TIMEOUT 90)
endif()

# Writes a one-partition mesh of points, or prints the points of one; the
# stitch tests below use it on both sides of scale_mesh.
add_executable(point_mesh point_mesh.cc)
target_link_libraries(point_mesh PRIVATE Ume)
target_compile_definitions(point_mesh PRIVATE ${COMMON_COMPILE_DEFINITIONS})

# Scale 2 stitches one copy of the mesh to its +x side, translated by the
# mesh's x extent.  The extent was taken as the largest x, which is right
# only when the smallest is 0: on [-1, 0] the copy landed on its original, and
# on [1, 3] it left a gap of 1.  Each test writes the mesh on the given
# coordinates, scales it, and matches every point of the result.  A copy
# keeps its original's y and z: scale_mesh once set only the stitched
# component of a copied point and left the other two at whatever the resize
# put there, which a mesh with every point at y = z = 0 could not show and
# the off_axis mesh does.
function(add_scale_mesh_stitch_test name expected)
  set(_basename ${_scale_mesh_dir}/${name})
  add_test(NAME scale_mesh_write_${name}_mesh
    COMMAND point_mesh write ${_basename}.00000.ume ${ARGN})
  set_tests_properties(scale_mesh_write_${name}_mesh PROPERTIES
    FIXTURES_SETUP scale_mesh_${name}_input)
  add_test(NAME scale_mesh_stitch_${name}
    COMMAND scale_mesh ${_basename} 2)
  set_tests_properties(scale_mesh_stitch_${name} PROPERTIES
    FIXTURES_REQUIRED scale_mesh_${name}_input
    FIXTURES_SETUP scale_mesh_${name}_output
    TIMEOUT 90)
  add_test(NAME scale_mesh_stitched_${name}_points
    COMMAND point_mesh print ${_basename}.00002.00000.ume)
  set_tests_properties(scale_mesh_stitched_${name}_points PROPERTIES
    FIXTURES_REQUIRED scale_mesh_${name}_output
    PASS_REGULAR_EXPRESSION "^${expected}$")
endfunction()

set(_p "[\r\n]+point")
add_scale_mesh_stitch_test(negative
  "point 0: -1 0 0${_p} 1: 0 0 0${_p} 2: 0 0 0${_p} 3: 1 0 0[\r\n]*"
  -1 0)
add_scale_mesh_stitch_test(positive
  "point 0: 1 0 0${_p} 1: 3 0 0${_p} 2: 3 0 0${_p} 3: 5 0 0[\r\n]*"
  1 3)
add_scale_mesh_stitch_test(off_axis
  "point 0: 0 1 2${_p} 1: 1 -3 0.5${_p} 2: 1 1 2${_p} 3: 2 -3 0.5[\r\n]*"
  0,1,2 1,-3,0.5)

# Writes a one-partition mesh with every entity, or prints the entity data of
# one; the tests below use it on both sides of scale_mesh.
add_executable(entity_mesh entity_mesh.cc)
target_link_libraries(entity_mesh PRIVATE Ume)
target_compile_definitions(entity_mesh PRIVATE ${COMMON_COMPILE_DEFINITIONS})

# The stitched copy of each entity has to carry its original's mask,
# communication type and subsets, and its index maps have to keep the -1
# sentinels.  scale_mesh copied only the index maps and the coordinates: the
# copied half had mask 0, which every kernel skips, and no subsets, and the
# offset turned each sentinel into the original's last entity: a boundary
# face's missing second zone, `m:f>z2 == -1`, became a zone of the original.
# The expected lines are entity_mesh's mesh doubled along x by hand.
set(_entity_lines
  "points lsize 4 size 4"
  "points mask: 1 -1 1 -1"
  "points comm_type: 1 2 1 2"
  "points subset corner lsize 2 elements: 1 3 mask: 3 3"
  "edges lsize 2 size 2"
  "edges mask: 1 1"
  "edges comm_type: 1 1"
  "faces lsize 6 size 6"
  "faces mask: 1 -1 1 1 -1 1"
  "faces comm_type: 1 2 3 1 2 3"
  "faces subset boundary lsize 4 elements: 1 2 4 5 mask: -1 -2 -1 -2"
  "sides lsize 4 size 4"
  "sides mask: 1 -1 1 -1"
  "sides comm_type: 1 2 1 2"
  "corners lsize 4 size 4"
  "corners mask: 1 -1 1 -1"
  "corners comm_type: 1 2 1 2"
  "zones lsize 4 size 4"
  "zones mask: 1 -1 1 -1"
  "zones comm_type: 1 2 1 2"
  "zones subset inner lsize 2 elements: 0 2 mask: 1 1"
  "zones subset all lsize 4 elements: 0 1 2 3 mask: 1 -1 1 -1"
  "iotas lsize 2 size 2"
  "iotas mask: 1 1"
  "iotas comm_type: 1 1"
  "m:e>p1: 0 2"
  "m:e>p2: 1 3"
  "m:f>z1: 0 0 1 2 2 3"
  "m:f>z2: 1 -1 -1 3 -1 -1"
  "m:s>z: 0 1 2 3"
  "m:s>p1: 0 1 2 3"
  "m:s>p2: 1 0 3 2"
  "m:s>e: 0 0 1 1"
  "m:s>f: 0 2 3 5"
  "m:s>c1: 0 1 2 3"
  "m:s>c2: 1 0 3 2"
  "m:s>s2: 1 0 3 2"
  "m:s>s3: 0 1 2 3"
  "m:s>s4: 1 0 3 2"
  "m:s>s5: 1 -1 3 -1"
  "m:c>p: 0 1 2 3"
  "m:c>z: 0 1 2 3"
  "m:a>z: 1 3"
  "m:a>f: 2 5"
  "m:a>p: 1 3"
  "m:a>e: 0 1"
  "m:a>s: 1 3")
string(JOIN "[\r\n]+" _entity_expected ${_entity_lines})
add_test(NAME scale_mesh_write_entity_mesh
  COMMAND entity_mesh write ${_scale_mesh_dir}/entities.00000.ume)
set_tests_properties(scale_mesh_write_entity_mesh PROPERTIES
  FIXTURES_SETUP scale_mesh_entity_input)
add_test(NAME scale_mesh_stitch_entities
  COMMAND scale_mesh ${_scale_mesh_dir}/entities 2)
set_tests_properties(scale_mesh_stitch_entities PROPERTIES
  FIXTURES_REQUIRED scale_mesh_entity_input
  FIXTURES_SETUP scale_mesh_entity_output
  PASS_REGULAR_EXPRESSION "Done\\."
  FAIL_REGULAR_EXPRESSION "Aborting\\."
  TIMEOUT 90)
add_test(NAME scale_mesh_stitched_entities
  COMMAND entity_mesh print ${_scale_mesh_dir}/entities.00002.00000.ume)
set_tests_properties(scale_mesh_stitched_entities PROPERTIES
  FIXTURES_REQUIRED scale_mesh_entity_output
  PASS_REGULAR_EXPRESSION "^${_entity_expected}[\r\n]*$")

# scale_mesh stitches a single-rank mesh without ghosts and has to refuse any
# other: a partition of a larger mesh, or a mesh with ghosts, whose copies
# would land among the local entities.  It went on and wrote a mesh whose
# copied half had no exchange lists.  scale_mesh_stitch_entities, the same
# mesh with one rank and no ghosts, is the control.
foreach(_refused IN ITEMS numpe ghost)
  if(_refused STREQUAL "numpe")
    set(_flag --numpe=2)
    set(_message "scale_mesh needs a single-rank mesh\\; this is rank 0 of 2")
  else()
    set(_flag --ghost-zone)
    set(_message "scale_mesh needs a mesh without ghosts\\; its zones \\[2, 3\\) are ghosts")
  endif()
  add_test(NAME scale_mesh_write_${_refused}_mesh
    COMMAND entity_mesh write ${_scale_mesh_dir}/${_refused}.00000.ume ${_flag})
  set_tests_properties(scale_mesh_write_${_refused}_mesh PROPERTIES
    FIXTURES_SETUP scale_mesh_${_refused}_input)
  add_test(NAME scale_mesh_refuses_${_refused}
    COMMAND scale_mesh ${_scale_mesh_dir}/${_refused} 2)
  set_tests_properties(scale_mesh_refuses_${_refused} PROPERTIES
    FIXTURES_REQUIRED scale_mesh_${_refused}_input
    PASS_REGULAR_EXPRESSION "${_message}\\.[\r\n]+Aborting\\."
    FAIL_REGULAR_EXPRESSION "Done\\."
    TIMEOUT 90)
endforeach()

# The scale factor has to be a whole positive power of 2.  atoi read "1x" as
# 1 and scaled by it; it read the most negative int as itself, which the
# power-of-2 test `scale & (scale - 1)` passed (by signed overflow), and the
# run wrote an unscaled mesh under a negative scale.  scale_mesh_intact_mesh,
# scale 1 on the same mesh, is the control.
foreach(_bad_scale IN ITEMS 1x -2147483648)
  add_test(NAME scale_mesh_bad_scale_${_bad_scale}
    COMMAND scale_mesh ${_scale_mesh_dir}/intact ${_bad_scale})
  set_tests_properties(scale_mesh_bad_scale_${_bad_scale} PROPERTIES
    FIXTURES_REQUIRED scale_mesh_intact_input
    PASS_REGULAR_EXPRESSION "Scale must be a power of 2, got \"${_bad_scale}\""
    FAIL_REGULAR_EXPRESSION "Done\\."
    TIMEOUT 90)
endforeach()

# Runs a driver through check_exit_status.cmake, which holds it to an exit
# status as well as to its output.
set(_exit_status_check
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_exit_status.cmake)

# ume_serial driver tests

# A rank file that ends early has to fail the run as one that cannot be
# opened does: ume_serial checked only the open, read the file as far as it
# went and computed on it, and exited 0.  The cut is scale_mesh_truncated_mesh's.
# The same mesh whole exits 0, the control.
set(_ume_serial_dir ${CMAKE_CURRENT_BINARY_DIR}/ume_serial_meshes)
file(MAKE_DIRECTORY ${_ume_serial_dir})
add_test(NAME ume_serial_write_truncated_mesh
  COMMAND write_zone_mask_mesh --drop-last=10
    ${_ume_serial_dir}/truncated.00000.ume 1 1 -1 -1)
add_test(NAME ume_serial_write_intact_mesh
  COMMAND write_zone_mask_mesh ${_ume_serial_dir}/intact.00000.ume 1 1 -1 -1)
set_tests_properties(
  ume_serial_write_truncated_mesh ume_serial_write_intact_mesh PROPERTIES
  FIXTURES_SETUP ume_serial_meshes)
add_test(NAME ume_serial_truncated_mesh
  COMMAND ${CMAKE_COMMAND}
    -DPROGRAM=$<TARGET_FILE:ume_serial>
    -DARGS=${_ume_serial_dir}/truncated.00000.ume
    -DSTATUS=1
    "-DEXPECT=Unable to read a mesh from file \"[^\"]*/truncated\\.00000\\.ume\"\\."
    ${_exit_status_check})
add_test(NAME ume_serial_intact_mesh
  COMMAND ${CMAKE_COMMAND}
    -DPROGRAM=$<TARGET_FILE:ume_serial>
    -DARGS=${_ume_serial_dir}/intact.00000.ume
    -DSTATUS=0
    "-DEXPECT=Reading: [^\r\n]*/intact\\.00000\\.ume"
    "-DREJECT=Unable to"
    ${_exit_status_check})
set_tests_properties(ume_serial_truncated_mesh ume_serial_intact_mesh
  PROPERTIES
  FIXTURES_REQUIRED ume_serial_meshes
  TIMEOUT 90)

# txt2bin driver tests

# txt2bin_point_mesh.txt is a one-point mesh in the text format.  A binary
# file that opens but cannot be written has to fail the conversion: an
# unchecked write left an empty file, which the verification read then took
# for a mesh (and a Debug build stopped at an assertion in it).  The zero
# file-size limit, with SIGXFSZ ignored, is scale_mesh_output_over_size_limit's.
# The check runs sh, which sets both and execs txt2bin: set on the cmake that
# runs the check instead, the limit would reach txt2bin but the ignored
# SIGXFSZ would not, since execute_process resets the child's signals.
# txt2bin_point_mesh, the same conversion with no limit, is the control.
set(_txt2bin_dir ${CMAKE_CURRENT_BINARY_DIR}/txt2bin_meshes)
file(MAKE_DIRECTORY ${_txt2bin_dir})
set(_txt2bin_input ${CMAKE_CURRENT_SOURCE_DIR}/txt2bin_point_mesh.txt)
add_test(NAME txt2bin_point_mesh
  COMMAND ${CMAKE_COMMAND}
    -DPROGRAM=$<TARGET_FILE:txt2bin>
    "-DARGS=${_txt2bin_input}$<SEMICOLON>${_txt2bin_dir}/point.ume"
    -DSTATUS=0
    "-DEXPECT=Copy verified"
    "-DREJECT=Unable to"
    ${_exit_status_check})
set_tests_properties(txt2bin_point_mesh PROPERTIES TIMEOUT 90)
if(UNIX)
  set(_size_limited_args
    -c "trap '' XFSZ && ulimit -f 0 && exec \"$0\" \"$@\""
    $<TARGET_FILE:txt2bin> ${_txt2bin_input}
    ${_txt2bin_dir}/size_limited.ume)
  list(JOIN _size_limited_args "$<SEMICOLON>" _size_limited_args)
  add_test(NAME txt2bin_output_over_size_limit
    COMMAND ${CMAKE_COMMAND}
      -DPROGRAM=sh
      "-DARGS=${_size_limited_args}"
      -DSTATUS=4
      "-DEXPECT=Unable to write a mesh to file \"[^\"]*/size_limited\\.ume\"\\."
      ${_exit_status_check})
  set_tests_properties(txt2bin_output_over_size_limit PROPERTIES TIMEOUT 90)
endif()

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
  # The lower-interior mesh cut short as in scale_mesh_truncated_mesh above.
  add_test(NAME ume_mpi_write_truncated_mesh
    COMMAND write_zone_mask_mesh --drop-last=10
      truncated.00000.ume 1 1 -1 -1
    WORKING_DIRECTORY ${_zone_mask_dir})
  set_tests_properties(
    ume_mpi_write_no_zone_mesh ume_mpi_write_lower_interior_mesh
    ume_mpi_write_long_basename_mesh ume_mpi_write_truncated_mesh PROPERTIES
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
  # A partition file that ends early has to take the job down as one that
  # cannot be opened does; unchecked, the rank ran on whatever Mesh::read had
  # reached.  ume_mpi_interior_zones_in_lower_half, on the same mesh whole, is
  # the control.
  add_test(NAME ume_mpi_truncated_mesh
    COMMAND ${CMAKE_COMMAND}
      -DARGS=truncated
      "-DEXPECT=Unable to read a mesh from file \"truncated\\.00000\\.ume\"\\..*Transport::abort: this rank could not read its mesh file"
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  # The count after -i has to be a non-negative integer.  atoi read "-1" as -1,
  # which became SIZE_MAX iterations of every timed loop, and "2x" as 2; both
  # have to take the job down instead.  A plain "2" runs to the end, the
  # control.
  add_test(NAME ume_mpi_negative_iteration_count
    COMMAND ${CMAKE_COMMAND}
      "-DARGS=lower-interior$<SEMICOLON>-i$<SEMICOLON>-1"
      "-DEXPECT=iteration count after -i .*got \"-1\".*Transport::abort: invalid iteration count"
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_trailing_garbage_iteration_count
    COMMAND ${CMAKE_COMMAND}
      "-DARGS=lower-interior$<SEMICOLON>-i$<SEMICOLON>2x"
      "-DEXPECT=iteration count after -i .*got \"2x\".*Transport::abort: invalid iteration count"
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  add_test(NAME ume_mpi_iteration_count
    COMMAND ${CMAKE_COMMAND}
      "-DARGS=lower-interior$<SEMICOLON>-i$<SEMICOLON>2"
      -DEXPECT_SUCCESS=ON
      "-DEXPECT=Checking gradient result.*Done\\."
      ${_ume_mpi_check}
    WORKING_DIRECTORY ${_zone_mask_dir})
  set_tests_properties(
    ume_mpi_no_interior_zone ume_mpi_interior_zones_in_lower_half
    ume_mpi_long_basename ume_mpi_truncated_mesh
    ume_mpi_negative_iteration_count ume_mpi_trailing_garbage_iteration_count
    ume_mpi_iteration_count PROPERTIES
    FIXTURES_REQUIRED ume_mpi_zone_mask_meshes
    TIMEOUT 90)
endif()
