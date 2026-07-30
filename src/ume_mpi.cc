/*
  Copyright (c) 2023, Triad National Security, LLC. All rights reserved.

  This is open source software; you can redistribute it and/or modify it under
  the terms of the BSD-3 License. If software is modified to produce derivative
  works, such modified software should be clearly marked, so as not to confuse
  it with the version available from LANL. Full text of the BSD-3 License can be
  found in the LICENSE.md file, and the full assertion of copyright in the
  NOTICE.md file.
*/

/*!
  \file ume_mpi.cc

  This is an example of MPI-based driver that reads one partition of an Ume
  binary mesh into each rank, and then performs (and tests) a gradient
  operation.

  Note that there must be as many *.ume files as there are MPI ranks, and they
  should have filenames of the form '<basename>.<pe>.ume', where <basename> is
  an arbitray string provided on the command line, and <pe> is a rank number
  with a printf format of "%05d" (zero-filled, five digits)

  ### Comparing a partitioned run against a serial one

  `--dump-grad` and `--compare` together answer the question a decomposition
  cannot answer about itself: does the mesh it produced still compute the same
  physics?  Everything else in the suite checks that a partition is *well
  formed* -- the audit, the cross-rank check, the round trip -- and a
  decomposition can be all of those and still be wrong.

  A one-rank partition is bit-identical to the mesh it came from, so
  `mpiexec -n 1 ume_mpi --dump-grad ref` *is* the serial answer, computed by
  this same binary along this same code path.  That removes the need for a
  second reference implementation, and for `ume_serial` to grow anything.

  The one thing a rank cannot do locally is agree with the reference about
  which zone to perturb: "the middle interior zone" names a different zone
  in every decomposition.  So the reference records the *serial* index of the
  zone it spiked, `--compare` reads it back, and each rank finds that zone
  through its `.gid` side-car -- which is also what lines the answers up zone
  for zone afterwards.

  The tolerance is relative and is not zero.  Partitioning genuinely reorders
  the accumulation in `gradzatz`, so the answers agree to roundoff rather than
  to the bit; `test_brick.cc`'s exact comparison would be the wrong test here.
*/

#include "Ume/Comm_MPI.hh"
#include "Ume/Global_Ids.hh"
#include "Ume/SOA_Idx_Mesh.hh"
#include "Ume/Timer.hh"
#include "Ume/face_area.hh"
#include "Ume/gradient.hh"
#include "Ume/process_mgmt.hh"
#include "Ume/renumbering.hh"
#include "Ume/utils.hh"
#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

using Mesh = Ume::SOA_Idx::Mesh;
using DBLV_T = typename Ume::DS_Types::DBLV_T;
using VEC3V_T = typename Ume::DS_Types::VEC3V_T;
using VEC3_T = typename Ume::DS_Types::VEC3_T;

namespace {

struct Options {
  std::string basename;
  std::size_t iterations{1};
  //! Write the zone gradients of this run out as a reference
  std::string dump_grad;
  //! Check the zone gradients of this run against a reference
  std::string compare;
  //! Allowed error, relative to the largest gradient in the reference
  double tol{1e-12};
};

//! One reference run's zone gradients, keyed by serial zone index
struct Reference {
  int spike_gid{-1};
  std::map<int, VEC3_T> zgrad;
  //! The largest gradient magnitude anywhere in the reference
  double scale{0.0};
};

void usage(std::ostream &os) {
  os << R"(Usage: ume_mpi [options] <basename>

Reads <basename>.NNNNN.ume on each rank and takes a zone-centered gradient.

Options:
  -i N               Repeat each timed kernel N times [1]
  --dump-grad FILE   Write this run's zone gradients to FILE as a reference.
                     Only meaningful on one rank, where the partition is the
                     serial mesh; requires <basename>.NNNNN.gid.
  --compare FILE     Spike the zone FILE says it spiked, and check this run's
                     zone gradients against it.  Requires the .gid side-car.
  --tol F            Error allowed by --compare, relative to the largest
                     gradient in the reference [1e-12]
  -h, --help         This message
)";
}

double norm(VEC3_T const &v) { return std::sqrt(dotprod(v, v)); }

template <typename T> bool parse(std::string_view text, T &out) {
  auto const last = text.data() + text.size();
  auto const [ptr, ec] = std::from_chars(text.data(), last, out);
  return ec == std::errc{} && ptr == last;
}

//! What reading the command line turned out to mean
/*! Nothing here may std::exit(): MPI is already up by the time the command
    line is read, and MPI_Finalize() is reached only by returning to main and
    running out through comm.stop(). */
enum class Parse { RUN, HELP, ERROR };

//! Says why, on anything but RUN
Parse parse_args(int argc, char *argv[], int numpe, Options &opt) {
  for (int i = 1; i < argc; ++i) {
    std::string_view const arg{argv[i]};
    //! The next argument, or nullptr having said that it is missing
    auto next = [&](char const *what) -> char const * {
      if (i + 1 >= argc) {
        std::cerr << "ume_mpi: " << what << " requires an argument\n";
        return nullptr;
      }
      return argv[++i];
    };

    if (arg == "-h" || arg == "--help") {
      usage(std::cout);
      return Parse::HELP;
    } else if (arg == "-i") {
      char const *const v = next("-i");
      if (!v)
        return Parse::ERROR;
      if (!parse(v, opt.iterations)) {
        std::cerr << "ume_mpi: -i wants a count\n";
        return Parse::ERROR;
      }
    } else if (arg == "--dump-grad") {
      char const *const v = next("--dump-grad");
      if (!v)
        return Parse::ERROR;
      opt.dump_grad = v;
    } else if (arg == "--compare") {
      char const *const v = next("--compare");
      if (!v)
        return Parse::ERROR;
      opt.compare = v;
    } else if (arg == "--tol") {
      char const *const v = next("--tol");
      if (!v)
        return Parse::ERROR;
      if (!parse(v, opt.tol)) {
        std::cerr << "ume_mpi: --tol wants a number\n";
        return Parse::ERROR;
      }
    } else if (arg.starts_with("-") && arg.size() > 1) {
      std::cerr << "ume_mpi: unknown option \"" << arg << "\"\n";
      return Parse::ERROR;
    } else if (opt.basename.empty()) {
      opt.basename = arg;
    } else {
      std::cerr << "ume_mpi: one mesh basename, not two\n";
      return Parse::ERROR;
    }
  }

  if (opt.basename.empty()) {
    usage(std::cerr);
    return Parse::ERROR;
  }
  if (!opt.dump_grad.empty() && !opt.compare.empty()) {
    std::cerr << "ume_mpi: --dump-grad and --compare are two different runs\n";
    return Parse::ERROR;
  }
  /* A reference has to be the whole answer.  On several ranks each one holds
     a part of it, and gathering them would be a second, unchecked code path
     for producing something we already get for free at -n 1. */
  if (!opt.dump_grad.empty() && numpe != 1) {
    std::cerr << "ume_mpi: --dump-grad wants one rank, so that the partition "
                 "it reads is the serial mesh; this run has "
              << numpe << '\n';
    return Parse::ERROR;
  }
  return Parse::RUN;
}

bool read_mesh(std::string const &basename, int const mype, Mesh &mesh) {
  std::string const fname = Ume::mesh_filename(basename, mype);
  std::ifstream is(fname);
  if (!is) {
    std::cerr << "Unable to open file \"" << fname << "\" for reading."
              << std::endl;
    return false;
  }
  mesh.read(is);
  is.close();
  if (!is) {
    std::cerr << "Unable to read a mesh from file \"" << fname << "\"."
              << std::endl;
    return false;
  }
  return true;
}

//! What this rank's zones were called in the mesh the partition came from
bool read_zone_gids(std::string const &basename, int const mype,
    Mesh const &mesh, std::vector<int> &zgid) {
  std::string const fname = Ume::Part::rank_file(basename, mype, "gid");
  std::ifstream is(fname);
  if (!is) {
    std::cerr << "Unable to open file \"" << fname << "\" for reading.  Write "
              << "it with `ume_partition --global-ids`." << std::endl;
    return false;
  }
  try {
    Ume::Part::Global_Ids const g = Ume::Part::read_global_ids(is);
    Ume::Part::Kind_Ids const &zones = g[Ume::Part::ZONE];
    /* A side-car from a different partitioning would name the wrong zones
       without ever being the wrong shape, so the sizes are worth checking. */
    if (g.mype != mype || zones.local != mesh.zones.local_size() ||
        zones.gid.size() != static_cast<std::size_t>(mesh.zones.size())) {
      std::cerr << "\"" << fname << "\" does not describe this mesh."
                << std::endl;
      return false;
    }
    zgid = zones.gid;
  } catch (std::exception const &e) {
    std::cerr << "\"" << fname << "\": " << e.what() << std::endl;
    return false;
  }
  return true;
}

bool write_reference(std::string const &path, Mesh const &mesh,
    std::vector<int> const &zgid, int const spike, VEC3V_T const &zgrad) {
  std::ofstream os(path);
  if (!os) {
    std::cerr << "Unable to open file \"" << path << "\" for writing."
              << std::endl;
    return false;
  }

  std::vector<int> real;
  for (int z = 0; z < mesh.zones.local_size(); ++z)
    if (mesh.zones.mask[static_cast<std::size_t>(z)] >= 1)
      real.push_back(z);

  os << "ume_grad 1\n";
  os << "spike " << zgid[static_cast<std::size_t>(spike)] << '\n';
  os << "zones " << real.size() << '\n';
  /* 17 significant digits is what a double round trips through decimal in,
     so the reference is exact and only the comparison has a tolerance. */
  os << std::scientific << std::setprecision(17);
  for (int const z : real) {
    VEC3_T const &g = zgrad[static_cast<std::size_t>(z)];
    os << zgid[static_cast<std::size_t>(z)] << ' ' << g[0] << ' ' << g[1] << ' '
       << g[2] << '\n';
  }
  return static_cast<bool>(os);
}

bool read_reference(std::string const &path, Reference &ref) {
  std::ifstream is(path);
  if (!is) {
    std::cerr << "Unable to open file \"" << path << "\" for reading."
              << std::endl;
    return false;
  }

  std::string tok;
  int version{}, nzones{};
  if (!(is >> tok >> version) || tok != "ume_grad" || version != 1) {
    std::cerr << "\"" << path << "\" is not a version 1 gradient reference."
              << std::endl;
    return false;
  }
  if (!(is >> tok >> ref.spike_gid) || tok != "spike" ||
      !(is >> tok >> nzones) || tok != "zones") {
    std::cerr << "\"" << path << "\" has a malformed header." << std::endl;
    return false;
  }

  for (int i = 0; i < nzones; ++i) {
    int gid{};
    VEC3_T g;
    if (!(is >> gid >> g[0] >> g[1] >> g[2])) {
      std::cerr << "\"" << path << "\" ends after " << i << " of " << nzones
                << " zones." << std::endl;
      return false;
    }
    ref.zgrad.emplace(gid, g);
    ref.scale = std::max(ref.scale, norm(g));
  }
  return true;
}

//! Check this rank's real zones against the reference; returns the mismatches
/*! Only the zones this rank holds are checked, because only those are the
    ones it computed.  That the ranks between them cover every zone of the
    serial mesh is what `check_partitions` establishes, so it is not
    re-established here. */
int compare_reference(Reference const &ref, Mesh const &mesh,
    std::vector<int> const &zgid, VEC3V_T const &zgrad, double const tol,
    std::ostream &os) {
  double const allowed = tol * ref.scale;
  int bad = 0, checked = 0;
  double worst = 0.0;

  for (int z = 0; z < mesh.zones.local_size(); ++z) {
    if (mesh.zones.mask[static_cast<std::size_t>(z)] < 1)
      continue;
    int const gid = zgid[static_cast<std::size_t>(z)];
    auto const it = ref.zgrad.find(gid);
    if (it == ref.zgrad.end()) {
      os << "PE" << mesh.mype << " zone " << z << " (global " << gid
         << ") is not in the reference\n";
      bad += 1;
      continue;
    }
    checked += 1;
    double const err = norm(zgrad[static_cast<std::size_t>(z)] - it->second);
    worst = std::max(worst, err);
    if (err > allowed) {
      os << "PE" << mesh.mype << " zone " << z << " (global " << gid
         << "): " << zgrad[static_cast<std::size_t>(z)] << " vs. " << it->second
         << ", off by " << err << '\n';
      bad += 1;
    }
  }

  os << "PE" << mesh.mype << ": " << checked << " zones compared, worst error "
     << worst << " of " << allowed << " allowed" << std::endl;
  return bad;
}

// Unused: its call in main is commented out.
[[maybe_unused]] bool test_point_gathscat(Mesh &mesh);
void check_gradzatz_diffs(Mesh const &mesh, int const &centered_zone_index,
    VEC3V_T const &zgrad, VEC3V_T const &zgrad_invert, VEC3V_T const &pgrad,
    VEC3V_T const &pgrad_invert);
void report_gradient_mismatches(int const mype,
    Ume::SOA_Idx::Entity const &entity, char const tag, std::string_view label,
    VEC3V_T const &grad, VEC3V_T const &grad_invert, double const tol);

} // namespace

int main(int argc, char *argv[]) {
  /* Initialize MPI and instantiate the MPI Transport. */
  Ume::Comm::MPI comm(&argc, &argv);

  /* Initialize Kokkos and the memory pool. This should happen right
   * after the call to MPI_Init for best performance -- and it has to happen
   * before the command line is read at all: Ume::initialize() takes argc by
   * reference and strips the --kokkos-* arguments, which parse_args would
   * otherwise reject as unknown options. */
  Ume::initialize(argc, argv);

  /* Nothing below may leave main without running out through the teardown at
     the bottom: MPI_Finalize() is reached only from comm.stop(), and a rank
     that skips it is a rank mpiexec reports as having died, even on --help.
     comm.abort() is the one legitimate early exit, and it is not one of
     these -- every check here fails identically on every rank, so no rank is
     left waiting on another. */
  auto teardown = [&comm](int status) {
    Ume::finalize();
    comm.stop();
    return status;
  };

  /* A missing basename is one of parse_args's own errors, which is why there
     is no separate argc check here. */
  Options opt;
  switch (parse_args(argc, argv, comm.numpe(), opt)) {
  case Parse::RUN:
    break;
  case Parse::HELP:
    return teardown(EXIT_SUCCESS);
  case Parse::ERROR:
    return teardown(EXIT_FAILURE);
  }
  size_t const ic = opt.iterations;
  bool const referencing = !opt.dump_grad.empty() || !opt.compare.empty();

  Reference ref;
  if (!opt.compare.empty() && !read_reference(opt.compare, ref)) {
    std::cerr << "Aborting." << std::endl;
    return teardown(EXIT_FAILURE);
  }

  /* Create a mesh instance and attach the communicator to the mesh. */
  Mesh mesh;
  mesh.comm = &comm;

  if (comm.pe() == 0)
    std::cout << "Initializing mesh..." << std::endl;

  /* Read the data file.  A rank that cannot read its partition has to take the
     whole job down: returning here would leave the other ranks waiting in the
     gathscat below, which hangs, or is killed by the launcher, rather than
     fails. */
  if (!read_mesh(opt.basename, comm.pe(), mesh)) {
    comm.abort("this rank could not read its mesh file");
    return EXIT_FAILURE; // not reached; abort() is not declared [[noreturn]]
  }

  /* The side-car is read only when something is going to use it, so an
     ordinary run still needs nothing but the .ume files. */
  std::vector<int> zgid;
  if (referencing && !read_zone_gids(opt.basename, comm.pe(), mesh, zgid)) {
    std::cerr << "Aborting." << std::endl;
    return EXIT_FAILURE;
  }

  /* This allows us to attach a debugger to a single rank specified in the
     UME_DEBUG_RANK environment variable. */
  Ume::debug_attach_point(comm.pe());

  /*
  if (test_point_gathscat(mesh)) {
    std::cout << comm.id() << ": test_point_gathscat PASS" << std::endl;
  } else {
    std::cout << comm.id() << ": test_point_gathscat FAIL" << std::endl;
  }
  */

  if (comm.pe() == 0)
    std::cout << "Creating zone field..." << std::endl;

  /* Find an interior local zone to set a value in. NOTE zone types
   * -1: exterior zone
   *  0: null zone or parallel ghost zone (halo)
   *  1: interior zone
   *
   * Take the middle of the interior zones rather than scanning upward from
   * the middle of the whole zone range: where the interior zones sit is a
   * property of the mesh generator, not something a driver can assume.
   * `gen_brick` appends every boundary zone after every real one, so on a
   * rank where boundary zones are at least half of the local zones the upper
   * half of the range is entirely mask == -1 and the scan runs off the end --
   * leaving czi == local_size(), which asserts in a Debug build and in a
   * Release build spikes a ghost zone, or writes past `zfield` entirely on a
   * rank that has no ghosts. */
  auto const &kztyp = mesh.zones.mask;
  int const zl = mesh.zones.local_size();
  auto const is_interior = [](short const mask) { return mask >= 1; };
  auto const num_interior =
      std::ranges::count_if(kztyp.begin(), kztyp.begin() + zl, is_interior);
  if (num_interior == 0) {
    /* One rank returning while the others wait in a gathscat would hang the
       job, so this has to take the whole job down.  The partitioning allows a
       rank with no zones, but this driver seeds its gradient in an interior
       zone of every rank. */
    comm.abort("this rank holds no interior zones: ume_mpi cannot run an "
               "over-decomposed mesh, since every rank needs an interior "
               "zone to seed its gradient in");
    return EXIT_FAILURE; // not reached; abort() is not declared [[noreturn]]
  }

  /* Under --compare the zone to spike is not this rank's to choose: the two
     runs have to perturb the same zone of the *serial* mesh, and only one
     rank of this run holds it.  Everywhere else czi is -1 and this rank
     simply contributes its share of a gradient it did not source. */
  int czi = -1;
  if (!opt.compare.empty()) {
    for (int z = 0; z < zl; ++z)
      if (zgid[static_cast<std::size_t>(z)] == ref.spike_gid &&
          is_interior(kztyp[z])) {
        czi = z;
        break;
      }
    if (czi >= 0)
      std::cout << "PE" << comm.pe() << " holds the reference's zone "
                << ref.spike_gid << " at index " << czi << std::endl;
  } else {
    /* The middle interior zone, located by counting rather than by building
       the list of them: on a production rank that list would be millions of
       ints, held only to read one of them. */
    for (int z = 0, seen = 0; z < zl; ++z) {
      if (is_interior(kztyp[z]) && seen++ == num_interior / 2) {
        czi = z;
        break;
      }
    }
    assert(czi >= 0);
  }

  /* Create a zone-field that is everywhere zero but in the centered
   * zone at index czi.  Ghost zones are left at zero even where they shadow
   * the spiked one: `gradzatz` reads the field only through corners with
   * `mask >= 1`, and a ghost zone has none. */
  DBLV_T zfield(mesh.zones.size(), 0.0);
  if (czi >= 0)
    zfield[czi] = 100000.0;

  if (comm.pe() == 0)
    std::cout << "Computing zone-centered gradient..." << std::endl;

  /* Do a parallel zone-centered gradient computation on that field.
   * Here we run the computation once to initialize mesh structures
   * and time only the second run. We then run the same computations
   * using inverted connectivities. */
  VEC3V_T pgrad, zgrad;
  Ume::Timer orig_time;
  Ume::gradzatz(mesh, zfield, zgrad, pgrad);
  orig_time.start();
  for (size_t i = 0; i < ic; i++) {
    Ume::gradzatz(mesh, zfield, zgrad, pgrad);
  }
  orig_time.stop();

  VEC3V_T pgrad_invert, zgrad_invert;
  Ume::Timer invert_time;
  Ume::gradzatz_invert(mesh, zfield, zgrad_invert, pgrad_invert);
  invert_time.start();
  for (size_t i = 0; i < ic; i++) {
    Ume::gradzatz_invert(mesh, zfield, zgrad_invert, pgrad_invert);
  }
  invert_time.stop();

  if (comm.pe() == 0) {
    std::cout << "Original algorithm took: " << orig_time.seconds() << "s\n";
    std::cout << "Inverted algorithm took: " << invert_time.seconds() << "s\n";
    std::cout << "Checking gradient result..." << std::endl;
  }

  /* Double check that the gradients are non-zero where we expect.  This is
     stated in terms of the spiked zone's neighborhood, so it has nothing to
     say on a rank that does not hold the spiked zone. */
  if (czi >= 0)
    check_gradzatz_diffs(mesh, czi, zgrad, zgrad_invert, pgrad, pgrad_invert);

  int status = EXIT_SUCCESS;
  if (!opt.dump_grad.empty()) {
    if (write_reference(opt.dump_grad, mesh, zgid, czi, zgrad))
      std::cout << "Wrote the reference gradients to \"" << opt.dump_grad
                << "\"" << std::endl;
    else
      status = EXIT_FAILURE;
  } else if (!opt.compare.empty()) {
    /* Every rank runs this and every rank returns its own verdict: mpiexec
       fails the job if any one of them does, and no rank leaves early, so
       nothing is left waiting in a collective. */
    if (compare_reference(ref, mesh, zgid, zgrad, opt.tol, std::cout) > 0)
      status = EXIT_FAILURE;
  }

  if (comm.pe() == 0)
    std::cout << "Computing face areas..." << std::endl;

  /* Do a face area computation on the mesh making use of the MPI comm
   * stencil. First create a result vector and initialize to an
   * impossible value. */
  DBLV_T face_area(mesh.faces.size(), -100000.0);
  Ume::Timer face_time;
  Ume::calc_face_area(mesh, face_area);
  face_time.start();
  for (size_t i = 0; i < ic; i++) {
    Ume::calc_face_area(mesh, face_area);
  }
  face_time.stop();

  if (comm.pe() == 0)
    std::cout << "Face area computation took: " << face_time.seconds() << "s\n";

  if (mesh.ivtag >= UME_VERSION_2) {
    if (comm.pe() == 0)
      std::cout << "Renumbering mesh entities..." << std::endl;

    if (mesh.dump_iotas) {
      Ume::Timer renumber_time;
      Ume::renumber_mesh(mesh);
      renumber_time.start();
      for (size_t i = 0; i < ic; i++) {
        Ume::renumber_mesh(mesh);
      }
      renumber_time.stop();

      if (comm.pe() == 0)
        std::cout << "Renumbering algorithm took: " << renumber_time.seconds()
                  << "s\n";
    } else {
      if (comm.pe() == 0)
        std::cout << "Iotas must be present in the mesh for renumbering."
                  << " Skipping renumbering..." << std::endl;
    }
  }

  if (comm.pe() == 0)
    std::cout << "Done." << std::endl;

  return teardown(status);
}

namespace {

bool test_point_gathscat(Mesh &mesh) {
  int const mype = mesh.comm->id();

  /* Initialize an integer field to the value 10 everywhere except:
   *   a copy = -1
   *   a source = the number of copies of this entity
   * A sum will set the source to zero, which should then get scattered
   * out to the copies. */
  int const background_val = 100;
  typename Ume::DS_Types::INTV_T int_field(mesh.points.size(), background_val);

  for (auto const &n : mesh.points.myCpys) {
    for (auto const &e : n.elements) {
      assert(int_field[e] == background_val);
      int_field[e] = -1;
    }
  }

  /* We're "pre-summing" the source counts so that we can check that
   * the source values in the int_field are set to the background
   * value. */
  std::map<int, int> src_counts;
  for (auto const &n : mesh.points.mySrcs) {
    for (auto const &e : n.elements) {
      auto [it, success] = src_counts.insert(std::make_pair(e, 0));
      it->second += 1;
    }
  }

  for (auto const &s : src_counts) {
    assert(int_field[s.first] == background_val);
    int_field[s.first] = s.second;
  }

  bool result = true;

  mesh.points.gathscat(Ume::Comm::Op::SUM, int_field);
  for (int i = 0; i < mesh.points.local_size(); ++i) {
    if (mesh.points.comm_type[i] == Ume::SOA_Idx::Entity::INTERNAL) {
      if (int_field[i] != background_val) {
        std::cout << "Rank " << mype << " expecting int_field[" << i
                  << "] == " << background_val << ", got " << int_field[i]
                  << '\n';
        result = false;
      }
    } else if (mesh.points.comm_type[i] == Ume::SOA_Idx::Entity::SOURCE ||
        mesh.points.comm_type[i] == Ume::SOA_Idx::Entity::COPY) {
      if (int_field[i] != 0) {
        std::cout << "Rank " << mype << " expecting int_field[" << i
                  << "] == 0, got " << int_field[i] << '\n';
        result = false;
      }
    }
  }

  return result;
}

void check_gradzatz_diffs(Mesh const &mesh, int const &centered_zone_index,
    VEC3V_T const &zgrad, VEC3V_T const &zgrad_invert, VEC3V_T const &pgrad,
    VEC3V_T const &pgrad_invert) {
  auto const &kztyp = mesh.zones.mask;

  constexpr double tol = 1e-6;
  report_gradient_mismatches(mesh.mype, mesh.zones, 'Z',
      "zgrad != zgrad_invert", zgrad, zgrad_invert, tol);
  report_gradient_mismatches(mesh.mype, mesh.points, 'P',
      "pgrad != pgrad_invert", pgrad, pgrad_invert, tol);

  auto const &z2pz = mesh.ds->caccess_intrr("m:z>pz");
  auto const &z2p = mesh.ds->caccess_intrr("m:z>p");
  auto const &kptyp = mesh.points.mask;
  std::vector<int> grad_zones;
  for (int z = 0; z < mesh.zones.size(); ++z) {
    if (z == centered_zone_index || kztyp[z] < 1)
      continue;
    if (zgrad[z] != 0.0)
      grad_zones.push_back(z);
  }

  std::vector<int> grad_points;
  for (int p = 0; p < mesh.points.size(); ++p) {
    if (kptyp[p] > 0 && pgrad[p] != 0.0)
      grad_points.push_back(p);
  }

  std::sort(grad_zones.begin(), grad_zones.end());
  std::sort(grad_points.begin(), grad_points.end());
  std::vector<int> diff;
  std::ranges::set_difference(
      grad_zones, z2pz[centered_zone_index], std::back_inserter(diff));
  if (!diff.empty()) {
    std::cout << "PE" << mesh.mype << " zone diff " << diff.size() << " found "
              << grad_zones.size() << " expected "
              << z2pz.size(centered_zone_index) << '\n';
  }

  diff.clear();
  std::ranges::set_difference(
      grad_points, z2p[centered_zone_index], std::back_inserter(diff));
  if (!diff.empty()) {
    std::cout << "PE" << mesh.mype << " pt diff " << diff.size() << " found "
              << grad_points.size() << " expected "
              << z2p.size(centered_zone_index) << '\n';
  }
}

/* Print each element of `entity` whose two gradients differ by more than `tol`
   in the L2 norm.  Any one component off by more than `tol` puts the norm past
   it, so a component-wise test would add nothing.  A NaN difference compares
   false and is not reported. */
void report_gradient_mismatches(int const mype,
    Ume::SOA_Idx::Entity const &entity, char const tag, std::string_view label,
    VEC3V_T const &grad, VEC3V_T const &grad_invert, double const tol) {
  if (grad == grad_invert)
    return;
  for (int i = 0; i < entity.size(); ++i) {
    if (Ume::vectormag(grad[i] - grad_invert[i]) > tol) {
      std::cout << "PE" << mype << ' ' << label << ' ' << tag << i << ' '
                << entity.mask[i] << ": " << grad[i] << " vs. "
                << grad_invert[i] << '\n';
    }
  }
}

} // namespace
