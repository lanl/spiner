// © (or copyright) 2026. Triad National Security, LLC. All rights reserved.
// Generative AI was used to assist with writing this file.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <ports-of-call/nqt_math.hpp>
#include <ports-of-call/portability.hpp>
#include <spiner/interpolation.hpp>

using FastGrid = Spiner::FastNonUniformGrid1D<double>;
using BinaryGrid = Spiner::NonUniformGrid1D<double>;
using duration = std::chrono::nanoseconds;

constexpr double fake_grid[110] = {
    0.0e+00, 8.9e-06, 1.8e-05, 4.5e-05, 8.9e-05, 1.8e-04, 4.5e-04, 8.9e-04,
    1.8e-03, 4.5e-03, 8.9e-03, 1.3e-02, 2.2e-02, 3.6e-02, 5.4e-02, 8.9e-02,
    1.3e-01, 2.2e-01, 3.6e-01, 5.4e-01, 7.2e-01, 8.9e-01, 1.1e+00, 1.3e+00,
    1.6e+00, 1.8e+00, 2.2e+00, 2.7e+00, 3.1e+00, 3.6e+00, 4.0e+00, 4.5e+00,
    4.9e+00, 5.4e+00, 5.8e+00, 6.3e+00, 6.7e+00, 7.2e+00, 7.4e+00, 7.6e+00,
    7.8e+00, 8.0e+00, 8.3e+00, 8.5e+00, 8.7e+00, 8.8e+00, 8.9e+00, 9.0e+00,
    9.2e+00, 9.4e+00, 9.6e+00, 9.8e+00, 1.0e+01, 1.1e+01, 1.2e+01, 1.3e+01,
    1.4e+01, 1.5e+01, 1.6e+01, 1.7e+01, 1.8e+01, 2.0e+01, 2.2e+01, 2.5e+01,
    2.7e+01, 2.9e+01, 3.1e+01, 3.4e+01, 3.6e+01, 4.0e+01, 4.5e+01, 4.9e+01,
    5.4e+01, 5.8e+01, 6.3e+01, 6.7e+01, 7.2e+01, 7.6e+01, 8.0e+01, 8.5e+01,
    8.9e+01, 9.8e+01, 1.1e+02, 1.2e+02, 1.3e+02, 1.6e+02, 1.8e+02, 2.0e+02,
    2.2e+02, 2.7e+02, 3.1e+02, 3.6e+02, 4.5e+02, 5.4e+02, 6.3e+02, 7.2e+02,
    8.0e+02, 8.9e+02, 1.1e+03, 1.3e+03, 1.8e+03, 2.7e+03, 3.6e+03, 5.4e+03,
    7.2e+03, 8.9e+03, 1.8e+04, 4.5e+04, 8.9e+04, 1.8e+05};

PORTABLE_INLINE_FUNCTION double nqtSinh(const double x) {
#ifdef SPINER_USE_PORTABLE_NQT
  return PortsOfCall::NQT::O1::Portable::sinh(x);
#else
  return PortsOfCall::NQT::O1::Aliased::sinh(x);
#endif
}

template <typename Grid>
double timeLookups(const Grid &grid, const int nqueries, long long &checksum) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#endif
  const auto start = std::chrono::high_resolution_clock::now();
  portableReduce(
      "Nonuniform grid lookup benchmark", 0, nqueries,
      PORTABLE_LAMBDA(const int i, long long &sum) {
        const double fraction = static_cast<double>(i) / (nqueries - 1);
        const double query = grid.min() + fraction * (grid.max() - grid.min());
        sum += grid.index(query);
      },
      checksum);
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#endif
  const auto stop = std::chrono::high_resolution_clock::now();
  return static_cast<double>(
             std::chrono::duration_cast<duration>(stop - start).count()) /
         nqueries;
}

int profileLookups(const char *name, const std::vector<double> &points,
                   const int nqueries,
                   const FastGrid::Settings &settings = {}) {
  BinaryGrid binary_host(points);
  FastGrid fast_host(points, settings);
  BinaryGrid binary = binary_host.getOnDevice();
  FastGrid fast = fast_host.getOnDevice();

  long long binary_checksum = 0;
  long long fast_checksum = 0;
  const double binary_time = timeLookups(binary, nqueries, binary_checksum);
  const double fast_time = timeLookups(fast, nqueries, fast_checksum);

  const std::size_t lookup_entries = fast_host.lookupSize();
  const std::size_t lookup_bytes = lookup_entries * sizeof(int);
  std::printf("%s %zu %d %zu %zu %.8e %.8e %.8e\n", name, points.size(),
              nqueries, lookup_entries, lookup_bytes, binary_time, fast_time,
              binary_time / fast_time);
  const int result = binary_checksum == fast_checksum ? 0 : 2;
  if (result != 0) {
    std::fprintf(stderr, "%s lookup checksums differ: %lld != %lld\n", name,
                 binary_checksum, fast_checksum);
  }

  fast.finalize();
  binary.finalize();
  fast_host.finalize();
  binary_host.finalize();
  return result;
}

int main(int argc, char *argv[]) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  int result = 0;
  {
    const int npoints = argc > 1 ? std::atoi(argv[1]) : 1024;
    const int nqueries = argc > 2 ? std::atoi(argv[2]) : 1000000;
    if (npoints < 2 || nqueries < 2) {
      std::fprintf(stderr, "Usage: %s [npoints>=2] [nqueries>=2]\n", argv[0]);
      result = 1;
    } else {
      constexpr double scale = 2.0;
      constexpr double transformed_min = -12.0;
      constexpr double transformed_max = 12.0;
      std::vector<double> points(npoints);
      for (int i = 0; i < npoints; ++i) {
        const double transformed =
            transformed_min + i * (transformed_max - transformed_min) /
                                  static_cast<double>(npoints - 1);
        points[i] = scale * nqtSinh(transformed);
      }

      const std::vector<double> fake_points(fake_grid, fake_grid + 110);
      std::printf("# table points queries lookup_entries lookup_bytes "
                  "binary_ns fast_ns speedup\n");
      result |= profileLookups(
          "signed_log", points, nqueries,
          FastGrid::Settings{.scale = -1,
                             .policy = FastGrid::Policy::RequireFast,
                             .max_lookup_ratio =
                                 FastGrid::DEFAULT_MAX_LOOKUP_RATIO});
      result |= profileLookups("fake", fake_points, nqueries);
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}
