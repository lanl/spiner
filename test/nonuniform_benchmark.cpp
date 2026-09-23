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

      BinaryGrid binary_host(points);
      FastGrid fast_host(
          points, FastGrid::Settings{.scale = scale,
                                     .policy = FastGrid::Policy::RequireFast,
                                     .max_lookup_ratio =
                                         FastGrid::DEFAULT_MAX_LOOKUP_RATIO});
      BinaryGrid binary = binary_host.getOnDevice();
      FastGrid fast = fast_host.getOnDevice();

      long long binary_checksum = 0;
      long long fast_checksum = 0;
      const double binary_time = timeLookups(binary, nqueries, binary_checksum);
      const double fast_time = timeLookups(fast, nqueries, fast_checksum);

      std::printf(
          "# points queries lookup_entries binary_ns fast_ns speedup\n");
      std::printf("%d %d %zu %.8e %.8e %.8e\n", npoints, nqueries,
                  fast_host.lookupSize(), binary_time, fast_time,
                  binary_time / fast_time);
      if (binary_checksum != fast_checksum) {
        std::fprintf(stderr, "Lookup checksums differ: %lld != %lld\n",
                     binary_checksum, fast_checksum);
        result = 2;
      }

      fast.finalize();
      binary.finalize();
      fast_host.finalize();
      binary_host.finalize();
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}
