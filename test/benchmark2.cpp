// © (or copyright) 2019-2021. Triad National Security, LLC. All rights
// reserved.  This program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
// operated by Triad National Security, LLC for the U.S.  Department of
// Energy/National Nuclear Security Administration. All rights in the
// program are reserved by Triad National Security, LLC, and the
// U.S. Department of Energy/National Nuclear Security
// Administration. The Government is granted for itself and others acting
// on its behalf a nonexclusive, paid-up, irrevocable worldwide license
// in this material to reproduce, prepare derivative works, distribute
// copies to the public, perform publicly and display publicly, and to
// permit others to do so.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/spiner_types.hpp>

using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Spiner::DBDeleter;

using duration = std::chrono::nanoseconds;

constexpr Real KX = 2;
constexpr Real KY = 3;
constexpr Real KZ = 4;

constexpr Real xmin = 0;
constexpr Real xmax = 1;

// ------------------------------------------------------------------------------------------------

void fence() {
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
}

// ------------------------------------------------------------------------------------------------

PORTABLE_INLINE_FUNCTION
Real testFunction(Real z, Real y, Real x) {
  return sin(2 * M_PI * KX * x) * sin(2 * M_PI * KY * y) *
         sin(2 * M_PI * KZ * z);
}

// ------------------------------------------------------------------------------------------------

template <typename Callable>
void run_tests(
    const std::vector<int> & nfine,
    const std::vector<RegularGrid1D> & gfine,
    const std::string & message,
    Callable callable) {
  printf("# %s\n", message.c_str());
  printf("#  %-15s   %-15s   %-15s\n", "nfine", "time/point (ns)", "L2 error");
  for (int ifine = 0; ifine < nfine.size(); ++ifine) {
    auto n = nfine[ifine];
    auto g = gfine[ifine];
    Real d3x = g.dx() * g.dx() * g.dx();
    Real del = g.dx() / 3.0;
    Real L2_error;

    fence();
    auto start = std::chrono::high_resolution_clock::now();
    portableReduce(
        "interp", 0, n, 0, n, 0, n,
        PORTABLE_LAMBDA(const int iz, const int iy, const int ix, Real &reduce) {
          Real zz = g.x(iz);
          Real yy = g.x(iy);
          Real xx = g.x(ix);
          for (int sz = -8; sz <= 8; ++sz) {
            for (int sy = -8; sy <= 8; ++sy) {
              for (int sx = -8; sx <= 8; ++sx) {
                Real z = zz + sz * del;
                Real y = yy + sy * del;
                Real x = xx + sx * del;
                Real f_interp = callable(z, y, x);
                Real f_true = testFunction(z, y, x);
                Real difference = f_interp - f_true;
                reduce += difference * difference * d3x;
              }
            }
          }
        },
        L2_error);
    fence();

    auto stop = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<duration>(stop - start);
    L2_error = std::sqrt(L2_error);
    Real time_per_point = time.count() / std::pow(static_cast<Real>(n), 3);
    printf("   %-15d   %-15.6e   %-15.6e\n", n, time_per_point, L2_error);
  }
}

// ------------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    // Parse arguments
    if (argc < 3) {
      printf("Usage: %s ncoarse [nfine...]\n"
             "where nfine is any number of fine grids, with nfine > ncoarse "
             "increasing\n",
             argv[0]);
    }
    std::vector<int> nfine;
    const int ncoarse = std::atoi(argv[1]);
    for (int i = 2; i < argc; ++i) {
      nfine.push_back(std::atoi(argv[i]));
      if (nfine.back() <= ncoarse) {
        std::cerr << "nfine must be greater than ncoarse" << std::endl;
        std::exit(1);
      }
      if ((i > 2) && nfine[i - 2] <= nfine[i - 2 - 1]) {
        std::cerr << "nfine must be increasing" << std::endl;
        std::exit(1);
      }
      if (nfine.back() % ncoarse != 0) {
        std::cerr << "nfine must be an integer multiple of ncoarse" << std::endl;
        std::exit(1);
      }
    }

    // Build grid and databox
    RegularGrid1D gcoarse(xmin, xmax, ncoarse);
    std::vector<RegularGrid1D> gfine;
    for (const auto &n : nfine) {
      gfine.push_back(RegularGrid1D(xmin, xmax, n));
    }

    std::cout << "# ncoarse = " << ncoarse << std::endl;
    std::unique_ptr<DataBox, DBDeleter> pdb(new DataBox(
        Spiner::AllocationTarget::Device, ncoarse, ncoarse, ncoarse));
    for (int d = 0; d < pdb->rank(); d++) {
      pdb->setRange(d, xmin, xmax, ncoarse);
    }
    auto db = *pdb;
    portableFor(
        "Filling databox", 0, ncoarse, 0, ncoarse, 0, ncoarse,
        PORTABLE_LAMBDA(const int iz, const int iy, const int ix) {
          Real z = gcoarse.x(iz);
          Real y = gcoarse.x(iy);
          Real x = gcoarse.x(ix);
          db(iz, iy, ix) = testFunction(z, y, x);
        });

    // Test interpToReal()
    run_tests(nfine, gfine, "interpToReal 3D (hand-coded method)",
            PORTABLE_LAMBDA(Real z, Real y, Real x){ return db.interpToReal(z, y, x); });

    // Test interpolate()
    run_tests(nfine, gfine, "interpolate 3D (recursive method)",
            PORTABLE_LAMBDA(Real z, Real y, Real x){ return db.interpolate(z, y, x); });
  }

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return 0;
}
