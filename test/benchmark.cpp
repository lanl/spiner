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

using Spiner::DataBox;
using Spiner::RegularGrid1D;

using duration = std::chrono::nanoseconds;

constexpr Real KX = 2;
constexpr Real KY = 3;
constexpr Real KZ = 4;

constexpr Real xmin = 0;
constexpr Real xmax = 1;

PORTABLE_INLINE_FUNCTION
Real testFunction(Real z, Real y, Real x) {
  return sin(2 * M_PI * KX * x) * sin(2 * M_PI * KY * y) *
         sin(2 * M_PI * KZ * z);
}

int main(int argc, char *argv[]) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    if (argc < 3) {
      printf("Usage: %s ncoarse [nfine...]\n"
             "where nfine is any number of fnie grids, with nfine > ncoarse increasing\n",
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
      if ((i > 2) && nfine[i-2] <= nfine[i-2-1]) {
        std::cerr << "nfine must be increasing" << std::endl;
        std::exit(1);
      }
    }
    
    RegularGrid1D gcoarse(xmin, xmax, ncoarse);
    std::vector<RegularGrid1D> gfine;
    for (const auto & n : nfine) {
      gfine.push_back(RegularGrid1D(xmin, xmax, n));
    }

    std::cout << "# ncoarse = " << ncoarse << std::endl;
    std::unique_ptr<DataBox, Spiner::DBDeleter> pdb(
        new DataBox(Spiner::AllocationTarget::Device, ncoarse, ncoarse, ncoarse));
    for (int d = 0; d < pdb->rank(); d++) {
      pdb->setRange(d, xmin, xmax, ncoarse);
    }
    auto db = *pdb;
    portableFor("Filling databox", 0, ncoarse, 0, ncoarse, 0, ncoarse,
                PORTABLE_LAMBDA(const int iz, const int iy, const int ix) {
                  Real z = gcoarse.x(iz);
                  Real y = gcoarse.x(iy);
                  Real x = gcoarse.x(ix);
                  db(iz, iy, ix) = testFunction(z, y, x);
                });

    std::cout << "# nfine\ttime/point (us)\tL2 error" << std::endl;
    for (int ifine = 0; ifine < nfine.size(); ++ifine) {
      auto n = nfine[ifine];
      auto g = gfine[ifine];
      Real d3x = g.dx()*g.dx()*g.dx();
      Real L2_error;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto start = std::chrono::high_resolution_clock::now();     
      portableReduce("interp", 0, n, 0, n, 0, n,
                     PORTABLE_LAMBDA(const int iz, const int iy, const int ix,
                                     Real &reduce) {
                       Real z = g.x(iz);
                       Real y = g.x(iy);
                       Real x = g.x(ix);
                       Real f_interp = db.interpToReal(z, y, x);
                       Real f_true = testFunction(z, y, x);
                       Real difference = f_interp - f_true;
                       reduce += difference * difference * d3x;
                     }, L2_error);
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
      auto stop = std::chrono::high_resolution_clock::now();
      auto time = std::chrono::duration_cast<duration>(stop - start);
      L2_error = std::sqrt(L2_error);
      printf("%d\t%.14e\t%.14e\n",
             n, time.count() / std::pow(static_cast<Real>(n), 3), L2_error);
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return 0;
}
