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
#include <fstream>
#include <iostream>
#include <string>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/spiner_types.hpp>

using Spiner::DataBox;
using Spiner::RegularGrid1D;

const std::string outname = "convergence.dat";

constexpr Real KX = 2;
constexpr Real KY = 3;
constexpr Real KZ = 4;

constexpr Real xmin = 0;
constexpr Real xmax = 1;

constexpr int NGRIDS = 4;
Real errors[NGRIDS];
constexpr Real NCOARSE[NGRIDS] = {8, 32, 128, 512};
constexpr int NFINE = 1024;

inline Real testFunction(Real z, Real y, Real x) {
  return sin(2 * M_PI * KX * x) * sin(2 * M_PI * KY * y) *
         sin(2 * M_PI * KZ * z);
}

int main() {
  RegularGrid1D gfine(xmin, xmax, NFINE);

  std::cout << "Convergence test for Spiner" << std::endl;
  std::cout << std::scientific;
  std::cout.precision(15);

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    for (int i = 0; i < NGRIDS; i++) {
      int n = NCOARSE[i];
      Real *data = (Real *)PORTABLE_MALLOC(sizeof(Real) * n * n * n);
      std::cout << "\tCoarse resolution = " << n << std::endl;
      DataBox db(data, n, n, n);
      for (int d = 0; d < db.rank(); d++) {
        db.setRange(d, xmin, xmax, n);
      }
      std::cout << "\t\tFilling DataBox" << std::endl;
      portableFor(
          "filling databox", 0, n, 0, n, 0, n,
          PORTABLE_LAMBDA(const int iz, const int iy, const int ix) {
            Real z = db.range(2).x(iz);
            Real y = db.range(1).x(iy);
            Real x = db.range(0).x(ix);
            db(iz, iy, ix) = testFunction(z, y, x);
          });
      std::cout << "\t\tInterpolating..." << std::endl;
      errors[i] = 0;
      portableReduce(
          "interpolating", 0, NFINE, 0, NFINE, 0, NFINE,
          PORTABLE_LAMBDA(const int iz, const int iy, const int ix,
                          Real &reduced) {
            Real z = gfine.x(iz);
            Real y = gfine.x(iy);
            Real x = gfine.x(ix);
            Real f_interp = db.interpToReal(z, y, x);
            Real f_true = testFunction(z, y, x);
            Real difference = f_interp - f_true;
            reduced += difference * difference;
          },
          errors[i]);
      errors[i] = sqrt(fabs(errors[i])) / (NFINE * NFINE * NFINE);
      std::cout << "\t\t...error = " << errors[i] << std::endl;
      PORTABLE_FREE(data);
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  std::cout << "Saving convergence results to: " << outname << std::endl;
  std::ofstream outfile;
  outfile.open(outname);
  outfile.precision(15);
  outfile << std::scientific;
  outfile << "# [0]:resolution, [1]:error" << std::endl;
  for (int i = 0; i < NGRIDS; i++) {
    outfile << NCOARSE[i] << "\t" << errors[i] << std::endl;
  }
  outfile.close();
  std::cout << "Done!" << std::endl;

  return 0;
}
