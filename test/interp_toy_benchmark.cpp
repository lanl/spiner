// © (or copyright) 2026. Triad National Security, LLC. All rights
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

// Generative AI was used to assist with modifications to this file.

// Compares Spiner's hand-unrolled DataBox::interpToReal against the
// generic interp_toy machinery on identical random query points.
//
// Usage: interp_toy_benchmark.bin [ncoarse [nidx [npoints [ntrials]]]]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interp_toy.hpp>
#include <spiner/spiner_types.hpp>

using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Spiner::interp::at;
using Spiner::interp::interpolate_with;
using Spiner::interp::linear;

constexpr Real xmin = 0;
constexpr Real xmax = 1;

PORTABLE_INLINE_FUNCTION Real testFunction(Real z, Real y, Real x) {
  return std::sin(2 * M_PI * 2 * x) * std::sin(2 * M_PI * 3 * y) *
         std::sin(2 * M_PI * 4 * z);
}

// Query points, stored in device memory.
struct Points {
  int n;
  Real *x, *y, *z;
  int *idx;
};

template <typename T>
T *toDevice(const std::vector<T> &host) {
  const std::size_t bytes = host.size() * sizeof(T);
  T *device = (T *)PORTABLE_MALLOC(bytes);
  portableCopyToDevice(device, host.data(), bytes);
  return device;
}

Points makePoints(const int n, const int nidx) {
  std::mt19937_64 rng(12345);
  std::uniform_real_distribution<Real> coord(xmin, xmax);
  std::uniform_int_distribution<int> index(0, nidx - 1);
  std::vector<Real> x(n), y(n), z(n);
  std::vector<int> idx(n);
  for (int i = 0; i < n; ++i) {
    x[i] = coord(rng);
    y[i] = coord(rng);
    z[i] = coord(rng);
    idx[i] = index(rng);
  }
  return {n, toDevice(x), toDevice(y), toDevice(z), toDevice(idx)};
}

void freePoints(Points &p) {
  PORTABLE_FREE(p.x);
  PORTABLE_FREE(p.y);
  PORTABLE_FREE(p.z);
  PORTABLE_FREE(p.idx);
}

inline void fence() {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#endif
}

// Sum kernel(i) over all points. Returns ns per point.
template <typename Kernel>
double timeKernel(const int n, const Kernel &kernel, Real &sum) {
  sum = 0;
  fence();
  const auto start = std::chrono::steady_clock::now();
  portableReduce(
      "interp_toy_benchmark", 0, n,
      PORTABLE_LAMBDA(const int i, Real &acc) { acc += kernel(i); }, sum);
  fence();
  const auto stop = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::nano>(stop - start).count() / n;
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t m = v.size() / 2;
  return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

// Times hand and toy kernels, alternating which runs first in each
// trial, and checks that they agree.
template <typename Hand, typename Toy>
void runCase(const char *label, const std::size_t table_bytes, const int n,
             const int ntrials, const Hand &hand, const Toy &toy) {
  Real hand_sum, toy_sum, sq_diff;
  timeKernel(n, hand, hand_sum); // warm up
  timeKernel(n, toy, toy_sum);
  timeKernel(
      n, PORTABLE_LAMBDA(const int i) {
        const Real d = hand(i) - toy(i);
        return d * d;
      },
      sq_diff);

  std::vector<double> th, tt;
  for (int t = 0; t < ntrials; ++t) {
    if (t % 2 == 0) {
      th.push_back(timeKernel(n, hand, hand_sum));
      tt.push_back(timeKernel(n, toy, toy_sum));
    } else {
      tt.push_back(timeKernel(n, toy, toy_sum));
      th.push_back(timeKernel(n, hand, hand_sum));
    }
  }
  const double mh = median(th);
  const double mt = median(tt);
  std::printf("%-12s %10.1f %9.3f %9.3f %9.3f %9.3f %8.3f %10.2e\n", label,
              table_bytes / 1024.0, *std::min_element(th.begin(), th.end()),
              mh, *std::min_element(tt.begin(), tt.end()), mt, mt / mh,
              std::sqrt(sq_diff / n));
}

int main(int argc, char *argv[]) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    const int ncoarse = argc > 1 ? std::atoi(argv[1]) : 32;
    const int nidx = argc > 2 ? std::atoi(argv[2]) : 8;
    const int npoints = argc > 3 ? std::atoi(argv[3]) : (1 << 22);
    const int ntrials = argc > 4 ? std::atoi(argv[4]) : 11;
    const int n = ncoarse;

    std::printf("# ncoarse = %d, nidx = %d, npoints = %d, ntrials = %d\n", n,
                nidx, npoints, ntrials);
    std::printf("# times in ns/point; ratio = toy median / hand median; "
                "rms = rms(hand - toy)\n");
    std::printf("%-12s %10s %9s %9s %9s %9s %8s %10s\n", "# case",
                "table_KiB", "hand_min", "hand_med", "toy_min", "toy_med",
                "ratio", "rms");

    Points p = makePoints(npoints, nidx);
    const RegularGrid1D g(xmin, xmax, n);

    { // 2D
      DataBox db(Spiner::AllocationTarget::Device, n, n);
      for (int d = 0; d < 2; ++d)
        db.setRange(d, xmin, xmax, n);
      portableFor(
          "fill 2D", 0, n, 0, n, PORTABLE_LAMBDA(const int iy, const int ix) {
            db(iy, ix) = testFunction(0.25, g.x(iy), g.x(ix));
          });
      runCase(
          "2D", db.sizeBytes(), npoints, ntrials,
          PORTABLE_LAMBDA(const int i) {
            return db.interpToReal(p.y[i], p.x[i]);
          },
          PORTABLE_LAMBDA(const int i) {
            return interpolate_with(db, linear(db.range(1), p.y[i]),
                                    linear(db.range(0), p.x[i]));
          });
      free(db);
    }

    { // 2D with a fastest-moving index
      DataBox db(Spiner::AllocationTarget::Device, n, n, nidx);
      for (int d = 1; d < 3; ++d)
        db.setRange(d, xmin, xmax, n);
      portableFor(
          "fill 2D+idx", 0, n, 0, n, 0, nidx,
          PORTABLE_LAMBDA(const int iy, const int ix, const int i) {
            db(iy, ix, i) = (1 + i) * testFunction(0.25, g.x(iy), g.x(ix));
          });
      runCase(
          "2D+idx", db.sizeBytes(), npoints, ntrials,
          PORTABLE_LAMBDA(const int i) {
            return db.interpToReal(p.y[i], p.x[i], p.idx[i]);
          },
          PORTABLE_LAMBDA(const int i) {
            return interpolate_with(db, linear(db.range(2), p.y[i]),
                                    linear(db.range(1), p.x[i]), at(p.idx[i]));
          });
      free(db);
    }

    { // 3D
      DataBox db(Spiner::AllocationTarget::Device, n, n, n);
      for (int d = 0; d < 3; ++d)
        db.setRange(d, xmin, xmax, n);
      portableFor(
          "fill 3D", 0, n, 0, n, 0, n,
          PORTABLE_LAMBDA(const int iz, const int iy, const int ix) {
            db(iz, iy, ix) = testFunction(g.x(iz), g.x(iy), g.x(ix));
          });
      runCase(
          "3D", db.sizeBytes(), npoints, ntrials,
          PORTABLE_LAMBDA(const int i) {
            return db.interpToReal(p.z[i], p.y[i], p.x[i]);
          },
          PORTABLE_LAMBDA(const int i) {
            return interpolate_with(db, linear(db.range(2), p.z[i]),
                                    linear(db.range(1), p.y[i]),
                                    linear(db.range(0), p.x[i]));
          });
      free(db);
    }

    { // 3D with a fastest-moving index
      DataBox db(Spiner::AllocationTarget::Device, n, n, n, nidx);
      for (int d = 1; d < 4; ++d)
        db.setRange(d, xmin, xmax, n);
      portableFor(
          "fill 3D+idx", 0, n, 0, n, 0, n, 0, nidx,
          PORTABLE_LAMBDA(const int iz, const int iy, const int ix,
                          const int i) {
            db(iz, iy, ix, i) =
                (1 + i) * testFunction(g.x(iz), g.x(iy), g.x(ix));
          });
      runCase(
          "3D+idx", db.sizeBytes(), npoints, ntrials,
          PORTABLE_LAMBDA(const int i) {
            return db.interpToReal(p.z[i], p.y[i], p.x[i], p.idx[i]);
          },
          PORTABLE_LAMBDA(const int i) {
            return interpolate_with(db, linear(db.range(3), p.z[i]),
                                    linear(db.range(2), p.y[i]),
                                    linear(db.range(1), p.x[i]),
                                    at(p.idx[i]));
          });
      free(db);
    }

    freePoints(p);
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}
