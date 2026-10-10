#ifndef SPINER_INTERPOLATION_STENCILS
#define SPINER_INTERPOLATION_STENCILS
//======================================================================
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
//======================================================================

// Generative AI was used to assist with modifications to this file.

#include <ports-of-call/portability.hpp>
#include <ports-of-call/robust_utils.hpp>

#include "fast_nonuniform_grid_1d.hpp"
#include "nonuniform_grid_1d.hpp"
#include "piecewise_grid_1d.hpp"
#include "regular_grid_1d.hpp"

namespace Spiner {
namespace InterpolationHelpers {
// named tuples for value + gradient data. They work with structured
// binding.
template <typename T>
struct Result1D {
  T value, d_dx1;
};
template <typename T>
struct Result2D {
  T value, d_dx2, d_dx1;
};
} // namespace InterpolationHelpers

namespace InterpolationStencils {
struct Linear {
  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  interpToRealWithGrads(const T x1, const GridArray_t &g,
                        const Array_t &a) noexcept {
    int ix;
    weights_t<T> w, dw;
    g[0].weightsWithGrad(x1, ix, w, dw);
    const T f0 = a(ix);
    const T f1 = a(ix + 1);
    return {w[0] * f0 + w[1] * f1, dw[1] * (f1 - f0)};
  }
  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result2D<T>
  interpToRealWithGrads(const T x2, const T x1, const GridArray_t &g,
                        const Array_t &a) noexcept {
    int ix1, ix2;
    weights_t<T> w1, w2, dw1, dw2;
    g[0].weightsWithGrad(x1, ix1, w1, dw1);
    g[1].weightsWithGrad(x2, ix2, w2, dw2);
    const T f00 = a(ix2, ix1);
    const T f01 = a(ix2, ix1 + 1);
    const T f10 = a(ix2 + 1, ix1);
    const T f11 = a(ix2 + 1, ix1 + 1);
    return {w2[0] * (w1[0] * f00 + w1[1] * f01) +
                w2[1] * (w1[0] * f10 + w1[1] * f11),
            dw2[1] * (w1[0] * (f10 - f00) + w1[1] * (f11 - f01)),
            dw1[1] * (w2[0] * (f01 - f00) + w2[1] * (f11 - f10))};
  }
}; // struct Linear

// Rational function interpolation as described by Kerley, 1977, LA-6903-MS
struct Rational {
  constexpr int NCONTROL = 4;
  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  interpToRealWithGrads(const T x, const GridArray_t &g,
                        const Array_t &a) noexcept {
    const int N = a.GetDim(1);
    PORTABLE_REQUIRE(
        N >= NCONTROL,
        "Four point stencil required for rational function interpolation.");
    int ix = std::min(N - 3, std::max(1, g[0].index(x)));
    const std::array<T, NCONTROL> xg = {g[0].x(i - 1), g[0].x(i), g[0].x(i + 1),
                                        g[0].x(i + 2)};
    const std::array<T, NCONTROL> f = {a(ix - 1), a(ix), a(ix + 1), a(ix + 2)};
    return interp1DHelper_(x, xg, f);
  }

  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result2D<T>
  interpToRealWithGrads(const T x2, const T x1, const GridArray_t &g,
                        const Array_t &a) noexcept {
    using PortsOfCall::Robust::ratio;
    const int N1 = a.GetDim(1);
    const int N2 = a.GetDim(2);
    PORTABLE_REQUIRE(
        N1 >= NCONTROL,
        "Four point stencil required for rational function interpolation.");
    PORTABLE_REQUIRE(
        N2 >= NCONTROL,
        "Four point stencil required for rational function interpolation.");
    int ix1 = std::min(N1 - 3, std::max(1, g[0].index(x1)));
    int ix2 = std::min(N2 - 3, std::max(1, g[1].index(x2)));

    std::array<T, NCONTROL> x1g, x2g;
    for (int c = 0; c < NCONTROL; ++c) {
      x1g[c] = g[0].x(ix1 - 1 + c);
      x2g[c] = g[1].x(ix2 - 1 + c);
    }

    // Kerley eqn 17
    // compute rj(x) = r_{ix2}(x1)
    std::array<T, NCONTROL> f;
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2, ix1 - 1 + c);
    }
    auto [rjx, drjdx] = interp1DHelper_(x1, x1g, f);
    // compute r_{j + 1}(x) = r_{ix2 + 1}(x1)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2 + 1, ix1 - 1 + c);
    }
    auto [rjp1x, drjp1dx] = interp1DHelper_(x1, x1g, f);
    // compute ri(y) = r_{ix1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix1 - 1 + c, ix2);
    }
    auto [riy, dridy] = interp1DHelper_(x2, x2g, f);
    // compute r_{i+1}(y) = r_{ix1+1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix1 - 1 + c, ix2+1);
    }
    auto [rip1y, drip1dy] = interp1DHelper_(x2, x2g, f);

    // compute slopes for Kerley eqn 18
    T dx = x1g[ix1 + 1] - x1g[ix1];
    T dy = x2g[ix2 + 1] - x2g[ix2];
    T qx = ratio(x1 - x1g[ix1], dx);
    T qy = ratio(x2 - x2g[ix2], dy);
    // Kerley eqn 18
    T fout = rjx * (1 - qy) + rjp1x * qy + riy * (1 - qx) + rip1y * qx
      - a(ix2, ix1) * (1 - qx) * (1 - qy) - a(ix2 + 1, ix1) * (1 - qx) * qy
      - a(ix2, ix1 + 1) * qx * (1 - qy) - a(ix2 + 1, ix1 + 1) *qx * qy;

    T dfdx2 = ratio(dx * (a(ix2, ix1) - a(ix2 + 1, ix1))
                    - (a(ix2, ix1) - a(ix2 + 1, ix1) - a(ix2, ix1 + 1) + a(ix2 + 1, ix1 + 1)) * (x1 - x1g[ix1])
                    - dx * rjx + dx * rjp1x + dy * (dx - x + x1g[ix1]) * dridy + dy * (x1 - x1g[ix1]) * drip1dy,
                    dx * dy);
    T dfdx1 = ratio(dy * (a(ix2, ix1) - a(ix2, ix1 + 1))
                    - (a(ix2, ix1) - a(ix2 + 1, ix1) - a(ix2, ix1 + 1) + a(ix2 + 1, ix1 + 1)) * (x2 - x2g[ix2])
                    - dy * diy + dy * rip1y + dx * (dy - x2 + x2g[ix2]) * drjdx + dx * (x2 - x2g[ix2]) * drjp1dx,
                    dx * dy);

    return {fout, dfdx2, dfdx1};
  }

  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  interp1DHelper_(const T x, const std::array<T, NCONTROL> &xg,
                  const std::array<T, NCONTROL> &f) {
    using PortsOfCall::Robust::ratio;
    T fout, dfout, C1, C2;
    std::array<T, NCONTROL - 1> dx, S;
    for (int i = 0; i < NCONTROL - 1; ++i) {
      dx[i] = xg[i + 1] - xg[i]; // Kerley eqn 1
      S[i] = ratio(f[i + 1] - f[i], dx[i]);
    }

    if (x <= x[0]) { // linear extrap off bottom
      dfout = ratio(S[0], dx[0]);
      fout = f[0] + (x - x[0]) * dfout;
    } else if (x >= x[3]) { // linear extrap off top
      dfout = ratio(S[3], dx[3]);
      fout = f[3] + (x - x[3]) * dfout;
    }
    if (x < x[1]) { // Kerley eqn 12
      C2 = ratio(S[1] - S[0], dx[1] + dx[0]);
      if (S[0] * (S[0] - dx[0] * C2) <= 0) C2 = ratio(S[0], dx[0]);
      fout = f[0] + (x - x[0]) * (S[0] - C2 * (x[1] - x));
      dfout = S[0] + c2 * (2 * x - x[0] - x[1]);
    } else if (x >= x[2]) { // Kerley eqn 13
      C1 = ratio(S[2] - S[1], dx[2] + dx[1]);
      fout = f[2] + (x - x[2]) * (S[2] - C1 * (x[3] - x));
      dfout = S[2] + C1 * (2 * x - x[2] - x[3]);
    } else { // Kerley eqn 11
      C1 = ratio(S[1] - S[0], dx[1] + dx[0]);
      C2 = ratio(S[2] - S[1], dx[2] + dx[1]);
      if (S[0] * (S[0] - dx[0] * C2) <= 0) C1 = ratio(S[1] - 2 * S[0], dx[1]);

      T mu1 = std::abs(C2 * (x[2] - x));
      T mu2 = std::abs(C1 * (x - x[1]));

      fout = f[1] + (x - x[1]) * (S[1] - ratio(C1 * mu1 + C2 * mu2, mu1 + mu2) *
                                             (x[2] - x));
      dfout = S[1] + ratio((x[1] - x[2]) * mu1 * mu2 *
                               (x * (mu1 + mu2) - x[2] * mu1 - x[1] * mu2),
                           (x - x[1]) * (x - x[2]) * (mu1 + mu2) * (mu1 + mu2));
    }

    return {fout, dfout};
  }
};
} // namespace InterpolationStencils
} // namespace Spiner

#endif // SPINER_INTERPOLATION_STENCILS
