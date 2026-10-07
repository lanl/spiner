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
    using PortsOfCall::Robust::ratio;
    T fout, dfout, C1, C2;

    const int N = a.GetDim(1);
    PORTABLE_REQUIRE(
        N >= NCONTROL,
        "Four point stencil required for rational function interpolation.");
    int ix = std::min(N - 3, std::max(1, g[0].index(x)));
    const std::array<T, NCONTROL> xg = {g[0].x(i - 1), g[0].x(i), g[0].x(i + 1),
                                        g[0].x(i + 2)};
    const std::array<T, NCONTROL> f = {a(ix - 1), a(ix), a(ix + 1), a(ix + 2)};
    std::array<T, NCONTROL - 1> dx, S;
    for (int i = 0; i < NCONTROL - 1; ++i) {
      dx[i] = xg[i + 1] - xg[i]; // Kerley eqn 1
      S[i] = ratio(f[i + 1] - f[i], dx[i]);
    }

    if (x <= x[0]) { // linear extrap off bottom
      dfout = ratio(S[0], dx[0]);
      fout = f[0]  + (x - x[0]) * dfout;
    } else if (x >= x[3]) { // linear extrap off top
      dfout = ratio(S[3], dx[3]);
      fout = f[3] + (x - x[3]) * dfout;
    } if (x < x[1]) { // Kerley eqn 12
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
