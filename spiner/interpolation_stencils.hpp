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

#include <algorithm>
#include <array>
#include <cmath>

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

// Rational function interpolation as described by Kerley, 1977, LA-6903-MS.
// Requires at least four points per interpolated axis. Outside the table, use
// linear secant extrapolation; at exact endpoints use the interior derivative.
struct Rational {
  static constexpr int NCONTROL = 4;
  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  interpToRealWithGrads(const T x, const GridArray_t &g,
                        const Array_t &a) noexcept {
    const int N = a.GetDim(1);
    PORTABLE_REQUIRE(
        N >= NCONTROL,
        "Four point stencil required for rational function interpolation.");
    int ix = std::min(N - 3, std::max(1, g[0].index(x)));
    const std::array<T, NCONTROL> xg = {g[0].x(ix - 1), g[0].x(ix),
                                        g[0].x(ix + 1), g[0].x(ix + 2)};
    const std::array<T, NCONTROL> f = {a(ix - 1), a(ix), a(ix + 1), a(ix + 2)};
    return interp1DHelper_(x, xg, f, ix == 1);
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
    // Indices of cell enclosing x1 and x2
    const int ix1 = std::min(N1 - 2, std::max(0, g[0].index(x1)));
    const int ix2 = std::min(N2 - 2, std::max(0, g[1].index(x2)));
    // indices of the first cell in the full four-point stencil
    const int start1 = std::min(N1 - 4, std::max(0, ix1 - 1));
    const int start2 = std::min(N2 - 4, std::max(0, ix2 - 1));

    std::array<T, NCONTROL> x1g, x2g;
    for (int c = 0; c < NCONTROL; ++c) {
      x1g[c] = g[0].x(start1 + c);
      x2g[c] = g[1].x(start2 + c);
    }

    // Kerley eqn 17
    // compute rj(x) = r_{ix2}(x1)
    std::array<T, NCONTROL> f;
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2, start1 + c);
    }
    auto [rjx, drjdx] = interp1DHelper_(x1, x1g, f, ix1 == 1);
    // compute r_{j + 1}(x) = r_{ix2 + 1}(x1)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2 + 1, start1 + c);
    }
    auto [rjp1x, drjp1dx] = interp1DHelper_(x1, x1g, f, ix1 == 1);
    // compute ri(y) = r_{ix1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(start2 + c, ix1);
    }
    auto [riy, dridy] = interp1DHelper_(x2, x2g, f, ix2 == 1);
    // compute r_{i+1}(y) = r_{ix1+1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(start2 + c, ix1 + 1);
    }
    auto [rip1y, drip1dy] = interp1DHelper_(x2, x2g, f, ix2 == 1);

    const T dx = g[0].x(ix1 + 1) - g[0].x(ix1);
    const T dy = g[1].x(ix2 + 1) - g[1].x(ix2);
    const T qx = ratio(x1 - g[0].x(ix1), dx);
    const T qy = ratio(x2 - g[1].x(ix2), dy);
    const T f00 = a(ix2, ix1), f10 = a(ix2, ix1 + 1);
    const T f01 = a(ix2 + 1, ix1), f11 = a(ix2 + 1, ix1 + 1);
    // Kerley eqn 18 and its partial derivatives before factoring:
    // clang-format off
    // const T fout = rjx * (1 - qy) + rjp1x * qy + riy * (1 - qx) + rip1y * qx -
    //                f00 * (1 - qx) * (1 - qy) - f01 * (1 - qx) * qy -
    //                f10 * qx * (1 - qy) - f11 * qx * qy;
    // const T dfdx1 =
    //     (1 - qy) * drjdx + qy * drjp1dx +
    //     ratio(rip1y - riy - (1 - qy) * (f10 - f00) - qy * (f11 - f01), dx);
    // const T dfdx2 =
    //     (1 - qx) * dridy + qx * drip1dy +
    //     ratio(rjp1x - rjx - (1 - qx) * (f01 - f00) - qx * (f11 - f10), dy);
    // clang-format on
    // Share the bilinear corner differences and edge corrections between
    // the value and both gradients. m is the mixed corner difference.
    const T p = f10 - f00;
    const T q = f01 - f00;
    const T m = (f11 - f01) - p;
    const T cx = rip1y - riy - p;
    const T cy = rjp1x - rjx - q;
    const T tx = cx - qy * m;
    const T ty = cy - qx * m;
    const T fout = rjx + riy - f00 + qx * cx + qy * ty;
    const T dfdx1 = drjdx + qy * (drjp1dx - drjdx) + ratio(tx, dx);
    const T dfdx2 = dridy + qx * (drip1dy - dridy) + ratio(ty, dy);

    return {fout, dfdx2, dfdx1};
  }

  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  interp1DHelper_(const T x, const std::array<T, NCONTROL> &xg,
                  const std::array<T, NCONTROL> &f,
                  const bool second_table_interval) noexcept {
    using PortsOfCall::Robust::ratio;
    T fout, dfout, C1, C2;
    std::array<T, NCONTROL - 1> dx, S;
    for (int i = 0; i < NCONTROL - 1; ++i) {
      dx[i] = xg[i + 1] - xg[i]; // Kerley eqn 1
      S[i] = ratio(f[i + 1] - f[i], dx[i]);
    }

    // Linearly extrapolate strictly outside the table.
    // Exact endpoints retain the interior derivative (a jump is possible).
    if (x < xg[0]) return {f[0] + (x - xg[0]) * S[0], S[0]};
    if (x > xg[3]) return {f[3] + (x - xg[3]) * S[2], S[2]};
    if (x < xg[1]) { // Kerley eqn 12
      C2 = ratio(S[1] - S[0], dx[1] + dx[0]);
      if (S[0] * (S[0] - dx[0] * C2) <= 0) C2 = ratio(S[0], dx[0]);
      // Original quadratic value and derivative:
      // fout = f[0] + (x - xg[0]) * (S[0] - C2 * (xg[1] - x));
      // dfout = S[0] + C2 * (2 * x - xg[0] - xg[1]);
      const T Q = x - xg[0], R = xg[1] - x;
      const T L = S[0] - R * C2;
      fout = f[0] + Q * L;
      dfout = L + Q * C2;
    } else if (x >= xg[2]) { // Kerley eqn 13
      C1 = ratio(S[2] - S[1], dx[2] + dx[1]);
      // Original quadratic value and derivative:
      // fout = f[2] + (x - xg[2]) * (S[2] - C1 * (xg[3] - x));
      // dfout = S[2] + C1 * (2 * x - xg[2] - xg[3]);
      const T Q = x - xg[2], R = xg[3] - x;
      const T L = S[2] - R * C1;
      fout = f[2] + Q * L;
      dfout = L + Q * C1;
    } else { // Kerley eqn 11
      C1 = ratio(S[1] - S[0], dx[1] + dx[0]);
      C2 = ratio(S[2] - S[1], dx[2] + dx[1]);
      // Appendix A, lines 50–51: propagate the lower endpoint constraint
      // only into the second interval of the full table.
      if (second_table_interval && S[0] * (S[0] - dx[0] * C1) <= 0)
        C1 = ratio(S[1] - 2 * S[0], dx[1]);

      // Appendix A, lines 54–59, including its exact-knot convention.
      const T Q = x - xg[1], R = xg[2] - x;
      const T mu1 = std::abs(C2 * R), mu2 = std::abs(C1 * Q);
      const T sum = mu1 + mu2;
      const T w = sum > 0 ? ratio(mu1, sum) : T(0);
      const T C = C2 + w * (C1 - C2);
      // Original Appendix A value and derivative:
      // fout = f[1] + Q * (S[1] - R * C);
      // dfout = S[1] + (Q - R) * C + dx[1] * (C - C2) * (1 - w);
      // Reuse the value's inner expression in the derivative.
      const T L = S[1] - R * C;
      fout = f[1] + Q * L;
      dfout = L + Q * C + dx[1] * (C - C2) * (1 - w);
    }

    return {fout, dfout};
  }
};
} // namespace InterpolationStencils
} // namespace Spiner

#endif // SPINER_INTERPOLATION_STENCILS
