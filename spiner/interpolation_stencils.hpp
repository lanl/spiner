#ifndef SPINER_INTERPOLATION_STENCILS_
#define SPINER_INTERPOLATION_STENCILS_
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
    // JMM: We don't use a deduction guide because the CUDA compiler
    // doesn't accept private scope deduction guides.
    const AxisGeometry<T> geometry(x, g[0], a.GetDim(1));
    const int start = geometry.start;
    const std::array<T, NCONTROL> f = {a(start), a(start + 1), a(start + 2),
                                       a(start + 3)};
    return evaluateLine_(geometry, f);
  }

  template <typename T, typename GridArray_t, typename Array_t>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result2D<T>
  interpToRealWithGrads(const T x2, const T x1, const GridArray_t &g,
                        const Array_t &a) noexcept {
    // JMM: We don't use a deduction guide because the CUDA compiler
    // doesn't accept private scope deduction guides.
    const AxisGeometry<T> gx(x1, g[0], a.GetDim(1));
    const AxisGeometry<T> gy(x2, g[1], a.GetDim(2));
    const int ix1 = gx.cell, ix2 = gy.cell;
    const int start1 = gx.start, start2 = gy.start;

    /*
      JMM: Note the notation here gets tricky because I followed
      Kerley's basic convetion and the 2nd-order interpolation needs a
      way of representing the "continuum" 1D interpolation objects
      based at grid lines. This doesn't match the rest of
      Spiner. Convention is:

      - r is the 1D interpolant.

      - i or j indicates index in x1 or x2 respectively for location
        on mesh the interpolant is based. It then interpolates in the
        ORTHOGNAL direction as a continuum object.

      - x or y specifies the continuum interpolant.

      E.g., riy specifies interpolant continuous in x2 based at point in x1 i
      rjp1x specifies interpolant continuous in x1 based at point j+1 in x2.

      This gets more confusing below when we actually compute the
      final formulae, as I (i.e., Claude) introduced some intermediate
      values that don't really have useful/meaningful names. I just
      picked something. Sorry.
    */

    // Kerley eqn 17
    // compute rj(x) = r_{ix2}(x1)
    std::array<T, NCONTROL> f;
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2, start1 + c);
    }
    auto [rjx, drjdx] = evaluateLine_(gx, f);
    // compute r_{j + 1}(x) = r_{ix2 + 1}(x1)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(ix2 + 1, start1 + c);
    }
    auto [rjp1x, drjp1dx] = evaluateLine_(gx, f);
    // compute ri(y) = r_{ix1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(start2 + c, ix1);
    }
    auto [riy, dridy] = evaluateLine_(gy, f);
    // compute r_{i+1}(y) = r_{ix1+1}(x2)
    for (int c = 0; c < NCONTROL; ++c) {
      f[c] = a(start2 + c, ix1 + 1);
    }
    auto [rip1y, drip1dy] = evaluateLine_(gy, f);

    const T inverse_dx = gx.inverse_cell_dx;
    const T inverse_dy = gy.inverse_cell_dx;
    // Originally qx = ratio(gx.cell_offset, dx), likewise for qy.
    const T qx = gx.cell_offset * inverse_dx; // x-slope denominator
    const T qy = gy.cell_offset * inverse_dy; // y-slope denominator
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
    const T p = f10 - f00;        // fip1j - fij, x-slope numerator at i,j
    const T q = f01 - f00;        // fijp1 - fij; y-slope numerator at i,j
    const T m = (f11 - f01) - p;  // x-slope numerator at i, j+1 minus p
    const T cx = rip1y - riy - p; // r(y) at i+1 - r(i) at i - p.
    const T cy = rjp1x - rjx - q; // r(x) at j+1 - r(x) at j - p.
    const T tx = cx - qy * m;     // build "cross terms" of interpolation ...
    const T ty = cy - qx * m;
    // fout
    rjx += riy - f00 + qx * cx + qy * ty;
    // Reuse the guarded reciprocals for the original ratio(tx, dx) and
    // ratio(ty, dy) gradient corrections.
    // dfdx1
    drjdx += qy * (drjp1dx - drjdx) + tx * inverse_dx;
    // dfdx2
    dridy += qx * (drip1dy - dridy) + ty * inverse_dy;

    return {rjx, dridy, drjdx};
  }

 private:
  enum class AxisMode {
    Interior,
    LowerQuadratic,
    UpperQuadratic,
    LowerExtrapolation,
    UpperExtrapolation
  };

  // Per-query scratch only: no allocations or ownership of grid memory.
  template <typename T>
  struct AxisGeometry {
    int cell, start;
    std::array<T, NCONTROL - 1> dx;
    std::array<T, NCONTROL - 1> inverse_dx;
    std::array<T, 2> inv_curv_width;
    T Q, R, cell_offset, inverse_cell_dx;
    AxisMode mode;
    bool second_table_interval;

    template <typename Grid_t>
    PORTABLE_FORCEINLINE_FUNCTION AxisGeometry(const T x, const Grid_t &grid,
                                               const int N) noexcept {
      PORTABLE_REQUIRE(
          N >= NCONTROL,
          "Four point stencil required for rational function interpolation.");
      cell = std::min(N - 2, std::max(0, grid.index(x)));
      start = std::min(N - 4, std::max(0, cell - 1));
      second_table_interval = cell == 1;
      std::array<T, NCONTROL> xg;
      for (int c = 0; c < NCONTROL; ++c)
        xg[c] = grid.x(start + c);
      for (int c = 0; c < NCONTROL - 1; ++c) {
        // Kerley eqn 1.
        dx[c] = xg[c + 1] - xg[c];
        inverse_dx[c] = grid.idx(start + c);
      }
      inv_curv_width[0] = grid.invCurvatureWidth(start);
      inv_curv_width[1] = grid.invCurvatureWidth(start + 1);
      const int k = cell - start;
      cell_offset = x - (k == 0 ? xg[0] : (k == 1 ? xg[1] : xg[2]));
      inverse_cell_dx =
          k == 0 ? inverse_dx[0] : (k == 1 ? inverse_dx[1] : inverse_dx[2]);

      // Q is anchored to the endpoint for extrapolation, to the interval's
      // left point otherwise.
      R = T(0);
      if (x < xg[0]) {
        mode = AxisMode::LowerExtrapolation;
        Q = x - xg[0];
      } else if (x > xg[3]) {
        mode = AxisMode::UpperExtrapolation;
        Q = x - xg[3];
      } else if (x < xg[1]) {
        mode = AxisMode::LowerQuadratic;
        Q = x - xg[0];
        R = xg[1] - x;
      } else if (x >= xg[2]) {
        mode = AxisMode::UpperQuadratic;
        Q = x - xg[2];
        R = xg[3] - x;
      } else {
        mode = AxisMode::Interior;
        Q = x - xg[1];
        R = xg[2] - x;
      }
    }
  };

  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static InterpolationHelpers::Result1D<T>
  evaluateLine_(const AxisGeometry<T> &geometry,
                const std::array<T, NCONTROL> &f) noexcept {
    const auto &dx = geometry.dx;
    const auto &inverse_dx = geometry.inverse_dx;
    const T Q = geometry.Q, R = geometry.R;
    T fout, dfout, C1, C2;
    std::array<T, NCONTROL - 1> S;
    for (int i = 0; i < NCONTROL - 1; ++i)
      S[i] = (f[i + 1] - f[i]) * inverse_dx[i];

    // Linearly extrapolate strictly outside the table.
    // Exact endpoints retain the interior derivative (a jump is possible).
    if (geometry.mode == AxisMode::LowerExtrapolation)
      return {f[0] + Q * S[0], S[0]};
    if (geometry.mode == AxisMode::UpperExtrapolation)
      return {f[3] + Q * S[2], S[2]};
    if (geometry.mode == AxisMode::LowerQuadratic) { // Kerley eqn 12
      C2 = (S[1] - S[0]) * geometry.inv_curv_width[0];
      // Original limited curvature: C2 = ratio(S[0], dx[0]).
      if (S[0] * (S[0] - dx[0] * C2) <= 0) C2 = S[0] * inverse_dx[0];
      // Original quadratic value and derivative:
      // fout = f[0] + (x - xg[0]) * (S[0] - C2 * (xg[1] - x));
      // dfout = S[0] + C2 * (2 * x - xg[0] - xg[1]);
      const T L = S[0] - R * C2;
      fout = f[0] + Q * L;
      dfout = L + Q * C2;
    } else if (geometry.mode == AxisMode::UpperQuadratic) { // Kerley eqn 13
      C1 = (S[2] - S[1]) * geometry.inv_curv_width[1];
      // Original quadratic value and derivative:
      // fout = f[2] + (x - xg[2]) * (S[2] - C1 * (xg[3] - x));
      // dfout = S[2] + C1 * (2 * x - xg[2] - xg[3]);
      const T L = S[2] - R * C1;
      fout = f[2] + Q * L;
      dfout = L + Q * C1;
    } else { // Kerley eqn 11
      C1 = (S[1] - S[0]) * geometry.inv_curv_width[0];
      C2 = (S[2] - S[1]) * geometry.inv_curv_width[1];
      // Appendix A, lines 50–51: propagate the lower endpoint constraint
      // only into the second interval of the full table.
      // Original limited curvature: C1 = ratio(S[1] - 2 * S[0], dx[1]).
      if (geometry.second_table_interval && S[0] * (S[0] - dx[0] * C1) <= 0)
        C1 = (S[1] - 2 * S[0]) * inverse_dx[1];

      // Appendix A, lines 54–59, including its exact-knot convention.
      const T mu1 = std::abs(C2 * R), mu2 = std::abs(C1 * Q);
      const T sum = mu1 + mu2;
      const T w = sum > 0 ? mu1 / sum : T(0);
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

#endif // SPINER_INTERPOLATION_STENCILS_
