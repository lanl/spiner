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

#include <ports-of-call/portability.hpp>
#include <spiner/regular_grid_1d.hpp>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinRel;

// test transform: reverse(forward(xlo)) and reverse(forward(xhi)) don't map
// 100% correctly and end up inside the bounds of [xlo, xhi], so long as xlo
// and xhi bracket xref = 0.5
struct TxNarrow {
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T forward(const T x) {
    constexpr T xref = 0.5;
    constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
    return (static_cast<T>(1) - eps) * (x - xref) + xref;
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T reverse(const T u) {
    constexpr T uref = 0.5;
    constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
    return (static_cast<T>(1) - eps) * (u - uref) + uref;
  }
};

// test transform: reverse(forward(xlo)) and reverse(forward(xhi)) don't map
// 100% correctly and end up outside the bounds of [xlo, xhi], so long as xlo
// and xhi bracket xref = 0.5
struct TxExpand {
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T forward(const T x) {
    constexpr T xref = 0.5;
    constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
    return (static_cast<T>(1) + eps) * (x - xref) + xref;
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T reverse(const T u) {
    constexpr T uref = 0.5;
    constexpr T eps = 10 * std::numeric_limits<T>::epsilon();
    return (static_cast<T>(1) + eps) * (u - uref) + uref;
  }
};

// test transform: monotonically decreasing instead of increasing
struct TxDecrease {
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T forward(const T x) {
    return 1 - x;
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T reverse(const T u) {
    return 1 - u;
  }
};

TEST_CASE("RegularGrid1D", "[RegularGrid1D]") {
  SECTION("A regular grid 1d emits appropriate metadata") {
    constexpr Real min = -1;
    constexpr Real max = 1;
    constexpr size_t N = 10;
    Spiner::RegularGrid1D<Real> g(min, max, N);
    REQUIRE(g.min() == min);
    REQUIRE(g.max() == max);
    REQUIRE(g.nPoints() == N);
  }
}

TEST_CASE("RegularGrid1D with production transformations",
          "[RegularGrid1D][transformations]") {
  constexpr double min = 1;
  constexpr double max = 1024;
  constexpr size_t N = 11;
  Spiner::RegularGrid1D<double, Spiner::TransformLinear> glin(min, max, N);
  Spiner::RegularGrid1D<double, Spiner::TransformLogarithmic> glog(min, max, N);
  REQUIRE(glin.min() == min);
  REQUIRE(glog.min() == min);
  REQUIRE(glin.max() == max);
  REQUIRE(glog.max() == max);
  REQUIRE(glin.nPoints() == N);
  REQUIRE(glog.nPoints() == N);

  // Check all fixed points (lin)
  for (std::size_t n = 0; n < N; ++n) {
    const double xx = 102.3 * double(n) + 1;
    CHECK_THAT(glin.x(n), WithinRel(xx, 1.0e-12));
  }

  // Check all fixed points (log)
  for (std::size_t n = 0; n < N; ++n) {
    const double xx = std::pow(double(2), n);
    CHECK_THAT(glog.x(n), WithinRel(xx, 1.0e-12));
  }

  // Check weights (lin)
  for (std::size_t n = 1; n < N; ++n) {
    const double xlin = 0.5 * (glin.x(n - 1) + glin.x(n));
    int ilin;
    Spiner::weights_t<double> wlin;
    glin.weights(xlin, ilin, wlin);
    CHECK(ilin == n - 1);
    CHECK_THAT(wlin.first, WithinRel(0.5));
    CHECK_THAT(wlin.second, WithinRel(0.5));
  }

  // Check weights (log)
  for (std::size_t n = 1; n < N; ++n) {
    const double xlog = std::sqrt(glog.x(n - 1) * glog.x(n));
    int ilog;
    Spiner::weights_t<double> wlog;
    glog.weights(xlog, ilog, wlog);
    CHECK(ilog == n - 1);
    CHECK_THAT(wlog.first, WithinRel(0.5));
    CHECK_THAT(wlog.second, WithinRel(0.5));
  }
}

TEST_CASE("RegularGrid1D with test transformations",
          "[RegularGrid1D][transformations]") {
  constexpr double min = 0.0;
  constexpr double max = 1.0;
  constexpr size_t N = 101;
  Spiner::RegularGrid1D<double, TxNarrow> gn(min, max, N);
  Spiner::RegularGrid1D<double, TxExpand> ge(min, max, N);
  Spiner::RegularGrid1D<double, TxDecrease> gd(min, max, N);

  // Basic tests (narrow)
  REQUIRE(gn.min() == min);
  REQUIRE(gn.max() == max);
  REQUIRE(gn.nPoints() == N);
  // TODO: Do we want the bounds to be exact?
  CHECK(gn.x(0) == min);
  CHECK(gn.x(N - 1) == max);

  // Basic tests (expand)
  REQUIRE(ge.min() == min);
  REQUIRE(ge.max() == max);
  REQUIRE(ge.nPoints() == N);
  // TODO: Do we want the bounds to be exact?
  CHECK(ge.x(0) == min);
  CHECK(ge.x(N - 1) == max);

  // Basic tests (decrease)
  REQUIRE(gd.min() == min);
  REQUIRE(gd.max() == max);
  REQUIRE(gd.nPoints() == N);
  // TODO: Do we want the bounds to be exact?
  CHECK(gd.x(0) == min);
  CHECK(gd.x(N - 1) == max);

  // Check all fixed points (narrow)
  for (std::size_t n = 0; n < N; ++n) {
    const double xx = 0.01 * double(n);
    CHECK_THAT(gn.x(n), WithinRel(xx, 1.0e-12));
  }

  // Check all fixed points (expand)
  for (std::size_t n = 0; n < N; ++n) {
    const double xx = 0.01 * double(n);
    CHECK_THAT(ge.x(n), WithinRel(xx, 1.0e-12));
  }

  // Check all fixed points (decrease)
  for (std::size_t n = 0; n < N; ++n) {
    const double xx = 0.01 * double(n);
    CHECK_THAT(gd.x(n), WithinRel(xx, 1.0e-12));
  }

  // Check weights (narrow)
  for (std::size_t n = 1; n < N; ++n) {
    const double x = 0.5 * (gn.x(n - 1) + gn.x(n));
    int index;
    Spiner::weights_t<double> weights;
    gn.weights(x, index, weights);
    CHECK(index == n - 1);
    CHECK_THAT(weights.first, WithinRel(0.5, 1.0e-12));
    CHECK_THAT(weights.second, WithinRel(0.5, 1.0e-12));
  }

  // Check weights (expand)
  for (std::size_t n = 1; n < N; ++n) {
    const double x = 0.5 * (ge.x(n - 1) + ge.x(n));
    int index;
    Spiner::weights_t<double> weights;
    ge.weights(x, index, weights);
    CHECK(index == n - 1);
    CHECK_THAT(weights.first, WithinRel(0.5, 1.0e-12));
    CHECK_THAT(weights.second, WithinRel(0.5, 1.0e-12));
  }

  // Check weights (decrease)
  for (std::size_t n = 1; n < N; ++n) {
    const double x = 0.5 * (gd.x(n - 1) + gd.x(n));
    int index;
    Spiner::weights_t<double> weights;
    gd.weights(x, index, weights);
    CHECK(index == n - 1);
    CHECK_THAT(weights.first, WithinRel(0.5, 1.0e-12));
    CHECK_THAT(weights.second, WithinRel(0.5, 1.0e-12));
  }
}
