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

TEST_CASE("RegularGrid1D with transformations", "[RegularGrid1D]") {
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
    for(std::size_t n = 0; n < N; ++n) {
        const double xx = 102.3 * double(n) + 1;
        const double uu = xx;
        CHECK_THAT(glin.u(n), WithinRel(uu, 1.0e-12));
        CHECK_THAT(glin.x(n), WithinRel(xx, 1.0e-12));
    }

    // Check all fixed points (log)
    for(std::size_t n = 0; n < N; ++n) {
        const double xx = std::pow(double(2), n);
        const double uu = std::log(xx);
        CHECK_THAT(glog.u(n), WithinRel(uu, 1.0e-12));
        CHECK_THAT(glog.x(n), WithinRel(xx, 1.0e-12));
    }

    // Check weights (lin)
    for(std::size_t n = 1; n < N; ++n) {
        const double xlin = 0.5 * (glin.x(n-1) + glin.x(n));
        int ilin;
        Spiner::weights_t<double> wlin;
        glin.weights(xlin, ilin, wlin);
        CHECK(ilin == n - 1);
        CHECK_THAT(wlin.first, WithinRel(0.5));
        CHECK_THAT(wlin.second, WithinRel(0.5));
    }

    // Check weights (log)
    for(std::size_t n = 1; n < N; ++n) {
        const double xlog = std::sqrt(glog.x(n-1) * glog.x(n));
        int ilog;
        Spiner::weights_t<double> wlog;
        glog.weights(xlog, ilog, wlog);
        CHECK(ilog == n - 1);
        CHECK_THAT(wlog.first, WithinRel(0.5));
        CHECK_THAT(wlog.second, WithinRel(0.5));
    }

}

