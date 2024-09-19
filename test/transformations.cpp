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

#include <spiner/transformations.hpp>

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

namespace {
    PORTABLE_FUNCTION double relerr(const double z, const double zz) {
        if (z == 0) {
            return std::abs(zz - z);
        } else {
            return std::abs(double(1) - zz/z);
        }
    }
};

TEST_CASE("transformation: linear", "") {
    // This one is almost too simple to meaningfully test, but we can at least
    // ensure that it compiles and does the trivially-obvious things.
    using Transform = Spiner::TransformLinear;

    // Test on CPU
    const double x0 = 3.14159;
    CHECK(Transform::forward(x0) == x0);
    CHECK(Transform::reverse(x0) == x0);
    const float x1 = 2.71828;
    CHECK(Transform::forward(x1) == x1);
    CHECK(Transform::reverse(x1) == x1);

    // Test on GPU (or CPU if no GPU available)
    const int num_threads = 10;
    double accum = 0;
    portableReduce(
        "parallel_reduce", 0, num_threads,
        PORTABLE_LAMBDA(const int n, double &reduce)
        {
            const double x = static_cast<double>(n);
            const double y = x;
            const double yy = Transform::forward(x);
            const double xx = Transform::reverse(y);
            reduce += 0.5 * (relerr(x, xx) + relerr(y, yy));
        },
        accum);
    CHECK(accum / num_threads == 0);
}

TEST_CASE("transformation: logarithmic", "") {
    using Transform = Spiner::TransformLogarithmic;

    // Test on CPU
    auto matcher = [](const double d) {
        //return WithinULP(d, 500);
        return WithinRel(d, 6.0e-12);
    };
    // Scan across most of the range
    for (int p = -300; p <= 300; ++p) {
        const double x = std::pow(double(10), p);
        const double y = std::log(x);
        // Basic checks
        CHECK_THAT(Transform::forward(x), matcher(y));
        CHECK_THAT(Transform::reverse(y), matcher(x));
        // Round trip
        CHECK_THAT(Transform::reverse(Transform::forward(x)), matcher(x));
    }
    // Special value
    CHECK(std::isfinite(Transform::forward(0)));
    CHECK(Transform::reverse(Transform::forward(0)) == 0);

    // Test on GPU (or CPU if no GPU available)
    const int num_threads = 101;
    double accum = 0;
    portableReduce(
        "parallel_reduce", 0, num_threads,
        PORTABLE_LAMBDA(const int n, double &reduce)
        {
            const int p = n - num_threads / 2;
            const double x = std::pow(double(10.0), p);
            const double y = std::log(x);
            const double yy = Transform::forward(x);
            const double xx = Transform::reverse(y);
            reduce += 0.5 * (relerr(x, xx) + relerr(y, yy));
        },
        accum);
    CHECK(accum / num_threads <= 6.0e-12);
}

