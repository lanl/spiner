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
    return std::abs(double(1) - zz / z);
  }
}
}; // namespace

TEST_CASE("transformation: linear", "[transformations]") {
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
      PORTABLE_LAMBDA(const int n, double &reduce) {
        const double x = static_cast<double>(n);
        const double y = x;
        const double yy = Transform::forward(x);
        const double xx = Transform::reverse(y);
        reduce += 0.5 * (relerr(x, xx) + relerr(y, yy));
      },
      accum);
  CHECK(accum / num_threads == 0);
}

TEST_CASE("transformation: logarithmic", "[transformations]") {
  using Transform = Spiner::TransformLogarithmic;

  // Test on CPU
  auto matcher = [](const double d) {
    // return WithinULP(d, 500);
    return WithinRel(d, 1.0e-12);
  };
  // Scan across most of the range
  for (int p = -307; p <= 307; ++p) {
    const double x = std::pow(double(10), p);
    const double y = std::log(x);
    // Basic checks
    // -- The "exact" calculation (y = log(x)) won't match the transformation
    //    (y = log(x + epsilon) for very small values of x, so skip checks.
    if (p >= -295) {
      CHECK_THAT(Transform::forward(x), matcher(y));
      CHECK_THAT(Transform::reverse(y), matcher(x));
    }
    // Round trip
    CHECK_THAT(Transform::reverse(Transform::forward(x)), matcher(x));
  }
  // Special value
  // -- Transform::forward(0) will infer the type to be an integer, and you
  //    will get a WILDLY incorrect answer for the round trip.
  CHECK(std::isfinite(Transform::forward(0.0)));
  CHECK(Transform::reverse(Transform::forward(static_cast<double>(0))) == 0);
  CHECK(Transform::reverse(Transform::forward(static_cast<float>(0))) == 0);

  // Test on GPU (or CPU if no GPU available)
  const int num_threads = 101;
  double accum = 0;
  portableReduce(
      "parallel_reduce", 0, num_threads,
      PORTABLE_LAMBDA(const int n, double &reduce) {
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
