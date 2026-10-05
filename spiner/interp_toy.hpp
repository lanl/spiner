#ifndef _SPINER_INTERP_TOY_HPP_
#define _SPINER_INTERP_TOY_HPP_
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

// Generic N-dimensional interpolation built from one rule per dimension.
//
// A rule is any object with a member
//
//   template <class Sample>
//   PORTABLE_INLINE_FUNCTION auto evaluate(const Sample &sample) const;
//
// where sample(i) returns the (already reduced) value of the dimensions
// to the right of this one, with this dimension fixed at index i. A
// rule may call sample any number of times at any valid indices and
// combine the results however it likes, including nonlinearly (e.g.,
// limiters or log-space interpolation).
//
// Example, for a rank-3 DataBox db and a query point (z, y, x), with
// the arguments ordered slowest to fastest as in db(iz, iy, ix):
//
//   using namespace Spiner::interp;
//   Real v = interpolate_with(db, linear(db.range(2), z),
//                             linear(db.range(1), y),
//                             linear(db.range(0), x));
//
// NOTE: Rules are ordered like operator(), slowest to fastest. DataBox
// grids are numbered the other way: range(0) is the fastest dimension.
// For a rank-N box, rule D pairs with db.range(N - 1 - D).

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_errors.hpp"
#include "regular_grid_1d.hpp"

namespace Spiner {
namespace interp {

// Fix one dimension at a particular array index.
struct FixedIndex {
  int index;

  template <class Sample>
  PORTABLE_INLINE_FUNCTION auto evaluate(const Sample &sample) const {
    return sample(index);
  }
};

PORTABLE_INLINE_FUNCTION constexpr FixedIndex at(const int index) {
  return {index};
}

// Interpolate one dimension using M samples.
//
// Reducer signature:
//   R(const std::array<V, M> &values)
// where V is the value type returned by the samples, i.e., the
// array's value type for the rightmost dimension, or the result type
// of the rule to the right otherwise. A generic lambda taking
// `const auto &values` works for any V.
//
// The query coordinate, grid positions, and any precomputed weights or
// limiter parameters the reducer needs are held in its captures. The
// reducer is called on whichever execution space evaluate runs on, so
// it must be device-callable when used in a kernel. A reducer created
// on the host and captured into a kernel must be a PORTABLE_LAMBDA.
//
// The rule's builder is responsible for keeping every entry of indices
// in bounds, e.g., by shifting a wide stencil inward near the edges of
// the grid.
template <std::size_t M, class F>
struct Interpolation {
  static_assert(M > 0, "An interpolation needs at least one sample");

  std::array<int, M> indices;
  F reduce;

  template <class Sample>
  PORTABLE_INLINE_FUNCTION auto evaluate(const Sample &sample) const {
    using V = std::remove_cvref_t<decltype(sample(indices[0]))>;
    std::array<V, M> values;
    for (std::size_t j = 0; j < M; ++j) {
      values[j] = sample(indices[j]);
    }
    return reduce(values);
  }
};

// Prepare an arbitrary interpolation rule.
template <std::size_t M, class F>
PORTABLE_INLINE_FUNCTION auto
make_interpolation(const std::array<int, M> &indices, F reduce) {
  return Interpolation<M, F>{indices, std::move(reduce)};
}

// Prepare linear interpolation on any Spiner grid, computing the
// weights only once. Query points outside the grid are handled
// according to the grid's weights() method.
template <class Grid>
PORTABLE_INLINE_FUNCTION auto linear(const Grid &grid,
                                     const typename Grid::ValueType x) {
  using T = typename Grid::ValueType;
  int ix;
  weights_t<T> w;
  grid.weights(x, ix, w);
  const T w0 = w.first;
  const T w1 = w.second;
  return make_interpolation(
      std::array<int, 2>{ix, ix + 1},
      [=](const auto &values) { return w0 * values[0] + w1 * values[1]; });
}

namespace detail {

// Debug-only check that sample index i of rule D is in bounds. Only
// applies to arrays that report their shape through DataBox-style
// rank() and dim(), where dim(1) is the fastest dimension.
template <int D, int N, class Array>
PORTABLE_INLINE_FUNCTION void require_index_in_bounds(const Array &array,
                                                      const int i) {
  if constexpr (requires { array.dim(1); }) {
    PORTABLE_REQUIRE(0 <= i && i < array.dim(N - D),
                     "Interpolation sample index out of bounds");
  }
}

// Evaluate rule D. access(i_D, ..., i_{N-1}) reads array at the
// indices already chosen by rules 0, ..., D-1 followed by its
// arguments. Each level fixes one more leading index, so no index
// storage or tuple machinery is needed.
template <int D, int N, class Array, class Access, class Rule, class... Rest>
PORTABLE_INLINE_FUNCTION auto evaluate(const Array &array, const Access &access,
                                       const Rule &rule, const Rest &...rest) {
  return rule.evaluate([&](const int i) {
    require_index_in_bounds<D, N>(array, i);
    if constexpr (sizeof...(Rest) == 0) {
      return access(i);
    } else {
      return evaluate<D + 1, N>(
          array, [&](const auto... js) { return access(i, js...); }, rest...);
    }
  });
}

} // namespace detail

// One rule per array dimension, ordered slowest to fastest (left to
// right, as in array(i, j, k, ...)). The rightmost dimension is
// reduced first. The result type is whatever the leftmost rule
// returns, which for the provided rules is the array's value type.
//
// Preconditions, checked in debug builds where the array exposes
// rank() and dim():
// - The array's runtime rank equals sizeof...(Rules).
// - All sample indices are valid.
// Each rule's numerical preconditions must also be satisfied.
//
// The array must be accessible from the calling execution space.
template <class Array, class... Rules>
  requires(sizeof...(Rules) > 0)
PORTABLE_INLINE_FUNCTION auto interpolate_with(const Array &array,
                                               const Rules &...rules) {
  constexpr int N = sizeof...(Rules);
  if constexpr (requires { array.rank(); }) {
    PORTABLE_REQUIRE(array.rank() == N,
                     "Number of interpolation rules must equal array rank");
  }
  return detail::evaluate<0, N>(
      array,
      [&](const auto... is) -> std::remove_cvref_t<decltype(array(is...))> {
        return array(is...);
      },
      rules...);
}

} // namespace interp
} // namespace Spiner

#endif // _SPINER_INTERP_TOY_HPP_
