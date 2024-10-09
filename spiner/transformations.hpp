#ifndef _SPINER_TRANSFORM_HPP_
#define _SPINER_TRANSFORM_HPP_
//======================================================================
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
//======================================================================

#include "ports-of-call/portability.hpp"

#include <cmath>
#include <limits>

namespace Spiner {

// Note on notation:
// -- "real" space is called x
// -- "transformed" space is called u
// -- u = forward(x)
// -- x = reverse(u)

// linear transformation (aka no-op): y = x
struct TransformLinear {
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T forward(const T x) {
    return x;
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T reverse(const T u) {
    return u;
  }
};

// logarithmic transformation: y = log(x + small)
struct TransformLogarithmic {
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T forward(const T x) {
    constexpr T eps = eps_f<T>();
    return std::log(std::abs(x) + eps);
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T reverse(const T u) {
    constexpr T eps = eps_r<T>();
    return std::exp(u) - eps;
  }
private:
  // When possible, we use asymetric epsilon values to ensure that
  // reverse(forward(0)) is exact.  In general, a performant calculation is
  // more important than getting this value exactly correct, so we require that
  // our epsilon values be constexpr.
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T eps_f() {
    return std::numeric_limits<T>::min();
  }
  template <typename T>
  PORTABLE_FORCEINLINE_FUNCTION static constexpr T eps_r() {
    // If std::exp and std::log were constexpr, we could explicitly calculate
    // the right epsilon value as a constexpr value:
    //     return std::exp(forward<T>(static_cast<T>(0)));
    // Unfortunately, we have to wait for C++26 for constexpr math.  We
    // hard-code certain known values, but default to a symmetric epsilon if
    // that's the best that we can do.
    if constexpr (std::is_same_v<std::decay_t<T>, float> &&
                  std::numeric_limits<T>::is_iec559) {
      return 1.175490707446280263444352135525329e-38f;
    } else if constexpr (std::is_same_v<std::decay_t<T>, double> &&
                         std::numeric_limits<T>::is_iec559) {
      return 2.225073858507262647230317031903882e-308;
    } else {
      return eps_f<T>();
    }
  }
};

// TODO: log_NQT and arcsinh_NQT, but these require adding a dependency on
//       https://github.com/lanl/not-quite-transcendental.  I may leave this for
//       Jonah ;-)

} // namespace Spiner
#endif // _SPINER_TRANSFORM_HPP_
