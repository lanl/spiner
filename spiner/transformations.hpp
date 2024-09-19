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

#include <limits>

namespace Spiner {

  // Note on notation:
  // -- "real" space is called x
  // -- "transformed" space is called u
  // -- u = forward(x)
  // -- x = reverse(u)

  // linear transformation (aka no-op): y = x
  struct TransformLinear{
    template<typename T>
    PORTABLE_INLINE_FUNCTION static T forward(const T x) {
      return x;
    }
    template<typename T>
    PORTABLE_INLINE_FUNCTION static T reverse(const T u) {
      return u;
    }
  };

  // logarithmic transformation: y = log(x + small)
  struct TransformLogarithmic{
    template<typename T>
    PORTABLE_INLINE_FUNCTION static T forward(const T x) {
      return std::log(x + std::numeric_limits<T>::denorm_min());
    }
    template<typename T>
    PORTABLE_INLINE_FUNCTION static T reverse(const T u) {
      return std::exp(u) - std::numeric_limits<T>::denorm_min();
    }
  };

  // TODO: log_NQT

  // TODO: arcsinh_NQT

} // namespace Spiner
#endif // _SPINER_TRANSFORM_HPP_
