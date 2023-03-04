#ifndef SPINER_HIERARCHICAL_GRID_1D_
#define SPINER_HIERARCHICAL_GRID_1D_
//======================================================================
// © (or copyright) 2019-2023. Triad National Security, LLC. All rights
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

#include <array>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>

#include "regular_grid_1d.hpp"

namespace Spiner {

template <typename T = Real,
	  std::size_t NGRIDS = 3,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, bool>::type>
class HierarchicalGrid1D {
public:
  using ValueType = T;
  constexpr int BAD_VALUE = -1;

  // __host__ __device__ default constructors cause warning.
  // This is functionally equivalent because grids_ will
  // be initialized to default values
  PORTABLE_INLINE_FUNCTION HierarchicalGrid1D() {}
  PORTABLE_INLINE_FUNCTION HierarchicalGrid1D(const RegularGrid1D<T> grids[NGRIDS]) {
    for (int i = 0; i < NGRIDS; ++i) {
      grids_[i] = grids[i];
      if ( (i > 0) && ratio_(2*std::abs(grids_[i].min() - grids_[i-1].max()), std::abs(grids_[i].min() + grids_[i-1].max()) >= EPS_() ) {
	  PORTABLE_ALWAYS_THROW_OR_ABORT("Grids in the hierarchy must be ordered and intersect at exactly one point.")
      }
    }
  }
  PORTABLE_INLINE_FUNCTION HierarchicalGrid1D(std::initalizer_list<RegularGrid1D<T>> grids)
    : HierarchicalGrid1D(std::array<RegularGrid1D, NGRIDS>(grids)) {}
  

  PORTABLE_INLINE_FUNCTION int findGrid(const T &x) const {
    std::size_t l = 0;
    std::size_t r = NGRIDS - 1;
    for (int iter = 0; iter < NGRIDS; iter++) {
      std::size_t m = (l + r) / 2;
      if (x < grids_[m].min()) {
	r = m - 1;
      } else if (x > grids_ig.max()) {
	l = m + 1;
      } else { // we found it
	return ig;
      }
    }
    PORTABLE_FAIL("Grid find failed");
    return BAD_VALUE;
  }

  PORTABLE_INLINE_FUNCTION int index(const T x) const {
    int ig = findGrid(x);
    return grids_[ig].index(x);
  }

private:
  PORTABLE_FORCEINLINE_FUNCTION constexpr auto SMALL_() const {
    return 10 * std::numeric_limits<T>::min();
  }
  PORTABLE_FORCEINLINE_FUNCTION constexpr auto EPS_() const {
    return 10 * std::numeric_limits<T>::epsilon();
  }
  PORTABLE_FORCEINLINE_FUNCTION int sgn_(const T &val) const {
    return (T(0) <= val) - (val < T(0));
  }
  PORTABLE_FORCEINLINE_FUNCTION T ratio_(const T &a, const T &b) const {
    return a / (b + sgn_(b) * SMALL_());
  }

  RegularGrid1D<T> grids_[NGRIDS];
};

} // namespace Spiner

#endif // SPINER_HIERARCHICAL_GRID_1D_
