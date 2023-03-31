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
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include "regular_grid_1d.hpp"

namespace Spiner {

template <typename T = Real,
	  int NGRIDS = 3,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, bool>::type>
class HierarchicalGrid1D {
public:
  using ValueType = T;
  static constexpr int BAD_VALUE = -1;

  // __host__ __device__ default constructors cause warning.
  // This is functionally equivalent because grids_ will
  // be initialized to default values
  PORTABLE_INLINE_FUNCTION HierarchicalGrid1D() {}
  PORTABLE_INLINE_FUNCTION HierarchicalGrid1D(const RegularGrid1D<T> grids[NGRIDS]) {
    int point_tot = 0;
    for (int i = 0; i < NGRIDS; ++i) {
      grids_[i] = grids[i];
      pointTotals_[i] = point_tot;
      point_tot += grids[i].nPoints();
      if ( (i > 0) && ratio_(2*std::abs(grids_[i].min() - grids_[i-1].max()),
                             std::abs(grids_[i].min() + grids_[i-1].max()) >= EPS_() )) {
        PORTABLE_ALWAYS_THROW_OR_ABORT("Grids in the hierarchy must be ordered "
                                       "and intersect at exactly one point.");
      }
    }
  }
  HierarchicalGrid1D(std::initializer_list<RegularGrid1D<T>> grids)
    : HierarchicalGrid1D(std::vector<RegularGrid1D<T>>(grids).data()) {}

  template<typename F>
  PORTABLE_INLINE_FUNCTION int findGrid(const F &direction) const {
    int l = 0;
    int r = NGRIDS - 1;
    for (int iter = 0; iter < NGRIDS; iter++) {
      int m = (l + r) / 2;
      int d = direction(m);
      if (d < 0) {
        if (m == 0) return 0;
        r = m - 1;
      } else if (d > 0) {
        if (m >= (NGRIDS - 1)) return NGRIDS - 1;
        l = m + 1;
      } else {
        return m;
      }
    }
    PORTABLE_ABORT("Grid find failed");
    return BAD_VALUE;
  }
  PORTABLE_INLINE_FUNCTION int findGridFromGlobalIdx(const int i) const {
    int ig = findGrid([&](const int m) {
      if (i < pointTotals_[m]) return -1;
      if ((m < NGRIDS - 1) && (i >= pointTotals_[m+1])) return 1;
      return 0;
    });
    return ig;
  }
  PORTABLE_INLINE_FUNCTION int findGridFromPosition(const T x) const {
    int ig = findGrid([&](const int m) {
      if (x < grids_[m].min()) return -1;
      if (x > grids_[m].max()) return 1;
      return 0;
    });
    return ig;
  }

  PORTABLE_INLINE_FUNCTION T x(const int i) const {
    int ig = findGridFromGlobalIdx(i);
    return grids_[ig].x(i - pointTotals_[ig]);
  }

  PORTABLE_INLINE_FUNCTION int index(const T x) const {
    int ig = findGridFromPosition(x);
    return pointTotals_[ig] + grids_[ig].index(x);
  }

  // Returns closest index and weights for interpolation
  PORTABLE_INLINE_FUNCTION void weights(const T &x, int &ix,
                                        weights_t<T> &w) const {
    int ig = findGridFromPosition(x);
    grids_[ig].weights(x, ix, w);
    ix += pointTotals_[ig];
  }

  PORTABLE_INLINE_FUNCTION
  bool operator==(const HierarchicalGrid1D<T,NGRIDS> &other) const {
    for (int ig = 0; ig < NGRIDS; ++ig) {
      if (grids_[ig] != other.grids_[ig]) return false;
    }
    return true;
  }
  PORTABLE_INLINE_FUNCTION
  bool operator!=(const HierarchicalGrid1D<T, NGRIDS> &other) const {
    return !(*this == other);
  }
  PORTABLE_INLINE_FUNCTION T min() const {
    return grids_[0].min();
  }
  PORTABLE_INLINE_FUNCTION T max() const {
    return grids_[NGRIDS-1].max();
  }
  PORTABLE_INLINE_FUNCTION size_t nPoints() const {
    return pointTotals_[NGRIDS-1] + grids_[NGRIDS-1].nPoints();
  }
  PORTABLE_INLINE_FUNCTION T dx(const int ig) const {
    assert( ig < NGRIDS );
    return grids_[ig].dx();
  }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    for (int ig = 0; ig < NGRIDS; ++ig) {
      if (grids_[ig].isnan()) return true;
    }
    return false;
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const { return !isnan(); }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                  "Spiner HDF5 only defined for these data types: float, double");
    herr_t status = 0;
    hid_t group = H5Gcreate(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int ngrids = static_cast<int>(NGRIDS);
    status = H5LTset_attribute_int(loc, name.c_str(), SP5::H1D::NGRIDS, &ngrids, 1);
    for (int i = 0; i < NGRIDS; ++i) {
      status += grids_[i].saveHDF(group, gridname_(i).c_str());
    }
    status += H5Gclose(group);
    return status;
  }

  inline herr_t loadHDF(hid_t, const std::string &name) const {
    static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                  "Spiner HDF5 only defined for these data types: float, double");
    herr_t status = 0;
    hid_t group = H5Gopen(loc, name.c_str(), H5P_DEFAULT);
    int ngrids;
    H5LTget_attribute_int(loc, name.c_str(), SP5::H1D::NGRIDS, &ngrids);
    assert( ngrids == NGRIDS );
    int point_tot = 0;
    for (int i = 0; i < NGRIDS; ++i) {
      status += grids_[i].loadHDF(group, gridname_(i).c_str());
      pointTotals_[i] = point_tot;
      point_tot += grids_[i].nPoints();
      if ( (i > 0) && ratio_(2*std::abs(grids_[i].min() - grids_[i-1].max()),
                             std::abs(grids_[i].min() + grids_[i-1].max()) >= EPS_() )) {
        PORTABLE_ALWAYS_THROW_OR_ABORT("Grids in the hierarchy must be ordered "
                                       "and intersect at exactly one point.");
      }
    }
    status += H5Gclose(group);
    return status;
  }
#endif

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
  inline std::string gridname_(int i) const {
    return SP5::H1D::GRID_FORMAT[0] + std::to_string(i + 1) +
           SP5::H1D::GRID_FORMAT[1];
  }

  RegularGrid1D<T> grids_[NGRIDS];
  int pointTotals_[NGRIDS];
};

} // namespace Spiner

#endif // SPINER_HIERARCHICAL_GRID_1D_
