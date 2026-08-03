#ifndef SPINER_HIERARCHICAL_GRID_1D_
#define SPINER_HIERARCHICAL_GRID_1D_
//======================================================================
// © (or copyright) 2019-2026. Triad National Security, LLC. All rights
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

#include <array>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include "regular_grid_1d.hpp"

namespace Spiner {

template <typename T = Real, int NGRIDSMAX = 5,
          typename =
              typename std::enable_if<std::is_arithmetic<T>::value, bool>::type>
class PiecewiseGrid1D {
 public:
  using ValueType = T;
  static constexpr int BAD_VALUE = -1;

  // __host__ __device__ default constructors cause warning.
  // This is functionally equivalent because grids_ will
  // be initialized to default values
  PORTABLE_INLINE_FUNCTION PiecewiseGrid1D() {}
  PORTABLE_INLINE_FUNCTION
  PiecewiseGrid1D(const PiecewiseGrid1D &) = default;
  PORTABLE_INLINE_FUNCTION PiecewiseGrid1D &
  operator=(const PiecewiseGrid1D &) = default;
  PiecewiseGrid1D(const std::vector<RegularGrid1D<T>> grids) {
    NGRIDS_ = grids.size();
    PORTABLE_ALWAYS_REQUIRE(
        NGRIDS_ <= NGRIDSMAX,
        "Total number of grids must be within maximum allowed");
    int point_tot = 0;
    for (int i = 0; i < NGRIDS_; ++i) {
      grids_[i] = grids[i];
      pointTotals_[i] = point_tot;
      point_tot += grids[i].nPoints();
      if (i > 0) {
        const Real diff = std::abs(grids_[i].min() - grids_[i - 1].max());
        const Real avg =
            0.5 * (std::abs(grids_[i].min()) + std::abs(grids_[i - 1].max()));
        if (ratio_(diff, avg) >= EPS_()) {
          PORTABLE_ALWAYS_THROW_OR_ABORT(
              "Grids must be ordered and intersect at exactly one point.");
        }
      }
    }
  }
  PiecewiseGrid1D(std::initializer_list<RegularGrid1D<T>> grids)
      : PiecewiseGrid1D(std::vector<RegularGrid1D<T>>(grids)) {}

  template <typename F>
  PORTABLE_INLINE_FUNCTION int findGrid(const F &direction) const {
    int l = 0;
    int r = NGRIDS_ - 1;
    for (int iter = 0; iter < NGRIDS_; iter++) {
      int m = (l + r) / 2;
      int d = direction(m);
      if (d < 0) {
        if (m == 0) return 0;
        r = m - 1;
      } else if (d > 0) {
        if (m >= (NGRIDS_ - 1)) return NGRIDS_ - 1;
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
      if ((m < NGRIDS_ - 1) && (i >= pointTotals_[m + 1])) return 1;
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

  PORTABLE_INLINE_FUNCTION DataStatus dataStatus() const {
    int status = static_cast<int>(DataStatus::Trivial);
    for (int i = 0; i < NGRIDS_; ++i) {
      int gs = static_cast<int>(grids_[i].dataStatus());
      PORTABLE_REQUIRE(
          !((status == static_cast<int>(DataStatus::AllocatedHost)) &&
            (gs == static_cast<int>(DataStatus::AllocatedDevice))),
          "Can't mix allocated host/allocated device!");
      status = std::max(gs, status);
    }
    return static_cast<DataStatus>(status);
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

  PORTABLE_INLINE_FUNCTION T min() const { return grids_[0].min(); }
  PORTABLE_INLINE_FUNCTION T max() const { return grids_[NGRIDS_ - 1].max(); }
  PORTABLE_INLINE_FUNCTION size_t nPoints() const {
    return pointTotals_[NGRIDS_ - 1] + grids_[NGRIDS_ - 1].nPoints();
  }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    for (int ig = 0; ig < NGRIDS_; ++ig) {
      if (grids_[ig].isnan()) return true;
    }
    return false;
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const {
    return NGRIDS_ > 0 && !isnan();
  }
  PORTABLE_INLINE_FUNCTION int nGrids() const { return NGRIDS_; }

  // Binary serialization is intended for transient communication between
  // compatible Spiner builds, not as a persistent or portable file format.
  std::size_t dynamicMemorySizeInBytes() const {
    std::size_t size = 0;
    for (int i = 0; i < NGRIDS_; ++i) {
      size += grids_[i].dynamicMemorySizeInBytes();
    }
    return size;
  }

  std::size_t serializedSizeInBytes() const {
    return sizeof(*this) + dynamicMemorySizeInBytes();
  }

  std::size_t dumpDynamicMemory(std::byte *dst) const {
    std::size_t offset = 0;
    for (int i = 0; i < NGRIDS_; ++i) {
      offset += grids_[i].dumpDynamicMemory(dst + offset);
    }
    return offset;
  }

  std::size_t serialize(std::byte *dst) const {
    std::memcpy(dst, this, sizeof(*this));
    return sizeof(*this) + dumpDynamicMemory(dst + sizeof(*this));
  }

  std::size_t setPointer(std::byte *src) {
    std::size_t offset = 0;
    for (int i = 0; i < NGRIDS_; ++i) {
      offset += grids_[i].setPointer(src + offset);
    }
    return offset;
  }

  std::size_t deSerialize(std::byte *src) {
    finalize(); // TODO(JMM): Maybe guard this finalize
    std::memcpy(this, src, sizeof(*this));
    PORTABLE_REQUIRE(0 <= NGRIDS_ && NGRIDS_ <= NGRIDSMAX,
                     "Invalid number of piecewise grids");
    return sizeof(*this) + setPointer(src + sizeof(*this));
  }

  PiecewiseGrid1D<T, NGRIDSMAX> getOnDevice() const {
    PiecewiseGrid1D<T, NGRIDSMAX> grid(*this);
    for (int i = 0; i < NGRIDS_; ++i) {
      grid.grids_[i] = grids_[i].getOnDevice();
    }
    return grid;
  }

  void finalize() {
    for (int i = 0; i < NGRIDS_; ++i) {
      grids_[i].finalize(); // TODO(JMM): Maybe guard this finalize
    }
    NGRIDS_ = 0;
  }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    herr_t status = 0;
    hid_t group =
        H5Gcreate(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int ngrids = static_cast<int>(NGRIDS_);
    status =
        H5LTset_attribute_int(loc, name.c_str(), SP5::H1D::NGRIDS, &ngrids, 1);
    for (int i = 0; i < NGRIDS_; ++i) {
      status += grids_[i].saveHDF(group, gridname_(i).c_str());
    }
    status += H5Gclose(group);
    return status;
  }

  inline herr_t loadHDF(hid_t loc, const std::string &name) {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    herr_t status = 0;
    int ngrids;
    int is_piecewise_grid =
        H5Aexists_by_name(loc, name.c_str(), SP5::H1D::NGRIDS, H5P_DEFAULT);
    if (is_piecewise_grid) {
      hid_t group = H5Gopen(loc, name.c_str(), H5P_DEFAULT);
      H5LTget_attribute_int(loc, name.c_str(), SP5::H1D::NGRIDS, &ngrids);
      NGRIDS_ = ngrids;
      PORTABLE_ALWAYS_REQUIRE(
          NGRIDS_ <= NGRIDSMAX,
          "Total number of grids must be within maximum allowed");
      int point_tot = 0;
      for (int i = 0; i < NGRIDS_; ++i) {
        status += grids_[i].loadHDF(group, gridname_(i).c_str());
        pointTotals_[i] = point_tot;
        point_tot += grids_[i].nPoints();
        if ((i > 0) &&
            ratio_(2 * std::abs(grids_[i].min() - grids_[i - 1].max()),
                   std::abs(grids_[i].min() + grids_[i - 1].max()) >= EPS_())) {
          PORTABLE_ALWAYS_THROW_OR_ABORT(
              "Grids must be ordered and intersect at exactly one point.");
        }
      }
      status += H5Gclose(group);
    } else {
      NGRIDS_ = 1;
      pointTotals_[0] = 0;
      status += grids_[0].loadHDF(loc, name);
    }
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

  RegularGrid1D<T> grids_[NGRIDSMAX];
  int pointTotals_[NGRIDSMAX]{};
  int NGRIDS_ = 0;
};

} // namespace Spiner

#endif // SPINER_HIERARCHICAL_GRID_1D_
