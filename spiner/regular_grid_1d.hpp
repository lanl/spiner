#ifndef SPINER_REGULAR_GRID_1D_
#define SPINER_REGULAR_GRID_1D_
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

#include <assert.h>
#include <cmath>
#include <limits>
#include <type_traits>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string>
#endif

#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_arrays.hpp"
#include "ports-of-call/portable_errors.hpp"
#include "sp5.hpp"
#include "spiner_types.hpp"

namespace Spiner {

// a poor-man's std::pair
template <typename T = Real>
struct weights_t {
  T first, second;
  PORTABLE_INLINE_FUNCTION Real &operator[](const int i) {
    assert(0 <= i && i <= 1);
    return i == 0 ? first : second;
  }
};

template <typename T = Real,
          typename std::enable_if<std::is_arithmetic<T>::value, bool>::type =
              true>
class RegularGrid1D {
 public:
  using ValueType = T;
  static constexpr T rNaN = std::numeric_limits<T>::signaling_NaN();
  static constexpr int iNaN = std::numeric_limits<int>::signaling_NaN();

  // Constructors
  PORTABLE_INLINE_FUNCTION RegularGrid1D()
      : min_(rNaN), max_(rNaN), dx_(rNaN), idx_(rNaN), N_(iNaN) {}
  PORTABLE_INLINE_FUNCTION RegularGrid1D(T min, T max, size_t N)
      : min_(min), max_(max), dx_((max - min) / ((T)(N - 1))), idx_(1 / dx_),
        N_(N) {
    PORTABLE_ALWAYS_REQUIRE(min_ < max_ && N_ > 0, "Valid grid");
  }

  // Assignment operator
  /*
  Default copy constructable
  PORTABLE_INLINE_FUNCTION RegularGrid1D &operator=(const RegularGrid1D &src) {
    if (this != &src) {
      min_ = src.min_;
      max_ = src.max_;
      dx_ = src.dx_;
      idx_ = src.idx_;
      N_ = src.N_;
    }
    return *this;
  }
  */

  // Forces x in the interval
  PORTABLE_INLINE_FUNCTION int bound(int ix) const {
#ifndef SPINER_DISABLE_BOUNDS_CHECKS
    if (ix < 0) ix = 0;
    if (ix >= (int)N_ - 1) ix = (int)N_ - 2; // Ensures ix+1 exists
#endif
    return ix;
  }

  // Gets real value at index
  PORTABLE_INLINE_FUNCTION T x(const int i) const { return i * dx_ + min_; }
  PORTABLE_INLINE_FUNCTION int index(const T x) const {
    return bound(idx_ * (x - min_));
  }

  // Returns closest index and weights for interpolation
  PORTABLE_INLINE_FUNCTION void weights(const T &x, int &ix,
                                        weights_t<T> &w) const {
    ix = index(x);
    const auto floor = static_cast<T>(ix) * dx_ + min_;
    w[1] = idx_ * (x - floor);
    w[0] = (1. - w[1]);
  }

  // 1D interpolation
  PORTABLE_INLINE_FUNCTION T operator()(const T &x,
                                        const PortableMDArray<T> &A) const {
    int ix;
    weights_t<T> w;
    weights(x, ix, w);
    return w[0] * A(ix) + w[1] * A(ix + 1);
  }

  // utitilies
  PORTABLE_INLINE_FUNCTION bool
  operator==(const RegularGrid1D<T> &other) const {
    return (min_ == other.min_ && max_ == other.max_ && dx_ == other.dx_ &&
            idx_ == other.idx_ && N_ == other.N_);
  }
  PORTABLE_INLINE_FUNCTION bool
  operator!=(const RegularGrid1D<T> &other) const {
    return !(*this == other);
  }
  PORTABLE_INLINE_FUNCTION T min() const { return min_; }
  PORTABLE_INLINE_FUNCTION T max() const { return max_; }
  PORTABLE_INLINE_FUNCTION T dx() const { return dx_; }
  PORTABLE_INLINE_FUNCTION size_t nPoints() const { return N_; }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    return (std::isnan(min_) || std::isnan(max_) || std::isnan(dx_) ||
            std::isnan(idx_) || std::isnan((T)N_));
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const { return !isnan(); }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    auto H5T_T =
        std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    herr_t status;
    T range[] = {min_, max_, dx_};
    hsize_t range_dims[] = {3};
    int n = static_cast<int>(N_);
    status = H5LTmake_dataset(loc, name.c_str(), SP5::RG1D::RANGE_RANK,
                              range_dims, H5T_T, range);
    status += H5LTset_attribute_int(loc, name.c_str(), SP5::RG1D::N, &n, 1);
    status += H5LTset_attribute_string(
        loc, name.c_str(), SP5::RG1D::RANGE_INFONAME, SP5::RG1D::RANGE_INFO);
    return status;
  }

  inline herr_t loadHDF(hid_t loc, const std::string &name) {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    auto H5T_T =
        std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    herr_t status;
    T range[3];
    int n;
    status = H5LTread_dataset(loc, name.c_str(), H5T_T, range);
    min_ = range[0];
    max_ = range[1];
    dx_ = range[2];
    idx_ = 1. / dx_;
    status += H5LTget_attribute_int(loc, name.c_str(), SP5::RG1D::N, &n);
    N_ = n;
    return status;
  }
#endif

 private:
  T min_, max_;
  T dx_, idx_;
  size_t N_;
};

} // namespace Spiner
#endif // SPINER_REGULAR_GRID_1D_
