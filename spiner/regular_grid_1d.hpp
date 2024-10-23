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
#include "transformations.hpp"

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

template <typename T = Real, typename Transform = TransformLinear,
          typename std::enable_if<std::is_arithmetic<T>::value, bool>::type =
              true>
class RegularGrid1D {
public:
  using ValueType = T;
  static constexpr T rNaN = std::numeric_limits<T>::signaling_NaN();
  static constexpr int iNaN = std::numeric_limits<int>::signaling_NaN();

  // Constructors
  PORTABLE_INLINE_FUNCTION RegularGrid1D()
      : umin_(rNaN), umax_(rNaN), du_(rNaN), inv_du_(rNaN), N_(iNaN) {}
  PORTABLE_INLINE_FUNCTION RegularGrid1D(T xmin, T xmax, size_t N)
      : xmin_(xmin)
      , xmax_(xmax)
      , umin_(Transform::forward(xmin))
      , umax_(Transform::forward(xmax))
      , du_((umax_ - umin_) / static_cast<T>(N - 1))
      , inv_du_(1 / du_)
      , N_(N)
  {
    // A transform could be monotonically decreasing, so there's no guarantee
    // that umin_ < umax_
    PORTABLE_ALWAYS_REQUIRE(xmin_ < xmax_ && N_ > 0, "Valid grid");
  }

  // Returns closest index and weights for interpolation
  PORTABLE_INLINE_FUNCTION void weights(const T &x, int &ix, weights_t<T> &w) const {
    const T u = Transform::forward(x);
    ix = index_u(u);
    const auto floor = static_cast<T>(ix) * du_ + umin_;
    w[1] = inv_du_ * (u - floor);
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

  // (in)equality comparison
  PORTABLE_INLINE_FUNCTION bool
  operator==(const RegularGrid1D<T, Transform> &other) const {
    return (umin_ == other.umin_ &&
            umax_ == other.umax_ &&
            du_ == other.du_ &&
            N_ == other.N_);
  }
  PORTABLE_INLINE_FUNCTION bool
  operator!=(const RegularGrid1D<T, Transform> &other) const {
    return !(*this == other);
  }

  // queries
  PORTABLE_INLINE_FUNCTION T min() const { return xmin_; }
  PORTABLE_INLINE_FUNCTION T max() const { return xmax_; }
  PORTABLE_INLINE_FUNCTION size_t nPoints() const { return N_; }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    return (std::isnan(xmin_) ||
            std::isnan(xmax_) ||
            std::isnan(umin_) ||
            std::isnan(umax_) ||
            std::isnan(du_) ||
            std::isnan(inv_du_) ||
            std::isnan((T)N_));
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const { return !isnan(); }
  // TODO: min() and x(0) won't necessarily match.
  //       max() and x(nPoints-1) won't necessarily match.
  //       Should we do anything about this?
  //       The easiest fix: if i == 0 return xmin_ and similar for xmax_
  // Translate between x coordinate and index
  PORTABLE_INLINE_FUNCTION T x(const int i) const {
    return Transform::reverse(u(i));
  }
  PORTABLE_INLINE_FUNCTION int index(const T x) const {
    return index_u(Transform::forward(x));
  }

  // HDF
#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    auto H5T_T =
        std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    herr_t status;
    T range[] = {umin_, umax_, du_};
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
    umin_ = range[0];
    umax_ = range[1];
    du_ = range[2];
    inv_du_ = 1. / du_;
    status += H5LTget_attribute_int(loc, name.c_str(), SP5::RG1D::N, &n);
    N_ = n;
    return status;
  }
#endif

private:
  // Forces x in the interval
  PORTABLE_INLINE_FUNCTION int bound(int ix) const {
#ifndef SPINER_DISABLE_BOUNDS_CHECKS
    if (ix < 0) ix = 0;
    if (ix >= (int)N_ - 1) ix = (int)N_ - 2; // Ensures ix+1 exists
#endif
    return ix;
  }

  // Translate between u (transformed variable) coordinate and index
  PORTABLE_INLINE_FUNCTION T u(const int i) const { return i * du_ + umin_; }
  PORTABLE_INLINE_FUNCTION int index_u(const T u) const {
    return bound(inv_du_ * (u - umin_));
  }

  T xmin_, xmax_;
  T umin_, umax_;
  T du_, inv_du_;
  size_t N_;
};

} // namespace Spiner
#endif // SPINER_REGULAR_GRID_1D_
