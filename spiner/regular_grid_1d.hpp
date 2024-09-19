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

// TODO: Need to expand testing

// TODO: Do transformations need state?
//    -- My thinking is that transformations generally won't need state, because it's things like
//       y = x or y = log(x).
//    -- Based on this, I would invoke the transform as Transform::forward(x) or
//       Transform::reverse(u).
//    -- But this means that transformations _cannot_ have state if that use-case comes up in the
//       future.
//    -- The other option would be to have the transformation object be passed into the
//       constructor.
//       -- We could default to Transform(), but that means all transformations _must_ have a
//          default constructor even if it makes no sense.
//       -- We could have no default, but then users are now required to update all their
//          constructors to add TransformLinear() as an argument, and we can't hide the addition of
//          the template.

template <typename T = Real,
          typename Transform = TransformLinear,
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
      : umin_(Transform::forward(xmin))
      , umax_(Transform::forward(xmax))
      , du_((umax_ - umin_) / ((T)(N - 1)))
      , inv_du_(1 / du_)
      , N_(N) {
    PORTABLE_ALWAYS_REQUIRE(umin_ < umax_ && N_ > 0, "Valid grid");
  }

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

  // Translate between x coordinate and index
  PORTABLE_INLINE_FUNCTION T x(const int i) const { return Transform::reverse(u(i)); }
  PORTABLE_INLINE_FUNCTION int index(const T x) const {
    return index_u(Transform::forward(x));
  }

  // Returns closest index and weights for interpolation
  PORTABLE_INLINE_FUNCTION void weights(const T &x, int &ix,
                                        weights_t<T> &w) const {
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

  // utitilies
  PORTABLE_INLINE_FUNCTION bool
  operator==(const RegularGrid1D<T> &other) const {
    return (umin_ == other.umin_ && umax_ == other.umax_ && du_ == other.du_ &&
            inv_du_ == other.inv_du_ && N_ == other.N_);
  }
  PORTABLE_INLINE_FUNCTION bool
  operator!=(const RegularGrid1D<T> &other) const {
    return !(*this == other);
  }
  // TODO: This way of constructing min() and max() implicitly asserts that the u-space
  //       representation is the "ground truth" and the x-space representation is merely derived
  //       from there.  We could change this and assert that the x-space representation is the
  //       "ground truth", but we would have to make some changes.  This arises because it's not
  //       guaranteed that every transformation is _exactly_ one-to-one reversible, so you could
  //       introduce a small gap between xmin and reverse(umin) or between forward(xmin) and umin
  //       (and similarly for max).
  // TODO: Should umin() and umax() not be publicly exposed?  If so, they shouldn't exist because
  //       they're only public getters for the private umin_ / umax_ variables.
  PORTABLE_INLINE_FUNCTION T umin() const { return umin_; }
  PORTABLE_INLINE_FUNCTION T umax() const { return umax_; }
  PORTABLE_INLINE_FUNCTION T min() const { return Transform::reverse(umin_); }
  PORTABLE_INLINE_FUNCTION T max() const { return Transform::reverse(umax_); }
  // TODO: Should there be a public du() function?
  // TODO: dx() is now ill-defined -- do we need it?
  //PORTABLE_INLINE_FUNCTION T dx() const { return dx_; }
  PORTABLE_INLINE_FUNCTION size_t nPoints() const { return N_; }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    return (std::isnan(umin_) || std::isnan(umax_) || std::isnan(du_) ||
            std::isnan(inv_du_) || std::isnan((T)N_));
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const { return !isnan(); }

  // TODO: I made very simple, naive changes to saveHDF and loadHDF, because it's not clear to me
  //       if something better should be done.
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
  T umin_, umax_;
  T du_, inv_du_;
  size_t N_;
};

} // namespace Spiner
#endif // SPINER_REGULAR_GRID_1D_
