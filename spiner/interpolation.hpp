#ifndef _SPINER_INTERP_
#define _SPINER_INTERP_
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

#include <assert.h>
#include <cmath>
#include <limits>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string>
#endif

#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_arrays.hpp"
#include "sp5.hpp"
#include "spiner_types.hpp"

namespace Spiner {
// TODO: be more careful about what this number should be
// sqrt machine epsilon or something
constexpr Real EPS = 10.0 * std::numeric_limits<Real>::epsilon();
constexpr Real rNaN = std::numeric_limits<Real>::signaling_NaN();
constexpr int iNaN = std::numeric_limits<int>::signaling_NaN();

// a poor-man's std::double
struct weights_t {
  Real first, second;
  PORTABLE_INLINE_FUNCTION Real &operator[](const int i) {
    assert(0 <= i && i <= 1);
    return i == 0 ? first : second;
  }
};

class RegularGrid1D {
 public:
  // Constructors
  PORTABLE_INLINE_FUNCTION RegularGrid1D()
      : min_(rNaN), max_(rNaN), dx_(rNaN), idx_(rNaN), N_(iNaN) {}
  PORTABLE_INLINE_FUNCTION RegularGrid1D(Real min, Real max, size_t N)
      : min_(min), max_(max), dx_((max - min) / ((Real)(N - 1))), idx_(1 / dx_),
        N_(N) {}

  // Assignment operator
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

  // Forces x in the interval
  PORTABLE_INLINE_FUNCTION int bound(int ix) const {
#ifndef SPINER_DISABLE_BOUNDS_CHECKS
    if (ix < 0) ix = 0;
    if (ix >= (int)N_ - 1) ix = (int)N_ - 2; // Ensures ix+1 exists
#endif
    return ix;
  }

  // Gets real value at index
  PORTABLE_INLINE_FUNCTION Real x(const int i) const { return i * dx_ + min_; }
  PORTABLE_INLINE_FUNCTION int index(const Real x) const {
    return bound(idx_ * (x - min_));
  }

  // Returns closest index and weights for interpolation
  PORTABLE_INLINE_FUNCTION void weights(Real x, int &ix, weights_t &w) const {
    ix = index(x);
    const Real floor = ix * dx_ + min_;
    w[1] = idx_ * (x - floor);
    w[0] = (1. - w[1]);
  }

  // 1D interpolation
  PORTABLE_INLINE_FUNCTION Real
  operator()(const Real x, const PortableMDArray<Real> &A) const {
    int ix;
    weights_t w;
    weights(x, ix, w);
    return w[0] * A(ix) + w[1] * A(ix + 1);
  }

  // utitilies
  PORTABLE_INLINE_FUNCTION bool operator==(const RegularGrid1D &other) const {
    return (min_ == other.min_ && max_ == other.max_ && dx_ == other.dx_ &&
            idx_ == other.idx_ && N_ == other.N_);
  }
  PORTABLE_INLINE_FUNCTION bool operator!=(const RegularGrid1D &other) const {
    return !(*this == other);
  }
  PORTABLE_INLINE_FUNCTION Real min() const { return min_; }
  PORTABLE_INLINE_FUNCTION Real max() const { return max_; }
  PORTABLE_INLINE_FUNCTION Real dx() const { return dx_; }
  PORTABLE_INLINE_FUNCTION Real nPoints() const { return N_; }
  PORTABLE_INLINE_FUNCTION bool isnan() const {
    return (std::isnan(min_) || std::isnan(max_) || std::isnan(dx_) ||
            std::isnan(idx_) || std::isnan((Real)N_));
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const { return !isnan(); }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    herr_t status;
    Real range[] = {min_, max_, dx_};
    hsize_t range_dims[] = {3};
    int n = static_cast<int>(N_);
    status = H5LTmake_dataset(loc, name.c_str(), SP5::RG1D::RANGE_RANK,
                              range_dims, H5T_REAL, range);
    status += H5LTset_attribute_int(loc, name.c_str(), SP5::RG1D::N, &n, 1);
    status += H5LTset_attribute_string(
        loc, name.c_str(), SP5::RG1D::RANGE_INFONAME, SP5::RG1D::RANGE_INFO);
    return status;
  }

  inline herr_t loadHDF(hid_t loc, const std::string &name) {
    herr_t status;
    Real range[3];
    int n;
    status = H5LTread_dataset(loc, name.c_str(), H5T_REAL, range);
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
  Real min_, max_;
  Real dx_, idx_;
  size_t N_;
};

} // namespace Spiner
#endif // _SPINER_INTERP_
