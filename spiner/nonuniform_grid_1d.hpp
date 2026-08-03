#ifndef SPINER_NONUNIFORM_GRID_1D_
#define SPINER_NONUNIFORM_GRID_1D_
//======================================================================
// © (or copyright) 2026. Triad National Security, LLC. All rights
// reserved.  This program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
// operated by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. All rights in the
// program are reserved by Triad National Security, LLC, and the
// U.S. Department of Energy/National Nuclear Security Administration.
//======================================================================

// Generative AI was used to assist with writing this file.

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <vector>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string>
#endif

#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_errors.hpp"
#include "regular_grid_1d.hpp"
#include "spiner_types.hpp"

namespace Spiner {

template <typename T = Real,
          typename std::enable_if<std::is_arithmetic<T>::value, bool>::type =
              true>
class NonUniformGrid1D {
 public:
  using ValueType = T;

  NonUniformGrid1D() = default;

  explicit NonUniformGrid1D(const std::vector<T> &points) {
    allocate_(points.data(), points.size());
  }
  NonUniformGrid1D(std::initializer_list<T> points)
      : NonUniformGrid1D(std::vector<T>(points)) {}

  // This constructor borrows caller-owned memory. The caller must
  // keep it alive and unchanged for this grid's lifetime, and is
  // responsible for ensuring that its coordinates are finite and
  // strictly increasing.  Validation is performed on the specified
  // PortsOfCall Execution Space. Code will segfault if the incorrect
  // space is specified.
  template <typename ExecutionSpace = PortsOfCall::Exec::Host>
  NonUniformGrid1D(
      T *points, const std::size_t n,
      const ExecutionSpace &exec_space = PortsOfCall::Exec::Host()) {
    if (n > 0) {
      n_ = n;
      data_ = points;
      status_ = DataStatus::Unmanaged;
    }
    validate_(exec_space);
  }

  PORTABLE_INLINE_FUNCTION T x(const int i) const { return data_[i]; }

  PORTABLE_INLINE_FUNCTION int index(const T value) const {
    if (value <= data_[0]) return 0;
    if (value >= data_[n_ - 1]) return static_cast<int>(n_) - 2;

    int lo = 0;
    int hi = n_ - 1;
    int iter = 0;
    while (hi - lo > 1) {
      // no-op if NDEBUG not set
      PORTABLE_REQUIRE(iter++ < n_ + 1,
                       "This binary search failed to converge/terminate.");
      const int mid = lo + (hi - lo) / 2;
      if (data_[mid] <= value) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  PORTABLE_INLINE_FUNCTION void weights(const T &value, int &ix,
                                        weights_t<T> &w) const {
    ix = index(value);
    const T dx = data_[ix + 1] - data_[ix];
    w[1] = (value - data_[ix]) / dx;
    w[0] = T(1) - w[1];
  }

  PORTABLE_INLINE_FUNCTION T min() const { return data_[0]; }
  PORTABLE_INLINE_FUNCTION T max() const { return data_[n_ - 1]; }
  PORTABLE_INLINE_FUNCTION std::size_t nPoints() const { return n_; }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const {
    return pointsWellFormed(data_, n_);
  }
  PORTABLE_INLINE_FUNCTION DataStatus dataStatus() const { return status_; }
  PORTABLE_INLINE_FUNCTION T *data() const { return data_; }

  std::size_t dynamicMemorySizeInBytes() const { return n_ * sizeof(T); }
  std::size_t serializedSizeInBytes() const {
    return sizeof(*this) + dynamicMemorySizeInBytes();
  }
  std::size_t dumpDynamicMemory(std::byte *dst) const {
    PORTABLE_REQUIRE(status_ != DataStatus::AllocatedDevice,
                     "Cannot dump device-resident grid memory");
    if (n_ > 0) std::memcpy(dst, data_, dynamicMemorySizeInBytes());
    return dynamicMemorySizeInBytes();
  }
  std::size_t serialize(std::byte *dst) const {
    PORTABLE_REQUIRE(status_ != DataStatus::AllocatedDevice,
                     "Cannot serialize device-resident grid memory");
    std::memcpy(dst, this, sizeof(*this));
    return sizeof(*this) + dumpDynamicMemory(dst + sizeof(*this));
  }
  std::size_t setPointer(std::byte *src) {
    data_ = n_ == 0 ? nullptr : reinterpret_cast<T *>(src);
    status_ = n_ == 0 ? DataStatus::Empty : DataStatus::Unmanaged;
    if (n_ > 0) validate_(PortsOfCall::Exec::Host{});
    return dynamicMemorySizeInBytes();
  }
  std::size_t deSerialize(std::byte *src) {
    PORTABLE_REQUIRE(
        (status_ == DataStatus::Empty || status_ == DataStatus::Unmanaged),
        "Must not de-serialize into an active grid.");
    std::memcpy(this, src, sizeof(*this));
    PORTABLE_REQUIRE(n_ == 0 || n_ >= 2, "Invalid nonuniform grid size");
    return sizeof(*this) + setPointer(src + sizeof(*this));
  }

  NonUniformGrid1D getOnDevice() const {
    PORTABLE_REQUIRE(status_ != DataStatus::AllocatedDevice,
                     "Cannot copy a device grid to device");
    NonUniformGrid1D grid;
    grid.n_ = n_;
    if (n_ > 0) {
      grid.data_ =
          static_cast<T *>(PORTABLE_MALLOC(dynamicMemorySizeInBytes()));
      PORTABLE_ALWAYS_REQUIRE(grid.data_ != nullptr, "Grid allocation failed");
      portableCopyToDevice(grid.data_, data_, dynamicMemorySizeInBytes());
      grid.status_ = DataStatus::AllocatedDevice;
    }
    return grid;
  }

  // Explicitly make an independent host-owned copy. Normal copy construction
  // and assignment remain shallow.
  void copy(const NonUniformGrid1D &other) {
    if (this == &other) return;
    PORTABLE_REQUIRE(other.status_ != DataStatus::AllocatedDevice,
                     "Cannot deep copy a device-resident grid to host");
    T *data = nullptr;
    if (other.n_ > 0) {
      data = static_cast<T *>(std::malloc(other.dynamicMemorySizeInBytes()));
      PORTABLE_ALWAYS_REQUIRE(data != nullptr, "Grid allocation failed");
      std::memcpy(data, other.data_, other.dynamicMemorySizeInBytes());
    }
    PORTABLE_REQUIRE(
        (status_ == DataStatus::Empty || status_ == DataStatus::Unmanaged),
        "Must not copy into an active grid.");
    n_ = other.n_;
    data_ = data;
    status_ = (n_ == 0) ? DataStatus::Empty : DataStatus::AllocatedHost;
  }

  void finalize() {
    // Note that finalizes are projections for Grid objects (unlike
    // databoxes) which prevents double-frees and means we can freely
    // call finalize whenever we need to set state.
    if (status_ != DataStatus::Unmanaged) {
      if (status_ == DataStatus::AllocatedHost) {
        std::free(data_);
      } else if (status_ == DataStatus::AllocatedDevice) {
        PORTABLE_FREE(data_);
      }
      data_ = nullptr;
      n_ = 0;
      status_ = DataStatus::Empty;
    }
  }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    PORTABLE_REQUIRE(status_ != DataStatus::AllocatedDevice,
                     "Cannot save device-resident grid memory");
    const auto h5_type =
        std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    const hsize_t dims[] = {static_cast<hsize_t>(n_)};
    return H5LTmake_dataset(loc, name.c_str(), 1, dims, h5_type, data_);
  }

  inline herr_t loadHDF(hid_t loc, const std::string &name) {
    static_assert(
        std::is_same<T, double>::value || std::is_same<T, float>::value,
        "Spiner HDF5 only defined for these data types: float, double");
    int rank = 0;
    herr_t status = H5LTget_dataset_ndims(loc, name.c_str(), &rank);
    PORTABLE_ALWAYS_REQUIRE(status >= 0 && rank == 1,
                            "Nonuniform grid HDF5 data must be rank one");
    hsize_t dims[1];
    H5T_class_t class_id;
    size_t type_size;
    status +=
        H5LTget_dataset_info(loc, name.c_str(), dims, &class_id, &type_size);
    finalize();
    n_ = static_cast<std::size_t>(dims[0]);
    PORTABLE_ALWAYS_REQUIRE(n_ >= 2, "Valid nonuniform grid");
    data_ = static_cast<T *>(std::malloc(dynamicMemorySizeInBytes()));
    PORTABLE_ALWAYS_REQUIRE(data_ != nullptr, "Grid allocation failed");
    status_ = DataStatus::AllocatedHost;
    const auto h5_type =
        std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;
    status += H5LTread_dataset(loc, name.c_str(), h5_type, data_);
    validate_(PortsOfCall::Exec::Host{});
    return status;
  }
#endif

  // JMM: Public because otherwise Cuda complains. Not part of the
  // public grid API, but since it's a static helper function, doesn't
  // hurt.
  PORTABLE_INLINE_FUNCTION
  static bool pointsWellFormed(const T *points, const std::size_t n) {
    if (points == nullptr || n < 2) return false;
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(points[i])) return false;
      if (i > 0 && !(points[i - 1] < points[i])) return false;
    }
    return true;
  }

 private:
  void allocate_(const T *src, const std::size_t n) {
    pointsWellFormed(src, n);
    n_ = n;
    data_ = static_cast<T *>(std::malloc(dynamicMemorySizeInBytes()));
    PORTABLE_ALWAYS_REQUIRE(data_ != nullptr, "Grid allocation failed");
    std::memcpy(data_, src, dynamicMemorySizeInBytes());
    status_ = DataStatus::AllocatedHost;
  }

  template <typename ExecutionSpace>
  void validate_(const ExecutionSpace &E) const {
    bool valid = false;
    T *points = data_; // TODO(JMM): Expose PORTABLE_CLASS_LAMBDA
    std::size_t n = n_;
    portableReduce(
        "Validate UnstructuredGrid1D", E, 0, 1,
        PORTABLE_LAMBDA(const int, bool &b) {
          b = pointsWellFormed(points, n);
        },
        valid);
    PORTABLE_ALWAYS_REQUIRE(valid, "Dataset is well formed");
  }

  std::size_t n_ = 0;
  T *data_ = nullptr;
  DataStatus status_ = DataStatus::Empty;
};

} // namespace Spiner

#endif // SPINER_NONUNIFORM_GRID_1D_
