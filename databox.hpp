#ifndef _SPINER_DATABOX_HPP_
#define _SPINER_DATABOX_HPP_
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

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <assert.h>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

#include "interpolation.hpp"
#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_arrays.hpp"
#include "sp5.hpp"
#include "spiner_types.hpp"

// TODO: get named indices working
// TODO: If asserts are too slow, remove them.

namespace Spiner {

enum class IndexType { Interpolated = 0, Named = 1, Indexed = 2 };
enum class DataStatus { Empty, Unmanaged, AllocatedHost, AllocatedDevice };
enum class AllocationTarget { Host, Device };
constexpr int MAXRANK = PortableMDArray<Real>::MAXDIM;

class DataBox {
 public:
  // Base constructor
  DataBox() = default;

  // Rank constructors w/ pointer
  // args should be ints.
  // example call
  // DataBox(data, nx3, nx2, nx1)
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION __attribute__((nothrow))
  DataBox(Real *data, Args... args)
      : rank_(sizeof...(args)), status_(DataStatus::Unmanaged), data_(data) {
    dataView_.NewPortableMDArray(data, std::forward<Args>(args)...);
    setAllIndexed_();
  }

  // Rank constructors w/o pointer
  // args should be ints.
  // example call
  // DataBox(data, nx3, nx2, nx1)
  template <typename... Args>
  inline __attribute__((nothrow)) DataBox(AllocationTarget t, Args... args)
      : rank_(sizeof...(args)) {
    allocate_(t, std::forward<Args>(args)...);
    dataView_.NewPortableMDArray(data_, std::forward<Args>(args)...);
    setAllIndexed_();
  }
  template <typename... Args>
  inline __attribute__((nothrow)) DataBox(Args... args)
      : DataBox(AllocationTarget::Host, std::forward<Args>(args)...) {}

  // Copy and move constructors. All shallow.
  inline __attribute__((nothrow)) DataBox(PortableMDArray<Real> A)
      : rank_(A.GetRank()), status_(DataStatus::Unmanaged), data_(A.data()),
        dataView_(A) {
    setAllIndexed_();
  }
  inline __attribute__((nothrow)) DataBox(PortableMDArray<Real> &A)
      : rank_(A.GetRank()), status_(DataStatus::Unmanaged), data_(A.data()),
        dataView_(A) {
    setAllIndexed_();
  }
  PORTABLE_INLINE_FUNCTION __attribute__((nothrow)) DataBox(const DataBox &src)
      : rank_(src.rank_), status_(src.status_), data_(src.data_) {
    setAllIndexed_();
    dataView_.InitWithShallowSlice(src.dataView_, 6, 0, src.dim(6));
    for (int i = 0; i < rank_; i++) {
      indices_[i] = src.indices_[i];
      grids_[i] = src.grids_[i];
    }
  }

  // Slice constructor
  PORTABLE_INLINE_FUNCTION __attribute__((nothrow))
  DataBox(const DataBox &b, const int dim, const int indx, const int nvar)
      : status_(DataStatus::Unmanaged), data_(b.data_) {
    dataView_.InitWithShallowSlice(b.dataView_, dim, indx, nvar);
    rank_ = dataView_.GetRank();
    for (int i = 0; i < rank_; i++) {
      indices_[i] = b.indices_[i];
      grids_[i] = b.grids_[i];
    }
  }

#ifdef SPINER_USE_HDF
  // HDF5 constructors
  inline DataBox(const std::string &filename) { loadHDF(filename); }
  inline DataBox(hid_t loc, const std::string &groupname) {
    loadHDF(loc, groupname);
  }
#endif // SPINER_USE_HDF

  // Read an array, shallow
  inline void setArray(PortableMDArray<Real> &A);

  // Re-allocates memory for DataBox either on host or device.
  // This is destructive. Memory is freed!
  template <typename... Args>
  inline void resize(AllocationTarget t, Args... args);
  template <typename... Args>
  inline void resize(Args... args) {
    resize(AllocationTarget::Host, std::forward<Args>(args)...);
  }

  // Index operators
  // examle calls:
  // Real x = db(n4, n3, n2, n1);
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION Real &operator()(Args... args) {
    return dataView_(std::forward<Args>(args)...);
  }
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION Real &operator()(Args... args) const {
    return dataView_(std::forward<Args>(args)...);
  }

  // Slice operation
  PORTABLE_INLINE_FUNCTION
  DataBox slice(const int dim, const int indx, const int nvar) const {
    return DataBox(*this, dim, indx, nvar);
  }
  PORTABLE_INLINE_FUNCTION DataBox slice(const int indx) const {
    return slice(rank_, indx, 1);
  }
  PORTABLE_INLINE_FUNCTION DataBox slice(const int ix2, const int ix1) const {
    // DataBox a(*this, rank_, ix2, 1);
    // return DataBox(a, a.rank_, ix1, 1);
    return slice(ix2).slice(ix1);
  }

  // Reshape the view of the array
  // example call:
  // db.reshape(n2,n1);
  // for DataBox db with total size shape n2*n1;
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION void reshape(Args... args) {
    dataView_.Reshape(std::forward<Args>(args)...);
  }

  // Interpolates whole DataBox to a real number,
  // x1 is fastest index. xN is slowest.
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline)) interpToReal(const Real x) const;
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline))
  interpToReal(const Real x2, const Real x1) const;
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline))
  interpToReal(const Real x3, const Real x2, const Real x1) const;
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline))
  interpToReal(const Real x3, const Real x2, const Real x1,
               const int idx) const;
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline))
  interpToReal(const Real x4, const Real x3, const Real x2,
               const Real x1) const;
  // Interpolates the whole databox to a real number,
  // with one intermediate, non-interpolatable index,
  // which is simply indexed into
  // JMM: Trust me---this is a common pattern
  PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
  __attribute__((always_inline))
  interpToReal(const Real x4, const Real x3, const Real x2, const int idx,
               const Real x1) const;
  // Interpolates SLOWEST indices of databox to a new
  // DataBox, interpolated at that slowest index.
  // WARNING: requires memory to be pre-allocated.
  // TODO: add 3d and higher interpFromDB if necessary
  PORTABLE_INLINE_FUNCTION void interpFromDB(const DataBox &db, const Real x);
  PORTABLE_INLINE_FUNCTION void interpFromDB(const DataBox &db, const Real x2,
                                             const Real x1);
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION DataBox interpToDB(Args... args) {
    DataBox db;
    db.interpFromDB(*this, std::forward<Args>(args)...);
    return db;
  }

  // Setters
  // NOTE: i ranges from 0 to N-1, where 0 is the FASTEST moving
  // index, and N-1 is the slowest
  // TODO: is this intuitive?
  inline void setIndexType(int i, IndexType t) {
    assert(0 <= i && i < rank_);
    indices_[i] = t;
  }
  inline void setRange(int i, RegularGrid1D g) {
    setIndexType(i, IndexType::Interpolated);
    grids_[i] = g;
  }
  inline void setRange(int i, Real xmin, Real xmax, int N) {
    setRange(i, RegularGrid1D(xmin, xmax, N));
  }

  // Reshapes from other databox, but does not allocate memory.
  // Does no checks that memory is available.
  // Optionally copies shape of source with ndims fewer slowest-moving
  // dimensions
  PORTABLE_INLINE_FUNCTION void copyShape(const DataBox &db,
                                          const int ndims = 0);
  // Returns new databox with same memory and metadata
  inline void copyMetadata(const DataBox &src);

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF() const { return saveHDF(SP5::DB::FILENAME); }
  inline herr_t saveHDF(const std::string &filename) const;
  inline herr_t saveHDF(hid_t loc, const std::string &groupname) const;
  inline herr_t loadHDF() { return loadHDF(SP5::DB::FILENAME); }
  inline herr_t loadHDF(const std::string &filename);
  inline herr_t loadHDF(hid_t loc, const std::string &groupname);
#endif

  // Reference accessors
  inline IndexType &indexType(const int i) { return indices_[i]; }
  inline RegularGrid1D &range(const int i) { return grids_[i]; }

  // Assignment and move, both perform shallow copy
  PORTABLE_INLINE_FUNCTION DataBox &operator=(const DataBox &other);
  inline void copy(const DataBox &src);

  // utility info
  inline DataStatus dataStatus() const { return status_; }
  inline bool isReference() { return status_ == DataStatus::Unmanaged; }
  inline bool ownsAllocatedMemory() {
    return (status_ != DataStatus::Unmanaged);
  }
  inline bool operator==(const DataBox &other) const;
  inline bool operator!=(const DataBox &other) const {
    return !(*this == other);
  }

  // call this to tell a databox it no longer owns its memory
  PORTABLE_INLINE_FUNCTION void makeShallow() {
    status_ = DataStatus::Unmanaged;
  }
  // call this to reset a databox to the default constructed state
  // this decouples it from other databoxes that might point at the same memory
  PORTABLE_INLINE_FUNCTION void reset() {
    data_ = nullptr;
    status_ = DataStatus::Empty;
    rank_ = 0;
  }

  PORTABLE_INLINE_FUNCTION Real *data() const { return data_; }
  PORTABLE_INLINE_FUNCTION Real min() const;
  PORTABLE_INLINE_FUNCTION Real max() const;
  PORTABLE_INLINE_FUNCTION int rank() const { return rank_; }
  PORTABLE_INLINE_FUNCTION int size() const { return dataView_.GetSize(); }
  PORTABLE_INLINE_FUNCTION int sizeBytes() const {
    return dataView_.GetSizeInBytes();
  }
  PORTABLE_INLINE_FUNCTION int dim(int i) const { return dataView_.GetDim(i); }
  PORTABLE_INLINE_FUNCTION void range(int i, Real &min, Real &max, Real &dx,
                                      int &N) const;
  PORTABLE_INLINE_FUNCTION RegularGrid1D range(int i) const;
  PORTABLE_INLINE_FUNCTION IndexType indexType(const int i) const {
    return indices_[i];
  }

  // TODO(JMM): Add more code for more portability strategies
  DataBox getOnDevice() const { // getOnDevice is always a deep copy
#ifdef PORTABILITY_STRATEGY_KOKKOS
    using HS = Kokkos::HostSpace;
    using DMS = Kokkos::DefaultExecutionSpace::memory_space;
    constexpr const bool execution_is_host{
        Kokkos::SpaceAccessibility<DMS, HS>::accessible};
    if (execution_is_host) {
      DataBox a;
      a.copy(*this); // a.copy handles setting allocation status
      return a;
    } else {
      using memUnmanaged = Kokkos::MemoryUnmanaged;
      using HostView_t = Kokkos::View<Real *, HS, memUnmanaged>;
      using DeviceView_t = Kokkos::View<Real *, memUnmanaged>;
      using Kokkos::deep_copy;
      Real *device_data = (Real *)PORTABLE_MALLOC(sizeBytes());
      DeviceView_t devView(device_data, dataView_.GetSize());
      HostView_t hostView(data_, dataView_.GetSize());
      deep_copy(devView, hostView);
      DataBox a{devView.data(), dim(6), dim(5), dim(4), dim(3), dim(2), dim(1)};
      a.copyShape(*this);
      a.status_ = DataStatus::AllocatedDevice;
      return a;
    }
#else  // no kokkos
    DataBox a;
    a.copy(*this); // a.copy handles allocation status.
    return a;
#endif // kokkos
  }

  // TODO(JMM): Potentially use this for device-free
  void finalize() {
    assert(ownsAllocatedMemory());
    if (status_ == DataStatus::AllocatedDevice) {
      PORTABLE_FREE(data_);
    } else if (status_ == DataStatus::AllocatedHost) {
      free(data_);
    }
    status_ = DataStatus::Empty;
  }

 private:
  int rank_ = 0; // after dataView_ b/c dataView_ should be initialized first
  DataStatus status_ = DataStatus::Empty;
  // when we manage our own data on host, it lives here
  Real *data_ = nullptr;           // points at data, managed or not
  PortableMDArray<Real> dataView_; // always used
  IndexType indices_[MAXRANK];
  RegularGrid1D grids_[MAXRANK];

  PORTABLE_INLINE_FUNCTION void setAllIndexed_();
  PORTABLE_INLINE_FUNCTION bool canInterpToReal_(const int interpOrder) const;
  inline std::string gridname_(int i) const {
    return SP5::DB::GRID_FORMAT[0] + std::to_string(i + 1) +
           SP5::DB::GRID_FORMAT[1];
  }

  int prod_() { return 1; }
  int prod_(int i) { return i; }
  template <typename... Rest>
  int prod_(int i, Rest... rest) {
    return i * prod_(std::forward<Rest>(rest)...);
  }
  template <typename... Args>
  inline void allocate_(AllocationTarget t, Args... args) {
    finalize();
    size_t size = sizeof(Real) * prod_(std::forward<Args>(args)...);
    if (t == AllocationTarget::Device) {
      data_ = (Real *)PORTABLE_MALLOC(size);
      status_ = DataStatus::AllocatedDevice;
    } else {
      data_ = (Real *)malloc(size);
      status_ = DataStatus::AllocatedHost;
    }
  }
};

// Read an array, shallow
inline void DataBox::setArray(PortableMDArray<Real> &A) {
  dataView_ = A;
  rank_ = A.GetRank();
  setAllIndexed_();
}

// Allocate memory of constructed DataBox
template <typename... Args>
inline void DataBox::resize(AllocationTarget t, Args... args) {
  assert(ownsAllocatedMemory());
  rank_ = sizeof...(args);
  allocate_(t, std::forward<Args>(args)...);
  setAllIndexed_();
  dataView_.NewPortableMDArray(data_, std::forward<Args>(args)...);
}

PORTABLE_INLINE_FUNCTION Real DataBox::interpToReal(const Real x) const {
  assert(canInterpToReal_(1));
  return grids_[0](x, dataView_);
}

PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
__attribute__((always_inline))
DataBox::interpToReal(const Real x2, const Real x1) const {
  assert(canInterpToReal_(2));
  int ix1, ix2;
  weights_t w1, w2;
  grids_[0].weights(x1, ix1, w1);
  grids_[1].weights(x2, ix2, w2);
  // TODO: prefectch corners for speed?
  // TODO: re-order access pattern?
  return (w2[0] *
              (w1[0] * dataView_(ix2, ix1) + w1[1] * dataView_(ix2, ix1 + 1)) +
          w2[1] * (w1[0] * dataView_(ix2 + 1, ix1) +
                   w1[1] * dataView_(ix2 + 1, ix1 + 1)));
}

PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
__attribute__((always_inline))
DataBox::interpToReal(const Real x3, const Real x2, const Real x1) const {
  assert(canInterpToReal_(3));
  int ix[3];
  weights_t w[3];
  grids_[0].weights(x1, ix[0], w[0]);
  grids_[1].weights(x2, ix[1], w[1]);
  grids_[2].weights(x3, ix[2], w[2]);
  // TODO: prefect corners for speed?
  // TODO: re-order access pattern?
  return (
      w[2][0] * (w[1][0] * (w[0][0] * dataView_(ix[2], ix[1], ix[0]) +
                            w[0][1] * dataView_(ix[2], ix[1], ix[0] + 1)) +
                 w[1][1] * (w[0][0] * dataView_(ix[2], ix[1] + 1, ix[0]) +
                            w[0][1] * dataView_(ix[2], ix[1] + 1, ix[0] + 1))) +
      w[2][1] *
          (w[1][0] * (w[0][0] * dataView_(ix[2] + 1, ix[1], ix[0]) +
                      w[0][1] * dataView_(ix[2] + 1, ix[1], ix[0] + 1)) +
           w[1][1] * (w[0][0] * dataView_(ix[2] + 1, ix[1] + 1, ix[0]) +
                      w[0][1] * dataView_(ix[2] + 1, ix[1] + 1, ix[0] + 1))));
}

PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
__attribute__((always_inline))
DataBox::interpToReal(const Real x3, const Real x2, const Real x1,
                      const int idx) const {
  assert(rank_ == 4);
  for (int r = 1; r < rank_; ++r) {
    assert(indices_[r] == IndexType::Interpolated);
    assert(grids_[r].isWellFormed());
  }
  int ix[3];
  weights_t w[3];
  grids_[1].weights(x1, ix[0], w[0]);
  grids_[2].weights(x2, ix[1], w[1]);
  grids_[3].weights(x3, ix[2], w[2]);
  // TODO: prefect corners for speed?
  // TODO: re-order access pattern?
  return (
      w[2][0] *
          (w[1][0] * (w[0][0] * dataView_(ix[2], ix[1], ix[0], idx) +
                      w[0][1] * dataView_(ix[2], ix[1], ix[0] + 1, idx)) +
           w[1][1] * (w[0][0] * dataView_(ix[2], ix[1] + 1, ix[0], idx) +
                      w[0][1] * dataView_(ix[2], ix[1] + 1, ix[0] + 1, idx))) +
      w[2][1] *
          (w[1][0] * (w[0][0] * dataView_(ix[2] + 1, ix[1], ix[0], idx) +
                      w[0][1] * dataView_(ix[2] + 1, ix[1], ix[0] + 1, idx)) +
           w[1][1] *
               (w[0][0] * dataView_(ix[2] + 1, ix[1] + 1, ix[0], idx) +
                w[0][1] * dataView_(ix[2] + 1, ix[1] + 1, ix[0] + 1, idx))));
}

PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
__attribute__((always_inline))
DataBox::interpToReal(const Real x4, const Real x3, const Real x2,
                      const Real x1) const {
  assert(canInterpToReal_(4));
  Real x[] = {x1, x2, x3, x4};
  int ix[4];
  weights_t w[4];
  for (int i = 0; i < 4; ++i) {
    grids_[i].weights(x[i], ix[i], w[i]);
  }
  // TODO(JMM): This is getty pretty gross. Should we automate?
  // Hand-written is probably faster, though.
  // Breaking line-limit to make this easier to read
  return (
      w[3][0] *
          (w[2][0] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3], ix[2], ix[1], ix[0]) +
                     w[0][1] * dataView_(ix[3], ix[2], ix[1], ix[0] + 1)) +
                w[1][1] *
                    (w[0][0] * dataView_(ix[3], ix[2], ix[1] + 1, ix[0]) +
                     w[0][1] * dataView_(ix[3], ix[2], ix[1] + 1, ix[0] + 1))) +
           w[2][1] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3], ix[2] + 1, ix[1], ix[0]) +
                     w[0][1] * dataView_(ix[3], ix[2] + 1, ix[1], ix[0] + 1)) +
                w[1][1] *
                    (w[0][0] * dataView_(ix[3], ix[2] + 1, ix[1] + 1, ix[0]) +
                     w[0][1] *
                         dataView_(ix[3], ix[2] + 1, ix[1] + 1, ix[0] + 1)))) +
      w[3][1] *
          (w[2][0] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3] + 1, ix[2], ix[1], ix[0]) +
                     w[0][1] * dataView_(ix[3] + 1, ix[2], ix[1], ix[0] + 1)) +
                w[1][1] *
                    (w[0][0] * dataView_(ix[3] + 1, ix[2], ix[1] + 1, ix[0]) +
                     w[0][1] *
                         dataView_(ix[3] + 1, ix[2], ix[1] + 1, ix[0] + 1))) +
           w[2][1] * (w[1][0] * (w[0][0] * dataView_(ix[3] + 1, ix[2] + 1,
                                                     ix[1], ix[0]) +
                                 w[0][1] * dataView_(ix[3] + 1, ix[2] + 1,
                                                     ix[1], ix[0] + 1)) +
                      w[1][1] * (w[0][0] * dataView_(ix[3] + 1, ix[2] + 1,
                                                     ix[1] + 1, ix[0]) +
                                 w[0][1] * dataView_(ix[3] + 1, ix[2] + 1,
                                                     ix[1] + 1, ix[0] + 1))))

  );
}

PORTABLE_INLINE_FUNCTION Real __attribute__((nothrow))
__attribute__((always_inline))
DataBox::interpToReal(const Real x4, const Real x3, const Real x2,
                      const int idx, const Real x1) const {
  assert(rank_ == 5);
  assert(indices_[0] == IndexType::Interpolated);
  assert(grids_[0].isWellFormed());
  for (int i = 2; i < 5; ++i) {
    assert(indices_[i] == IndexType::Interpolated);
    assert(grids_[i].isWellFormed());
  }
  Real x[] = {x1, x2, x3, x4};
  int ix[4];
  weights_t w[4];
  grids_[0].weights(x[0], ix[0], w[0]);
  for (int i = 1; i < 4; ++i) {
    grids_[i + 1].weights(x[i], ix[i], w[i]);
  }
  // TODO(JMM): This is getty pretty gross. Should we automate?
  // Hand-written is probably faster, though.
  // Breaking line-limit to make this easier to read
  return (
      w[3][0] *
          (w[2][0] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3], ix[2], ix[1], idx, ix[0]) +
                     w[0][1] * dataView_(ix[3], ix[2], ix[1], idx, ix[0] + 1)) +
                w[1][1] *
                    (w[0][0] * dataView_(ix[3], ix[2], ix[1] + 1, idx, ix[0]) +
                     w[0][1] *
                         dataView_(ix[3], ix[2], ix[1] + 1, idx, ix[0] + 1))) +
           w[2][1] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3], ix[2] + 1, ix[1], idx, ix[0]) +
                     w[0][1] *
                         dataView_(ix[3], ix[2] + 1, ix[1], idx, ix[0] + 1)) +
                w[1][1] * (w[0][0] * dataView_(ix[3], ix[2] + 1, ix[1] + 1, idx,
                                               ix[0]) +
                           w[0][1] * dataView_(ix[3], ix[2] + 1, ix[1] + 1, idx,
                                               ix[0] + 1)))) +
      w[3][1] *
          (w[2][0] *
               (w[1][0] *
                    (w[0][0] * dataView_(ix[3] + 1, ix[2], ix[1], idx, ix[0]) +
                     w[0][1] *
                         dataView_(ix[3] + 1, ix[2], ix[1], idx, ix[0] + 1)) +
                w[1][1] * (w[0][0] * dataView_(ix[3] + 1, ix[2], ix[1] + 1, idx,
                                               ix[0]) +
                           w[0][1] * dataView_(ix[3] + 1, ix[2], ix[1] + 1, idx,
                                               ix[0] + 1))) +
           w[2][1] *
               (w[1][0] * (w[0][0] * dataView_(ix[3] + 1, ix[2] + 1, ix[1], idx,
                                               ix[0]) +
                           w[0][1] * dataView_(ix[3] + 1, ix[2] + 1, ix[1], idx,
                                               ix[0] + 1)) +
                w[1][1] * (w[0][0] * dataView_(ix[3] + 1, ix[2] + 1, ix[1] + 1,
                                               idx, ix[0]) +
                           w[0][1] * dataView_(ix[3] + 1, ix[2] + 1, ix[1] + 1,
                                               idx, ix[0] + 1))))

  );
}

PORTABLE_INLINE_FUNCTION void DataBox::interpFromDB(const DataBox &db,
                                                    const Real x) {
  assert(db.indices_[db.rank_ - 1] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 1].isWellFormed());
  assert(size() == (db.size() / db.dim(db.rank_)));

  int ix;
  weights_t w;
  copyShape(db, 1);

  db.grids_[db.rank_ - 1].weights(x, ix, w);
  DataBox lower(db.slice(ix)), upper(db.slice(ix + 1));
  // lower = db.slice(ix);
  // upper = db.slice(ix+1);
  for (int i = 0; i < size(); i++) {
    dataView_(i) = w[0] * lower(i) + w[1] * upper(i);
  }
}

PORTABLE_INLINE_FUNCTION void
DataBox::interpFromDB(const DataBox &db, const Real x2, const Real x1) {
  assert(db.rank_ >= 2);
  assert(db.indices_[db.rank_ - 1] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 1].isWellFormed());
  assert(db.indices_[db.rank_ - 2] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 2].isWellFormed());
  assert(size() == (db.size() / (db.dim(db.rank_) * db.dim(db.rank_ - 1))));

  int ix2, ix1;
  weights_t w2, w1;
  copyShape(db, 2);

  db.grids_[db.rank_ - 2].weights(x1, ix1, w1);
  db.grids_[db.rank_ - 1].weights(x2, ix2, w2);
  DataBox corners[2][2]{{db.slice(ix2, ix1), db.slice(ix2 + 1, ix1)},
                        {db.slice(ix2, ix1 + 1), db.slice(ix2 + 1, ix1 + 1)}};
  //    copyShape(db,2);
  //
  //    db.grids_[db.rank_-2].weights(x1, ix1, w1);
  //    db.grids_[db.rank_-1].weights(x2, ix2, w2);
  // corners[0][0] = db.slice(ix2,   ix1   );
  // corners[1][0] = db.slice(ix2,   ix1+1 );
  // corners[0][1] = db.slice(ix2+1, ix1   );
  // corners[1][1] = db.slice(ix2+1, ix1+1 );
  /*
  for (int i = 0; i < size(); i++) {
    dataView_(i) = (   w2[0]*w1[0]*corners[0][0](i)
                     + w2[0]*w1[1]*corners[1][0](i)
                     + w2[1]*w1[0]*corners[0][1](i)
                     + w2[1]*w1[1]*corners[1][1](i));
  }
  */
  for (int i = 0; i < size(); i++) {
    dataView_(i) =
        (w2[0] * (w1[0] * corners[0][0](i) + w1[1] * corners[1][0](i)) +
         w2[1] * (w1[0] * corners[0][1](i) + w1[1] * corners[1][1](i)));
  }
}

// Reshapes from other databox, but does not allocate memory.
// Does no checks that memory is available.
// Optionally copies shape of source with ndims fewer slowest-moving dimensions
PORTABLE_INLINE_FUNCTION void DataBox::copyShape(const DataBox &db,
                                                 const int ndims) {
  rank_ = db.rank_ - ndims;
  int dims[MAXRANK];
  for (int i = 0; i < MAXRANK; i++)
    dims[i] = 1;
  setAllIndexed_(); // TODO: can remove
  for (int i = rank_ - 1; i >= 0; i--)
    dims[i] = db.dim(i + 1);
  reshape(dims[5], dims[4], dims[3], dims[2], dims[1], dims[0]);
  rank_ = db.rank_ - ndims;
  for (int i = 0; i < rank_; i++) {
    indices_[i] = db.indices_[i];
    grids_[i] = db.grids_[i];
  }
}
// reallocates and then copies shape from other databox
// everything but the actual copy in a deep copy
inline void DataBox::copyMetadata(const DataBox &src) {
  AllocationTarget t =
      (src.status_ == DataStatus::AllocatedDevice ? AllocationTarget::Device
                                                  : AllocationTarget::Host);
  resize(t, src.dim(6), src.dim(5), src.dim(4), src.dim(3), src.dim(2),
         src.dim(1));
  rank_ = src.rank_; // resize sets rank
  for (int i = 0; i < rank_; i++) {
    grids_[i] = src.grids_[i];
    indices_[i] = src.indices_[i];
  }
}

#ifdef SPINER_USE_HDF
inline herr_t DataBox::saveHDF(const std::string &filename) const {
  herr_t status;
  hid_t file;

  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = saveHDF(file, SP5::DB::GRPNAME);
  status += H5Fclose(file);
  return status;
}

inline herr_t DataBox::saveHDF(hid_t loc, const std::string &groupname) const {
  hid_t group, grids;
  herr_t status = 0;

  std::vector<int> dims_int(rank_);
  for (int i = 0; i < rank_; i++)
    dims_int[i] = dim(i + 1);

  std::vector<hsize_t> dims_hsize_t(rank_);
  for (int i = 0; i < rank_; i++)
    dims_hsize_t[i] = dim(rank_ - i);

  std::vector<int> index_types(rank_);
  for (int i = 0; i < rank_; i++) {
    index_types[i] = static_cast<int>(indices_[i]);
  }

  // Greate group
  group =
      H5Gcreate(loc, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Record rank as an attribute
  status += H5LTset_attribute_int(loc, groupname.c_str(), SP5::DB::RANKNAME,
                                  &rank_, 1);

  // Index types
  status += H5LTset_attribute_int(loc, groupname.c_str(), SP5::DB::IDXSNAME,
                                  index_types.data(), rank_);
  status += H5LTset_attribute_string(loc, groupname.c_str(),
                                     SP5::DB::IDXINFONAME, SP5::DB::IDXINFO);

  // Dimensions of the PortableMDArray, set as an attribute
  status += H5LTset_attribute_int(loc, groupname.c_str(), SP5::DB::DIMSNAME,
                                  dims_int.data(), rank_);

  // Save the PortableMDArray, in a human-readable shape
  status += H5LTmake_dataset(group, SP5::DB::DSETNAME, rank_,
                             dims_hsize_t.data(), H5T_REAL, dataView_.data());

  // Grids
  grids = H5Gcreate(group, SP5::DB::GRIDNAME, H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  for (int i = 0; i < rank_; i++) {
    if (indices_[i] == IndexType::Interpolated) {

      status += grids_[i].saveHDF(grids, gridname_(i).c_str());
    }
  }
  status += H5Gclose(grids);
  status += H5Gclose(group);
  return status;
}

inline herr_t DataBox::loadHDF(const std::string &filename) {
  herr_t status;
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  status = loadHDF(file, SP5::DB::GRPNAME);
  return status;
}

inline herr_t DataBox::loadHDF(hid_t loc, const std::string &groupname) {
  hid_t group, grids;
  herr_t status = 0;
  std::vector<int> index_types;
  std::vector<int> dims(6, 1);

  // Open group
  group = H5Gopen(loc, groupname.c_str(), H5P_DEFAULT);
  // Get rank
  status +=
      H5LTget_attribute_int(loc, groupname.c_str(), SP5::DB::RANKNAME, &rank_);
  // Resize metadata fields
  setAllIndexed_();

  // Get dimensions
  dims.resize(rank_);
  status += H5LTget_attribute_int(loc, groupname.c_str(), SP5::DB::DIMSNAME,
                                  dims.data());

  // Allocate PortableMDArray and read it in to buffer
  allocate_(AllocationTarget::Host, dims[5], dims[4], dims[3], dims[2], dims[1],
            dims[0]);
  dataView_.NewPortableMDArray(data_, dims[5], dims[4], dims[3], dims[2],
                               dims[1], dims[0]);
  status +=
      H5LTread_dataset(group, SP5::DB::DSETNAME, H5T_REAL, dataView_.data());

  // Get index types
  index_types.resize(rank_);
  status += H5LTget_attribute_int(loc, groupname.c_str(), SP5::DB::IDXSNAME,
                                  index_types.data());
  for (int i = 0; i < rank_; i++) {
    indices_[i] = static_cast<IndexType>(index_types[i]);
  }

  // Get grids
  grids = H5Gopen(group, SP5::DB::GRIDNAME, H5P_DEFAULT);
  for (int i = 0; i < rank_; i++) {
    if (indices_[i] == IndexType::Interpolated) {
      status += grids_[i].loadHDF(grids, gridname_(i).c_str());
    }
  }
  status += H5Gclose(grids);

  status += H5Gclose(group);
  return status;
}
#endif // SPINER_USE_HDF

// Performs shallow copy by default
PORTABLE_INLINE_FUNCTION DataBox &DataBox::operator=(const DataBox &src) {
  if (this != &src) {
    rank_ = src.rank_;
    status_ = src.status_;
    data_ = src.data_;
    dataView_.InitWithShallowSlice(src.dataView_, 6, 0, src.dim(6));
    for (int i = 0; i < rank_; i++) {
      indices_[i] = src.indices_[i];
      grids_[i] = src.grids_[i];
    }
  }
  return *this;
}

// Performs a deep copy
inline void DataBox::copy(const DataBox &src) {
  copyMetadata(src);
  for (int i = 0; i < src.size(); i++)
    dataView_(i) = src(i);
}

inline bool DataBox::operator==(const DataBox &other) const {
  if (rank_ != other.rank_) return false;
  for (int i = 0; i < rank_; i++) {
    if (indices_[i] != other.indices_[i]) return false;
    if (indices_[i] == IndexType::Interpolated &&
        grids_[i] != other.grids_[i]) {
      return false;
    }
  }
  return dataView_ == other.dataView_;
}

// TODO: should this be std::reduce?
PORTABLE_INLINE_FUNCTION Real DataBox::min() const {
  Real min = std::numeric_limits<Real>::infinity();
  for (int i = 0; i < size(); i++) {
    min = std::min(min, dataView_(i));
  }
  return min;
}

PORTABLE_INLINE_FUNCTION Real DataBox::max() const {
  Real max = -std::numeric_limits<Real>::infinity();
  for (int i = 0; i < size(); i++) {
    max = std::max(max, dataView_(i));
  }
  return max;
}

PORTABLE_INLINE_FUNCTION void DataBox::range(int i, Real &min, Real &max,
                                             Real &dx, int &N) const {
  assert(0 <= i && i < rank_);
  assert(indices_[i] == IndexType::Interpolated);
  min = grids_[i].min();
  max = grids_[i].max();
  dx = grids_[i].dx();
  N = grids_[i].nPoints();
}

PORTABLE_INLINE_FUNCTION RegularGrid1D DataBox::range(int i) const {
  assert(0 <= i && i < rank_);
  assert(indices_[i] == IndexType::Interpolated);
  return grids_[i];
}

PORTABLE_INLINE_FUNCTION void DataBox::setAllIndexed_() {
  for (int i = 0; i < rank_; i++) {
    indices_[i] = IndexType::Indexed;
  }
}

PORTABLE_INLINE_FUNCTION bool
DataBox::canInterpToReal_(const int interpOrder) const {
  if (rank_ != interpOrder) return false;
  for (int i = 0; i < rank_; i++) {
    if (indices_[i] != IndexType::Interpolated) return false;
    if (!(grids_[i].isWellFormed())) return false;
  }
  return true;
}

inline DataBox getOnDeviceDataBox(const DataBox &a_host) {
  return a_host.getOnDevice();
}
inline void free(DataBox &db) { db.finalize(); }

struct DBDeleter {
  template <typename T>
  void operator()(T *ptr) {
    ptr->finalize();
    delete ptr;
  }
};
} // namespace Spiner
#endif // _SPINER_DATABOX_HPP_
