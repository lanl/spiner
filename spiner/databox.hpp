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
enum class DataStatus { Empty = 0, Unmanaged = 1, AllocatedHost = 2, AllocatedDevice = 3 };
enum class AllocationTarget { Host, Device };

template <typename T = Real, typename Grid_t = RegularGrid1D<T>,
          typename Concept =
              typename std::enable_if<std::is_arithmetic<T>::value, bool>::type>
class DataBox {
 public:
  using ValueType = T;
  using GridType = Grid_t;
  static constexpr int MAXRANK = PortableMDArray<T>::MAXDIM;
  static constexpr T EPS = 10.0 * std::numeric_limits<T>::epsilon();

  // Base constructor
  DataBox() = default;

  // Rank constructors w/ pointer
  // args should be ints.
  // example call
  // DataBox(data, nx3, nx2, nx1)
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION DataBox(T *data, Args... args) noexcept
      : rank_(sizeof...(args)), status_(DataStatus::Unmanaged), data_(data) {
    dataView_.NewPortableMDArray(data, std::forward<Args>(args)...);
    setAllIndexed_();
  }

  // Rank constructors w/o pointer
  // args should be ints.
  // example call
  // DataBox(data, nx3, nx2, nx1)
  template <typename... Args>
  inline DataBox(AllocationTarget t, Args... args) noexcept
      : rank_(sizeof...(args)) {
    allocate_(t, std::forward<Args>(args)...);
    dataView_.NewPortableMDArray(data_, std::forward<Args>(args)...);
    setAllIndexed_();
  }
  template <typename... Args>
  inline DataBox(Args... args) noexcept
      : DataBox(AllocationTarget::Host, std::forward<Args>(args)...) {}

  // Copy and move constructors. All shallow.
  inline DataBox(PortableMDArray<T> A) noexcept
      : rank_(A.GetRank()), status_(DataStatus::Unmanaged), data_(A.data()),
        dataView_(A) {
    setAllIndexed_();
  }
  inline DataBox(PortableMDArray<T> &A) noexcept
      : rank_(A.GetRank()), status_(DataStatus::Unmanaged), data_(A.data()),
        dataView_(A) {
    setAllIndexed_();
  }
  PORTABLE_INLINE_FUNCTION
  DataBox(const DataBox<T, Grid_t, Concept> &src) noexcept
      : rank_(src.rank_), status_(src.status_), data_(src.data_) {
    setAllIndexed_();
    dataView_.InitWithShallowSlice(src.dataView_, 6, 0, src.dim(6));
    for (int i = 0; i < rank_; i++) {
      indices_[i] = src.indices_[i];
      grids_[i] = src.grids_[i];
    }
  }

  // Slice constructor
  PORTABLE_INLINE_FUNCTION
  DataBox(const DataBox<T, Grid_t, Concept> &b, const int dim, const int indx,
          const int nvar) noexcept
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
  inline void setArray(PortableMDArray<T> &A);

  // Re-allocates memory for DataBox either on host or device.
  // This is destructive. Memory is freed!
  template <typename... Args>
  inline void resize(AllocationTarget t, Args... args) {
    assert(ownsAllocatedMemory());
    rank_ = sizeof...(args);
    allocate_(t, std::forward<Args>(args)...);
    setAllIndexed_();
    dataView_.NewPortableMDArray(data_, std::forward<Args>(args)...);
  }
  template <typename... Args>
  inline void resize(Args... args) {
    resize(AllocationTarget::Host, std::forward<Args>(args)...);
  }

  // Index operators
  // examle calls:
  // T x = db(n4, n3, n2, n1);
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION T &operator()(Args... args) {
    return dataView_(std::forward<Args>(args)...);
  }
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION T &operator()(Args... args) const {
    return dataView_(std::forward<Args>(args)...);
  }

  // Slice operation
  PORTABLE_INLINE_FUNCTION
  DataBox<T, Grid_t, Concept> slice(const int dim, const int indx,
                                    const int nvar) const {
    return DataBox(*this, dim, indx, nvar);
  }
  PORTABLE_INLINE_FUNCTION DataBox<T, Grid_t, Concept>
  slice(const int indx) const {
    return slice(rank_, indx, 1);
  }
  PORTABLE_INLINE_FUNCTION DataBox<T, Grid_t, Concept>
  slice(const int ix2, const int ix1) const {
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
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x) const noexcept;
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x2,
                                               const T x1) const noexcept;
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x3, const T x2,
                                               const T x1) const noexcept;
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x3, const T x2,
                                               const T x1,
                                               const int idx) const noexcept;
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x4, const T x3,
                                               const T x2,
                                               const T x1) const noexcept;
  // Interpolates the whole databox to a real number,
  // with one intermediate, non-interpolatable index,
  // which is simply indexed into
  // JMM: Trust me---this is a common pattern
  PORTABLE_FORCEINLINE_FUNCTION T interpToReal(const T x4, const T x3,
                                               const T x2, const int idx,
                                               const T x1) const noexcept;
  // Interpolates SLOWEST indices of databox to a new
  // DataBox, interpolated at that slowest index.
  // WARNING: requires memory to be pre-allocated.
  // TODO: add 3d and higher interpFromDB if necessary
  PORTABLE_INLINE_FUNCTION void
  interpFromDB(const DataBox<T, Grid_t, Concept> &db, const T x);
  PORTABLE_INLINE_FUNCTION void
  interpFromDB(const DataBox<T, Grid_t, Concept> &db, const T x2, const T x1);
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION DataBox<T, Grid_t, Concept>
  interpToDB(Args... args) {
    DataBox<T, Grid_t, Concept> db;
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
  inline void setRange(int i, Grid_t g) {
    setIndexType(i, IndexType::Interpolated);
    grids_[i] = g;
  }
  template <typename... Args>
  inline void setRange(int i, Args &&...args) {
    setRange(i, Grid_t(std::forward<Args>(args)...));
  }

  // Reshapes from other databox, but does not allocate memory.
  // Does no checks that memory is available.
  // Optionally copies shape of source with ndims fewer slowest-moving
  // dimensions
  PORTABLE_INLINE_FUNCTION void copyShape(const DataBox<T, Grid_t, Concept> &db,
                                          const int ndims = 0);
  // Returns new databox with same memory and metadata
  inline void copyMetadata(const DataBox<T, Grid_t, Concept> &src);

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
  inline Grid_t &range(const int i) { return grids_[i]; }

  // Assignment and move, both perform shallow copy
  PORTABLE_INLINE_FUNCTION DataBox<T, Grid_t, Concept> &
  operator=(const DataBox<T, Grid_t, Concept> &other);
  inline void copy(const DataBox<T, Grid_t, Concept> &src);

  // utility info
  inline DataStatus dataStatus() const { return status_; }
  inline bool isReference() { return status_ == DataStatus::Unmanaged; }
  inline bool ownsAllocatedMemory() {
    return (status_ != DataStatus::Unmanaged);
  }
  inline bool operator==(const DataBox<T, Grid_t, Concept> &other) const;
  inline bool operator!=(const DataBox<T, Grid_t, Concept> &other) const {
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

  PORTABLE_INLINE_FUNCTION T *data() const { return data_; }
  PORTABLE_INLINE_FUNCTION T min() const;
  PORTABLE_INLINE_FUNCTION T max() const;
  PORTABLE_INLINE_FUNCTION int rank() const { return rank_; }
  PORTABLE_INLINE_FUNCTION int size() const { return dataView_.GetSize(); }
  PORTABLE_INLINE_FUNCTION int sizeBytes() const {
    return dataView_.GetSizeInBytes();
  }
  PORTABLE_INLINE_FUNCTION int dim(int i) const { return dataView_.GetDim(i); }
  PORTABLE_INLINE_FUNCTION Grid_t range(int i) const;
  PORTABLE_INLINE_FUNCTION IndexType indexType(const int i) const {
    return indices_[i];
  }

  DataBox<T, Grid_t, Concept>
  getOnDevice() const { // getOnDevice is always a deep copy
    if (size() == 0 || status_ == DataStatus::Empty) { // edge case for unallocated
      DataBox<T, Grid_t, Concept> a;
      return a;
    }
    // create device memory (host memory if no device)
    T *device_data = (T *)PORTABLE_MALLOC(sizeBytes());
    // copy to device
    portableCopyToDevice(device_data, data_, sizeBytes());
    // create new databox of size size
    DataBox<T, Grid_t, Concept> a{device_data, dim(6), dim(5), dim(4),
                                  dim(3),      dim(2), dim(1)};
    a.copyShape(*this);
    // set correct allocation status of the new databox
    // note this is ALWAYS device, even if host==device.
    a.status_ = DataStatus::AllocatedDevice;
    return a;
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
  T *data_ = nullptr;           // points at data, managed or not
  PortableMDArray<T> dataView_; // always used
  IndexType indices_[MAXRANK];
  Grid_t grids_[MAXRANK];

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
    size_t size = sizeof(T) * prod_(std::forward<Args>(args)...);
    if (t == AllocationTarget::Device) {
      data_ = (T *)PORTABLE_MALLOC(size);
      status_ = DataStatus::AllocatedDevice;
    } else {
      data_ = (T *)malloc(size);
      status_ = DataStatus::AllocatedHost;
    }
  }
};

// Read an array, shallow
template <typename T, typename Grid_t, typename Concept>
inline void DataBox<T, Grid_t, Concept>::setArray(PortableMDArray<T> &A) {
  dataView_ = A;
  rank_ = A.GetRank();
  setAllIndexed_();
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION T
DataBox<T, Grid_t, Concept>::interpToReal(const T x) const noexcept {
  assert(canInterpToReal_(1));
  int ix;
  weights_t<T> w;
  grids_[0].weights(x, ix, w);
  return w[0] * dataView_(ix) + w[1] * dataView_(ix + 1);
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_FORCEINLINE_FUNCTION T DataBox<T, Grid_t, Concept>::interpToReal(
    const T x2, const T x1) const noexcept {
  assert(canInterpToReal_(2));
  int ix1, ix2;
  weights_t<T> w1, w2;
  grids_[0].weights(x1, ix1, w1);
  grids_[1].weights(x2, ix2, w2);
  // TODO: prefectch corners for speed?
  // TODO: re-order access pattern?
  return (w2[0] *
              (w1[0] * dataView_(ix2, ix1) + w1[1] * dataView_(ix2, ix1 + 1)) +
          w2[1] * (w1[0] * dataView_(ix2 + 1, ix1) +
                   w1[1] * dataView_(ix2 + 1, ix1 + 1)));
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_FORCEINLINE_FUNCTION T DataBox<T, Grid_t, Concept>::interpToReal(
    const T x3, const T x2, const T x1) const noexcept {
  assert(canInterpToReal_(3));
  int ix[3];
  weights_t<T> w[3];
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

template <typename T, typename Grid_t, typename Concept>
PORTABLE_FORCEINLINE_FUNCTION T DataBox<T, Grid_t, Concept>::interpToReal(
    const T x3, const T x2, const T x1, const int idx) const noexcept {
  assert(rank_ == 4);
  for (int r = 1; r < rank_; ++r) {
    assert(indices_[r] == IndexType::Interpolated);
    assert(grids_[r].isWellFormed());
  }
  int ix[3];
  weights_t<T> w[3];
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

// DH: this is a large function to force an inline, perhaps just make it a
// suggestion to the compiler?
template <typename T, typename Grid_t, typename Concept>
PORTABLE_FORCEINLINE_FUNCTION T DataBox<T, Grid_t, Concept>::interpToReal(
    const T x4, const T x3, const T x2, const T x1) const noexcept {
  assert(canInterpToReal_(4));
  T x[] = {x1, x2, x3, x4};
  int ix[4];
  weights_t<T> w[4];
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

template <typename T, typename Grid_t, typename Concept>
PORTABLE_FORCEINLINE_FUNCTION T DataBox<T, Grid_t, Concept>::interpToReal(
    const T x4, const T x3, const T x2, const int idx,
    const T x1) const noexcept {
  assert(rank_ == 5);
  assert(indices_[0] == IndexType::Interpolated);
  assert(grids_[0].isWellFormed());
  for (int i = 2; i < 5; ++i) {
    assert(indices_[i] == IndexType::Interpolated);
    assert(grids_[i].isWellFormed());
  }
  T x[] = {x1, x2, x3, x4};
  int ix[4];
  weights_t<T> w[4];
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

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION void
DataBox<T, Grid_t, Concept>::interpFromDB(const DataBox<T, Grid_t, Concept> &db,
                                          const T x) {
  assert(db.indices_[db.rank_ - 1] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 1].isWellFormed());
  assert(size() == (db.size() / db.dim(db.rank_)));

  int ix;
  weights_t<T> w;
  copyShape(db, 1);

  db.grids_[db.rank_ - 1].weights(x, ix, w);
  DataBox<T, Grid_t, Concept> lower(db.slice(ix)), upper(db.slice(ix + 1));
  // lower = db.slice(ix);
  // upper = db.slice(ix+1);
  for (int i = 0; i < size(); i++) {
    dataView_(i) = w[0] * lower(i) + w[1] * upper(i);
  }
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION void
DataBox<T, Grid_t, Concept>::interpFromDB(const DataBox<T, Grid_t, Concept> &db,
                                          const T x2, const T x1) {
  assert(db.rank_ >= 2);
  assert(db.indices_[db.rank_ - 1] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 1].isWellFormed());
  assert(db.indices_[db.rank_ - 2] == IndexType::Interpolated);
  assert(db.grids_[db.rank_ - 2].isWellFormed());
  assert(size() == (db.size() / (db.dim(db.rank_) * db.dim(db.rank_ - 1))));

  int ix2, ix1;
  weights_t<T> w2, w1;
  copyShape(db, 2);

  db.grids_[db.rank_ - 2].weights(x1, ix1, w1);
  db.grids_[db.rank_ - 1].weights(x2, ix2, w2);
  DataBox<T, Grid_t, Concept> corners[2][2]{
      {db.slice(ix2, ix1), db.slice(ix2 + 1, ix1)},
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
template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION void
DataBox<T, Grid_t, Concept>::copyShape(const DataBox<T, Grid_t, Concept> &db,
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
template <typename T, typename Grid_t, typename Concept>
inline void DataBox<T, Grid_t, Concept>::copyMetadata(
    const DataBox<T, Grid_t, Concept> &src) {
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
template <typename T, typename Grid_t, typename Concept>
inline herr_t
DataBox<T, Grid_t, Concept>::saveHDF(const std::string &filename) const {
  herr_t status;
  hid_t file;

  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status = saveHDF(file, SP5::DB::GRPNAME);
  status += H5Fclose(file);
  return status;
}

template <typename T, typename Grid_t, typename Concept>
inline herr_t
DataBox<T, Grid_t, Concept>::saveHDF(hid_t loc,
                                     const std::string &groupname) const {
  hid_t group, grids;
  herr_t status = 0;
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "Spiner HDF5 only defined for these data types: float, double");
  // Runtime because HDF5 is doing something evil under the hood with these
  // macros
  auto H5T_T =
      std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

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
                             dims_hsize_t.data(), H5T_T, dataView_.data());

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

template <typename T, typename Grid_t, typename Concept>
inline herr_t
DataBox<T, Grid_t, Concept>::loadHDF(const std::string &filename) {
  herr_t status;
  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  status = loadHDF(file, SP5::DB::GRPNAME);
  return status;
}

template <typename T, typename Grid_t, typename Concept>
inline herr_t
DataBox<T, Grid_t, Concept>::loadHDF(hid_t loc, const std::string &groupname) {
  hid_t group, grids;
  herr_t status = 0;
  std::vector<int> index_types;
  std::vector<int> dims(6, 1);
  static_assert(std::is_same<T, double>::value || std::is_same<T, float>::value,
                "Spiner HDF5 only defined for these data types: float, double");
  // Runtime because HDF5 is doing something evil under the hood with these
  // macros
  auto H5T_T =
      std::is_same<T, double>::value ? H5T_NATIVE_DOUBLE : H5T_NATIVE_FLOAT;

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
  status += H5LTread_dataset(group, SP5::DB::DSETNAME, H5T_T, dataView_.data());

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
template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION DataBox<T, Grid_t, Concept> &
DataBox<T, Grid_t, Concept>::operator=(const DataBox<T, Grid_t, Concept> &src) {
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
template <typename T, typename Grid_t, typename Concept>
inline void
DataBox<T, Grid_t, Concept>::copy(const DataBox<T, Grid_t, Concept> &src) {
  copyMetadata(src);
  for (int i = 0; i < src.size(); i++)
    dataView_(i) = src(i);
}

template <typename T, typename Grid_t, typename Concept>
inline bool DataBox<T, Grid_t, Concept>::operator==(
    const DataBox<T, Grid_t, Concept> &other) const {
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
template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION T DataBox<T, Grid_t, Concept>::min() const {
  T min = std::numeric_limits<T>::infinity();
  for (int i = 0; i < size(); i++) {
    min = std::min(min, dataView_(i));
  }
  return min;
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION T DataBox<T, Grid_t, Concept>::max() const {
  T max = -std::numeric_limits<T>::infinity();
  for (int i = 0; i < size(); i++) {
    max = std::max(max, dataView_(i));
  }
  return max;
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION Grid_t
DataBox<T, Grid_t, Concept>::range(int i) const {
  assert(0 <= i && i < rank_);
  assert(indices_[i] == IndexType::Interpolated);
  return grids_[i];
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION void DataBox<T, Grid_t, Concept>::setAllIndexed_() {
  for (int i = 0; i < rank_; i++) {
    indices_[i] = IndexType::Indexed;
  }
}

template <typename T, typename Grid_t, typename Concept>
PORTABLE_INLINE_FUNCTION bool
DataBox<T, Grid_t, Concept>::canInterpToReal_(const int interpOrder) const {
  if (rank_ != interpOrder) return false;
  for (int i = 0; i < rank_; i++) {
    if (indices_[i] != IndexType::Interpolated) return false;
    if (!(grids_[i].isWellFormed())) return false;
  }
  return true;
}

template <typename T, typename Grid_t, typename Concept>
inline DataBox<T, Grid_t, Concept>
getOnDeviceDataBox(const DataBox<T, Grid_t, Concept> &a_host) {
  return a_host.getOnDevice();
}
template <typename T, typename Grid_t, typename Concept>
inline void free(DataBox<T, Grid_t, Concept> &db) {
  db.finalize();
}

struct DBDeleter {
  template <typename T>
  void operator()(T *ptr) {
    ptr->finalize();
    delete ptr;
  }
};
} // namespace Spiner
#endif // _SPINER_DATABOX_HPP_
