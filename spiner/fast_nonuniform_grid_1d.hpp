#ifndef SPINER_FAST_NONUNIFORM_GRID_1D_
#define SPINER_FAST_NONUNIFORM_GRID_1D_
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

#include <algorithm>
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

#include "nonuniform_grid_1d.hpp"
#include "ports-of-call/nqt_math.hpp"
#include "ports-of-call/portability.hpp"
#include "ports-of-call/portable_errors.hpp"
#include "ports-of-call/robust_utils.hpp"
#include "regular_grid_1d.hpp"
#include "spiner_types.hpp"

namespace Spiner {

template <typename T = Real>
class FastNonUniformGrid1D {
  static_assert(std::is_same<T, double>::value,
                "FastNonUniformGrid1D is currently defined only for double");

 public:
  using ValueType = T;
  static constexpr int DEFAULT_MAX_LOOKUP_RATIO = 32;

  enum class Policy { Automatic = 0, RequireFast = 1, ForceBinary = 2 };

  struct Settings {
    // A negative scale infers the transition scale from the coordinates.
    T scale = T(-1);
    Policy policy = Policy::Automatic;
    int max_lookup_ratio = DEFAULT_MAX_LOOKUP_RATIO;
  };

  FastNonUniformGrid1D() = default;

  FastNonUniformGrid1D(const std::vector<T> &points,
                       const Settings &settings = {})
      : coordinates_(points), requested_policy_(settings.policy),
        max_lookup_ratio_(settings.max_lookup_ratio) {
    scale_ = resolveScale_(settings.scale);
    validateSettings_();
    validateScaleRange_();
    configureLookup_();
  }

  FastNonUniformGrid1D(std::initializer_list<T> points,
                       const Settings &settings = {})
      : FastNonUniformGrid1D(std::vector<T>(points), settings) {}

  PORTABLE_INLINE_FUNCTION T x(const int i) const { return coordinates_.x(i); }

  PORTABLE_INLINE_FUNCTION int index(const T value) const {
    if (!usesFastLookup()) return coordinates_.index(value);
    if (value <= min()) return 0;
    if (value >= max()) return nPoints() - 2;

    const T transformed = transform_(value / scale_);
    int lookup_index =
        std::max(0, std::min(lookup_grid_.index(transformed),
                             lookupSize() - 1));
    int coordinate_index = lookup_[lookup_index];
    // the lookup table is essentially discrete interpolation. This
    // off-by-one check interpolates to the correct index.
    if (coordinate_index < nPoints() - 2 &&
        value >= coordinates_.x(coordinate_index + 1)) {
      ++coordinate_index;
    }
    return coordinate_index;
  }

  PORTABLE_INLINE_FUNCTION void weights(const T &value, int &ix,
                                        weights_t<T> &w) const {
    ix = index(value);
    const T dx = coordinates_.x(ix + 1) - coordinates_.x(ix);
    w[1] = (value - coordinates_.x(ix)) / dx;
    w[0] = T(1) - w[1];
  }

  PORTABLE_INLINE_FUNCTION T min() const { return coordinates_.min(); }
  PORTABLE_INLINE_FUNCTION T max() const { return coordinates_.max(); }
  PORTABLE_INLINE_FUNCTION std::size_t nPoints() const {
    return coordinates_.nPoints();
  }
  PORTABLE_INLINE_FUNCTION bool isWellFormed() const {
    const bool lookup_well_formed =
        !usesFastLookup() || lookup_grid_.isWellFormed();
    return coordinates_.isWellFormed() && std::isfinite(scale_) && scale_ > 0 &&
           max_lookup_ratio_ > 0 && lookup_well_formed;
  }
  PORTABLE_INLINE_FUNCTION DataStatus dataStatus() const {
    return coordinates_.dataStatus();
  }
  PORTABLE_INLINE_FUNCTION const T *data() const { return coordinates_.data(); }
  PORTABLE_INLINE_FUNCTION T scale() const { return scale_; }
  Settings settings() const {
    return {scale_, requested_policy_, max_lookup_ratio_};
  }
  PORTABLE_INLINE_FUNCTION int lookupSize() const {
    return usesFastLookup() ? lookup_grid_.nPoints() - 1 : 0;
  }
  PORTABLE_INLINE_FUNCTION int maxLookupRatio() const {
    return max_lookup_ratio_;
  }
  PORTABLE_INLINE_FUNCTION Policy requestedPolicy() const {
    return requested_policy_;
  }
  PORTABLE_INLINE_FUNCTION bool usesFastLookup() const {
    return lookup_status_ != DataStatus::Empty;
  }

  void reconfigureLookup(const Settings &settings) {
    PORTABLE_ALWAYS_REQUIRE(
        coordinates_.dataStatus() == DataStatus::AllocatedHost,
        "Lookup reconfiguration requires a host-owned grid");
    const T scale = resolveScale_(settings.scale);
    validateSettings_(scale, settings.policy, settings.max_lookup_ratio);
    validateScaleRange_(scale);
    const bool scale_changed = scale != scale_;
    scale_ = scale;
    requested_policy_ = settings.policy;
    max_lookup_ratio_ = settings.max_lookup_ratio;

    if (requested_policy_ == Policy::ForceBinary) {
      releaseLookup_();
      return;
    }

    const int max_entries = maxLookupEntries_();
    if (!scale_changed && usesFastLookup() && lookupSize() <= max_entries) {
      return;
    }
    releaseLookup_();
    const bool success = buildLookup_(max_entries);
    if (requested_policy_ == Policy::RequireFast) {
      PORTABLE_ALWAYS_REQUIRE(
          success, "Fast lookup table exceeds its limit or is invalid");
    }
  }

  void reconfigureLookup(const Policy policy, const int max_lookup_ratio) {
    Settings updated = settings();
    updated.policy = policy;
    updated.max_lookup_ratio = max_lookup_ratio;
    reconfigureLookup(updated);
  }

  std::size_t dynamicMemorySizeInBytes() const {
    return coordinates_.dynamicMemorySizeInBytes() + lookupSize() * sizeof(int);
  }
  std::size_t serializedSizeInBytes() const {
    return sizeof(*this) + dynamicMemorySizeInBytes();
  }
  std::size_t dumpDynamicMemory(std::byte *dst) const {
    PORTABLE_REQUIRE(dataStatus() != DataStatus::AllocatedDevice &&
                         lookup_status_ != DataStatus::AllocatedDevice,
                     "Cannot dump device-resident fast grid memory");
    std::size_t offset = coordinates_.dumpDynamicMemory(dst);
    if (usesFastLookup()) {
      std::memcpy(dst + offset, lookup_, lookupSize() * sizeof(int));
      offset += lookupSize() * sizeof(int);
    }
    return offset;
  }
  std::size_t serialize(std::byte *dst) const {
    PORTABLE_REQUIRE(dataStatus() != DataStatus::AllocatedDevice,
                     "Cannot serialize device-resident fast grid memory");
    std::memcpy(dst, this, sizeof(*this));
    return sizeof(*this) + dumpDynamicMemory(dst + sizeof(*this));
  }
  std::size_t setPointer(std::byte *src) {
    std::size_t offset = coordinates_.setPointer(src);
    if (usesFastLookup()) {
      const std::size_t lookup_size = lookupSize();
      lookup_ = reinterpret_cast<int *>(src + offset);
      lookup_status_ = DataStatus::Unmanaged;
      offset += lookup_size * sizeof(int);
    } else {
      lookup_ = nullptr;
      lookup_status_ = DataStatus::Empty;
    }
    return offset;
  }
  std::size_t deSerialize(std::byte *src) {
    PORTABLE_REQUIRE((dataStatus() == DataStatus::Empty ||
                      dataStatus() == DataStatus::Unmanaged) &&
                         (lookup_status_ == DataStatus::Empty ||
                          lookup_status_ == DataStatus::Unmanaged),
                     "Must not de-serialize into an active fast grid");
    std::memcpy(this, src, sizeof(*this));
    validateSettings_();
    const std::size_t offset = sizeof(*this) + setPointer(src + sizeof(*this));
    validateScaleRange_();
    return offset;
  }

  FastNonUniformGrid1D getOnDevice() const {
    PORTABLE_REQUIRE(dataStatus() != DataStatus::AllocatedDevice,
                     "Cannot copy a device fast grid to device");
    FastNonUniformGrid1D grid;
    grid.coordinates_ = coordinates_.getOnDevice();
    grid.lookup_grid_ = lookup_grid_;
    grid.scale_ = scale_;
    grid.requested_policy_ = requested_policy_;
    grid.max_lookup_ratio_ = max_lookup_ratio_;
    if (usesFastLookup()) {
      grid.lookup_ =
          static_cast<int *>(PORTABLE_MALLOC(lookupSize() * sizeof(int)));
      PORTABLE_ALWAYS_REQUIRE(grid.lookup_ != nullptr,
                              "Lookup allocation failed");
      portableCopyToDevice(grid.lookup_, lookup_, lookupSize() * sizeof(int));
      grid.lookup_status_ = DataStatus::AllocatedDevice;
    }
    return grid;
  }

  void copy(const FastNonUniformGrid1D &other) {
    if (this == &other) return;
    PORTABLE_REQUIRE((dataStatus() == DataStatus::Empty ||
                      dataStatus() == DataStatus::Unmanaged) &&
                         (lookup_status_ == DataStatus::Empty ||
                          lookup_status_ == DataStatus::Unmanaged),
                     "Must not copy into an active fast grid");
    PORTABLE_REQUIRE(other.dataStatus() != DataStatus::AllocatedDevice &&
                         other.lookup_status_ != DataStatus::AllocatedDevice,
                     "Cannot deep copy a device-resident fast grid to host");

    coordinates_.copy(other.coordinates_);
    lookup_grid_ = other.lookup_grid_;
    scale_ = other.scale_;
    requested_policy_ = other.requested_policy_;
    max_lookup_ratio_ = other.max_lookup_ratio_;
    if (other.usesFastLookup()) {
      lookup_ =
          static_cast<int *>(std::malloc(other.lookupSize() * sizeof(int)));
      PORTABLE_ALWAYS_REQUIRE(lookup_ != nullptr, "Lookup allocation failed");
      std::memcpy(lookup_, other.lookup_, other.lookupSize() * sizeof(int));
      lookup_status_ = DataStatus::AllocatedHost;
    } else {
      lookup_ = nullptr;
      lookup_status_ = DataStatus::Empty;
    }
  }

  void finalize() {
    coordinates_.finalize();
    if (lookup_status_ != DataStatus::Unmanaged) releaseLookup_();
  }

#ifdef SPINER_USE_HDF
  inline herr_t saveHDF(hid_t loc, const std::string &name) const {
    PORTABLE_REQUIRE(dataStatus() != DataStatus::AllocatedDevice,
                     "Cannot save device-resident fast grid memory");
    herr_t status = 0;
    hid_t group =
        H5Gcreate(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status += coordinates_.saveHDF(group, SP5::FNG1D::COORDINATES);
    const int policy = static_cast<int>(requested_policy_);
    PORTABLE_REQUIRE(max_lookup_ratio_ <= std::numeric_limits<int>::max(),
                     "Maximum lookup ratio cannot be represented in HDF5");
    status += H5LTset_attribute_double(loc, name.c_str(), SP5::FNG1D::SCALE,
                                       &scale_, 1);
    status += H5LTset_attribute_int(loc, name.c_str(),
                                    SP5::FNG1D::LOOKUP_POLICY, &policy, 1);
    status += H5LTset_attribute_int(
        loc, name.c_str(), SP5::FNG1D::MAX_LOOKUP_RATIO, &max_lookup_ratio_, 1);
    status += H5Gclose(group);
    return status;
  }

  inline herr_t loadHDF(hid_t loc, const std::string &name) {
    finalize();
    if (lookup_status_ == DataStatus::Unmanaged) clearLookupMetadata_();

    herr_t status = 0;
    hid_t group = H5Gopen(loc, name.c_str(), H5P_DEFAULT);
    status += coordinates_.loadHDF(group, SP5::FNG1D::COORDINATES);
    int policy = 0;
    int ratio = 0;
    status +=
        H5LTget_attribute_double(loc, name.c_str(), SP5::FNG1D::SCALE, &scale_);
    status += H5LTget_attribute_int(loc, name.c_str(),
                                    SP5::FNG1D::LOOKUP_POLICY, &policy);
    status += H5LTget_attribute_int(loc, name.c_str(),
                                    SP5::FNG1D::MAX_LOOKUP_RATIO, &ratio);
    status += H5Gclose(group);
    PORTABLE_ALWAYS_REQUIRE(ratio > 0, "Maximum lookup ratio must be positive");
    requested_policy_ = static_cast<Policy>(policy);
    max_lookup_ratio_ = ratio;
    validateSettings_();
    validateScaleRange_();
    configureLookup_();
    return status;
  }
#endif

 private:
  PORTABLE_INLINE_FUNCTION static T transform_(const T value) {
#ifdef SPINER_USE_PORTABLE_NQT
    return PortsOfCall::NQT::O1::Portable::asinh(value);
#else
    return PortsOfCall::NQT::O1::Aliased::asinh(value);
#endif
  }

  // necessary because the enum might be set by static cast
  static bool validPolicy_(const Policy policy) {
    return policy == Policy::Automatic || policy == Policy::RequireFast ||
           policy == Policy::ForceBinary;
  }

  T resolveScale_(const T requested_scale) const {
    PORTABLE_ALWAYS_REQUIRE(std::isfinite(requested_scale) &&
                                requested_scale != T(0),
                            "Fast grid scale must be finite and nonzero");
    if (requested_scale > T(0)) return requested_scale;

    T inferred_scale = std::numeric_limits<T>::infinity();
    for (int i = 0; i < nPoints(); ++i) {
      const T magnitude = std::abs(coordinates_.x(i));
      if (magnitude >= PortsOfCall::Robust::SMALL<T>()) {
        inferred_scale = std::min(inferred_scale, magnitude);
      }
    }
    PORTABLE_ALWAYS_REQUIRE(
        std::isfinite(inferred_scale),
        "Cannot infer a fast grid scale from these coordinates");
    return inferred_scale;
  }

  void validateSettings_() const {
    validateSettings_(scale_, requested_policy_, max_lookup_ratio_);
  }

  static void validateSettings_(const T scale, const Policy policy,
                                const int max_lookup_ratio) {
    PORTABLE_ALWAYS_REQUIRE(std::isfinite(scale) && scale > 0,
                            "Fast grid scale must be finite and positive");
    PORTABLE_ALWAYS_REQUIRE(max_lookup_ratio > 0,
                            "Maximum lookup ratio must be positive");
    PORTABLE_ALWAYS_REQUIRE(validPolicy_(policy), "Invalid fast lookup policy");
  }

  void validateScaleRange_(const T scale) const {
    for (int i = 0; i < nPoints(); ++i) {
      PORTABLE_ALWAYS_REQUIRE(
          std::isfinite(coordinates_.x(i) / scale),
          "Fast grid coordinate range is too large for its scale");
    }
  }

  void validateScaleRange_() const { validateScaleRange_(scale_); }

  std::size_t maxLookupEntries_() const {
    // guard against overflow... probably not necessary?
    if (nPoints() > std::numeric_limits<int>::max() / max_lookup_ratio_) {
      return std::numeric_limits<int>::max();
    }
    return nPoints() * max_lookup_ratio_;
  }

  void configureLookup_() {
    if (requested_policy_ == Policy::ForceBinary) return;
    const bool success = buildLookup_(maxLookupEntries_());
    if (requested_policy_ == Policy::RequireFast) {
      PORTABLE_ALWAYS_REQUIRE(
          success, "Fast lookup table exceeds its limit or is invalid");
    }
  }

  bool buildLookup_(const int max_entries) {
    const T transformed_min = transform_(min() / scale_);
    const T transformed_max = transform_(max() / scale_);
    if (!(std::isfinite(transformed_min) && std::isfinite(transformed_max) &&
          (transformed_min < transformed_max))) {
      return false;
    }

    // compute min spacing between coordinates
    T min_spacing = std::numeric_limits<T>::infinity();
    T previous = transformed_min;
    for (int i = 1; i < nPoints(); ++i) {
      const T current = transform_(coordinates_.x(i) / scale_);
      const T spacing = current - previous;
      if (!std::isfinite(current) || !(spacing > 0)) return false;
      min_spacing = std::min(min_spacing, spacing);
      previous = current;
    }

    // check our lookup table fits. required is nonfinite if we risk overflow
    const T slightly_smaller_spacing = std::nextafter(min_spacing, T(0));
    const T span = transformed_max - transformed_min;
    const T required = std::ceil(span / slightly_smaller_spacing);
    const int max_cells = std::min(
        max_entries, std::numeric_limits<int>::max() - 1);
    if (!std::isfinite(required) || required < 1 ||
        required > static_cast<T>(max_cells)) {
      return false;
    }

    // One additional cell provides margin against rounding in the lookup grid.
    const int cells = required;
    if (cells == max_cells) return false;
    const int conservative_cells = cells + 1;
    const RegularGrid1D<T> candidate_grid(transformed_min, transformed_max,
                                          conservative_cells + 1);
    // Sanity check to ensure lookup table doesn't predict wrong point
    // in coordinates due to roundoff
    if (!lookupLayoutIsValid_(candidate_grid, conservative_cells)) {
      return false;
    }

    int *candidate =
        static_cast<int *>(std::malloc(conservative_cells * sizeof(int)));
    if (candidate == nullptr) return false;
    fillLookup_(candidate, candidate_grid, conservative_cells);
    lookup_ = candidate;
    lookup_grid_ = candidate_grid;
    lookup_status_ = DataStatus::AllocatedHost;
    return true;
  }

  bool lookupLayoutIsValid_(const RegularGrid1D<T> &grid,
                            const int cells) const {
    int source_index = 0;
    // for each lookup index, find source index
    for (int j = 0; j < cells; ++j) {
      const T left = grid.x(j);
      const T right = grid.x(j + 1);
      while (source_index < nPoints() - 2 &&
             transform_(coordinates_.x(source_index + 1) / scale_) < left) {
        ++source_index;
      }
      // and ensure there's only 1 source index per lookup index
      if (source_index + 2 <= nPoints() - 2 &&
          transform_(coordinates_.x(source_index + 2) / scale_) < right) {
        return false;
      }
    }
    return true;
  }

  void fillLookup_(int *lookup, const RegularGrid1D<T> &grid,
                   const std::size_t cells) const {
    int source_index = 0;
    for (std::size_t j = 0; j < cells; ++j) {
      const T left = grid.x(j);
      while (source_index < nPoints() - 2 &&
             transform_(coordinates_.x(source_index + 1) / scale_) < left) {
        ++source_index;
      }
      lookup[j] = source_index;
    }
  }

  void clearLookupMetadata_() {
    lookup_ = nullptr;
    lookup_grid_ = RegularGrid1D<T>{};
    lookup_status_ = DataStatus::Empty;
  }

  void releaseLookup_() {
    if (lookup_status_ == DataStatus::AllocatedHost) {
      std::free(lookup_);
    } else if (lookup_status_ == DataStatus::AllocatedDevice) {
      PORTABLE_FREE(lookup_);
    }
    clearLookupMetadata_();
  }

  NonUniformGrid1D<T> coordinates_;
  RegularGrid1D<T> lookup_grid_;
  int *lookup_ = nullptr;
  T scale_ = std::numeric_limits<T>::signaling_NaN();
  Policy requested_policy_ = Policy::Automatic;
  int max_lookup_ratio_ = 0;
  DataStatus lookup_status_ = DataStatus::Empty;
};

} // namespace Spiner

#endif // SPINER_FAST_NONUNIFORM_GRID_1D_
