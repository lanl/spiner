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

// Generative AI was used to assist with modifications to this file.

#include <algorithm> // std::min, std::max
#include <array>
#include <cmath> // sqrt
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>
#include <spiner/interp_toy.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/spiner_types.hpp>

#ifdef SPINER_USE_HDF
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

using DataBox = Spiner::DataBox<Real>;
using Spiner::IndexType;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using NonUniformGrid1D = Spiner::NonUniformGrid1D<Real>;

// Because this uses NQT only valid specifically for doubles
using FastNonUniformGrid1D = Spiner::FastNonUniformGrid1D<double>;
using FastGridSettings = FastNonUniformGrid1D::Settings;
using FastGridPolicy = FastNonUniformGrid1D::Policy;

using Spiner::DBDeleter;
const Real EPSTEST = std::sqrt(DataBox::EPS);
template <int N>
using PiecewiseGrid1D = Spiner::PiecewiseGrid1D<Real, N>;
template <int N>
using PiecewiseDB = Spiner::DataBox<Real, PiecewiseGrid1D<N>>;
using Spiner::DataStatus;
using NonUniformDB = Spiner::DataBox<Real, NonUniformGrid1D>;
using FastNonUniformDB = Spiner::DataBox<double, FastNonUniformGrid1D>;

PORTABLE_INLINE_FUNCTION double nqtSinhForFastGridTest(const double x) {
#ifdef SPINER_USE_PORTABLE_NQT
  return PortsOfCall::NQT::O1::Portable::sinh(x);
#else
  return PortsOfCall::NQT::O1::Aliased::sinh(x);
#endif
}

PORTABLE_INLINE_FUNCTION Real linearFunction(Real z, Real y, Real x) {
  return x + y + z;
}
PORTABLE_INLINE_FUNCTION Real linearFunction(Real a, Real z, Real y, Real x) {
  return x + y + z + a;
}
PORTABLE_INLINE_FUNCTION Real linearFunction(Real b, Real a, Real z, Real y,
                                             Real x) {
  return x + y + z + a + b;
}

SCENARIO("PortableMDArrays can be allocated from a pointer",
         "[PortableMDArray]") {
  constexpr int N = 2;
  constexpr int M = 3;
  std::vector<int> data(N * M);
  PortableMDArray<int> a;
  int tot = 0;
  for (int i = 0; i < N * M; i++) {
    data[i] = tot;
    tot++;
  }
  a.NewPortableMDArray(data.data(), M, N);

  SECTION("Shape should be NxM") {
    REQUIRE(a.GetDim1() == N);
    REQUIRE(a.GetDim2() == M);
  }

  SECTION("Stride is as set by initialized pointer") {
    int tot = 0;
    for (int j = 0; j < M; j++) {
      for (int i = 0; i < N; i++) {
        REQUIRE(a(j, i) == tot);
        tot++;
      }
    }
  }

  SECTION("Identical slices of the same data should compare equal") {
    PortableMDArray<int> aslc1, aslc2;
    aslc1.InitWithShallowSlice(a, 1, 0, 2);
    aslc2.InitWithShallowSlice(a, 1, 0, 2);
    REQUIRE(aslc1 == aslc2);
  }
}

TEST_CASE("RegularGrid1D", "[RegularGrid1D]") {
  SECTION("A regular grid 1d emits appropriate metadata") {
    constexpr Real min = -1;
    constexpr Real max = 1;
    constexpr size_t N = 10;
    RegularGrid1D g(min, max, N);
    REQUIRE(g.min() == min);
    REQUIRE(g.max() == max);
    REQUIRE(g.nPoints() == N);
  }

  SECTION("A regular grid can be serialized and deserialized") {
    RegularGrid1D grid(-1.0, 2.0, 7);
    std::vector<std::byte> serialized(grid.serializedSizeInBytes());

    const std::size_t written = grid.serialize(serialized.data());
    REQUIRE(written == serialized.size());

    RegularGrid1D restored;
    const std::size_t consumed = restored.deSerialize(serialized.data());
    REQUIRE(consumed == written);
    REQUIRE(restored.min() == grid.min());
    REQUIRE(restored.max() == grid.max());
    REQUIRE(restored.nPoints() == grid.nPoints());
    REQUIRE(restored.setPointer(serialized.data()) == 0);
    const auto device_grid = restored.getOnDevice();
    REQUIRE(device_grid.min() == grid.min());
    REQUIRE(device_grid.max() == grid.max());
    REQUIRE(device_grid.nPoints() == grid.nPoints());
    restored.finalize();
    grid.finalize();
  }

  SECTION("Shallow and deep copies of a regular grid are plain copies") {
    RegularGrid1D grid(-1.0, 2.0, 7);
    RegularGrid1D shallow, deep;
    shallow.shallowCopy(grid);
    deep.deepCopy(grid);
    for (const RegularGrid1D &g : {shallow, deep}) {
      REQUIRE(g.dataStatus() == DataStatus::Trivial);
      REQUIRE(g.min() == grid.min());
      REQUIRE(g.max() == grid.max());
      REQUIRE(g.nPoints() == grid.nPoints());
    }
  }
}

TEST_CASE("NonUniformGrid1D", "[NonUniformGrid1D]") {
  const std::vector<Real> points = {-1.0, -0.5, 0.25, 2.0};

  SECTION("An owning grid maps coordinates and computes local weights") {
    NonUniformGrid1D grid(points);
    REQUIRE(grid.dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(grid.nPoints() == points.size());
    REQUIRE(grid.min() == points.front());
    REQUIRE(grid.max() == points.back());
    REQUIRE(grid.x(2) == points[2]);
    REQUIRE(grid.index(-2.0) == 0);
    REQUIRE(grid.index(0.0) == 1);
    REQUIRE(grid.index(2.0) == 2);

    int ix;
    Spiner::weights_t<Real> w;
    grid.weights(0.0, ix, w);
    REQUIRE(ix == 1);
    REQUIRE(std::abs(w[0] - 1.0 / 3.0) <= EPSTEST);
    REQUIRE(std::abs(w[1] - 2.0 / 3.0) <= EPSTEST);

    grid.weights(2.0, ix, w);
    REQUIRE(ix == 2);
    REQUIRE(std::abs(w[0]) <= EPSTEST);
    REQUIRE(std::abs(w[1] - 1.0) <= EPSTEST);
    grid.finalize();
    grid.finalize();
    REQUIRE(grid.dataStatus() == DataStatus::Empty);
  }

  SECTION("A borrowed grid never frees its caller-owned coordinates") {
    std::vector<Real> borrowed_points = points;
    NonUniformGrid1D grid(borrowed_points.data(), borrowed_points.size());
    REQUIRE(grid.dataStatus() == DataStatus::Unmanaged);
    grid.finalize();
    REQUIRE(std::abs(borrowed_points[2] - 0.25) <= EPSTEST);
  }

  SECTION("deepCopy creates independent host-owned coordinates") {
    NonUniformGrid1D source(points);
    NonUniformGrid1D copy;
    copy.deepCopy(source);
    REQUIRE(copy.dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(copy.data() != source.data());
    REQUIRE(copy.nPoints() == source.nPoints());
    for (std::size_t i = 0; i < source.nPoints(); ++i)
      REQUIRE(copy.x(i) == source.x(i));

    source.data()[1] = -0.25;
    REQUIRE(copy.x(1) == -0.5);
    copy.finalize();
    source.finalize();
  }

  SECTION("A non uniform grid can be deep copied.") {
    NonUniformGrid1D source(points);
    NonUniformGrid1D copy;
    copy.deepCopy(source);
    REQUIRE(copy.data() != source.data());
    REQUIRE(source.x(1) == -0.5);
    copy.finalize();
    source.finalize();
  }

  SECTION("shallowCopy creates a non-owning handle") {
    NonUniformGrid1D source(points);
    NonUniformGrid1D shallow;
    shallow.shallowCopy(source);
    REQUIRE(shallow.dataStatus() == DataStatus::Unmanaged);
    REQUIRE(source.dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(shallow.data() == source.data());
    REQUIRE(shallow.nPoints() == source.nPoints());

    // An ordinary copy, by contrast, copies the ownership status.
    NonUniformGrid1D alias = source;
    REQUIRE(alias.dataStatus() == DataStatus::AllocatedHost);

    // Finalizing the non-owning handle leaves the owner intact.
    shallow.finalize();
    REQUIRE(source.dataStatus() == DataStatus::AllocatedHost);
    for (std::size_t i = 0; i < points.size(); ++i)
      REQUIRE(source.x(i) == points[i]);
    source.finalize();
  }

  SECTION("shallowCopy of an empty grid is empty") {
    NonUniformGrid1D source;
    NonUniformGrid1D shallow;
    shallow.shallowCopy(source);
    REQUIRE(shallow.dataStatus() == DataStatus::Empty);
  }

  SECTION("getOnDevice creates device-owned coordinates") {
    NonUniformGrid1D host_grid(points);
    NonUniformGrid1D grid = host_grid.getOnDevice();
    REQUIRE(grid.dataStatus() == DataStatus::AllocatedDevice);

    Real value = 0;
    portableReduce(
        "Interpolate with a directly allocated nonuniform grid", 0, 1,
        PORTABLE_LAMBDA(const int, Real &result) {
          int ix;
          Spiner::weights_t<Real> w;
          grid.weights(0.0, ix, w);
          result += w[0] * grid.x(ix) + w[1] * grid.x(ix + 1);
        },
        value);
    REQUIRE(std::abs(value) <= EPSTEST);
    grid.finalize();
    host_grid.finalize();
  }

  SECTION("A grid serializes its owned coordinates and relocates them") {
    NonUniformGrid1D grid(points);
    std::vector<std::byte> serialized(grid.serializedSizeInBytes());
    REQUIRE(grid.serialize(serialized.data()) == serialized.size());

    NonUniformGrid1D restored;
    REQUIRE(restored.deSerialize(serialized.data()) == serialized.size());
    REQUIRE(restored.nPoints() == grid.nPoints());
    for (std::size_t i = 0; i < grid.nPoints(); ++i)
      REQUIRE(restored.x(i) == grid.x(i));
    REQUIRE(restored.dataStatus() == DataStatus::Unmanaged);
    REQUIRE(reinterpret_cast<const std::byte *>(restored.data()) >=
            serialized.data());
    restored.finalize();
    grid.finalize();
  }

  SECTION(
      "A DataBox delegates serialization and device copying to every axis") {
    constexpr int N = 4;
    constexpr int RANK = 3;
    NonUniformDB db(N, N, N);
    db.setRange(0, points);
    db.setRange(1, std::vector<Real>{-2.0, -0.25, 0.5, 3.0});
    db.setRange(2, std::vector<Real>{-3.0, -1.0, 0.5, 4.0});
    for (int k = 0; k < N; ++k) {
      for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
          db(k, j, i) = db.range(0).x(i) + 2.0 * db.range(1).x(j) +
                        3.0 * db.range(2).x(k);
        }
      }
    }

    NonUniformDB device = db.getOnDevice();
    for (int i = 0; i < RANK; ++i) {
      REQUIRE(device.range(i).data() != db.range(i).data());
      REQUIRE(device.range(i).dataStatus() == DataStatus::AllocatedDevice);
    }
    Real value = 0;
    portableReduce(
        "Interpolate with nonuniform grids", 0, 1,
        PORTABLE_LAMBDA(const int, Real &result) {
          result += device.interpToReal(Real(0), Real(0), Real(0));
        },
        value);
    REQUIRE(std::abs(value) <= EPSTEST);
    // finalize frees only values. Grids are freed by finalizeGrids.
    device.finalize();
    for (int i = 0; i < RANK; ++i)
      REQUIRE(device.range(i).dataStatus() == DataStatus::AllocatedDevice);
    device.finalizeGrids();
    for (int i = 0; i < RANK; ++i)
      REQUIRE(device.range(i).dataStatus() == DataStatus::Empty);

    db.setIndexType(1, IndexType::Indexed);
    db.setIndexType(2, IndexType::Named);
    std::size_t expected = sizeof(db) + db.sizeBytes();
    for (int i = 0; i < RANK; ++i)
      expected += db.range(i).dynamicMemorySizeInBytes();
    REQUIRE(db.serializedSizeInBytes() == expected);

    std::vector<std::byte> serialized(db.serializedSizeInBytes());
    REQUIRE(db.serialize(serialized.data()) == serialized.size());
    NonUniformDB restored;
    REQUIRE(restored.deSerialize(serialized.data()) == serialized.size());
    REQUIRE(restored.indexType(0) == IndexType::Interpolated);
    REQUIRE(restored.indexType(1) == IndexType::Indexed);
    REQUIRE(restored.indexType(2) == IndexType::Named);
    const std::byte *begin = serialized.data();
    const std::byte *end = begin + serialized.size();
    for (int i = 0; i < RANK; ++i) {
      const std::byte *grid_data =
          reinterpret_cast<const std::byte *>(restored.range(i).data());
      REQUIRE(restored.range(i).dataStatus() == DataStatus::Unmanaged);
      REQUIRE(grid_data >= begin);
      REQUIRE(grid_data < end);
      REQUIRE(restored.range(i).nPoints() == db.range(i).nPoints());
      for (std::size_t j = 0; j < db.range(i).nPoints(); ++j)
        REQUIRE(restored.range(i).x(j) == db.range(i).x(j));
    }

    db.finalizeGrids();
    db.finalize();
    for (int i = 0; i < RANK; ++i)
      REQUIRE(db.range(i).dataStatus() == DataStatus::Empty);
  }

  SECTION("A DataBox can copy values without deep-copying grids") {
    constexpr int N = 4;
    NonUniformDB db(N);
    db.setRange(0, points);
    for (int i = 0; i < N; ++i)
      db(i) = static_cast<Real>(i);

    NonUniformDB device = db.getOnDevice(false);
    REQUIRE(device.dataStatus() == DataStatus::AllocatedDevice);
    REQUIRE(&device.range(0) != &db.range(0));
    REQUIRE(device.range(0).data() == db.range(0).data());
    REQUIRE(device.range(0).dataStatus() == DataStatus::Unmanaged);

    Real sum = 0;
    portableReduce(
        "Read a DataBox copied without grids", 0, N,
        PORTABLE_LAMBDA(const int i, Real &result) { result += device(i); },
        sum);
    REQUIRE(sum == 6.0);

    // The shallow grid alias does not finalize the host-owned coordinates.
    device.finalizeGrids();
    REQUIRE(db.range(0).dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(db.range(0).x(1) == points[1]);
    device.finalize();
    db.finalizeGrids();
    db.finalize();
  }

  SECTION("A DataBox separates ownership of values and grids") {
    constexpr int N = 4;
    NonUniformDB db(N, N);
    db.setRange(0, points);
    db.setRange(1, points);

    SECTION("Slices hold non-owning grids") {
      NonUniformDB slice = db.slice(1);
      REQUIRE(slice.dataStatus() == DataStatus::Unmanaged);
      REQUIRE(slice.range(0).dataStatus() == DataStatus::Unmanaged);
      REQUIRE(slice.range(0).data() == db.range(0).data());
      slice.finalizeGrids();
      REQUIRE(db.range(0).dataStatus() == DataStatus::AllocatedHost);
    }

    SECTION("Copies and assignment are handles with the same ownership") {
      NonUniformDB copy(db);
      NonUniformDB assigned;
      assigned = db;
      for (const NonUniformDB *handle : {&copy, &assigned}) {
        REQUIRE(handle->data() == db.data());
        for (int i = 0; i < db.rank(); ++i) {
          REQUIRE(handle->range(i).data() == db.range(i).data());
          REQUIRE(handle->range(i).dataStatus() == DataStatus::AllocatedHost);
        }
      }
    }

    SECTION("deepCopy copies values and grids") {
      for (int i = 0; i < db.size(); ++i)
        db(i) = static_cast<Real>(i);
      NonUniformDB copy;
      copy.deepCopy(db);
      REQUIRE(copy.data() != db.data());
      for (int i = 0; i < db.size(); ++i)
        REQUIRE(copy(i) == db(i));
      for (int i = 0; i < db.rank(); ++i) {
        REQUIRE(copy.indexType(i) == db.indexType(i));
        REQUIRE(copy.range(i).dataStatus() == DataStatus::AllocatedHost);
        REQUIRE(copy.range(i).data() != db.range(i).data());
      }
      // The copy is fully independent of the original.
      copy.finalizeGrids();
      copy.finalize();
      REQUIRE(db.range(0).x(1) == points[1]);
    }

    SECTION("finalize leaves grids and finalizeGrids frees them") {
      NonUniformDB other(N);
      other.setRange(0, points);
      other.finalize();
      REQUIRE(other.dataStatus() == DataStatus::Empty);
      REQUIRE(other.range(0).dataStatus() == DataStatus::AllocatedHost);
      other.finalizeGrids();
      REQUIRE(other.range(0).dataStatus() == DataStatus::Empty);
    }

    db.finalizeGrids();
    db.finalize();
  }
}

TEST_CASE("FastNonUniformGrid1D", "[FastNonUniformGrid1D]") {
  using FastReferenceGrid = Spiner::NonUniformGrid1D<double>;
  const std::vector<double> points = {-100.0, -20.0, -4.0, -1.0, 0.0,
                                      1.0,    4.0,   20.0, 100.0};

  SECTION("Fast lookup exactly matches binary lookup") {
    FastReferenceGrid binary(points);
    FastNonUniformGrid1D fast(points);
    REQUIRE(fast.usesFastLookup());
    REQUIRE(fast.maxLookupRatio() == 32);
    REQUIRE(fast.lookupSize() <= fast.maxLookupRatio() * fast.nPoints());
    REQUIRE(fast.requestedPolicy() == FastGridPolicy::Automatic);
    REQUIRE(fast.scale() == 1.0);
    REQUIRE(fast.settings().scale == 1.0);
    REQUIRE(fast.dataStatus() == DataStatus::AllocatedHost);

    std::vector<double> queries = {-200.0, -100.0, -99.0, -20.0, -3.0,
                                   -1.0,   -0.5,   0.0,   0.5,   1.0,
                                   3.0,    20.0,   99.0,  100.0, 200.0};
    for (const double point : points) {
      queries.push_back(
          std::nextafter(point, -std::numeric_limits<double>::infinity()));
      queries.push_back(
          std::nextafter(point, std::numeric_limits<double>::infinity()));
    }
    for (int i = 0; i <= 2000; ++i)
      queries.push_back(-150.0 + 300.0 * i / 2000.0);

    for (const double query : queries) {
      REQUIRE(fast.index(query) == binary.index(query));
      int fast_index, binary_index;
      Spiner::weights_t<double> fast_weights, binary_weights;
      fast.weights(query, fast_index, fast_weights);
      binary.weights(query, binary_index, binary_weights);
      REQUIRE(fast_index == binary_index);
      REQUIRE(fast_weights[0] == binary_weights[0]);
      REQUIRE(fast_weights[1] == binary_weights[1]);
    }

    fast.finalize();
    binary.finalize();
  }

  SECTION("Signed-log grids remain exact across scales and sizes") {
    constexpr double transformed_min = -8.0;
    constexpr double transformed_max = 8.0;
    for (const double scale : {0.1, 2.0, 10.0}) {
      for (const int npoints : {3, 17, 65}) {
        std::vector<double> signed_log_points(npoints);
        for (int i = 0; i < npoints; ++i) {
          const double transformed =
              transformed_min + i * (transformed_max - transformed_min) /
                                    static_cast<double>(npoints - 1);
          signed_log_points[i] = scale * nqtSinhForFastGridTest(transformed);
        }
        FastReferenceGrid binary(signed_log_points);
        FastNonUniformGrid1D fast(
            signed_log_points,
            FastGridSettings{.scale = scale,
                             .policy = FastGridPolicy::RequireFast,
                             .max_lookup_ratio = 8});
        REQUIRE(fast.usesFastLookup());
        for (int i = 0; i <= 2000; ++i) {
          const double query =
              1.1 * signed_log_points.front() +
              i * 1.1 * (signed_log_points.back() - signed_log_points.front()) /
                  2000.0;
          REQUIRE(fast.index(query) == binary.index(query));
        }
        fast.finalize();
        binary.finalize();
      }
    }
  }

  SECTION("Policies and host reconfiguration control the lookup table") {
    const std::vector<double> uneven = {0.0, 0.01, 0.02, 1.0};
    FastNonUniformGrid1D automatic(
        uneven, FastGridSettings{.scale = 10.0, .max_lookup_ratio = 8});
    REQUIRE_FALSE(automatic.usesFastLookup());

    FastGridSettings settings = automatic.settings();
    settings.max_lookup_ratio = 32;
    automatic.reconfigureLookup(settings);
    REQUIRE(automatic.usesFastLookup());
    const std::size_t lookup_size = automatic.lookupSize();

    settings.policy = FastGridPolicy::RequireFast;
    automatic.reconfigureLookup(settings);
    REQUIRE(automatic.lookupSize() == lookup_size);
    REQUIRE(automatic.requestedPolicy() == FastGridPolicy::RequireFast);

    settings.policy = FastGridPolicy::Automatic;
    settings.max_lookup_ratio = 8;
    automatic.reconfigureLookup(settings);
    REQUIRE_FALSE(automatic.usesFastLookup());
    REQUIRE(automatic.lookupSize() == 0);

    settings.policy = FastGridPolicy::ForceBinary;
    settings.max_lookup_ratio = 32;
    automatic.reconfigureLookup(settings);
    REQUIRE_FALSE(automatic.usesFastLookup());
    REQUIRE(automatic.requestedPolicy() == FastGridPolicy::ForceBinary);

    settings.policy = FastGridPolicy::Automatic;
    settings.scale = -1.0;
    automatic.reconfigureLookup(settings);
    REQUIRE(automatic.usesFastLookup());
    REQUIRE(automatic.scale() == 0.01);
    automatic.finalize();

    FastNonUniformGrid1D binary(
        points, FastGridSettings{.policy = FastGridPolicy::ForceBinary,
                                 .max_lookup_ratio = 8});
    REQUIRE_FALSE(binary.usesFastLookup());
    REQUIRE(binary.index(3.0) == 5);
    binary.finalize();
  }

  SECTION("Copy, serialization, and device transfer preserve active mode") {
    FastNonUniformGrid1D source(points);
    FastNonUniformGrid1D copy;
    copy.deepCopy(source);
    REQUIRE(copy.usesFastLookup());
    REQUIRE(copy.data() != source.data());
    REQUIRE(copy.lookupSize() == source.lookupSize());
    REQUIRE(copy.index(3.0) == source.index(3.0));

    std::vector<std::byte> serialized(source.serializedSizeInBytes());
    REQUIRE(source.serialize(serialized.data()) == serialized.size());
    FastNonUniformGrid1D restored;
    REQUIRE(restored.deSerialize(serialized.data()) == serialized.size());
    REQUIRE(restored.dataStatus() == DataStatus::Unmanaged);
    REQUIRE(restored.usesFastLookup());
    REQUIRE(restored.index(-3.0) == source.index(-3.0));

    FastNonUniformGrid1D device = source.getOnDevice();
    REQUIRE(device.dataStatus() == DataStatus::AllocatedDevice);
    double error = 0;
    portableReduce(
        "Fast nonuniform device lookup", 0, 1,
        PORTABLE_LAMBDA(const int, double &local_error) {
          const double query = 3.0;
          int ix;
          Spiner::weights_t<double> weights;
          device.weights(query, ix, weights);
          const double reconstructed =
              weights[0] * device.x(ix) + weights[1] * device.x(ix + 1);
          local_error += std::abs(reconstructed - query);
        },
        error);
    REQUIRE(error <= EPSTEST);

    device.finalize();

    FastNonUniformGrid1D binary_host(
        points, FastGridSettings{.policy = FastGridPolicy::ForceBinary,
                                 .max_lookup_ratio = 8});
    FastNonUniformGrid1D binary_device = binary_host.getOnDevice();
    REQUIRE_FALSE(binary_device.usesFastLookup());
    error = 0;
    portableReduce(
        "Fast nonuniform device fallback", 0, 1,
        PORTABLE_LAMBDA(const int, double &local_error) {
          local_error += binary_device.index(3.0) == 5 ? 0.0 : 1.0;
        },
        error);
    REQUIRE(error == 0.0);
    binary_device.finalize();
    binary_host.finalize();

    restored.finalize();
    copy.finalize();
    source.finalize();
  }

  SECTION("shallowCopy creates a non-owning handle to coordinates and table") {
    FastNonUniformGrid1D source(points);
    REQUIRE(source.usesFastLookup());
    FastNonUniformGrid1D shallow;
    shallow.shallowCopy(source);
    REQUIRE(shallow.dataStatus() == DataStatus::Unmanaged);
    REQUIRE(shallow.data() == source.data());
    REQUIRE(shallow.usesFastLookup());
    REQUIRE(shallow.lookupSize() == source.lookupSize());
    REQUIRE(shallow.index(3.0) == source.index(3.0));

    shallow.finalize();
    REQUIRE(source.dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(source.usesFastLookup());
    REQUIRE(source.index(3.0) == 5);
    source.finalize();
  }

  SECTION("DataBox interpolation uses the accelerated grid") {
    const int npoints = static_cast<int>(points.size());
    FastNonUniformDB db(npoints, npoints);
    db.setRange(0, points, FastGridSettings{});
    db.setRange(1, points, FastGridSettings{});
    for (int j = 0; j < npoints; ++j)
      for (int i = 0; i < npoints; ++i)
        db(j, i) = 2.0 * db.range(0).x(i) + 3.0 * db.range(1).x(j) - 1.0;
    REQUIRE(db.range(0).usesFastLookup());
    REQUIRE(db.range(1).usesFastLookup());
    REQUIRE(std::abs(db.interpToReal(-3.0, 3.0) - (-4.0)) <= EPSTEST);

    FastNonUniformDB device = db.getOnDevice();
    double result = 0;
    portableReduce(
        "Fast nonuniform DataBox interpolation", 0, 1,
        PORTABLE_LAMBDA(const int, double &value) {
          value += device.interpToReal(-3.0, 3.0);
        },
        result);
    REQUIRE(std::abs(result - (-4.0)) <= EPSTEST);
    device.finalizeGrids();
    device.finalize();
    db.finalizeGrids();
    db.finalize();
  }
}

TEST_CASE("PiecewiseGrid1D", "[PiecewiseGrid1D]") {
  SECTION("A default piecewise grid has a valid empty lifecycle") {
    PiecewiseGrid1D<3> grid;
    REQUIRE(grid.nGrids() == 0);
    std::vector<std::byte> serialized(grid.serializedSizeInBytes());
    REQUIRE(grid.serialize(serialized.data()) == serialized.size());

    PiecewiseGrid1D<3> restored;
    REQUIRE(restored.deSerialize(serialized.data()) == serialized.size());
    REQUIRE(restored.nGrids() == 0);
    restored.finalize();
    grid.finalize();
  }

  GIVEN("Some regular grid 1Ds") {
    RegularGrid1D g1(0, 0.25, 3);
    RegularGrid1D g2(0.25, 0.75, 11);
    RegularGrid1D g3(0.75, 1, 7);
    THEN("We can construct a piecewise grid object") {
      PiecewiseGrid1D<3> h = {{g1, g2, g3}};
      AND_THEN("We can find each grid based on physical position") {
        REQUIRE(h.findGridFromPosition(0.1) == 0);
        REQUIRE(h.findGridFromPosition(0.3) == 1);
        REQUIRE(h.findGridFromPosition(0.8) == 2);
        // extrapolation
        REQUIRE(h.findGridFromPosition(-5) == 0);
        REQUIRE(h.findGridFromPosition(5) == 2);
      }
      AND_THEN("We can find each grid based on global index") {
        REQUIRE(h.findGridFromGlobalIdx(-1) == 0);
        REQUIRE(h.findGridFromGlobalIdx(0) == 0);
        REQUIRE(h.findGridFromGlobalIdx(2) == 0);
        REQUIRE(h.findGridFromGlobalIdx(3) == 1);
        REQUIRE(h.findGridFromGlobalIdx(4) == 1);
        REQUIRE(h.findGridFromGlobalIdx(13) == 1);
        REQUIRE(h.findGridFromGlobalIdx(14) == 2);
        REQUIRE(h.findGridFromGlobalIdx(20) == 2);
        REQUIRE(h.findGridFromGlobalIdx(21) == 2);
      }
      AND_THEN("We can get x from global index") {
        REQUIRE(std::abs(h.x(2) - 0.25) < EPSTEST);
        REQUIRE(std::abs(h.x(3) - 0.25) < EPSTEST);
      }
      AND_THEN("We can global index from x") {
        REQUIRE(h.index(-1) == 0);
        REQUIRE(h.index(0.25 - 1e-3) == 1);
        REQUIRE(h.index(0.2501) == 3);
        REQUIRE(h.index(100) == 19);
      }
      AND_THEN("We can compute weights") {
        Spiner::weights_t<Real> w;
        int ix;
        h.weights(0.8751, ix, w);
        REQUIRE(ix == 17);
        REQUIRE(std::abs(w[1] - 0.0024) < EPSTEST);
        REQUIRE(std::abs(w[0] - (1 - 0.0024)) < EPSTEST);
      }
      AND_THEN("We can serialize and deserialize the nested grids") {
        const std::size_t expected = sizeof(h) + h.dynamicMemorySizeInBytes();
        REQUIRE(h.serializedSizeInBytes() == expected);
        std::vector<std::byte> serialized(expected);
        REQUIRE(h.serialize(serialized.data()) == expected);

        PiecewiseGrid1D<3> restored;
        REQUIRE(restored.deSerialize(serialized.data()) == expected);
        REQUIRE(restored.nPoints() == h.nPoints());
        for (std::size_t i = 0; i < h.nPoints(); ++i)
          REQUIRE(restored.x(i) == h.x(i));
        REQUIRE(restored.nPoints() == h.nPoints());
        REQUIRE(restored.index(0.8) == h.index(0.8));

        auto device = h.getOnDevice();
        REQUIRE(device.nPoints() == h.nPoints());
        for (std::size_t i = 0; i < h.nPoints(); ++i)
          REQUIRE(device.x(i) == h.x(i));
        device.finalize();
        restored.finalize();
      }
    }
  }
}

TEST_CASE("DataBox Basics", "[DataBox]") {

  SECTION("DataBoxes are initialized with correct rank") {
    DataBox db(2);
    DataBox db4(5, 4, 2, 2);
    REQUIRE(db.rank() == 1);
    REQUIRE(db4.rank() == 4);

    db.finalize(); // free data
    db4.finalize();
  }

  SECTION("A DataBox can be written to and read from") {

    constexpr int M = 3;
    constexpr int N = 2;

    DataBox db(M, N);
    int tot = 0;
    for (int j = 0; j < M; j++) {
      for (int i = 0; i < N; i++) {
        db(j, i) = tot;
        tot++;
      }
    }
    tot = 0;
    for (int j = 0; j < M; j++) {
      for (int i = 0; i < N; i++) {
        REQUIRE(db(j, i) == tot);
        tot++;
      }
    }

    SECTION("DataBox min and max can be correctly computed") {
      REQUIRE(db.max() == tot - 1);
      REQUIRE(db.min() == 0);
    }

    SECTION("DataBox metadata can be copied") {
      DataBox dbCopy;
      dbCopy.copyMetadata(db);
      REQUIRE(dbCopy.rank() == db.rank());
      for (int i = 0; i < db.rank(); i++) {
        REQUIRE(dbCopy.dim(i + 1) == db.dim(i + 1));
        REQUIRE(dbCopy.indexType(i) == db.indexType(i));
      }
      SECTION("DataBoxes can be resized") {
        dbCopy.resize(5, 4, 3);
        REQUIRE(dbCopy.rank() == 3);
        REQUIRE(dbCopy.dim(1) == 3);
        REQUIRE(dbCopy.dim(2) == 4);
        REQUIRE(dbCopy.dim(3) == 5);
      }
      dbCopy.finalize(); // re-allocations require free
    }

    SECTION("DataBoxes can be shallow copied") {
      DataBox db2(db);
      REQUIRE(&(db2(0)) == &(db(0)));
      db2 = db;
      REQUIRE(&(db2(0)) == &(db(0)));
    }

    SECTION("DataBoxes can be deep copied") {
      DataBox db2;
      db2.deepCopy(db);
      REQUIRE(&(db2(0)) != &(db(0)));
      db2.finalize(); // deep copies require free
    }

    SECTION("DataBoxes can be sliced in 2D") {
      DataBox dbslc = db.slice(0);
      DataBox dbslc2(db, 1, 0, 2);

      REQUIRE(dbslc2.rank() == 1);
      REQUIRE(dbslc2.dim(dbslc2.rank()) == 2);

      SECTION("DataBox slices are correctly indexed") {
        int tot = 0;
        for (int i = 0; i < dbslc.dim(dbslc.rank()); i++) {
          REQUIRE(dbslc(i) == tot);
          tot++;
        }
      }

      SECTION("DataBox slices are shallow") {
        REQUIRE(&(dbslc(0)) == &(db(0)));
      }
    }
    db.finalize(); // free original data
  }
}

TEST_CASE("DataBox interpolation", "[DataBox]") {

  GIVEN("A four-dimensional data box filled with a linear function") {
    constexpr int NCOARSE = 5;
    constexpr int NFINE = 20;
    constexpr int RANK = 4;
    DataBox db(Spiner::AllocationTarget::Device, NCOARSE, NCOARSE, NCOARSE,
               NCOARSE);

    constexpr Real xmin = 0;
    constexpr Real xmax = 1;

    for (int i = 0; i < RANK; i++)
      db.setRange(i, xmin, xmax, NCOARSE);

    portableFor(
        "Fill 4D databox", 0, NCOARSE, 0, NCOARSE, 0, NCOARSE, 0, NCOARSE,
        PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                        const int ix) {
          RegularGrid1D grid(xmin, xmax, NCOARSE);
          Real a = grid.x(ia);
          Real z = grid.x(iz);
          Real y = grid.x(iy);
          Real x = grid.x(ix);
          db(ia, iz, iy, ix) = linearFunction(a, z, y, x);
        });
    THEN("interpToReal in 4D is exact for linear functions") {
      Real error = 0;
      portableReduce(
          "Interpolate 4D databox", 0, NFINE, 0, NFINE, 0, NFINE, 0, NFINE,
          PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                          const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            Real a = grid.x(ia);
            Real z = grid.x(iz);
            Real y = grid.x(iy);
            Real x = grid.x(ix);
            Real f_true = linearFunction(a, z, y, x);
            Real difference = db.interpToReal(a, z, y, x) - f_true;
            accumulate += (difference * difference);
          },
          error);
      REQUIRE(error <= EPSTEST);
    }
    THEN("interpToReal in 3D with one index is exact for linear functions") {
      Real error = 0;
      portableReduce(
          "Interpolate + index 4D databox", 0, NFINE, 0, NFINE, 0, NFINE, 0,
          NCOARSE,
          PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                          const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            RegularGrid1D grid_coarse(xmin, xmax, NCOARSE);
            Real a = grid.x(ia);
            Real z = grid.x(iz);
            Real y = grid.x(iy);
            Real x = grid_coarse.x(ix);
            Real f_true = linearFunction(a, z, y, x);
            Real difference = db.interpToReal(a, z, y, ix) - f_true;
            accumulate += (difference * difference);
          },
          error);
      REQUIRE(error <= EPSTEST);
    }
    free(db);
  }

  GIVEN("A data box filled with a linear function and some indexing") {
    constexpr int NCOARSE = 5;
    constexpr int NIDX = 5;
    constexpr int NFINE = 20;
    DataBox db(Spiner::AllocationTarget::Device, NCOARSE, NCOARSE, NCOARSE,
               NIDX, NCOARSE);

    constexpr Real xmin = 0;
    constexpr Real xmax = 1;

    db.setRange(0, xmin, xmax, NCOARSE);
    db.setRange(2, xmin, xmax, NCOARSE);
    db.setRange(3, xmin, xmax, NCOARSE);
    db.setRange(4, xmin, xmax, NCOARSE);
    portableFor(
        "Fill 5D databox", 0, NCOARSE, 0, NCOARSE, 0, NCOARSE, 0, NIDX, 0,
        NCOARSE,
        PORTABLE_LAMBDA(const int ib, const int ia, const int iz, const int iy,
                        const int ix) {
          RegularGrid1D grid(xmin, xmax, NCOARSE);
          Real b = grid.x(ib);
          Real a = grid.x(ia);
          Real z = grid.x(iz);
          Real y = grid.x(iy);
          Real x = grid.x(ix);
          db(ib, ia, iz, iy, ix) = linearFunction(b, a, z, y, x);
        });
    THEN("interpToReal in 4D with one non-interpolated index is exact for "
         "linear functions") {
      Real error = 0;
      portableReduce(
          "Interpolate 5D databox", 0, NFINE, 0, NFINE, 0, NFINE, 0, NIDX, 0,
          NFINE,
          PORTABLE_LAMBDA(const int ib, const int ia, const int iz,
                          const int iy, const int ix, Real &accumulate) {
            RegularGrid1D grid1(xmin, xmax, NFINE);
            RegularGrid1D grid2(xmin, xmax, NIDX);
            Real b = grid1.x(ib);
            Real a = grid1.x(ia);
            Real z = grid1.x(iz);
            Real y = grid2.x(iy);
            Real x = grid1.x(ix);
            Real f_true = linearFunction(b, a, z, y, x);
            Real difference = db.interpToReal(b, a, z, iy, x) - f_true;
            accumulate += (difference * difference);
          },
          error);
      REQUIRE(error <= EPSTEST);
    }
    free(db);
  }

  constexpr int NFINE = 100;
  constexpr int RANK = 3;
  constexpr int NZ = 8;
  constexpr int NY = 10;
  constexpr int NX = 12;
  DataBox db(NZ, NY, NX);

  constexpr Real xmin = 0;
  constexpr Real xmax = 1;
  constexpr Real ymin = -0.5;
  constexpr Real ymax = 0.5;
  constexpr Real zmin = -1;
  constexpr Real zmax = 0;

  std::array<RegularGrid1D, RANK> grids = {RegularGrid1D(xmin, xmax, NX),
                                           RegularGrid1D(ymin, ymax, NY),
                                           RegularGrid1D(zmin, zmax, NZ)};
  std::array<RegularGrid1D, RANK> fine_grids = {
      RegularGrid1D(xmin, xmax, NFINE), RegularGrid1D(ymin, ymax, NFINE),
      RegularGrid1D(zmin, zmax, NFINE)};

  for (int i = 0; i < RANK; i++)
    db.setRange(i, grids[i]);

  for (int iz = 0; iz < NZ; iz++) {
    Real z = grids[2].x(iz);
    for (int iy = 0; iy < NY; iy++) {
      Real y = grids[1].x(iy);
      for (int ix = 0; ix < NX; ix++) {
        Real x = grids[0].x(ix);
        db(iz, iy, ix) = linearFunction(z, y, x);
      }
    }
  }

  SECTION("interpToReal in 3D is exact for linear functions") {
    Real error = 0;
    for (int iz = 0; iz < NFINE; iz++) {
      Real z = fine_grids[2].x(iz);
      for (int iy = 0; iy < NFINE; iy++) {
        Real y = fine_grids[1].x(iy);
        for (int ix = 0; ix < NFINE; ix++) {
          Real x = fine_grids[0].x(ix);
          Real f_true = linearFunction(z, y, x);
          Real difference = db.interpToReal(z, y, x) - f_true;
          error += (difference * difference);
        }
      }
    }
    error = sqrt(error);
    REQUIRE(error <= EPSTEST);
  }

  SECTION("interpToReal in 3D with one non-interpolated index") {
    Real error = 0;
    for (int iz = 0; iz < NFINE; iz++) {
      Real z = fine_grids[2].x(iz);
      for (int iy = 0; iy < NFINE; iy++) {
        Real y = fine_grids[1].x(iy);
        for (int ix = 0; ix < NX; ix++) {
          Real x = grids[0].x(ix);
          Real f_true = linearFunction(z, y, x);
          Real difference = db.interpToReal(z, y, ix) - f_true;
          error += (difference * difference);
        }
      }
    }
    error = sqrt(error);
    REQUIRE(error <= EPSTEST);
  }

  SECTION("interpFromDB 3D->2D") {
    constexpr Real z = (zmax + zmin) / 2.;

    SECTION("Slicing relevant for interpFromDB in slowest index works") {
      int iz = grids[RANK - 1].index(z);
      DataBox lower = db.slice(iz);
      DataBox upper = db.slice(iz + 1);

      Real error = 0;
      for (int iy = 0; iy < NY; iy++) {
        for (int ix = 0; ix < NX; ix++) {
          Real difference = lower(iy, ix) - db(iz, iy, ix);
          error += difference * difference;
          difference = upper(iy, ix) - db(iz + 1, iy, ix);
          error += difference * difference;
        }
      }
      error = sqrt(0.5 * error);
      REQUIRE(error <= EPSTEST);
    }

    // The destination must already exist with the shape of the
    // faster dimensions of the source.
    DataBox db2d(NY, NX);
    Real *const storage = db2d.data();
    db2d.interpFromDB(db, z);

    // Only values are written. Storage, shape, and metadata are the
    // caller's.
    REQUIRE(db2d.data() == storage);
    REQUIRE(db2d.dataStatus() == DataStatus::AllocatedHost);
    REQUIRE(db2d.rank() == 2);
    REQUIRE(db2d.dim(1) == NX);
    REQUIRE(db2d.dim(2) == NY);
    for (int i = 0; i < db2d.rank(); i++) {
      REQUIRE(db2d.indexType(i) == IndexType::Indexed);
    }

    Real error = 0;
    for (int iy = 0; iy < NY; iy++) {
      Real y = grids[1].x(iy);
      for (int ix = 0; ix < NX; ix++) {
        Real x = grids[0].x(ix);
        Real f_true = linearFunction(z, y, x);
        Real difference = db2d(iy, ix) - f_true;
        error += (difference * difference);
      }
    }
    error = sqrt(error);
    REQUIRE(error <= EPSTEST);

    SECTION("interpToReal 2D after the caller sets the grids") {
      for (int i = 0; i < db2d.rank(); i++) {
        db2d.setRange(i, grids[i]);
      }
      Real error = 0;
      for (int iy = 0; iy < NFINE; iy++) {
        Real y = fine_grids[1].x(iy);
        for (int ix = 0; ix < NFINE; ix++) {
          Real x = fine_grids[0].x(ix);
          Real f_true = linearFunction(z, y, x);
          Real difference = db2d.interpToReal(y, x) - f_true;
          error += (difference * difference);
        }
      }
      error = sqrt(error);
      REQUIRE(error <= EPSTEST);
    }

    SECTION("Repeated fills reuse the same destination") {
      constexpr Real z2 = zmin + 0.25 * (zmax - zmin);
      db2d.interpFromDB(db, z2);
      REQUIRE(db2d.data() == storage);
      Real error = 0;
      for (int iy = 0; iy < NY; iy++) {
        Real y = grids[1].x(iy);
        for (int ix = 0; ix < NX; ix++) {
          Real x = grids[0].x(ix);
          Real difference = db2d(iy, ix) - linearFunction(z2, y, x);
          error += (difference * difference);
        }
      }
      error = sqrt(error);
      REQUIRE(error <= EPSTEST);
    }
    free(db2d);
  }

  SECTION("interpFromDB 3D->1D") {
    constexpr Real z = (zmax + zmin) / 2.;
    constexpr Real y = (ymax + ymin) / 2.;

    SECTION("Slicing in 2D works") {
      int iz = grids[RANK - 1].index(z);
      int iy = grids[RANK - 2].index(y);
      DataBox corner = db.slice(iz, iy);
      Real error = 0;
      for (int ix = 0; ix < NX; ix++) {
        error += (corner(ix) - db(iz, iy, ix)) * (corner(ix) - db(iz, iy, ix));
      }
      error = sqrt(error);
      REQUIRE(error <= EPSTEST);
    }

    DataBox db1d(NX);
    Real *const storage = db1d.data();
    db1d.interpFromDB(db, z, y);
    REQUIRE(db1d.data() == storage);
    REQUIRE(db1d.rank() == 1);
    REQUIRE(db1d.dim(1) == NX);
    REQUIRE(db1d.indexType(0) == IndexType::Indexed);

    Real error = 0;
    for (int ix = 0; ix < NX; ix++) {
      Real x = grids[0].x(ix);
      Real f_true = linearFunction(z, y, x);
      Real difference = db1d(ix) - f_true;
      error += difference * difference;
    }
    error = sqrt(error);
    REQUIRE(error <= EPSTEST);
    free(db1d);
  }

  SECTION("interpFromDB fills a device DataBox inside a kernel") {
    constexpr Real z = (zmax + zmin) / 2.;
    DataBox src = db.getOnDevice();
    DataBox dst(Spiner::AllocationTarget::Device, NY, NX);

    portableFor(
        "interpFromDB on device", 0, 1, PORTABLE_LAMBDA(const int) {
          // Captures are const. Copies are shallow, so filling a copy
          // fills dst.
          DataBox fill = dst;
          fill.interpFromDB(src, z);
        });

    Real error = 0;
    portableReduce(
        "Check interpFromDB on device", 0, NY, 0, NX,
        PORTABLE_LAMBDA(const int iy, const int ix, Real &accumulate) {
          RegularGrid1D ygrid(ymin, ymax, NY);
          RegularGrid1D xgrid(xmin, xmax, NX);
          Real difference =
              dst(iy, ix) - linearFunction(z, ygrid.x(iy), xgrid.x(ix));
          accumulate += difference * difference;
        },
        error);
    REQUIRE(std::sqrt(error) <= EPSTEST);
    free(dst);
    free(src);
  }

  free(db); // free databox
}

TEST_CASE("interpFromDB with memory-owning grids",
          "[DataBox][NonUniformGrid1D]") {
  constexpr int NZ = 3;
  constexpr int NY = 4;
  constexpr int NX = 5;
  const std::vector<Real> zpoints = {-1.0, -0.25, 0.0};
  const std::vector<Real> ypoints = {-0.5, -0.1, 0.2, 0.5};
  const std::vector<Real> xpoints = {0.0, 0.1, 0.3, 0.6, 1.0};

  NonUniformDB src(NZ, NY, NX);
  src.setRange(0, xpoints);
  src.setRange(1, ypoints);
  src.setRange(2, zpoints);
  for (int iz = 0; iz < NZ; iz++) {
    for (int iy = 0; iy < NY; iy++) {
      for (int ix = 0; ix < NX; ix++) {
        src(iz, iy, ix) = linearFunction(zpoints[iz], ypoints[iy], xpoints[ix]);
      }
    }
  }

  // The destination owns its own grids, set up by the caller.
  NonUniformDB dst(NY, NX);
  dst.setRange(0, xpoints);
  dst.setRange(1, ypoints);
  const Real *const xdata = dst.range(0).data();
  const Real *const ydata = dst.range(1).data();

  constexpr Real z = -0.5;
  dst.interpFromDB(src, z);

  // interpFromDB must not alias the source grids into the destination.
  REQUIRE(dst.range(0).data() == xdata);
  REQUIRE(dst.range(1).data() == ydata);
  REQUIRE(dst.range(0).data() != src.range(0).data());
  REQUIRE(dst.range(1).data() != src.range(1).data());

  Real error = 0;
  for (int iy = 0; iy < NY; iy++) {
    for (int ix = 0; ix < NX; ix++) {
      Real difference =
          dst(iy, ix) - linearFunction(z, ypoints[iy], xpoints[ix]);
      error += difference * difference;
    }
  }
  REQUIRE(std::sqrt(error) <= EPSTEST);

  // Freeing the destination must leave the source grids intact.
  dst.finalizeGrids();
  free(dst);
  for (int i = 0; i < src.rank(); i++) {
    REQUIRE(src.range(i).dataStatus() == DataStatus::AllocatedHost);
  }
  for (int ix = 0; ix < NX; ix++) {
    REQUIRE(src.range(0).x(ix) == xpoints[ix]);
  }
  constexpr Real y = 0.1;
  constexpr Real x = 0.45;
  REQUIRE(std::abs(src.interpToReal(z, y, x) - linearFunction(z, y, x)) <=
          EPSTEST);
  src.finalizeGrids();
  free(src);
}

TEST_CASE("DataBox Interpolation with piecewise grids",
          "[DataBox][PiecewiseGrid1D]") {
  GIVEN("A piecewise grid") {
    constexpr int NGRIDS = 2;
    constexpr Real xmin = 0;
    constexpr Real xmax = 1;

    RegularGrid1D g1(xmin, 0.35 * (xmax - xmin), 3);
    RegularGrid1D g2(0.35 * (xmax - xmin), xmax, 4);
    PiecewiseGrid1D<NGRIDS> g = {{g1, g2}};

    const int NCOARSE = g.nPoints();

    THEN("The piecewise grid contains a number of points equal the sum of "
         "the points of the individual grids") {
      REQUIRE(g.nPoints() == g1.nPoints() + g2.nPoints());
    }

    WHEN("We construct and fill a 3D DataBox based on this grid") {
      constexpr int RANK = 3;
      PiecewiseDB<NGRIDS> dbh(Spiner::AllocationTarget::Host, NCOARSE, NCOARSE,
                              NCOARSE);
      for (int i = 0; i < RANK; ++i) {
        dbh.setRange(i, g);
      }
      for (int iz = 0; iz < NCOARSE; ++iz) {
        for (int iy = 0; iy < NCOARSE; ++iy) {
          for (int ix = 0; ix < NCOARSE; ++ix) {
            Real x = g.x(ix);
            Real y = g.x(iy);
            Real z = g.x(iz);
            dbh(iz, iy, ix) = linearFunction(z, y, x);
          }
        }
      }
      auto db = dbh.getOnDevice();

      THEN("We can interpolate it to a finer grid and get the right answer") {
        Real error = 0;
        constexpr int NFINE = 21;
        portableReduce(
            "Interpolate 3D databox", 0, NFINE, 0, NFINE, 0, NFINE,
            PORTABLE_LAMBDA(const int iz, const int iy, const int ix,
                            Real &accumulate) {
              RegularGrid1D gfine(xmin, xmax, NFINE);
              Real x = gfine.x(ix);
              Real y = gfine.x(iy);
              Real z = gfine.x(iz);
              Real f_true = linearFunction(z, y, x);
              Real difference = db.interpToReal(z, y, x) - f_true;
              accumulate += (difference * difference);
            },
            error);
        REQUIRE(error <= EPSTEST);
      }

      // cleanup
      free(db);
      free(dbh);
    }

    WHEN("We construct a 3D databox based on this grid, where the slowest "
         "moving index is not interpolatable") {
      constexpr int NSLOW = 3;
      constexpr int RANK = 3;
      PiecewiseDB<NGRIDS> dbh(Spiner::AllocationTarget::Host, NSLOW, NCOARSE,
                              NCOARSE);
      for (int i = 0; i < RANK - 1; ++i) {
        dbh.setRange(i, g);
      }
      dbh.setIndexType(RANK - 1, Spiner::IndexType::Indexed);
      for (int iz = 0; iz < NSLOW; ++iz) {
        for (int iy = 0; iy < NCOARSE; ++iy) {
          for (int ix = 0; ix < NCOARSE; ++ix) {
            Real x = g.x(ix);
            Real y = g.x(iy);
            dbh(iz, iy, ix) = linearFunction(static_cast<Real>(iz), y, x);
          }
        }
      }
      auto db = dbh.getOnDevice();

      THEN("We can do mixed slice interpolation operations") {
        constexpr int NFINE = 21;
        for (int iz = 0; iz < NSLOW; ++iz) {
          auto slc = db.slice(iz);
          Real error = 0;
          portableReduce(
              "Interpolate 2D databox", 0, NFINE, 0, NFINE,
              PORTABLE_LAMBDA(const int iy, const int ix, Real &accumulate) {
                RegularGrid1D gfine(xmin, xmax, NFINE);
                Real x = gfine.x(ix);
                Real y = gfine.x(iy);
                Real f_true = linearFunction(iz, y, x);
                Real difference = slc.interpToReal(y, x) - f_true;
                accumulate += (difference * difference);
              },
              error);
          REQUIRE(error <= EPSTEST);
        }
      }
      free(db);
      free(dbh);
    }
  }
}

SCENARIO("Serializing and deserializing a DataBox",
         "[DataBox][PiecewiseGrid1D][Serialize]") {
  GIVEN("A piecewise grid") {
    constexpr int NGRIDS = 2;
    constexpr Real xmin = 0;
    constexpr Real xmax = 1;

    RegularGrid1D g1(xmin, 0.35 * (xmax - xmin), 3);
    RegularGrid1D g2(0.35 * (xmax - xmin), xmax, 4);
    PiecewiseGrid1D<NGRIDS> g = {{g1, g2}};

    const int NCOARSE = g.nPoints();

    THEN("The piecewise grid contains a number of points equal the sum of "
         "the points of the individual grids") {
      REQUIRE(g.nPoints() == g1.nPoints() + g2.nPoints());
    }

    WHEN("We construct and fill a 3D DataBox based on this grid") {
      constexpr int RANK = 3;
      PiecewiseDB<NGRIDS> dbh(Spiner::AllocationTarget::Host, NCOARSE, NCOARSE,
                              NCOARSE);
      for (int i = 0; i < RANK; ++i) {
        dbh.setRange(i, g);
      }
      for (int iz = 0; iz < NCOARSE; ++iz) {
        for (int iy = 0; iy < NCOARSE; ++iy) {
          for (int ix = 0; ix < NCOARSE; ++ix) {
            Real x = g.x(ix);
            Real y = g.x(iy);
            Real z = g.x(iz);
            dbh(iz, iy, ix) = linearFunction(z, y, x);
          }
        }
      }
      WHEN("We serialize the DataBox") {
        std::size_t serial_size = dbh.serializedSizeInBytes();
        REQUIRE(serial_size == (sizeof(dbh) + dbh.sizeBytes() +
                                RANK * g.dynamicMemorySizeInBytes()));

        std::byte *db_serial = (std::byte *)malloc(serial_size * sizeof(char));
        std::size_t write_offst = dbh.serialize(db_serial);
        REQUIRE(write_offst == serial_size);

        THEN("We can initialize a new databox based on the serialized one") {
          PiecewiseDB<NGRIDS> dbh2;
          std::size_t read_offst = dbh2.deSerialize(db_serial);
          REQUIRE(read_offst == write_offst);

          AND_THEN("They do not point to the same memory") {
            // checks DataBox pointer
            REQUIRE(dbh2.data() != dbh.data());
            // checks accessor agrees
            REQUIRE(&dbh2(0) != &dbh(0));
          }

          WHEN("We initialize a THIRD databox on the serialized one") {
            PiecewiseDB<NGRIDS> dbh3;
            std::size_t read_offst3 = dbh3.deSerialize(db_serial);
            REQUIRE(read_offst3 == write_offst);
            THEN("The second and third databoxes DO point at the same memory") {
              REQUIRE(dbh2.data() == dbh3.data());
              REQUIRE(&dbh3(0) == &dbh2(0));
              AND_THEN("But they are separate objects") {
                REQUIRE(&dbh2 != &dbh3);
              }
            }
          }

          AND_THEN("The shape is correct") {
            REQUIRE(dbh2.rank() == dbh.rank());
            REQUIRE(dbh2.size() == dbh.size());
            for (int d = 1; d <= 3; ++d) {
              REQUIRE(dbh2.dim(d) == dbh.dim(d));
            }
          }

          AND_THEN("The contents are correct") {
            for (int i = 0; i < dbh.size(); ++i) {
              REQUIRE(dbh(i) == dbh2(i));
            }
          }

          AND_THEN("The grid metadata is correct") {
            for (int i = 0; i < RANK; ++i) {
              REQUIRE(dbh2.indexType(i) == IndexType::Interpolated);
              REQUIRE(dbh2.range(i).nPoints() == dbh.range(i).nPoints());
              for (std::size_t j = 0; j < dbh.range(i).nPoints(); ++j)
                REQUIRE(dbh2.range(i).x(j) == dbh.range(i).x(j));
            }
          }
        }

        // cleanup
        free(db_serial);
      }

      // cleanup
      free(dbh);
    }
  }
}

DataBox MakeFilledDB(int N, int &tot) {
  DataBox db(N, N, N);
  tot = 0;
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        db(k, j, i) = tot++;
      }
    }
  }
  db.setRange(0, 0, 1, 10);
  return db;
}
SCENARIO("Moving a databox", "[DataBox]") {
  WHEN("A databox is asigned and the original one goes out of scope") {
    constexpr int N = 2;
    int tot;
    DataBox db;
    {
      auto db2 = MakeFilledDB(N, tot);
      db = db2;
    }
    THEN("The databox status is correct") {
      REQUIRE(db.dataStatus() == Spiner::DataStatus::AllocatedHost);
      AND_THEN("The data is present") {
        REQUIRE(db(N - 1, N - 1, N - 1) == tot - 1);
      }
    }
    free(db);
  }
}

SCENARIO("Decoupling two databoxes", "[DataBox][reset]") {
  GIVEN("A databox") {
    constexpr int N = 5;
    DataBox db(N);
    for (int i = 0; i < db.size(); ++i) {
      db(i) = i;
    }
    WHEN("Another databox is copied from it") {
      DataBox db2 = db;
      THEN("The new databox can be decoupled from the original") {
        db2.reset();
        db2.resize(N);
        for (int i = 0; i < db.size(); ++i) {
          db2(i) = db(i) + 1;
        }
        AND_THEN("The original databox is unchanged") {
          for (int i = 0; i < db.size(); ++i) {
            REQUIRE(std::abs(db(i) - i) <= EPSTEST);
          }
        }
        free(db2);
      }
    }
    free(db);
  }
}

SCENARIO("Allocating a DataBox on device", "[Databox][Constructor]") {
  GIVEN("A databox is allocated on device") {
    constexpr int N = 2;
    constexpr Real factor = 1.275; // something arbitrary
    DataBox db_dev(Spiner::AllocationTarget::Device, N, N, N);
    WHEN("It is set to a given value") {
      portableFor(
          "Fill the databox", 0, N, 0, N, 0, N,
          PORTABLE_LAMBDA(int k, int j, int i) { db_dev(k, j, i) = factor; });
      THEN("That value can be recovered") {
        Real sum = 0;
        portableReduce(
            "Sum up the databox", 0, N, 0, N, 0, N,
            PORTABLE_LAMBDA(int k, int j, int i, Real &val) {
              val += db_dev(k, j, i);
            },
            sum);
        REQUIRE(std::abs(sum - factor * N * N * N) <= EPSTEST);
      }
    }
    free(db_dev);
  }
}

SCENARIO("Copying a DataBox to device", "[DataBox][GetOnDevice]") {
  GIVEN("A databox allocated on the host") {
    constexpr int N = 2;
    constexpr Real factor = 1.275;
    DataBox db_host(N, N, N);
    for (int i = 0; i < db_host.size(); ++i) {
      db_host(i) = factor;
    }
    WHEN("It is copied to device") {
      DataBox db_dev = db_host.getOnDevice();
      THEN("It can be read on device") {
        Real sum = 0;
        portableReduce(
            "Sum up the databox", 0, N, 0, N, 0, N,
            PORTABLE_LAMBDA(int k, int j, int i, Real &val) {
              val += db_dev(k, j, i);
            },
            sum);
        REQUIRE(std::abs(sum - factor * N * N * N) <= EPSTEST);
      }
      printf("free db_dev\n");
      free(db_dev);
    }
    printf("free db_host\n");
    free(db_host);
  }
  GIVEN("An empty databox") {
    DataBox db;
    WHEN("We copy it to device") {
      DataBox db2 = db.getOnDevice();
      THEN("The new object is still empty") {
        REQUIRE(db.dataStatus() == Spiner::DataStatus::Empty);
      }
    }
  }
}

SCENARIO("Using unique pointers to garbage collect DataBox",
         "[DataBox][GarbageCollection]") {
  constexpr int N = 1000;
  GIVEN("A databox allocated on device with a unique pointer") {
    std::unique_ptr<DataBox, DBDeleter<>> pdb(
        new DataBox(Spiner::AllocationTarget::Device, N));
    THEN("We can access it") {
      auto db = *pdb; // shallow copy
      portableFor(
          "Just do something", 0, N,
          PORTABLE_LAMBDA(int i) { db(i) = 2.0 * i; });
    }
  }
  GIVEN("A databox with memory-owning grids and a unique pointer") {
    // DBDeleter<true> also finalizes the grids. A leak checker such as
    // LeakSanitizer reports the coordinates if it does not.
    std::unique_ptr<NonUniformDB, DBDeleter<true>> pdb(new NonUniformDB(4));
    pdb->setRange(0, std::vector<Real>{-1.0, -0.5, 0.25, 2.0});
    THEN("Its grids own their memory") {
      REQUIRE(pdb->range(0).dataStatus() == DataStatus::AllocatedHost);
    }
  }
}

#if SPINER_USE_HDF
TEST_CASE("NonUniformGrid1D HDF5", "[NonUniformGrid1D][HDF5]") {
  const std::vector<Real> points = {-1.0, -0.5, 0.25, 2.0};
  const std::string filename = "nonuniform_grid_test.sp5";
  const std::string grid_name = "grid";
  NonUniformGrid1D grid(points);

  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  herr_t status = grid.saveHDF(file, grid_name);
  status += H5Fclose(file);
  REQUIRE(status == H5_SUCCESS);

  NonUniformGrid1D loaded;
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  status = loaded.loadHDF(file, grid_name);
  status += H5Fclose(file);
  REQUIRE(status == H5_SUCCESS);
  REQUIRE(loaded.nPoints() == grid.nPoints());
  for (std::size_t i = 0; i < grid.nPoints(); ++i)
    REQUIRE(loaded.x(i) == grid.x(i));

  constexpr int N = 4;
  NonUniformDB db(N);
  db.setRange(0, points);
  for (int i = 0; i < N; ++i)
    db(i) = db.range(0).x(i);
  const std::string db_filename = "nonuniform_databox_test.sp5";
  REQUIRE(db.saveHDF(db_filename) == H5_SUCCESS);
  NonUniformDB loaded_db;
  REQUIRE(loaded_db.loadHDF(db_filename) == H5_SUCCESS);
  REQUIRE(loaded_db.range(0).nPoints() == db.range(0).nPoints());
  for (std::size_t i = 0; i < db.range(0).nPoints(); ++i)
    REQUIRE(loaded_db.range(0).x(i) == db.range(0).x(i));
  REQUIRE(std::abs(loaded_db.interpToReal(0.0)) <= EPSTEST);

  loaded_db.finalizeGrids();
  loaded_db.finalize();
  db.finalizeGrids();
  db.finalize();
  loaded.finalize();
  grid.finalize();
}

TEST_CASE("FastNonUniformGrid1D HDF5", "[FastNonUniformGrid1D][HDF5]") {
  const std::vector<double> points = {-100.0, -20.0, -4.0, -1.0, 0.0,
                                      1.0,    4.0,   20.0, 100.0};
  const std::string filename = "fast_nonuniform_grid_test.sp5";
  const std::string grid_name = "grid";
  FastNonUniformGrid1D grid(
      points, FastGridSettings{.scale = 2.0,
                               .policy = FastGridPolicy::RequireFast,
                               .max_lookup_ratio = 8});
  REQUIRE(grid.usesFastLookup());

  hid_t file =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  herr_t status = grid.saveHDF(file, grid_name);
  status += H5Fclose(file);
  REQUIRE(status == H5_SUCCESS);

  FastNonUniformGrid1D loaded;
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  status = loaded.loadHDF(file, grid_name);
  status += H5Fclose(file);
  REQUIRE(status == H5_SUCCESS);
  REQUIRE(loaded.usesFastLookup());
  REQUIRE(loaded.scale() == grid.scale());
  REQUIRE(loaded.maxLookupRatio() == grid.maxLookupRatio());
  REQUIRE(loaded.requestedPolicy() == grid.requestedPolicy());
  REQUIRE(loaded.index(3.0) == grid.index(3.0));

  FastGridSettings settings = loaded.settings();
  settings.policy = FastGridPolicy::ForceBinary;
  loaded.reconfigureLookup(settings);
  REQUIRE_FALSE(loaded.usesFastLookup());
  REQUIRE(loaded.index(3.0) == grid.index(3.0));

  loaded.finalize();
  grid.finalize();
}

SCENARIO("PiecewiseGrid HDF5", "[PiecewiseGrid1D][HDF5]") {
  GIVEN("A piecewise grid") {
    RegularGrid1D g1(0, 0.25, 3);
    RegularGrid1D g2(0.25, 0.75, 11);
    RegularGrid1D g3(0.75, 1, 7);
    PiecewiseGrid1D<3> piecewise_grid = {{g1, g2, g3}};
    THEN("We can save it to file") {
      const std::string filename = "piecewise_test.sp5";
      const std::string grid_name = "grid";
      herr_t status;
      hid_t file;
      file =
          H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      status = piecewise_grid.saveHDF(file, grid_name.c_str());
      status += H5Fclose(file);
      REQUIRE(status == H5_SUCCESS);

      AND_THEN("We can read it from file and get the same information out") {
        PiecewiseGrid1D<3> loaded_grid;
        herr_t status;
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        status = loaded_grid.loadHDF(file, grid_name.c_str());
        status += H5Fclose(file);
        REQUIRE(status == H5_SUCCESS);

        REQUIRE(loaded_grid.nPoints() == piecewise_grid.nPoints());
        for (std::size_t i = 0; i < piecewise_grid.nPoints(); ++i)
          REQUIRE(loaded_grid.x(i) == piecewise_grid.x(i));
      }
    }
    GIVEN("A single regular grid") {
      RegularGrid1D g1(0, 0.25, 3);
      WHEN("We save it to file") {
        const std::string filename = "backwards_compatibility_test.sp5";
        const std::string grid_name = "grid";
        herr_t status;
        hid_t file;
        file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
        status = g1.saveHDF(file, grid_name.c_str());
        status += H5Fclose(file);
        REQUIRE(status == H5_SUCCESS);

        THEN("We can read the file back out with a piecewise grid") {
          PiecewiseGrid1D<3> loaded_grid;
          herr_t status;
          hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
          status = loaded_grid.loadHDF(file, grid_name.c_str());
          status += H5Fclose(file);
          REQUIRE(status == H5_SUCCESS);

          REQUIRE(loaded_grid.min() == g1.min());
          REQUIRE(loaded_grid.max() == g1.max());
          REQUIRE(loaded_grid.nPoints() == g1.nPoints());
          REQUIRE(loaded_grid.nGrids() == 1);
        }
      }
    }
  }
}

SCENARIO("DataBox HDF5", "[DataBox][HDF5]") {
  constexpr int N = 2;
  herr_t status;

  DataBox db(N, N, N);
  int tot = 0;
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        db(k, j, i) = tot++;
      }
    }
  }
  db.setRange(0, 0, 1, 10);

  GIVEN("DataBox can be saved to HDF5") {
    status = db.saveHDF();
    REQUIRE(status == H5_SUCCESS);

    WHEN("DataBox can be loaded from HDF5") {
      DataBox db2;
      status = db2.loadHDF();
      REQUIRE(status == H5_SUCCESS);

      THEN("Metadata read in is consistent") {
        REQUIRE(db2.rank() == db.rank());
        for (int i = 0; i < db.rank(); i++) {
          REQUIRE(db.indexType(i) == db2.indexType(i));
          REQUIRE(db.dim(i + 1) == db2.dim(i + 1));
          if (db.indexType(i) == IndexType::Interpolated) {
            REQUIRE(db.range(i).nPoints() == db2.range(i).nPoints());
            for (std::size_t j = 0; j < db.range(i).nPoints(); ++j)
              REQUIRE(db.range(i).x(j) == db2.range(i).x(j));
          }
        }
        AND_THEN("Data itself is consistent") {
          for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
              for (int i = 0; i < N; i++) {
                REQUIRE(db(k, j, i) == db2(k, j, i));
              }
            }
          }
        }
      }
      free(db2);
    }
  }
  free(db);
}
#endif

#ifdef PORTABILITY_STRATEGY_KOKKOS
SCENARIO("Kokkos functionality: interpolation", "[DataBox],[Kokkos]") {
  constexpr int NFINE = 100;
  constexpr int RANK = 3;
  constexpr int NZ = 8;
  constexpr int NY = 9;
  constexpr int NX = 12;
  DataBox db(NZ, NY, NX);

  constexpr Real xmin = 0;
  constexpr Real xmax = 1;
  constexpr Real ymin = -0.5;
  constexpr Real ymax = 0.5;
  constexpr Real zmin = -1;
  constexpr Real zmax = 0;

  std::array<RegularGrid1D, RANK> grids = {RegularGrid1D(xmin, xmax, NX),
                                           RegularGrid1D(ymin, ymax, NY),
                                           RegularGrid1D(zmin, zmax, NZ)};

  Kokkos::View<RegularGrid1D *> fine_grids("fine grids", RANK);
  auto fine_grids_h = Kokkos::create_mirror_view(fine_grids);
  fine_grids_h[0] = RegularGrid1D(xmin, xmax, NFINE);
  fine_grids_h[1] = RegularGrid1D(ymin, ymax, NFINE);
  fine_grids_h[2] = RegularGrid1D(zmin, zmax, NFINE);
  Kokkos::deep_copy(fine_grids, fine_grids_h);

  for (int i = 0; i < RANK; i++)
    db.setRange(i, grids[i]);

  for (int iz = 0; iz < NZ; iz++) {
    Real z = grids[2].x(iz);
    for (int iy = 0; iy < NY; iy++) {
      Real y = grids[1].x(iy);
      for (int ix = 0; ix < NX; ix++) {
        Real x = grids[0].x(ix);
        db(iz, iy, ix) = linearFunction(z, y, x);
      }
    }
  }

  using DeviceView_t = Kokkos::View<Real *, Kokkos::MemoryUnmanaged>;
  using HostView_t =
      Kokkos::View<Real *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  Real *device_data = (Real *)PORTABLE_MALLOC(db.sizeBytes());
  DeviceView_t deviceView(device_data, db.size());
  HostView_t hostView(db.data(), db.size());
  Kokkos::deep_copy(deviceView, hostView);
  DataBox db_dev(device_data, NZ, NY, NX);
  db_dev.copyShape(db);

  Real error = 0;
  using Policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
  Kokkos::parallel_reduce(
      Policy3D({0, 0, 0}, {NFINE, NFINE, NFINE}),
      PORTABLE_LAMBDA(const int iz, const int iy, const int ix, Real &update) {
        DataBox db2 = db_dev; // checks that copying works on device
        const Real z = fine_grids[2].x(iz);
        const Real y = fine_grids[1].x(iy);
        const Real x = fine_grids[0].x(ix);
        const Real f_true = linearFunction(z, y, x);
        const Real difference = db2.interpToReal(z, y, x) - f_true;
        update += difference * difference;
      },
      error);
  error = sqrt(error);
  REQUIRE(error <= EPSTEST);

  PORTABLE_FREE(device_data);
  free(db);
}
#endif

// A three-point quadratic rule, optionally limited to the range of the
// two samples bracketing x. Exercises interp_toy with a nonlinear
// reducer and a stencil that must be shifted inward at the upper edge.
template <class Grid>
PORTABLE_INLINE_FUNCTION auto limitedQuadratic(const Grid &grid, const Real x,
                                               const bool limit) {
  int ix;
  Spiner::weights_t<Real> w;
  grid.weights(x, ix, w);
  const int n = static_cast<int>(grid.nPoints());
  const int s = (ix + 2 < n) ? ix : n - 3;
  const int k = ix - s; // offset of the left bracketing sample
  const Real x0 = grid.x(s);
  const Real x1 = grid.x(s + 1);
  const Real x2 = grid.x(s + 2);
  const Real l0 = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
  const Real l1 = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
  const Real l2 = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
  return Spiner::interp::make_interpolation(
      std::array<int, 3>{s, s + 1, s + 2}, [=](const auto &v) {
        auto q = l0 * v[0] + l1 * v[1] + l2 * v[2];
        if (limit) {
          const auto lo = v[k] < v[k + 1] ? v[k] : v[k + 1];
          const auto hi = v[k] < v[k + 1] ? v[k + 1] : v[k];
          q = q < lo ? lo : (q > hi ? hi : q);
        }
        return q;
      });
}

TEST_CASE("interp_toy interpolate_with", "[DataBox][interp_toy]") {
  using Spiner::interp::at;
  using Spiner::interp::interpolate_with;
  using Spiner::interp::linear;

  constexpr Real xmin = 0;
  constexpr Real xmax = 1;

  GIVEN("A four-dimensional data box filled with a linear function") {
    constexpr int NCOARSE = 5;
    constexpr int NFINE = 12;
    constexpr int RANK = 4;
    DataBox db(Spiner::AllocationTarget::Device, NCOARSE, NCOARSE, NCOARSE,
               NCOARSE);
    for (int i = 0; i < RANK; i++)
      db.setRange(i, xmin, xmax, NCOARSE);

    portableFor(
        "Fill 4D databox", 0, NCOARSE, 0, NCOARSE, 0, NCOARSE, 0, NCOARSE,
        PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                        const int ix) {
          RegularGrid1D grid(xmin, xmax, NCOARSE);
          db(ia, iz, iy, ix) =
              linearFunction(grid.x(ia), grid.x(iz), grid.x(iy), grid.x(ix));
        });

    THEN("Fully linear interpolation is exact and matches interpToReal") {
      Real error = 0;
      Real mismatch = 0;
      portableReduce(
          "interpolate_with 4D", 0, NFINE, 0, NFINE, 0, NFINE, 0, NFINE,
          PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                          const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            const Real a = grid.x(ia);
            const Real z = grid.x(iz);
            const Real y = grid.x(iy);
            const Real x = grid.x(ix);
            const Real v = interpolate_with(
                db, linear(db.range(3), a), linear(db.range(2), z),
                linear(db.range(1), y), linear(db.range(0), x));
            const Real difference = v - linearFunction(a, z, y, x);
            accumulate += difference * difference;
          },
          error);
      portableReduce(
          "interpolate_with vs interpToReal 4D", 0, NFINE, 0, NFINE, 0, NFINE,
          0, NFINE,
          PORTABLE_LAMBDA(const int ia, const int iz, const int iy,
                          const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            const Real a = grid.x(ia);
            const Real z = grid.x(iz);
            const Real y = grid.x(iy);
            const Real x = grid.x(ix);
            const Real v = interpolate_with(
                db, linear(db.range(3), a), linear(db.range(2), z),
                linear(db.range(1), y), linear(db.range(0), x));
            const Real difference = v - db.interpToReal(a, z, y, x);
            accumulate += difference * difference;
          },
          mismatch);
      REQUIRE(error <= EPSTEST);
      REQUIRE(mismatch <= EPSTEST);
    }

    THEN("A fixed index may sit in any dimension") {
      Real error = 0;
      portableReduce(
          "interpolate_with 4D with indices", 0, NCOARSE, 0, NFINE, 0, NFINE,
          0, NFINE,
          PORTABLE_LAMBDA(const int i, const int iz, const int iy,
                          const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            RegularGrid1D coarse(xmin, xmax, NCOARSE);
            const Real c = coarse.x(i);
            const Real z = grid.x(iz);
            const Real y = grid.x(iy);
            const Real x = grid.x(ix);
            const auto lz = linear(db.range(2), z);
            const auto ly = linear(db.range(1), y);
            const auto lx = linear(db.range(0), x);
            const Real slowest =
                interpolate_with(db, at(i), lz, ly, lx) -
                linearFunction(c, z, y, x);
            const Real middle =
                interpolate_with(db, linear(db.range(3), z), lz, at(i), lx) -
                linearFunction(z, z, c, x);
            const Real fastest =
                interpolate_with(db, linear(db.range(3), z), lz, ly, at(i)) -
                db.interpToReal(z, z, y, i);
            accumulate +=
                slowest * slowest + middle * middle + fastest * fastest;
          },
          error);
      REQUIRE(error <= EPSTEST);
    }
    free(db);
  }

  GIVEN("A 2D data box holding a function quadratic in x") {
    constexpr int NX = 11;
    constexpr int NY = 4;
    constexpr int NFINE = 37;
    DataBox db(Spiner::AllocationTarget::Device, NY, NX);
    db.setRange(0, xmin, xmax, NX);
    db.setRange(1, xmin, xmax, NY);
    portableFor(
        "Fill 2D databox", 0, NY, 0, NX,
        PORTABLE_LAMBDA(const int iy, const int ix) {
          RegularGrid1D gx(xmin, xmax, NX);
          RegularGrid1D gy(xmin, xmax, NY);
          const Real x = gx.x(ix);
          db(iy, ix) = gy.x(iy) + x * x;
        });

    THEN("Mixing linear and limited quadratic rules is exact") {
      Real error = 0;
      portableReduce(
          "interpolate_with limited quadratic", 0, NFINE, 0, NFINE,
          PORTABLE_LAMBDA(const int iy, const int ix, Real &accumulate) {
            RegularGrid1D grid(xmin, xmax, NFINE);
            const Real y = grid.x(iy);
            const Real x = grid.x(ix);
            const Real v =
                interpolate_with(db, linear(db.range(1), y),
                                 limitedQuadratic(db.range(0), x, true));
            const Real difference = v - (y + x * x);
            accumulate += difference * difference;
          },
          error);
      REQUIRE(error <= EPSTEST);
    }
    free(db);
  }

  GIVEN("A 1D data box with a spike") {
    constexpr int NX = 11;
    DataBox db(NX);
    db.setRange(0, xmin, xmax, NX);
    for (int ix = 0; ix < NX; ix++) {
      db(ix) = (ix == 5) ? 1 : 0;
    }
    // Stencil {3, 4, 5} sees the spike; x is bracketed by two zeros.
    const Real x = 0.35;
    THEN("The unlimited quadratic undershoots") {
      REQUIRE(interpolate_with(db, limitedQuadratic(db.range(0), x, false)) <
              -0.1);
    }
    THEN("The limiter keeps the result within the bracketing samples") {
      REQUIRE(interpolate_with(db, limitedQuadratic(db.range(0), x, true)) ==
              0);
    }
    free(db);
  }

  SECTION("The result type follows the data type") {
    using FloatDB = Spiner::DataBox<float>;
    using Spiner::RegularGrid1D;
    STATIC_REQUIRE(
        std::is_same_v<decltype(interpolate_with(
                           std::declval<const FloatDB &>(),
                           linear(std::declval<RegularGrid1D<float>>(), 0.f),
                           at(0))),
                       float>);
  }
}

int main(int argc, char *argv[]) {

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  int result;
  {
    result = Catch::Session().run(argc, argv);
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}
