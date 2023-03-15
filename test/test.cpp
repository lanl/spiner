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

#include <algorithm> // std::min, std::max
#include <array>
#include <cmath> // sqrt
#include <memory>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/spiner_types.hpp>

#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

using DataBox = Spiner::DataBox<Real>;
using Spiner::IndexType;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Spiner::DBDeleter;
const Real EPSTEST = std::sqrt(DataBox::EPS);

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

TEST_CASE("PortableMDArrays can be allocated from a pointer",
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
    REQUIRE(g.dx() == (max - min) / ((Real)(N - 1)));
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
      REQUIRE(dbCopy != db);

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
      db2.copy(db);
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
        REQUIRE(dbslc == dbslc2);
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

  SECTION("interpFromDB 3D->2D") {
    constexpr Real z = (zmax + zmin) / 2.;

    SECTION("Slicing relevant for interpToDB in slowest index works") {
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

    DataBox db2d;
    db2d.resize(db.size() / db.dim(db.rank()));
    db2d.interpFromDB(db, z);

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

    SECTION("interpToReal 2D") {
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

    DataBox db1d;
    db1d.resize(db.size() / (db.dim(db.rank()) * db.dim(db.rank() - 1)));
    db1d.interpFromDB(db, z, y);
    REQUIRE(db1d.rank() == 1);
    REQUIRE(db1d.dim(1) == NX);

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
  free(db); // free databox
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
      free(db_dev);
    }
    free(db_host);
  }
}

SCENARIO("Using unique pointers to garbage collect DataBox",
         "[DataBox][GarbageCollection]") {
  constexpr int N = 1000;
  GIVEN("A databox allocated on device with a unique pointer") {
    std::unique_ptr<DataBox, DBDeleter> pdb(
        new DataBox(Spiner::AllocationTarget::Device, N));
    THEN("We can access it") {
      auto db = *pdb; // shallow copy
      portableFor(
          "Just do something", 0, N,
          PORTABLE_LAMBDA(int i) { db(i) = 2.0 * i; });
    }
  }
}

#if SPINER_USE_HDF
SCENARIO("DataBox HDF5", "[DataBox],[HDF5]") {
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
            REQUIRE(db.range(i) == db2.range(i));
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

int main(int argc, char *argv[]) {

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  int result;
  { result = Catch::Session().run(argc, argv); }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}
