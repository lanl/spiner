# Changelog

## Current develop

### Added (new features/APIs/variables/...)
- [[MR 148]](https://github.com/lanl/spiner/pull/148) Add `FastNonUniformGrid1D` with exact O(1) NQT-asinh lookup acceleration, configurable binary fallback, and host-side reconfiguration.
- [[MR 146]](https://github.com/lanl/spiner/pull/146) Add UnstructuredGrid1D. Note this bumps ports-of-call dependency and the next spiner release should depend on the newer ports-of-call.
- [[MR 144]](https://github.com/lanl/spiner/pull/144) Add the ability for grid objects to own memory. This extends the public API but the real changes will come later.
- [[MR 103]](https://github.com/lanl/spiner/pull/103) Add config_summary to cmake

### Fixed (Repair bugs, etc)
- [[MR 150]](https://github.com/lanl/spiner/pull/150) Cleaned up memory model for grids that own data
- [[MR 141]](https://github.com/lanl/spiner/pull/142) fix issue with range() method on device

### Changed (changing behavior/API/variables/...)
- [[MR 150]](https://github.com/lanl/spiner/pull/150) Replace the grid `copy()` method with `shallowCopy()`, which makes a non-owning (`Unmanaged`) handle, and `deepCopy()`, the previous `copy()`. Rename `DataBox::copy()` to `DataBox::deepCopy()`, which now deep copies grids too.
- [[MR 150]](https://github.com/lanl/spiner/pull/150) `DataBox::finalize()` and `free()` no longer finalize grids. Use the new `DataBox::finalizeGrids()` to free grid-owned memory. `DBDeleter` is now `DBDeleter<bool FinalizeGrids = false>`, so write `DBDeleter<>` in place of `DBDeleter`. Slices, `copyShape()`, and `getOnDevice(false)` now hold non-owning grids, and `DataBox` assignment no longer deep copies grids.
- [[MR 150]](https://github.com/lanl/spiner/pull/150) `DataBox::interpFromDB` now requires a destination that already has storage and the right shape, and writes only values, leaving the destination's grids and index types unchanged.
- [[MR 145]](https://github.com/lanl/spiner/pull/145) Remove comparator and dx methods from public API

### Infrastructure (changes irrelevant to downstream codes)
- [[MR 141]](https://github.com/lanl/spiner/pull/142) update contribution rules

### Deprecated (soon to be removed behavior/API/variables/...)
- [[MR 150]](https://github.com/lanl/spiner/pull/150) Removed interpToDB

## Release 1.7.0
Date: 12/15/2026

This is the start of changelog

© 2021-2026. Triad National Security, LLC. All rights reserved.  This
program was produced under U.S. Government contract 89233218CNA000001
for Los Alamos National Laboratory (LANL), which is operated by Triad
National Security, LLC for the U.S.  Department of Energy/National
Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of
Energy/National Nuclear Security Administration. The Government is
granted for itself and others acting on its behalf a nonexclusive,
paid-up, irrevocable worldwide license in this material to reproduce,
prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.
