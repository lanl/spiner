# Pointer-Owning Non-Uniform Grids

Generative AI was used to assist with this document.

## Summary

Introduce `NonUniformGrid1D<T>`, storing an ordered coordinate array and using binary search for interpolation. Each `DataBox<T, NonUniformGrid1D<T>>` axis owns or references its own coordinates. Existing regular and piecewise grids remain unchanged, and mixed grid classes within one `DataBox` remain out of scope.

## Public Interface and Behavior

- Add `spiner/nonuniform_grid_1d.hpp` with the grid interface expected by `DataBox`: default, owning vector/list, and unmanaged-pointer construction; coordinate lookup; interpolation weights; lifecycle; serialization; device copying; and HDF5 I/O.
- Require at least two finite, strictly increasing coordinates. Find intervals with a device-portable `O(log N)` binary search and preserve regular-grid extrapolation by clamping the interval to `[0, N-2]`.
- Copies retain Spiner's shallow pointer semantics. They do not copy or reference-count coordinate storage, so callers finalize an owned allocation exactly once.

## Ownership and Integration

- Store point count, pointer, and `DataStatus`. Vector/list construction owns host memory; pointer construction is unmanaged; `getOnDevice()` deep-copies coordinates; and idempotent `finalize()` releases only owned host/device storage.
- Serialize `[NonUniformGrid1D header][coordinate values]`; relocation makes the coordinate pointer unmanaged.
- Use `DataBox<T, NonUniformGrid1D<T>>` without `DataBox` special cases: its existing lifecycle delegation applies independently to every axis.
- Persist coordinates as a one-dimensional HDF5 dataset and document ownership, search, and extrapolation behavior.

## Test Plan

- Test lookup, interpolation, extrapolation, invalid coordinate sequences, ownership, copies, and idempotent cleanup.
- Test standalone and multi-axis `DataBox` serialization, device copies, portable-kernel interpolation, and HDF5 round trips.
- Run existing regular-grid, piecewise-grid, DataBox lifecycle, serialization, HDF5, and accelerator regressions.

## Assumptions

- The public type is `Spiner::NonUniformGrid1D<T>`.
- Borrowed coordinates remain valid and unchanged for the grid lifetime.
- Owning construction uses host storage; `getOnDevice()` provides host-to-device
  copies.
- Mixing regular, piecewise, and non-uniform classes in one `DataBox` is deferred.
