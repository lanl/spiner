# Grid Resource Lifecycle and Serialization Plan

## Goal

Extend the grid API so that `DataBox` can manage grids that may eventually own
dynamic host or device data. Add serialization, pointer relocation, device-copy,
and finalization operations to the current grid types, and thread those
operations through `PiecewiseGrid1D` and `DataBox`.

The current grids contain only inline data, so most initial implementations are
trivial. The design must nevertheless exercise the complete ownership path so a
future dynamically allocated grid does not require another `DataBox` API
change.

## Resource-lifecycle policy

Apply grid serialization, device copying, pointer relocation, and finalization
to every grid slot in `[0, rank)`, regardless of its `IndexType`.

`IndexType` controls interpolation semantics, not resource lifecycle. Processing
all grid slots:

- keeps the lifecycle loops simple and consistent;
- prevents an inactive grid from retaining an unrelocated host pointer in a
  device-side `DataBox`;
- prevents leaks when a dimension changes from interpolated to indexed or
  named;
- avoids making resource ownership depend on mutable interpolation metadata.

This requires every default-constructed grid to be a valid empty object.
`PiecewiseGrid1D` must therefore initialize `NGRIDS_` to zero and initialize any
other state required by its empty operations.

The following invariants should hold for every supported grid:

- default construction creates a valid empty grid;
- empty-grid serialization and device copying are valid;
- `setPointer()` may consume zero bytes;
- `finalize()` is idempotent and safe on empty or unmanaged storage;
- `getOnDevice()` performs a deep copy of any owned allocation;
- deserialized storage is unmanaged and refers into the supplied byte buffer.

## Public grid API

Add the following interface to `RegularGrid1D` and `PiecewiseGrid1D`:

```cpp
std::size_t serializedSizeInBytes() const;
std::size_t serialize(char *dst) const;
std::size_t deSerialize(char *src);
std::size_t setPointer(char *src);
GridType getOnDevice() const;
void finalize();
```

Preserve the existing `deSerialize` capitalization used by `DataBox` unless a
separate project-wide API rename is undertaken.

### Serialization contract

Each serializable object uses this recursive representation:

```text
[object bytes][owned payload or recursively serialized child objects]
```

- `serializedSizeInBytes()` includes `sizeof(*this)` and all subordinate
  serialized storage.
- `serialize()` writes the complete representation and returns the number of
  bytes written.
- `deSerialize()` copies the object bytes, calls `setPointer()` on the remainder,
  and returns the total number of bytes consumed.
- `setPointer()` assumes the object's inline metadata has already been restored.
  It relocates owned pointers into the supplied buffer and processes subordinate
  records, returning the number of payload bytes consumed.

Offsets returned by every method are accumulated by callers. Tests must verify
that calculated, written, and consumed sizes agree.

## `RegularGrid1D` implementation

Implement the new operations in `spiner/regular_grid_1d.hpp`.

Because all current state is inline:

- `serializedSizeInBytes()` returns `sizeof(*this)`;
- `serialize()` copies `sizeof(*this)` bytes to `dst`;
- `deSerialize()` copies `sizeof(*this)` bytes from `src`, then calls
  `setPointer(src + sizeof(*this))`;
- `setPointer()` returns zero;
- `getOnDevice()` returns a value copy of `*this`;
- `finalize()` is a no-op.

Keep host-only memory-management routines undecorated unless the build requires
a portability annotation. Interpolation and access methods remain portable.

## `PiecewiseGrid1D` implementation

Implement the new operations in `spiner/piecewise_grid_1d.hpp`.

1. Make the default constructor establish an empty state, including
   `NGRIDS_ == 0`.
2. Serialize the `PiecewiseGrid1D` object followed by a complete serialized
   record for each child in `[0, NGRIDS_)`.
3. In `setPointer()`, deserialize each child record in sequence. Calling the
   child's `deSerialize()` is appropriate because each child record starts with
   that child's object bytes.
4. In `getOnDevice()`, copy inline metadata and replace every active child with
   `grids_[i].getOnDevice()`.
5. In `finalize()`, call `finalize()` on every child in `[0, NGRIDS_)`, then
   leave the piecewise grid in a valid empty/finalized state.

This recursive implementation is intentionally nontrivial even though
`RegularGrid1D` currently owns no allocation. It ensures that a future
pointer-owning child grid is handled automatically.

Review empty-grid query behavior separately from lifecycle behavior. Methods
such as `min()`, `max()`, and `nPoints()` may continue to require a nonempty
grid, while serialization, device copying, and finalization must support an
empty grid.

## `DataBox` implementation

Update the existing lifecycle methods in `spiner/databox.hpp`.

### Serialized representation

Use the following layout:

```text
[DataBox bytes][DataBox value data][serialized grid 0]...[serialized grid rank-1]
```

The grid records duplicate the inline grid metadata present in the raw
`DataBox` image. This small cost provides a uniform, independently testable grid
serialization contract and permits each grid to restore its own future
allocation without `DataBox` knowing its representation.

### Size and serialization

- Extend `serializedSizeInBytes()` by summing
  `grids_[i].serializedSizeInBytes()` for every `i` in `[0, rank)`.
- Keep the prohibition against serializing device-resident `DataBox` values.
- After writing the existing `DataBox` header and value buffer, call
  `grids_[i].serialize()` for every dimension and accumulate each returned
  offset.

### Deserialization and pointer relocation

- Preserve the existing precondition that an owning, active `DataBox` cannot be
  overwritten by deserialization.
- Copy the `DataBox` object bytes first so its dimensions, rank, and grid
  metadata are available.
- Have `DataBox::setPointer()` relocate the value pointer and rebuild the
  `PortableMDArray`, then call `grids_[i].deSerialize()` for every dimension.
- Mark the resulting value storage unmanaged, as today.
- Ensure each deserialized grid likewise regards memory inside the serialization
  buffer as unmanaged.

The byte buffer must outlive every `DataBox` deserialized from it.

### Device copy

Retain the current deep copy of the value buffer. After copying shape and index
metadata, replace every grid in `[0, rank)` with the result of that grid's
`getOnDevice()`.

Avoid leaving a shallow host copy of any grid allocation in the returned
device-side `DataBox`, including for indexed or named dimensions.

### Finalization

Call `finalize()` on every grid in `[0, rank)` before freeing the `DataBox`
value allocation. Then clear the `DataBox` pointer/view state consistently with
its empty status.

Grid finalization must be safe when `DataBox::finalize()` is reached during
allocation or resize of a default/empty grid. Preserve the current explicit
manual-lifetime model; introducing RAII ownership is outside this change.

### Range replacement and index changes

Since every slot participates in lifecycle management, stale grids will still
be finalized by the owning `DataBox`. Nevertheless, review these mutation
paths:

- replacing a grid with `setRange()`;
- changing an interpolated dimension to indexed or named;
- `reset()`, shallow copy, slicing, assignment, and `copyShape()`.

Document and test which object remains responsible for finalizing shallow grid
allocations. Do not silently introduce a second owner through ordinary copy or
assignment. A future owning grid will likely need an ownership-status field
parallel to `DataStatus`; its normal copies should be shallow/unmanaged, while
`getOnDevice()` is explicitly deep.

## HDF5 behavior

Do not change the HDF5 representation as part of this work.

HDF5 represents meaningful interpolation metadata rather than the complete
in-memory object graph, so `DataBox::saveHDF()` and `loadHDF()` may continue to
visit only interpolated dimensions. They already delegate grid persistence to
the grid's `saveHDF()` and `loadHDF()` methods. A future grid that owns data must
extend those methods independently.

Run the existing HDF5 round-trip tests as regressions.

## Serialization compatibility documentation

Cross-version serialization compatibility is explicitly not a requirement. No
version marker or compatibility reader is needed.

Update `doc/sphinx/src/databox.rst` and add a concise comment next to the
serialization implementation stating:

> DataBox serialization is intended for transient communication between
> processes using the same Spiner version and a compatible architecture. The
> binary representation is not guaranteed to be compatible across Spiner
> versions, architectures, compilers, or build configurations. It should not be
> used as a persistent storage format or a web interchange format. Use HDF5 for
> persistent, portable storage.

Also document:

- the revised `DataBox` and nested-grid layout;
- that grid records are emitted for every dimension in `[0, rank)`, including
  indexed and named dimensions;
- that deserialized value and grid storage refer into an externally owned
  buffer;
- that the buffer must remain alive and suitably aligned for the lifetime of
  the deserialized object;
- that size calculation, serialization, and deserialization must use a
  compatible Spiner build.

The existing architecture and alignment limitations remain unchanged.

## Testing plan

### 1. `RegularGrid1D` unit tests

- Round-trip a populated regular grid.
- Verify calculated, written, and consumed sizes are all `sizeof(grid)`.
- Verify equality, `x()`, `index()`, and `weights()` after deserialization.
- Exercise `setPointer()`, `getOnDevice()`, and repeated `finalize()` calls.
- Exercise the same lifecycle methods on a default-constructed grid.

### 2. `PiecewiseGrid1D` unit tests

- Verify a default piecewise grid is safely empty and supports all lifecycle
  operations.
- Round-trip a multi-segment piecewise grid.
- Verify its serialized size includes its own header and every child record.
- Verify equality, point count, coordinate lookup, indexing, weights, and
  interpolation behavior after deserialization.
- Verify device copying and repeated finalization recursively process children.

### 3. Existing `DataBox` serialization tests

Extend the serialization scenario in `test/test.cpp`:

- update the expected size to include one serialized grid record per dimension;
- verify calculated, written, and consumed sizes agree;
- verify rank, dimensions, values, index types, and all grid ranges after
  deserialization;
- retain the test that independently deserialized `DataBox` objects refer to
  the same serialized buffer while remaining distinct C++ objects;
- add a mixture of interpolated, indexed, and named dimensions and verify all
  grid records are processed;
- test empty/default grid slots explicitly.

### 4. Device integration tests

Extend the `DataBox::getOnDevice()` tests:

- assign grids before copying;
- perform interpolation inside a portable kernel;
- verify both the values and grid data are usable on the device;
- include a non-interpolated dimension whose grid lifecycle method must still
  be invoked;
- finalize the device and host objects independently.

Run these tests with host-only portability and with an actual accelerator
backend where available.

### 5. Test-only owning grid

Add an `OwningTestGrid1D` in the test code. It should implement the grid
interface and own a small dynamically allocated coordinate or coefficient
array. Use it as `DataBox<T, OwningTestGrid1D>` to verify behavior that trivial
production grids cannot expose:

- serialization includes the owned payload;
- deserialization relocates the grid pointer into the serialized buffer;
- `DataBox::setPointer()` reaches every grid slot;
- indexed and named dimensions do not retain stale original pointers;
- `getOnDevice()` deep-copies grid storage;
- interpolation on the device reads the device allocation;
- finalizing a `DataBox` releases both its grid allocations and value
  allocation;
- empty/unmanaged finalization is harmless;
- shallow copies do not become independent owners;
- calculated, written, and consumed offsets remain identical at every nesting
  level.

This test double is the primary proof that all delegation paths have been
threaded correctly.

### 6. Regression matrix

Build and run:

- CPU-only tests;
- Kokkos host tests;
- CUDA or another device backend where available;
- HDF5 disabled;
- HDF5 enabled, including existing regular and piecewise HDF5 round trips.

Run formatting and the repository's normal compile/test checks after the
implementation.

## Suggested implementation order

1. Establish valid empty-grid invariants.
2. Add and test the trivial `RegularGrid1D` API.
3. Add and test recursive `PiecewiseGrid1D` behavior.
4. Thread all operations through every `DataBox` dimension.
5. Add the test-only owning grid and integration coverage.
6. Update serialization documentation and source comments.
7. Run the full build matrix and address ownership or portability failures.

## Acceptance criteria

- Both production grid types expose the complete lifecycle API.
- Every `DataBox` grid slot in `[0, rank)` participates in serialization,
  deserialization, pointer relocation, device copying, and finalization.
- Static grids retain their current interpolation and HDF5 behavior.
- A test-only dynamically allocated grid works end-to-end through `DataBox`.
- Serialization byte counts agree at every layer.
- Device tests demonstrate that no grid retains a host-only pointer.
- The documentation clearly disclaims cross-version and cross-architecture
  compatibility and recommends HDF5 for persistent storage.

Generative AI (OpenAI Codex) was used to assist with modifications to this plan.
