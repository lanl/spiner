# Grid spacing and reciprocal API proposal

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

Status: implemented and validated on CPU; GPU validation unavailable.

## API and numerical contract

Add portable, const, inline accessors to every grid type:

```cpp
T dx(int i) const;                 // interval width
T idx(int i) const;                // ratio(T(1), dx(i))
T invCurvatureWidth(int i) const;  // ratio(T(1), dx(i) + dx(i+1))
```

The first two accept 0 <= i < nPoints()-1; the last accepts
0 <= i < nPoints()-2. These are global interval indices, with no implicit
clamping. The stencil already clamps its window. Keep bounds checks consistent
with project conventions. The reciprocal contract uses PortsOfCall's guarded
ratio, matching the current rational stencil. Preserve addition of adjacent
spacings rather than substituting x(i+2)-x(i), which can round differently.
Nonuniform interval widths are x(i+1)-x(i); regular grids use their stored dx_.
Public accessors follow existing naming conventions; private helpers end in _.

## Implementation by grid

- NonUniformGrid1D owns three additional arrays: N-1 spacings, N-1 inverse
  spacings, and N-2 inverse two-cell spacings. Build eagerly after coordinate
  validation for owned grids. One contiguous allocation containing coordinates
  followed by these arrays simplifies ownership and host/device transfers.
  For N >= 2 this contains 4N-4 values, versus N today.
- FastNonUniformGrid1D forwards accessors to coordinates_. Its lookup structure
  remains independent; do not duplicate the caches in the wrapper.
- RegularGrid1D returns stored scalar quantities, preparing guarded reciprocals
  at construction/loading rather than in accessors. Use the same guarded idx_ scalar for indexing, weights, and the idx() accessor.
  Stored dx_ can differ slightly from differences of rounded x(i) values;
  include this numerical change in validation. No rational-stencil specialization.
- PiecewiseGrid1D delegates intervals and within-subgrid curvature widths to
  subgrids, and precomputes curvature reciprocals for spans crossing joins.
  Resolve join interval indexing explicitly and manage any added storage through
  the same copy, device transfer, serialization, and finalization paths.

Accessors only return stored quantities or delegate to another grid. Do not add
on-demand calculation, missing-cache fallback branches, or optional cache layouts.

## Ownership, mutation, and persistence

Keep the existing allocation states and explicit finalize()/shallow/deep copy
semantics. Cached pointers always refer to the same memory space as coordinates.
Update constructors, HDF loading, deepCopy(), getOnDevice(), shallowCopy(),
finalize(), dynamic memory sizing, dumping, serialization, and setPointer()
together. Audit DataBox and fast-grid consumers of those sizes and offsets.
Never read device coordinates from host code to build a cache.

The pointer constructor remains allocation-free and unmanaged, but requires a
caller-owned buffer large enough for coordinates AND all caches: 4N-4 values
for N >= 2. Use the same contiguous layout for owned and unmanaged grids:

```text
coordinates[N] | dx[N-1] | idx[N-1] | invCurvatureWidth[N-2]
```

Add an explicit boolean constructor option, `initialize_caches`, alongside the
execution-space argument. Require callers to choose it rather than silently
preserving the old coordinates-only buffer contract:

- true: validate coordinates, then populate the cache regions in caller memory.
- false: validate coordinates and existing cache contents without modifying them.

Both paths execute in the supplied execution space, including device buffers.
The caller guarantees sufficient buffer capacity and lifetime; the constructor
cannot infer allocation size from a pointer. Update all pointer-constructor call
sites and document the new buffer contract. Validate caches against the specified
formulas with a documented floating-point comparison tolerance appropriate to T.
Shallow views borrow the entire layout and never free it.

Raw serialization always includes caches; deserialization fixes all pointers
within the borrowed buffer without recalculating or validating cached quantities. This format is
already documented for compatible builds only. HDF continues to store just
coordinates and rebuilds caches on host when loading.

Assume coordinates do not change during use. Leave data() and existing pointer
const correctness unchanged. Do not add mutation enforcement, invalidation,
refresh mechanisms, or alias tracking.

## Follow-up: one regular-grid reciprocal

Removed guarded_idx_. Construction and HDF loading now initialize idx_ with
`ratio(T(1), dx_)`; all reciprocal consumers share it. This intentionally applies
the denominator guard to indexing and linear weights as well as rational geometry.
Validation: six-thread builds passed; Debug and Release each passed 28/28
tests, and the HDF-enabled suite passed 31/31. git diff --check passed.

## Follow-up: comparison helper naming

Renamed cacheEqual_ to approxEqual_ and updated its callers in
NonUniformGrid1D, which supplies cache validation for FastNonUniformGrid1D.
The trailing underscore and comparison behavior are unchanged.

## Follow-up: restore serialized caches directly

Binary serialization includes the entire cache payload. Removed cache validation
from NonUniformGrid1D::setPointer(), so deserialization (including the fast-grid
wrapper) only restores the payload and rebinds storage without recomputing cache
formulas. Unmanaged constructor validation remains explicitly available. Added a
byte-for-byte buffer preservation check to the existing grid round-trip tests.
HDF loading remains a separate coordinate-based file path.

## Follow-up: shared robust comparison tolerance

Use PortsOfCall::Robust::EPS<T>() in approxEqual_ instead of a local epsilon
multiplier. The current dependency defines EPS<T>() as ten machine epsilons;
this tightens the previous sixteen-epsilon relative tolerance. Exact equality
and finite-value checks remain unchanged.
