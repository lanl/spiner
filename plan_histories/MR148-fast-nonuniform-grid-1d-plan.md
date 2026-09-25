# Fast NQT-Accelerated Nonuniform Grid

Generative AI was used to assist with this document.

## Summary

Add `FastNonUniformGrid1D<double>`, a distinct grid type that owns an
authoritative `NonUniformGrid1D<double>` and optionally augments it with an
exact O(1) lookup table uniformly spaced in corrected NQT-asinh coordinates.

The accelerator is selected at construction, HDF5 load, or explicit host-side
reconfiguration. When it cannot satisfy the configured memory limit, automatic
mode uses the underlying grid's binary search.

## Public API and Configuration

- Add `FastNonUniformGridPolicy` with:
  - `Automatic`: build the LUT when it fits; otherwise use binary search.
  - `RequireFast`: require the LUT to fit and be valid, otherwise terminate
    with a clear requirement failure.
  - `ForceBinary`: omit or discard the LUT.
- Provide vector and initializer-list constructors:
  - `(points, scale)` uses `Automatic` and an 8x LUT-entry cap.
  - `(points, scale, policy, max_lookup_ratio)` provides full control.
- Add `reconfigureLookup(policy, max_lookup_ratio)`:
  - Operates only on a host-owned grid.
  - Builds a missing LUT, retains a compatible existing LUT, or discards it
    according to the new policy and cap.
  - `ForceBinary` immediately frees the LUT while retaining the underlying
    coordinates.
  - Requires a positive cap even when selecting `ForceBinary`, so the setting
    remains usable in later reconfiguration.
  - Must be called only on the canonical owner when no shallow alias depends on
    its LUT, matching existing explicit-lifecycle conventions.
- Require a finite positive scale. Transform coordinates and queries as
  `NQT::O1::<variant>::asinh(x / scale)`.
- Constrain the grid template to `double`; single-precision builds cannot
  instantiate this type.
- Preserve the standard grid API: `x`, `index`, `weights`, `min`, `max`,
  `nPoints`, `isWellFormed`, `dataStatus`, and lifecycle/serialization
  operations.
- Provide read-only inspection through `const data()`, `scale()`,
  `lookupSize()`, `maxLookupRatio()`, `requestedPolicy()`, and
  `usesFastLookup()`. Do not expose mutable coordinates or the mutable wrapped
  grid.
- Add `SPINER_USE_PORTABLE_NQT`, defaulting to `OFF`. `OFF` selects
  `O1::Aliased::asinh`; `ON` selects `O1::Portable::asinh`. Export the selection
  through the header-only CMake target.
- Update the Ports-of-Call dependency from `v3.1.0` to corrected release
  `v3.1.1`.

## Implementation

- Store the underlying `NonUniformGrid1D<double>`, an inline
  `RegularGrid1D<double>` describing transformed cells, an optional owned `int`
  LUT, construction settings, and LUT allocation status.
- Build or rebuild the accelerator on the host:
  1. Transform the original coordinates and determine their minimum positive
     transformed spacing.
  2. Choose a uniform cell width strictly no larger than that spacing, using a
     conservative `nextafter` adjustment.
  3. Create enough cells to cover the transformed coordinate range.
  4. Verify that every cell crosses at most one original boundary, enlarging
     the table within the cap if floating-point binning requires it.
  5. Store the original interval index at each cell's lower edge without
     retaining the transformed coordinate array.
- Limit the LUT to `max_lookup_ratio * nPoints()` entries with overflow-safe
  size calculations. Invalid, collapsed, non-finite, or oversized transformed
  layouts trigger binary fallback in `Automatic` and failure in `RequireFast`.
- Reconfiguration behavior:
  - Reuse the current LUT when it remains valid under the new settings.
  - Discard it when selecting `ForceBinary` or when a reduced automatic cap no
    longer permits it.
  - Attempt construction when moving from binary fallback to `Automatic` or
    `RequireFast`, or when a larger cap now permits acceleration.
  - Reject reconfiguration of device-resident, unmanaged, or empty grids.
- Implement accelerated lookup by:
  1. Preserving existing endpoint/extrapolation clamps.
  2. Transforming the query and finding its uniform auxiliary cell.
  3. Reading one candidate interval from the LUT.
  4. Performing at most one comparison with the next physical coordinate to
     obtain the exact interval.
- Compute interpolation weights from the accelerated index and authoritative
  physical coordinates. Do not delegate fast-mode weights to the wrapped
  grid's binary-search implementation.
- Follow existing pointer semantics: ordinary copies are shallow and require
  one explicit finalization; do not add RAII copy/move behavior.
- Deep-copy both coordinate and LUT allocations in `copy()` and
  `getOnDevice()`. Extend `finalize`, allocation status, dynamic-memory sizing,
  binary serialization, relocation, and deserialization to cover both arrays
  and both active modes.
- Serialize binary data as the inline wrapper followed by underlying
  coordinates and optional LUT indices.
- For HDF5, store the authoritative coordinates plus scale, ratio, and
  requested policy. Rebuild the derived LUT on load rather than persisting it,
  then allow callers to change or discard it with `reconfigureLookup()`.
- Export the new header through `interpolation.hpp`, document ownership,
  reconfiguration, memory limits, transform, fallback, and build configuration,
  and add a changelog entry.

## Tests and Benchmark

- Compare `index`, `weights`, and `DataBox` interpolation against
  `NonUniformGrid1D` at every coordinate, immediately around boundaries,
  throughout cells, across zero, in both linear and logarithmic regions, and
  outside the domain.
- Exercise successful automatic acceleration, automatic fallback from
  pathological spacing, required acceleration, and forced binary lookup.
- Test reconfiguration from fast to binary, binary to fast after increasing
  the cap, fast to automatic fallback after reducing the cap, and compatible
  settings that retain the existing LUT.
- Verify that reconfiguration rejects non-host-owned states and document/test
  the shallow-alias ownership precondition where practical.
- Test default/custom scales and ratios, read-only metadata, shallow-copy
  conventions, explicit deep copy, idempotent owned cleanup, binary
  serialization/relocation, HDF5 reconstruction and post-load LUT disposal, and
  multi-axis `DataBox` integration.
- Run fast and fallback grids in portable device kernels and verify host/device
  results and allocation states.
- Configure and test both aliased-default and portable-NQT builds.
- Add a repeatable nonuniform-grid benchmark that generates signed-log-spaced
  coordinates, compares binary and accelerated lookup/interpolation on
  host/device backends, fences timing correctly, and reports time per query plus
  result agreement. Do not impose a machine-dependent speed threshold.

## Assumptions

- Ports-of-Call `v3.1.1` makes both O1 `asinh` variants monotone while
  preserving their existing namespace and signatures.
- Fast-grid construction is owning-only; unmanaged coordinate construction is
  intentionally excluded.
- Full reconfiguration is host-owned only and is performed before creating
  device copies.
- Fallback mode is fixed between explicit reconfiguration events.
- Accelerated and fallback modes exactly reproduce `NonUniformGrid1D` interval,
  interpolation, and extrapolation semantics.
