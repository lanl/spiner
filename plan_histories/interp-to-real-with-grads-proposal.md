# Interpolation with gradients proposal

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

Status: implemented and validated on 2026-10-05.

## Interface

Add plain aggregate result types containing `T value, d_dx1` for 1D and
`T value, d_dx2, d_dx1` for 2D. Add portable, const, noexcept
`interpToRealWithGrads(x1)` and `interpToRealWithGrads(x2, x1)` methods.
Structured bindings follow the existing slowest-to-fastest argument order.
Avoid requiring standard-library tuple operations in device code.
Initially cover fully interpolated rank-1 and rank-2 boxes. An indexed trailing
dimension can be supported separately if desired.

An alternative interface returns T and writes gradients to references:
`T interpToRealWithGrads(T x1, T &d_dx1) const noexcept` and
`T interpToRealWithGrads(T x2, T x1, T &d_dx2, T &d_dx1) const noexcept`.
Choose one interface: mixing these with the aggregate overloads makes a
two-lvalue call ambiguous between 1D reference output and 2D aggregate input.

## Calculation

Use the same stencil selection and value arithmetic as interpToReal. Fetch
each corner once and analytically differentiate the interpolation weights.
For 1D, df/dx1 = (f1-f0)/h1. For 2D, write fij with i indexing x2 and j
indexing x1, and t1, t2 denoting the upper weights:

- df/dx1 = ((1-t2)(f01-f00) + t2(f11-f10))/h1.
- df/dx2 = ((1-t1)(f10-f00) + t1(f11-f01))/h2.

Prefer adding a grid helper `weightsWithGrad(x, ix, w, dw)` returning the
existing weights and their derivatives from a single lookup. Regular grids
use their stored inverse spacing; nonuniform and fast nonuniform grids use
the selected physical coordinate interval. Piecewise grids delegate to the
selected subgrid and adjust the resulting index exactly as weights does.
This preserves the derivative of the actual weight formula and avoids
roundoff from reconstructing regular-grid spacing by subtracting coordinates.
Fast lookup transformations do not alter the physical interpolation variable.

## Semantics and portability

- Preserve current validation and cell selection, including tie breaking.
- At knots return the derivative in the cell selected by weights; a unique
  mathematical derivative need not exist there.
- Preserve existing extrapolation behavior: with normal bounds checks the
  boundary cell is selected but weights extrapolate. Gradients differentiate
  that extrapolation. Do not promise safe out-of-range behavior when bounds
  checks are disabled.
- Derivatives are with respect to supplied grid coordinates. Callers using
  logarithmic coordinates apply any further chain rule themselves.
- Use portable annotations and scalar locals only; no allocations, ownership
  changes, or new DataBox/Grid copy or move semantics.
- Data and grid coordinates must be accessible in the calling execution space.
- Target floating-point T for useful derivatives, without changing existing
  DataBox template constraints as part of this proposal.

## Implementation and verification plan

1. Add aggregate results and grid derivative helpers, then implement the two
   DataBox overloads while leaving existing value-only calls unchanged.
2. Test constant/affine 1D functions and a + b*x1 + c*x2 + d*x1*x2 in 2D,
   using unequal spacings to catch axis swaps. Check returned value agreement.
3. Cover all built-in grid types, knots, piecewise transitions, and supported
   extrapolation. Compare against finite differences strictly inside cells.
4. Exercise supported host/device execution with matching data and grid memory
   spaces. Build with six threads and run the applicable interpolation tests.

## Implementation results

- Added nested aggregate result types and both fully interpolated DataBox
  overloads, preserving coordinate order and value arithmetic.
- Added weightsWithGrad to all four built-in grid types. The regular,
  nonuniform, and fast nonuniform helpers call weights once and compute the
  derivative from the selected interval without a second lookup. Piecewise
  grids delegate to the selected subgrid.
- Preserved the pre-existing const weights_t accessor change in the workspace.
- Added 675 assertions covering all four grid types: constant/affine 1D data,
  bilinear 2D data with a cross term, unequal axis spacings, endpoints,
  extrapolation, changing slopes at knots and piecewise transitions, finite
  differences inside cells, and portable structured bindings on deep device
  copies with explicit cleanup.
- `cmake --build build-debug --parallel 6` succeeded.
- `build-debug/test/test.bin '[Gradients]'` passed all 675 assertions.
- `ctest --test-dir build-debug --output-on-failure -j 6` passed all 21 tests.
- Formatted changed C++ regions with clang-format 19; git diff --check passed.
- The available debug configuration uses Kokkos Serial, with CUDA disabled.
  Actual GPU compilation and execution were not verified.
