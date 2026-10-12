# Rational stencil optimization candidates

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

Status: recommendations only; no interpolation code changed.

User direction: production rational interpolation will predominantly use
FastNonUniformGrid1D. Regular-grid specialization is excluded. Prioritize a common
spacing API and grid-owned precomputed geometry for nonuniform grids.

## Review plan and evidence

1. Inspect the current stencil, RegularGrid1D, Robust::ratio, and DataBox dispatch.
   Completed.
2. Compare source-level work with the archived throughput baseline. Completed.
3. Rank candidates, identify numerical constraints, and specify experiments.
   Completed below. Speedups have not been measured for these candidates.

The stencil SHA256 matches the archived baseline:
`4793b5a65a905fbce45e2dcedc9d597c57637c659c2b4dd02e92140001f487be`.
Pinned Kokkos Serial results were approximately 6.1–6.2 ns/query in 1D and
22.8–23.2 ns/query in 2D. The approximately fourfold increase is consistent with
four 1D edge evaluations, but does not establish a hardware bottleneck.

On an ordinary interior path with nonzero blend denominator and inactive limiter:

- Each 1D evaluation expresses three slope divisions, two curvature divisions,
  and one blend division: six calls to `ratio`.
- The 2D implementation expresses four such evaluations plus four divisions for
  normalized coordinates and partial derivatives: 28 calls to `ratio`.
- The 2D gathers express 16 sample reads plus four corner reads, but the four
  edges contain only 12 distinct samples. The four corners occur in both an
  x-directed and a y-directed edge, and again in the subtraction term.

These are source-level counts, not assembly instruction counts. Inlining,
common-subexpression elimination, vectorization, and load reuse can reduce them.
Inspect optimized code for the actual backend before attributing time to them.

## First experiments: local algebraic simplifications

These require no grid API changes, caches, or regular-grid specialization and
should precede the larger structural candidates below.

### A. Factor the final 2D value and gradients together

The current Eq. (18) expression repeatedly forms corner weights and differences.
Define two corner differences and a mixed difference, then reuse the gradient
corrections in the value calculation:

```cpp
const T p = f10 - f00;
const T q = f01 - f00;
const T m = (f11 - f01) - p;
const T cx = rip1y - riy - p;
const T cy = rjp1x - rjx - q;
const T tx = cx - qy * m;
const T ty = cy - qx * m;

const T fout = rjx + riy - f00 + qx * cx + qy * ty;
const T dfdx1 = drjdx + qy * (drjp1dx - drjdx) + ratio(tx, dx);
const T dfdx2 = dridy + qx * (drip1dy - dridy) + ratio(ty, dy);
```

This is the same polynomial blend and its derivatives for arbitrary qx/qy,
including extrapolation. It does not change the four edge interpolants or their
limiters. The mixed difference m is the cross coefficient of the bilinear
corner interpolant; tx and ty are its corrections to the edge differences.

Counting just the combined final value/gradient arithmetic, excluding qx/qy
construction, the 1D edge evaluations, loads, and internal `ratio` guards:

| Form | Multiplications | Adds/subtracts | Guarded divisions |
| --- | ---: | ---: | ---: |
| Current, sharing 1-qx and 1-qy | 20 | 23 | 2 |
| Factored above | 6 | 20 | 2 |

These counts are source-level scalar operations, with FMA counted as its two
constituent operations. Compiler reassociation/CSE may change the realized
savings. The full 2D cost includes four 1D evaluations, so the operation-count
reduction is not a claim of a proportional end-to-end speedup.

The expressions are algebraically axis-symmetric, but this particular evaluation
order is not bitwise symmetric under transposition. Floating-point cancellation,
overflow, and rounding can differ; retain scale-aware comparison and transpose
checks. A smaller first patch can use `a + t*(b-a)` for just the four two-term
edge blends, reducing their multiplications from eight to four, before adopting
the joint factoring above.

### B. Reuse the 1D value's inner expression in the derivative

After computing C, the interior helper already forms `S[1] - R*C` for its value.
Name and reuse it:

```cpp
const T L = S[1] - R * C;
fout = f[1] + Q * L;
dfout = L + Q * C + dx[1] * (C - C2) * (1 - w);
```

This removes one add/subtract operation from the combined interior value and
derivative evaluation. It is a small change, not a major 1D speedup hypothesis.
For the lower/upper boundary quadratic, use its local Q and R and the same
pattern: `value = f_left + Q*L`, `gradient = L + Q*C`. Relative to the current
`S + C*(2*x - x_left - x_right)`, this saves three arithmetic operations when Q,
R, and L are reused from the value calculation. It also avoids forming `2*x`,
though it is not a blanket accuracy guarantee.

Optionally expose `deltaC = w*(C1-C2)` when forming `C = C2+deltaC` and use deltaC
instead of recomputing `C-C2` in the derivative. This saves another subtraction
but is a separate numerical change: after rounding, `(C2+deltaC)-C2` need not
return deltaC. Measure error as well as performance.

### C. Let FMA opportunities follow the simpler expressions

The factored blends and `f_left + Q*L` naturally allow fused multiply-add
instructions where the compiler/backend permits contraction. First inspect the
existing optimized code; do not add explicit FMA calls or global fast-math flags
without evidence. These change rounding, and a forced function call or extra
instructions can be counterproductive on some backends. Keep robust denominator
handling and knot/limiter semantics unchanged.

Validation so far: 500 exact rational-arithmetic comparisons of the current and
factored 2D value/gradient expressions passed, with signed values and weights
inside and outside [0,1]. The 1D derivative reuse identity was checked in the same
script. This establishes algebraic equivalence, not floating-point or throughput
results. No stencil implementation was changed.

## 1. Prepare geometry once per axis

Split selection/geometry from interpolation of the four values. A small local
geometry record can contain the window start, enclosing cell, local interval,
three spacings, two adjacent-spacing sums, query offsets, and boundary context.
In 2D, prepare one record for x and one for y, then reuse each for two edges and
for Eq. (18)'s weights/derivatives. Use the already-gathered coordinates for the
enclosing cell instead of calling `grid.x` again.

First preserve division operations and reuse denominator preparation, offsets,
and branch classification. This is the lowest-risk structural experiment and
applies to nonuniform grids as well as regular grids. Some of this work may
already be shared by the compiler.

Next try storing guarded reciprocal spacings and reciprocal spacing sums in
that record, multiplying each numerator by the corresponding reciprocal.
In general 2D interior geometry needs five such reciprocals per axis. Together
with four data-dependent blend divisions, this gives a conceptual 14 divisions
instead of 28, excluding limiter work. It does not improve general-grid 1D
division counts by itself because each denominator is ordinarily used once.

Reciprocal multiplication changes floating-point rounding. Preserve `ratio`'s
small-denominator treatment and compare errors before adopting it. Sharing a
guarded denominator while retaining division is a separate, less invasive change.

## 2. Expose spacing for all grids; cache nonuniform spacing in the grid

Recommended in conjunction with shared per-axis geometry. Give every grid type
portable accessors, provisionally `spacing(i)` and `inverseSpacing(i)`, for the
interval between global point indices i and i+1, with `0 <= i < nPoints()-1`.
These take an already-selected interval index and must not repeat position lookup.
The rational stencil consumes the same API for every grid type, with no
regular-grid-specific interpolation implementation.

- NonUniformGrid1D: precompute N-1 spacings and N-1 inverse spacings when building
  an owning grid, storing them alongside its coordinates.
- FastNonUniformGrid1D: delegate these accessors to its existing
  `coordinates_` NonUniformGrid1D member. Avoid a second copy of the same caches
  in the fast grid; coordinate storage owns them and fast lookup remains separate.
- RegularGrid1D: expose its existing stored spacing and inverse spacing through
  the common interface, without new allocations.
- PiecewiseGrid1D: implement the same global-interval contract, including joins.
  Account explicitly for duplicated join coordinates/zero-width intervals in
  its current indexing; do not blindly return a neighboring subgrid's spacing
  for such an interval or change indexing as part of this optimization.

Start with the two arrays the user suggested. Their added storage is
`2*(N-1)*sizeof(T)` bytes per nonuniform axis (about 16*(N-1) bytes for double).
They are axis geometry, not per-table-cell data, so all rows/columns and repeated
queries can reuse them. Count the extra coordinate-array reads as well as saved
arithmetic: precomputation may increase memory traffic for large scattered grids.

Using cached inverse spacings removes three geometric divisions from each 1D
interior evaluation, leaving two curvature divisions and one blend division.
In 2D, reusing prepared axis geometry also lets the final weights and derivative
terms use the cached reciprocals. Preparing the two adjacent-spacing-sum
reciprocals once per axis would leave a conceptual eight per-query divisions:
four sum reciprocals and four blends, excluding limiter work. These are arithmetic
counts, not predicted speedups or compiler instruction counts.

As a separate measured extension, cache `1/(spacing(i)+spacing(i+1))` for the
N-2 adjacent interval pairs. This targets the remaining curvature denominators
and costs another `(N-2)*sizeof(T)` bytes per axis. Do not substitute
`inverseSpacing(i)+inverseSpacing(i+1)`, which is mathematically different.
With these additional caches, only the data-dependent blend divisions remain
on the ordinary interior path. They cannot be precomputed from grid geometry.

### Ownership and lifecycle requirements

Keep coordinates and their derived arrays under one consistent grid ownership
and memory-space policy. A single allocation partitioned into array spans is a
candidate; multiple allocations are also possible if every lifecycle operation
handles all of them. Preserve only the existing Trivial, AllocatedHost,
AllocatedDevice, Empty, and Unmanaged states.

Audit construction, allocation failure handling, deepCopy, shallowCopy,
getOnDevice, finalize, byte-count accounting, serialization/deserialization,
metadata/device transfers, and HDF loading. Deep copies must copy or rebuild all
arrays in the destination memory space. Shallow aliases borrow all arrays;
finalization frees each owned allocation once. Ordinary copies retain the
existing pointer semantics. Do not add Grid/DataBox RAII copy or move behavior.

An unmanaged coordinate buffer must not silently acquire an independently owned
cache while still reporting Unmanaged. Provide a coherent externally supplied
geometry buffer, or compute spacing on demand when no cache exists. In-place
unmanaged deserialization needs the same explicit decision: serialize cache
arrays with a documented layout/version, or retain the cache-free fallback and
rebuild on an explicit owning copy. Keep coordinate-only HDF data compatible by
rebuilding derived arrays after load. Never read device coordinates on the host
while building caches.

Cache invalidation must be explicit because NonUniformGrid1D currently exposes
mutable `data()` and shallow aliases. Coordinate edits require rebuilding the
spacing arrays, and for FastNonUniformGrid1D also its lookup data. Define whether
edits require reconstructing the grid or an explicit refresh operation; document
alias lifetimes if rebuilding replaces allocations. Changes to tabulated function
values alone do not invalidate geometry caches.

### Numerical contract and validation

Specify whether `inverseSpacing` returns a mathematical reciprocal or a reciprocal
with the same small-denominator guard as the stencil's `ratio`. Do not silently
replace guarded division with an unguarded cached reciprocal. A mathematically
named raw reciprocal accessor can coexist with a guarded stencil calculation,
provided extreme-spacing cases retain an appropriate fallback. Resolve this API
contract before implementation, including overflow/subnormal and zero-width
interval behavior.

Reciprocal multiplication changes rounding even with identical guards. Regular
stored spacing can differ in its last bits from reconstructed coordinate
differences; define and test that behavior for the common interface. Exercise
host/device copies, unmanaged views, repeated finalization, serialization, HDF
loads when enabled, and mutation/rebuild behavior in addition to numerical tests.

## 3. Handle extrapolation before gathering a full stencil

Currently 1D gathers four samples and computes three slopes before discovering
that only one endpoint secant is needed. Strictly outside the table, gather just
the first or last pair and return the linear value and constant derivative.
Similarly, boundary quadratics need three samples and two slopes, not four and
three. Branch placement should be measured against the common interior workload.

There is a stronger 2D simplification. If x is strictly outside its table, its
two horizontal edge interpolants are linear. Their contribution in Eq. (18)
cancels the bilinear subtraction, leaving:

```text
R = (1-qx)*ri(y) + qx*rip1(y)
dR/dx = (rip1(y)-ri(y))/dx
dR/dy = (1-qx)*ri'(y) + qx*rip1'(y)
```

Evaluate only the two vertical interpolants. The y-outside case is symmetric.
When both coordinates are outside, evaluate bilinear extrapolation from four
corner samples. This mainly benefits the boundary/extrapolation distribution,
not the default interior benchmark.

Use strict outside comparisons. At an exact endpoint the current policy returns
an interior derivative; that case must not enter the secant shortcut. Algebraic
cancellation changes floating-point evaluation order, so verify tolerances.

## 4. Reuse the 12 unique 2D values

Gather the four enclosing-cell corners once and use them for both directions
and the Eq. (18) subtraction. Gather the remaining eight edge samples once.
Do not load a full 4x4 block: four of its points are unused by this method.

This can remove redundant addressing/loads if the compiler has not already done
so. Scalar locals or a small fixed-size structure may expose reuse, but retaining
all samples can increase register pressure. On GPUs compare registers, spills,
and throughput; on CPUs check optimized loads. Merely replacing `std::array` with
raw arrays or hand-unrolling its four-element loops is not a demonstrated win.

## 5. Separate the interior evaluation path from boundary logic

Once geometry is prepared, two evaluations on the same axis can share the
interior/boundary/extrapolation decision. Keep the usual interior evaluation
compact and move uncommon cases into separate inline helpers as an experiment.
Retain the data-dependent limiter and blend-degeneracy handling where required.

Do not assume branchless code is faster. Computing all branches can add divisions
and register pressure, and predication can evaluate invalid denominators.
Index-based classification must reproduce coordinate comparisons near knots;
test exact coordinates and their adjacent representable values.

## 6. Optional prepared coefficients for immutable tables

For applications that repeatedly query an unchanged table, precompute secants
and curvature coefficients. This removes substantial per-query work, but adds
memory loads/storage and setup cost. In 2D, axis-specific coefficient storage
can be significant; measure larger, shuffled tables before assuming a benefit.

This should be an explicit prepared-interpolator API, not a hidden DataBox cache.
Shallow aliases can mutate table values without an automatic invalidation signal.
The caller must explicitly rebuild prepared data after changes, select its
memory space, and free it. Preserve DataBox/Grid pointer semantics and avoid
adding RAII copy/move behavior. This is a larger design change than the first
five candidates and should be considered separately.

## Recommended sequence and acceptance criteria

1. First extend the benchmark to FastNonUniformGrid1D with representative
   production coordinate distributions and report whether fast lookup is active
   or has fallen back to binary search. Retain NonUniformGrid1D as a comparison.
   The archived regular-grid timings remain historical context, not the primary
   production baseline. Establish fresh measurements before optimizing.
2. Inspect that Release kernel's divisions, loads, branches, and vectorization.
   Do not infer cache or arithmetic bottlenecks from traffic-control subtraction.
   Once the production workload is established, keep it fixed between candidates.
3. Try the local algebraic simplifications above, starting with the final 2D
   blend and 1D shared intermediates. Check numerical error before timing.
4. Try shared per-axis geometry. Measure denominator reuse and reciprocal
   multiplication separately so numerical changes can be isolated.
5. Add the common spacing API and grid-owned nonuniform caches. Benchmark basic
   spacing caches and the optional adjacent-spacing-sum cache separately.
6. Try extrapolation shortcuts, then explicit corner reuse and branch layout.
7. Consider prepared coefficients only with a concrete repeated-table use case.

For each candidate, run existing Rational tests and independent benchmark
validation. Cover regular/nonuniform grids, signed and mixed curvatures, zero
curvature, the second-interval limiter, exact/near knots, both extrapolation
sides, rectangular tables, tiny spacings, extreme amplitudes, and coordinate
offsets. Numerical rearrangements should use scale-aware tolerances and error
comparisons to the existing implementation/reference, not just aggregate checksums.

Benchmark FastNonUniformGrid1D in 1D/2D across shuffled, ordered, and boundary
workloads, plus small and large tables; report fast-lookup status and cache bytes. Compare baseline and candidate on the same pinned CPU/compiler
settings; require changes larger than observed noise. CPU and GPU improvements
are separate claims. Verify device compilation/execution and register use when
a GPU is available. Keep all temporary data local to the execution space and
perform any coefficient preparation/copies explicitly outside timed evaluation.
