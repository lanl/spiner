# Rational interpolation performance baseline proposal

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

Date: 2026-10-10.
Status: approved and implemented on 2026-10-10; original proposal retained below.

## Objective and scope

Measure sustained throughput of many independent calls to rational interpolation
in 1D and 2D. Keep inputs, timing methodology, and output format stable so that
future optimizations can be compared against this baseline.

Add `test/rational_benchmark.cpp`, built as `rational_benchmark.bin` through
`test/CMakeLists.txt` and linked like the existing benchmark targets. Retain the
existing `benchmark.bin`: it measures 3D linear interpolation and includes
analytic function evaluation and error reduction in its timing, so it is not
directly comparable to the proposed measurements.

Benchmark the public API through
`DataBox<Real, RegularGrid1D<Real>, InterpolationStencils::Rational>`:

- 1D: `interpToRealWithGrads(x)`, consuming value and derivative.
- 2D: `interpToRealWithGrads(y, x)`, consuming value and both derivatives.

Use the Linear policy with the same value-and-gradient API, tables, queries,
and timing harness as a control. Do not use `interpToReal` as the rational entry
point: that method does not dispatch through the gradient stencil policy.
The primary metric is rational throughput; the Linear comparison supplies
context, not a substitute for before/after rational measurements.

## Workloads

Start with regular grids to keep lookup cost predictable and focus on the scheme.
Nonuniform grids and direct helper-only benchmarks can be added later as separate
cases. Do not mix their results into the initial baseline.

Proposed default configuration:

| Parameter | Default |
| --- | --- |
| 1D table | 1,024 points |
| 2D table | 256 x 192 points, x fastest |
| Domain | x in [0, 1], y in [0, 2] |
| Queries per pass | 1,048,576 |
| Timed passes per trial | 20 |
| Trials | 7 |
| Warmup passes per case | 3 |
| Input seed | Fixed, recorded in output |
| Policies | Rational and Linear |
| Default distribution | Interior, shuffled |

This produces 20,971,520 interpolation calls per trial, per dimension and policy.
Dimensions are independent workloads: a 2D query count is the total number of
points, not the number of points on each axis. Counts and byte sizes use checked
64-bit arithmetic; each launch must respect the backend's supported index range.

Fill tables once, outside timing, with smooth nonlinear data such as:

```text
f(x)   = 2 + sin(6*pi*x) + 0.2*cos(14*pi*x) + x*x
F(x,y) = 2 + sin(6*pi*x)*cos(2*pi*y) + 0.2*sin(4*pi*x*y) + x*x - 0.3*y*y
```

The 2D data are nonseparable. Both tables have varying signed curvature and avoid
making a constant or linear special case the main performance result. All
trigonometric evaluations happen during setup. Store the same generated values
for both policies, and record an input hash so numerical library differences
across machines are visible.

Generate query arrays on the host with a specified integer PRNG and fixed seed,
then copy them once to the execution space. Use deterministic integer-to-fraction
mapping rather than implementation-dependent standard-library distributions.

Provide three selectable query distributions:

1. **Interior, shuffled (default):** uniformly select eligible cells, then choose
   fractional coordinates in [0.1, 0.9]. Exclude the first and last cells of each
   axis. In 2D, generate x and y independently, not along a diagonal. Shuffle the
   query order reproducibly to measure scattered table access.
2. **Interior, ordered:** use exactly the same interior query multiset, sorted by
   cell (row-major in 2D). This isolates the benefit of spatial locality.
3. **Boundary/extrapolation mix:** in 1D allocate 25% of queries to each of the
   first cell, last cell, below-table region, and above-table region. Outside
   distances range from 0.1 to 1 boundary-cell width. In 2D apply this mix
   independently to both coordinates, yielding all 16 region combinations.
   Report separately because these branches cost less than interior interpolation.

Exact knots remain correctness-test cases rather than a significant fraction of
the throughput workload. Linear extrapolation is the established secant policy.

Allow arbitrary table sizes of at least four points per axis. Suggested optional
size sweeps are 1D: 64 / 1,024 / 65,536, and 2D: 32x24 / 256x192 / 2048x1536.
Report table bytes and query/output bytes; do not label sizes as cache-resident
without reference to the measured machine. Defaults give a compact first run;
large tables and alternate distributions are explicit additional runs.

## Timing and preventing unused work

Use `portableFor` with one interpolation per query. Write every returned component
to preallocated output arrays: two arrays in 1D, three in 2D. This keeps the full
value-and-gradient calculation observable and provides many independent queries
for CPU vectorization and GPU parallelism.

For each case:

1. Allocate and populate the table, queries, and outputs outside timing.
2. Perform three untimed warmup passes and fence.
3. For each timed pass, fence, start `std::chrono::steady_clock`, launch the
   interpolation kernel, fence, and stop the timer.
4. After each pass, outside its timed interval, consume all output components in
   a checksum reduction and accumulate the checksums into observable host state.
   Print the accumulated results after the trials. Consume every pass, not only
   the final overwritten buffer, so repeated work cannot be discarded.
5. Sum the interpolation-kernel durations within each trial and divide by the
   exact number of interpolation calls. Preserve all seven trial measurements.

The checksum reduction runs in the same execution space as the outputs. Transfer
only its small result to the host. Its execution, transfer, and synchronization
are outside the measured interval. Input and output buffers are reused, so this
is a warmed, steady-state workload; the checksum also touches outputs between
passes. Use the same procedure for every optimization comparison.

Timing includes query loads, grid lookup, sample gathering, interpolation,
gradient evaluation, output stores, kernel launch, and completion synchronization.
It excludes allocation, initialization, query generation, data transfers,
correctness checks, checksum reduction, and printing. This is amortized batch
throughput, not single-call latency or isolated arithmetic cost.

Record a companion load/store-only kernel measurement for each dimension to
show the cost of query/output traffic and launch overhead. Do not subtract this
time from interpolation results: costs overlap and the difference is not a
reliable arithmetic-only measurement. Warn when pass durations are too short
for stable timing; increase queries per pass rather than repeating an identical
interpolation inside each thread, where compiler hoisting could distort results.

Alternate policy order between trials to reduce ordering bias. Use Release
builds and fixed backend/thread settings; avoid simultaneous benchmark jobs.

## Correctness and memory spaces

Before timing, compare a deterministic sample of Rational results against a
small host reference based on the paper's Appendix A and Eq. (18), with Spiner's
linear extrapolation policy. Use explicit relative and absolute tolerances
appropriate to `Real`. Check both values and gradients and reject nonfinite
results. Compare Linear with its own reference, not with Rational.

Do not use agreement with the analytic generating function as an exact numerical
oracle: interpolation has approximation error. The existing unit tests remain
the more extensive correctness gate. Checksums are diagnostic and need not be
bitwise identical across parallel reduction orders.

Host setup and validation use host allocations. Timed kernels access only
execution-space-accessible tables, query arrays, and outputs. Copies, fences,
and deallocation are explicit. Preserve DataBox shallow-copy semantics; free
each owning allocation once and finalize owned grids separately where applicable.
Exercise the harness with a small Serial run and with a GPU backend when one is
available. Do not claim GPU results from Kokkos Serial.

## Output and reproducibility

Emit machine-readable CSV with one row per trial plus summary rows. Include:

- Policy, dimension, grid type, nx/ny, distribution, seed, and input hash.
- Queries per pass, passes, trial, total calls, table/query/output bytes.
- Elapsed interpolation time, ns/call, and million calls/second.
- Separate checksums for value and each derivative.
- Median, minimum, and interquartile range across trials in summary rows.
- Matched Linear/Rational throughput ratio, clearly labeled.

Record compiler/version, effective compiler flags, build type, `Real` precision,
backend, Kokkos/Ports-of-Call versions, thread configuration, source commit and
dirty status, and a user-supplied hardware label in accompanying metadata.
Retain raw trial rows so medians and variability can be recomputed.

Example proposed interface:

```sh
cmake --build build-release --target rational_benchmark.bin --parallel 6
./build-release/test/rational_benchmark.bin \
  --dim both --nx1d 1024 --nx2d 256 --ny2d 192 \
  --queries 1048576 --passes 20 --trials 7 --warmup 3 \
  --distribution interior-shuffled --policy both --seed 12345 \
  --output rational-baseline.csv --hardware-label workstation
```

Validate arguments before allocations, including point counts, positive query
and pass counts, and allocation-size overflow. Support small counts for smoke
tests. Performance runs are manual, with no noisy timing thresholds in CTest.

## Implementation and acceptance plan

1. Add the executable and CLI with documented defaults and current copyright/AI
   notices in each changed file.
2. Implement deterministic datasets and queries, memory setup, and untimed
   numerical validation.
3. Add shared timing machinery with compile-time dimension/policy dispatch,
   warmups, pass consumption, the traffic control, and CSV/metadata output.
4. Build with six threads; run the existing tests and tiny benchmark cases for
   both dimensions, policies, and all three distributions. Check count accounting,
   invalid arguments, checksums, finite timings, and explicit cleanup.
5. Run the default Release baseline twice on the same idle machine and compare
   medians and variability. If they are unstable, increase work per pass or
   address the environment before recording a baseline.
6. Save raw measurements and environment metadata with the exact baseline revision.
   For future optimizations, rerun baseline and candidate on the same machine
   with identical options, and require correctness before comparing performance.

Review decision: accept or adjust the default sizes/counts, value-plus-gradient
scope, and the separate executable before implementation.

## Implementation results

Completed the approved plan in `test/rational_benchmark.cpp`, with usage in
`test/rational_benchmark.md` and build metadata generated by
`test/rational_benchmark_metadata.cmake`. All runtime controls, the three query
distributions, both dimensions/policies, the traffic control, reference validation,
per-pass output consumption, trial/summary CSV, and metadata are implemented.
The defaults match the proposal.

Implementation details: the CLI spells the boundary workload `boundary`. Counts
per region differ by at most one when query counts are not divisible by 4 or 16.
The independent host reference differentiates the rational quotient directly.
It validates up to 257 queries per policy; checksums consume all results on every
pass. CTest includes three tiny workload runs and an invalid-query-count test,
without performance thresholds. Compile-command export supplies effective flags
for Makefile/Ninja builds. Source hashes identify the uncommitted baseline exactly.

Validation completed:

- Six-thread Debug and Release builds succeeded.
- Debug and Release test suites: 24/24 passed in each configuration.
- Checked CSV call counts, byte counts, timing units, and row counts.
- Repeated seed/input hashes and checksums matched. Ordered and shuffled query
  multisets gave matching checksums within floating-point tolerance.
- Eleven invalid CLI cases rejected correctly. Four minimum-size cases covered
  each dimension and policy with one query and four-point grids.
- Default Release measurements were repeated; raw trials and metadata are archived
  in `rational-interpolation-benchmark-baseline.md`, including CPU-affinity runs.
- Tested backend: Kokkos Serial, double precision, GNU 15.2.0 on Apple aarch64.
  GPU execution was not available and has not been validated.
