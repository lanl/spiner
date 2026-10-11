# Rational interpolation benchmark

Copyright 2026 Triad National Security, LLC. All rights reserved.
Generative AI was used to assist with this document.

From the repository root, configure a fresh CPU build with CMake 3.24 or newer
and a C++20 compiler, then build with six threads. To match the Kokkos Serial
baseline:

```sh
cmake -S . -B build-rational -DCMAKE_BUILD_TYPE=Release \
  -DSPINER_BUILD_TESTS=ON -DSPINER_TEST_USE_KOKKOS=ON \
  -DSPINER_TEST_USE_KOKKOS_CUDA=OFF \
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF
cmake --build build-rational --target rational_benchmark.bin --parallel 6
./build-rational/test/rational_benchmark.bin --hardware-label my-machine
```

CUDA is enabled by default in Spiner's test configuration, so explicitly disable
it for this CPU build. CMake finds installed dependencies or downloads
Ports-of-Call, Catch2, and Kokkos as needed. For a CPU build without Kokkos, use
`-DSPINER_TEST_USE_KOKKOS=OFF` in a separate fresh directory. This uses the native
Ports-of-Call serial backend and should be recorded as a different configuration
from the Kokkos baseline.

If a Release directory is already configured with those options:

```sh
cmake --build build-release --target rational_benchmark.bin --parallel 6
./build-release/test/rational_benchmark.bin --hardware-label my-machine
```

Use `--help` for options. Defaults run both dimensions and both Rational and
Linear policies on regular grids, with 1,048,576 queries per pass, 20 passes per
trial, seven trials, and three warmup passes. Grid sizes are 1,024 points in 1D
and 256 x 192 in 2D. Each dimension/policy performs 20,971,520 calls per trial.
The `traffic` rows measure a control kernel that reads query coordinates and
writes output arrays, without looking up grid cells, reading the interpolation
table, or performing interpolation. Specifically:

- In 1D, it reads each query's x coordinate and writes x to both output arrays.
- In 2D, it reads each query's x and y coordinates and writes x, x, and y to the
  three output arrays.

These writes use the same output layout as the interpolation kernels, but their
values are coordinate copies, not interpolated values or derivatives. The timing
includes query/output memory traffic, kernel launch, and completion synchronization.
For `traffic` rows, `ns/call` means nanoseconds per processed query; no interpolation
call occurs. Likewise, its checksum columns sum the coordinate copies.

Use `traffic` to gauge the overhead of processing the query and output buffers.
It omits interpolation-table accesses, so it is not a complete model of the memory
cost of interpolation. Do not subtract it from Rational or Linear timings to
estimate arithmetic cost: memory access and computation can overlap.

For a quick correctness and CLI check:

```sh
./build-release/test/rational_benchmark.bin \
  --queries 256 --passes 1 --trials 1 --warmup 0 \
  --nx1d 7 --nx2d 7 --ny2d 4 --distribution boundary \
  --output smoke.csv
```

Runtime options include `--dim 1|2|both`, `--policy rational|linear|both`,
`--nx1d`, `--nx2d`, `--ny2d`, `--queries`, `--passes`, `--trials`, `--warmup`,
`--seed`, `--output`, `--grid regular|fast-nonuniform`, and `--hardware-label`. The policy and dimension inside
kernels remain compile-time specializations. Sizes must fit the library/backend
integer indexing and available memory.

The three `--distribution` choices are:

- `interior-shuffled` (default): uniform eligible cell selection and fractional
  positions in [0.1, 0.9), excluding the first and last interval on each axis.
- `interior-ordered`: the same point multiset, sorted by cell, x fastest in 2D.
- `boundary`: approximately equal counts in the first/last interval and below/
  above the domain, with all 16 combinations in 2D. Extrapolation distances are
  [0.1, 1) boundary-cell widths. Counts differ by at most one per region combination
  when the query count is not divisible by 4 (1D) or 16 (2D).

Queries use SplitMix64 with explicit fraction mapping and bounded rejection
sampling. Tables contain nonlinear functions with mixed curvatures. Generation,
allocation, transfers, and a separate full-table host reference check are untimed.
Host validation samples up to 257 points per policy; checksums inspect all outputs
on every pass. Tables and queries are identical between policies.

Each pass times a single `portableFor` kernel using `steady_clock` and completion
fences. Each query writes its value and all gradients. Immediately afterward, an
untimed device reduction consumes every component and brings back small checksum
results. Timing therefore includes query loads, grid lookup, interpolation,
output stores, launch, and synchronization, but excludes validation and checksum
work. It measures warmed batch throughput, not isolated arithmetic or single-call
latency. Checksums touch output buffers between timed passes. Do not subtract the
traffic-control time from interpolation time.

`rational-baseline.csv` (or `--output`) contains raw trial and summary rows.
`seconds` in a summary is the median trial time; `min_ns_per_call` and
`iqr_ns_per_call` describe variability. Checksum columns are meaningful only in
trial rows. `rational_over_linear_throughput` is populated on Rational summary
rows when both policies run: values below one mean Rational is slower. Other
unused summary/trial fields are zero. Repeated checksum values are expected;
parallel reduction orders need not give identical last bits.

The accompanying `.metadata.txt` records UTC time, build revision, tracked dirty
status, hashes of relevant source files (including untracked benchmark sources),
compiler flags/compile command where available, dependency/backend information,
precision, seed/options, and hardware label. Kokkos configuration includes backend
concurrency. Input hashes cover native table/query bytes, so ordered and shuffled
runs have different hashes despite sharing a query multiset. Hashes are diagnostic
across architectures, not a portable serialization format.

Run baseline and candidate on the same otherwise idle machine, with identical
compiler/backend settings and options. Keep raw CSV and metadata together. The
short-pass warning suggests increasing queries to amortize timing/launch noise.
Smoke tests are registered in CTest with no performance thresholds. GPU results
require an actual GPU backend; Kokkos Serial results are CPU measurements.

## Fast nonuniform grids

Use `--grid fast-nonuniform` for the production-oriented spacing-cache benchmark.
The regular grid remains the default for comparison with the original baseline.
Nonuniform coordinates are `extent * (t + t*t)/2`, where `t=i/(N-1)` and extent
is 1 for x and 2 for y. Interval widths vary smoothly by approximately a factor
of three. Fast lookup is required; the benchmark does not silently fall back to
binary search. Query fractions and extrapolation distances use local cell widths.
Grid type is selected outside timing, with separate compiled kernel instantiations.

Each policy owns a separate device copy of its grids, explicitly freed along
with the table. Host construction, cache preparation, device copies, and freeing
are outside timing. CSV grid labels and source metadata identify this workload.
Compare ordered/shuffled queries for small and large tables with identical
compiler and execution settings; saved divisions alone do not establish a speedup.

```sh
./build-release/test/rational_benchmark.bin --grid fast-nonuniform \
  --distribution interior-shuffled --output rational-fast.csv
```
