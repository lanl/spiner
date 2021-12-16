Spiner
===

![Build Status](https://github.com/LANL/spiner/actions/workflows/tests.yml/badge.svg)

Performance portable utilities for representing and interpolating
tabulated data. Named for [Brent
Spiner](https://en.wikipedia.org/wiki/Brent_Spiner).

## Performance portability

`Spiner` is compatible with code on CPU, GPU, and everything in between. Read the `ports-of-call` [README](ports-of-call/README.md) for more details.

## Building and Installation

`Spiner` is self-contained. Simply clone it as
```bash
git clone git@gitlab.lanl.gov:jonahm/spiner.git
```
To build and run unit tests,
```bash
mkdir bin
cmake -DBUILD_TESTING=ON ..
make -j
make test
```
and to do convergence testing,
```bash
make convergence
```
after building.

To install,
```bash
make install
```
after configuring and building.

### Build options

- `SPINER_USE_HDF` enables or disables HDF5. Default is `OFF`
- `SPINER_USE_KOKKOS` enables or disables Kokkos. Default is `OFF`.
- `SPINER_USE_CUDA` enables or disables Cuda. Requires Kokkos. Default is `OFF`.
- `BUILD_TESTING` enables or disables tests. Default is `OFF`. If
  this is disabled, then configuration *only* prepares for install and
  provides targets for in-tree builds, as no build step is necessary.
- `SPINER_HDF5_INSTALL_DIR` a hint for cmake about where you may have stashed HDF5.
- `SPINER_KOKKOS_INSTALL_DIR` a hint for cmake about where you may have stashed Kokkos.

## Including spiner in your project

You can build `spiner` in-line with your project, or pre-install
it. It's header-only and the include directories should have the
expected structure. If you build inline, add the following targets to your `cmake`:
```cmake
target_link_libraries(my_project PRIVATE spiner::flags spiner::libs)
```
`spiner::flags` contains compile and include flags, to be included at
compile lines. `spiner::libs` contains linker flags. Since `spiner` is
header-only, `spiner::libs` only contains link flags for dependencies,
such as `HDF5` or `Kokkos`.

For access to both, use
```cmake
target_link_libraries(my_project PRIVATE spiner::spiner)
```

## Dependencies

`Spiner` has no dependencies for the `databox` tool. Simply include it in your project under the `spiner` directory. It is header-only and requires only a few files:

- `spiner/databox.hpp`
- `spiner/interpolation.hpp`
- `spiner/spiner_types.hpp`
- `spiner/sp5.hpp`
- `ports-of-call/portability.hpp`
- `ports-of-call/portable_arrays.hpp`

To use the build system (rather than simply cloning and including the files) requires `cmake`.

The testing tooling requires a few different pieces:

- Unit testing requires [Catch2](https://github.com/catchorg/Catch2),
  which is downloaded automatically if needed.
- Convergence testing requires the scientific python stack, including:
  - python3
  - numpy
  - matplotlib

### HDF5

`Spiner` supports reading and writing DataBox objects into a custom HDF5 format called `SP5`. 
To enable this, compile with the appropriate `HDF5` linking and the flag `-DSPINER_USE_HDF`.
If you use the cmake build system, just configure with `-DSPINER_USE_HDF=ON`.

### CUDA and Kokkos

`Spiner` uses the `ports-of-call` code to optionally support
compilation with CUDA, Kokkos, or none of the above. If `Kokkos` is
discoverable by cmake (for example if you installed it with `spack`),
then the build system should find it automatically. Otherwise you can
specify a location for Kokkos with `SPINER_KOKKOS_INSTALL_DIR`. 

The following spack install was tested with a V100 GPU:
```bash
spack install kokkos-nvcc-wrapper
spack install kokkos~shared+cuda+cuda_lambda+cuda_relocatable_device_code+wrapper cuda_arch=70
```
and then the following cmake configuration line
```C++
cmake -DSPINER_USE_KOKKOS=ON -DSPINER_USE_CUDA=ON -DBUILD_TESTING=ON -DCMAKE_CXX_COMPILER=nvcc_wrapper ..
```
builds the tests for CUDA.

### Clang-Format

Clang-format version 12 is required for committing, and a github
workflow is used to check that code meets format requirements. We
provide a make target in the build system. After configuration, simply
type
```bash
make format
```
to format the code.

Other versions of `clang-format` may work. If you would like to try,
please examine the diff and see if the formatting appears
stable. Otherwise, you may need to upgrade your version of
`clang-format`.

In general, we recommend formatting regularly so that the format calls
do not pollute the diffs. If a format call necessarily pollutes the
diff, do it as a separate commit.

## Features

- Spiner supports interpolation in arbitrary dimensions, and it's fast in 3d and fewer.
- Spiner supports interpolation onto "subtables"

## Interpolation

Interpolation is linear. Here's an example of interpolation in 3D (2D
slice shown). Convergence is second-order, as expected.

![convergence plot](figs/convergence.png)

## Contributors

`Spiner` was primarily developed by Jonah Miller in collaboration with
- Chad Meyer
- Daniel Holladay
- Josh Dolence

The underlying portability machinery in `ports-of-call` was primarily developed by Chad Meyer in collaboration with
- Jonah Miller
- Daniel Holladay
- Josh Dolence

Continuous integration and build system support has been provided by
- Jonah Miller
- Karen Tsai
- Christopher Mauney

## Copyright

© (or copyright) 2019-2021. Triad National Security, LLC. All rights
reserved.  This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S.  Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the
U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting
on its behalf a nonexclusive, paid-up, irrevocable worldwide license
in this material to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.

This program is open source under the BSD-3 License.  Redistribution
and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
