# Ports of call

Contains portability utilities used in Spiner. This functionality is
header-only.

We define a few portability macros which are useful:

1. `PORTABLE_FUNCTION`: decorators necessary for compiling a kernel function
2. `PORTABLE_INLINE_FUNCTION`: ditto, but for when functions ought to be inlined
3. `_WITH_KOKKOS_`: Defined if Kokkos is enabled.
4. `_WITH_CUDA_`: Defined when Cuda is enabled
5. `Real`: a typedef to double (default) or float (if you define
`SINGLE_PRECISION_ENABLED`)
6. `PORTABLE_MALLOC()`, `PORTABLE_FREE()`: A macro or wrapper for
   kokkos_malloc or cudaMalloc, or raw malloc.

At compile time, you define `PORTABILITY_STRATEGY_{KOKKOS,CUDA,NONE}` (if you don't
define it, it defaults to NONE).

There are to be two headers in this folder, for different use cases.

### portability.hpp

- Provides certain macros for decorating functions, lambdas, etc.  Also
provides device malloc/device free operations. This is the principle
path for turning on/off Kokkos. 
- Provides loop abstractions that can be leveraged by the code.  

### portable_array.hpp:

Provides an implementation of a multidimensional array with round
parentheses access. Could be a Kokkos view or something that is not
portable

## Copyright

© (or copyright) 2019. Triad National Security, LLC. All rights
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
