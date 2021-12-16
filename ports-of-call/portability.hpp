#ifndef _PORTABILITY_HPP_
#define _PORTABILITY_HPP_

// ========================================================================================
// © (or copyright) 2019-2021. Triad National Security, LLC. All rights
// reserved.  This program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
// operated by Triad National Security, LLC for the U.S.  Department of
// Energy/National Nuclear Security Administration. All rights in the
// program are reserved by Triad National Security, LLC, and the
// U.S. Department of Energy/National Nuclear Security
// Administration. The Government is granted for itself and others acting
// on its behalf a nonexclusive, paid-up, irrevocable worldwide license
// in this material to reproduce, prepare derivative works, distribute
// copies to the public, perform publicly and display publicly, and to
// permit others to do so.
// ========================================================================================

#include <string>

#ifdef PORTABILITY_STRATEGY_KOKKOS
#include "Kokkos_Core.hpp"
#define PORTABLE_FUNCTION KOKKOS_FUNCTION
#define PORTABLE_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#define PORTABLE_FORCEINLINE_FUNCTION KOKKOS_FORCEINLINE_FUNCTION
#define PORTABLE_LAMBDA KOKKOS_LAMBDA
#define _WITH_KOKKOS_
// The following is a malloc on default memory space
#define PORTABLE_MALLOC(size) Kokkos::kokkos_malloc<>(size)
#define PORTABLE_FREE(ptr) Kokkos::kokkos_free<>(ptr)
// Do we want to include additional terms here (for memory spaces, etc.)?
#else
#ifdef PORTABILITY_STRATEGY_CUDA
#include "cuda.h"
#define PORTABLE_FUNCTION __host__ __device__
#define PORTABLE_INLINE_FUNCTION __host__ __device__ inline
#define PORTABLE_FORCEINLINE_FUNCTION                                          \
  __host__ __device__ inline __attribute__((always_inline))
#define PORTABLE_LAMBDA [=] __host__ __device__
void *PORTABLE_MALLOC(size_t size) {
  void *devPtr cudaError_t e = cudaMalloc(devPtr, size);
  return devPtr;
}
void PORTABLE_FREE(void *ptr) { cudaError_t e = cudaFree(ptr); }
#define _WITH_CUDA_
// It is worth noting here that we will not define
// _WITH_CUDA_ when we are doing KOKKOS (even with the
// CUDA backend)  Rely on KOKKOS_HAVE_CUDA in that case
#else
#define PORTABLE_FUNCTION
#define PORTABLE_INLINE_FUNCTION inline
#define PORTABLE_FORCEINLINE_FUNCTION inline __attribute__((always_inline))
#define PORTABLE_LAMBDA [=]
#define PORTABLE_MALLOC(size) malloc(size)
#define PORTABLE_FREE(ptr) free(ptr)
#endif
#endif

#ifndef SINGLE_PRECISION_ENABLED
#define SINGLE_PRECISION_ENABLED 0
#endif

#if SINGLE_PRECISION_ENABLED
typedef float Real;
#else
typedef double Real;
#endif

template <typename Function>
void portableFor(const char *name, int start, int stop, Function function) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using policy = Kokkos::RangePolicy<>;
  Kokkos::parallel_for(name, policy(start, stop), function);
#else
  for (int i = start; i < stop; i++) {
    function(i);
  }
#endif
}

template <typename Function>
void portableFor(const char *name, int starty, int stopy, int startx, int stopx,
                 Function function) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
  Kokkos::parallel_for(name, Policy2D({starty, startx}, {stopy, stopx}),
                       function);
#else
  for (int iy = starty; iy < stopy; iy++) {
    for (int ix = startx; ix < stopx; ix++) {
      function(iy, ix);
    }
  }
#endif
}

template <typename Function>
void portableFor(const char *name, int startz, int stopz, int starty, int stopy,
                 int startx, int stopx, Function function) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
  Kokkos::parallel_for(
      name, Policy3D({startz, starty, startx}, {stopz, stopy, stopx}),
      function);
#else
  for (int iz = startz; iz < stopz; iz++) {
    for (int iy = starty; iy < stopy; iy++) {
      for (int ix = startx; ix < stopx; ix++) {
        function(iz, iy, ix);
      }
    }
  }
#endif
}

template <typename Function>
void portableFor(const char *name, int starta, int stopa, int startz, int stopz,
                 int starty, int stopy, int startx, int stopx,
                 Function function) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy4D = Kokkos::MDRangePolicy<Kokkos::Rank<4>>;
  Kokkos::parallel_for(
      name,
      Policy4D({starta, startz, starty, startx}, {stopa, stopz, stopy, stopx}),
      function);
#else
  for (int ia = starta; ia < stopa; ia++) {
    for (int iz = startz; iz < stopz; iz++) {
      for (int iy = starty; iy < stopy; iy++) {
        for (int ix = startx; ix < stopx; ix++) {
          function(ia, iz, iy, ix);
        }
      }
    }
  }
#endif
}

template <typename Function>
void portableFor(const char *name, int startb, int stopb, int starta, int stopa,
                 int startz, int stopz, int starty, int stopy, int startx,
                 int stopx, Function function) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy5D = Kokkos::MDRangePolicy<Kokkos::Rank<5>>;
  Kokkos::parallel_for(name,
                       Policy5D({startb, starta, startz, starty, startx},
                                {stopb, stopa, stopz, stopy, stopx}),
                       function);
#else
  for (int ib = startb; ib < stopb; ib++) {
    for (int ia = starta; ia < stopa; ia++) {
      for (int iz = startz; iz < stopz; iz++) {
        for (int iy = starty; iy < stopy; iy++) {
          for (int ix = startx; ix < stopx; ix++) {
            function(ib, ia, iz, iy, ix);
          }
        }
      }
    }
  }
#endif
}

template <typename Function, typename T>
void portableReduce(const char *name, int startz, int stopz, int starty,
                    int stopy, int startx, int stopx, Function function,
                    T &reduced) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy3D = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
  Kokkos::parallel_reduce(
      name, Policy3D({startz, starty, startx}, {stopz, stopy, stopx}), function,
      reduced);
#else
  for (int iz = startz; iz < stopz; iz++) {
    for (int iy = starty; iy < stopy; iy++) {
      for (int ix = startx; ix < stopx; ix++) {
        function(iz, iy, ix, reduced);
      }
    }
  }
#endif
}

template <typename Function, typename T>
void portableReduce(const char *name, int starta, int stopa, int startz,
                    int stopz, int starty, int stopy, int startx, int stopx,
                    Function function, T &reduced) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy4D = Kokkos::MDRangePolicy<Kokkos::Rank<4>>;
  Kokkos::parallel_reduce(
      name,
      Policy4D({starta, startz, starty, startx}, {stopa, stopz, stopy, stopx}),
      function, reduced);
#else
  for (int ia = starta; ia < stopa; ia++) {
    for (int iz = startz; iz < stopz; iz++) {
      for (int iy = starty; iy < stopy; iy++) {
        for (int ix = startx; ix < stopx; ix++) {
          function(ia, iz, iy, ix, reduced);
        }
      }
    }
  }
#endif
}

template <typename Function, typename T>
void portableReduce(const char *name, int startb, int stopb, int starta,
                    int stopa, int startz, int stopz, int starty, int stopy,
                    int startx, int stopx, Function function, T &reduced) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  using Policy5D = Kokkos::MDRangePolicy<Kokkos::Rank<5>>;
  Kokkos::parallel_reduce(name,
                          Policy5D({startb, starta, startz, starty, startx},
                                   {stopb, stopa, stopz, stopy, stopx}),
                          function, reduced);
#else
  for (int ib = startb; ib < stopb; ib++) {
    for (int ia = starta; ia < stopa; ia++) {
      for (int iz = startz; iz < stopz; iz++) {
        for (int iy = starty; iy < stopy; iy++) {
          for (int ix = startx; ix < stopx; ix++) {
            function(ib, ia, iz, iy, ix, reduced);
          }
        }
      }
    }
  }
#endif
}

#endif // PORTABILITY_HPP
