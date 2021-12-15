#ifndef _PORTABLE_ARRAYS_HPP_
#define _PORTABLE_ARRAYS_HPP_
//========================================================================================
// Portable MDAarray, based on AthenaArray in Athena++
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//
// © (or copyright) 2019-2021. Triad National Security, LLC. All rights
// reserved.  This program was produced under U.S. Government contract
// 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
// is operated by Triad National Security, LLC for the U.S.
// Department of Energy/National Nuclear Security Administration. All
// rights in the program are reserved by Triad National Security, LLC,
// and the U.S. Department of Energy/National Nuclear Security
// Administration. The Government is granted for itself and others
// acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works,
// distribute copies to the public, perform publicly and display
// publicly, and to permit others to do so.
// ========================================================================================

//  The operator() is overloaded, e.g. elements of a 4D array of size
//  [N4xN3xN2xN1] are accessed as:  A(n,k,j,i) = A[i + N1*(j + N2*(k + N3*n))]
//  NOTE THE TRAILING INDEX INSIDE THE PARENTHESES IS INDEXED FASTEST

#include "portability.hpp"
#include <algorithm>
#include <assert.h>
#include <cstddef> // size_t
#include <cstring> // memset()
#include <utility> // swap()

template <typename T>
class PortableMDArray {
 public:
  static constexpr int MAXDIM = 6;

  // ctors
  // default ctor: simply set null PortableMDArray
  PORTABLE_FUNCTION
  PortableMDArray() noexcept
      : pdata_(nullptr), nx1_(0), nx2_(0), nx3_(0), nx4_(0), nx5_(0), nx6_(0) {}
  PORTABLE_FUNCTION PortableMDArray(T *data, int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(1), nx3_(1), nx4_(1), nx5_(1), nx6_(1) {}
  PORTABLE_FUNCTION
  PortableMDArray(T *data, int nx2, int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(nx2), nx3_(1), nx4_(1), nx5_(1), nx6_(1) {
  }
  PORTABLE_FUNCTION
  PortableMDArray(T *data, int nx3, int nx2, int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(1), nx5_(1),
        nx6_(1) {}
  PORTABLE_FUNCTION
  PortableMDArray(T *data, int nx4, int nx3, int nx2, int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(1),
        nx6_(1) {}
  PORTABLE_FUNCTION
  PortableMDArray(T *data, int nx5, int nx4, int nx3, int nx2, int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(nx5),
        nx6_(1) {}
  PORTABLE_FUNCTION
  PortableMDArray(T *data, int nx6, int nx5, int nx4, int nx3, int nx2,
                  int nx1) noexcept
      : pdata_(data), nx1_(nx1), nx2_(nx2), nx3_(nx3), nx4_(nx4), nx5_(nx5),
        nx6_(nx6) {}

  // define copy constructor and overload assignment operator so both do deep
  // copies.
  PortableMDArray(const PortableMDArray<T> &t) noexcept;
  PortableMDArray<T> &operator=(const PortableMDArray<T> &t) noexcept;

  // public functions to allocate/deallocate memory for 1D-5D data
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx1) noexcept;
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx2, int nx1) noexcept;
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx3, int nx2,
                                            int nx1) noexcept;
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx4, int nx3, int nx2,
                                            int nx1) noexcept;
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx5, int nx4, int nx3,
                                            int nx2, int nx1) noexcept;
  PORTABLE_FUNCTION void NewPortableMDArray(T *data, int nx6, int nx5, int nx4,
                                            int nx3, int nx2, int nx1) noexcept;

  // public function to swap underlying data pointers of two equally-sized
  // arrays
  void SwapPortableMDArray(PortableMDArray<T> &array2);

  // functions to get array dimensions
  PORTABLE_FORCEINLINE_FUNCTION int GetDim1() const { return nx1_; }
  PORTABLE_FORCEINLINE_FUNCTION int GetDim2() const { return nx2_; }
  PORTABLE_FORCEINLINE_FUNCTION int GetDim3() const { return nx3_; }
  PORTABLE_FORCEINLINE_FUNCTION int GetDim4() const { return nx4_; }
  PORTABLE_FORCEINLINE_FUNCTION int GetDim5() const { return nx5_; }
  PORTABLE_FORCEINLINE_FUNCTION int GetDim6() const { return nx6_; }
  PORTABLE_INLINE_FUNCTION int GetDim(size_t i) const {
    // TODO: remove if performance cirtical
    assert(0 < i && i <= 6 && "PortableMDArrays are max 6D");
    switch (i) {
    case 1:
      return GetDim1();
    case 2:
      return GetDim2();
    case 3:
      return GetDim3();
    case 4:
      return GetDim4();
    case 5:
      return GetDim5();
    case 6:
      return GetDim6();
    }
    return -1;
  }

  // a function to get the total size of the array
  PORTABLE_FORCEINLINE_FUNCTION int GetSize() const {
    return nx1_ * nx2_ * nx3_ * nx4_ * nx5_ * nx6_;
  }
  PORTABLE_FORCEINLINE_FUNCTION std::size_t GetSizeInBytes() const {
    return nx1_ * nx2_ * nx3_ * nx4_ * nx5_ * nx6_ * sizeof(T);
  }

  PORTABLE_INLINE_FUNCTION size_t GetRank() const {
    for (int i = 6; i >= 1; i--) {
      if (GetDim(i) > 1) return i;
    }
    return 0;
  }

  PORTABLE_INLINE_FUNCTION void Reshape(int nx6, int nx5, int nx4, int nx3,
                                        int nx2, int nx1) {
    assert(nx6 * nx5 * nx4 * nx3 * nx2 * nx1 == GetSize());
    nx1_ = nx1;
    nx2_ = nx2;
    nx3_ = nx3;
    nx4_ = nx4;
    nx5_ = nx5;
    nx6_ = nx6;
  }
  PORTABLE_INLINE_FUNCTION void Reshape(int nx5, int nx4, int nx3, int nx2,
                                        int nx1) {
    Reshape(1, nx5, nx4, nx3, nx2, nx1);
  }
  PORTABLE_INLINE_FUNCTION void Reshape(int nx4, int nx3, int nx2, int nx1) {
    Reshape(1, 1, nx4, nx3, nx2, nx1);
  }
  PORTABLE_INLINE_FUNCTION void Reshape(int nx3, int nx2, int nx1) {
    Reshape(1, 1, 1, nx3, nx2, nx1);
  }
  PORTABLE_INLINE_FUNCTION void Reshape(int nx2, int nx1) {
    Reshape(1, 1, 1, 1, nx2, nx1);
  }
  PORTABLE_INLINE_FUNCTION void Reshape(int nx1) {
    Reshape(1, 1, 1, 1, 1, nx1);
  }

  PORTABLE_FORCEINLINE_FUNCTION bool IsShallowSlice() { return true; }
  PORTABLE_FORCEINLINE_FUNCTION bool IsEmpty() { return GetSize() < 1; }
  // "getter" function to access private data member
  // TODO(felker): Replace this unrestricted "getter" with a limited, safer
  // alternative.
  // TODO(felker): Rename function. Conflicts with "PortableMDArray<> data"
  // OutputData member.
  PORTABLE_FORCEINLINE_FUNCTION T *data() { return pdata_; }
  PORTABLE_FORCEINLINE_FUNCTION const T *data() const { return pdata_; }
  PORTABLE_FORCEINLINE_FUNCTION T *begin() { return pdata_; }
  PORTABLE_FORCEINLINE_FUNCTION T *end() { return pdata_ + GetSize(); }

  // overload "function call" operator() to access 1d-5d data
  // provides Fortran-like syntax for multidimensional arrays vs. "subscript"
  // operator[]

  // "non-const variants" called for "PortableMDArray<T>()" provide read/write
  // access via returning by reference, enabling assignment on returned l-value,
  // e.g.: a(3) = 3.0;
  PORTABLE_FORCEINLINE_FUNCTION T &operator()() { return pdata_[0]; }

  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n) { return pdata_[n]; }
  // "const variants" called for "const PortableMDArray<T>" returns T by value,
  // since T is typically a built-in type (versus "const T &" to avoid copying
  // for general types)
  PORTABLE_FORCEINLINE_FUNCTION T &operator()() const { return pdata_[0]; }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n) const {
    return pdata_[n];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int i) {
    return pdata_[i + nx1_ * n];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int i) const {
    return pdata_[i + nx1_ * n];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int j,
                                              const int i) {
    return pdata_[i + nx1_ * (j + nx2_ * n)];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int j,
                                              const int i) const {
    return pdata_[i + nx1_ * (j + nx2_ * n)];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int k,
                                              const int j, const int i) {
    return pdata_[i + nx1_ * (j + nx2_ * (k + nx3_ * n))];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int n, const int k,
                                              const int j, const int i) const {
    return pdata_[i + nx1_ * (j + nx2_ * (k + nx3_ * n))];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &
  operator()(const int m, const int n, const int k, const int j, const int i) {
    return pdata_[i + nx1_ * (j + nx2_ * (k + nx3_ * (n + nx4_ * m)))];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int m, const int n,
                                              const int k, const int j,
                                              const int i) const {
    return pdata_[i + nx1_ * (j + nx2_ * (k + nx3_ * (n + nx4_ * m)))];
  }
  // int l?, int o?
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int p, const int m,
                                              const int n, const int k,
                                              const int j, const int i) {
    return pdata_[i +
                  nx1_ * (j + nx2_ * (k + nx3_ * (n + nx4_ * (m + nx5_ * p))))];
  }
  PORTABLE_FORCEINLINE_FUNCTION T &operator()(const int p, const int m,
                                              const int n, const int k,
                                              const int j, const int i) const {
    return pdata_[i +
                  nx1_ * (j + nx2_ * (k + nx3_ * (n + nx4_ * (m + nx5_ * p))))];
  }

  PortableMDArray<T> &operator*=(T scale) {
    std::transform(pdata_, pdata_ + GetSize(), pdata_,
                   [scale](T val) { return scale * val; });
    return *this;
  }

  PortableMDArray<T> &operator+=(const PortableMDArray<T> &other) {
    assert(GetSize() == other.GetSize());
    std::transform(pdata_, pdata_ + GetSize(), other.pdata_, pdata_,
                   std::plus<T>());
    return *this;
  }

  PortableMDArray<T> &operator-=(const PortableMDArray<T> &other) {
    assert(GetSize() == other.GetSize());
    std::transform(pdata_, pdata_ + GetSize(), other.pdata_, pdata_,
                   std::minus<T>());
    return *this;
  }

  // Checks that arrays point to same data with same shape
  // note this POINTER equivalence, not data equivalence
  bool operator==(const PortableMDArray<T> &other) const;
  bool operator!=(const PortableMDArray<T> &other) const {
    return !(*this == other);
  }

  // (deferred) initialize an array with slice from another array
  PORTABLE_FUNCTION
  void InitWithShallowSlice(const PortableMDArray<T> &src, const int dim,
                            const int indx, const int nvar);

 private:
  T *pdata_;
  int nx1_, nx2_, nx3_, nx4_, nx5_, nx6_;
};

// copy constructor (does a shallow copy)

template <typename T>
PortableMDArray<T>::PortableMDArray(const PortableMDArray<T> &src) noexcept {
  nx1_ = src.nx1_;
  nx2_ = src.nx2_;
  nx3_ = src.nx3_;
  nx4_ = src.nx4_;
  nx5_ = src.nx5_;
  nx6_ = src.nx6_;
  if (src.pdata_) pdata_ = src.pdata_;
}

// shallow copy assignment operator

template <typename T>
PortableMDArray<T> &
PortableMDArray<T>::operator=(const PortableMDArray<T> &src) noexcept {
  if (this != &src) {
    nx1_ = src.nx1_;
    nx2_ = src.nx2_;
    nx3_ = src.nx3_;
    nx4_ = src.nx4_;
    nx5_ = src.nx5_;
    nx6_ = src.nx6_;
    pdata_ = src.pdata_;
  }
  return *this;
}

// Checks that arrays point to same data with same shape
// note this POINTER equivalence, not data equivalence
template <typename T>
bool PortableMDArray<T>::operator==(const PortableMDArray<T> &rhs) const {
  return (pdata_ == rhs.pdata_ && nx1_ == rhs.nx1_ && nx2_ == rhs.nx2_ &&
          nx3_ == rhs.nx3_ && nx4_ == rhs.nx4_ && nx5_ == rhs.nx5_ &&
          nx6_ == rhs.nx6_);
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::InitWithShallowSlice()
//  \brief shallow copy of nvar elements in dimension dim of an array, starting
//  at index=indx. Copies pointer to data, but not data itself.

//  Shallow slice is only able to address the "nvar" range in "dim", and all
//  entries of the src array for d<dim (cannot access any nx4=2, etc. entries if
//  dim=3 for example)

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::InitWithShallowSlice(const PortableMDArray<T> &src,
                                         const int dim, const int indx,
                                         const int nvar) {
  pdata_ = src.pdata_;
  if (dim == 6) {
    nx6_ = nvar;
    nx5_ = src.nx5_;
    nx4_ = src.nx4_;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx * (nx1_ * nx2_ * nx3_ * nx4_ * nx5_);
  } else if (dim == 5) {
    nx6_ = 1;
    nx5_ = nvar;
    nx4_ = src.nx4_;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx * (nx1_ * nx2_ * nx3_ * nx4_);
  } else if (dim == 4) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = nvar;
    nx3_ = src.nx3_;
    nx2_ = src.nx2_;
    nx1_ = src.nx1_;
    pdata_ += indx * (nx1_ * nx2_ * nx3_);
  } else if (dim == 3) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = nvar;     // nx3
    nx2_ = src.nx2_; // nx2
    nx1_ = src.nx1_; // nx1
    pdata_ += indx * (nx1_ * nx2_);
  } else if (dim == 2) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = 1;
    nx2_ = nvar;
    nx1_ = src.nx1_;
    pdata_ += indx * (nx1_);
  } else if (dim == 1) {
    nx6_ = 1;
    nx5_ = 1;
    nx4_ = 1;
    nx3_ = 1;
    nx2_ = 1;
    nx1_ = nvar;
    pdata_ += indx;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief allocate new 1D array with elements initialized to zero.

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = 1;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief 2d data allocation

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx2, int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = 1;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief 3d data allocation

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx3, int nx2,
                                       int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = 1;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief 4d data allocation

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx4, int nx3, int nx2,
                                       int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = 1;
  nx6_ = 1;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief 5d data allocation

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx5, int nx4, int nx3,
                                       int nx2, int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = nx5;
  nx6_ = 1;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::NewPortableMDArray()
//  \brief 6d data allocation

template <typename T>
PORTABLE_FUNCTION void
PortableMDArray<T>::NewPortableMDArray(T *data, int nx6, int nx5, int nx4,
                                       int nx3, int nx2, int nx1) noexcept {
  nx1_ = nx1;
  nx2_ = nx2;
  nx3_ = nx3;
  nx4_ = nx4;
  nx5_ = nx5;
  nx6_ = nx6;
  pdata_ = data;
}

//----------------------------------------------------------------------------------------
//! \fn PortableMDArray::SwapPortableMDArray()
//  \brief  swap pdata_ pointers of two equally sized PortableMDArrays (shallow
//  swap)
// Does not allocate memory for either array

template <typename T>
void PortableMDArray<T>::SwapPortableMDArray(PortableMDArray<T> &array2) {
  std::swap(pdata_, array2.pdata_);
  return;
}

#endif // _PORTABLE_ARRAYS_HPP_
