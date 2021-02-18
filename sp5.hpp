#ifndef _SPINER_SP5_HPP_
#define _SPINER_SP5_HPP_
//======================================================================
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
//======================================================================

// This file contains strings defined for use accross the SP5 data
// format

namespace SP5 {

  namespace DB {
    constexpr char FILENAME[] = "databox.sp5";
    constexpr char GRPNAME[]  = "databox";
    constexpr char DSETNAME[] = "data";
    constexpr char RANKNAME[] = "rank";
    constexpr char DIMSNAME[] = "dims";
    constexpr char IDXSNAME[] = "index_types";
    constexpr char IDXINFONAME[] = "index_types_info";
    constexpr char IDXINFO[]  = "Interpolated:0\nNamed:1\nIndexed:2";
    constexpr char GRIDNAME[] = "grids";
    const std::string GRID_FORMAT[] = {"grid_[","]"};
  }

  namespace RG1D {
    constexpr char RANGE_NAME[] = "range";
    constexpr char N[] = "npoints";
    constexpr char RANGE_INFONAME[] = "range columns";
    constexpr char RANGE_INFO[] = "[0]:min [1]:max [2]:dx";
    constexpr int RANGE_RANK = 1;
  }
  
}

#endif // _SPINER_SP5_HPP_
