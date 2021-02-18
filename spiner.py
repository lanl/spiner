#!/usr/bin/env python

#======================================================================
# © (or copyright) 2019-2021. Triad National Security, LLC. All rights
# reserved.  This program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
# operated by Triad National Security, LLC for the U.S.  Department of
# Energy/National Nuclear Security Administration. All rights in the
# program are reserved by Triad National Security, LLC, and the
# U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting
# on its behalf a nonexclusive, paid-up, irrevocable worldwide license
# in this material to reproduce, prepare derivative works, distribute
# copies to the public, perform publicly and display publicly, and to
# permit others to do so.
#======================================================================

from enum import IntEnum

class IndexType(IntEnum):
    Interpolated = 0
    Named = 1
    Indexed = 2

class DataBox:
    def __init__(self,
                 indices=[],
                 grids={},
                 data=None):
        self._indices=indices
        self._grids=grids
        self._data=data

    def __getitem__(self,key):
        return self._data[key]

    def shape(self):
        return self._data.shape

    def grids(self,i):
        return self._grids[i]

    def indices(self,i):
        return self._indices

    def data(self):
        return self._data

    @classmethod
    def fromHDF(cls,loc):
        import h5py
        import numpy as np
        indices = loc.attrs['index_types']
        grids = {}
        for i in range(indices.shape[0]):
            if indices[i] == IndexType.Interpolated:
                gridname = 'grids/grid_[{}]'.format(i+1)
                xmin,xmax,dx = loc[gridname][()]
                nx = loc[gridname].attrs['npoints'][0]
                grids[i] = np.linspace(xmin,xmax,nx)
        data = loc['data'][()]
        return cls(indices,grids,data)
