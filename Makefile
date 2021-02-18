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

USE_HDF5=no

CC=g++
INC= spiner_types.hpp \
	databox.hpp \
	interpolation.hpp
GCH=$(INC:.hpp=.hpp.gch)

H5CC=h5c++
H5_DIR=/usr/local/hdf5-parallel
H5_INCLUDE=-I${H5_DIR}/include
H5_LIB=-L${H5_DIR}/lib
H5_LINK=-lhdf5_hl -lhdf5 -Wl,-rpath=${H5_DIR}/lib/

CATCH_INCLUDE=-ICatch2/single_include/catch2
PORTS_INCLUDE=-Iports-of-call
BASE_FLAGS=-std=c++14 -Wall -fdiagnostics-color=always -g

ifeq ($(USE_HDF5),yes)
	CC=${H5CC}
	INCLUDE_FLAGS=${CATCH_INCLUDE} ${PORTS_INCLUDE} ${H5_INCLUDE}
	LIB_FLAGS=${H5_LIB}
	LINK_FLAGS=${H5_LINK}
	CFLAGS=${BASE_FLAGS} ${INCLUDE_FLAGS} -DSPINER_USE_HDF
else
	INCLUDE_FLAGS=${CATCH_INCLUDE} ${PORTS_INCLUDE}
	LIB_FLAGS=
	LINK_FLAGS=
	CFLAGS=${BASE_FLAGS} ${INCLUDE_FLAGS}
endif

LFLAGS=${LIB_FLAGS} ${LINK_FLAGS}

default: test

all: test convergence

test: test.bin
	./test.bin

convergence: convergence.bin
	./convergence.bin
	python plot_convergence.py

%.bin: %.o
	$(CC) ${CFLAGS} ${LFLAGS} -g -o $@ $^

%.o: %.cpp ${INC}
	$(CC) ${CFLAGS} -c -o $@ $<

.PHONY: default all test convergence clean info

clean:
	$(RM) test.bin convergence.bin
	$(RM) test.o convergence.o databox.o
	$(RM) ${GCH}
	$(RM) convergence.dat
	$(RM) databox.sp5
