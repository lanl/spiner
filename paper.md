---
title: 'Spiner: Performance Portable Routines for Generic, Tabulated, Multi-Dimensional Data'
tags:
  - C++
  - Performance portability
  - GPUs
  - Numerical methods
  - Interpolation
  - Tabulated data
authors:
  - name: Jonah M. Miller^[jonahm\@lanl.gov]
    orcid: 0000-0001-6432-7860
    affiliation: "1, 2"

affiliations:
  - name: CCS-2, Computational Physics and Methods, Los Alamos National Laboratory, Los Alamos, NM
    index: 1
  - name: Center for Theoretical Astrophysics, Los Alamos National Laboratory, Los Alamos, NM
    index: 2
date: 16 March 2022
bibliography: paper.bib

# Summary

We present `Spiner`, a new, performance-portable library for working
with tabulated data. Spiner provides efficient routines for
multi-dimensional interpolation and indexing on CPUs and GPUs,
including interwoven interpolation and indexing access patterns, as
needed for radiation transport. Importantly, Spiner defines a data
format, based on HDF5, that couples the tabulated data to the
information required to interpolate it, which Spiner can read and move
to GPU.

