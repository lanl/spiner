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
  - name: Daniel Holladay
    affiliation: "2, 3"
  - name: Chad D. Meyer
    affiliation: 4
  - name: Joshua C. Dolence
    affiliation: "1, 2"
  - name: Sriram Swaminarayan
    affiliation: 3
  - name: Christopher M. Mauney
    affiliation: "2, 5"
  - name: Karen Tsai
    affiliation: 3
affiliations:
  - name: CCS-2, Computational Physics and Methods, Los Alamos National Laboratory, Los Alamos, NM
    index: 1
  - name: Center for Theoretical Astrophysics, Los Alamos National Laboratory, Los Alamos, NM
    index: 2
  - name: CCS-7, Applied Computer Science, Los Alamos National Laboratory, Los ALamos, NM
    index: 3
  - name: XCP-4, Continuum Models and Numerical Methods, Los Alamos National Laboratory, Los ALamos, NM
    index: 4
  - name: HPC-ENV, HPC Environments, Los Alamo National Laboratory, Los Alamos, NM
    index: 5
date: 16 March 2022
bibliography: paper.bib
---

# Summary

We present `Spiner`, a new, performance-portable library for working
with tabulated data. Spiner provides efficient routines for
multi-dimensional interpolation and indexing on CPUs and GPUs,
including interwoven interpolation and indexing access patterns, as
needed for radiation transport. Importantly, Spiner defines a data
format, based on HDF5, that couples the tabulated data to the
information required to interpolate it, which Spiner can read and move
to GPU.

# Statement of Need

As Moore's law comes to an end, more and more performance comes from
specialized hardware, such as GPUs. A key tool in the toolbox for many
scientific codes is tabulated data. Fluid and continuum dynamics codes
often encapsulate the equation of state as data tabulated in density
and temperature, for example as published in the Sesame database
[@sesame] or the stellar collapse database
[@stellarcollapsetables,@stellarcollapseweb]. Radiation transport,
such as that performed by [@fornax] and [@nubhlight] uses emissivity
and absorption opacity on tables such as those computed in
[@SullivanWeak]. Spiner is now used in the open-source and on-going
Singularity-EOS [@singularityeos], Singularity-Opac
[@singularityopac], and Phoebus [@phoebus] projects, which have
separate code papers in-prep.
