---
title: PICLas Documentation
author: 
  - Institute for Aerodynamics and Gas Dynamics (IAG)
  - Institute for Space Systems (IRS)
  - University of Stuttgart, Germany
date: \today
documentclass: scrreprt
lang: en-US
papersize: a4
colorlinks: yes
toc: yes
header-includes:
  - \input{header}

---

\hypertarget{introduction}{}

# Introduction

 [**PICLas**](http://https://github.com/piclas-framework/piclas)  is a three-dimensional simulation
 framework for Particle-in-Cell and Direct Simulation Monte Carlo methods that can be coupled for
 the simulation of collisional plasma flows.
 It features a high-order discontinuous 
 Galerkin (DG) simulation module for the solution of the time-dependent Maxwell 
 equations on unstructured hexahedral elements in three space dimensions. 
 The code was specifically designed for very high order accurate simulations on massively parallel 
 systems. 
 It is licensed under GPLv3, written in Fortran and parallelized with MPI. Implemented features are
 
 * Coupled Particle-in-Cell with Direct Simulation Monte Carlo methods
 * Arbitrary order nodal polynomial tensor product basis using Gauss or Gauss Lobatto collocation points for electrostatic and electromagnetic solvers
 * Matching high order curved mesh generation from external mesh formats (CGNS, GMSH) or simple analytic blocks via the open source preprocessor [**HOPR**](http://hopr-project.org)
 * Nonconforming interfaces based on the mortar approach (electromagnetic solver)
 * Non-reflecting boundary conditions via CFS-PMLs (electromagnetic solver)
 * Automatic domain decomposition for parallel simulations based on a space filling curve
 * High order low-storage explicit Runge-Kutta time integration
 * I/O using the HDF5 library optimized for massively parallel jobs
 
## How this documentation is organized

This user guide is organized to both guide the first steps as well as provide a complete overview of 
the simulation code's features from a user and a developer point of view.

* Chapter \ref{chap:gettingstarted} contains step by step instructions from obtaining the source 
  code up to running a first simulation and visualizing the simulation results. In addition, it 
  provides an overview of the whole simulation framework and the currently implemented features.

* Chapter \ref{chap:codeoptions} lists the most important compiler and runtime options.

* Chapter \ref{chap:workflow} is meant as a complete user guide with a detailed description how to 
  use and apply the features of **PICLas** from a user's point of view. This includes setting up 
  solver settings, initial and boundary conditions, the mesh interface, parallel execution and the 
  currently available post processing capabilities.

* Chapter \ref{chap:toolsoverview} lists all tools contained in the **PICLas** repository, including 
  **POSTI** post-processing tools. 

* Simulation tutorials are contained in Chapter \ref{chap:tutorials}.

* The unit test system used to test key routines with CTest is described in Chapter \ref{chap:unittest}.
