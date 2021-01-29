---
title: PICLas
logo: ../logo.png
subtitle: Developer Guide
version: Version 2.0
author: 
  - University of Stuttgart, Germany
  - Institute for Aerodynamics and Gas Dynamics (IAG)
  - Institute for Space Systems (IRS)
institute: IAG IRS
date: \today
documentclass: scrreprt
lang: en-US
papersize: a4
fontsize: 11pt
geometry: "left=2.0cm,right=2.0cm,top=3.5cm,bottom=2.5cm"
colorlinks: yes
toc: yes
header-includes:
  - \input{header}
bibliography: ../references.bib
csl: ../ieee.csl
link-citations: true
---

\hypertarget{introduction}{}

# Introduction

 [**PICLas**](http://https://github.com/piclas-framework/piclas)  is a three-dimensional simulation
 framework for Particle-in-Cell, Direct Simulation Monte Carlo and other particle methods that can be coupled for
 the simulation of collisional plasma flows.
 It features a high-order discontinuous 
 Galerkin (DG) simulation module for the solution of the time-dependent Maxwell 
 equations on unstructured hexahedral elements in three space dimensions. 
 The code was specifically designed for very high order accurate simulations on massively parallel 
 systems. 
 It is licensed under GPLv3, written in Fortran and parallelized with MPI. Implemented features are
 
 * Coupled Particle-in-Cell with Direct Simulation Monte Carlo methods
 * Arbitrary order nodal polynomial tensor product basis using Gauss or Gauss Lobatto collocation 
   points for electrostatic and electromagnetic solvers
 * Matching high order curved mesh generation from external mesh formats (CGNS, GMSH) or 
   simple analytic blocks via the open source preprocessor [**HOPR**](http://hopr-project.org) [@Hindenlang2015]
 * Non-conforming interfaces [@Sonntag2017] based on the mortar approach [@Kopriva2001;@Bui2012] (electromagnetic solver)
 * Non-reflecting boundary conditions via CFS-PMLs [@Copplestone2017] (electromagnetic solver)
 * Automatic domain decomposition for parallel simulations based on a space filling curve
 * High order low-storage explicit Runge-Kutta time integration [@Carpenter1994]
 * I/O using the [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) library optimized for massively parallel jobs

## How this documentation is organized

This guide is organized to guide the first implementation steps as well as provide a complete overview of 
the simulation code's features from a developer's point of view.

Preliminary Table of Contents

1. Gitlab Workflow
   1. Issues & Milestones
   2. Release & Deploy
2. Style Guide
3. MPI Implementation
4. Regression Testing
5. Compiler Options

* The first Chapter \ref{chap:git_workflow} shall give an overview over the development workflow within the Gitlab environment, and the necessary steps to create a release, deploy the update to the Collaborative Numerics Group and GitHub.
* The second Chapter \ref{chap:style_guide} describes the rules and guidelines regarding code development such as how the header of functions and subroutines look like.
* Chapter \ref{chap:compiler_options}
