# User Guide

PICLas comes with a vast variety of models and methods. Originally being centered around Particle-in-Cell (PIC) and Direct
Simulation Monte Carlo (DSMC) methods, PICLas has been extended to other particle-based methods, namely Bhatnagar-Gross-Krook (BGK)
and Fokker-Planck (FP) models.
Each of these models, some of which can be combined, offer distinctive features such as

* Coupled Particle-in-Cell with Direct Simulation Monte Carlo methods
* Particle-based Bhatnagar-Gross-Krook (Ellipsoidal Statistical, Shakov, Unified) and Fokker-Planck (Cubic, Ellipsoidal
  Statistical) models for continuum gas flows
* Line-by-line radiation solver and photon Monte Carlo radiation transport
* Raytracing and photo-ionization module
* Variable order nodal polynomial tensor product basis using Gauss or Gauss Lobatto collocation points for electrostatic and
  electromagnetic solvers
* Matching high order curved mesh generation from external mesh formats (CGNS, GMSH) or
  simple analytic blocks via the open source preprocessor [PyHOPE](https://github.com/hopr-framework/PyHOPE) {cite}`Kopper2025`
* Non-conforming interfaces {cite}`Sonntag2017` based on the mortar approach {cite}`Kopriva2001,Bui2012` (electromagnetic solver)
* Non-reflecting boundary conditions via CFS-PMLs {cite}`Copplestone2017` (electromagnetic solver)
* Automatic domain decomposition for parallel simulations based on a space filling curve
* High order low-storage explicit Runge-Kutta time integration {cite}`Carpenter1994`
* I/O using the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) library optimized for massively parallel jobs

This user guide is organized to both guide the first steps as well as provide a complete overview of
the simulation code's features from a user and a developer point of view.

* Chapter {ref}`userguide/installation:Installation` contains step by step instructions from obtaining the source
  code up to running a first simulation and visualizing the simulation results. In addition, it
  provides an overview of the whole simulation framework and the currently implemented features.
* Chapter {ref}`userguide/meshing:Mesh Generation` describes the preprocessing step of creating mesh files via the in-house tool
  [PyHOPE](https://github.com/hopr-framework/PyHOPE) that also handles mesh formats created with external mesh generators
* Chapter {ref}`userguide/workflow:Workflow` outlines the workflow and the visualization of results produced with **PICLas**.
* Chapter {ref}`userguide/features-and-models/index:Features & Models` shall serve as a reference for the models and features implemented in **PICLas**.
* Chapter {ref}`userguide/visu_output:Visualization & Output` presents the options and parameters for the output of particle data, field and flow variables.
* Chapter {ref}`userguide/tools:Tools` lists tools within the **PICLas** repository, including the post-processing tools.
* Simulation tutorials are contained in Chapter {ref}`userguide/tutorials/index:Tutorials`.
* Cluster-specific user guidelines are given in Chapter {ref}`userguide/cluster_guide:Cluster Guidelines`.