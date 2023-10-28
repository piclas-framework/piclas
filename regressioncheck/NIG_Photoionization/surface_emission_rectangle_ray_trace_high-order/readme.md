# Photoionization: Surface Emission via SEE for ray tracing with high-order refinement
* **Comparing**: RadiationSurfState.h5 and RadiationVolState.h5 with reference files, the total number of real electrons in the system with a numerical ref. solution
+ hopr mesh is built on-the-fly
  - 0.) single-core hopr run (pre-external)
+ Ray tracing model from which surface emission is calculated
  - 1.) reference files are created by single-core run (pre-external)
* Particle emission due to secondary electron emission from a surface
* No deposition, no interpolation 
* Comparison of the number of emitted electrons with the reference solution
* Different MPF and number of MPI ranks are tested to yield the same result
  - 2.) single + multi node ray tracing runs and compare volume+surface output to .h5 (post-external)
  - 3.) single + multi node plasma restart runs and compare the analytical or reference electron density over time (post-external)
* Note: Because the volume is exactly 1 cubic metre, the calculated electron density is exactly the number of real electrons in the system
