# Photoionization: Surface Emission via SEE for ray tracing
* **Comparing**: RadiationSurfState.h5 and RadiationVolState.h5 with reference files, the total number of real electrons in the system with a numerical ref. solution
* Ray tracing model from which surface emission is calculated
* Particle emission due to secondary electron emission from a surface
* No deposition, no interpolation 
* Comparison of the number of emitted electrons with the reference solution
* Different MPF and number of MPI ranks are tested to yield the same result
* Note: Because the volume is exactly 1 cubic metre, the calculated electron density is exactly the number of real electrons in the system
