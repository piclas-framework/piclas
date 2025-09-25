# Photoionization: Surface Emission via SEE for ray tracing
* **Comparing**: RadiationSurfState.h5, RadiationVolState.h5 and reference_NodeSourceExtGlobal.h5 with reference files, the total number of real electrons in the system with a numerical ref. solution
  * reference_NodeSourceExtGlobal.h5 contains the charge by the electron holes that are created on the VDL surface
* Ray tracing model from which surface emission is calculated
* Particle emission due to secondary electron emission from a surface
* No interpolation, but deposition is required to test the VDL model
* Comparison of the number of emitted electrons with the reference solution
* Different MPF and number of MPI ranks are tested to yield the same result
* Note: Because the volume is exactly 1 cubic metre, the calculated electron density is exactly the number of real electrons in the system
