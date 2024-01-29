# Photoionization in the volume (rectangle) for ray tracing with high-order refinement 
- **Comparing**: the total number of real electrons in the system with a numerical ref. solution
* ray tracing + volume ionization reactions
* reference density in Electrons_ref.csv calculated with the old model and 1x1x1 emission region for volume see and 1e-3 J
* Particle emission due to photoionization of $`H_{2}`$ in a volume only (no surface emission)
* no deposition, no interpolation 
* comparison of the number of emitted electrons with the reference solution (1 MPF and MPI=1)
* different MPF and number of MPI ranks are tested to yield the same result
* Note: Because the volume is exactly 1 cubic metre, the calculated electron density is exactly the number of real electrons in the system
