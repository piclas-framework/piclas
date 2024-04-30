# Photoionization in the volume (rectangle) for ray tracing with high-order refinement and bilinear tracking
- **Comparing**: the total number of real electrons in the system with a numerical ref. solution
* convert cubit mesh with hopr to .h5
* ray tracing + volume ionization reactions
* reference density in Electrons_ref.csv calculated with the old model and 1x1x1 emission region for volume see and 1e-3 J
  * the size of the domain in this example is 1mm x 1mm x 1.33 mm, which is much smaller than the other reggies with 1m3 simulation
    size, however, the resulting density must be the same as the irradiation is the same, just on a smaller scale
  * the MPF is therefore chosen much smaller: MPF=0.01
* Particle emission due to photoionization of $`H_{2}`$ in a volume only (no surface emission)
* no deposition, no interpolation 
* comparison of the number of emitted electrons with the reference solution (1 MPF and MPI=1)
* different MPF and number of MPI ranks are tested to yield the same result
* Note: Because the volume is exactly 1 cubic metre, the calculated electron density is exactly the number of real electrons in the system
