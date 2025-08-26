# 3D dielectric surface charge deposition via inner BCs (=100)

    * Dielectric sphere in sphere example (linear mesh NGeo=1)
    * Particle removal in dielectric regions during initial emission
    * Inner BC (100) which are treated as normal inner faces in the HDG solver
    * Particles impacting on the dielectric surface are swapped to no particles (removed)
    * Their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
    * Checks dielectric interface detection for large number of MPI ranks
    * PartStateBoundary output and DSMC SurfState output
