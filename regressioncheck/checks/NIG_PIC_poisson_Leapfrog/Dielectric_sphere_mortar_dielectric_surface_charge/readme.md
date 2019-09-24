# 3D dielectric surface charge deposition via inner BCs (=100)
Tests
    * Dielectric sphere in sphere example (linear mesh NGeo=1)
    * Particle removal in dielectric regions during initial emission
    * Inner BC (100) which are treated as normal inner faces in the HDG solver
    * Particles impacting on the dielectric surface are swapped to heavy particle and slowed down due to extremely cold surface (accommodation=1)
    * Their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
