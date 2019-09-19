# 2D dielectric surface charge deposition via inner BCs (=100)
Tests
    * particle removal in dielectric regions during initial emission
    * inner BC (100) which are treated as normal inner faces in the HDG solver
    * particles impacting on the dielectric surface are remove (via swap species procedures) and
        their charge is deposited (via cell_volweight_mean) to the surface nodes (cell vertices)
