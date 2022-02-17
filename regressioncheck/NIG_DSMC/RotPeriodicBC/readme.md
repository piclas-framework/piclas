# Rotationally Periodic BCs
- Boundary conditions for particles that perform a 90 degree periodic transformation
- Checks if no particles are moved outside of the simulation domain
- Run different number of MPI ranks (1,2,7,15,25) to test single- and multi-node functionality
- For multi-node (PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE)
  - A warning should be written to std out when the halo region if a proc merely reaches a rot BC but not any further(hence not
  finding a corresponding rotationally periodic side)
  - Debugging information should be written to TestRotatingWall_LostRotPeriodicSides_000.00000000000000000.h5 file
    (cannot be tested automatically unfortunately) when setting CalcMeshInfo=T

