# Rotationally Periodic BCs
- Boundary conditions for particles that perform a 90 degree periodic transformation
- Checks if no particles are moved outside of the simulation domain
- Run different number of MPI ranks (1,2,7,15,25) to test single- and multi-node functionality
- For multi-node (PICLAS_SHARED_MEMORY = OMPI_COMM_TYPE_CORE)
  - A warning should be written to std out when the halo region if a proc merely reaches a rot BC but not any further(hence not
  finding a corresponding rotationally periodic side)
  - Debugging information should be written to TestRotatingWall_LostRotPeriodicSides_000.00000000000000000.h5 file
    (cannot be tested automatically unfortunately) when setting CalcMeshInfo=T
- for debugging, use 2 procs and smaller mesh via

        DEFVAR=(INT):    i01 = 2   ! no. elems in left and right block
        DEFVAR=(INT):    i02 = 2   ! no. elems in upper block (should be twice the value of i01)

        DEFVAR=(INT):    ir1 = 3   ! no. elems in r for first ring
        DEFVAR=(REAL):   r01 = 3.5 ! middle square dim
        DEFVAR=(REAL):   r02 = 7.0 ! middle square dim
        DEFVAR=(REAL):   s0  = 0.2857142857142857 ! middle square dim


        DEFVAR=(INT):    iz1 = 3
        DEFVAR=(INT):    iz2 = 6
        DEFVAR=(REAL):   lz = 1    ! length of domain in z

        DEFVAR=(REAL):   f1 = 1.    ! stretching factor in first ring

