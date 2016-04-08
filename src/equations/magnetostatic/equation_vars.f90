MODULE MOD_Equation_Vars
!===================================================================================================================================
! Contains the constant Advection Velocity Vector used for the linear scalar advection equation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL              :: Pi
REAL              :: IniWavenumber(3) ! wavenumbers in 3 directions (sinus periodic with exactfunc=6) 
INTEGER           :: IniExactFunc
REAL              :: IniCenter(3)
REAL              :: IniAmplitude 
REAL              :: IniHalfwidth 

REAL,ALLOCATABLE  :: chitens(:,:,:,:,:,:) !diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitensInv(:,:,:,:,:,:) ! inverse of diffusion 3x3 tensor on each gausspoint
REAL,ALLOCATABLE  :: chitens_face(:,:,:,:,:) !diffusion 3x3 tensor on each face gausspoint

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames(3)=(/ CHARACTER(LEN=255) :: 'MagneticFieldX'           , &
                                                                                     'MagneticFieldY', &
                                                                                     'MagneticFieldZ'/)

LOGICAL           :: EquationInitIsDone=.FALSE.
REAL              :: eps0 
REAL              :: mu0, smu0 
REAL              :: c
REAL              :: c2
REAL              :: c2_inv
REAL              :: c_inv
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
REAL,ALLOCATABLE  :: B(:,:,:,:,:)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
