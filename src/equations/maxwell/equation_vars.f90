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
LOGICAL           :: CentralFlux                             ! flag for central or upwind flux
REAL              :: c_corr
REAL              :: c_corr2    !c_corr^2
REAL              :: c_corr_c   !c_corr*c
REAL              :: c_corr_c2  !c_corr*c^2
REAL              :: eta_c      !(c_corr -1 )*c
REAL              :: fDamping
INTEGER           :: IniExactFunc
INTEGER           :: BCType(6)=-999
INTEGER           :: BoundaryCondition(6,2)
LOGICAL           :: EquationInitIsDone=.FALSE.
LOGICAL           :: DoParabolicDamping
REAL              :: c
REAL              :: c_inv
REAL              :: c2      ! c^2
REAL              :: c2_inv
REAL              :: eps0 
REAL              :: mu0 
REAL              :: smu0
REAL              :: DipoleOmega ! electric dipole angular frequency
REAL              :: tPulse
INTEGER           :: alpha_shape
REAL              :: shapeFuncPrefix
REAL              :: rCutoff
! planar wave and gaussian beam
REAL              :: E_0(1:3)                               !> electric field vector of wave
REAL              :: BeamEta                                !> impedance factor (2*impedance): BeamEta=2.*SQRT(mu0/eps0)
REAL              :: BeamWaveNumber                         !> wave number: BeamWaveNumber= 2*pi/WaveLength
REAL              :: BeamOmegaW                             !> angular frequency: BeamOmegaW = WaveNumber*c
INTEGER           :: BeamIdir1,BeamIdir2,BeamIdir3          !> wave beam direction aux variables
REAL              :: WaveVector(1:3)                        !> wave vector
REAL              :: WaveLength                             !> wave length
REAL,DIMENSION(3) :: WaveBasePoint                          !> wave base point || origin
REAL              :: tFWHM                                  !> time for full wave half maximum
REAL              :: Beam_a0                                !> value to scale max. electric field
REAL              :: BeamAmpFac                             !> decide if pulse maxima is scaled by intensity or a_0 parameter
REAL              :: tDelay                                 !> delay time filter for gaussian beam
REAL              :: I_0                                    !> max. intensity
REAL              :: sigma_t                                !> sigma_t can be used instead of tFWHM
REAL              :: omega_0, omega_0_2inv                  !> spot size and inv of spot size
REAL              :: TEScale                                !> scaling of input TE-wave strength
INTEGER           :: TERotation                             !> left or right rotating TE wave
REAL              :: TEFrequency                            !> frequency of TE wave
LOGICAL           :: TEPulse                                !> Flag for pulsed or continuous wave
LOGICAL           :: DoExactFlux                            !> Flag to switch emission to flux superposition at certain positions
INTEGER           :: FluxDir                                !> direction of flux
! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE  :: nBCByType(:)
INTEGER,ALLOCATABLE  :: BCSideID(:,:)
! can specify BC state
CHARACTER(LEN=255):: BCStateFile

CHARACTER(LEN=255),DIMENSION(8),PARAMETER :: StrVarNames(8)=(/ CHARACTER(LEN=255) :: 'ElectricFieldX', &
                                                                                     'ElectricFieldY', &
                                                                                     'ElectricFieldZ', &
                                                                                     'MagneticFieldX', &
                                                                                     'MagneticFieldY', &
                                                                                     'MagneticFieldZ', &
                                                                                     'Phi'           , &
                                                                                     'Psi           ' /)
!===================================================================================================================================
END MODULE MOD_Equation_Vars
