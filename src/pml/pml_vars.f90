MODULE MOD_PML_Vars
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! PML region damping factor
LOGICAL                               :: DoPML
REAL,ALLOCATABLE                      :: PMLzeta(:,:,:,:,:)   ! damping factor in xyz
INTEGER                               :: nPMLElems
INTEGER,ALLOCATABLE                   :: PMLtoElem(:)    
INTEGER,ALLOCATABLE                   :: ElemtoPML(:)    
INTEGER                               :: nTotalPML
! PML auxiliary variables P_t=E & Q_t=B
REAL,ALLOCATABLE                      :: U2(:,:,:,:,:)     ! P=U2(1:3) and Q=U2(4:6)
REAL,ALLOCATABLE                      :: U2t(:,:,:,:,:)    ! P=U2(1:3) and Q=U2(4:6)
! Probes
TYPE ProbeVars
  REAL, ALLOCATABLE                   :: Distance(:,:,:,:,:)        ! probe distance to volume gauss points
  REAL                                :: Coordinates(3,3)   ! probe coordinates for 3 probes
  INTEGER                             :: Element(3)         ! iElem of the probes
  INTEGER                             :: iElemMinLoc(3,4)   ! i,j,k,iElem of the probes
END TYPE
TYPE (ProbeVars)                      :: Probes
INTEGER                               :: PMLzetaShape ! shape functions for particle deposition and PML damping coefficient
INTEGER                               :: PMLwriteZeta ! output zeta field for debug
INTEGER                               :: PMLspread
REAL,DIMENSION(6)                     :: xyzPhysicalMinMax ! physical boundary coordinates, outside = PML region
REAL                                  :: PMLzeta0             ! damping constant for PML region
LOGICAL                               :: PMLzetaNorm
!===================================================================================================================================
END MODULE MOD_PML_Vars
