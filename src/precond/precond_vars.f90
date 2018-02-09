#include "boltzplatz.h"

MODULE MOD_Precond_Vars
!===================================================================================================================================
! Contains global variables used by the Timedisc modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE      :: invP(:,:,:)  !inverse of block Jacobian for each element (1:nDOF_elem,1:nDOFelem,1:PP_nElems)
INTEGER               :: PrecondType
INTEGER               :: DebugMatrix
LOGICAL               :: doVol,doSurf
LOGICAL               :: PrecondInitIsDone
REAL,ALLOCATABLE      :: invBJ(:,:,:,:,:,:)  !inverse of block Jacobian for each DOF (1:PP_nVar,1:PP_nVar,0:N,0:N,0:N,1:PP_nElems)
REAL,ALLOCATABLE      :: invJ(:,:,:,:,:)  ! inverse of Jacobian (1:PP_nVar,0:N,0:N,0:N,1:PP_nElems)
INTEGER               :: PrecondMethod
INTEGER,ALLOCATABLE   :: NeighborElemID(:,:)  ! 1:6, 1:PP_nElems
INTEGER,ALLOCATABLE   :: MatPatternElemID(:,:) ! number of row, ElemID, row equals diag-elemid
INTEGER,ALLOCATABLE   :: MatPattern(:,:)       ! row, locSideID
INTEGER,ALLOCATABLE   :: Lower(:,:)            ! index range for MatPattern, containing the lower triangular matrix
INTEGER,ALLOCATABLE   :: Upper(:,:)            ! index range for MatPattern, containing the upper triangular matrix
REAL,ALLOCATABLE      :: ProcJacobian(:,:,:,:)
INTEGER               :: BlockSize
INTEGER               :: nBlockSize
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: invXi,invEta,invZeta
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: dRdXi,dRdEta,dRdZeta
LOGICAL               :: BuildNvecisDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_Precond_Vars
