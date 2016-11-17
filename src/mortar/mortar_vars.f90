MODULE MOD_Mortar_Vars
!===================================================================================================================================
!> Variables used for mortars: mortar interpolation and projection matrices
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
!> 1D-Mortar Operator: interpolation full interval 0: [-1,1] to left interval 1: [-1,0] and right intervall 2: [0,1]
REAL,ALLOCATABLE,TARGET :: M_0_1(:,:),M_0_2(:,:)
!> 1D-Mortar Operator: projection left interval 1: [-1,0] and right intervall 2: [0,1] to full intervall 0: [-1,1]
REAL,ALLOCATABLE,TARGET :: M_1_0(:,:),M_2_0(:,:)
LOGICAL                 :: MortarInitIsDone=.FALSE. !< marks whether mortar init routines are complete
!===================================================================================================================================
END MODULE MOD_Mortar_Vars
