#if USE_QDS_DG
MODULE MOD_QDS_Equation_Vars
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
INTEGER,PARAMETER :: QDSnVar=40          ! number of QDS variables: 40
INTEGER           :: QDSIniExactFunc     ! initial condition: exact function
REAL,ALLOCATABLE  :: U_Face_old(:,:,:,:) ! auxiliary face array for first order absorbing BC
!===================================================================================================================================
END MODULE MOD_QDS_Equation_Vars
#endif /*USE_QDS_DG*/
