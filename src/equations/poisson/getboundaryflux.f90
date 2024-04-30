!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_GetBoundaryFlux
!===================================================================================================================================
! Fills the boundary part of the flux list (from nInnerSidess+1 to nSides)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
!INTERFACE GetBoundaryFlux
  !MODULE PROCEDURE GetBoundaryFlux
!END INTERFACE

INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

!PUBLIC::GetBoundaryFlux
PUBLIC:: InitBC,FinalizeBC
!===================================================================================================================================

CONTAINS


SUBROUTINE InitBC()
!===================================================================================================================================
! Initialize boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: BCData,nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone)) CALL abort(__STAMP__,&
    "InitBC not ready to be called or already called.")
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  ! should not be required || example for MaxBCState
  !IF((locType.NE.22).AND.locType.NE.3) MaxBCState = MAX(MaxBCState,locState)
  !IF((locType.EQ.4).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No temperature (refstate) defined for BC_TYPE',locType)
  !IF((locType.EQ.23).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow Mach number in refstate (x,Ma,x,x,x) defined for BC_TYPE',locType)
  !IF((locType.EQ.24).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow pressure in refstate defined for BC_TYPE',locType)
  !IF((locType.EQ.25).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow pressure in refstate defined for BC_TYPE',locType)
  !IF((locType.EQ.27).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No inflow refstate (dens,v1,v2,v3,pressure) in refstate defined for BC_TYPE',locType)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

! Sanity check for BCs
!IF(MaxBCState.GT.nRefState)&
!  CALL abort(__STAMP__,&
!    'ERROR: Boundary RefState not defined! (MaxBCState,nRefState):',MaxBCState,REAL(nRefState))


! Initialize State File Boundary condition
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
  IF(locType.EQ.20)THEN
    ! Allocate buffer array to store temp data for all BC sides
    ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
    BCData=0.
    !CALL ReadBCFlow(BCStateFile)
    !CALL abort(__STAMP__,&
    !     'no BC defined in maxwell/getboundaryflux.f90!')
    EXIT
  END IF
END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC
!SUBROUTINE GetBoundaryFlux(F_Face,BCType,BCState,xGP_Face,normal,t,tDeriv,U_Face) ! old
!SUBROUTINE GetBoundaryFlux(t,tDeriv, Flux, U_Minus, NormVec, TangVec1, TangVec2, BCFace_xGP)
!!===================================================================================================================================
!! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!! BCType: 1...periodic, 2...exact BC
!! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!!              SUBROUTINE CalcSurfInt
!! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!!===================================================================================================================================
!! MODULES
!USE MOD_Globals,ONLY:Abort
!USE MOD_PreProc
!USE MOD_Mesh_Vars    ,  ONLY:nBCSides
!USE MOD_Equation,ONLY:ExactFunc
!USE MOD_Equation_Vars,ONLY:IniExactFunc
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!INTEGER, INTENT(IN)                  :: BCType,BCState
!!REAL,INTENT(IN)                   :: U_Face(PP_nVar,0:PP_N,0:PP_N)
!!REAL,INTENT(IN)                      :: xGP_Face(3,0:PP_N,0:PP_N)
!!REAL,INTENT(IN)                      :: normal(3,0:PP_N,0:PP_N)
!REAL,INTENT(IN)                      :: t
!INTEGER,INTENT(IN)                   :: tDeriv
!REAL,INTENT(IN)                      :: U_Minus(     PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
!REAL,INTENT(IN)                      :: NormVec(           3,0:PP_N,0:PP_N,1:nBCSides)
!REAL,INTENT(IN)                      :: TangVec1(          3,0:PP_N,0:PP_N,1:nBCSides)
!REAL,INTENT(IN)                      :: TangVec2(          3,0:PP_N,0:PP_N,1:nBCSides)
!REAL,INTENT(IN)                      :: BCFace_xGP(        3,0:PP_N,0:PP_N,1:nBCSides)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!REAL,INTENT(OUT)                     :: Flux(        PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                              :: p,q
!INTEGER                              :: BCType,BCState,nBCLoc
!REAL                                 :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
!!===================================================================================================================================
!SELECT CASE(BCType)
!CASE(1) !Periodic already filled!
!CASE DEFAULT ! unknown BCType
  !CALL abort(&
!__STAMP__&
!,'no BC defined in navierstokes/getboundaryflux.f90!',999,999.)
!END SELECT ! BCType
!END SUBROUTINE GetBoundaryFlux


SUBROUTINE FinalizeBC()
!===================================================================================================================================
! Initialize boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC
END MODULE MOD_GetBoundaryFlux
