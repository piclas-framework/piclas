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
! Contains FillBoundary (which depends on the considered equation)
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
INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

PUBLIC::GetBoundaryFlux
PUBLIC:: InitBC,FinalizeBC
!===================================================================================================================================

CONTAINS



!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars_FV     ,ONLY: EquationInitIsDone_FV
#if USE_HDG
USE MOD_Equation_Vars   ,ONLY:nBCByType,BCSideID
#else
USE MOD_Equation_Vars_FV,ONLY:nBCByType,BCSideID
#endif
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,BoundaryType_FV,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState, locType_FV, locState_FV
INTEGER :: MaxBCState,MaxBCStateGlobal
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone_FV))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locType_FV =BoundaryType_FV(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  locState_FV=BoundaryType_FV(BC(iSide),BC_STATE)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

! Initialize State File Boundary condition
DO i=1,nBCs
  locType =BoundaryType_FV(i,BC_TYPE)
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

!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(t,tDeriv,Flux,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals         ,ONLY: Abort
USE MOD_Mesh_Vars       ,ONLY: nBCSides,nBCs,BoundaryType_FV
#if USE_HDG
USE MOD_Equation_Vars   ,ONLY: nBCByType,BCSideID,E
USE MOD_Mesh_Vars       ,ONLY: SideToElem
#else
USE MOD_Equation_Vars_FV,ONLY: nBCByType,BCSideID
#endif
USE MOD_Equation_Vars_FV,ONLY: IniExactFunc_FV,RefState_FV
USE MOD_Gradient_Vars   ,ONLY: Grad_dx_master
USE MOD_Riemann
USE MOD_Equation_FV     ,ONLY: ExactFunc_FV
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: tDeriv      ! deriv
REAL,INTENT(IN)                      :: UPrim_master(     PP_nVar_FV,0:0,0:0,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,0:0,0:0,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:0,0:0,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:0,0:0,1:nBCSides)
REAL,INTENT(IN)                      :: Face_xGP(        3,0:0,0:0,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux( PP_nVar_FV,0:0,0:0,1:nBCSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,SideID,ElemID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: UPrim_boundary(PP_nVar_FV,0:0,0:0), GradSide
INTEGER                              :: p,q
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType_FV(iBC,BC_TYPE)
  BCState=BoundaryType_FV(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  
  CASE(2) !Exact function or RefState_FV
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      ElemID=SideToElem(S2E_ELEM_ID,SideID)
      IF(BCState.EQ.0) THEN
        CALL ExactFunc_FV(IniExactFunc_FV,t,0,Face_xGP(:,0,0,SideID),UPrim_boundary(:,0,0))
      ELSE
        UPrim_boundary(:,0,0)=RefState_FV(:,BCState)
      END IF
#ifdef drift_diffusion
      GradSide=(UPrim_master(1,0,0,SideID)-UPrim_boundary(1,0,0))/(2*SQRT((Grad_dx_master(1,SideID))**2 &
                                       +(Grad_dx_master(2,SideID))**2 &
                                       +(Grad_dx_master(3,SideID))**2))
#endif
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID),&
                   GradSide,E(1:3,0,0,0,ElemID),E(1:3,0,0,0,ElemID))
    END DO
  
  CASE(3) !von Neumann
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      UPrim_boundary=UPrim_master(:,:,:,SideID)
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID), &
                   GradSide=0.,E_L=E(1:3,0,0,0,ElemID),E_R=E(1:3,0,0,0,ElemID))
    END DO

  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in DVM/getboundaryflux.f90!')
  END SELECT ! BCType
END DO

END SUBROUTINE GetBoundaryFlux

SUBROUTINE FinalizeBC()
!===================================================================================================================================
! Initialize boundary conditions
!===================================================================================================================================
! MODULES
#if USE_HDG
USE MOD_Equation_Vars   ,ONLY:nBCByType,BCSideID
#else
USE MOD_Equation_Vars_FV,ONLY:nBCByType,BCSideID
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC


END MODULE MOD_GetBoundaryFlux
