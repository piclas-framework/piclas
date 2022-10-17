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
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

! Initialize State File Boundary condition
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
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
!> Computes the boundary state for the different boundary conditions.
!==================================================================================================================================
SUBROUTINE GetBoundaryState(SideID,t,Nloc,UPrim_boundary,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BoundaryType,BC
USE MOD_Equation    ,ONLY: ExactFunc
USE MOD_DistFunc     ,ONLY: MaxwellDistribution, MaxwellScattering, MacroValuesFromDistribution
USE MOD_Equation_Vars,ONLY: IniExactFunc, RefStatePrim, DVMSpeciesData
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: SideID                                          !< ID of current side
REAL,INTENT(IN)         :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)      :: Nloc    !< polynomial degree
REAL,INTENT(IN)         :: UPrim_master(  PP_nVar,0:Nloc,0:Nloc) !< inner surface solution
REAL,INTENT(IN)         :: NormVec(                 3,0:Nloc,0:Nloc) !< normal surface vectors
REAL,INTENT(IN)         :: TangVec1(                3,0:Nloc,0:Nloc) !< tangent surface vectors 1
REAL,INTENT(IN)         :: TangVec2(                3,0:Nloc,0:Nloc) !< tangent surface vectors 2
REAL,INTENT(IN)         :: Face_xGP(                3,0:Nloc,0:Nloc) !< positions of surface flux points
REAL,INTENT(OUT)        :: UPrim_boundary(PP_nVar,0:Nloc,0:Nloc) !< resulting boundary state

! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
INTEGER                 :: BCType, BCState
REAL                    :: MacroVal(8), tau
!===================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = BoundaryType(BC(SideID),BC_STATE)

SELECT CASE(BCType)

CASE(1) !Periodic already filled!
  CALL Abort(__STAMP__, &
      "GetBoundaryState called for periodic side!")

CASE(2) !Exact function or refstate
  IF(BCState.EQ.0) THEN
    DO q=0,Nloc; DO p=0,Nloc
      CALL ExactFunc(IniExactFunc,t,0,Face_xGP(:,p,q),UPrim_boundary(:,p,q))
    END DO; END DO
  ELSE
    DO q=0,Nloc; DO p=0,Nloc
      CALL MaxwellDistribution(RefStatePrim(:,BCState),UPrim_boundary(:,p,q))
    END DO; END DO
  END IF

CASE(4) ! maxwell scattering
  MacroVal(:) = RefStatePrim(:,BCState)
  DO q=0,Nloc; DO p=0,Nloc
    CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
    CALL MaxwellScattering(UPrim_boundary(:,p,q),UPrim_master(:,p,q),NormVec(:,p,q),1,t) ! t=tDeriv here
  END DO; END DO

CASE(5) !constant static pressure+temperature inlet
  DO q=0,Nloc; DO p=0,Nloc
    CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q),t,tau,1)
    MacroVal(1)=RefStatePrim(1,BCState)
    MacroVal(5)=RefStatePrim(5,BCState)
    CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
  END DO; END DO

CASE(6) !constant static pressure outlet
  DO q=0,Nloc; DO p=0,Nloc
    CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q),t,tau,1)
    MacroVal(5)=RefStatePrim(5,BCState)*RefStatePrim(1,BCState)/MacroVal(1) !to get the pressure given by refstate
    CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
  END DO; END DO

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in linearscalaradvection/getboundaryflux.f90!')

END SELECT ! BCType

END SUBROUTINE GetBoundaryState

!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(t,tDeriv,Flux,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
USE MOD_Riemann
USE MOD_TimeDisc_Vars,ONLY : dt
USE MOD_DistFunc     ,ONLY: MaxwellDistribution, MacroValuesFromDistribution
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: tDeriv      ! deriv
REAL,INTENT(IN)                      :: UPrim_master(     PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: Face_xGP(        3,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux( PP_nVar,0:PP_N,0:PP_N,1:nBCSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: UPrim_boundary(PP_nVar,0:PP_N,0:PP_N)
REAL                                 :: MacroVal(8), tau
INTEGER                              :: p,q
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!

  CASE(2) !Exact function or refstate
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL GetBoundaryState(SideID,t,PP_N,UPrim_boundary,UPrim_master(:,:,:,SideID),NormVec(:,:,:,SideID) &
        ,TangVec1(:,:,:,SideID),TangVec2(:,:,:,SideID),Face_xGP(:,:,:,SideID))
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(3) ! specular reflection
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N
        DO p=0,PP_N
          CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q,SideID),dt,tau,1)
          MacroVal(2:4) = MacroVal(2:4) - 2.*DOT_PRODUCT(NormVec(1:3,p,q,SideID),MacroVal(2:4))*NormVec(1:3,p,q,SideID)
          CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
        END DO ! p
      END DO ! q
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  ! CASE(4,5,6) ! diffusive or constant static pressure in/outlet
  !   CALL GetBoundaryState(SideID,dt,PP_N,UPrim_boundary,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
  !   CALL Riemann(Flux,UPrim_master,UPrim_boundary,NormVec)

  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryflux.f90!')
  END SELECT ! BCType
END DO

END SUBROUTINE GetBoundaryFlux

!==================================================================================================================================
!> Computes the gradient at a boundary for FV subcells.
!==================================================================================================================================
SUBROUTINE GetBoundaryFVgradient(SideID,t,gradU,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP,sdx_Face)
! MODULES
USE MOD_PreProc
USE MOD_Globals       ,ONLY: Abort
USE MOD_Mesh_Vars     ,ONLY: BoundaryType,BC
USE MOD_Equation     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc, RefStatePrim
USE MOD_TimeDisc_Vars, ONLY : dt
USE MOD_DistFunc      ,ONLY: MaxwellDistribution, MacroValuesFromDistribution, MaxwellScattering
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(IN)   :: t
REAL,INTENT(IN)   :: UPrim_master(PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(OUT)  :: gradU       (PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: NormVec (              3,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: TangVec1(              3,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: TangVec2(              3,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: Face_xGP(              3,0:PP_N,0:PP_N)
REAL,INTENT(IN)   :: sdx_Face(                0:PP_N,0:PP_N,3)

!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar), MacroVal(8), tau
!==================================================================================================================================
BCType  = Boundarytype(BC(SideID),BC_TYPE)
BCState = Boundarytype(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  DO q=0,PP_N; DO p=0,PP_N
    IF(BCState.EQ.0) THEN ! Determine the exact BC state
      CALL ExactFunc(IniExactFunc,t,0,Face_xGP(:,p,q),UPrim_boundary)
    ELSE
      CALL MaxwellDistribution(RefStatePrim(:,BCState),UPrim_boundary)
    END IF
    gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
  END DO; END DO

CASE(3) ! specular reflection
  DO q=0,PP_N
    DO p=0,PP_N
      CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q),dt,tau,2)
      MacroVal(2:4) = MacroVal(2:4) - 2.*DOT_PRODUCT(NormVec(1:3,p,q),MacroVal(2:4))*NormVec(1:3,p,q)
      CALL MaxwellDistribution(MacroVal,UPrim_boundary)
      gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
    END DO ! p
  END DO ! q

CASE(4) ! diffusive
  MacroVal(:) = RefStatePrim(:,BCState)
  DO q=0,PP_N
    DO p=0,PP_N
      CALL MaxwellDistribution(MacroVal,UPrim_boundary)
      CALL MaxwellScattering(UPrim_boundary,UPrim_master(:,p,q),NormVec(:,p,q),2,dt)
      gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
    END DO ! p
  END DO ! q

CASE(5) !constant static pressure+temperature inlet
  DO q=0,PP_N; DO p=0,PP_N
    CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q),dt,tau,2)
    MacroVal(1)=RefStatePrim(1,BCState)
    MacroVal(5)=RefStatePrim(5,BCState)
    CALL MaxwellDistribution(MacroVal,UPrim_boundary)
    gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
  END DO; END DO

CASE(6) !constant static pressure outlet
  DO q=0,PP_N; DO p=0,PP_N
    CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q),dt,tau,2)
    MacroVal(5)=RefStatePrim(5,BCState)*RefStatePrim(1,BCState)/MacroVal(1)
    CALL MaxwellDistribution(MacroVal,UPrim_boundary)
    gradU(:,p,q) = (UPrim_master(:,p,q) - UPrim_boundary) * sdx_Face(p,q,3)
  END DO; END DO

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
       'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryFVgradient

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
