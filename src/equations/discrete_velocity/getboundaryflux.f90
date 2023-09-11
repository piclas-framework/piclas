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
USE MOD_Equation_Vars_FV  ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars_FV  ,ONLY: nBCByType,BCSideID
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
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(t,tDeriv,Flux,UPrim_master,NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Equation_Vars_FV,ONLY: nBCByType,BCSideID,IniExactFunc_FV,RefState
USE MOD_Riemann
USE MOD_TimeDisc_Vars,ONLY : dt
USE MOD_Equation_FV  ,ONLY: ExactFunc_FV
USE MOD_DistFunc     ,ONLY: MaxwellDistribution, MaxwellScattering, MacroValuesFromDistribution
USE MOD_Equation_Vars_FV,ONLY: DVMSpeciesData,DVMnVelos,DVMVelos,DVMVeloDisc,DVMVeloMax,DVMVeloMin!,DVMBGKModel,DVMWeights,DVMDim,Pi
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: tDeriv      ! deriv
REAL,INTENT(IN)                      :: UPrim_master(     PP_nVar_FV,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:PP_N,0:PP_N,1:nBCSides)
REAL,INTENT(IN)                      :: Face_xGP(        3,0:PP_N,0:PP_N,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux( PP_nVar_FV,0:PP_N,0:PP_N,1:nBCSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: UPrim_boundary(PP_nVar_FV,0:PP_N,0:PP_N)
INTEGER                              :: p,q
INTEGER                              :: iVel,jVel,kVel,upos, upos_sp
REAL                                 :: MacroVal(14), tau, vnormal!, vwall, Sin, Sout, weight, MovTerm, WallDensity
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
      IF(BCState.EQ.0) THEN
        DO q=0,PP_N; DO p=0,PP_N
          CALL ExactFunc_FV(IniExactFunc_FV,t,0,Face_xGP(:,p,q,SideID),UPrim_boundary(:,p,q))
        END DO; END DO
      ELSE
        DO q=0,PP_N; DO p=0,PP_N
          CALL MaxwellDistribution(RefState(:,BCState),UPrim_boundary(:,p,q))
        END DO; END DO
      END IF
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(3) !specular reflection
    IF(DVMVeloDisc.EQ.2.AND.ANY((DVMVeloMin+DVMVeloMax).NE.0.)) THEN
      CALL abort(__STAMP__,'Specular reflection only implemented for zero-centered velocity grid')
    END IF
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        IF (BCState.NE.0) CALL abort(__STAMP__,'DVM specular bc with moving wall not working (yet?)') !THEN
        !   IF (DVMVeloDisc.NE.3) CALL abort(__STAMP__,'DVM specular bc error: moving wall only with Gauss-Hermite disc')
        !   CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q,SideID),dt/2.,tau,1)
        !   Sin=0.
        !   Sout=0.
        !   DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
        !     upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
        !     weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
        !     vnormal = DVMVelos(iVel,1)*NormVec(1,p,q,SideID) + DVMVelos(jVel,2)*NormVec(2,p,q,SideID) + DVMVelos(kVel,3)*NormVec(3,p,q,SideID)
        !     vwall = DVMVelos(iVel,1)*RefState(2,BCState) + DVMVelos(jVel,2)*RefState(3,BCState) + DVMVelos(kVel,3)*RefState(4,BCState)
        !     IF (vnormal.LT.0.) THEN !inflow
        !       Sin = Sin + 2.*weight*vwall*EXP(-(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)/DVMSpeciesData%R_S/MacroVal(5)/2.) &
        !       /DVMSpeciesData%R_S/MacroVal(5)/(2.*Pi*DVMSpeciesData%R_S*MacroVal(5))**(DVMDim/2.)
        !     ELSE IF (vnormal.GT.0.) THEN
        !       Sout = Sout + weight*UPrim_master(upos,p,q,SideID)
        !     ELSE
        !       Sout = Sout + 2.*weight*UPrim_master(upos,p,q,SideID)
        !     END IF
        !   END DO; END DO; END DO
        !   WallDensity = Sout/(1-Sin)
        ! END IF
        DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
          upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
          vnormal = DVMVelos(iVel,1)*NormVec(1,p,q,SideID) + DVMVelos(jVel,2)*NormVec(2,p,q,SideID) + DVMVelos(kVel,3)*NormVec(3,p,q,SideID)
          ! vwall = DVMVelos(iVel,1)*RefState(2,BCState) + DVMVelos(jVel,2)*RefState(3,BCState) + DVMVelos(kVel,3)*RefState(4,BCState)
          IF (ABS(ABS(NormVec(1,p,q,SideID)) - 1.).LE.1e-6) THEN !x-perpendicular boundary
            upos_sp=(DVMnVelos(1)+1-iVel)+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
          ELSE IF (ABS(ABS(NormVec(2,p,q,SideID)) - 1.).LE.1e-6) THEN !y-perpendicular boundary
            upos_sp=iVel+(DVMnVelos(2)-jVel)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
          ELSE IF (ABS(ABS(NormVec(3,p,q,SideID)) - 1.).LE.1e-6) THEN !z-perpendicular boundary
            upos_sp=iVel+(jVel-1)*DVMnVelos(1)+(DVMnVelos(3)-kVel)*DVMnVelos(1)*DVMnVelos(2)
          ELSE
            CALL abort(__STAMP__,'Specular reflection only implemented for boundaries perpendicular to velocity grid')
          END IF
          IF (vnormal.LT.0.) THEN !inflow
            ! IF (BCState.NE.0) THEN
            !   MovTerm = 2.*WallDensity*EXP(-(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)/DVMSpeciesData%R_S/MacroVal(5)/2.) &
            !   /DVMSpeciesData%R_S/MacroVal(5)/(2.*Pi*DVMSpeciesData%R_S*MacroVal(5))**(DVMDim/2.)
            ! ELSE
            !   MovTerm = 0.
            ! END IF
            UPrim_boundary(upos,p,q)=UPrim_master(upos_sp,p,q,SideID)! + MovTerm
            IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
              UPrim_boundary(PP_nVar_FV/2+upos,p,q)=UPrim_master(PP_nVar_FV/2+upos_sp,p,q,SideID)! + MovTerm
            END IF
          ELSE
            UPrim_boundary(upos,p,q)=UPrim_master(upos,p,q,SideID)
            IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
              UPrim_boundary(PP_nVar_FV/2+upos,p,q)=UPrim_master(PP_nVar_FV/2+upos,p,q,SideID)
            END IF
          END IF
        END DO; END DO; END DO
      END DO; END DO
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(4,14) ! maxwell scattering
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      MacroVal(:) = RefState(:,BCState)
      DO q=0,PP_N; DO p=0,PP_N
        CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
        CALL MaxwellScattering(UPrim_boundary(:,p,q),UPrim_master(:,p,q,SideID),NormVec(:,p,q,SideID),1,dt/2.)
      END DO; END DO
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(5) !constant static pressure+temperature inlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q,SideID),dt/2.,tau,1)
        MacroVal(1)=RefState(1,BCState)
        MacroVal(5)=RefState(5,BCState)
        CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
      END DO; END DO
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(6) !constant static pressure outlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
        CALL MacroValuesFromDistribution(MacroVal,UPrim_master(:,p,q,SideID),dt/2.,tau,1)
        MacroVal(5)=RefState(5,BCState)*RefState(1,BCState)/MacroVal(1) !to get the pressure given by refstate
        CALL MaxwellDistribution(MacroVal,UPrim_boundary(:,p,q))
      END DO; END DO
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
    END DO

  CASE(7) !open outlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      UPrim_boundary=UPrim_master(:,:,:,SideID)
      CALL Riemann(Flux(:,:,:,SideID),UPrim_master(:,:,:,SideID),UPrim_boundary,NormVec(:,:,:,SideID))
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
USE MOD_Equation_Vars_FV,ONLY: BCData,nBCByType,BCSideID
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
