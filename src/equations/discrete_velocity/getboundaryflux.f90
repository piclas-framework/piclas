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
USE MOD_Globals
USE MOD_Equation_Vars_FV  ,ONLY: EquationInitIsDone_FV
#if USE_HDG
USE MOD_Equation_Vars     ,ONLY:nBCByType,BCSideID
#else
USE MOD_Equation_Vars_FV  ,ONLY:nBCByType,BCSideID
#endif
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone_FV))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitBC not ready to be called or already called.")
END IF

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
SUBROUTINE GetBoundaryFlux(t,Flux,UPrim_master,NormVec,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: nBCSides,nBCs
USE MOD_Mesh_Vars_FV ,ONLY: BoundaryType_FV
#if USE_HDG
USE MOD_Equation_Vars   ,ONLY: nBCByType,BCSideID
#else
USE MOD_Equation_Vars_FV,ONLY: nBCByType,BCSideID
#endif
USE MOD_Equation_Vars_FV,ONLY: IniExactFunc_FV,RefState_FV
USE MOD_Riemann
USE MOD_TimeDisc_Vars,ONLY : dt
USE MOD_Equation_FV  ,ONLY: ExactFunc_FV
USE MOD_DistFunc     ,ONLY: MaxwellDistribution, MaxwellScatteringDVM, MacroValuesFromDistribution, MoleculeRelaxEnergy
USE MOD_Equation_Vars_FV,ONLY: DVMDim,DVMSpecData,DVMVeloDisc,DVMnSpecies, DVMMethod, DVMnMacro, BCTempGrad, DVMnSpecTot
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
REAL,INTENT(IN)                      :: UPrim_master(     PP_nVar_FV,1:nBCSides)
REAL,INTENT(IN)                      :: NormVec(           3,1:nBCSides)
REAL,INTENT(IN)                      :: Face_xGP(        3,1:nBCSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux( PP_nVar_FV,1:nBCSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,SideID
INTEGER                              :: BCType,BCState,nBCLoc
REAL                                 :: UPrim_boundary(PP_nVar_FV)
INTEGER                              :: iVel,jVel,kVel,upos, upos_sp, iSpec, vFirstID, vLastID, BCorient
REAL                                 :: MacroVal(DVMnMacro), vnormal
REAL                                 :: MacroValInside(DVMnMacro,DVMnSpecTot),rho,Pr,tau,prefac,relaxFac
REAL                                 :: Erot(DVMnSpecTot), ErelaxTrans, ErelaxRot(DVMnSpecies)
REAL                                 :: Evib(DVMnSpecTot), ErelaxVib(DVMnSpecies)
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType_FV(iBC,BC_TYPE)
  BCState=BoundaryType_FV(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!

  CASE(2) !Exact function or refstate
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      IF(BCState.EQ.0) THEN
        CALL ExactFunc_FV(IniExactFunc_FV,t,Face_xGP(:,SideID),UPrim_boundary(:))
      ELSE
        vFirstID=1
        vLastID=0
        DO iSpec=1,DVMnSpecies
          vLastID = vLastID + DVMSpecData(iSpec)%nVar
          CALL MaxwellDistribution(RefState_FV(:,iSpec,BCState),UPrim_boundary(vFirstID:vLastID),iSpec)
          vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
        END DO
      END IF
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(3,13) !specular reflection
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      IF (BCState.NE.0) CALL abort(__STAMP__,'DVM specular bc with moving wall not impemented')
      IF      (ALMOSTZERO(ABS(ABS(NormVec(1,SideID)) - 1.))) THEN !x-perpendicular boundary
        BCorient=1
      ELSE IF (ALMOSTZERO(ABS(ABS(NormVec(2,SideID)) - 1.))) THEN !y-perpendicular boundary
        BCorient=2
      ELSE IF (ALMOSTZERO(ABS(ABS(NormVec(3,SideID)) - 1.))) THEN !z-perpendicular boundary
        BCorient=3
      ELSE
        CALL abort(__STAMP__,'Specular reflection only implemented for boundaries perpendicular to velocity grid')
      END IF
      vFirstID=0
      DO iSpec=1,DVMnSpecies
        ASSOCIATE(Sp => DVMSpecData(iSpec))
        IF (BCType.EQ.13.AND.Sp%InterID.NE.4) THEN
           ! absorb ions for plasma sheath
          UPrim_boundary(vFirstID+1:vFirstID+Sp%nVar) = 0.
          vFirstID = vFirstID + Sp%nVar
          CYCLE
        END IF
        IF(DVMVeloDisc(iSpec).EQ.2.AND.(Sp%VeloMin(BCorient)+Sp%VeloMax(BCorient)).NE.0.) THEN
          print*, iSpec, BCorient, Sp%VeloMin(BCorient), Sp%VeloMax(BCorient)
          CALL abort(__STAMP__,'Specular reflection only implemented for zero-centered velocity grid')
        END IF
        DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
          upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
          vnormal = Sp%Velos(iVel,1)*NormVec(1,SideID) + Sp%Velos(jVel,2)*NormVec(2,SideID) + Sp%Velos(kVel,3)*NormVec(3,SideID)
          SELECT CASE(BCorient)
          CASE(1) !x-perpendicular boundary
            upos_sp=(Sp%nVelos(1)+1-iVel)+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
          CASE(2) !y-perpendicular boundary
            upos_sp=iVel+(Sp%nVelos(2)-jVel)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
          CASE(3) !z-perpendicular boundary
            upos_sp=iVel+(jVel-1)*Sp%nVelos(1)+(Sp%nVelos(3)-kVel)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
          END SELECT
          IF (vnormal.LT.0.) THEN !inflow
            UPrim_boundary(upos)=UPrim_master(upos_sp,SideID)
            IF (DVMDim.LT.3) THEN
              UPrim_boundary(Sp%nVarReduced+upos)=UPrim_master(Sp%nVarReduced+upos_sp,SideID)
            END IF
            IF (Sp%Xi_Rot.GT.0.) THEN
              UPrim_boundary(Sp%nVarErotStart+upos)=UPrim_master(Sp%nVarErotStart+upos_sp,SideID)
            END IF
            IF (Sp%T_Vib.GT.0.) THEN
              UPrim_boundary(Sp%nVarEvibStart+upos)=UPrim_master(Sp%nVarEvibStart+upos_sp,SideID)
            END IF
          ELSE
            UPrim_boundary(upos)=UPrim_master(upos,SideID)
            IF (DVMDim.LT.3) THEN
              UPrim_boundary(Sp%nVarReduced+upos)=UPrim_master(Sp%nVarReduced+upos,SideID)
            END IF
            IF (Sp%Xi_Rot.GT.0.) THEN
              UPrim_boundary(Sp%nVarErotStart+upos)=UPrim_master(Sp%nVarErotStart+upos,SideID)
            END IF
            IF (Sp%T_Vib.GT.0.) THEN
              UPrim_boundary(Sp%nVarEvibStart+upos)=UPrim_master(Sp%nVarEvibStart+upos,SideID)
            END IF
          END IF
        END DO; END DO; END DO
        vFirstID = vFirstID + Sp%nVar
        END ASSOCIATE
      END DO ! iSpec
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(4,14,24,25) ! maxwell scattering
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      IF (DVMMethod.GT.0) THEN
        CALL MacroValuesFromDistribution(MacroValInside,UPrim_master(:,SideID),dt/2.,tau,1,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot,Evib=Evib)
        CALL MoleculeRelaxEnergy(ErelaxTrans,ErelaxRot,ErelaxVib,MacroValInside(5,DVMnSpecTot),Erot(1:DVMnSpecies),Evib(1:DVMnSpecies),Pr)
        IF (dt.EQ.0.) THEN
          prefac = 1.
        ELSE
          SELECT CASE(DVMMethod)
          CASE(1)
            prefac = 0.
            IF (tau.GT.0.) THEN
              relaxFac = dt/tau/2.
              IF (CHECKEXP(relaxFac)) THEN
                prefac = 2.*tau*(1.-EXP(-relaxFac))/dt ! f from f2~
              END IF
            END IF
          CASE(2)
            prefac = 2.*tau/(2.*tau+dt/2.)
          END SELECT
        END IF
      END IF
      vFirstID=1
      vLastID=0
      DO iSpec=1,DVMnSpecies
        vLastID = vLastID + DVMSpecData(iSpec)%nVar
        MacroVal(:) = RefState_FV(:,iSpec,BCState)
        IF (BCType.EQ.24) MacroVal(5) = MacroVal(5)+Face_xGP(1,SideID)*BCTempGrad
        ! IF (BCType.EQ.25) MacroVal(5) = MacroVal(5)+(9.-Face_xGP(1,SideID))*2.*BCTempGrad
        CALL MaxwellDistribution(MacroVal,UPrim_boundary(vFirstID:vLastID),iSpec)
        CALL MaxwellScatteringDVM(iSpec,UPrim_boundary(vFirstID:vLastID),UPrim_master(vFirstID:vLastID,SideID),NormVec(:,SideID),prefac, &
                MacroValInside(:,DVMnSpecTot),MacroValInside(1,iSpec),rho,Pr,ErelaxTrans,ErelaxRot(iSpec),ErelaxVib(iSpec))
        vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
      END DO
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(5) !constant static pressure+temperature inlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL MacroValuesFromDistribution(MacroValInside,UPrim_master(:,SideID),dt/2.,tau,1)
      vFirstID=1
      vLastID=0
      DO iSpec=1,DVMnSpecies
        vLastID = vLastID + DVMSpecData(iSpec)%nVar
        MacroValInside(1,iSpec)=RefState_FV(1,iSpec,BCState)
        MacroValInside(5,iSpec)=RefState_FV(5,iSpec,BCState)
        CALL MaxwellDistribution(MacroValInside(:,iSpec),UPrim_boundary(vFirstID:vLastID),iSpec)
        vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
      END DO
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(6) !constant static pressure outlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      CALL MacroValuesFromDistribution(MacroValInside,UPrim_master(:,SideID),dt/2.,tau,1)
      vFirstID=1
      vLastID=0
      DO iSpec=1,DVMnSpecies
        vLastID = vLastID + DVMSpecData(iSpec)%nVar
        MacroValInside(5,iSpec)=RefState_FV(5,iSpec,BCState)*RefState_FV(1,iSpec,BCState)/MacroValInside(1,iSpec) !to get the pressure given by refstate
        CALL MaxwellDistribution(MacroValInside(:,iSpec),UPrim_boundary(vFirstID:vLastID),iSpec)
        vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
      END DO
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(7) !open outlet
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      UPrim_boundary=UPrim_master(:,SideID)
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
    END DO

  CASE(8) !dirichlet zero
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      UPrim_boundary=0.0
      CALL Riemann(Flux(:,SideID),UPrim_master(:,SideID),UPrim_boundary,NormVec(:,SideID))
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
#if !(USE_HDG)
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
#if !(USE_HDG)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
#endif
END SUBROUTINE FinalizeBC


END MODULE MOD_GetBoundaryFlux
