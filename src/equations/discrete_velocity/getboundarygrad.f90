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

#if (USE_FV)
MODULE MOD_GetBoundaryGrad
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

PUBLIC:: GetBoundaryGrad
!===================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes the gradient at a boundary for Finite Volumes reconstruction (2nd order version).
!==================================================================================================================================
SUBROUTINE GetBoundaryGrad(SideID,gradU,gradUinside,UPrim_master,NormVec,Face_xGP,output)
! MODULES
USE MOD_PreProc
USE MOD_Globals          ,ONLY: Abort
USE MOD_Mesh_Vars        ,ONLY: BC
USE MOD_Mesh_Vars_FV     ,ONLY: BoundaryType_FV
USE MOD_Equation_FV      ,ONLY: ExactFunc_FV
USE MOD_Equation_Vars_FV ,ONLY: IniExactFunc_FV, RefState_FV
USE MOD_Equation_Vars_FV ,ONLY: DVMVeloDisc, DVMSpecData, DVMnSpecies, DVMDim, DVMMethod, DVMnMacro, BCTempGrad
USE MOD_TimeDisc_Vars    ,ONLY: dt, time
USE MOD_DistFunc         ,ONLY: MaxwellDistribution, MacroValuesFromDistribution, MaxwellScatteringDVM, MoleculeRelaxEnergy
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN):: SideID
REAL,INTENT(OUT)  :: gradU       (PP_nVar_FV)
REAL,INTENT(IN)   :: gradUinside (PP_nVar_FV)
REAL,INTENT(IN)   :: UPrim_master(PP_nVar_FV)
REAL,INTENT(IN)   :: NormVec (3)
REAL,INTENT(IN)   :: Face_xGP(3)
LOGICAL,INTENT(IN):: output
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BCType,BCState
REAL    :: UPrim_boundary(1:PP_nVar_FV), fplus(1:PP_nVar_FV)
REAL    :: MacroVal(DVMnMacro), tau, vnormal, prefac, MacroValInside(DVMnMacro,DVMnSpecies+1),rho,Pr,relaxFac
INTEGER :: iVel, jVel, kVel, upos, upos_sp, iSpec, vFirstID, vLastID, BCorient
REAL    :: Erot(DVMnSpecies+1), ErelaxTrans, ErelaxRot(DVMnSpecies)
REAL    :: Evib(DVMnSpecies+1), Erelaxvib(DVMnSpecies)
!==================================================================================================================================
BCType  = BoundaryType_FV(BC(SideID),BC_TYPE)
BCState = BoundaryType_FV(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  IF(BCState.EQ.0) THEN ! Determine the exact BC state
    CALL ExactFunc_FV(IniExactFunc_FV,time,Face_xGP,UPrim_boundary)
  ELSE
    vFirstID=1
    vLastID=0
    DO iSpec=1,DVMnSpecies
      vLastID = vLastID + DVMSpecData(iSpec)%nVar
      CALL MaxwellDistribution(RefState_FV(:,iSpec,BCState),UPrim_boundary(vFirstID:vLastID),iSpec)
      vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
    END DO
  END IF
  gradU = UPrim_master - UPrim_boundary


CASE(3) ! specular reflection
  UPrim_boundary(:)=UPrim_master(:)-gradUinside(:)
  IF (BCState.NE.0) CALL abort(__STAMP__,'DVM specular bc with moving wall not implemented')
  IF      (ABS(ABS(NormVec(1)) - 1.).LE.1e-6) THEN !x-perpendicular boundary
    BCorient=1
  ELSE IF (ABS(ABS(NormVec(2)) - 1.).LE.1e-6) THEN !y-perpendicular boundary
    BCorient=2
  ELSE IF (ABS(ABS(NormVec(3)) - 1.).LE.1e-6) THEN !z-perpendicular boundary
    BCorient=3
  ELSE
    CALL abort(__STAMP__,'Specular reflection only implemented for boundaries perpendicular to velocity grid')
  END IF
  vFirstID=0
  DO iSpec=1,DVMnSpecies
    ASSOCIATE(Sp => DVMSpecData(iSpec))
    IF(DVMVeloDisc(iSpec).EQ.2.AND.(Sp%VeloMin(BCorient)+Sp%VeloMax(BCorient)).NE.0.) THEN
      CALL abort(__STAMP__,'Specular reflection only implemented for zero-centered velocity grid')
    END IF
    DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
      upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      vnormal = Sp%Velos(iVel,1)*NormVec(1) + Sp%Velos(jVel,2)*NormVec(2) + Sp%Velos(kVel,3)*NormVec(3)
      SELECT CASE(BCorient)
      CASE(1) !x-perpendicular boundary
        upos_sp=(Sp%nVelos(1)+1-iVel)+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      CASE(2) !y-perpendicular boundary
        upos_sp=iVel+(Sp%nVelos(2)-jVel)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      CASE(3) !z-perpendicular boundary
        upos_sp=iVel+(jVel-1)*Sp%nVelos(1)+(Sp%nVelos(3)-kVel)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
      END SELECT
      IF (vnormal.LT.0.) THEN !inflow
        gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos_sp))
        IF (DVMDim.LT.3) THEN
          gradU(Sp%nVarReduced+upos) = 2.*(UPrim_master(Sp%nVarReduced+upos) - UPrim_boundary(Sp%nVarReduced+upos_sp))
        END IF
        IF (Sp%Xi_Rot.GT.0.) THEN
          gradU(Sp%nVarErotStart+upos)=2.*(UPrim_master(Sp%nVarErotStart+upos) - UPrim_boundary(Sp%nVarErotStart+upos_sp))
        END IF
        IF (Sp%T_Vib.GT.0.) THEN
          gradU(Sp%nVarEvibStart+upos)=2.*(UPrim_master(Sp%nVarEvibStart+upos) - UPrim_boundary(Sp%nVarEvibStart+upos_sp))
        END IF
      ELSE
        gradU(upos) = 2.*gradUinside(upos)
        IF (DVMDim.LT.3) THEN
          gradU(Sp%nVarReduced+upos) = 2.*gradUinside(Sp%nVarReduced+upos)
        END IF
        IF (Sp%Xi_Rot.GT.0.) THEN
          gradU(Sp%nVarErotStart+upos) = 2.*gradUinside(Sp%nVarErotStart+upos)
        END IF
        IF (Sp%T_Vib.GT.0.) THEN
          gradU(Sp%nVarEvibStart+upos) = 2.*gradUinside(Sp%nVarEvibStart+upos)
        END IF
      END IF
    END DO; END DO; END DO
    END ASSOCIATE
  END DO


CASE(4,24,25) ! diffusive order 2 (see Baranger et al. 2019, MCS)
  fplus(:)=UPrim_master(:)-gradUinside(:)
  IF (output) THEN
    CALL MacroValuesFromDistribution(MacroValInside,fplus,dt,tau,1,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot,Evib=Evib)
  ELSE
    CALL MacroValuesFromDistribution(MacroValInside,fplus,dt/2.,tau,2,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot,Evib=Evib)
  END IF
  CALL MoleculeRelaxEnergy(ErelaxTrans,ErelaxRot,ErelaxVib,MacroValInside(5,DVMnSpecies+1),Erot(1:DVMnSpecies),Evib(1:DVMnSpecies),Pr)
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
  vFirstID=1
  vLastID=0
  DO iSpec = 1,DVMnSpecies
    ASSOCIATE(Sp => DVMSpecData(iSpec))
    vLastID = vLastID + Sp%nVar
    MacroVal(:) = RefState_FV(:,iSpec,BCState)
    IF (BCType.EQ.24) MacroVal(5) = MacroVal(5)+Face_xGP(1)*BCTempGrad
    ! IF (BCType.EQ.25) MacroVal(5) = MacroVal(5)+(9.-Face_xGP(1))*2.*BCTempGrad
    CALL MaxwellDistribution(MacroVal,UPrim_boundary(vFirstID:vLastID),iSpec)
    CALL MaxwellScatteringDVM(iSpec,UPrim_boundary(vFirstID:vLastID),fplus(vFirstID:vLastID),NormVec,prefac, &
                              MacroValInside(:,DVMnSpecies+1),MacroValInside(1,iSpec),rho,Pr,ErelaxTrans,ErelaxRot(iSpec),ErelaxVib(iSpec))
    DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
      upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID-1
      vnormal = Sp%Velos(iVel,1)*NormVec(1) &
              + Sp%Velos(jVel,2)*NormVec(2) &
              + Sp%Velos(kVel,3)*NormVec(3)
      IF (vnormal.LT.0.) THEN
        gradU(upos) = 2.*(UPrim_master(upos) - UPrim_boundary(upos))
        IF (DVMDim.LT.3) THEN
          gradU(Sp%nVarReduced+upos) = 2.*(UPrim_master(Sp%nVarReduced+upos) - UPrim_boundary(Sp%nVarReduced+upos))
        END IF
        IF (Sp%Xi_Rot.GT.0.) THEN
          gradU(Sp%nVarErotStart+upos)=2.*(UPrim_master(Sp%nVarErotStart+upos) - UPrim_boundary(Sp%nVarErotStart+upos))
        END IF
        IF (Sp%T_Vib.GT.0.) THEN
          gradU(Sp%nVarEvibStart+upos)=2.*(UPrim_master(Sp%nVarEvibStart+upos) - UPrim_boundary(Sp%nVarEvibStart+upos))
        END IF
      ELSE
        gradU(upos) = 2.*gradUinside(upos)
        IF (DVMDim.LT.3) THEN
          gradU(Sp%nVarReduced+upos) = 2.*gradUinside(Sp%nVarReduced+upos)
        END IF
        IF (Sp%Xi_Rot.GT.0.) THEN
          gradU(Sp%nVarErotStart+upos) = 2.*gradUinside(Sp%nVarErotStart+upos)
        END IF
        IF (Sp%T_Vib.GT.0.) THEN
          gradU(Sp%nVarEvibStart+upos) = 2.*gradUinside(Sp%nVarEvibStart+upos)
        END IF
      END IF
    END DO; END DO; END DO
    vFirstID = vFirstID + Sp%nVar
    END ASSOCIATE
  END DO

CASE(14) ! diffusive order 1
  IF (output) THEN
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt,tau,1,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot,Evib=Evib)
  ELSE
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt/2.,tau,2,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot,Evib=Evib)
  END IF
  CALL MoleculeRelaxEnergy(ErelaxTrans,ErelaxRot,ErelaxVib,MacroValInside(5,DVMnSpecies+1),Erot(1:DVMnSpecies),Evib(1:DVMnSpecies),Pr)
  IF (dt.EQ.0.) THEN
    prefac = 1.
  ELSE
    SELECT CASE(DVMMethod)
    CASE(1)
      prefac = 2.*tau*(1.-EXP(-dt/tau/2.))/dt ! f from f2~
    CASE(2)
      prefac = 2.*tau/(2.*tau+dt/2.)
    END SELECT
  END IF
  vFirstID=1
  vLastID=0
  DO iSPec = 1,DVMnSpecies
    vLastID = vLastID + DVMSpecData(iSpec)%nVar
    MacroVal(:) = RefState_FV(:,iSpec,BCState)
    CALL MaxwellDistribution(MacroVal,UPrim_boundary(vFirstID:vLastID),iSpec)
    CALL MaxwellScatteringDVM(iSpec,UPrim_boundary(vFirstID:vLastID),UPrim_master(vFirstID:vLastID),NormVec,prefac, &
                              MacroValInside(:,DVMnSpecies+1),MacroValInside(1,iSpec),rho,Pr,ErelaxTrans,ErelaxRot(iSpec),ErelaxVib(iSpec))
    vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
  END DO
  gradU = 2.*(UPrim_master-UPrim_boundary)

CASE(5) ! constant pressure and temperature
  IF (output) THEN
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt,tau,1)
  ELSE
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt/2.,tau,2)
  END IF
  vFirstID=1
  vLastID=0
  DO iSpec = 1,DVMnSpecies
    vLastID = vLastID + DVMSpecData(iSpec)%nVar
    MacroValInside(1,iSpec)=RefState_FV(1,iSpec,BCState)
    MacroValInside(5,iSpec)=RefState_FV(5,iSpec,BCState)
    CALL MaxwellDistribution(MacroValInside(:,iSpec),UPrim_boundary(vFirstID:vLastID),iSpec)
    vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
  END DO
  gradU = UPrim_master - UPrim_boundary

CASE(6) ! constant pressure
  IF (output) THEN
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt,tau,1)
  ELSE
    CALL MacroValuesFromDistribution(MacroValInside,UPrim_master,dt/2.,tau,2)
  END IF
  vFirstID=1
  vLastID=0
  DO iSpec = 1,DVMnSpecies
    vLastID = vLastID + DVMSpecData(iSpec)%nVar
    CALL MacroValuesFromDistribution(MacroVal,UPrim_master(vFirstID:vLastID),dt/2.,tau,1)
    MacroValInside(5,iSpec)=RefState_FV(5,iSpec,BCState)*RefState_FV(1,iSpec,BCState)/MacroValInside(1,iSpec) !to get the pressure given by refstate
    CALL MaxwellDistribution(MacroValInside(:,iSpec),UPrim_boundary(vFirstID:vLastID),iSpec)
    vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
  END DO
  gradU = UPrim_master - UPrim_boundary

CASE(7) ! open outlet
  gradU = 0.


CASE(8) ! dirichlet zero
  gradU = UPrim_master

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,&
        'no BC defined in DVM/getboundarygrad.f90!')
END SELECT ! BCType

END SUBROUTINE GetBoundaryGrad

END MODULE MOD_GetBoundaryGrad
#endif /*USE_FV*/
