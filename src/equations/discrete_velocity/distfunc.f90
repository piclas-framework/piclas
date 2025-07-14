!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"

!==================================================================================================================================
!> Routines providing distribution function management for the exponential integration method with discrete velocities
!==================================================================================================================================
MODULE MOD_DistFunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: MacroValuesFromDistribution, MaxwellDistribution, GradDistribution
PUBLIC:: TargetDistribution, MoleculeRelaxEnergy
PUBLIC:: MaxwellScatteringDVM, RescaleU, RescaleInit, ForceStep, IntegrateFluxValues
!==================================================================================================================================

CONTAINS

SUBROUTINE MacroValuesFromDistribution(MacroVal,U,tDeriv,tau,tilde,charge,MassDensity,PrandtlNumber,Erot)
!===================================================================================================================================
! Calculates the moments of the distribution function
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMMethod, DVMBGKModel, DVMDim, DVMnSpecies, DVMColl, DVMnMacro
USE MOD_PreProc
USE MOD_Globals                  ,ONLY: abort
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: U(PP_nVar_FV), tDeriv
INTEGER,INTENT(IN)              :: tilde
REAL, INTENT(OUT)               :: MacroVal(DVMnMacro,DVMnSpecies+1), tau
REAL, INTENT(OUT), OPTIONAL     :: charge, MassDensity, PrandtlNumber, Erot(DVMnSpecies+1)
! REAL, INTENT(OUT), OPTIONAL     :: skewness(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rhoTotal, dens(DVMnSpecies+1), rhoU(3,DVMnSpecies+1)
REAL                            :: densEtr(DVMnSpecies+1), densErot(DVMnSpecies+1)
REAL                            :: uVelo(3,DVMnSpecies+1), cVel(3,DVMnSpecies+1), cMag2(DVMnSpecies+1)
REAL                            :: mu(DVMnSpecies+1), thermalcond(DVMnSpecies+1), cP
REAL                            :: PressTens(6,DVMnSpecies+1), Heatflux(3,DVMnSpecies+1)
REAL                            :: weight, prefac, Phi, Prandtl, PrandtlCorr1, PrandtlCorr2, relaxFac
INTEGER                         :: iVel,jVel,kVel,upos,iSpec,jSpec,vFirstID,total
!===================================================================================================================================
rhoTotal = 0.
tau = 0.
dens = 0.
rhoU = 0.
densEtr = 0.
densErot = 0.
IF (PRESENT(Erot)) Erot = 0.
PressTens = 0.
Heatflux = 0.
MacroVal = 0.
cP = 0.
mu = 0.
thermalcond = 0.
PrandtlCorr1 = 0.
PrandtlCorr2 = 0.
total = DVMnSpecies+1
IF (PRESENT(charge)) charge = 0.

vFirstID = 0
DO iSpec=1, DVMnSpecies
  ASSOCIATE(Sp    => DVMSpecData(iSpec))

  ! IF (PRESENT(skewness)) skewness = 0.
  DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
    upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
    weight = Sp%Weights(iVel,1)*Sp%Weights(jVel,2)*Sp%Weights(kVel,3)
    dens(iSpec) = dens(iSpec) + weight*U(upos)
    rhoU(1,iSpec) = rhoU(1,iSpec) + weight*Sp%Mass*Sp%Velos(iVel,1)*U(upos)
    rhoU(2,iSpec) = rhoU(2,iSpec) + weight*Sp%Mass*Sp%Velos(jVel,2)*U(upos)
    rhoU(3,iSpec) = rhoU(3,iSpec) + weight*Sp%Mass*Sp%Velos(kVel,3)*U(upos)
    IF (DVMDim.LT.3) THEN
      densEtr(iSpec) = densEtr(iSpec) + weight*0.5*Sp%Mass*((Sp%Velos(iVel,1)**2.+Sp%Velos(jVel,2)**2.+Sp%Velos(kVel,3)**2.)*U(upos) &
                                            +U(Sp%nVarReduced+upos))
    ELSE
      densEtr(iSpec) = densEtr(iSpec) + weight*0.5*Sp%Mass*(Sp%Velos(iVel,1)**2.+Sp%Velos(jVel,2)**2.+Sp%Velos(kVel,3)**2.)*U(upos)
    END IF
    IF (Sp%Xi_Rot.GT.0) THEN
      densErot(iSpec) = densErot(iSpec) + weight*U(Sp%nVarErotStart+upos)
    END IF
  END DO; END DO; END DO
  uVelo(1:3,iSpec) = rhoU(1:3,iSpec)/dens(iSpec)/Sp%Mass

  MacroVal(1,iSpec) = dens(iSpec)
  MacroVal(2:4,iSpec) = uVelo(1:3,iSpec)
  MacroVal(5,iSpec) = (densEtr(iSpec) - 0.5*(DOT_PRODUCT(rhoU(:,iSpec),rhoU(:,iSpec)))/dens(iSpec)/Sp%Mass) &
                      *2./3./dens(iSpec)/BoltzmannConst
  IF (.NOT.(PRESENT(charge)).AND.MacroVal(5,iSpec).LE.0) THEN
    CALL abort(__STAMP__,'DVM negative temperature! Species n°',IntInfoOpt=iSpec)
  ELSE IF (.NOT.(PRESENT(charge))) THEN
    mu(iSpec) = Sp%mu_Ref*(MacroVal(5,iSpec)/Sp%T_Ref)**(Sp%omegaVHS+0.5)
    IF((Sp%InterID.EQ.2).OR.(Sp%InterID.EQ.20)) THEN ! inner DOF
      ! Istomin et. al., "Eucken correction in high-temperature gases with electronic excitation", J. Chem. Phys. 140,
      ! 184311 (2014)
      thermalcond(iSpec) = 0.25 * (15. + 2. * Sp%Xi_Rot * 1.328) * mu(iSpec) * BoltzmannConst / Sp%Mass
    ELSE ! atoms
      thermalcond(iSpec) = 0.25 * 15. * mu(iSpec) * BoltzmannConst / Sp%Mass
    END IF
  END IF

  rhoTotal = rhoTotal + Sp%Mass*dens(iSpec)
  dens(total) = dens(total) + dens(iSpec)
  rhoU(1:3,total) = rhoU(1:3,total) + rhoU(1:3,iSpec)
  densEtr(total) = densEtr(total) + densEtr(iSpec)
  densErot(total) = densErot(total) + densErot(iSpec)
  IF (PRESENT(charge)) charge = charge + Sp%Charge*dens(iSpec)

  vFirstID = vFirstID + Sp%nVar
  END ASSOCIATE
END DO
IF (PRESENT(Erot)) Erot(:) = densErot(:)/dens(:)
IF (PRESENT(charge)) RETURN

uVelo(1:3,total) = rhoU(1:3,total)/rhoTotal
MacroVal(1,total) = dens(total)
IF (PRESENT(MassDensity)) MassDensity = rhoTotal
MacroVal(2:4,total) = uVelo(1:3,total)
MacroVal(5,total) = (densEtr(total) - 0.5*(DOT_PRODUCT(rhoU(:,total),rhoU(:,total)))/rhoTotal)*2./3./dens(total)/BoltzmannConst
IF (MacroVal(5,total).LE.0) CALL abort(__STAMP__,'DVM negative total temperature!')

vFirstID = 0
DO iSpec=1, DVMnSpecies
  ASSOCIATE(Sp    => DVMSpecData(iSpec))

  DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
    upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
    weight = Sp%Weights(iVel,1)*Sp%Weights(jVel,2)*Sp%Weights(kVel,3)*Sp%Mass

    cVel(1,iSpec) = Sp%Velos(iVel,1) - uVelo(1,iSpec)
    cVel(2,iSpec) = Sp%Velos(jVel,2) - uVelo(2,iSpec)
    cVel(3,iSpec) = Sp%Velos(kVel,3) - uVelo(3,iSpec)
    cMag2(iSpec)  = DOT_PRODUCT(cVel(:,iSpec),cVel(:,iSpec))

    cVel(1,total) = Sp%Velos(iVel,1) - uVelo(1,total)
    cVel(2,total) = Sp%Velos(jVel,2) - uVelo(2,total)
    cVel(3,total) = Sp%Velos(kVel,3) - uVelo(3,total)
    cMag2(total)  = DOT_PRODUCT(cVel(:,total),cVel(:,total))

    PressTens(1,iSpec) = PressTens(1,iSpec) + weight*cVel(1,iSpec)*cVel(1,iSpec)*U(upos)
    PressTens(1,total) = PressTens(1,total) + weight*cVel(1,total)*cVel(1,total)*U(upos)
    IF (DVMDim.LT.3) THEN
      IF (DVMDim.EQ.1) THEN
        PressTens(2,iSpec) = PressTens(2,iSpec) + weight*0.5*U(Sp%nVarReduced+upos)
        PressTens(3,iSpec) = PressTens(3,iSpec) + weight*0.5*U(Sp%nVarReduced+upos)
        PressTens(2,total) = PressTens(2,total) + weight*0.5*U(Sp%nVarReduced+upos)
        PressTens(3,total) = PressTens(3,total) + weight*0.5*U(Sp%nVarReduced+upos)
      ELSE ! DVMDim.EQ.2
        PressTens(4,iSpec) = PressTens(4,iSpec) + weight*cVel(1,iSpec)*cVel(2,iSpec)*U(upos)
        PressTens(2,iSpec) = PressTens(2,iSpec) + weight*cVel(2,iSpec)*cVel(2,iSpec)*U(upos)
        PressTens(3,iSpec) = PressTens(3,iSpec) + weight*U(Sp%nVarReduced+upos)
        PressTens(4,total) = PressTens(4,total) + weight*cVel(1,total)*cVel(2,total)*U(upos)
        PressTens(2,total) = PressTens(2,total) + weight*cVel(2,total)*cVel(2,total)*U(upos)
        PressTens(3,total) = PressTens(3,total) + weight*U(Sp%nVarReduced+upos)
      END IF
      Heatflux(1,iSpec) = Heatflux(1,iSpec) + weight*0.5*cVel(1,iSpec)*(U(upos)*cMag2(iSpec)+U(Sp%nVarReduced+upos))
      Heatflux(2,iSpec) = Heatflux(2,iSpec) + weight*0.5*cVel(2,iSpec)*(U(upos)*cMag2(iSpec)+U(Sp%nVarReduced+upos))
      Heatflux(3,iSpec) = Heatflux(3,iSpec) + weight*0.5*cVel(3,iSpec)*(U(upos)*cMag2(iSpec)+U(Sp%nVarReduced+upos))
      Heatflux(1,total) = Heatflux(1,total) + weight*0.5*cVel(1,total)*(U(upos)*cMag2(total)+U(Sp%nVarReduced+upos))
      Heatflux(2,total) = Heatflux(2,total) + weight*0.5*cVel(2,total)*(U(upos)*cMag2(total)+U(Sp%nVarReduced+upos))
      Heatflux(3,total) = Heatflux(3,total) + weight*0.5*cVel(3,total)*(U(upos)*cMag2(total)+U(Sp%nVarReduced+upos))
    ELSE
      PressTens(2,iSpec) = PressTens(2,iSpec) + weight*cVel(2,iSpec)*cVel(2,iSpec)*U(upos)
      PressTens(3,iSpec) = PressTens(3,iSpec) + weight*cVel(3,iSpec)*cVel(3,iSpec)*U(upos)
      PressTens(4,iSpec) = PressTens(4,iSpec) + weight*cVel(1,iSpec)*cVel(2,iSpec)*U(upos)
      PressTens(5,iSpec) = PressTens(5,iSpec) + weight*cVel(1,iSpec)*cVel(3,iSpec)*U(upos)
      PressTens(6,iSpec) = PressTens(6,iSpec) + weight*cVel(2,iSpec)*cVel(3,iSpec)*U(upos)
      Heatflux(1:3,iSpec) = Heatflux(1:3,iSpec) + weight*0.5*cVel(1:3,iSpec)*(U(upos)*cMag2(iSpec))
      PressTens(2,total) = PressTens(2,total) + weight*cVel(2,total)*cVel(2,total)*U(upos)
      PressTens(3,total) = PressTens(3,total) + weight*cVel(3,total)*cVel(3,total)*U(upos)
      PressTens(4,total) = PressTens(4,total) + weight*cVel(1,total)*cVel(2,total)*U(upos)
      PressTens(5,total) = PressTens(5,total) + weight*cVel(1,total)*cVel(3,total)*U(upos)
      PressTens(6,total) = PressTens(6,total) + weight*cVel(2,total)*cVel(3,total)*U(upos)
      Heatflux(1:3,total) = Heatflux(1:3,total) + weight*0.5*cVel(1:3,total)*(U(upos)*cMag2(total))
      ! IF (PRESENT(skewness)) THEN
      !   skewness(1) = skewness(1) + weight*cVel(1)*(U(upos)*cVel(1)*cVel(1))
      !   skewness(2) = skewness(2) + weight*cVel(2)*(U(upos)*cVel(2)*cVel(2))
      !   skewness(3) = skewness(3) + weight*cVel(3)*(U(upos)*cVel(3)*cVel(3))
      ! END IF
    END IF
  END DO; END DO; END DO

  Macroval(6:11,iSpec)  = PressTens(1:6,iSpec)
  MacroVal(12:14,iSpec) = Heatflux(1:3,iSpec)

  IF (DVMnSpecies.GT.1) THEN
    ! Wilke's mixing rules
    Phi = 0.
    DO jSpec=1, DVMnSpecies
        Phi = Phi + dens(jSpec) * (1.+SQRT(mu(iSpec)/mu(jSpec)) &
        * (DVMSpecData(jSpec)%Mass/Sp%Mass)**(0.25) )**(2.0) &
        / ( SQRT(8.0 * (1.0 + Sp%Mass/DVMSpecData(jSpec)%Mass)) )
    END DO
    mu(total) = mu(total) + dens(iSpec)*mu(iSpec)/Phi
    thermalcond(total) = thermalcond(total) + dens(iSpec)*thermalcond(iSpec)/Phi
    cP = cP + ((5.+Sp%Xi_Rot)/2.) * BoltzmannConst * dens(iSpec)/rhoTotal
    PrandtlCorr1 = PrandtlCorr1 + (5.+Sp%Xi_Rot)*dens(iSpec)/Sp%Mass
    PrandtlCorr2 = PrandtlCorr2 + (5.+Sp%Xi_Rot)*dens(iSpec)
  ELSE
    mu(total) = mu(1)
    thermalcond(total) = thermalcond(1)
  END IF
  vFirstID = vFirstID + Sp%nVar
  END ASSOCIATE
END DO

IF (DVMnSpecies.GT.1) THEN
  Prandtl = cP*mu(total)/thermalcond(total)*PrandtlCorr1*rhoTotal/PrandtlCorr2/dens(total) !Pr = alpha * cP * mu / K
ELSE
  Prandtl = 2.*(DVMSpecData(1)%Xi_Rot + 5.)/(2.*DVMSpecData(1)%Xi_Rot + 15.)
END IF
tau = mu(total)/BoltzmannConst/dens(total)/MacroVal(5,total)
IF (DVMBGKModel.EQ.1) tau = tau/Prandtl !ESBGK
IF (DVMBGKModel.EQ.6) tau = tau*2./3. !Grad13BGK

IF (PRESENT(PrandtlNumber)) PrandtlNumber = Prandtl

Macroval(6:11,total)  = PressTens(1:6,total)
MacroVal(12:14,total) = Heatflux(1:3,total)

IF (DVMColl.AND.DVMMethod.GT.0.AND.tDeriv.GT.0.) THEN
  IF(tilde.EQ.1) THEN ! higher moments from f~
    SELECT CASE(DVMMethod)
    CASE(1) !EDDVM
      relaxFac = tDeriv/tau
      IF (CHECKEXP(relaxFac)) THEN
        prefac = (1.-EXP(-relaxFac))/relaxFac
      ELSE
        prefac = 0.
      END IF
      SELECT CASE(DVMBGKModel)
      CASE(1) !ESBGK
        Macroval(6,:)   = (Macroval(6,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        Macroval(7,:)   = (Macroval(7,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        Macroval(8,:)   = (Macroval(8,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        MacroVal(9:11,:)  = Macroval(9:11,:)*prefac/(1./Prandtl+prefac*(1.-1./Prandtl))
        MacroVal(12:14,:) = MacroVal(12:14,:)*prefac
      CASE(2,5) !Shakhov/SN
        MacroVal(6,:)   = Macroval(6,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        MacroVal(7,:)   = Macroval(7,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        MacroVal(8,:)   = Macroval(8,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        MacroVal(9:11,:)  = MacroVal(9:11,:)*prefac
        MacroVal(12:14,:) = MacroVal(12:14,:)*prefac/(Prandtl+prefac*(1-Prandtl))
      CASE(3,4) !Maxwell
        MacroVal(6,:)   = Macroval(6,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        MacroVal(7,:)   = Macroval(7,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        MacroVal(8,:)   = Macroval(8,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst
        Macroval(9:14,:)  = prefac*Macroval(9:14,:) ! non-eq moments should be zero anyway
      CASE(6) !Double moment distributions
        Macroval(6,:)   = (Macroval(6,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst*2./3.) &
        /(2./3.+prefac/3.)
        Macroval(7,:)   = (Macroval(7,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst*2./3.) &
        /(2./3.+prefac/3.)
        Macroval(8,:)   = (Macroval(8,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst*2./3.) &
        /(2./3.+prefac/3.)
        MacroVal(9:11,:)  = Macroval(9:11,:)*prefac/(2./3.+prefac/3.)
        MacroVal(12:14,:) = MacroVal(12:14,:)*prefac/(2.*Prandtl/3.+prefac*(1.-2.*Prandtl/3.))
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM-BGKCollModel does not exist')
      END SELECT
    CASE(2) !DUGKS
      SELECT CASE(DVMBGKModel)
      CASE(1) !ESBGK
        prefac = (2.*tau)/(2.*tau+tDeriv)
        Macroval(6,:)   = (Macroval(6,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        Macroval(7,:)   = (Macroval(7,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        Macroval(8,:)   = (Macroval(8,:)*prefac+(1.-prefac)*MacroVal(5,:)*dens(:)*BoltzmannConst/Prandtl) &
                          /(1./Prandtl+prefac*(1.-1./Prandtl))
        MacroVal(9:11,:)  = Macroval(9:11,:)*prefac/(1./Prandtl+prefac*(1.-1./Prandtl))
        MacroVal(12:14,:) = MacroVal(12:14,:)*prefac
      CASE(2,5) !Shakhov/SN
        Macroval(6,:)   = Macroval(6,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(7,:)   = Macroval(7,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(8,:)   = Macroval(8,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(9:11,:)  = Macroval(9:11,:)*2.*tau/(2.*tau+tDeriv)
        MacroVal(12:14,:) = MacroVal(12:14,:)*2.*tau/(2.*tau+tDeriv*Prandtl)
      CASE(3,4) !Maxwell
        Macroval(6,:)   = Macroval(6,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(7,:)   = Macroval(7,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(8,:)   = Macroval(8,:)*2.*tau/(2.*tau+tDeriv) + MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+tDeriv)
        Macroval(9:14,:)  = Macroval(9:14,:)*2.*tau/(2.*tau+tDeriv) ! non-eq moments should be zero anyway
      CASE(6) !Double moment distributions
        Macroval(6,:)   = Macroval(6,:)*2.*tau/(2.*tau+2.*tDeriv/3.) &
                        + 2.*MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+2.*tDeriv/3.)/3.
        Macroval(7,:)   = Macroval(7,:)*2.*tau/(2.*tau+2.*tDeriv/3.) &
                        + 2.*MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+2.*tDeriv/3.)/3.
        Macroval(8,:)   = Macroval(8,:)*2.*tau/(2.*tau+2.*tDeriv/3.) &
                        + 2.*MacroVal(5,:)*dens(:)*BoltzmannConst*tDeriv/(2.*tau+2.*tDeriv/3.)/3.
        Macroval(9:11,:)  = Macroval(9:11,:)*2.*tau/(2.*tau+2.*tDeriv/3.)
        MacroVal(12:14,:) = MacroVal(12:14,:)*2.*tau/(2.*tau+2.*tDeriv*Prandtl/3.)
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM-BGKCollModel does not exist')
      END SELECT
    CASE DEFAULT
      CALL abort(__STAMP__,'DVM-Method does not exist')
    END SELECT
  END IF
END IF

END SUBROUTINE MacroValuesFromDistribution

SUBROUTINE MoleculeRelaxEnergy(ErelaxTrans,ErelaxRot,TempTrans,Erot,Prandtl)!,Zrot
!===================================================================================================================================
! Calculate relaxation energies (translational and rotational) for molecules
! (Hild & Pfeiffer 2024, J. Comp. Phys. 514, 113226)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMnSpecies, DVMnInnerE
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: TempTrans, Erot(DVMnSpecies), Prandtl
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: ErelaxTrans,ErelaxRot(DVMnSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec
REAL                            :: EdiffRot
!===================================================================================================================================
ErelaxTrans = TempTrans*BoltzmannConst*3./2.
ErelaxRot = 0.

IF (DVMnInnerE.EQ.0) RETURN

DO iSpec = 1, DVMnSpecies
  ASSOCIATE(Sp => DVMSpecData(iSpec))
  IF (Sp%Xi_Rot.GT.0) THEN
    EdiffRot = 30.*(TempTrans*BoltzmannConst*Sp%Xi_Rot/2.-Erot(iSpec))/(4.-2.*Sp%omegaVHS)/(6.-2.*Sp%omegaVHS)/Prandtl/Sp%Z_Rot
    ! change this to include actual multispecies collfreq
    ErelaxRot(iSpec) = Erot(iSpec) + EdiffRot
    ErelaxTrans = ErelaxTrans - EdiffRot
  END IF
  END ASSOCIATE
END DO

END SUBROUTINE MoleculeRelaxEnergy

SUBROUTINE MaxwellDistributionCons(MacroVal,fMaxwell,iSpec,densSpec)
!===================================================================================================================================
! conservative maxwell distribution (cf Mieussens 2000)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fMaxwell(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, Temp, uVelo(DVMDim), vMag, weight, gM
REAL,DIMENSION(2+DVMDim)        :: alpha, psi, rhovec
REAL                            :: J(2+DVMDim,2+DVMDim), B(2+DVMDim)
INTEGER                         :: iVel,jVel,kVel, upos, countz, IPIV(2+DVMDim), info_dgesv
!===================================================================================================================================
IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
  ! molecules with rotational DOF
  CALL abort(__STAMP__,'Conservative DVM not implemented for species with rotational DOF')
END IF

rho = densSpec*DVMSpecData(iSpec)%Mass
uVelo(1:DVMDim) = MacroVal(2:1+DVMDim)
Temp = MacroVal(5)

! vector of conservative variables
rhovec(1) = rho
rhovec(2:1+DVMDim)=rho*uVelo(:)
rhovec(2+DVMDim)=rho*(FLOAT(DVMDim)*DVMSpecData(iSpec)%R_S*Temp+DOT_PRODUCT(uVelo,uVelo))/2.

alpha(1) = LOG(rho/(2.*Pi*DVMSpecData(iSpec)%R_S*Temp)**(DVMDim/2.))-DOT_PRODUCT(uVelo,uVelo)/2./DVMSpecData(iSpec)%R_S/Temp
alpha(2:1+DVMDim) = uVelo(1:DVMDim)/DVMSpecData(iSpec)%R_S/Temp
alpha(2+DVMDim) = -1/DVMSpecData(iSpec)%R_S/Temp

! init counter
countz=0

!Newton algorithm to find alpha so that <psi exp(alpha.psi)> = rhovec
DO WHILE (countz.LT.1000)
  countz=countz+1
  J = 0.

  DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
    weight = DVMSpecData(iSpec)%Weights(iVel,1)*DVMSpecData(iSpec)%Weights(jVel,2)*DVMSpecData(iSpec)%Weights(kVel,3)
    vMag = DVMSpecData(iSpec)%Velos(iVel,1)**2 + DVMSpecData(iSpec)%Velos(jVel,2)**2 + DVMSpecData(iSpec)%Velos(kVel,3)**2
    psi(1)=1
    psi(2)=DVMSpecData(iSpec)%Velos(iVel,1)
    IF (DVMDim.GT.1) psi(3)=DVMSpecData(iSpec)%Velos(jVel,2)
    IF (DVMDim.GT.2) psi(4)=DVMSpecData(iSpec)%Velos(kVel,3)
    psi(2+DVMDim) = vMag/2.
    gM = EXP(DOT_PRODUCT(alpha,psi))

    ! J: derivative of -B = <psi exp(alpha.psi)> - rhovec
    ! J = <psi x psi exp(alpha.psi)>
    J(:,1) = J(:,1)+weight*gM*psi(:)
    J(:,2) = J(:,2)+weight*gM*psi(2)*psi(:)
    IF (DVMDim.GT.1) J(:,3) = J(:,3)+weight*gM*psi(3)*psi(:)
    IF (DVMDim.GT.2) J(:,4) = J(:,4)+weight*gM*psi(4)*psi(:)
    J(:,2+DVMDim) = J(:,2+DVMDim)+weight*gM*psi(2+DVMDim)*psi(:)
  END DO; END DO; END DO

  B(:) = rhovec(:) - J(:,1)

  ! solve JX = B
  CALL DGESV(2+DVMDim,1,J,2+DVMDim,IPIV,B,2+DVMDim,info_dgesv)

  IF(info_dgesv.NE.0) CALL abort(__STAMP__,'Newton DGESV fail')
  IF (.NOT.(ANY(ABS(B(:)).GT.(1e-5*ABS(alpha)+1e-12)))) EXIT
  ! IF (.NOT.(ANY(ABS(B(:,1)).GT.(1e-6*ABS(alpha)+1e-16)))) EXIT

  ! update alpha with Newton increment X (stored in B)
  alpha = alpha + B
END DO
IF (countz.GE.1000) print*, 'Newton max iter reached'

alpha = alpha + B
DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  weight = DVMSpecData(iSpec)%Weights(iVel,1)*DVMSpecData(iSpec)%Weights(jVel,2)*DVMSpecData(iSpec)%Weights(kVel,3)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
  vMag = DVMSpecData(iSpec)%Velos(iVel,1)**2 + DVMSpecData(iSpec)%Velos(jVel,2)**2 + DVMSpecData(iSpec)%Velos(kVel,3)**2
  psi(1)=1
  psi(2)=DVMSpecData(iSpec)%Velos(iVel,1)
  IF (DVMDim.GT.1) psi(3)=DVMSpecData(iSpec)%Velos(jVel,2)
  IF (DVMDim.GT.2) psi(4)=DVMSpecData(iSpec)%Velos(kVel,3)
  psi(2+DVMDim) = vMag/2.
  fMaxwell(upos) = EXP(DOT_PRODUCT(alpha,psi))/DVMSpecData(iSpec)%Mass
  IF (DVMDim.LT.3) THEN
    fMaxwell(DVMSpecData(iSpec)%nVarReduced+upos) = fMaxwell(upos)*DVMSpecData(iSpec)%R_S*Temp*(3.-DVMDim)
  END IF
END DO; END DO; END DO

END SUBROUTINE MaxwellDistributionCons

SUBROUTINE MaxwellDistribution(MacroVal,fMaxwell,iSpec,densSpec,ErelaxRot)
!===================================================================================================================================
! Maxwell distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_PreProc
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fMaxwell(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro)
REAL, INTENT(IN), OPTIONAL       :: densSpec,ErelaxRot
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: dens, Temp, uVelo(3), cVel(3), cMag2
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
IF (PRESENT(densSpec)) THEN
  dens = densSpec
ELSE
  dens = MacroVal(1)
END IF
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
  cVel(1) = DVMSpecData(iSpec)%Velos(iVel,1) - uVelo(1)
  cVel(2) = DVMSpecData(iSpec)%Velos(jVel,2) - uVelo(2)
  cVel(3) = DVMSpecData(iSpec)%Velos(kVel,3) - uVelo(3)
  cMag2 = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  fMaxwell(upos) = dens/((2.*Pi*DVMSpecData(iSpec)%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpecData(iSpec)%R_S*Temp))
  IF (DVMDim.LT.3) THEN
    fMaxwell(DVMSpecData(iSpec)%nVarReduced+upos) = fMaxwell(upos)*DVMSpecData(iSpec)%R_S*Temp*(3.-DVMDim)
  END IF
  IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
    ! molecules with rotational DOF
    IF (PRESENT(ErelaxRot)) THEN
      fMaxwell(DVMSpecData(iSpec)%nVarErotStart+upos) = fMaxwell(upos)*ErelaxRot
    ELSE
      fMaxwell(DVMSpecData(iSpec)%nVarErotStart+upos) = fMaxwell(upos)*Temp*BoltzmannConst*DVMSpecData(iSpec)%Xi_Rot/2.
    END IF
  END IF
END DO; END DO; END DO

END SUBROUTINE MaxwellDistribution

SUBROUTINE ShakhovDistribution(MacroVal,fShakhov,iSpec,densSpec,rhoTotal,Prandtl)
!===================================================================================================================================
! Shakhov distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_PreProc
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fShakhov(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec,rhoTotal,Prandtl
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: dens, Temp, uVelo(3), cVel(3), cMag2, gM, q(3), ShakhFac1, ShakhFac2, densTotal
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
dens = densSpec
densTotal = MacroVal(1)
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)
q(1:3) = MacroVal(12:14)

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
  cVel(1) = DVMSpecData(iSpec)%Velos(iVel,1) - uVelo(1)
  cVel(2) = DVMSpecData(iSpec)%Velos(jVel,2) - uVelo(2)
  cVel(3) = DVMSpecData(iSpec)%Velos(kVel,3) - uVelo(3)
  cMag2 = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  gM = dens/((2.*Pi*DVMSpecData(iSpec)%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpecData(iSpec)%R_S*Temp))
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*rhoTotal*(DVMSpecData(iSpec)%R_S*Temp)**2)
  ShakhFac2 = cMag2/(DVMSpecData(iSpec)%R_S*Temp)
  fShakhov(upos) = gM*(1.+(1.-Prandtl)*ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMDim.LT.3) THEN
    fShakhov(DVMSpecData(iSpec)%nVarReduced+upos) = gM*DVMSpecData(iSpec)%R_S*Temp*(DVMSpecData(iSpec)%Xi_Rot+3.-DVMDim &
                                + (1.-Prandtl)*ShakhFac1*((ShakhFac2-DVMDim)*(DVMSpecData(iSpec)%Xi_Rot+3.-DVMDim)&
                                - 2.*DVMSpecData(iSpec)%Xi_Rot))
  END IF
END DO; END DO; END DO

END SUBROUTINE ShakhovDistribution

SUBROUTINE ESBGKDistribution(MacroVal,fESBGK,iSpec,densSpec,Prandtl,ErelaxTrans,ErelaxRot)
!===================================================================================================================================
! ESBGK distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_Basis                    ,ONLY: INV33
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fESBGK(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec,Prandtl,ErelaxTrans,ErelaxRot
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: dens,Temp,uVelo(3),cVel(3),pressTens(3,3),pressProduct,ilambda(3,3),ldet,hfac,densTotal,TtransRel
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
dens             = densSpec
densTotal        = MacroVal(1)
uVelo(1:3)       = MacroVal(2:4)
Temp             = MacroVal(5)
pressTens(1,1)   = MacroVal(6)
pressTens(2,2)   = MacroVal(7)
pressTens(3,3)   = MacroVal(8)
pressTens(1,2:3) = MacroVal(9:10)
pressTens(2:3,1) = MacroVal(9:10)
pressTens(2,3)   = MacroVal(11)
pressTens(3,2)   = MacroVal(11)
TtransRel = 2.*ErelaxTrans/BoltzmannConst/3.

pressTens = (1.-1./Prandtl)*pressTens*TtransRel/Temp/densTotal/DVMSpecData(iSpec)%Mass
pressTens(1,1) = pressTens(1,1)+DVMSpecData(iSpec)%R_S*TtransRel/Prandtl
pressTens(2,2) = pressTens(2,2)+DVMSpecData(iSpec)%R_S*TtransRel/Prandtl
pressTens(3,3) = pressTens(3,3)+DVMSpecData(iSpec)%R_S*TtransRel/Prandtl

CALL INV33(pressTens,ilambda,ldet)
IF (ldet.LE.0.) THEN
  CALL abort(__STAMP__,'ESBGK matrix not positive-definite')
ELSE
  ! determinant from reduced pressure tensor (size D*D)
  IF (DVMDim.LE.2) ldet = ldet/pressTens(3,3)
  IF (DVMDim.LE.1) ldet = ldet/pressTens(2,2)
END IF

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)

  cVel(1) = DVMSpecData(iSpec)%Velos(iVel,1) - uVelo(1)
  cVel(2) = DVMSpecData(iSpec)%Velos(jVel,2) - uVelo(2)
  cVel(3) = DVMSpecData(iSpec)%Velos(kVel,3) - uVelo(3)

  pressProduct = cVel(1)*DOT_PRODUCT(ilambda(:,1),cVel) &
                + cVel(2)*DOT_PRODUCT(ilambda(:,2),cVel) &
                + cVel(3)*DOT_PRODUCT(ilambda(:,3),cVel)

  fESBGK(upos) = dens/sqrt(ldet*(2*Pi)**DVMDim)*EXP(-pressProduct/2.)
  IF (DVMDim.LT.3) THEN
    hfac = 0.
    IF (DVMDim.LE.2) hfac = hfac + pressTens(3,3)
    IF (DVMDim.LE.1) hfac = hfac + pressTens(2,2)
    fESBGK(DVMSpecData(iSpec)%nVarReduced+upos) = fESBGK(upos)*hfac
  END IF
  IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
    ! molecules with rotational DOF
    fESBGK(DVMSpecData(iSpec)%nVarErotStart+upos) = fESBGK(upos)*ErelaxRot
  END IF
END DO; END DO; END DO

END SUBROUTINE ESBGKDistribution

SUBROUTINE GradDistribution(MacroVal,fGrad,iSpec,densSpec,rhoTotal,ErelaxRot)
!===================================================================================================================================
! Grad 13 moments distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fGrad(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro)
REAL, INTENT(IN), OPTIONAL       :: densSpec,rhoTotal,ErelaxRot
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: dens,Temp,uVelo(3),cVel(3),cMag2,gM,q(3),pressTens(3,3),ShakhFac1,ShakhFac2,pressFac,pressProduct
REAL                            :: pressProduct2, densTotal, rho
INTEGER                         :: iVel,jVel,kVel,upos
!===================================================================================================================================
IF (PRESENT(densSpec).AND.PRESENT(rhoTotal)) THEN
  dens = densSpec
  rho = rhoTotal
ELSE
  dens = MacroVal(1)
  rho = MacroVal(1)*DVMSpecData(iSpec)%Mass
END IF
densTotal        = MacroVal(1)
uVelo(1:3)       = MacroVal(2:4)
Temp             = MacroVal(5)
pressTens(1,1)   = MacroVal(6)
pressTens(2,2)   = MacroVal(7)
pressTens(3,3)   = MacroVal(8)
pressTens(1,2:3) = MacroVal(9:10)
pressTens(2:3,1) = MacroVal(9:10)
pressTens(2,3)   = MacroVal(11)
pressTens(3,2)   = MacroVal(11)
q(1:3)           = MacroVal(12:14)

! here the traceless pressure tensor is used (init with MacroVal(6:8)=0 for Tx=Ty=Tz)
IF (ABS(SUM(MacroVal(6:8))).GT.1.e-12*(densTotal*BoltzmannConst*Temp)) CALL abort(__STAMP__, &
                                  'Diagonal entries of the pressure tensor should add up to zero',0,SUM(MacroVal(6:8)))

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)

  cVel(1) = DVMSpecData(iSpec)%Velos(iVel,1) - uVelo(1)
  cVel(2) = DVMSpecData(iSpec)%Velos(jVel,2) - uVelo(2)
  cVel(3) = DVMSpecData(iSpec)%Velos(kVel,3) - uVelo(3)
  cMag2 = DOT_PRODUCT(cVel,cVel)

  pressProduct = cVel(1)*DOT_PRODUCT(pressTens(:,1),cVel) &
               + cVel(2)*DOT_PRODUCT(pressTens(:,2),cVel) &
               + cVel(3)*DOT_PRODUCT(pressTens(:,3),cVel)

  IF (DVMDim.LT.3) THEN
    pressProduct2 = pressProduct + (5.-DVMDim)*DVMSpecData(iSpec)%R_S*Temp*pressTens(3,3)
    pressProduct = pressProduct + (3.-DVMDim)*DVMSpecData(iSpec)%R_S*Temp*pressTens(3,3)
  END IF

  pressFac = rho*(DVMSpecData(iSpec)%R_S*Temp)**2
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*pressFac)
  ShakhFac2 = cMag2/(DVMSpecData(iSpec)%R_S*Temp)

  gM = dens/((2.*Pi*DVMSpecData(iSpec)%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpecData(iSpec)%R_S*Temp))
  fGrad(upos) = gM*(1.+0.5*pressProduct/pressFac+ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMDim.LT.3) THEN
    fGrad(DVMSpecData(iSpec)%nVarReduced+upos) = gM*DVMSpecData(iSpec)%R_S*Temp*(3.-DVMDim) &
                  *(1 + 0.5*pressProduct2/pressFac + ShakhFac1*(ShakhFac2-DVMDim))
  END IF
  IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
    ! molecules with rotational DOF
    IF (PRESENT(ErelaxRot)) THEN
      fGrad(DVMSpecData(iSpec)%nVarErotStart+upos) = fGrad(upos)*ErelaxRot
    ELSE
      fGrad(DVMSpecData(iSpec)%nVarErotStart+upos) = fGrad(upos)*Temp*BoltzmannConst*DVMSpecData(iSpec)%Xi_Rot/2.
    END IF
  END IF
END DO; END DO; END DO

END SUBROUTINE GradDistribution

SUBROUTINE GradDistributionPrandtl(MacroVal,fGrad,iSpec,densSpec,rhoTotal,Prandtl)
!===================================================================================================================================
! Grad 13 moments distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMnMacro
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars             ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fGrad(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec,rhoTotal,Prandtl
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: MacroValPrandtl(DVMnMacro)
!===================================================================================================================================
IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
  ! molecules with rotational DOF
  CALL abort(__STAMP__,'Grad13 BGK not implemented for species with rotational DOF')
END IF

MacroValPrandtl(1:5)   = MacroVal(1:5)

MacroValPrandtl(6)     = (MacroVal(6)-MacroVal(1)*BoltzmannConst*MacroVal(5))/3.
MacroValPrandtl(7)     = (MacroVal(7)-MacroVal(1)*BoltzmannConst*MacroVal(5))/3.
MacroValPrandtl(8)     = (MacroVal(8)-MacroVal(1)*BoltzmannConst*MacroVal(5))/3.
MacroValPrandtl(9:11)  = MacroVal(9:11)/3.

MacroValPrandtl(12:14) = (1.-2.*Prandtl/3.)*MacroVal(12:14)

CALL GradDistribution(MacroValPrandtl,fGrad,iSpec,densSpec,rhoTotal)

END SUBROUTINE GradDistributionPrandtl

SUBROUTINE SkewNormalDistribution(MacroVal,fSkew,iSpec,densSpec,rhoTotal,Prandtl)
!===================================================================================================================================
! Skew-normal distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMDim, Pi, DVMnMacro
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fSkew(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec,rhoTotal,Prandtl
! REAL, INTENT(IN),OPTIONAL        :: skewness(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: dens, Temp, uVelo(3), cVel(3), cMag2, q(3)
INTEGER                         :: iVel,jVel,kVel, upos
REAL                            :: skew(1:3), delta(1:3), alpha(1:3), ksi(1:3), omega(1:3), Phi(1:3), max_skew
!===================================================================================================================================
IF (DVMSpecData(iSpec)%InterID.EQ.2.OR.DVMSpecData(iSpec)%InterID.EQ.20) THEN
  ! molecules with rotational DOF
  CALL abort(__STAMP__,'Skew normal BGK not implemented for species with rotational DOF')
END IF

dens = densSpec
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)
q(1:3) = MacroVal(12:14)

max_skew = 0.5 * (4. - Pi) * (2. / (Pi - 2.)) ** 1.5
skew = (1-Prandtl)*2.*(q/rhoTotal)*(DVMSpecData(iSpec)%R_S*Temp)**(-3./2.)
! skew = (1-DVMSpecData(iSpec)%Prandtl)*(skewness/rho)*(DVMSpecData(iSpec)%R_S*Temp)**(-3./2.)

skew = SIGN(1.,skew)*MIN(ABS(skew),0.9*max_skew)

delta = (SIGN(1.,skew)*(2*ABS(skew)/(4-Pi))**(1./3.))/SQRT(2./Pi*(1+(2*ABS(skew)/(4-Pi))**(2./3.)))
alpha = delta/SQRT(1.-delta*delta)
omega = SQRT(DVMSpecData(iSpec)%R_S*Temp/(1.-2*delta*delta/Pi))
ksi = uVelo - SQRT(2/Pi)*omega*delta

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
  cVel(1) = (DVMSpecData(iSpec)%Velos(iVel,1) - ksi(1))/omega(1)
  cVel(2) = (DVMSpecData(iSpec)%Velos(jVel,2) - ksi(2))/omega(2)
  cVel(3) = (DVMSpecData(iSpec)%Velos(kVel,3) - ksi(3))/omega(3)
  cMag2 = DOT_PRODUCT(cVel,cVel)
  Phi = 1.+ERF(alpha*cVel/sqrt(2.))

  fSkew(upos) = dens*Phi(1)*Phi(2)*Phi(3)*EXP(-cMag2/2.)/PRODUCT(omega(1:DVMDim))/(2.*Pi)**(DVMDim/2.)

  IF (DVMDim.LT.3) THEN
    fSkew(DVMSpecData(iSpec)%nVarReduced+upos) = fSkew(upos)*DVMSpecData(iSpec)%R_S*Temp*(3.-DVMDim)
  END IF

END DO; END DO; END DO

END SUBROUTINE SkewNormalDistribution

SUBROUTINE TargetDistribution(MacroVal,fTarget,iSpec,densSpec,rho,Pr,ErelaxTrans,ErelaxRot)
!===================================================================================================================================
! Target distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMBGKModel,DVMSpecData,DVMnSpecies,DVMnMacro
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
REAL,INTENT(OUT)                 :: fTarget(DVMSpecData(iSpec)%nVar)
REAL, INTENT(IN)                 :: MacroVal(DVMnMacro),densSpec,rho,Pr
REAL, INTENT(IN)                 :: ErelaxTrans,ErelaxRot
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: MassDensity
!===================================================================================================================================

SELECT CASE(DVMBGKModel)
  CASE(1)
    CALL ESBGKDistribution(MacroVal,fTarget,iSpec,densSpec,Pr,ErelaxTrans,ErelaxRot)
  CASE(2)
    CALL ShakhovDistribution(MacroVal,fTarget,iSpec,densSpec,rho,Pr)
  CASE(3)
    CALL MaxwellDistribution(MacroVal,fTarget,iSpec,densSpec)
  CASE(4)
    CALL MaxwellDistributionCons(MacroVal,fTarget,iSpec,densSpec)
  CASE(5)
    CALL SkewNormalDistribution(MacroVal,fTarget,iSpec,densSpec,rho,Pr)
  CASE(6)
    CALL GradDistributionPrandtl(MacroVal,fTarget,iSpec,densSpec,rho,Pr)
  CASE DEFAULT
    CALL abort(__STAMP__,'DVM BGK Model not implemented.')
END SELECT

END SUBROUTINE TargetDistribution

SUBROUTINE MaxwellScatteringDVM(iSpec,fBoundary,U,NormVec,prefac,MacroVal,densSpec,MassDensity,PrandtlNumber,ErelaxTrans,ErelaxRot)
!===================================================================================================================================
! Gets accurate density for the half maxwellian at diffusive boundaries
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMnMacro, DVMMethod
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: iSpec
REAL,INTENT(INOUT)                   :: fBoundary(DVMSpecData(iSpec)%nVar)
REAL,INTENT(IN)                      :: U(DVMSpecData(iSpec)%nVar),NormVec(3),prefac,MacroVal(DVMnMacro),densSpec,MassDensity,PrandtlNumber
REAL,INTENT(IN)                      :: ErelaxTrans,ErelaxRot
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Sout, Sin, weight, vnormal
REAL                                 :: fRelaxed(DVMSpecData(iSpec)%nVar)
INTEGER                              :: iVel, jVel, kVel, upos
!===================================================================================================================================

IF (DVMMethod.GT.0) THEN
  CALL TargetDistribution(MacroVal, fRelaxed, iSpec, densSpec, MassDensity, PrandtlNumber,ErelaxTrans,ErelaxRot)
  fRelaxed = U*prefac + fRelaxed*(1.-prefac)
ELSE
  ! first order method
  fRelaxed = U
END IF

Sin = 0.
Sout = 0.

DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
  upos= iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)
  vnormal = DVMSpecData(iSpec)%Velos(iVel,1)*NormVec(1) &
          + DVMSpecData(iSpec)%Velos(jVel,2)*NormVec(2) &
          + DVMSpecData(iSpec)%Velos(kVel,3)*NormVec(3)
  weight = DVMSpecData(iSpec)%Weights(iVel,1)*DVMSpecData(iSpec)%Weights(jVel,2)*DVMSpecData(iSpec)%Weights(kVel,3)
  IF (vnormal.GT.0.) THEN !outflow
    Sout = Sout + weight*vnormal*fRelaxed(upos)
  ELSE !inflow
    Sin = Sin - weight*vnormal*fBoundary(upos)
  END IF
END DO; END DO; END DO

fBoundary = fBoundary * (Sout/Sin)
! no additional rescaling needed because it is an equilibrium distribution

END SUBROUTINE MaxwellScatteringDVM

SUBROUTINE RescaleU(tilde,tDeriv)
!===================================================================================================================================
! Rescales distribution function for EDDVM/DUGKS
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY : DVMMomentSave, DVMMethod, DVMSpecData, DVMnSpecies, DVMnMacro, DVMInnerESave, DVMnInnerE
USE MOD_Globals,        ONLY :abort
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY : nElems
USE MOD_FV_Vars,        ONLY : U_FV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)              :: tDeriv
INTEGER, INTENT(IN)           :: tilde
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: MacroVal(DVMnMacro,DVMnSpecies+1), tau, prefac, relaxFac, Pr, rho
REAL                            :: ErelaxTrans, Erot(DVMnSpecies+1), ErelaxRot(DVMnSpecies)
INTEGER                         :: i,j,k,iElem,iSpec,vFirstID,vLastID
REAL,ALLOCATABLE                :: fTarget(:)
!===================================================================================================================================
Erot = 0.
DO iElem =1, nElems
  ! DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
  DO k=0, 0; DO j=0, 0; DO i=0, 0
    SELECT CASE (tilde)
      CASE(1) ! f~  -----> f2^    (tDeriv=dt)
        CALL MacroValuesFromDistribution(MacroVal,U_FV(:,i,j,k,iElem),tDeriv,tau,tilde,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot)
        DVMMomentSave(1:DVMnMacro,:,iElem) = MacroVal(1:DVMnMacro,:)
        DVMMomentSave(DVMnMacro+1,:,iElem) = tau
        DVMMomentSave(DVMnMacro+2,:,iElem) = rho
        DVMMomentSave(DVMnMacro+3,:,iElem) = Pr
        IF (DVMnInnerE.GT.0) DVMInnerESave(1,1:DVMnSpecies+1,iElem) = Erot(1:DVMnSpecies+1)
        SELECT CASE(DVMMethod)
        CASE(0) !First order
          relaxFac = tDeriv/tau
          IF (CHECKEXP(relaxFac)) THEN
            prefac = EXP(-relaxFac)
          ELSE
            prefac = 0.
          END IF
        CASE(1) !EDDVM
          relaxFac = tDeriv/tau/2.
          IF (CHECKEXP(3.*relaxFac)) THEN
            prefac = (EXP(-relaxFac)-EXP(-3.*relaxFac))/(1.-EXP(-relaxFac))/2.
          ELSE
            prefac = 0.
          END IF
        CASE(2) !DUGKS
          prefac = (2.*tau-tDeriv/2.)/(2.*tau+tDeriv)
        END SELECT
      CASE(2) ! f2^ -----> f^     (tDeriv=dt/2)
        MacroVal(1:DVMnMacro,:) = DVMMomentSave(1:DVMnMacro,:,iElem)
        tau = DVMMomentSave(DVMnMacro+1,DVMnSpecies+1,iElem)
        rho = DVMMomentSave(DVMnMacro+2,DVMnSpecies+1,iElem)
        Pr = DVMMomentSave(DVMnMacro+3,DVMnSpecies+1,iElem)
        IF (DVMnInnerE.GT.0) Erot(1:DVMnSpecies+1) = DVMInnerESave(1,1:DVMnSpecies+1,iElem)
        SELECT CASE(DVMMethod)
        CASE(1)
          relaxFac = tDeriv/tau
          IF (CHECKEXP(2.*relaxFac)) THEN
            prefac = 2.*(EXP(-relaxFac)-EXP(-2.*relaxFac))/(1.-EXP(-2.*relaxFac))
          ELSE
            prefac = 0.
          END IF
        CASE(2)
          prefac = (4./3.)-(1./3.)*(2.*tau+2.*tDeriv)/(2.*tau-tDeriv)
        END SELECT
    END SELECT
    CALL MoleculeRelaxEnergy(ErelaxTrans, ErelaxRot, MacroVal(5,DVMnSpecies+1), Erot, Pr)
    vFirstID=1
    vLastID=0
    DO iSpec=1,DVMnSpecies
      vLastID = vLastID + DVMSpecData(iSpec)%nVar
      ALLOCATE(fTarget(DVMSpecData(iSpec)%nVar))
      CALL TargetDistribution(MacroVal(:,DVMnSpecies+1), fTarget, iSpec, MacroVal(1,iSpec), rho, Pr, ErelaxTrans, Erelaxrot(iSpec))
      U_FV(vFirstID:vLastID,i,j,k,iElem) = U_FV(vFirstID:vLastID,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
      DEALLOCATE(fTarget)
      vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
    END DO
  END DO; END DO; END DO
END DO
END SUBROUTINE RescaleU

SUBROUTINE RescaleInit(tDeriv)
!===================================================================================================================================
! Initial rescale (f->ftilde) for initialization with non equilibrium flow
! TODO: Should also be used for restart
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY: DVMMethod, DVMSpecData, DVMnSpecies, DVMnMacro
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY : nElems
USE MOD_FV_Vars,        ONLY : U_FV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)              :: tDeriv
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: MacroVal(DVMnMacro,DVMnSpecies+1), tau, prefac, rho, Pr
INTEGER                         :: i,j,k,iElem,iSpec,vFirstID,vLastID
REAL, ALLOCATABLE               :: fTarget(:)
REAL                            :: ErelaxTrans, Erot(DVMnSpecies+1), ErelaxRot(DVMnSpecies)
!===================================================================================================================================
SWRITE(UNIT_stdOut,*) 'INITIAL DISTRIBUTION FUNCTION RESCALE'
vFirstID=1
vLastID=0
DO iElem =1, nElems
  DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
    CALL MacroValuesFromDistribution(MacroVal,U_FV(:,i,j,k,iElem),0.,tau,1,MassDensity=rho,PrandtlNumber=Pr,Erot=Erot) ! tDeriv=0 to get heatflux from original distribution
    SELECT CASE (DVMMethod)
      CASE(1)
        prefac = (tDeriv/tau)/(1. - (EXP(-tDeriv/tau)))
      CASE(2)
        prefac = (2.*tau+tDeriv)/(2.*tau)
    END SELECT
    CALL MoleculeRelaxEnergy(ErelaxTrans, ErelaxRot, MacroVal(5,DVMnSpecies+1), Erot(1:DVMnSpecies), Pr)
    DO iSpec=1,DVMnSpecies
      vLastID = vLastID + DVMSpecData(iSpec)%nVar
      ALLOCATE(fTarget(DVMSpecData(iSpec)%nVar))
      CALL TargetDistribution(MacroVal(:,DVMnSpecies+1), fTarget, iSpec, MacroVal(1,iSpec), rho, Pr, ErelaxTrans, Erelaxrot(iSpec))
      U_FV(vFirstID:vLastID,i,j,k,iElem) = U_FV(vFirstID:vLastID,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
      DEALLOCATE(fTarget)
      vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
    END DO
  END DO; END DO; END DO
END DO
END SUBROUTINE RescaleInit

SUBROUTINE ForceStep(tDeriv,ploesma)
!===================================================================================================================================
! Calculates force term (to add in 2 parts (Strang splitting) for 2nd order accuracy)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY: DVMSpecData, DVMnSpecies, DVMDim, DVMForce, DVMnMacro
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY : nElems
USE MOD_FV_Vars,        ONLY : U_FV
#if USE_HDG
USE MOD_Equation_Vars,  ONLY: E
USE MOD_Interpolation_Vars,ONLY: wGP
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)              :: tDeriv
LOGICAL, OPTIONAL             :: ploesma
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: MacroVal(DVMnMacro,DVMnSpecies+1), tau, forceTerm, forceTerm2, cVel(3), velodiff, gamma
INTEGER                         :: i,j,k,iElem,iVel,jVel,kVel,upos,iSpec,vFirstID ,upos1,upos2
REAL, ALLOCATABLE               :: fTarget(:)
REAL                            :: Eloc(3)
!===================================================================================================================================
DO iElem =1, nElems
#if USE_HDG
  IF (PRESENT(ploesma)) THEN
    Eloc = 0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Eloc(1:3) = Eloc(1:3) + wGP(i)*wGP(j)*wGP(k)*E(1:3,i,j,k,iElem)/((PP_N+1.)**3) !need jacobi here for noncartesian
    END DO; END DO; END DO
  END IF
#endif
  ! DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
  DO k=0, 0; DO j=0, 0; DO i=0, 0
    ! CALL MacroValuesFromDistribution(MacroVal,U_FV(:,i,j,k,iElem),tDeriv,tau,1)
    ! SELECT CASE (DVMBGKModel)
    !   CASE(1)
    !     CALL MaxwellDistribution(MacroVal,fTarget)
    !   CASE(2)
    !     CALL ShakhovDistribution(MacroVal,fTarget)
    !   CASE DEFAULT
    !     CALL abort(__STAMP__,'DVM BGK Model not implemented.')
    !   END SELECT
    ! gamma = tau*(1.-EXP(-tDeriv/tau))/tDeriv
    vFirstID = 0
    DO iSpec=1,DVMnSpecies
      ASSOCIATE(Sp => DVMSpecData(iSpec))
      ! ALLOCATE(fTarget(Sp%nVar))
      ! CALL MaxwellDistribution(MacroVal(1:DVMnMacro,iSpec),fTarget,iSpec) !species-specific equilibrium approximation (bad idea?)

      DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
        upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2)
        !equilibrium approximation
        ! cVel(1) = Sp%Velos(iVel,1) - MacroVal(2,iSpec)
        ! cVel(2) = Sp%Velos(jVel,2) - MacroVal(3,iSpec)
        ! cVel(3) = Sp%Velos(kVel,3) - MacroVal(4,iSpec)
        ! IF (PRESENT(ploesma)) THEN
        !   forceTerm = (Sp%Charge/Sp%Mass) &
        !             * DOT_PRODUCT(Eloc,cVel)/(Sp%R_S*MacroVal(5,iSpec)) * fTarget(upos)
        ! ELSE
        !   forceTerm = DOT_PRODUCT(DVMForce,cVel)/(Sp%R_S*MacroVal(5,iSpec)) * fTarget(upos)
        ! END IF
        ! non equilibrium version
        IF (iVel.EQ.1) THEN
          upos1=upos
          upos2 = upos + 1 !iVel+1+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2)
          velodiff=Sp%Velos(iVel+1,1)-Sp%Velos(iVel,1)
        ELSE IF (iVel.EQ.Sp%nVelos(1)) THEN
          upos1 = upos - 1 !iVel-1+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2)
          upos2=upos
          velodiff=Sp%Velos(iVel,1)-Sp%Velos(iVel-1,1)
        ELSE
          upos1 = upos - 1 !iVel-1+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2)
          upos2 = upos + 1 !iVel+1+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2)
          velodiff=Sp%Velos(iVel+1,1)-Sp%Velos(iVel-1,1)
        END IF
        IF (PRESENT(ploesma)) THEN
          forceTerm = - (Sp%Charge/Sp%Mass)*Eloc(1)*(U_FV(upos2+vFirstID,i,j,k,iElem)-U_FV(upos1+vFirstID,i,j,k,iElem))/velodiff
        ELSE
          forceTerm = - DVMForce(1)*(U_FV(upos2+vFirstID,i,j,k,iElem)-U_FV(upos1+vFirstID,i,j,k,iElem))/velodiff
        END IF
        ! forceTerm = - DVMForce(1)*(gamma*(U(upos2,i,j,k,iElem)-U(upos1,i,j,k,iElem)) &
        !                        +(1-gamma)*(fTarget(upos2)-fTarget(upos1)))/velodiff

        U_FV(upos+vFirstID,i,j,k,iElem) = U_FV(upos+vFirstID,i,j,k,iElem) + forceTerm*tDeriv/2. !t/2 for strang splitting
        IF (DVMDim.LT.3) THEN
          IF (PRESENT(ploesma)) THEN
            forceTerm2 = - (Sp%Charge/Sp%Mass)*Eloc(1)*(U_FV(Sp%nVarReduced+upos2+vFirstID,i,j,k,iElem)-U_FV(Sp%nVarReduced+upos1+vFirstID,i,j,k,iElem))/velodiff
          ELSE
            forceTerm2 = - DVMForce(1)*(U_FV(Sp%nVarReduced+upos2+vFirstID,i,j,k,iElem)-U_FV(Sp%nVarReduced+upos1+vFirstID,i,j,k,iElem))/velodiff
          END IF
          U_FV(Sp%nVarReduced+upos+vFirstID,i,j,k,iElem) = U_FV(Sp%nVarReduced+upos+vFirstID,i,j,k,iElem) + forceTerm2*tDeriv/2.
        END IF
      END DO; END DO; END DO
      vFirstID = vFirstID + Sp%nVar
      ! DEALLOCATE(fTarget)
      END ASSOCIATE
    END DO
  END DO; END DO; END DO
END DO
END SUBROUTINE ForceStep


SUBROUTINE IntegrateFluxValues(MacroVal,U)
!===================================================================================================================================
! Calculates the surface macro values from distribution fluxes
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpecData, DVMnSpecies, DVMDim
USE MOD_PreProc
USE MOD_Globals               ,ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: U(PP_nVar_FV)
REAL, INTENT(OUT)               :: MacroVal(5)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, rhoU(3), densE, weight
INTEGER                         :: iVel,jVel,kVel, upos, iSpec,vFirstID
!===================================================================================================================================
MacroVal = 0.

vFirstID = 0
DO iSpec=1, DVMnSpecies
  rho = 0.
  rhoU = 0.
  densE = 0.
  ASSOCIATE(Sp    => DVMSpecData(iSpec))
  DO kVel=1, Sp%nVelos(3);   DO jVel=1, Sp%nVelos(2);   DO iVel=1, Sp%nVelos(1)
    upos= iVel+(jVel-1)*Sp%nVelos(1)+(kVel-1)*Sp%nVelos(1)*Sp%nVelos(2) + vFirstID
    weight = Sp%Weights(iVel,1)*Sp%Weights(jVel,2)*Sp%Weights(kVel,3)*Sp%Mass
    rho = rho + weight*U(upos)
    rhoU(1) = rhoU(1) + weight*Sp%Velos(iVel,1)*U(upos)
    rhoU(2) = rhoU(2) + weight*Sp%Velos(jVel,2)*U(upos)
    rhoU(3) = rhoU(3) + weight*Sp%Velos(kVel,3)*U(upos)
    IF (DVMDim.LT.3) THEN
      densE = densE + weight*0.5*((Sp%Velos(iVel,1)**2.+Sp%Velos(jVel,2)**2.+Sp%Velos(kVel,3)**2.)*U(upos)+U(Sp%nVarReduced+upos))
    ELSE
      densE = densE + weight*0.5*(Sp%Velos(iVel,1)**2.+Sp%Velos(jVel,2)**2.+Sp%Velos(kVel,3)**2.)*U(upos)
    END IF
  END DO; END DO; END DO

  MacroVal(1) = MacroVal(1) + rho ! mass flow
  MacroVal(2:4) = MacroVal(2:4) + rhoU(1:3) ! force per area
  MacroVal(5) = MacroVal(5) + densE ! heat flux

  vFirstID = vFirstID + Sp%nVar
  END ASSOCIATE
END DO

END SUBROUTINE IntegrateFluxValues

END MODULE MOD_DistFunc