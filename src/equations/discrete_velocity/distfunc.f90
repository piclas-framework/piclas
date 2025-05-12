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
PUBLIC:: TargetDistribution
PUBLIC:: MaxwellScattering, RescaleU, RescaleInit, ForceStep, IntegrateFluxValues
!==================================================================================================================================

CONTAINS

SUBROUTINE MacroValuesFromDistribution(MacroVal,U,tDeriv,tau,tilde)
!===================================================================================================================================
! Calculates the moments of the distribution function
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMWeights, DVMSpeciesData, DVMMethod, DVMBGKModel, DVMDim
USE MOD_PreProc
USE MOD_Globals               ,ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: U(PP_nVar_FV), tDeriv
INTEGER,INTENT(IN)              :: tilde
REAL, INTENT(OUT)               :: MacroVal(14), tau
! REAL, INTENT(OUT), OPTIONAL     :: skewness(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, rhoU(3), rhoE, PressTens(6), Heatflux(3), uVelo(3), cV, cVel(3),cMag2, mu, weight, prefac
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = 0.
rhoU = 0.
rhoE = 0.
PressTens = 0.
Heatflux = 0.
MacroVal = 0.
! IF (PRESENT(skewness)) skewness = 0.
DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  rho = rho + weight*U(upos)
  rhoU(1) = rhoU(1) + weight*DVMVelos(iVel,1)*U(upos)
  rhoU(2) = rhoU(2) + weight*DVMVelos(jVel,2)*U(upos)
  rhoU(3) = rhoU(3) + weight*DVMVelos(kVel,3)*U(upos)
  IF (DVMDim.LT.3) THEN
    rhoE = rhoE + weight*0.5*((DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)+U(PP_nVar_FV/2+upos))
  ELSE
    rhoE = rhoE + weight*0.5*(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)
  END IF
END DO; END DO; END DO

uVelo = rhoU/rho
cV = DVMSpeciesData%R_S*rho*(3.+DVMSpeciesData%Internal_DOF)/2.

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag2 = DOT_PRODUCT(cVel,cVel)
  PressTens(1) = PressTens(1) + weight*cVel(1)*cVel(1)*U(upos)
  IF (DVMDim.LT.3) THEN
    IF (DVMDim.EQ.1) THEN
      PressTens(2) = PressTens(2) + weight*0.5*U(PP_nVar_FV/2+upos)
      PressTens(3) = PressTens(3) + weight*0.5*U(PP_nVar_FV/2+upos)
    ELSE ! DVMDim.EQ.2
      PressTens(4) = PressTens(4) + weight*cVel(1)*cVel(2)*U(upos)
      PressTens(2) = PressTens(2) + weight*cVel(2)*cVel(2)*U(upos)
      PressTens(3) = PressTens(3) + weight*U(PP_nVar_FV/2+upos)
    END IF
    Heatflux(1) = Heatflux(1) + weight*0.5*cVel(1)*(U(upos)*cMag2+U(PP_nVar_FV/2+upos))
    Heatflux(2) = Heatflux(2) + weight*0.5*cVel(2)*(U(upos)*cMag2+U(PP_nVar_FV/2+upos))
    Heatflux(3) = Heatflux(3) + weight*0.5*cVel(3)*(U(upos)*cMag2+U(PP_nVar_FV/2+upos))
  ELSE
    PressTens(2) = PressTens(2) + weight*cVel(2)*cVel(2)*U(upos)
    PressTens(3) = PressTens(3) + weight*cVel(3)*cVel(3)*U(upos)
    PressTens(4) = PressTens(4) + weight*cVel(1)*cVel(2)*U(upos)
    PressTens(5) = PressTens(5) + weight*cVel(1)*cVel(3)*U(upos)
    PressTens(6) = PressTens(6) + weight*cVel(2)*cVel(3)*U(upos)
    Heatflux(1:3) = Heatflux(1:3) + weight*0.5*cVel(1:3)*(U(upos)*cMag2)
    ! IF (PRESENT(skewness)) THEN
    !   skewness(1) = skewness(1) + weight*cVel(1)*(U(upos)*cVel(1)*cVel(1))
    !   skewness(2) = skewness(2) + weight*cVel(2)*(U(upos)*cVel(2)*cVel(2))
    !   skewness(3) = skewness(3) + weight*cVel(3)*(U(upos)*cVel(3)*cVel(3))
    ! END IF
  END IF
END DO; END DO; END DO

MacroVal(1) = rho
MacroVal(2:4) = uVelo
MacroVal(5) = (rhoE - 0.5*(DOT_PRODUCT(rhoU,rhoU))/rho)/cV
IF (MacroVal(5).LE.0) CALL abort(__STAMP__,'DVM negative temperature!')

mu = DVMSpeciesData%mu_Ref*(MacroVal(5)/DVMSpeciesData%T_Ref)**(DVMSpeciesData%omegaVHS+0.5)
tau = mu/(DVMSpeciesData%R_S*MacroVal(1)*MacroVal(5))
IF (DVMBGKModel.EQ.1) tau = tau/DVMSpeciesData%Prandtl !ESBGK
IF (DVMBGKModel.EQ.7) tau = tau*2./3.

Macroval(6:11)  = PressTens(1:6)
MacroVal(12:14) = Heatflux(1:3)
IF (tDeriv.GT.0.) THEN
  SELECT CASE (tilde)
  CASE(1) ! higher moments from f~
    SELECT CASE(DVMMethod)
    CASE(1) !EDDVM
      prefac = (1.-EXP(-tDeriv/tau))/(tDeriv/tau)
      SELECT CASE(DVMBGKModel)
      CASE(1,4) !ESBGK
        Macroval(6:8)   = (Macroval(6:8)*prefac+(1.-prefac)*MacroVal(5)*DVMSpeciesData%R_S*rho/DVMSpeciesData%Prandtl) &
                          /(1./DVMSpeciesData%Prandtl+prefac*(1.-1./DVMSpeciesData%Prandtl))
        MacroVal(9:11)  = Macroval(9:11)*prefac/(1./DVMSpeciesData%Prandtl+prefac*(1.-1./DVMSpeciesData%Prandtl))
        MacroVal(12:14) = MacroVal(12:14)*prefac
      CASE(2,6) !Shakhov/SN
        MacroVal(6:8)   = Macroval(6:8)*prefac+(1.-prefac)*MacroVal(5)*DVMSpeciesData%R_S*rho
        MacroVal(9:11)  = MacroVal(9:11)*prefac
        MacroVal(12:14) = MacroVal(12:14)*prefac/(DVMSpeciesData%Prandtl+prefac*(1-DVMSpeciesData%Prandtl))
      CASE(3,5) !Maxwell
        MacroVal(6:8)   = Macroval(6:8)*prefac+(1.-prefac)*MacroVal(5)*DVMSpeciesData%R_S*rho
        Macroval(9:14)  = prefac*Macroval(9:14) ! non-eq moments should be zero anyway
      CASE(7) !Double moment distributions
        Macroval(6:8)   = (Macroval(6:8)*prefac+(1.-prefac)*MacroVal(5)*DVMSpeciesData%R_S*rho*2./3.) &
        /(2./3.+prefac/3.)
        MacroVal(9:11)  = Macroval(9:11)*prefac/(2./3.+prefac/3.)
        MacroVal(12:14) = MacroVal(12:14)*prefac/(2.*DVMSpeciesData%Prandtl/3.+prefac*(1.-2.*DVMSpeciesData%Prandtl/3.))
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM-BGKCollModel does not exist')
      END SELECT
    CASE(2) !DUGKS
      SELECT CASE(DVMBGKModel)
      CASE(1,4) !ESBGK
        prefac = (2.*tau)/(2.*tau+tDeriv)
        Macroval(6:8)   = (Macroval(6:8)*prefac+(1.-prefac)*MacroVal(5)*DVMSpeciesData%R_S*rho/DVMSpeciesData%Prandtl) &
                          /(1./DVMSpeciesData%Prandtl+prefac*(1.-1./DVMSpeciesData%Prandtl))
        MacroVal(9:11)  = Macroval(9:11)*prefac/(1./DVMSpeciesData%Prandtl+prefac*(1.-1./DVMSpeciesData%Prandtl))
        MacroVal(12:14) = MacroVal(12:14)*prefac
      CASE(2,6) !Shakhov/SN
        Macroval(6:8)   = Macroval(6:8)*2.*tau/(2.*tau+tDeriv) + MacroVal(5)*DVMSpeciesData%R_S*rho*tDeriv/(2.*tau+tDeriv)
        Macroval(9:11)  = Macroval(9:11)*2.*tau/(2.*tau+tDeriv)
        MacroVal(12:14) = MacroVal(12:14)*2.*tau/(2.*tau+tDeriv*DVMSpeciesData%Prandtl)
      CASE(3,5) !Maxwell
        Macroval(6:8)   = Macroval(6:8)*2.*tau/(2.*tau+tDeriv) + MacroVal(5)*DVMSpeciesData%R_S*rho*tDeriv/(2.*tau+tDeriv)
        Macroval(9:14)  = Macroval(9:14)*2.*tau/(2.*tau+tDeriv) ! non-eq moments should be zero anyway
      CASE(7) !Double moment distributions
        Macroval(6:8)   = Macroval(6:8)*2.*tau/(2.*tau++2.*tDeriv/3.) &
                        + 2.*MacroVal(5)*DVMSpeciesData%R_S*rho*tDeriv/(2.*tau+2.*tDeriv/3.)/3.
        Macroval(9:11)  = Macroval(9:11)*2.*tau/(2.*tau+2.*tDeriv/3.)
        MacroVal(12:14) = MacroVal(12:14)*2.*tau/(2.*tau+2.*tDeriv*DVMSpeciesData%Prandtl/3.)
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM-BGKCollModel does not exist')
      END SELECT
    END SELECT
  CASE(2) ! higher moments from f^
    MacroVal(6:14) = 0. ! will get copied from earlier f~ macroval
  CASE DEFAULT
    CALL abort(__STAMP__,'DVM-Method does not exist')
  END SELECT
END IF

END SUBROUTINE MacroValuesFromDistribution

SUBROUTINE MaxwellDistributionCons(MacroVal,fMaxwell)
!===================================================================================================================================
! conservative maxwell distribution (cf Mieussens 2000)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi, DVMWeights
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fMaxwell(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, Temp, uVelo(DVMDim), vMag, weight, gM
REAL,DIMENSION(2+DVMDim)        :: alpha, psi, rhovec
REAL                            :: J(2+DVMDim,2+DVMDim), B(2+DVMDim)
INTEGER                         :: iVel,jVel,kVel, upos, countz, IPIV(2+DVMDim), info_dgesv
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:DVMDim) = MacroVal(2:1+DVMDim)
Temp = MacroVal(5)

! vector of conservative variables
rhovec(1) = rho
rhovec(2:1+DVMDim)=rho*uVelo(:)
rhovec(2+DVMDim)=rho*(FLOAT(DVMDim)*DVMSpeciesData%R_S*Temp+DOT_PRODUCT(uVelo,uVelo))/2.

alpha(1) = LOG(rho/(2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))-DOT_PRODUCT(uVelo,uVelo)/2./DVMSpeciesData%R_S/Temp
alpha(2:1+DVMDim) = uVelo(1:DVMDim)/DVMSpeciesData%R_S/Temp
alpha(2+DVMDim) = -1/DVMSpeciesData%R_S/Temp

! init counter
countz=0

!Newton algorithm to find alpha so that <psi exp(alpha.psi)> = rhovec
DO WHILE (countz.LT.1000)
  countz=countz+1
  J = 0.

  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
    vMag = DVMVelos(iVel,1)**2 + DVMVelos(jVel,2)**2 + DVMVelos(kVel,3)**2
    psi(1)=1
    psi(2)=DVMVelos(iVel,1)
    IF (DVMDim.GT.1) psi(3)=DVMVelos(jVel,2)
    IF (DVMDim.GT.2) psi(4)=DVMVelos(kVel,3)
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
DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  vMag = DVMVelos(iVel,1)**2 + DVMVelos(jVel,2)**2 + DVMVelos(kVel,3)**2
  psi(1)=1
  psi(2)=DVMVelos(iVel,1)
  IF (DVMDim.GT.1) psi(3)=DVMVelos(jVel,2)
  IF (DVMDim.GT.2) psi(4)=DVMVelos(kVel,3)
  psi(2+DVMDim) = vMag/2.
  fMaxwell(upos) = EXP(DOT_PRODUCT(alpha,psi))
  IF (DVMDim.LT.3) THEN
    fMaxwell(PP_nVar_FV/2+upos) = fMaxwell(upos)*DVMSpeciesData%R_S*Temp*(DVMSpeciesData%Internal_DOF+3.-DVMDim)
  END IF
END DO; END DO; END DO

END SUBROUTINE MaxwellDistributionCons

SUBROUTINE MaxwellDistribution(MacroVal,fMaxwell)
!===================================================================================================================================
! Maxwell distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fMaxwell(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, Temp, uVelo(3), cVel(3), cMag2
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag2 = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  fMaxwell(upos) = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpeciesData%R_S*Temp))
  IF (DVMDim.LT.3) THEN
    fMaxwell(PP_nVar_FV/2+upos) = fMaxwell(upos)*DVMSpeciesData%R_S*Temp*(DVMSpeciesData%Internal_DOF+3.-DVMDim)
  END IF
END DO; END DO; END DO

END SUBROUTINE MaxwellDistribution

SUBROUTINE ShakhovDistribution(MacroVal,fShakhov)
!===================================================================================================================================
! Shakhov distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fShakhov(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, Temp, uVelo(3), cVel(3), cMag2, gM, q(3), ShakhFac1, ShakhFac2
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)
q(1:3) = MacroVal(12:14)

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag2 = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  gM = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpeciesData%R_S*Temp))
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*rho*DVMSpeciesData%R_S*DVMSpeciesData%R_S*Temp*Temp)
  ShakhFac2 = cMag2/(DVMSpeciesData%R_S*Temp)
  fShakhov(upos) = gM*(1.+(1.-DVMSpeciesData%Prandtl)*ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMDim.LT.3) THEN
    fShakhov(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*(DVMSpeciesData%Internal_DOF+3.-DVMDim &
                                + (1.-DVMSpeciesData%Prandtl)*ShakhFac1*((ShakhFac2-DVMDim)*(DVMSpeciesData%Internal_DOF+3.-DVMDim)&
                                - 2.*DVMSpeciesData%Internal_DOF))
  END IF
END DO; END DO; END DO

END SUBROUTINE ShakhovDistribution

SUBROUTINE ESBGKDistribution(MacroVal,fESBGK)
!===================================================================================================================================
! ESBGK distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
USE MOD_Basis                    ,ONLY: INV33
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fESBGK(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho,Temp,uVelo(3),cVel(3),pressTens(3,3),pressProduct,ilambda(3,3),ldet,hfac
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho              = MacroVal(1)
uVelo(1:3)       = MacroVal(2:4)
Temp             = MacroVal(5)
pressTens(1,1)   = MacroVal(6)
pressTens(2,2)   = MacroVal(7)
pressTens(3,3)   = MacroVal(8)
pressTens(1,2:3) = MacroVal(9:10)
pressTens(2:3,1) = MacroVal(9:10)
pressTens(2,3)   = MacroVal(11)
pressTens(3,2)   = MacroVal(11)

pressTens = (1.-1./DVMSpeciesData%Prandtl)*pressTens/rho
pressTens(1,1) = pressTens(1,1)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl
pressTens(2,2) = pressTens(2,2)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl
pressTens(3,3) = pressTens(3,3)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl

CALL INV33(pressTens,ilambda,ldet)
IF (ldet.LE.0.) THEN
  CALL abort(__STAMP__,'ESBGK matrix not positive-definite')
ELSE
  ! determinant from reduced pressure tensor (size D*D)
  IF (DVMDim.LE.2) ldet = ldet/pressTens(3,3)
  IF (DVMDim.LE.1) ldet = ldet/pressTens(2,2)
END IF

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)

  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)

  pressProduct = cVel(1)*DOT_PRODUCT(ilambda(:,1),cVel) &
                + cVel(2)*DOT_PRODUCT(ilambda(:,2),cVel) &
                + cVel(3)*DOT_PRODUCT(ilambda(:,3),cVel)

  fESBGK(upos) = rho/sqrt(ldet*(2*Pi)**DVMDim)*EXP(-pressProduct/2.)
  IF ((DVMSpeciesData%Internal_DOF .GT.0.0).OR.(DVMDim.LT.3)) THEN
    hfac = DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF
    IF (DVMDim.LE.2) hfac = hfac + pressTens(3,3)
    IF (DVMDim.LE.1) hfac = hfac + pressTens(2,2)
    fESBGK(PP_nVar_FV/2+upos) = fESBGK(upos)*hfac
  END IF
END DO; END DO; END DO

END SUBROUTINE ESBGKDistribution

SUBROUTINE ESBGKDistributionCons(MacroVal,fESBGK)
!===================================================================================================================================
! Conservative ESBGK distribution from macro values (cf Andries et al. 2000)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi, DVMWeights
USE MOD_Basis                    ,ONLY: INV33
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fESBGK(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                              :: rho,Temp,uVelo(DVMDim),pressTens(3,3),pressVec(2*DVMDim),weight,gM,hfac
REAL                                              :: iLambda(3,3),ldet, pressProduct(DVMDim)
REAL,DIMENSION(1+2*DVMDim+DVMDim*(DVMDim-1)/2)    :: rhovec,alpha,psi,B
REAL                                              :: J(1+2*DVMDim+DVMDim*(DVMDim-1)/2,1+2*DVMDim+DVMDim*(DVMDim-1)/2)
INTEGER                                           :: iVel,jVel,kVel,upos,countz,IPIV(1+2*DVMDim+DVMDim*(DVMDim-1)/2),info_dgesv
INTEGER                                           :: iMat,jMat
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:DVMDim) = MacroVal(2:1+DVMDim)
Temp = MacroVal(5)
pressTens(1,1)   = MacroVal(6)
pressTens(2,2)   = MacroVal(7)
pressTens(3,3)   = MacroVal(8)
pressTens(1,2:3) = MacroVal(9:10)
pressTens(2:3,1) = MacroVal(9:10)
pressTens(2,3)   = MacroVal(11)
pressTens(3,2)   = MacroVal(11)

pressTens = (1.-1./DVMSpeciesData%Prandtl)*pressTens/rho
pressTens(1,1) = pressTens(1,1)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl
pressTens(2,2) = pressTens(2,2)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl
pressTens(3,3) = pressTens(3,3)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl

CALL INV33(pressTens,ilambda,ldet)
IF (ldet.LE.0.) THEN
  CALL abort(__STAMP__,'ESBGK matrix not positive-definite')
ELSE
  ! determinant from reduced pressure tensor (size D*D)
  IF (DVMDim.LE.2) ldet = ldet/pressTens(3,3)
  IF (DVMDim.LE.1) ldet = ldet/pressTens(2,2)
END IF

! pressure tensor to vector
pressVec = 0.
pressVec(1:DVMDim)   = MacroVal(6:6+DVMDim-1)
pressVec(DVMDim+1:DVMDim+DVMDim*(DVMDim-1)/2) = MacroVal(9:9+DVMDim*(DVMDim-1)/2-1)

pressVec = (1.-1./DVMSpeciesData%Prandtl)*pressVec/rho
pressVec(1:DVMDim) = pressVec(1:DVMDim)+DVMSpeciesData%R_S*Temp/DVMSpeciesData%Prandtl

! vector of conservative variables + cross terms
rhovec(1) = rho
rhovec(2:1+DVMDim)=rho*uVelo(:)
rhovec(DVMDim+2:2*DVMDim+1)=rho*(uVelo(1:DVMDim)*uVelo(1:DVMDim) + pressVec(1:DVMDim))
IF (DVMDim.GT.1) rhovec(2*DVMDim+2) = rho*(uVelo(1)*uVelo(2) + pressVec(DVMDim+1))
IF (DVMDim.GT.2) THEN
  rhovec(2*DVMDim+3) = rho*(uVelo(1)*uVelo(3) + pressVec(DVMDim+2))
  rhovec(2*DVMDim+4) = rho*(uVelo(2)*uVelo(3) + pressVec(DVMDim+3))
END IF

pressProduct(1:DVMDim) = MATMUL(iLambda(1:DVMDim,1:DVMDim),uVelo)

alpha = 0.
alpha(1) = LOG(rho/sqrt(ldet*(2*Pi)**DVMDim))-DOT_PRODUCT(pressProduct,uVelo)/2.
alpha(2:1+DVMDim) = pressProduct(1:DVMDim)
alpha(2+DVMDim) = -iLambda(1,1)/2.
IF (DVMDim.GT.1) THEN
  alpha(3+DVMDim)= -iLambda(2,2)/2.
  alpha(2+2*DVMDim) = -iLambda(1,2)/2.
  IF (DVMDim.GT.2) THEN
    alpha(4+DVMDim) = -iLambda(3,3)/2.
    alpha(3+2*DVMDim) = -iLambda(1,3)/2.
    alpha(4+2*DVMDim) = -iLambda(2,3)/2.
  END IF
END IF

! init counter
countz=0

!Newton algorithm to find alpha so that <psi exp(alpha.psi)> = rhovec
DO WHILE (countz.LT.1000)
  countz=countz+1
  J = 0.

  DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
    weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
    psi(1)=1
    psi(2)=DVMVelos(iVel,1)
    psi(2+DVMDim) = DVMVelos(iVel,1)**2
    IF (DVMDim.GT.1) THEN
      psi(3)=DVMVelos(jVel,2)
      psi(3+DVMDim)=DVMVelos(jVel,2)**2
      psi(2+2*DVMDim) = DVMVelos(iVel,1)*DVMVelos(jVel,2)
      IF (DVMDim.GT.2) THEN
        psi(4)=DVMVelos(kVel,3)
        psi(4+DVMDim)=DVMVelos(kVel,3)**2
        psi(3+2*DVMDim) = DVMVelos(iVel,1)*DVMVelos(kVel,3)
        psi(4+2*DVMDim) = DVMVelos(jVel,2)*DVMVelos(kVel,3)
      END IF
    END IF
    gM = EXP(DOT_PRODUCT(alpha,psi))

    ! J: derivative of -B = <psi exp(alpha.psi)> - rhovec
    ! J = <psi x psi exp(alpha.psi)>
    DO iMat=1,1+2*DVMDim+DVMDim*(DVMDim-1)/2
      DO jMat=1,1+2*DVMDim+DVMDim*(DVMDim-1)/2
        J(iMat,jMat) = J(iMat,jMat) + weight*gM*psi(iMat)*psi(jMat)
      END DO
    END DO
  END DO; END DO; END DO

  B(:) = rhovec(:) - J(:,1)

  ! solve JX = B
  CALL DGESV(1+2*DVMDim+DVMDim*(DVMDim-1)/2,1,J,1+2*DVMDim+DVMDim*(DVMDim-1)/2,IPIV,B,1+2*DVMDim+DVMDim*(DVMDim-1)/2,info_dgesv)

  IF(info_dgesv.NE.0) CALL abort(__STAMP__,'Newton DGESV fail')
  IF (.NOT.(ANY(ABS(B(:)).GT.(1e-5*ABS(alpha)+1e-12)))) EXIT
  ! IF (.NOT.(ANY(ABS(B(:,1)).GT.(1e-6*ABS(alpha)+1e-16)))) EXIT

  ! update alpha with Newton increment X (stored in B)
  alpha = alpha + B
END DO
IF (countz.GE.1000) print*, 'Newton max iter reached'

alpha = alpha + B
DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  psi(1)=1
  psi(2)=DVMVelos(iVel,1)
  psi(2+DVMDim) = DVMVelos(iVel,1)**2
  IF (DVMDim.GT.1) THEN
    psi(3)=DVMVelos(jVel,2)
    psi(3+DVMDim)=DVMVelos(jVel,2)**2
    psi(2+2*DVMDim) = DVMVelos(iVel,1)*DVMVelos(jVel,2)
    IF (DVMDim.GT.2) THEN
      psi(4)=DVMVelos(kVel,3)
      psi(4+DVMDim)=DVMVelos(kVel,3)**2
      psi(3+2*DVMDim) = DVMVelos(iVel,1)*DVMVelos(kVel,3)
      psi(4+2*DVMDim) = DVMVelos(jVel,2)*DVMVelos(kVel,3)
    END IF
  END IF
  fESBGK(upos) = EXP(DOT_PRODUCT(alpha,psi))
  IF ((DVMSpeciesData%Internal_DOF .GT.0.0).OR.(DVMDim.LT.3)) THEN
    hfac = DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF
    IF (DVMDim.LE.2) hfac = hfac + pressTens(3,3)
    IF (DVMDim.LE.1) hfac = hfac + pressTens(2,2)
    fESBGK(PP_nVar_FV/2+upos) = fESBGK(upos)*hfac
  END IF
END DO; END DO; END DO

END SUBROUTINE ESBGKDistributionCons

SUBROUTINE GradDistribution(MacroVal,fGrad)
!===================================================================================================================================
! Grad 13 moments distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fGrad(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho,Temp,uVelo(3),cVel(3),cMag2,gM,q(3),pressTens(3,3),ShakhFac1,ShakhFac2,pressFac,pressProduct
REAL                            :: pressProduct2
INTEGER                         :: iVel,jVel,kVel,upos
!===================================================================================================================================
rho              = MacroVal(1)
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
IF (ABS(SUM(MacroVal(6:8))).GT.1.e-12*(rho*DVMSpeciesData%R_S*Temp)) CALL abort(__STAMP__, &
                                  'Diagonal entries of the pressure tensor should add up to zero',0,SUM(MacroVal(6:8)))

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)

  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag2 = DOT_PRODUCT(cVel,cVel)

  pressProduct = cVel(1)*DOT_PRODUCT(pressTens(:,1),cVel) &
               + cVel(2)*DOT_PRODUCT(pressTens(:,2),cVel) &
               + cVel(3)*DOT_PRODUCT(pressTens(:,3),cVel)

  IF (DVMDim.LT.3) THEN
    pressProduct2 = pressProduct + (5.-DVMDim)*DVMSpeciesData%R_S*Temp*pressTens(3,3)
    pressProduct = pressProduct + (3.-DVMDim)*DVMSpeciesData%R_S*Temp*pressTens(3,3)
  END IF

  pressFac = rho*DVMSpeciesData%R_S*DVMSpeciesData%R_S*Temp*Temp
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*pressFac)
  ShakhFac2 = cMag2/(DVMSpeciesData%R_S*Temp)

  gM = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag2/(2.*DVMSpeciesData%R_S*Temp))
  fGrad(upos) = gM*(1.+0.5*pressProduct/pressFac+ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMDim.LT.3) THEN
    fGrad(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*(DVMSpeciesData%Internal_DOF+3.-DVMDim) &
                  *(1 + 0.5*pressProduct2/pressFac + ShakhFac1*(ShakhFac2-DVMDim))
  END IF
END DO; END DO; END DO

END SUBROUTINE GradDistribution

SUBROUTINE GradDistributionPrandtl(MacroVal,fGrad)
!===================================================================================================================================
! Grad 13 moments distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSpeciesData
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fGrad(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: MacroValPrandtl(14)
!===================================================================================================================================
MacroValPrandtl(1:5)   = MacroVal(1:5)

MacroValPrandtl(6)     = (MacroVal(6)-MacroVal(1)*DVMSpeciesData%R_S*MacroVal(5))/3.
MacroValPrandtl(7)     = (MacroVal(7)-MacroVal(1)*DVMSpeciesData%R_S*MacroVal(5))/3.
MacroValPrandtl(8)     = (MacroVal(8)-MacroVal(1)*DVMSpeciesData%R_S*MacroVal(5))/3.
MacroValPrandtl(9:11)  = MacroVal(9:11)/3.

MacroValPrandtl(12:14) = (1.-2.*DVMSpeciesData%Prandtl/3.)*MacroVal(12:14)

CALL GradDistribution(MacroValPrandtl,fGrad)

END SUBROUTINE GradDistributionPrandtl

SUBROUTINE SkewNormalDistribution(MacroVal,fSkew)
!===================================================================================================================================
! Skew-normal distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)                 :: fSkew(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
! REAL, INTENT(IN),OPTIONAL        :: skewness(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, Temp, uVelo(3), cVel(3), cMag2, q(3)
INTEGER                         :: iVel,jVel,kVel, upos
REAL                            :: skew(1:3), delta(1:3), alpha(1:3), ksi(1:3), omega(1:3), Phi(1:3), max_skew
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)
q(1:3) = MacroVal(12:14)

max_skew = 0.5 * (4. - Pi) * (2. / (Pi - 2.)) ** 1.5
skew = (1-DVMSpeciesData%Prandtl)*2.*(q/rho)*(DVMSpeciesData%R_S*Temp)**(-3./2.)
! skew = (1-DVMSpeciesData%Prandtl)*(skewness/rho)*(DVMSpeciesData%R_S*Temp)**(-3./2.)

skew = SIGN(1.,skew)*MIN(ABS(skew),0.9*max_skew)

delta = (SIGN(1.,skew)*(2*ABS(skew)/(4-Pi))**(1./3.))/SQRT(2./Pi*(1+(2*ABS(skew)/(4-Pi))**(2./3.)))
alpha = delta/SQRT(1.-delta*delta)
omega = SQRT(DVMSpeciesData%R_S*Temp/(1.-2*delta*delta/Pi))
ksi = uVelo - SQRT(2/Pi)*omega*delta

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  cVel(1) = (DVMVelos(iVel,1) - ksi(1))/omega(1)
  cVel(2) = (DVMVelos(jVel,2) - ksi(2))/omega(2)
  cVel(3) = (DVMVelos(kVel,3) - ksi(3))/omega(3)
  cMag2 = DOT_PRODUCT(cVel,cVel)
  Phi = 1.+ERF(alpha*cVel/sqrt(2.))

  fSkew(upos) = rho*Phi(1)*Phi(2)*Phi(3)*EXP(-cMag2/2.)/PRODUCT(omega(1:DVMDim))/(2.*Pi)**(DVMDim/2.)

  IF (DVMDim.LT.3) THEN
    fSkew(PP_nVar_FV/2+upos) = fSkew(upos)*DVMSpeciesData%R_S*Temp*(DVMSpeciesData%Internal_DOF+3.-DVMDim)
  END IF

END DO; END DO; END DO

END SUBROUTINE SkewNormalDistribution

SUBROUTINE TargetDistribution(MacroVal,fTarget)
!===================================================================================================================================
! Target distribution from macro values
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMBGKModel
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: fTarget(PP_nVar_FV)
REAL, INTENT(IN)                 :: MacroVal(14)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SELECT CASE(DVMBGKModel)
  CASE(1)
    CALL ESBGKDistribution(MacroVal,fTarget)
  CASE(2)
    CALL ShakhovDistribution(MacroVal,fTarget)
  CASE(3)
    CALL MaxwellDistribution(MacroVal,fTarget)
  CASE(4)
    CALL ESBGKDistributionCons(MacroVal,fTarget)
  CASE(5)
    CALL MaxwellDistributionCons(MacroVal,fTarget)
  CASE(6)
    CALL SkewNormalDistribution(MacroVal,fTarget)
  CASE(7)
    CALL GradDistributionPrandtl(MacroVal,fTarget)
  CASE DEFAULT
    CALL abort(__STAMP__,'DVM BGK Model not implemented.')
END SELECT

END SUBROUTINE TargetDistribution

SUBROUTINE MaxwellScattering(fBoundary,U,NormVec,tilde,tDeriv)
!===================================================================================================================================
! Gets accurate density for the half maxwellian at diffusive boundaries
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMWeights, DVMMethod
USE MOD_PreProc
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: fBoundary(PP_nVar_FV)
REAL,INTENT(IN)                      :: U(PP_nVar_FV), NormVec(3), tDeriv
INTEGER,INTENT(IN)                   :: tilde
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: Sout, Sin, weight, tau, prefac
REAL                                 :: fTarget(PP_nVar_FV), Urelaxed(PP_nVar_FV), MacroVal(14), vnormal
INTEGER                              :: iVel, jVel, kVel, upos
!===================================================================================================================================
CALL MacroValuesFromDistribution(MacroVal,U,tDeriv,tau,tilde)
IF (tDeriv.EQ.0.) THEN
  prefac = 1.
ELSE
SELECT CASE(tilde)
  CASE(1)
    SELECT CASE(DVMMethod)
    CASE(1)
      prefac = tau*(1.-EXP(-tDeriv/tau))/tDeriv ! f from f2~
    CASE(2)
      prefac = 2.*tau/(2.*tau+tDeriv)
    END SELECT
  CASE(2)
    SELECT CASE(DVMMethod)
    CASE(1)
      prefac = 1 !tau*(EXP(tDeriv/tau)-1.)/tDeriv ! f from f2^ (currently f=f2^: no relaxation to f in the boundary grad calculation)
    CASE(2)
      prefac = 2.*tau/(2.*tau-tDeriv)
    END SELECT
END SELECT
END IF

CALL TargetDistribution(MacroVal, fTarget)

Urelaxed = U*prefac + ftarget*(1.-prefac)

Sin = 0.
Sout = 0.

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  vnormal = DVMVelos(iVel,1)*NormVec(1) + DVMVelos(jVel,2)*NormVec(2) + DVMVelos(kVel,3)*NormVec(3)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  IF (vnormal.GT.0.) THEN !outflow
    Sout = Sout + weight*vnormal*Urelaxed(upos)
  ELSE !inflow
    Sin = Sin - weight*vnormal*fBoundary(upos)
  END IF
END DO; END DO; END DO

fBoundary = fBoundary * (Sout/Sin)
! no additional rescaling needed because it is an equilibrium distribution

END SUBROUTINE MaxwellScattering

SUBROUTINE RescaleU(tilde,tDeriv)
!===================================================================================================================================
! Rescales distribution function for EDDVM/DUGKS
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY : DVMMomentSave, DVMMethod, DVMSpeciesData
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
REAL                            :: MacroVal(14), tau, fTarget(PP_nVar_FV), prefac
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================
DO iElem =1, nElems
  DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
    SELECT CASE (tilde)
      CASE(1) ! f~  -----> f2^    (tDeriv=dt)
        CALL MacroValuesFromDistribution(MacroVal(:),U_FV(:,i,j,k,iElem),tDeriv,tau,tilde)
        DVMMomentSave(1:14,iElem) = MacroVal(1:14)
        DVMMomentSave(15,iElem) = tau
        SELECT CASE(DVMMethod)
        CASE(1)
          prefac = (EXP(-tDeriv/tau/2.)-EXP(-3.*tDeriv/tau/2.))/(1.-EXP(-tDeriv/tau/2.))/2.
        CASE(2)
          prefac = (2.*tau-tDeriv/2.)/(2.*tau+tDeriv)
        END SELECT
      CASE(2) ! f2^ -----> f^     (tDeriv=dt/2)
        MacroVal(1:14) = DVMMomentSave(1:14,iElem)
        tau = DVMMomentSave(15,iElem)
        SELECT CASE(DVMMethod)
        CASE(1)
          prefac = 2.*(EXP(-tDeriv/tau)-EXP(-2.*tDeriv/tau))/(1.-EXP(-2.*tDeriv/tau))
        CASE(2)
          prefac = (4./3.)-(1./3.)*(2.*tau+2.*tDeriv)/(2.*tau-tDeriv)
        END SELECT
    END SELECT
    ! IF (MacroVal(5).LE.0) print*, iElem, i,j,k
    CALL TargetDistribution(MacroVal, fTarget)
    U_FV(:,i,j,k,iElem) = U_FV(:,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
  END DO; END DO; END DO
END DO
END SUBROUTINE RescaleU

SUBROUTINE RescaleInit(tDeriv)
!===================================================================================================================================
! Initial rescale (f->ftilde) for initialization with non equilibrium flow
! TODO: Should also be used for restart
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY: DVMMethod
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
REAL                            :: MacroVal(14), tau, fTarget(PP_nVar_FV), prefac
INTEGER                         :: i,j,k,iElem
!===================================================================================================================================
SWRITE(UNIT_stdOut,*) 'INITIAL DISTRIBUTION FUNCTION RESCALE'
DO iElem =1, nElems
  DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
    CALL MacroValuesFromDistribution(MacroVal(:),U_FV(:,i,j,k,iElem),0.,tau,1) ! tDeriv=0 to get heatflux from original distribution
    CALL TargetDistribution(MacroVal, fTarget)
    SELECT CASE (DVMMethod)
    CASE(1)
      prefac = (tDeriv/tau)/(1. - (EXP(-tDeriv/tau)))
    CASE(2)
      prefac = (2.*tau+tDeriv)/(2.*tau)
    END SELECT
    U_FV(:,i,j,k,iElem) = U_FV(:,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
  END DO; END DO; END DO
END DO
END SUBROUTINE RescaleInit

SUBROUTINE ForceStep(tDeriv)
!===================================================================================================================================
! Calculates force term (to add in 2 parts (Strang splitting) for 2nd order accuracy)
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMForce !, DVMBGKModel
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
REAL                            :: MacroVal(14), tau, fTarget(PP_nVar_FV), forceTerm, cVel(3)!, velodiff, gamma
INTEGER                         :: i,j,k,iElem,iVel,jVel,kVel,upos !,upos1,upos2
!===================================================================================================================================
DO iElem =1, nElems
  DO k=0, PP_N; DO j=0, PP_N; DO i=0, PP_N
    CALL MacroValuesFromDistribution(MacroVal(:),U_FV(:,i,j,k,iElem),tDeriv,tau,1)
    ! SELECT CASE (DVMBGKModel)
    !   CASE(1)
    !     CALL MaxwellDistribution(MacroVal,fTarget)
    !   CASE(2)
    !     CALL ShakhovDistribution(MacroVal,fTarget)
    !   CASE DEFAULT
    !     CALL abort(__STAMP__,'DVM BGK Model not implemented.')
    !   END SELECT
    ! gamma = tau*(1.-EXP(-tDeriv/tau))/tDeriv

    CALL MaxwellDistribution(MacroVal,fTarget) !equilibrium approximation

    DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
      upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)

      !equilibrium approximation
      cVel(1) = DVMVelos(iVel,1) - MacroVal(2)
      cVel(2) = DVMVelos(jVel,2) - MacroVal(3)
      cVel(3) = DVMVelos(kVel,3) - MacroVal(4)
      forceTerm = DOT_PRODUCT(DVMForce,cVel)/(DVMSpeciesData%R_S*MacroVal(5)) * fTarget(upos)

      ! non equilibrium version
      ! IF (iVel.EQ.1) THEN
      !   upos1=upos
      !   upos2 = iVel+1+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
      !   velodiff=DVMVelos(iVel+1,1)-DVMVelos(iVel,1)
      ! ELSE IF (iVel.EQ.DVMnVelos(1)) THEN
      !   upos1 = iVel-1+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
      !   upos2=upos
      !   velodiff=DVMVelos(iVel,1)-DVMVelos(iVel-1,1)
      ! ELSE
      !   upos1 = iVel-1+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
      !   upos2 = iVel+1+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
      !   velodiff=DVMVelos(iVel+1,1)-DVMVelos(iVel-1,1)
      ! END IF
      ! forceTerm = - DVMForce(1)*(gamma*(U(upos2,i,j,k,iElem)-U(upos1,i,j,k,iElem)) &
      !                        +(1-gamma)*(fTarget(upos2)-fTarget(upos1)))/velodiff

      U_FV(upos,i,j,k,iElem) = U_FV(upos,i,j,k,iElem) + forceTerm*tDeriv/2 !t/2 for strang splitting
    END DO; END DO; END DO
  END DO; END DO; END DO
END DO
END SUBROUTINE ForceStep


SUBROUTINE IntegrateFluxValues(MacroVal,U)
!===================================================================================================================================
! Calculates the surface macro values from distribution fluxes
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMWeights, DVMDim
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
REAL                            :: rho, rhoU(3), rhoE, weight
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = 0.
rhoU = 0.
rhoE = 0.
MacroVal = 0.

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  rho = rho + weight*U(upos)
  rhoU(1) = rhoU(1) + weight*DVMVelos(iVel,1)*U(upos)
  rhoU(2) = rhoU(2) + weight*DVMVelos(jVel,2)*U(upos)
  rhoU(3) = rhoU(3) + weight*DVMVelos(kVel,3)*U(upos)
  IF (DVMDim.LT.3) THEN
    rhoE = rhoE + weight*0.5*((DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)+U(PP_nVar_FV/2+upos))
  ELSE
    rhoE = rhoE + weight*0.5*(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)
  END IF
END DO; END DO; END DO

MacroVal(1) = rho ! mass flow
MacroVal(2:4) = rhoU ! force per area
MacroVal(5) = rhoE ! heat flux

END SUBROUTINE IntegrateFluxValues

END MODULE MOD_DistFunc