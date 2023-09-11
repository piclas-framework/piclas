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

PUBLIC:: MacroValuesFromDistribution
PUBLIC:: MaxwellDistribution, MaxwellDistributionCons, ShakhovDistribution, ESBGKDistribution, GradDistribution
PUBLIC:: MaxwellScattering, RescaleU, RescaleInit, ForceStep
!==================================================================================================================================

CONTAINS

SUBROUTINE MacroValuesFromDistribution(MacroVal,U,tDeriv,tau,tilde)
!===================================================================================================================================
! Calculates the moments of the distribution function
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMWeights, DVMSpeciesData, DVMMethod, DVMBGKModel
USE MOD_PreProc
USE MOD_Globals               ,ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: U(PP_nVar_FV), tDeriv
INTEGER,INTENT(IN)              :: tilde
REAL, INTENT(OUT)               :: MacroVal(14), tau
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: rho, rhoU(3), rhoE, PressTens(6), Heatflux(3), uVelo(3), cV, cVel(3),cMag, mu, weight
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = 0.
rhoU = 0.
rhoE = 0.
PressTens = 0.
Heatflux = 0.
MacroVal = 0.
DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  rho = rho + weight*U(upos)
  rhoU(1) = rhoU(1) + weight*DVMVelos(iVel,1)*U(upos)
  rhoU(2) = rhoU(2) + weight*DVMVelos(jVel,2)*U(upos)
  rhoU(3) = rhoU(3) + weight*DVMVelos(kVel,3)*U(upos)
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    rhoE = rhoE + weight*0.5*((DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)+U(PP_nVar_FV/2+upos))
  ELSE
    rhoE = rhoE + weight*0.5*(DVMVelos(iVel,1)**2.+DVMVelos(jVel,2)**2.+DVMVelos(kVel,3)**2.)*U(upos)
  END IF
END DO; END DO; END DO

uVelo = rhoU/rho
! cV = (3.+DVMSpeciesData%Internal_DOF)/2.*DVMSpeciesData%R_S*rho
cV = 3./2.*DVMSpeciesData%R_S*rho

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag = DOT_PRODUCT(cVel,cVel)
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    Heatflux(1) = Heatflux(1) + weight*0.5*cVel(1)*(U(upos)*cMag+U(PP_nVar_FV/2+upos))
    Heatflux(2) = Heatflux(2) + weight*0.5*cVel(2)*(U(upos)*cMag+U(PP_nVar_FV/2+upos))
    Heatflux(3) = Heatflux(3) + weight*0.5*cVel(3)*(U(upos)*cMag+U(PP_nVar_FV/2+upos))
  ELSE
    PressTens(1) = PressTens(1) + weight*cVel(1)*cVel(1)*U(upos)
    PressTens(2) = PressTens(2) + weight*cVel(1)*cVel(2)*U(upos)
    PressTens(3) = PressTens(3) + weight*cVel(1)*cVel(3)*U(upos)
    PressTens(4) = PressTens(4) + weight*cVel(2)*cVel(2)*U(upos)
    PressTens(5) = PressTens(5) + weight*cVel(2)*cVel(3)*U(upos)
    PressTens(6) = PressTens(6) + weight*cVel(3)*cVel(3)*U(upos)
    Heatflux(1) = Heatflux(1) + weight*0.5*cVel(1)*(U(upos)*cMag)
    Heatflux(2) = Heatflux(2) + weight*0.5*cVel(2)*(U(upos)*cMag)
    Heatflux(3) = Heatflux(3) + weight*0.5*cVel(3)*(U(upos)*cMag)
  END IF
END DO; END DO; END DO

MacroVal(1) = rho
MacroVal(2:4) = uVelo
MacroVal(5) = (rhoE - 0.5*(DOT_PRODUCT(rhoU,rhoU))/rho)/cV
IF (MacroVal(5).LE.0) CALL abort(__STAMP__,'DVM negative temperature!')

mu = DVMSpeciesData%mu_Ref*(MacroVal(5)/DVMSpeciesData%T_Ref)**(DVMSpeciesData%omegaVHS+0.5)
tau = mu/(DVMSpeciesData%R_S*MacroVal(1)*MacroVal(5))
IF (DVMBGKModel.EQ.1) tau = tau/DVMSpeciesData%Prandtl !ESBGK

IF (tDeriv.EQ.0.) THEN
  Macroval(6:11)  = PressTens(1:6)
  MacroVal(12:14) = Heatflux(1:3)
ELSE
  SELECT CASE (tilde)
    CASE(1) ! higher moments from f~
      SELECT CASE(DVMMethod)
      CASE(1)
        Macroval(6:11)  = PressTens(1:6)*(1.-EXP(-tDeriv/tau))/(tDeriv/tau)
        MacroVal(12:14) = Heatflux(1:3)*(1.-EXP(-tDeriv*DVMSpeciesData%Prandtl/tau))/(tDeriv*DVMSpeciesData%Prandtl/tau)
      CASE(2)
        Macroval(6:11)  = PressTens(1:6)*2.*tau/(2.*tau+tDeriv)
        MacroVal(12:14) = Heatflux(1:3)*2.*tau/(2.*tau+tDeriv*DVMSpeciesData%Prandtl)
      END SELECT
    CASE(2) ! higher moments from f^
      ! MacroVal(12:14) = Heatflux(1:3)*(1.-EXP(-DVMSpeciesData%Prandtl*tDeriv/tau)) &
      !                            /(EXP(-DVMSpeciesData%Prandtl*tDeriv/tau)*tDeriv*DVMSpeciesData%Prandtl/tau)
      MacroVal(6:14) = 0. !will get copied from earlier f~ macroval
  END SELECT
END IF

END SUBROUTINE

SUBROUTINE MaxwellDistributionCons(MacroVal,fMaxwell)
!===================================================================================================================================
! conservative maxwell
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
REAL                            :: rho, Temp, uVelo(DVMDim), vMag, gM, weight
REAL,DIMENSION(2+DVMDim)        :: alpha, psi, rhovec
REAL                            :: J(2+DVMDim,2+DVMDim), B(2+DVMDim,1)
INTEGER                         :: iVel,jVel,kVel, upos, countz, IPIV(2+DVMDim), info_dgesv
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:DVMDim) = MacroVal(2:1+DVMDim)
Temp = MacroVal(5)
rhovec(1) = rho
rhovec(2:1+DVMDim)=rho*uVelo(:)
rhovec(2+DVMDim)=rho*(3*DVMSpeciesData%R_S*Temp+DOT_PRODUCT(uVelo,uVelo))/2.
countz=0

alpha(1) = LOG(rho/(2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))-DOT_PRODUCT(uVelo,uVelo)/2./DVMSpeciesData%R_S/Temp
alpha(2:1+DVMDim) = uVelo(1:DVMDim)/DVMSpeciesData%R_S/Temp
alpha(2+DVMDim) = -1/DVMSpeciesData%R_S/Temp

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
    J(:,1) = J(:,1)+weight*gM*psi(:)
    J(:,2) = J(:,2)+weight*gM*psi(2)*psi(:)
    IF (DVMDim.GT.1) J(:,3) = J(:,3)+weight*gM*psi(3)*psi(:)
    IF (DVMDim.GT.2) J(:,4) = J(:,4)+weight*gM*psi(4)*psi(:)
    J(:,2+DVMDim) = J(:,2+DVMDim)+weight*gM*psi(2+DVMDim)*psi(:)
  END DO; END DO; END DO

  B(:,1) = rhovec(:) - J(:,1)

  CALL DGESV(2+DVMDim,1,J,2+DVMDim,IPIV,B,2+DVMDim,info_dgesv)

  IF(info_dgesv.NE.0) CALL abort(__STAMP__,'Newton DGESV fail')
  IF (.NOT.(ANY(ABS(B(:,1)).GT.(1e-5*ABS(alpha)+1e-12)))) EXIT
  ! IF (.NOT.(ANY(ABS(B(:,1)).GT.(1e-6*ABS(alpha)+1e-16)))) EXIT
  alpha = alpha + B(:,1)
END DO
IF (countz.GE.1000) print*, 'Newton max iter reached'

J=0.

alpha = alpha + B(:,1)
DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  weight = DVMWeights(iVel,1)*DVMWeights(jVel,2)*DVMWeights(kVel,3)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  vMag = DVMVelos(iVel,1)**2 + DVMVelos(jVel,2)**2 + DVMVelos(kVel,3)**2
  psi(1)=1
  psi(2)=DVMVelos(iVel,1)
  IF (DVMDim.GT.1) psi(3)=DVMVelos(jVel,2)
  IF (DVMDim.GT.2) psi(4)=DVMVelos(kVel,3)
  psi(2+DVMDim) = vMag/2.
  gM = EXP(DOT_PRODUCT(alpha,psi))
  fMaxwell(upos)= gM
  J(:,1) = J(:,1)+weight*gM*psi(:)
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    fMaxwell(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF
  END IF
END DO; END DO; END DO

END SUBROUTINE

SUBROUTINE MaxwellDistribution(MacroVal,fMaxwell)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
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
REAL                            :: rho, Temp, uVelo(3), cVel(3), cMag, gM
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
  cMag = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  gM = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag/(2.*DVMSpeciesData%R_S*Temp))
  fMaxwell(upos)= gM
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    fMaxwell(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF
  END IF
END DO; END DO; END DO

END SUBROUTINE

SUBROUTINE ShakhovDistribution(MacroVal,fShakhov)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
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
REAL                            :: rho, Temp, uVelo(3), cVel(3), cMag, gM, q(3), Prandtl, ShakhFac1, ShakhFac2
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho = MacroVal(1)
uVelo(1:3) = MacroVal(2:4)
Temp = MacroVal(5)
q(1:3) = MacroVal(12:14)
! Prandtl = 2.*(DVMSpeciesData%Internal_DOF + 5.)/(2.*DVMSpeciesData%Internal_DOF + 15.)
Prandtl = 2./3.

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)
  cMag = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)
  gM = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag/(2.*DVMSpeciesData%R_S*Temp))
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*rho*DVMSpeciesData%R_S*DVMSpeciesData%R_S*Temp*Temp)
  ShakhFac2 = cMag/(DVMSpeciesData%R_S*Temp)
  fShakhov(upos) = gm*(1.+(1.-Prandtl)*ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    fShakhov(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF*(1. &
        +(1.-Prandtl)*ShakhFac1*((ShakhFac2-DVMDim)-2.*(DVMSpeciesData%Internal_DOF-3.+DVMDim)))
  END IF
END DO; END DO; END DO

END SUBROUTINE

SUBROUTINE ESBGKDistribution(MacroVal,fESBGK)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
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
REAL                            :: rho,Temp,uVelo(3),cVel(3),pressTens(3,3),pressProduct,ilambda(3,3),ldet
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho              = MacroVal(1)
uVelo(1:3)       = MacroVal(2:4)
Temp             = MacroVal(5)
pressTens        = 0.
pressTens(2:3,1) = MacroVal(7:8)
pressTens(1,2)   = MacroVal(7)
pressTens(3,2)   = MacroVal(10)
pressTens(1,3)   = MacroVal(8)
pressTens(2,3)   = MacroVal(10)

pressTens = (1.-1./DVMSpeciesData%Prandtl)*pressTens/rho
pressTens(1,1) = DVMSpeciesData%R_S*Temp
pressTens(2,2) = DVMSpeciesData%R_S*Temp
pressTens(3,3) = DVMSpeciesData%R_S*Temp

CALL INV33(pressTens,ilambda,ldet)
IF (ldet.LE.0.) CALL abort(__STAMP__,'DVM negative temperature!')

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)

  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)

  pressProduct = cVel(1)*DOT_PRODUCT(ilambda(:,1),cVel) &
               + cVel(2)*DOT_PRODUCT(ilambda(:,2),cVel) &
               + cVel(3)*DOT_PRODUCT(ilambda(:,3),cVel)

  ! print*, ldet, pressProduct

  fESBGK(upos) = rho/sqrt(ldet*(2*Pi)**3)*EXP(-pressProduct/2.)
  IF ((DVMSpeciesData%Internal_DOF .GT.0.0).OR.(DVMDim.LT.3)) THEN
    CALL abort(__STAMP__,'DVM ESBGK model implemented only for 3D monatomic')
  END IF
END DO; END DO; END DO

END SUBROUTINE

SUBROUTINE GradDistribution(MacroVal,fGrad)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMSpeciesData, DVMDim, Pi
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
REAL                            :: rho,Temp,uVelo(3),cVel(3),cMag,gM,q(3),pressTens(3,3),ShakhFac1,ShakhFac2,pressFac,pressProduct
INTEGER                         :: iVel,jVel,kVel, upos
!===================================================================================================================================
rho              = MacroVal(1)
uVelo(1:3)       = MacroVal(2:4)
Temp             = MacroVal(5)
q(1:3)           = MacroVal(12:14)
pressTens        = 0.
pressTens(2:3,1) = MacroVal(7:8)
pressTens(1,2)   = MacroVal(7)
pressTens(3,2)   = MacroVal(10)
pressTens(1,3)   = MacroVal(8)
pressTens(2,3)   = MacroVal(10)

DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
  upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)

  cVel(1) = DVMVelos(iVel,1) - uVelo(1)
  cVel(2) = DVMVelos(jVel,2) - uVelo(2)
  cVel(3) = DVMVelos(kVel,3) - uVelo(3)

  pressProduct = cVel(1)*DOT_PRODUCT(pressTens(:,1),cVel) &
               + cVel(2)*DOT_PRODUCT(pressTens(:,2),cVel) &
               + cVel(3)*DOT_PRODUCT(pressTens(:,3),cVel)

  cMag = cVel(1)*cVel(1) + cVel(2)*cVel(2)+ cVel(3)*cVel(3)

  pressFac = rho*DVMSpeciesData%R_S*DVMSpeciesData%R_S*Temp*Temp
  ShakhFac1 = DOT_PRODUCT(q,cVel)/(5.*pressFac)
  ShakhFac2 = cMag/(DVMSpeciesData%R_S*Temp)

  gM = rho/((2.*Pi*DVMSpeciesData%R_S*Temp)**(DVMDim/2.))*exp(-cMag/(2.*DVMSpeciesData%R_S*Temp))
  fGrad(upos) = gM*(1.+0.5*pressProduct/pressFac+ShakhFac1*(ShakhFac2-2.-DVMDim))
  IF (DVMSpeciesData%Internal_DOF .GT.0.0) THEN
    fGrad(PP_nVar_FV/2+upos) = gM*DVMSpeciesData%R_S*Temp*DVMSpeciesData%Internal_DOF*(1. &
        +ShakhFac1*((ShakhFac2-DVMDim)-2.*(DVMSpeciesData%Internal_DOF-3.+DVMDim)))
  END IF
END DO; END DO; END DO

END SUBROUTINE

SUBROUTINE MaxwellScattering(fBoundary,U,NormVec,tilde,tDeriv)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMnVelos, DVMVelos, DVMWeights, DVMBGKModel, DVMMethod
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
SELECT CASE (DVMBGKModel)
  CASE(1)
    CALL ESBGKDistribution(MacroVal,fTarget)
  CASE(2)
    CALL ShakhovDistribution(MacroVal,fTarget)
  CASE(3)
    CALL MaxwellDistribution(MacroVal,fTarget)
  CASE(4)
    CALL MaxwellDistributionCons(MacroVal,fTarget)
  CASE DEFAULT
    CALL abort(__STAMP__,'DVM BGK Model not implemented.')
END SELECT

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

END SUBROUTINE

SUBROUTINE RescaleU(tilde,tDeriv)
!===================================================================================================================================
! Specifies all the initial conditions. The state in conservative variables is returned.
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY : DVMBGKModel, DVMMomentSave, DVMMethod
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
    CALL MacroValuesFromDistribution(MacroVal(:),U_FV(:,i,j,k,iElem),tDeriv,tau,tilde)
    SELECT CASE (tilde)
      CASE(1) ! f~  -----> f2^    (tDeriv=dt)
        DVMMomentSave(1:9,iElem) = MacroVal(6:14)
        SELECT CASE(DVMMethod)
        CASE(1)
          prefac = (EXP(-tDeriv/tau/2.)-EXP(-3.*tDeriv/tau/2.))/(1.-EXP(-tDeriv/tau/2.))/2.
        CASE(2)
          prefac = (2.*tau-tDeriv/2.)/(2.*tau+tDeriv)
        END SELECT
      CASE(2) ! f2^ -----> f^     (tDeriv=dt/2)
        MacroVal(6:14) = DVMMomentSave(1:9,iElem)
        SELECT CASE(DVMMethod)
        CASE(1)
          prefac = 2.*(EXP(-tDeriv/tau)-EXP(-2.*tDeriv/tau))/(1.-EXP(-2.*tDeriv/tau))
        CASE(2)
          prefac = (4./3.)-(1./3.)*(2.*tau+2.*tDeriv)/(2.*tau-tDeriv)
        END SELECT
    END SELECT
    ! IF (MacroVal(5).LE.0) print*, iElem, i,j,k
    SELECT CASE (DVMBGKModel)
      CASE(1)
        CALL ESBGKDistribution(MacroVal,fTarget)
      CASE(2)
        CALL ShakhovDistribution(MacroVal,fTarget)
      CASE(3)
        CALL MaxwellDistribution(MacroVal,fTarget)
      CASE(4)
        CALL MaxwellDistributionCons(MacroVal,fTarget)
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM BGK Model not implemented.')
    END SELECT
    U_FV(:,i,j,k,iElem) = U_FV(:,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
  END DO; END DO; END DO
END DO
END SUBROUTINE

SUBROUTINE RescaleInit(tDeriv)
!===================================================================================================================================
! Initial rescale for initialization with non equilibrium flow
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV,  ONLY: DVMBGKModel, DVMMethod
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
    SELECT CASE (DVMBGKModel)
      CASE(1)
        CALL ESBGKDistribution(MacroVal,fTarget)
      CASE(2)
        CALL ShakhovDistribution(MacroVal,fTarget)
      CASE(3)
        CALL MaxwellDistribution(MacroVal,fTarget)
      CASE(4)
        CALL MaxwellDistributionCons(MacroVal,fTarget)
      CASE DEFAULT
        CALL abort(__STAMP__,'DVM BGK Model not implemented.')
    END SELECT
    SELECT CASE (DVMMethod)
    CASE(1)
      prefac = (tDeriv/tau)/(1. - (EXP(-tDeriv/tau)))
    CASE(2)
      prefac = (2.*tau+tDeriv)/(2.*tau)
    END SELECT
    U_FV(:,i,j,k,iElem) = U_FV(:,i,j,k,iElem)*prefac + fTarget(:)*(1.-prefac)
  END DO; END DO; END DO
END DO
END SUBROUTINE

SUBROUTINE ForceStep(tDeriv)
!===================================================================================================================================
! Initial rescale for initialization with non equilibrium flow
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
END SUBROUTINE

END MODULE MOD_DistFunc
