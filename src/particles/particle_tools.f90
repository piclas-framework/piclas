!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_part_tools
!===================================================================================================================================
! Contains tools for particle related operations. This routine is uses MOD_Particle_Boundary_Tools, but not vice versa!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

INTERFACE VeloFromDistribution
  MODULE PROCEDURE VeloFromDistribution
END INTERFACE

INTERFACE LIQUIDEVAP
  MODULE PROCEDURE LIQUIDEVAP
END INTERFACE

INTERFACE LIQUIDREFL
  MODULE PROCEDURE LIQUIDREFL
END INTERFACE

INTERFACE ALPHALIQUID
  MODULE PROCEDURE ALPHALIQUID
END INTERFACE

INTERFACE BETALIQUID
  MODULE PROCEDURE BETALIQUID
END INTERFACE

INTERFACE TSURUTACONDENSCOEFF
  MODULE PROCEDURE TSURUTACONDENSCOEFF
END INTERFACE

INTERFACE isChargedParticle
  MODULE PROCEDURE isChargedParticle
END INTERFACE

INTERFACE isPushParticle
  MODULE PROCEDURE isPushParticle
END INTERFACE

INTERFACE isDepositParticle
  MODULE PROCEDURE isDepositParticle
END INTERFACE

INTERFACE isInterpolateParticle
  MODULE PROCEDURE isInterpolateParticle
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: LIQUIDEVAP,LIQUIDREFL,ALPHALIQUID,BETALIQUID,TSURUTACONDENSCOEFF
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, VeloFromDistribution, GetParticleWeight, isChargedParticle
PUBLIC :: isPushParticle, isDepositParticle, isInterpolateParticle
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars        ,ONLY: PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem, PartPressureCell
USE MOD_Particle_Vars        ,ONLY: KeepWallParticles, PartState, VarTimeStep
USE MOD_DSMC_Vars            ,ONLY: useDSMC, CollInf
USE MOD_Particle_VarTimeStep ,ONLY: CalcVarTimeStep
#if USE_MPI
USE MOD_MPI_Vars             ,ONLY: OffSetElemMPI
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers   ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: counter1,i,n
INTEGER            :: ElemID
#if !USE_MPI
INTEGER            :: OffSetElemMPI(0) = 0            !> Dummy offset for single-thread mode
#endif
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF(PDM%maxParticleNumber.EQ.0) RETURN
counter1 = 1
IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell) THEN
  PEM%pNumber(:) = 0
  IF (KeepWallParticles) PEM%wNumber(:) = 0
END IF

n = PDM%ParticleVecLength !PDM%maxParticleNumber
PDM%ParticleVecLength = 0
PDM%insideParticleNumber = 0
IF (doParticleMerge) vMPF_SpecNumElem = 0

IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell) THEN
  DO i=1,n
    IF (.NOT.PDM%ParticleInside(i)) THEN
      IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
      PDM%nextFreePosition(counter1) = i
      counter1 = counter1 + 1
    ELSE
      ElemID = PEM%Element(i) - OffSetElemMPI(MyRank)
      IF (PEM%pNumber(ElemID).EQ.0) THEN
        PEM%pStart(ElemID) = i                     ! Start of Linked List for Particles in Elem
      ELSE
        PEM%pNext(PEM%pEnd(ElemID)) = i            ! Next Particle of same Elem (Linked List)
      END IF
      PEM%pEnd(ElemID) = i
      PEM%pNumber(ElemID) = &                      ! Number of Particles in Element
          PEM%pNumber(ElemID) + 1
      IF (VarTimeStep%UseVariableTimeStep) THEN
        VarTimeStep%ParticleTimeStep(i) = CalcVarTimeStep(PartState(1,i),PartState(2,i),ElemID)
      END IF
      IF (KeepWallParticles) THEN
        IF (PDM%ParticleAtWall(i)) THEN
          PEM%wNumber(ElemID) = PEM%wNumber(ElemID) + 1
        END IF
      END IF
      PDM%ParticleVecLength = i
      IF(doParticleMerge) vMPF_SpecNumElem(ElemID,PartSpecies(i)) = vMPF_SpecNumElem(ElemID,PartSpecies(i)) + 1
    END IF
  END DO
ELSE ! no DSMC
  DO i=1,n
    IF (.NOT.PDM%ParticleInside(i)) THEN
      PDM%nextFreePosition(counter1) = i
      counter1 = counter1 + 1
    ELSE
      PDM%ParticleVecLength = i
    END IF
  END DO
ENDIF
PDM%insideParticleNumber = PDM%ParticleVecLength - counter1+1
PDM%CurrentNextFreePosition = 0
DO i = n+1,PDM%maxParticleNumber
  IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
  PDM%nextFreePosition(counter1) = i
  counter1 = counter1 + 1
END DO
PDM%nextFreePosition(counter1:PDM%MaxParticleNumber)=0 ! exists if MaxParticleNumber is reached!!!
IF (counter1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber)=0

#if USE_LOADBALANCE
CALL LBPauseTime(LB_UNFP,tLBStart)
#endif /*USE_LOADBALANCE*/

  RETURN
END SUBROUTINE UpdateNextFreePosition

FUNCTION DiceUnitVector()
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals_Vars ,ONLY: Pi
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                     :: DiceUnitVector(3)
REAL                     :: iRan, bVec, aVec
!===================================================================================================================================
CALL RANDOM_NUMBER(iRan)
bVec              = 1. - 2.*iRan
aVec              = SQRT(1. - bVec**2.)
DiceUnitVector(3) = bVec
CALL RANDOM_NUMBER(iRan)
bVec              = Pi *2. * iRan
DiceUnitVector(1) = aVec * COS(bVec)
DiceUnitVector(2) = aVec * SIN(bVec)

END FUNCTION DiceUnitVector


FUNCTION VeloFromDistribution(distribution,specID,Tempergy)
!===================================================================================================================================
!> calculation of velocityvector (Vx,Vy,Vz) sampled from given distribution function
!>  liquid_evap: normal direction to surface with ARM from shifted evaporation rayleigh, tangential from normal distribution
!>  liquid_refl: normal direction to surface with ARM from shifted reflection rayleigh, tangential from normal distribution
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals                 ,ONLY: Abort,UNIT_stdOut
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_Particle_Vars           ,ONLY: Species
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: distribution !< specifying keyword for velocity distribution
INTEGER,INTENT(IN)          :: specID       !< input species
REAL,INTENT(IN)             :: Tempergy         !< input temperature [K] or energy [J] or velocity [m/s]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, PARAMETER :: xmin=0., xmax=5.
REAL            :: VeloFromDistribution(1:3)
REAL            :: alpha, beta
REAL            :: y1, f, ymax, i, binsize
REAL            :: sigma, val(1:2)
REAL            :: Velo1, Velo2, Velosq
REAL            :: RandVal(2)
!===================================================================================================================================
!-- set velocities
SELECT CASE(TRIM(distribution))
CASE('deltadistribution')
  !Velosq = 2
  !DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    !CALL RANDOM_NUMBER(RandVal)
    !Velo1 = 2.*RandVal(1) - 1.
    !Velo2 = 2.*RandVal(2) - 1.
    !Velosq = Velo1**2 + Velo2**2
  !END DO
  !VeloFromDistribution(1) = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
  !VeloFromDistribution(2) = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
  !CALL RANDOM_NUMBER(RandVal)
  !VeloFromDistribution(3) = SQRT(-2*LOG(RandVal(1)))

  ! Get random vector
  VeloFromDistribution = DiceUnitVector()
  ! Mirror z-component of velocity (particles are emitted from surface!)
  VeloFromDistribution(3) = ABS(VeloFromDistribution(3))
  ! Set magnitude
  VeloFromDistribution = Tempergy*VeloFromDistribution

CASE('liquid_evap','liquid_refl')
  ! sample normal direction with ARM from given, shifted rayleigh distribution function
  sigma = SQRT(BoltzmannConst*Tempergy/Species(SpecID)%MassIC)
  alpha=ALPHALIQUID(specID,Tempergy)
  beta=BETALIQUID(specID,Tempergy)
  ! define ymax used in ARM
  IF (beta.GE.1 .AND. TRIM(distribution).EQ.'liquid_evap') THEN
    i = xmin
    binsize = (xmax-xmin)/100.
    ymax = LIQUIDEVAP(beta,i,1.)
    DO WHILE (i.LE.xmax) !
      val(1)=LIQUIDEVAP(beta,i,1.)
      val(2)=LIQUIDEVAP(beta,i+binsize,1.)
      IF (val(2).GT.val(1)) THEN
        ymax = val(2)
      END IF
      i=i+binsize
    END DO
    ymax = ymax*1.1
  ELSE
    SELECT CASE(TRIM(distribution))
    CASE('liquid_evap')
      ymax = 0.7
    CASE('liquid_refl')
      ymax = 0.9
    END SELECT
  END IF
  ! do ARM loop
  y1=1
  f=0
  DO WHILE (y1-f.GT.0)
    CALL RANDOM_NUMBER(RandVal)
    Velo1=xmin+RandVal(1)*(xmax-xmin)
    SELECT CASE(TRIM(distribution))
    CASE('liquid_evap')
      f=LIQUIDEVAP(beta,Velo1,1.)
    CASE('liquid_refl')
      f=LIQUIDREFL(alpha,beta,Velo1,1.)
    END SELECT
    y1=ymax*RandVal(2)
  END DO
  VeloFromDistribution(3) = sigma*Velo1
  ! build tangential velocities from gauss (normal) distribution
  Velosq = 2
  DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    CALL RANDOM_NUMBER(RandVal)
    Velo1 = 2.*RandVal(1) - 1.
    Velo2 = 2.*RandVal(2) - 1.
    Velosq = Velo1**2 + Velo2**2
  END DO
  VeloFromDistribution(1) = Velo1*SQRT(-2.*LOG(Velosq)/Velosq)*sigma
  VeloFromDistribution(2) = Velo2*SQRT(-2.*LOG(Velosq)/Velosq)*sigma
CASE DEFAULT
  WRITE (UNIT_stdOut,'(A)') "distribution =", distribution
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')
END SELECT

END FUNCTION VeloFromDistribution


PURE REAL FUNCTION GetParticleWeight(iPart)
!===================================================================================================================================
!> Determines the appropriate particle weighting for the axisymmetric case with radial weighting and the variable time step. For
!> radial weighting, the radial factor is multiplied by the regular weighting factor. If only a variable time step is used, at the
!> moment, the regular weighting factor is not included.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Particle_Vars           ,ONLY: usevMPF, VarTimeStep, PartMPF
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  IF (VarTimeStep%UseVariableTimeStep) THEN
    GetParticleWeight = PartMPF(iPart) * VarTimeStep%ParticleTimeStep(iPart)
  ELSE
    GetParticleWeight = PartMPF(iPart)
  END IF
ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
  GetParticleWeight = VarTimeStep%ParticleTimeStep(iPart)
ELSE
  GetParticleWeight = 1.
END IF

END FUNCTION GetParticleWeight


PURE REAL FUNCTION LIQUIDEVAP(beta,x,sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: beta,x,sigma
REAL            :: betaLoc
!===================================================================================================================================
betaLoc = beta
IF (betaLoc.GE.2.) betaLoc = 2. - 1e-10
IF (betaLoc.LT.0.) betaLoc = 0.

liquidEvap=(1-betaLoc*exp(-0.5*(x/sigma)**2))/(1-betaLoc/2)  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
IF (liquidEvap.LT.0.) liquidEvap = 0.
END FUNCTION


REAL FUNCTION LIQUIDREFL(alpha,beta,x,sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: alpha,beta,x,sigma
REAL            :: betaLoc, alphaLoc
!===================================================================================================================================
betaLoc = beta
IF (betaLoc.GE.2.) betaLoc = 2. - 1e-10
IF (betaLoc.LT.0.) betaLoc = 0.
alphaLoc = alpha
IF (alphaLoc.GT.1.) alphaLoc = 1.
IF (alphaLoc.LT.0.) alphaLoc = 0.

if (alphaLoc.GE.1.) then
  if (betaLoc.LE.0) then
    liquidRefl = x/sigma**2  *  exp(-0.5*(x/sigma)**2)
  else
    liquidRefl = (betaLoc*exp(-0.5*(x/sigma)**2))/(1.-(1.-betaLoc/2.))  *   x/sigma**2  *  exp(-0.5*(x/sigma)**2)
  end if
else
  liquidRefl = (1.-alphaLoc+alphaLoc*betaLoc*exp(-0.5*(x/sigma)**2))/(1.-alphaLoc*(1.-betaLoc/2.)) &
             * x/sigma**2 * exp(-0.5*(x/sigma)**2)
end if

IF (liquidRefl.LT.0.) liquidRefl = 0.
END FUNCTION


PURE FUNCTION ALPHALIQUID(specID,temp) RESULT(alpha)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SpecSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE (SpecSurf(specID)%condensCase)
CASE (1)
  alpha = SpecSurf(specID)%liquidAlpha
CASE (2)
  alpha = exp(-((4-BETALIQUID(specID,temp))/(2*(2-BETALIQUID(specID,temp)))-1))
END SELECT
IF (alpha.GT.1.) alpha = 1.
IF (alpha.LT.0.) alpha = 0.
END FUNCTION


PURE FUNCTION BETALIQUID(specID,temp) RESULT(beta)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars ,ONLY: SpecSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE (SpecSurf(specID)%condensCase)
CASE (1)
  beta = SpecSurf(specID)%liquidBeta
CASE (2)
  beta = SpecSurf(specID)%liquidBetaCoeff(1)*temp**5 &
       + SpecSurf(specID)%liquidBetaCoeff(2)*temp**4 &
       + SpecSurf(specID)%liquidBetaCoeff(3)*temp**3 &
       + SpecSurf(specID)%liquidBetaCoeff(4)*temp**2 &
       + SpecSurf(specID)%liquidBetaCoeff(5)*temp    &
       + SpecSurf(specID)%liquidBetaCoeff(6)
END SELECT
IF (beta.GE.2.) beta = 2. - 1e-10
IF (beta.LT.0.) beta=0.
END FUNCTION


FUNCTION TSURUTACONDENSCOEFF(SpecID,normalVelo,temp) RESULT(sigma)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars  ,ONLY: BoltzmannConst
USE MOD_Particle_Vars ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: specID
REAL,INTENT(IN)    :: normalVelo,temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: sigma
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
sigma = ALPHALIQUID(specID,temp)*(1-BETALIQUID(specID,temp)*exp(-normalVelo**2*Species(specID)%MassIC/(2*Boltzmannconst*temp)))
IF (sigma.LT.0.) sigma = 0.
IF (sigma.GT.1.) sigma = 1.
END FUNCTION

PURE FUNCTION isChargedParticle(iPart)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if particle has charge unequal to zero and return T/F logical.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iPart
LOGICAL             :: isChargedParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0)THEN
  isChargedParticle = .TRUE.
ELSE
  isChargedParticle = .FALSE.
END IF ! ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0
END FUNCTION isChargedParticle


PURE FUNCTION isPushParticle(iPart)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if particle is to be evolved in time by the particle pusher (time integration).
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==300) /*FP-Flow*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#elif (PP_TimeDiscMethod==400) /*BGK*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#else /*all other methods, mainly PIC*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iPart
LOGICAL             :: isPushParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0)THEN
  isPushParticle = .TRUE.
ELSE
  isPushParticle = .FALSE.
END IF ! ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0
END FUNCTION isPushParticle


PURE FUNCTION isDepositParticle(iPart)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if particle is to be deposited on the grid (particle-to-grid coupling).
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==300) /*FP-Flow*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#elif (PP_TimeDiscMethod==400) /*BGK*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#else /*all other methods, mainly PIC*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iPart
LOGICAL             :: isDepositParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0)THEN
  isDepositParticle = .TRUE.
ELSE
  isDepositParticle = .FALSE.
END IF ! ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0
END FUNCTION isDepositParticle


PURE FUNCTION isInterpolateParticle(iPart)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if particle is to be interpolated (field-to-particle coupling), which is required for calculating the acceleration, e.g.,
! due to Lorentz forces at the position of the particle.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==300) /*FP-Flow*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#elif (PP_TimeDiscMethod==400) /*BGK*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species ! Change this when required
#else /*all other methods, mainly PIC*/
USE MOD_Particle_Vars ,ONLY: PartSpecies,Species
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iPart
LOGICAL             :: isInterpolateParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0)THEN
  isInterpolateParticle = .TRUE.
ELSE
  isInterpolateParticle = .FALSE.
END IF ! ABS(Species(PartSpecies(iPart))%ChargeIC).GT.0.0
END FUNCTION isInterpolateParticle


END MODULE MOD_part_tools
