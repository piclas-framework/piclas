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
! Contains tools for particles
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE UpdateNextFreePosition
  MODULE PROCEDURE UpdateNextFreePosition
END INTERFACE

INTERFACE VELOFROMDISTRIBUTION
  MODULE PROCEDURE VELOFROMDISTRIBUTION
END INTERFACE

INTERFACE CreateParticle
  MODULE PROCEDURE CreateParticle
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, VELOFROMDISTRIBUTION, GetParticleWeight, CreateParticle
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition()
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars        ,ONLY: PDM,PEM, PartSpecies, doParticleMerge, vMPF_SpecNumElem, PartPressureCell
USE MOD_Particle_Vars        ,ONLY: KeepWallParticles, PartState, VarTimeStep, usevMPF
USE MOD_DSMC_Vars            ,ONLY: useDSMC, CollInf
USE MOD_Particle_VarTimeStep ,ONLY: CalcVarTimeStep
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
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF(PDM%maxParticleNumber.EQ.0) RETURN
counter1 = 1
IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell.OR.usevMPF) THEN
  PEM%pNumber(:) = 0
  IF (KeepWallParticles) PEM%wNumber(:) = 0
END IF
n = PDM%ParticleVecLength !PDM%maxParticleNumber
PDM%ParticleVecLength = 0
PDM%insideParticleNumber = 0
IF (doParticleMerge) vMPF_SpecNumElem = 0
IF (useDSMC.OR.doParticleMerge.OR.PartPressureCell.OR.usevMPF) THEN
  DO i=1,n
    IF (.NOT.PDM%ParticleInside(i)) THEN
      IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
      PDM%nextFreePosition(counter1) = i
      counter1 = counter1 + 1
    ELSE
      IF (PEM%pNumber(PEM%Element(i)).EQ.0) THEN
        PEM%pStart(PEM%Element(i)) = i                    ! Start of Linked List for Particles in Elem
      ELSE
        PEM%pNext(PEM%pEnd(PEM%Element(i))) = i ! Next Particle of same Elem (Linked List)
      END IF
      PEM%pEnd(PEM%Element(i)) = i
      PEM%pNumber(PEM%Element(i)) = &                      ! Number of Particles in Element
          PEM%pNumber(PEM%Element(i)) + 1
      IF (VarTimeStep%UseVariableTimeStep) THEN
        VarTimeStep%ParticleTimeStep(i) = CalcVarTimeStep(PartState(i,1),PartState(i,2),PEM%Element(i))
      END IF
      IF (KeepWallParticles) THEN
        IF (PDM%ParticleAtWall(i)) THEN
          PEM%wNumber(PEM%Element(i)) = PEM%wNumber(PEM%Element(i)) + 1
        END IF
      END IF
      PDM%ParticleVecLength = i
      IF(doParticleMerge) vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) = vMPF_SpecNumElem(PEM%Element(i),PartSpecies(i)) + 1
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


FUNCTION VELOFROMDISTRIBUTION(distribution,specID,Tempergy)
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
USE MOD_Particle_Boundary_Tools ,ONLY: LIQUIDEVAP,LIQUIDREFL,ALPHALIQUID,BETALIQUID
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
REAL            :: veloFromDistribution(1:3)
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
  !veloFromDistribution(1) = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
  !veloFromDistribution(2) = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
  !CALL RANDOM_NUMBER(RandVal)
  !veloFromDistribution(3) = SQRT(-2*LOG(RandVal(1)))

  ! Get random vector
  veloFromDistribution = DiceUnitVector()
  ! Mirror z-component of velocity (particles are emitted from surface!)
  veloFromDistribution(3) = ABS(veloFromDistribution(3))
  ! Set magnitude
  veloFromDistribution = Tempergy*veloFromDistribution

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
  veloFromDistribution(3) = sigma*Velo1
  ! build tangential velocities from gauss (normal) distribution
  Velosq = 2
  DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    CALL RANDOM_NUMBER(RandVal)
    Velo1 = 2.*RandVal(1) - 1.
    Velo2 = 2.*RandVal(2) - 1.
    Velosq = Velo1**2 + Velo2**2
  END DO
  veloFromDistribution(1) = Velo1*SQRT(-2.*LOG(Velosq)/Velosq)*sigma
  veloFromDistribution(2) = Velo2*SQRT(-2.*LOG(Velosq)/Velosq)*sigma
CASE DEFAULT
  WRITE (UNIT_stdOut,'(A)') "distribution =", distribution
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')
END SELECT

END FUNCTION VELOFROMDISTRIBUTION


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

SUBROUTINE CreateParticle(Species,Pos,ElemID,Velocity,RotEnergy,VibEnergy,ElecEnergy,NewPartID)
!===================================================================================================================================
!> creates a single particle at correct array position and assign properties
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: PDM, PEM, PartState, LastPartPos, PartSpecies
USE MOD_DSMC_Vars     ,ONLY: useDSMC, CollisMode, DSMC, PartStateIntEn     ! , RadialWeighting
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN)           :: Species
REAL, INTENT(IN)              :: Pos(1:3)
INTEGER, INTENT(IN)           :: ElemID
REAL, INTENT(IN)              :: Velocity(1:3)
REAL, INTENT(IN)              :: RotEnergy
REAL, INTENT(IN)              :: VibEnergy
REAL, INTENT(IN)              :: ElecEnergy
INTEGER, INTENT(OUT),OPTIONAL :: NewPartID
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER :: newParticleID
!===================================================================================================================================


!IPWRITE(UNIT_stdOut,*) 'NEW PARTICLE!'

newParticleID = PDM%nextFreePosition(PDM%CurrentNextFreePosition+1)
IF (newParticleID .EQ. 0) CALL abort(&
__STAMP__&
,'ERROR in CreateParticle: newParticleID.EQ.0 - maximum nbr of particles reached?')
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1

PartSpecies(newParticleID) = Species
LastPartPos(newParticleID,1:3)=Pos(1:3)
PartState(newParticleID,1:3) = Pos(1:3)
PartState(newParticleID,4:6) = Velocity(1:3)

IF (useDSMC.AND.(CollisMode.GT.1)) THEN
  PartStateIntEn(newParticleID, 1) = VibEnergy
  PartStateIntEn(newParticleID, 2) = RotEnergy
  IF (DSMC%ElectronicModel) THEN
    PartStateIntEn(newParticleID, 3) = ElecEnergy
  ENDIF
END IF

PDM%ParticleInside(newParticleID) = .TRUE.
PDM%dtFracPush(newParticleID)     = .FALSE.
PDM%IsNewPart(newParticleID)      = .FALSE.   ! ??????? correct ????
PEM%Element(newParticleID)        = ElemID
PEM%lastElement(newParticleID)    = ElemID

! ?????? necessary?
! IF (VarTimeStep%UseVariableTimeStep) THEN
!   VarTimeStep%ParticleTimeStep(newParticleID) &
!     = CalcVarTimeStep(PartState(newParticleID,1),PartState(newParticleID,2),PEM%Element(newParticleID))
! END IF
! IF (RadialWeighting%DoRadialWeighting) THEN
!   PartMPF(newParticleID) = CalcRadWeightMPF(PartState(newParticleID,2), 1,newParticleID)
! END IF
IF (PRESENT(NewPartID)) NewPartID=newParticleID

END SUBROUTINE CreateParticle

END MODULE MOD_part_tools
