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

MODULE MOD_part_tools
!===================================================================================================================================
! Contains tools for particle related operations. This routine uses MOD_Particle_Boundary_Tools, but not vice versa!
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

INTERFACE DiceUnitVector
  MODULE PROCEDURE DiceUnitVector
END INTERFACE

INTERFACE isChargedParticle
  MODULE PROCEDURE isChargedParticle
END INTERFACE

INTERFACE isPushParticle
  MODULE PROCEDURE isPushParticle
END INTERFACE

INTERFACE InRotRefFrameCheck
  MODULE PROCEDURE InRotRefFrameCheck
END INTERFACE

INTERFACE isDepositParticle
  MODULE PROCEDURE isDepositParticle
END INTERFACE

INTERFACE isInterpolateParticle
  MODULE PROCEDURE isInterpolateParticle
END INTERFACE

INTERFACE StoreLostParticleProperties
  MODULE PROCEDURE StoreLostParticleProperties
END INTERFACE

INTERFACE InterpolateEmissionDistribution2D
  MODULE PROCEDURE InterpolateEmissionDistribution2D
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: UpdateNextFreePosition, DiceUnitVector, VeloFromDistribution, GetParticleWeight, CalcRadWeightMPF, isChargedParticle
PUBLIC :: isPushParticle, isDepositParticle, isInterpolateParticle, StoreLostParticleProperties, BuildTransGaussNums
PUBLIC :: CalcXiElec,ParticleOnProc,  CalcERot_particle, CalcEVib_particle, CalcEElec_particle, CalcVelocity_maxwell_particle
PUBLIC :: InitializeParticleMaxwell
PUBLIC :: InterpolateEmissionDistribution2D
PUBLIC :: MergeCells,InRotRefFrameCheck
!===================================================================================================================================

CONTAINS

SUBROUTINE UpdateNextFreePosition(WithOutMPIParts)
!===================================================================================================================================
! Updates next free position
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars            ,ONLY: useDSMC,CollInf
USE MOD_Particle_Vars        ,ONLY: PDM,PEM,PartSpecies
USE MOD_Particle_Vars        ,ONLY: PartState,PartTimeStep,usevMPF,UseVarTimeStep
USE MOD_Particle_TimeStep    ,ONLY: GetParticleTimeStep
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers   ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Particle_MPI_Vars    ,ONLY: PartTargetProc
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, OPTIONAL         :: WithOutMPIParts
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: counter,i
INTEGER            :: ElemID
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

IF(PDM%maxParticleNumber.EQ.0) RETURN

IF (useDSMC.OR.usevMPF) THEN
  PEM%pNumber(:) = 0
END IF

PDM%ParticleVecLengthOld = PDM%ParticleVecLength
counter                  = 0
PDM%ParticleVecLength    = 0

! Check size of PDM%ParticleInside array vs. PDM%ParticleVecLength. During particle splitting, the max particle number might be
! exceeded, which may lead to an out-of-bounds here
IF(usevMPF)THEN
  IF(PDM%ParticleVecLengthOld.GT.SIZE(PDM%ParticleInside))THEN
    IPWRITE(UNIT_StdOut,*) "PDM%ParticleVecLength    :", PDM%ParticleVecLengthOld
    IPWRITE(UNIT_StdOut,*) "SIZE(PDM%ParticleInside) :", SIZE(PDM%ParticleInside)
    CALL abort(__STAMP__,'PDM%ParticleVecLength exceeds allocated arrays. Possible vMPF overflow during particle split or '//&
                         'restart file with too many particles. Increase Part-maxParticleNumber or use more processors')
  END IF ! PDM%ParticleVecLengthOld.GT.SIZE(PDM%ParticleInside)
END IF ! usevMPF

IF (useDSMC.OR.usevMPF) THEN
  DO i = 1,PDM%ParticleVecLengthOld
    IF (.NOT.PDM%ParticleInside(i)) THEN
      IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
      counter = counter + 1
      PDM%nextFreePosition(counter) = i
#if USE_MPI
    ELSE IF (PRESENT(WithOutMPIParts)) THEN
      ! Particle is to be sent to another proc
      IF (PartTargetProc(i).NE.-1) THEN
        IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
        counter = counter + 1
        PDM%nextFreePosition(counter) = i
      ! Particle will stay on the current proc
      ELSE
        ElemID = PEM%LocalElemID(i)
        ! Start of linked list for particles in elem
        IF (PEM%pNumber(ElemID).EQ.0) THEN
          PEM%pStart(ElemID)          = i
        ! Next particle of same elem (linked list)
        ELSE
          PEM%pNext(PEM%pEnd(ElemID)) = i
        END IF
        PEM%pEnd(   ElemID)   = i
        PEM%pNumber(ElemID)   = PEM%pNumber(ElemID) + 1
        PDM%ParticleVecLength = i
        IF(UseVarTimeStep) PartTimeStep(i) = GetParticleTimeStep(PartState(1,i),PartState(2,i),ElemID)
      END IF
#endif
    ELSE
      ! Sanity check corrupted particle list (some or all entries of a particle become zero, including the species ID)
      IF(PartSpecies(i).LE.0) CALL abort(__STAMP__,'Species ID is zero for ipart=',IntInfoOpt=i)
      ElemID = PEM%LocalElemID(i)
      ! Start of linked list for particles in elem
      IF (PEM%pNumber(ElemID).EQ.0) THEN
        PEM%pStart(ElemID)          = i
      ! Next particle of same elem (linked list)
      ELSE
        PEM%pNext(PEM%pEnd(ElemID)) = i
      END IF
      PEM%pEnd(   ElemID)   = i
      PEM%pNumber(ElemID)   = PEM%pNumber(ElemID) + 1
      PDM%ParticleVecLength = i
      IF(UseVarTimeStep) PartTimeStep(i) = GetParticleTimeStep(PartState(1,i),PartState(2,i),ElemID)
    END IF
  END DO
! no DSMC
ELSE
  DO i = 1,PDM%ParticleVecLengthOld
    IF (.NOT.PDM%ParticleInside(i)) THEN
      counter = counter + 1
      PDM%nextFreePosition(counter) = i
#if USE_MPI
    ELSE IF (PRESENT(WithOutMPIParts)) THEN
      ! Particle is to be sent to another proc
      IF (PartTargetProc(i).NE.-1) THEN
        IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
        counter = counter + 1
        PDM%nextFreePosition(counter) = i
      ! Particle will stay on the current proc
      ELSE
        PDM%ParticleVecLength = i
      END IF
#endif
    ELSE
      PDM%ParticleVecLength = i
    END IF
  END DO
ENDIF
PDM%CurrentNextFreePosition = 0

! Positions after ParticleVecLength in freePosition
DO i = PDM%ParticleVecLengthOld+1,PDM%maxParticleNumber
  IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(i) = 0
  counter = counter + 1
  PDM%nextFreePosition(counter) = i
END DO

! Set nextFreePosition for occupied slots to zero
PDM%nextFreePosition(counter+1:PDM%maxParticleNumber) = 0
! If maxParticleNumber are inside, counter is greater than maxParticleNumber
IF (counter+1.GT.PDM%MaxParticleNumber) PDM%nextFreePosition(PDM%MaxParticleNumber) = 0

#if USE_LOADBALANCE
CALL LBPauseTime(LB_UNFP,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE UpdateNextFreePosition


SUBROUTINE StoreLostParticleProperties(iPart,ElemID,UsePartState_opt,PartMissingType_opt)
!----------------------------------------------------------------------------------------------------------------------------------!
! Store information of a lost particle (during restart and during the simulation)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Globals                ,ONLY: abort,myrank
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,PartSpecies,Species,PartState,LastPartPos
USE MOD_Particle_Tracking_Vars ,ONLY: PartStateLost,PartLostDataSize,PartStateLostVecLength
USE MOD_TimeDisc_Vars          ,ONLY: time
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: iPart
INTEGER,INTENT(IN)          :: ElemID ! Global element index
LOGICAL,INTENT(IN),OPTIONAL :: UsePartState_opt
INTEGER,INTENT(IN),OPTIONAL :: PartMissingType_opt ! 0: lost, 1: missing & found once, >1: missing & multiply found
INTEGER                     :: dims(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF
! Temporary arrays
REAL, ALLOCATABLE    :: PartStateLost_tmp(:,:)   ! (1:11,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,MPF,time,ElemID,iPart
!                                                !                 2nd index: 1 to number of lost particles
INTEGER              :: ALLOCSTAT
LOGICAL              :: UsePartState_loc
!===================================================================================================================================
UsePartState_loc = .FALSE.
IF (PRESENT(UsePartState_opt)) UsePartState_loc = UsePartState_opt

! Set macro particle factor
IF (usevMPF) THEN
  MPF = PartMPF(iPart)
ELSE
  MPF = Species(PartSpecies(iPart))%MacroParticleFactor
END IF

dims = SHAPE(PartStateLost)

ASSOCIATE( iMax => PartStateLostVecLength )
  ! Increase maximum number of boundary-impact particles
  iMax = iMax + 1

  ! Check if array maximum is reached.
  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
  IF(iMax.GT.dims(2))THEN

    ! --- PartStateLost ---
    ALLOCATE(PartStateLost_tmp(1:PartLostDataSize,1:dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateLost_tmp array!')
    ! Save old data
    PartStateLost_tmp(1:PartLostDataSize,1:dims(2)) = PartStateLost(1:PartLostDataSize,1:dims(2))

    ! Re-allocate PartStateLost to twice the size
    DEALLOCATE(PartStateLost)
    ALLOCATE(PartStateLost(1:PartLostDataSize,1:2*dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateLost array!')
    PartStateLost(1:PartLostDataSize,        1:  dims(2)) = PartStateLost_tmp(1:PartLostDataSize,1:dims(2))
    PartStateLost(1:PartLostDataSize,dims(2)+1:2*dims(2)) = 0.

  END IF

  ! 1-3: Particle position (last valid position)
  IF(UsePartState_loc)THEN
    PartStateLost(1:3,iMax) = PartState(1:3,iPart)
  ELSE
    PartStateLost(1:3,iMax) = LastPartPos(1:3,iPart)
  END IF ! UsePartState_loc
  ! 4-6: Particle velocity
  PartStateLost(4:6  ,iMax) = PartState(4:6,iPart)
  ! 7: SpeciesID
  PartStateLost(7    ,iMax) = REAL(PartSpecies(iPart))
  ! 8: Macro particle factor
  PartStateLost(8    ,iMax) = MPF
  ! 9: time of loss
  PartStateLost(9    ,iMax) = time
  ! 10: Global element ID
  PartStateLost(10   ,iMax) = REAL(ElemID)
  ! 11: Particle ID
  PartStateLost(11   ,iMax) = REAL(iPart)
  ! 12-14: Particle position (position of loss)
  PartStateLost(12:14,iMax) = PartState(1:3,iPart)
  ! 15: myrank
  PartStateLost(15,iMax) = myrank
  ! 16: missing type, i.e., 0: lost, 1: missing & found once, >1: missing & multiply found
  IF(PRESENT(PartMissingType_opt))THEN ! when particles go missing during restart (maybe they are found on other procs or lost)
    PartStateLost(16,iMax) = PartMissingType_opt
  ELSE ! simply lost during the simulation
    PartStateLost(16,iMax) = 0
  END IF ! PRESENT(PartMissingType_opt)
END ASSOCIATE

END SUBROUTINE StoreLostParticleProperties


FUNCTION DiceUnitVector()
!===================================================================================================================================
!> Calculates random unit vector
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
REAL                     :: rRan, cos_scatAngle, sin_scatAngle, rotAngle
!===================================================================================================================================
CALL RANDOM_NUMBER(rRan)

cos_scatAngle     = 2.*rRan-1.
sin_scatAngle     = SQRT(1. - cos_scatAngle ** 2.)
DiceUnitVector(1) = cos_scatAngle

CALL RANDOM_NUMBER(rRan)
rotAngle          = 2. * Pi * rRan

DiceUnitVector(2) = sin_scatAngle * COS(rotAngle)
DiceUnitVector(3) = sin_scatAngle * SIN(rotAngle)

END FUNCTION DiceUnitVector


!===================================================================================================================================
!> Calculation of velocity vector (Vx,Vy,Vz) sampled from given distribution function
!===================================================================================================================================
FUNCTION VeloFromDistribution(distribution,Tempergy,iNewPart,ProductSpecNbr)
! MODULES
USE MOD_Globals           ,ONLY: Abort,UNIT_stdOut,VECNORM
USE MOD_Globals_Vars      ,ONLY: eV2Joule,ElectronMass,c
USE MOD_SurfaceModel_Vars ,ONLY: BackupVeloABS
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: distribution   !< Specifying keyword for velocity distribution
REAL,INTENT(IN)             :: Tempergy       !< Input temperature in [K] or energy in [J] or [eV] or velocity in [m/s]
INTEGER,INTENT(IN)          :: iNewPart       !< The i-th particle that is inserted (only required for some distributions)
INTEGER,INTENT(IN)          :: ProductSpecNbr !< Total number of particles that are inserted (only required for some distributions)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: VeloFromDistribution(1:3) !< Velocity vector created from specific velocity distribution function
REAL            :: VeloABS                   !< Absolute velocity of the velocity vector
REAL            :: RandVal                   !< Pseudo random number
LOGICAL         :: ARM                       !< Acceptance rejection method
REAL            :: PDF                       !< Probability density function
REAL            :: eps,eps2                  !< kinetic electron energy [eV]
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!-- set velocities
SELECT CASE(TRIM(distribution))

CASE('deltadistribution')

  ! Get random vector
  VeloFromDistribution = DiceUnitVector()
  ! Mirror z-component of velocity (particles are emitted from surface!)
  VeloFromDistribution(3) = ABS(VeloFromDistribution(3))
  ! Set magnitude
  VeloFromDistribution = Tempergy*VeloFromDistribution ! Tempergy here is [m/s]

CASE('uniform-energy')

  ! Get random vector
  VeloFromDistribution = DiceUnitVector()
  ! Mirror z-component of velocity (particles are emitted from surface!)
  VeloFromDistribution(3) = ABS(VeloFromDistribution(3))
  ! Set uniform energy distribution. Note that Tempergy here is [eV], which is converted from eV to J
  CALL RANDOM_NUMBER(RandVal) ! random value between 0 and 1.
  VeloABS = SQRT(2.0 * RandVal * Tempergy * eV2Joule / ElectronMass)
  ! Set magnitude
  VeloFromDistribution = VeloABS*VeloFromDistribution

CASE('Morozov2004') ! Secondary electron emission (SEE) due to electron bombardment on dielectric surfaces

  IF(ProductSpecNbr.EQ.1)THEN ! 1 SEE

    ! ARM for energy distribution
    ARM = .TRUE.
    DO WHILE(ARM)
      CALL RANDOM_NUMBER(RandVal) ! random x-coordinate
      PDF = 4.0*RandVal*(1.0-RandVal)
      eps = RandVal ! RandVal is eps/eps_p (relative energy as compared with the incident electron energy)
      CALL RANDOM_NUMBER(RandVal) ! random y-coordinate
      IF (RandVal.LT.PDF) ARM = .FALSE.
    END DO
    VeloABS = SQRT(2.0 * eps * Tempergy * eV2Joule / ElectronMass) ! Tempergy here is [eV], which is converted from eV to J

  ELSE ! 2 SEE

    ! ARM for energy distribution
    IF(iNewPart.EQ.1)THEN ! 1st call
      ARM = .TRUE.
      DO WHILE(ARM)
        ! Pick 1st electron energy
        CALL RANDOM_NUMBER(RandVal)
        PDF = 4.0*RandVal*(1.0-RandVal)
        eps = RandVal ! RandVal is eps/eps_p (relative energy as compared with the incident electron energy)
        CALL RANDOM_NUMBER(RandVal)
        IF (RandVal.LT.PDF)THEN
          ! Pick 2nd electron energy
          CALL RANDOM_NUMBER(RandVal)
          PDF = 4.0*RandVal*(1.0-RandVal)
          eps2 = RandVal ! RandVal is eps/eps_p (relative energy as compared with the incident electron energy)
          CALL RANDOM_NUMBER(RandVal)
          IF(RandVal.LT.PDF)THEN
            ARM = .FALSE. ! success, skip this loop and skip the outer loop
            ! eV to J: store 2nd electron velocity for next function call
            BackupVeloABS = SQRT(2.0 * eps2 * Tempergy * eV2Joule / ElectronMass) ! Tempergy here is [eV] (converted from eV to J)
            IF(eps+eps2.GT.1.0) ARM = .TRUE. ! start again for both energies
          END IF
        END IF
      END DO
      VeloABS = SQRT(2.0 * eps * Tempergy * eV2Joule / ElectronMass) ! Tempergy here is [eV], which is converted from eV to J

    ELSE ! 2nd call of this function, velocity already known from last call and stored in "BackupVeloABS"

      VeloABS = BackupVeloABS

    END IF ! iNewPart.EQ.1

  END IF ! ProductSpecNbr.EQ.1

  ! Get random vector
  VeloFromDistribution = DiceUnitVector()
  ! Mirror z-component of velocity (particles are emitted from surface!)
  VeloFromDistribution(3) = ABS(VeloFromDistribution(3))
  ! Set magnitude
  VeloFromDistribution = VeloABS*VeloFromDistribution ! VeloABS is [m/s]

CASE DEFAULT

  CALL abort(__STAMP__,'Unknown velocity dsitribution: ['//TRIM(distribution)//']')

END SELECT

! Sanity check: is the newly created particle faster than c
IF(VECNORM(VeloFromDistribution).GT.c) CALL abort(__STAMP__,'VeloFromDistribution: Particle is faster than the speed of light: ',&
    RealInfoOpt=VeloABS)

END FUNCTION VeloFromDistribution


PPURE REAL FUNCTION GetParticleWeight(iPart)
!===================================================================================================================================
!> Determines the appropriate particle weighting for the axisymmetric case with radial weighting and the variable time step. For
!> radial weighting, the radial factor is multiplied by the regular weighting factor. If only a variable time step is used, at the
!> moment, the regular weighting factor is not included.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Particle_Vars           ,ONLY: usevMPF, UseVarTimeStep, PartTimeStep, PartMPF
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  IF (UseVarTimeStep) THEN
    GetParticleWeight = PartMPF(iPart) * PartTimeStep(iPart)
  ELSE
    GetParticleWeight = PartMPF(iPart)
  END IF
ELSE IF (UseVarTimeStep) THEN
  GetParticleWeight = PartTimeStep(iPart)
ELSE
  GetParticleWeight = 1.
END IF

END FUNCTION GetParticleWeight


PPURE REAL FUNCTION CalcRadWeightMPF(yPos, iSpec, iPart)
!===================================================================================================================================
!> Determines the weighting factor when using an additional radial weighting for axisymmetric simulations. Linear increase from the
!> rotational axis (y=0) to the outer domain boundary (y=ymax).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
USE MOD_Particle_Vars           ,ONLY: Species, PEM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemMidPoint_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: yPos
INTEGER, INTENT(IN)             :: iSpec
INTEGER, OPTIONAL,INTENT(IN)    :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                 :: yPosIn
!===================================================================================================================================

IF(RadialWeighting%CellLocalWeighting.AND.PRESENT(iPart)) THEN
  yPosIn = ElemMidPoint_Shared(2,PEM%CNElemID(iPart))
ELSE
  yPosIn = yPos
END IF

CalcRadWeightMPF = (1. + yPosIn/GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))*Species(iSpec)%MacroParticleFactor

RETURN

END FUNCTION CalcRadWeightMPF


PPURE FUNCTION isChargedParticle(iPart)
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


PPURE FUNCTION isPushParticle(iPart)
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

PPURE FUNCTION InRotRefFrameCheck(iPart)
!----------------------------------------------------------------------------------------------------------------------------------!
! Check if particle is in a rotating frame of reference region.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars ,ONLY: PartState,RotRefFramRegion,RotRefFrameAxis,nRefFrameRegions
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: iPart
LOGICAL             :: InRotRefFrameCheck
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iRegion
!===================================================================================================================================

IF(nRefFrameRegions.GT.0) THEN
  InRotRefFrameCheck = .FALSE.
  DO iRegion = 1, nRefFrameRegions
    IF((PartState(RotRefFrameAxis,iPart).GT.RotRefFramRegion(1,iRegion)).AND. &
       (PartState(RotRefFrameAxis,iPart).LT.RotRefFramRegion(2,iRegion))) THEN
      InRotRefFrameCheck = .TRUE.
      EXIT
    END IF
  END DO
ELSE
  InRotRefFrameCheck = .TRUE.
END IF

END FUNCTION InRotRefFrameCheck

PPURE FUNCTION isDepositParticle(iPart)
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


PPURE FUNCTION isInterpolateParticle(iPart)
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


SUBROUTINE BuildTransGaussNums(nPart, iRanPart)
!===================================================================================================================================
!> Builds random Gauss numbers with a zero mean and a variance of one
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: sumiRan(3), varianceiRan(3)
INTEGER                        :: iLoop,I
!===================================================================================================================================
sumiRan(1:3) = 0.0
varianceiRan(1:3) = 0.0
DO iLoop = 1, nPart
  iRanPart(1,iLoop) = rnor()
  iRanPart(2,iLoop) = rnor()
  iRanPart(3,iLoop) = rnor()
  sumiRan(1:3) = sumiRan(1:3) + iRanPart(1:3,iLoop)
END DO
sumiRan(1:3) = sumiRan(1:3)/nPart
DO iLoop = 1, nPart
  iRanPart(1:3,iLoop) = iRanPart(1:3,iLoop)-sumiRan(1:3)
  varianceiRan(1:3) = varianceiRan(1:3) + iRanPart(1:3,iLoop)*iRanPart(1:3,iLoop)
END DO
varianceiRan(1:3) = SQRT(varianceiRan(1:3)/nPart)

DO iLoop = 1, nPart
  DO I = 1, 3
    IF(varianceiRan(I).GT.0)THEN ! Catch division by zero
      iRanPart(I,iLoop) = iRanPart(I,iLoop)/varianceiRan(I)
    ELSE
      iRanPart(I,iLoop) = 0.
    END IF ! varianceiRan(I).GT.0
  END DO ! I = 1, 3
END DO

END SUBROUTINE BuildTransGaussNums


PPURE REAL FUNCTION CalcXiElec(Telec, iSpec)
!===================================================================================================================================
!> Calculation of the electronic degree of freedom for a given temperature and species
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: Telec  !
INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: iQua
REAL                            :: TempRatio, SumOne, SumTwo
!===================================================================================================================================

IF (Telec.LE.0.0) THEN
  CalcXiElec = 0.0
  RETURN
END IF

SumOne = 0.0
SumTwo = 0.0

DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant-1
  TempRatio = SpecDSMC(iSpec)%ElectronicState(2,iQua)/Telec
  IF(CHECKEXP(TempRatio)) THEN
    SumOne = SumOne + SpecDSMC(iSpec)%ElectronicState(1,iQua)*SpecDSMC(iSpec)%ElectronicState(2,iQua)*EXP(-TempRatio)
    SumTwo = SumTwo + SpecDSMC(iSpec)%ElectronicState(1,iQua)*EXP(-TempRatio)
  END IF
END DO

IF((SumOne.GT.0.0).AND.(SumTwo.GT.0.0)) THEN
  CalcXiElec = 2. * SumOne / (SumTwo * Telec)
ELSE
  CalcXiElec = 0.0
END IF

RETURN

END FUNCTION CalcXiElec


!===================================================================================================================================
!> Check whether particle host element ID is on the current proc
!===================================================================================================================================
PPURE LOGICAL FUNCTION ParticleOnProc(PartID) RESULT(L)
! MODULES
USE MOD_Preproc
USE MOD_Particle_Vars ,ONLY: PEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: PartID ! Particle index
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
L = (PEM%LocalElemID(PartID).GE.1).AND.(PEM%LocalElemID(PartID).LE.PP_nElems)
END FUNCTION ParticleOnProc


FUNCTION CalcVelocity_maxwell_particle(iSpec,Temp)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species
USE Ziggurat,                   ONLY : rnor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec
REAL, INTENT(IN)                :: Temp(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: CalcVelocity_maxwell_particle(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CalcVelocity_maxwell_particle(1:3) = 0.0

IF(Temp(1).GT.0.0) CalcVelocity_maxwell_particle(1) = rnor()*SQRT(BoltzmannConst*Temp(1)/Species(iSpec)%MassIC)
IF(Temp(2).GT.0.0) CalcVelocity_maxwell_particle(2) = rnor()*SQRT(BoltzmannConst*Temp(2)/Species(iSpec)%MassIC)
IF(Temp(3).GT.0.0) CalcVelocity_maxwell_particle(3) = rnor()*SQRT(BoltzmannConst*Temp(3)/Species(iSpec)%MassIC)

END FUNCTION CalcVelocity_maxwell_particle


REAL FUNCTION CalcEVib_particle(iSpec,TempVib,iPart)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC, PolyatomMolDSMC, VibQuantsPar, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: TempVib
INTEGER, INTENT(IN),OPTIONAL  :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: iRan
INTEGER                   :: iQuant, iDOF, iPolyatMole
LOGICAL                   :: SetVibQuant
!===================================================================================================================================

IF(PRESENT(iPart)) THEN
  SetVibQuant = .TRUE.
ELSE
  SetVibQuant = .FALSE.
END IF

IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
  ! set vibrational energy
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  IF(SetVibQuant) THEN
    IF(ALLOCATED(VibQuantsPar(iPart)%Quants)) DEALLOCATE(VibQuantsPar(iPart)%Quants)
    ALLOCATE(VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
  END IF
  CalcEVib_particle = 0.0
  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TempVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
    DO WHILE (iQuant.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*TempVib/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
    END DO
    CalcEVib_particle = CalcEVib_particle &
                                + (iQuant + DSMC%GammaQuant)*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
    IF(SetVibQuant) VibQuantsPar(iPart)%Quants(iDOF)=iQuant
  END DO
ELSE
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT(-LOG(iRan)*TempVib/SpecDSMC(iSpec)%CharaTVib)
  DO WHILE (iQuant.GE.SpecDSMC(iSpec)%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TempVib/SpecDSMC(iSpec)%CharaTVib)
  END DO
  CalcEVib_particle = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
END IF

RETURN

END FUNCTION CalcEVib_particle


REAL FUNCTION CalcERot_particle(iSpec,TempRot)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iSpec
REAL, INTENT(IN)          :: TempRot
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                      :: PartStateTempVar, NormProb, iRan2
!===================================================================================================================================

CalcERot_particle = 0.

IF (SpecDSMC(iSpec)%Xi_Rot.EQ.2) THEN
  CALL RANDOM_NUMBER(iRan2)
  CalcERot_particle = -BoltzmannConst*TempRot*LOG(iRan2)
ELSE IF (SpecDSMC(iSpec)%Xi_Rot.EQ.3) THEN
  CALL RANDOM_NUMBER(iRan2)
  PartStateTempVar = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
  NormProb = SQRT(PartStateTempVar)*EXP(-PartStateTempVar)/(SQRT(0.5)*EXP(-0.5))
  CALL RANDOM_NUMBER(iRan2)
  DO WHILE (iRan2.GE.NormProb)
    CALL RANDOM_NUMBER(iRan2)
    PartStateTempVar = iRan2*10 !the distribution function has only non-negligible  values betwenn 0 and 10
    NormProb = SQRT(PartStateTempVar)*EXP(-PartStateTempVar)/(SQRT(0.5)*EXP(-0.5))
    CALL RANDOM_NUMBER(iRan2)
  END DO
  CalcERot_particle = PartStateTempVar*BoltzmannConst*TempRot
END IF

RETURN

END FUNCTION CalcERot_particle


REAL FUNCTION CalcEElec_particle(iSpec,TempElec,iPart)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars      ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars         ,ONLY: SpecDSMC, DSMC, ElectronicDistriPart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: TempElec
INTEGER, INTENT(IN),OPTIONAL  :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iQua
REAL                      :: iRan, ElectronicPartition, ElectronicPartitionTemp, tmpExp
!===================================================================================================================================
ElectronicPartition  = 0.
CalcEElec_particle = 0.

IF(.NOT.PRESENT(iPart).AND.DSMC%ElectronicModel.EQ.2) THEN
  CALL abort(__STAMP__,'ERROR: Calculation of electronic energy using ElectronicModel = 2 requires the input of particle index!')
END IF

SELECT CASE(DSMC%ElectronicModel)
CASE(1,4)
  ElectronicPartitionTemp = 0.
  IF(TempElec.GT.0.0) THEN
    ! calculate sum over all energy levels == partition function for temperature Telec
    DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TempElec)
      IF ( ElectronicPartitionTemp .GT. ElectronicPartition ) THEN
        ElectronicPartition = ElectronicPartitionTemp
      END IF
    END DO
    ElectronicPartitionTemp = 0.
    ! select level
    CALL RANDOM_NUMBER(iRan)
    DO WHILE ( iRan .GE. ElectronicPartitionTemp / ElectronicPartition )
      CALL RANDOM_NUMBER(iRan)
      iQua = int( ( SpecDSMC(iSpec)%MaxElecQuant ) * iRan)
      ElectronicPartitionTemp = SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iQua)/TempElec)
      CALL RANDOM_NUMBER(iRan)
    END DO
  ELSE
    iQua = 0
  END IF
  CalcEElec_particle = BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
CASE(2)
  IF(ALLOCATED(ElectronicDistriPart(iPart)%DistriFunc)) DEALLOCATE(ElectronicDistriPart(iPart)%DistriFunc)
  ALLOCATE(ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(iSpec)%MaxElecQuant))
  CalcEElec_particle = 0.0
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TempElec
    IF (CHECKEXP(tmpExp)) &
      ElectronicPartition = ElectronicPartition + SpecDSMC(iSpec)%ElectronicState(1,iQua) * EXP(-tmpExp)
  END DO
  DO iQua = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
    tmpExp = SpecDSMC(iSpec)%ElectronicState(2,iQua) / TempElec
    IF (CHECKEXP(tmpExp)) THEN
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = SpecDSMC(iSpec)%ElectronicState(1,iQua)*EXP(-tmpExp)/ElectronicPartition
    ELSE
      ElectronicDistriPart(iPart)%DistriFunc(iQua+1) = 0.0
    END IF
    CalcEElec_particle = CalcEElec_particle + &
        ElectronicDistriPart(iPart)%DistriFunc(iQua+1) * BoltzmannConst * SpecDSMC(iSpec)%ElectronicState(2,iQua)
  END DO
CASE(3)
  ! Initialize in ground state
  CalcEElec_particle = 0.
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR: Unknown electronic relaxation model: ',IntInfoOpt=DSMC%ElectronicModel)
END SELECT

RETURN

END FUNCTION CalcEElec_particle


SUBROUTINE MergeCells()
!===================================================================================================================================
!> Routine for virtual merging of neighbouring cells. 
!> Currently, the merging is only done via the number of particles within the cells.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,        ONLY: VirtMergedCells, PEM, MinPartNumCellMerge, VirtualCellMergeSpread, MaxNumOfMergedCells
USE MOD_Mesh_Vars,            ONLY: nElems,offsetElem
USE MOD_Particle_Mesh_Vars,   ONLY: ElemToElemMapping,ElemToElemInfo, ElemVolume_Shared
USE MOD_Mesh_Tools,           ONLY: GetCNElemID, GetGlobalElemID
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nPart, CNElemID, GlobalElemID, iNbElem, GlobNbElem, LocNBElem, nPartMerged, MasterCellID
INTEGER                        :: CNNbElem, iElem, iOldElem, currentCellCount
INTEGER, ALLOCATABLE           :: tempCellID(:)
LOGICAL                        :: AllowBackMerge
!===================================================================================================================================
!Nullify every value
DO iElem = 1, nElems
  VirtMergedCells(iElem)%isMerged = .FALSE.
  VirtMergedCells(iElem)%MasterCell = 0  
  VirtMergedCells(iElem)%MergedVolume = 0.0
  IF (VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    DEALLOCATE(VirtMergedCells(iElem)%MergedCellID)
    VirtMergedCells(iElem)%NumOfMergedCells=0
  END IF
END DO

!Loop over all cells and neighbouring cells to merge them
ElemLoop: DO iElem = 1, nElems
  IF(VirtMergedCells(iElem)%isMerged) CYCLE
  nPart = PEM%pNumber(iElem)
  IF (nPart.LE.MinPartNumCellMerge) THEN
    GlobalElemID = iElem + offSetElem
    CNElemID = GetCNElemID(GlobalElemID)
    nPartMerged = nPart
    AllowBackMerge = .TRUE.
    NBElemLoop: DO iNbElem = 1,ElemToElemMapping(2,CNElemID)
      GlobNbElem = GetGlobalElemID(ElemToElemInfo(ElemToElemMapping(1,CNElemID)+iNbElem))
      LocNBElem = GlobNbElem-offSetElem
      CNNbElem = GetCNElemID(GlobNbElem)
      IF ((LocNBElem.LT.1).OR.(LocNBElem.GT.nElems)) CYCLE
      IF(VirtMergedCells(LocNBElem)%isMerged.AND.AllowBackMerge) THEN
        IF(VirtualCellMergeSpread.GT.1) THEN
          IF (VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
            MasterCellID = VirtMergedCells(LocNBElem)%MasterCell-offSetElem
            IF(VirtMergedCells(MasterCellID)%NumOfMergedCells.GE.(MaxNumOfMergedCells-1)) CYCLE NBElemLoop   
            ALLOCATE(tempCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells))
            tempCellID = VirtMergedCells(MasterCellID)%MergedCellID
            DEALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID)
            VirtMergedCells(MasterCellID)%NumOfMergedCells = VirtMergedCells(MasterCellID)%NumOfMergedCells + 1 + VirtMergedCells(iElem)%NumOfMergedCells
            ALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells))
            VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1-VirtMergedCells(iElem)%NumOfMergedCells) = &
              tempCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1-VirtMergedCells(iElem)%NumOfMergedCells)
            oldMasterElemLoop: DO iOldElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
              currentCellCount = VirtMergedCells(MasterCellID)%NumOfMergedCells-VirtMergedCells(iElem)%NumOfMergedCells -1 + iOldElem
              VirtMergedCells(MasterCellID)%MergedCellID(currentCellCount) = VirtMergedCells(iElem)%MergedCellID(iOldElem)
              VirtMergedCells(VirtMergedCells(iElem)%MergedCellID(iOldElem))%MasterCell = VirtMergedCells(LocNBElem)%MasterCell
            END DO oldMasterElemLoop
            VirtMergedCells(MasterCellID)%MergedCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells) = iElem
            VirtMergedCells(iElem)%MasterCell = VirtMergedCells(LocNBElem)%MasterCell
            VirtMergedCells(iElem)%isMerged = .TRUE.
            VirtMergedCells(MasterCellID)%MergedVolume=VirtMergedCells(MasterCellID)%MergedVolume+VirtMergedCells(iElem)%MergedVolume
            VirtMergedCells(iElem)%MergedVolume = 0.0
            VirtMergedCells(iElem)%NumOfMergedCells = 0
            DEALLOCATE(VirtMergedCells(iElem)%MergedCellID)
            DEALLOCATE(tempCellID)
            CYCLE ElemLoop
          ELSE
            MasterCellID = VirtMergedCells(LocNBElem)%MasterCell-offSetElem
            IF(VirtMergedCells(MasterCellID)%NumOfMergedCells.GE.(MaxNumOfMergedCells-1)) CYCLE NBElemLoop
            ALLOCATE(tempCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells))
            tempCellID = VirtMergedCells(MasterCellID)%MergedCellID
            DEALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID)
            VirtMergedCells(MasterCellID)%NumOfMergedCells = VirtMergedCells(MasterCellID)%NumOfMergedCells + 1
            ALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells))
            VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1) = &
              tempCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1)
            VirtMergedCells(MasterCellID)%MergedCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells) = iElem
            VirtMergedCells(MasterCellID)%MergedVolume=VirtMergedCells(MasterCellID)%MergedVolume+ElemVolume_Shared(CNElemID) 
            VirtMergedCells(iElem)%MasterCell = VirtMergedCells(LocNBElem)%MasterCell
            VirtMergedCells(iElem)%isMerged = .TRUE.            
            DEALLOCATE(tempCellID)
            CYCLE ElemLoop
          END IF
        ELSE
          CYCLE NBElemLoop
        END IF
      ELSE IF ((VirtMergedCells(LocNBElem)%NumOfMergedCells.GT.0).AND.AllowBackMerge) THEN
        IF(VirtualCellMergeSpread.GT.0) THEN
          IF (VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
            MasterCellID = LocNBElem
            IF(VirtMergedCells(MasterCellID)%NumOfMergedCells.GE.(MaxNumOfMergedCells-1)) CYCLE NBElemLoop
            ALLOCATE(tempCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells))
            tempCellID = VirtMergedCells(MasterCellID)%MergedCellID
            DEALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID)
            VirtMergedCells(MasterCellID)%NumOfMergedCells = VirtMergedCells(MasterCellID)%NumOfMergedCells + 1 + VirtMergedCells(iElem)%NumOfMergedCells
            ALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells))
            VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1-VirtMergedCells(iElem)%NumOfMergedCells) = &
              tempCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1-VirtMergedCells(iElem)%NumOfMergedCells)
            oldMasterElemLoop2: DO iOldElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
              currentCellCount = VirtMergedCells(MasterCellID)%NumOfMergedCells-VirtMergedCells(iElem)%NumOfMergedCells -1 + iOldElem
              VirtMergedCells(MasterCellID)%MergedCellID(currentCellCount) = VirtMergedCells(iElem)%MergedCellID(iOldElem)
              VirtMergedCells(VirtMergedCells(iElem)%MergedCellID(iOldElem))%MasterCell = MasterCellID + offSetElem
            END DO oldMasterElemLoop2
            VirtMergedCells(MasterCellID)%MergedCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells) = iElem
            VirtMergedCells(iElem)%MasterCell = MasterCellID + offSetElem
            VirtMergedCells(iElem)%isMerged = .TRUE.
            VirtMergedCells(MasterCellID)%MergedVolume=VirtMergedCells(MasterCellID)%MergedVolume+VirtMergedCells(iElem)%MergedVolume
            VirtMergedCells(iElem)%MergedVolume = 0.0
            VirtMergedCells(iElem)%NumOfMergedCells = 0
            DEALLOCATE(VirtMergedCells(iElem)%MergedCellID)
            DEALLOCATE(tempCellID)
            CYCLE ElemLoop
          ELSE
            MasterCellID = LocNBElem
            IF(VirtMergedCells(MasterCellID)%NumOfMergedCells.GE.(MaxNumOfMergedCells-1)) CYCLE NBElemLoop
            ALLOCATE(tempCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells))
            tempCellID = VirtMergedCells(MasterCellID)%MergedCellID
            DEALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID)
            VirtMergedCells(MasterCellID)%NumOfMergedCells = VirtMergedCells(MasterCellID)%NumOfMergedCells + 1
            ALLOCATE(VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells))
            VirtMergedCells(MasterCellID)%MergedCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1) = &
              tempCellID(1:VirtMergedCells(MasterCellID)%NumOfMergedCells-1)
            VirtMergedCells(MasterCellID)%MergedCellID(VirtMergedCells(MasterCellID)%NumOfMergedCells) = iElem
            VirtMergedCells(MasterCellID)%MergedVolume=VirtMergedCells(MasterCellID)%MergedVolume+ElemVolume_Shared(CNElemID) 
            VirtMergedCells(iElem)%MasterCell = MasterCellID + offSetElem
            VirtMergedCells(iElem)%isMerged = .TRUE.
            DEALLOCATE(tempCellID)
            CYCLE ElemLoop
          END IF
        ELSE
          CYCLE NBElemLoop
        END IF
      ELSE IF (VirtMergedCells(LocNBElem)%isMerged.AND.(.NOT.AllowBackMerge)) THEN
        CYCLE NBElemLoop
      ELSE IF ((VirtMergedCells(LocNBElem)%NumOfMergedCells.GT.0).AND.(.NOT.AllowBackMerge)) THEN
        CYCLE NBElemLoop
      ELSE
        IF(VirtualCellMergeSpread.LT.3) AllowBackMerge = .FALSE.
        IF (VirtMergedCells(iElem)%NumOfMergedCells.EQ.0) THEN
          VirtMergedCells(iElem)%MergedVolume = VirtMergedCells(iElem)%MergedVolume + ElemVolume_Shared(CNElemID)
          VirtMergedCells(iElem)%NumOfMergedCells = VirtMergedCells(iElem)%NumOfMergedCells + 1
          ALLOCATE(VirtMergedCells(iElem)%MergedCellID(VirtMergedCells(iElem)%NumOfMergedCells))
          VirtMergedCells(iElem)%MergedCellID(VirtMergedCells(iElem)%NumOfMergedCells) = LocNBElem
          VirtMergedCells(iElem)%MergedVolume = VirtMergedCells(iElem)%MergedVolume + ElemVolume_Shared(CNNbElem)
          VirtMergedCells(iElem)%MasterCell = iElem + offSetElem  
          VirtMergedCells(LocNBElem)%MasterCell = iElem + offSetElem  
          VirtMergedCells(LocNBElem)%isMerged = .TRUE.
        ELSE
          IF(VirtMergedCells(iElem)%NumOfMergedCells.GE.(MaxNumOfMergedCells-1)) CYCLE ElemLoop
          ALLOCATE(tempCellID(VirtMergedCells(iElem)%NumOfMergedCells))
          tempCellID = VirtMergedCells(iElem)%MergedCellID
          DEALLOCATE(VirtMergedCells(iElem)%MergedCellID)
          VirtMergedCells(iElem)%NumOfMergedCells = VirtMergedCells(iElem)%NumOfMergedCells + 1
          ALLOCATE(VirtMergedCells(iElem)%MergedCellID(VirtMergedCells(iElem)%NumOfMergedCells))
          VirtMergedCells(iElem)%MergedCellID(1:VirtMergedCells(iElem)%NumOfMergedCells-1) = &
            tempCellID(1:VirtMergedCells(iElem)%NumOfMergedCells-1)
          VirtMergedCells(iElem)%MergedCellID(VirtMergedCells(iElem)%NumOfMergedCells) = LocNBElem
          VirtMergedCells(iElem)%MergedVolume = VirtMergedCells(iElem)%MergedVolume + ElemVolume_Shared(CNNbElem)
          VirtMergedCells(LocNBElem)%MasterCell = iElem + offSetElem  
          VirtMergedCells(LocNBElem)%isMerged = .TRUE.
          DEALLOCATE(tempCellID)
        END IF
        nPartMerged = nPartMerged + PEM%pNumber(LocNBElem)
        IF (nPartMerged.GT.MinPartNumCellMerge) CYCLE ElemLoop
      END IF
    END DO NBElemLoop    
  END IF
END DO ElemLoop

END SUBROUTINE MergeCells



SUBROUTINE InitializeParticleMaxwell(iPart,iSpec,iElem,Mode,iInit)
!===================================================================================================================================
!> Initialize a particle from a given macroscopic result, requires the macroscopic velocity, translational and internal temperatures
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: offSetElem
USE MOD_Particle_Vars           ,ONLY: PDM, PartSpecies, PartState, PEM, UseVarTimeStep, PartTimeStep, PartMPF, Species
USE MOD_DSMC_Vars               ,ONLY: DSMC, PartStateIntEn, CollisMode, SpecDSMC, RadialWeighting, AmbipolElecVelo
USE MOD_Restart_Vars            ,ONLY: MacroRestartValues
USE MOD_Particle_TimeStep       ,ONLY: GetParticleTimeStep
USE MOD_Particle_Emission_Vars  ,ONLY: EmissionDistributionDim
!USE MOD_part_tools              ,ONLY: CalcRadWeightMPF, CalcEElec_particle, CalcEVib_particle, CalcERot_particle
!USE MOD_part_tools              ,ONLY: CalcVelocity_maxwell_particle
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart, iSpec, iElem
INTEGER, INTENT(IN)             :: Mode                  !< 1: Macroscopic restart (data for each element)
                                                         !< 2: Emission distribution (equidistant data from .h5 file)
INTEGER, INTENT(IN), OPTIONAL   :: iInit                 !< particle emission initialization ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=100) :: hilf        !< auxiliary variable fr error output
INTEGER            :: iSpecAD     !< ambipolar diffusion species ID
REAL               :: v(3),vAD(3) !< velocity vector (vAD corresponds to ambipolar diffusion)
REAL               :: T(3),TAD(3) !< temperature vector (TAD corresponds to ambipolar diffusion)
REAL               :: Tvib        !< vibrational temperature
REAL               :: Trot        !< rotational temperature
REAL               :: Telec       !< electronic temperature
!===================================================================================================================================

SELECT CASE(Mode)
CASE(1) ! Macroscopic restart (data for each element)
  v       = MacroRestartValues(iElem,iSpec,1:3)
  T       = MacroRestartValues(iElem,iSpec,4:6)
  IF(DSMC%DoAmbipolarDiff)THEN
    iSpecAD = DSMC%AmbiDiffElecSpec
    vAD     = MacroRestartValues(iElem,iSpecAD,1:3)
    TAD     = MacroRestartValues(iElem,iSpecAD,4:6)
  END IF ! DSMC%DoAmbipolarDiff
  Tvib    = MacroRestartValues(iElem,iSpec,DSMC_TVIB)
  Trot    = MacroRestartValues(iElem,iSpec,DSMC_TROT)
  Telec   = MacroRestartValues(iElem,iSpec,DSMC_TELEC)
CASE(2) ! Emission distribution (equidistant data from .h5 file)
  ! Sanity check
  hilf=' is not implemented in InitializeParticleMaxwell() in combination with EmissionDistribution yet!'
  IF(DSMC%DoAmbipolarDiff) CALL abort(__STAMP__,'DSMC%DoAmbipolarDiff=T'//TRIM(hilf))
  IF(UseVarTimeStep) CALL abort(__STAMP__,'UseVarTimeStep=T'//TRIM(hilf))
  IF(RadialWeighting%DoRadialWeighting) CALL abort(__STAMP__,'RadialWeighting%DoRadialWeighting=T'//TRIM(hilf))
  ! Check dimensionality of data
  SELECT CASE(EmissionDistributionDim)
  CASE(1)
    CALL abort(__STAMP__,'EmissionDistributionDim=1 is not implemented')
  CASE(2)
    ! Density field from .h5 file that is interpolated to the element origin (bilinear interpolation)
    T(1:3) = InterpolateEmissionDistribution2D(iSpec,iInit,PartState(1:3,iPart),dimLower=4,dimUpper=4,transformation=.FALSE.)
    v(1:3) = InterpolateEmissionDistribution2D(iSpec,iInit,PartState(1:3,iPart),dimLower=5,dimUpper=6,transformation=.TRUE.)
  CASE(3)
    CALL abort(__STAMP__,'EmissionDistributionDim=3 is not implemented')
  END SELECT
  Tvib    = T(1)
  Trot    = T(1)
  Telec   = T(1)
CASE DEFAULT
  CALL abort(__STAMP__,'InitializeParticleMaxwell: Mode option is unknown. Mode=',IntInfoOpt=Mode)
END SELECT

! 1) Set particle velocity from macroscopic bulk velocity and translational temperature in the cell
PartState(4:6,iPart) = CalcVelocity_maxwell_particle(iSpec,T(1:3)) + v(1:3)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(iSpec)%ChargeIC.GT.0.0) THEN
    IF (ALLOCATED(AmbipolElecVelo(iPart)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(iPart)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(iPart)%ElecVelo(3))
    AmbipolElecVelo(iPart)%ElecVelo(1:3) = CalcVelocity_maxwell_particle(iSpecAD, TAD(1:3) ) + vAD(1:3)
  END IF
END IF
! 2) Set internal energies (rotational, vibrational, electronic)
IF(CollisMode.GT.1) THEN
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    PartStateIntEn(1,iPart) = CalcEVib_particle(iSpec,Tvib,iPart)
    PartStateIntEn(2,iPart) = CalcERot_particle(iSpec,Trot)
  ELSE
    PartStateIntEn(1:2,iPart) = 0.0
  END IF
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((SpecDSMC(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      PartStateIntEn(3,iPart) = CalcEElec_particle(iSpec,Telec,iPart)
    ELSE
      PartStateIntEn(3,iPart) = 0.0
    END IF
  END IF
END IF

! 3) Set the species and element number
PartSpecies(iPart) = iSpec
PEM%GlobalElemID(iPart) = iElem+offSetElem
PEM%LastGlobalElemID(iPart) = iElem+offSetElem
PDM%ParticleInside(iPart) = .TRUE.

! 4) Set particle time step and weights (if required)
IF (UseVarTimeStep) THEN
  PartTimeStep(iPart) = GetParticleTimeStep(PartState(1,iPart),PartState(2,iPart),iElem)
END IF
IF (RadialWeighting%DoRadialWeighting) THEN
  PartMPF(iPart) = CalcRadWeightMPF(PartState(2,iPart),iSpec,iPart)
END IF

END SUBROUTINE InitializeParticleMaxwell


PPURE FUNCTION InterpolateEmissionDistribution2D(iSpec,iInit,Pos,dimLower,dimUpper,transformation)
!===================================================================================================================================
!> This function returns a scalar property or a vector property (transformation from cylindrical coordinates to Cartesian required!)
!> Interpolates the variable external field to the r- and z-position via bilinear interpolation
!>
!>   1.2 --------- 2.2
!>    |             |
!>  r |             |
!>    |             |
!>   1.1 --------- 2.1
!>           z
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
!USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistribution !< data in 2D cylindrical coordinates
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionDelta
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionMin
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionMax,EmissionDistributionNum
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: iSpec, iInit
REAL,INTENT(IN)          :: Pos(1:3)       !< coordinate
LOGICAL,INTENT(IN)       :: transformation !< transform from f(r,z) to f(x,y,z)
INTEGER,INTENT(IN)       :: dimLower,dimUpper !< lower and upper dimension of variables in EmissionDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateEmissionDistribution2D(dimLower:dimLower+2)  !< 3D Cartesian vector
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos,jPos                                !< index in array (equidistant subdivision assumed)
REAL                     :: f(dimLower:dimUpper)
REAL                     :: r,delta,mat(2,2),dx(1:2),dy(1:2),vec(1:2)
INTEGER                  :: idx1,idx2,idx3,idx4,i
!===================================================================================================================================
ASSOCIATE(&
      x => Pos(1) ,&
      y => Pos(2) ,&
      z => Pos(3) ,&
      EmissionDistribution => Species(iSpec)%Init(iInit)%EmissionDistribution&
      )
  r = SQRT(x*x + y*y)

  IF(r.GT.EmissionDistributionMax(1))THEN
    InterpolateEmissionDistribution2D = 0.
  ELSEIF(r.LT.EmissionDistributionMin(1))THEN
    InterpolateEmissionDistribution2D = 0.
  ELSEIF(z.GT.EmissionDistributionMax(2))THEN
    InterpolateEmissionDistribution2D = 0.
  ELSEIF(z.LT.EmissionDistributionMin(2))THEN
    InterpolateEmissionDistribution2D = 0.
  ELSE

    ! Get index in r and z
    iPos = INT((r-EmissionDistribution(1,1))/EmissionDistributionDelta(1)) + 1 ! dr = EmissionDistributionDelta(1)
    jPos = INT((z-EmissionDistribution(2,1))/EmissionDistributionDelta(2)) + 1 ! dz = EmissionDistributionDelta(2)
    ! Catch problem when r or z are exactly at the upper boundary and INT() does not round to the lower integer (do not add +1 in
    ! this case)
    iPos = MIN(iPos, EmissionDistributionNum(2) - 1 )
    jPos = MIN(jPos, EmissionDistributionNum(1) - 1 )

    ! Shift all points by Nz = EmissionDistributionNum(1)
    ASSOCIATE( Nz => EmissionDistributionNum(1) )
      ! 1.1
      idx1 = (iPos-1)*Nz + jPos
      ! 2.1
      idx2 = (iPos-1)*Nz + jPos + 1
      ! 1.2
      idx3 = iPos*Nz + jPos
      ! 2.2
      idx4 = iPos*Nz + jPos + 1
    END ASSOCIATE

    ! Interpolate
    delta = EmissionDistributionDelta(1)*EmissionDistributionDelta(2)
    delta = 1./delta

    dx(1) = EmissionDistribution(1,idx4)-r
    dx(2) = EmissionDistributionDelta(1) - dx(1)

    dy(1) = EmissionDistribution(2,idx4)-z
    dy(2) = EmissionDistributionDelta(2) - dy(1)

    DO i = dimLower, dimUpper
      mat(1,1) = EmissionDistribution(i,idx1)
      mat(2,1) = EmissionDistribution(i,idx2)
      mat(1,2) = EmissionDistribution(i,idx3)
      mat(2,2) = EmissionDistribution(i,idx4)

      vec(1) = dx(1)
      vec(2) = dx(2)
      f(i) = delta * DOT_PRODUCT( vec, MATMUL(mat, (/dy(1),dy(2)/) ) )
    END DO ! i =  dimLower, dimUpper

    ! Transform from Br, Bz to Bx, By, Bz
    IF(transformation)THEN
      r=1./r
      InterpolateEmissionDistribution2D(dimLower)   = f(dimLower)*x*r
      InterpolateEmissionDistribution2D(dimLower+1) = f(dimLower)*y*r
      InterpolateEmissionDistribution2D(dimLower+2) = f(dimUpper)
    ELSE
      InterpolateEmissionDistribution2D(dimLower:dimUpper) = f(dimLower:dimUpper)
      InterpolateEmissionDistribution2D(dimLower+1) = InterpolateEmissionDistribution2D(dimLower)
      InterpolateEmissionDistribution2D(dimLower+2) = InterpolateEmissionDistribution2D(dimLower)
    END IF ! transformation
  END IF ! r.GT.EmissionDistributionMax(1)
END ASSOCIATE

END FUNCTION InterpolateEmissionDistribution2D

END MODULE MOD_part_tools
