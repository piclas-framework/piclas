!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

MODULE MOD_Particle_Boundary_Tools
!===================================================================================================================================
! Tools used for boundary interactions
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
PUBLIC :: CalcWallSample
PUBLIC :: SampleImpactProperties
PUBLIC :: StoreBoundaryParticleProperties
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcWallSample(PartID,SurfSideID,p,q,SampleType,PartTrajectory_opt,SurfaceNormal_opt)
!===================================================================================================================================
!> Sample the energy of particles before and after a wall interaction for the determination of macroscopic properties such as heat
!> flux and force per area
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Globals                   ,ONLY: abort,DOTPRODUCT
USE MOD_DSMC_Vars                 ,ONLY: SpecDSMC,useDSMC,PartStateIntEn,RadialWeighting
USE MOD_DSMC_Vars                 ,ONLY: CollisMode,DSMC,AmbipolElecVelo
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallState,CalcSurfaceImpact
USE MOD_part_tools                ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: PartID,SurfSideID,p,q
CHARACTER(*),INTENT(IN)            :: SampleType
REAL,INTENT(IN),OPTIONAL           :: PartTrajectory_opt(1:3)
REAL,INTENT(IN),OPTIONAL           :: SurfaceNormal_opt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: ETrans, ETransAmbi, ERot, EVib, EElec, MomArray(1:3), MassIC, MPF
INTEGER         :: ETransID, ERotID, EVibID, EElecID, SpecID
!===================================================================================================================================
MomArray(:)=0.
EVib = 0.
ERot = 0.
EElec = 0.

SpecID = PartSpecies(PartID)

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  MPF = GetParticleWeight(PartID)
ELSE
  MPF = GetParticleWeight(PartID)*Species(SpecID)%MacroParticleFactor
END IF
MassIC = Species(SpecID)%MassIC

! Calculate the translational energy
ETrans = 0.5 * Species(SpecID)%MassIC * DOTPRODUCT(PartState(4:6,PartID))
IF (DSMC%DoAmbipolarDiff) THEN
  ! Add the translational energy of electron "attached" to the ion
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    ETransAmbi = 0.5 * Species(DSMC%AmbiDiffElecSpec)%MassIC * DOTPRODUCT(AmbipolElecVelo(PartID)%ElecVelo(1:3))
    ! Save the electron energy to sample it later in SampleImpactProperties
    ETrans = ETrans + ETransAmbi
  END IF
END IF
! Depending on whether the routine is called before (old) or after (new) a surface interaction, the momentum is added or removed
! from the sampling array. Additionally, the correct indices are set for the sampling array.
SELECT CASE (TRIM(SampleType))
CASE ('old')
  MomArray(1:3)   = MassIC * PartState(4:6,PartID) * MPF
  ETransID = SAMPWALL_ETRANSOLD
  ERotID   = SAMPWALL_EROTOLD
  EVibID   = SAMPWALL_EVIBOLD
  EElecID  = SAMPWALL_EELECOLD
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
      MomArray(1:3) = MomArray(1:3) + Species(DSMC%AmbiDiffElecSpec)%MassIC * AmbipolElecVelo(PartID)%ElecVelo(1:3) * MPF
    END IF
  END IF
  ! Species-specific simulation particle impact counter
  SampWallState(SAMPWALL_NVARS+SpecID,p,q,SurfSideID) = SampWallState(SAMPWALL_NVARS+SpecID,p,q,SurfSideID) + 1
  ! Sampling of species-specific impact energies and angles
  IF(CalcSurfaceImpact) THEN
    IF (useDSMC) THEN
      IF (CollisMode.GT.1) THEN
        EVib = PartStateIntEn(1,PartID)
        ERot = PartStateIntEn(2,PartID)
        IF(DSMC%ElectronicModel) THEN
          EElec = PartStateIntEn(3,PartID)
        END IF
      END IF
    END IF
    CALL SampleImpactProperties(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,EElec,PartTrajectory_opt,SurfaceNormal_opt,p,q)
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
        CALL SampleImpactProperties(SurfSideID,DSMC%AmbiDiffElecSpec,MPF,ETransAmbi,0.,0.,0.,PartTrajectory_opt,SurfaceNormal_opt,p,q)
      END IF
    END IF
  END IF
  ! Sample the time step for the correct determination of the heat flux
  IF (VarTimeStep%UseVariableTimeStep) THEN
    SampWallState(SAMPWALL_NVARS+nSpecies+1,p,q,SurfSideID) = SampWallState(SAMPWALL_NVARS+nSpecies+1,p,q,SurfSideID) &
                                                              + VarTimeStep%ParticleTimeStep(PartID)
  END IF
CASE ('new')
  ! must be old_velocity-new_velocity
  MomArray(1:3)   = -MassIC * PartState(4:6,PartID) * MPF
  ETransID = SAMPWALL_ETRANSNEW
  ERotID   = SAMPWALL_EROTNEW
  EVibID   = SAMPWALL_EVIBNEW
  EElecID  = SAMPWALL_EELECNEW
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
      MomArray(1:3) = MomArray(1:3) - Species(DSMC%AmbiDiffElecSpec)%MassIC * AmbipolElecVelo(PartID)%ElecVelo(1:3) * MPF
    END IF
  END IF
CASE DEFAULT
  CALL abort(&
    __STAMP__&
    ,'ERROR in CalcWallSample: wrong SampleType specified. Possible types -> ( old , new )')
END SELECT
!----  Sampling force at walls (correct sign is set above)
SampWallState(SAMPWALL_DELTA_MOMENTUMX,p,q,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMX,p,q,SurfSideID) + MomArray(1)
SampWallState(SAMPWALL_DELTA_MOMENTUMY,p,q,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMY,p,q,SurfSideID) + MomArray(2)
SampWallState(SAMPWALL_DELTA_MOMENTUMZ,p,q,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMZ,p,q,SurfSideID) + MomArray(3)
!----  Sampling the energy (translation) accommodation at walls
SampWallState(ETransID ,p,q,SurfSideID) = SampWallState(ETransID ,p,q,SurfSideID) + ETrans * MPF
IF (useDSMC) THEN
  IF (CollisMode.GT.1) THEN
    IF ((SpecDSMC(SpecID)%InterID.EQ.2).OR.SpecDSMC(SpecID)%InterID.EQ.20) THEN
      !----  Sampling the internal (rotational) energy accommodation at walls
      SampWallState(ERotID ,p,q,SurfSideID) = SampWallState(ERotID ,p,q,SurfSideID) + PartStateIntEn(2,PartID) * MPF
      !----  Sampling for internal (vibrational) energy accommodation at walls
      SampWallState(EVibID ,p,q,SurfSideID) = SampWallState(EVibID ,p,q,SurfSideID) + PartStateIntEn(1,PartID) * MPF
    END IF
    IF(DSMC%ElectronicModel) THEN
      !----  Sampling for internal (electronic) energy accommodation at walls
      SampWallState(EElecID ,p,q,SurfSideID) = SampWallState(EElecID ,p,q,SurfSideID) + PartStateIntEn(3,PartID) * MPF
    END IF
  END IF
END IF

END SUBROUTINE CalcWallSample


SUBROUTINE SampleImpactProperties(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,EElec,PartTrajectory,SurfaceNormal,p,q)
!===================================================================================================================================
!> Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!>
!===================================================================================================================================
!USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_Particle_Boundary_Vars ,ONLY: SampWallImpactEnergy,SampWallImpactVector
USE MOD_Particle_Boundary_Vars ,ONLY: SampWallImpactAngle ,SampWallImpactNumber
USE MOD_Globals_Vars           ,ONLY: PI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: SurfSideID          !< Surface ID
INTEGER,INTENT(IN) :: SpecID              !< Species ID
REAL,INTENT(IN)    :: MPF                 !< Particle macro particle factor
REAL,INTENT(IN)    :: ETrans              !< Translational energy of impacting particle
REAL,INTENT(IN)    :: ERot                !< Rotational energy of impacting particle
REAL,INTENT(IN)    :: EVib                !< Vibrational energy of impacting particle
REAL,INTENT(IN)    :: EElec               !< Electronic energy of impacting particle
REAL,INTENT(IN)    :: PartTrajectory(1:3) !< Particle trajectory vector (normalized)
REAL,INTENT(IN)    :: SurfaceNormal(1:3)  !< Surface normal vector (normalized)
INTEGER,INTENT(IN) :: p                 !< Surface sub-faces
INTEGER,INTENT(IN) :: q                 !< Surface sub-faces
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

!----- Sampling of impact energy for each species (trans, rot, vib)
SampWallImpactEnergy(SpecID,1,p,q,SurfSideID) = SampWallImpactEnergy(SpecID,1,p,q,SurfSideID) + ETrans * MPF
SampWallImpactEnergy(SpecID,2,p,q,SurfSideID) = SampWallImpactEnergy(SpecID,2,p,q,SurfSideID) + ERot   * MPF
SampWallImpactEnergy(SpecID,3,p,q,SurfSideID) = SampWallImpactEnergy(SpecID,3,p,q,SurfSideID) + EVib   * MPF
SampWallImpactEnergy(SpecID,4,p,q,SurfSideID) = SampWallImpactEnergy(SpecID,4,p,q,SurfSideID) + EElec  * MPF

!----- Sampling of impact vector ,SurfSideIDfor each species (x,y,z)
SampWallImpactVector(SpecID,1,p,q,SurfSideID) = SampWallImpactVector(SpecID,1,p,q,SurfSideID) + PartTrajectory(1) * MPF
SampWallImpactVector(SpecID,2,p,q,SurfSideID) = SampWallImpactVector(SpecID,2,p,q,SurfSideID) + PartTrajectory(2) * MPF
SampWallImpactVector(SpecID,3,p,q,SurfSideID) = SampWallImpactVector(SpecID,3,p,q,SurfSideID) + PartTrajectory(3) * MPF

!----- Sampling of impact angle for each species
SampWallImpactAngle(SpecID,p,q,SurfSideID) = SampWallImpactAngle(SpecID,p,q,SurfSideID) + &
    (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal)))) * MPF

!----- Sampling of impact number for each species
SampWallImpactNumber(SpecID,p,q,SurfSideID) = SampWallImpactNumber(SpecID,p,q,SurfSideID) + MPF

END SUBROUTINE SampleImpactProperties


SUBROUTINE StoreBoundaryParticleProperties(iPart,SpecID,PartPos,PartTrajectory,SurfaceNormal,mode,usevMPF_optIN)
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity and species to PartDataBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,Species,PartState
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength
USE MOD_TimeDisc_Vars          ,ONLY: time
USE MOD_Globals_Vars           ,ONLY: PI
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iPart
INTEGER,INTENT(IN) :: SpecID ! The species ID is required as it might not yet be set during emission
REAL,INTENT(IN)    :: PartPos(1:3)
REAL,INTENT(IN)    :: PartTrajectory(1:3)
REAL,INTENT(IN)    :: SurfaceNormal(1:3)
INTEGER,INTENT(IN) :: mode ! 1: particle impacts on BC (species is stored as positive value)
                             ! 2: particles is emitted from the BC into the simulation domain (species is stored as negative value)
LOGICAL,INTENT(IN),OPTIONAL :: usevMPF_optIN ! For setting MPF for cases when PartMPF(iPart) might not yet be set during emission
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF
INTEGER              :: dims(2)
! Temporary arrays
REAL, ALLOCATABLE    :: PartStateBoundary_tmp(:,:) ! (1:10,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle
!                                                  !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: ALLOCSTAT
!===================================================================================================================================
IF(PRESENT(usevMPF_optIN))THEN
  IF(usevMPF_optIN)THEN
    CALL abort(&
    __STAMP__&
    ,'StoreBoundaryParticleProperties: usevMPF_optIN cannot be true!')
  ELSE
    MPF = Species(SpecID)%MacroParticleFactor
  END IF ! usevMPF_optIN
ELSE
  IF (usevMPF) THEN
    MPF = PartMPF(iPart)
  ELSE
    MPF = Species(SpecID)%MacroParticleFactor
  END IF
END IF ! PRESENT(MPF_optIN)

dims = SHAPE(PartStateBoundary)

ASSOCIATE( iMax => PartStateBoundaryVecLength )
  ! Increase maximum number of boundary-impact particles
  iMax = iMax + 1

  ! Check if array maximum is reached.
  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
  IF(iMax.GT.dims(2))THEN

    ! --- PartStateBoundary ---
    ALLOCATE(PartStateBoundary_tmp(1:10,1:dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:10,1:dims(2)) = PartStateBoundary(1:10,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    ALLOCATE(PartStateBoundary(1:10,1:2*dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:10,        1:  dims(2)) = PartStateBoundary_tmp(1:10,1:dims(2))
    PartStateBoundary(1:10,dims(2)+1:2*dims(2)) = 0.

  END IF

  PartStateBoundary(1:3,iMax) = PartPos
  PartStateBoundary(4:6,iMax) = PartState(4:6,iPart)
  IF(mode.EQ.1)THEN
    PartStateBoundary(7  ,iMax) = REAL(SpecID)
  ELSEIF(mode.EQ.2)THEN
    PartStateBoundary(7  ,iMax) = -REAL(SpecID)
  ELSE
    CALL abort(&
    __STAMP__&
    ,'StoreBoundaryParticleProperties: mode must be either 1 or 2! mode=',IntInfoOpt=mode)
  END IF ! mode.EQ.1
  PartStateBoundary(8  ,iMax) = MPF
  PartStateBoundary(9  ,iMax) = time
  PartStateBoundary(10 ,iMax) = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal))))
END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties

END MODULE MOD_Particle_Boundary_Tools