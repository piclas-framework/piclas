!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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
PUBLIC :: StoreBoundaryParticleProperties
PUBLIC :: GetRadialDistance2D,GetRacetrackDistance2D
PUBLIC :: PointToSegmentDist2D
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcWallSample(PartID,SurfSideID,SampleType,SurfaceNormal_opt,PartPosImpact_opt)
!===================================================================================================================================
!> Sample the energy of particles before and after a wall interaction for the determination of macroscopic properties such as heat
!> flux and force per area
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Globals                   ,ONLY: abort,DOTPRODUCT
USE MOD_DSMC_Vars                 ,ONLY: useDSMC,PartIntEn, SpecDSMC
USE MOD_DSMC_Vars                 ,ONLY: CollisMode,DSMC
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallState,CalcSurfaceImpact
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: CalcTorque, SWITorqueCoefficientX, SWITorqueCoefficientY, SWITorqueCoefficientZ
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcSurfOutputPerGroup
#if USE_HDG
USE MOD_Particle_Boundary_Vars    ,ONLY: DoVirtualDielectricLayer
USE MOD_Particle_Vars             ,ONLY: ResetVDLSpecID
#endif/*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: PartID,SurfSideID
CHARACTER(*),INTENT(IN)            :: SampleType
REAL,INTENT(IN),OPTIONAL           :: SurfaceNormal_opt(1:3),PartPosImpact_opt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: ETrans, ETransAmbi, ERot, EVib, EElec, MomArray(1:3), MassIC, PartWeight, TorqueArray(1:3)
INTEGER         :: ETransID, ERotID, EVibID, EElecID, SpecID, SubP, SubQ
!===================================================================================================================================
MomArray(:)=0.
EVib = 0.
ERot = 0.
EElec = 0.

SubP = TrackInfo%p
SubQ = TrackInfo%q

SpecID = PartSpecies(PartID)
MassIC = Species(SpecID)%MassIC
#if USE_HDG
! Check particle index for VDL particles and reset to original species index
IF(DoVirtualDielectricLayer) SpecID = ResetVDLSpecID(PartID)
#endif/*USE_HDG*/

! Particle weight for the sampling of rates (heat flux, force, torque, impacts per second):
! Using only the actual particle weight without time step (which is included in GetParticleWeight function)
! The surface sampling is normalized with the global sampling time, analogous to the surface flux, where the particle time step is
! not considered, as the flux between elements stays constant. Note the difference to the species-specific time step, where for the
! respective species the elapsed time corresponds to TimeSample * Species(SpecID)%TimeStepFactor, thus including it in the weight here
IF(usevMPF) THEN
  PartWeight = PartMPF(PartID)
ELSE
  PartWeight = Species(SpecID)%MacroParticleFactor
END IF
IF(VarTimeStep%UseSpeciesSpecific) PartWeight = PartWeight / Species(SpecID)%TimeStepFactor

! Calculate the translational energy
ETrans = 0.5 * Species(SpecID)%MassIC * DOTPRODUCT(PartState(4:6,PartID))
IF (DSMC%DoAmbipolarDiff) THEN
  ! Add the translational energy of electron "attached" to the ion
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    ETransAmbi = 0.5 * Species(DSMC%AmbiDiffElecSpec)%MassIC * DOTPRODUCT(PartIntEn(PartID)%ElecVelo(1:3))
    ! Save the electron energy to sample it later in SampleImpactProperties
    ETrans = ETrans + ETransAmbi
  END IF
END IF
! Depending on whether the routine is called before (old) or after (new) a surface interaction, the momentum is added or removed
! from the sampling array. Additionally, the correct indices are set for the sampling array.
SELECT CASE (TRIM(SampleType))
CASE ('old')
  MomArray(1:3)   = MassIC * PartState(4:6,PartID) * PartWeight
  ETransID = SAMPWALL_ETRANSOLD
  ERotID   = SAMPWALL_EROTOLD
  EVibID   = SAMPWALL_EVIBOLD
  EElecID  = SAMPWALL_EELECOLD
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
      MomArray(1:3) = MomArray(1:3) + Species(DSMC%AmbiDiffElecSpec)%MassIC * PartIntEn(PartID)%ElecVelo(1:3) * PartWeight
    END IF
  END IF
  ! Species-specific simulation particle impact counter
  SampWallState(SAMPWALL_NVARS+SpecID,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_NVARS+SpecID,SubP,SubQ,SurfSideID) + 1
  ! Sampling of species-specific impact energies and angles
  IF(CalcSurfaceImpact) THEN
    IF (useDSMC) THEN
      IF (CollisMode.GT.1) THEN
        IF((Species(SpecID)%InterID.EQ.2).OR.(Species(SpecID)%InterID.EQ.20)) THEN
          EVib = PartIntEn(PartID)%EVib(1)
          ERot = PartIntEn(PartID)%ERot(1)
        END IF
        IF(DSMC%ElectronicModel.GT.0) THEN
          IF((Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized)) EElec = PartIntEn(PartID)%EElec(1)
        END IF
      END IF
    END IF
    CALL SampleImpactProperties(SurfSideID,SpecID,PartWeight,ETrans,EVib,ERot,EElec,TrackInfo%PartTrajectory,SurfaceNormal_opt)
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
        CALL SampleImpactProperties(SurfSideID,DSMC%AmbiDiffElecSpec,PartWeight,ETransAmbi,0.,0.,0.,TrackInfo%PartTrajectory,SurfaceNormal_opt)
      END IF
    END IF
  END IF
CASE ('new')
  ! must be old_velocity-new_velocity
  MomArray(1:3)   = -MassIC * PartState(4:6,PartID) * PartWeight
  ETransID = SAMPWALL_ETRANSNEW
  ERotID   = SAMPWALL_EROTNEW
  EVibID   = SAMPWALL_EVIBNEW
  EElecID  = SAMPWALL_EELECNEW
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
      MomArray(1:3) = MomArray(1:3) - Species(DSMC%AmbiDiffElecSpec)%MassIC * PartIntEn(PartID)%ElecVelo(1:3) * PartWeight
    END IF
  END IF
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR in CalcWallSample: wrong SampleType specified. Possible types -> ( old , new )')
END SELECT
!----  Sampling force at walls (correct sign is set above)
SampWallState(SAMPWALL_DELTA_MOMENTUMX,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMX,SubP,SubQ,SurfSideID) + MomArray(1)
SampWallState(SAMPWALL_DELTA_MOMENTUMY,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMY,SubP,SubQ,SurfSideID) + MomArray(2)
SampWallState(SAMPWALL_DELTA_MOMENTUMZ,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMZ,SubP,SubQ,SurfSideID) + MomArray(3)
!----  Sampling the energy (translation) accommodation at walls
SampWallState(ETransID ,SubP,SubQ,SurfSideID) = SampWallState(ETransID ,SubP,SubQ,SurfSideID) + ETrans * PartWeight
!----  Sampling torque
IF(CalcTorque) THEN
  TorqueArray(1) = PartPosImpact_opt(2) * MomArray(3) - PartPosImpact_opt(3) * MomArray(2)
  TorqueArray(2) = PartPosImpact_opt(3) * MomArray(1) - PartPosImpact_opt(1) * MomArray(3)
  TorqueArray(3) = PartPosImpact_opt(1) * MomArray(2) - PartPosImpact_opt(2) * MomArray(1)
  SampWallState(SWITorqueCoefficientX,SubP,SubQ,SurfSideID) = SampWallState(SWITorqueCoefficientX,SubP,SubQ,SurfSideID) + TorqueArray(1)
  SampWallState(SWITorqueCoefficientY,SubP,SubQ,SurfSideID) = SampWallState(SWITorqueCoefficientY,SubP,SubQ,SurfSideID) + TorqueArray(2)
  SampWallState(SWITorqueCoefficientZ,SubP,SubQ,SurfSideID) = SampWallState(SWITorqueCoefficientZ,SubP,SubQ,SurfSideID) + TorqueArray(3)
END IF
IF (useDSMC) THEN
  IF (CollisMode.GT.1) THEN
    IF ((Species(SpecID)%InterID.EQ.2).OR.Species(SpecID)%InterID.EQ.20) THEN
      !----  Sampling the internal (rotational) energy accommodation at walls
      SampWallState(ERotID ,SubP,SubQ,SurfSideID) = SampWallState(ERotID ,SubP,SubQ,SurfSideID) + PartIntEn(PartID)%ERot(1) * PartWeight
      !----  Sampling for internal (vibrational) energy accommodation at walls
      SampWallState(EVibID ,SubP,SubQ,SurfSideID) = SampWallState(EVibID ,SubP,SubQ,SurfSideID) + PartIntEn(PartID)%EVib(1) * PartWeight
    END IF
    IF(DSMC%ElectronicModel.GT.0) THEN
      !----  Sampling for internal (electronic) energy accommodation at walls
      IF((Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized)) &
        SampWallState(EElecID ,SubP,SubQ,SurfSideID) = SampWallState(EElecID ,SubP,SubQ,SurfSideID) + PartIntEn(PartID)%EElec(1) * PartWeight
    END IF
  END IF
END IF
!---- Sampling of integral group output for SurfaceAnalyze.csv
IF(CalcSurfOutputPerGroup) THEN
  CALL SampleSurfaceGroupProperties(SurfSideID,PartID,SpecID,SampleType,TorqueArray,ETrans,PartWeight)
END IF

END SUBROUTINE CalcWallSample


SUBROUTINE SampleImpactProperties(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,EElec,PartTrajectory,SurfaceNormal)
!===================================================================================================================================
!> Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!===================================================================================================================================
USE MOD_Particle_Boundary_Vars ,ONLY: SampWallImpactEnergy,SampWallImpactVector
USE MOD_Particle_Boundary_Vars ,ONLY: SampWallImpactAngle ,SampWallImpactNumber
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_Particle_Tracking_Vars ,ONLY: TrackInfo
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SubP, SubQ
!-----------------------------------------------------------------------------------------------------------------------------------

SubP = TrackInfo%p
SubQ = TrackInfo%q

!----- Sampling of impact energy for each species (trans, rot, vib)
SampWallImpactEnergy(SpecID,1,SubP,SubQ,SurfSideID) = SampWallImpactEnergy(SpecID,1,SubP,SubQ,SurfSideID) + ETrans * MPF
SampWallImpactEnergy(SpecID,2,SubP,SubQ,SurfSideID) = SampWallImpactEnergy(SpecID,2,SubP,SubQ,SurfSideID) + ERot   * MPF
SampWallImpactEnergy(SpecID,3,SubP,SubQ,SurfSideID) = SampWallImpactEnergy(SpecID,3,SubP,SubQ,SurfSideID) + EVib   * MPF
SampWallImpactEnergy(SpecID,4,SubP,SubQ,SurfSideID) = SampWallImpactEnergy(SpecID,4,SubP,SubQ,SurfSideID) + EElec  * MPF

!----- Sampling of impact vector ,SurfSideIDfor each species (x,y,z)
SampWallImpactVector(SpecID,1,SubP,SubQ,SurfSideID) = SampWallImpactVector(SpecID,1,SubP,SubQ,SurfSideID) + PartTrajectory(1) * MPF
SampWallImpactVector(SpecID,2,SubP,SubQ,SurfSideID) = SampWallImpactVector(SpecID,2,SubP,SubQ,SurfSideID) + PartTrajectory(2) * MPF
SampWallImpactVector(SpecID,3,SubP,SubQ,SurfSideID) = SampWallImpactVector(SpecID,3,SubP,SubQ,SurfSideID) + PartTrajectory(3) * MPF

!----- Sampling of impact angle for each species
SampWallImpactAngle(SpecID,SubP,SubQ,SurfSideID) = SampWallImpactAngle(SpecID,SubP,SubQ,SurfSideID) + &
    (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal)))) * MPF

!----- Sampling of impact number for each species
SampWallImpactNumber(SpecID,SubP,SubQ,SurfSideID) = SampWallImpactNumber(SpecID,SubP,SubQ,SurfSideID) + MPF

END SUBROUTINE SampleImpactProperties


!----------------------------------------------------------------------------------------------------------------------------------!
!> Save particle position, velocity and species to PartDataBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE StoreBoundaryParticleProperties(iPart,SpecID,PartPos,PartTrajectory,SurfaceNormal,iPartBound,mode,MPF_optIN,Velo_optIN)
! MODULES
USE MOD_Globals
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,Species,PartState
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundaryMemoryLimit,PartStateBoundaryMemory,PartStateBoundaryResizeCounter
USE MOD_TimeDisc_Vars          ,ONLY: time
USE MOD_Globals_Vars           ,ONLY: PI, Joule2eV
USE MOD_Array_Operations       ,ONLY: ChangeSizeArray
USE MOD_Particle_Analyze_Pure  ,ONLY: CalcEkinPart2
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iPart
INTEGER,INTENT(IN) :: SpecID !> The species ID is required as it might not yet be set during emission
REAL,INTENT(IN)    :: PartPos(1:3)
REAL,INTENT(IN)    :: PartTrajectory(1:3)
REAL,INTENT(IN)    :: SurfaceNormal(1:3)
INTEGER,INTENT(IN) :: iPartBound  !> Part-BoundaryX on which the impact occurs
INTEGER,INTENT(IN) :: mode !> 1: particle impacts on BC (species is stored as positive value)
                           !> 2: particles is emitted from the BC into the simulation domain (species is stored as negative value)
REAL,INTENT(IN),OPTIONAL :: MPF_optIN !> Supply the MPF in special cases
REAL,INTENT(IN),OPTIONAL :: Velo_optIN(1:3) !> Supply the velocity in special cases
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF,PartStateBoundaryMemUsage
INTEGER              :: dims(2)
!===================================================================================================================================
IF(PRESENT(MPF_optIN))THEN
  MPF = MPF_optIN
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
  ! Check if array maximum is reached and increase size if it does
  IF(iMax.GT.dims(2))THEN
    ! Utilizing routine using MOVE_ALLOC and increase array by 20%
    CALL ChangeSizeArray(PartStateBoundary,dims(2),CEILING(dims(2)*1.2),0.)
    ! Increment counter
    PartStateBoundaryResizeCounter = PartStateBoundaryResizeCounter + 1
    ! Check if the memory limit per core is reached
    PartStateBoundaryMemUsage = PartStateBoundaryMemory*1.2**(PartStateBoundaryResizeCounter-1)
    IF (PartStateBoundaryMemUsage.GT.PartStateBoundaryMemoryLimit) THEN
      ! Output warning
      IPWRITE(UNIT_StdOut,'(I0,A,I0,A,ES8.2,A,ES8.2,A,I0,A)') ' Warning: PartStateBoundary has been resized ',&
        PartStateBoundaryResizeCounter,' times to [',&
        PartStateBoundaryMemUsage,'] GB, which has passed the process limit of [',&
        PartStateBoundaryMemoryLimit, '] GB (',dims(2),' particles)'
    END IF
  END IF

  PartStateBoundary(1:3,iMax) = PartPos
  IF(PRESENT(Velo_optIN))THEN
    PartStateBoundary(4:6,iMax) = Velo_optIN
  ELSE
    PartStateBoundary(4:6,iMax) = PartState(4:6,iPart)
  END IF ! PRESENT(Velo_optIN)
  ! Mode 1: store normal species ID, mode 2: store negative species ID (special analysis of emitted particles in/from volume/surface)
  IF(mode.EQ.1)THEN
    PartStateBoundary(7  ,iMax) = REAL(SpecID)
  ELSEIF(mode.EQ.2)THEN
    PartStateBoundary(7  ,iMax) = -REAL(SpecID)
  ELSE
    CALL abort(__STAMP__,'StoreBoundaryParticleProperties: mode must be either 1 or 2! mode=',IntInfoOpt=mode)
  END IF ! mode.EQ.1
  IF(PartStateBoundary(7,iMax).EQ.0) CALL abort(__STAMP__,'Error in StoreBoundaryParticleProperties. SpecID is zero')
  ! Calculate kinetic energy (or set to zero for photon debug output)
  IF(SpecID.EQ.999)THEN
    PartStateBoundary(8  ,iMax) = 0.
  ELSE
    PartStateBoundary(8  ,iMax) = CalcEkinPart2(PartStateBoundary(4:6,iMax),SpecID,1.0) * Joule2eV
  END IF
  PartStateBoundary(9  ,iMax) = MPF
  PartStateBoundary(10 ,iMax) = time
  PartStateBoundary(11 ,iMax) = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal))))
  PartStateBoundary(12 ,iMax) = REAL(iPartBound)

END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties


!===================================================================================================================================
!> Determines the minimum and maximum radial distance from a side's bounding box to a given origin on a surface.
!>
!>    corner 4 (min,max) -------- corner 3 (max,max)
!>           |                           |
!>           |                           |
!>           |                           |
!>    corner 1 (min,min) -------- corner 2 (max,min)
!>
!===================================================================================================================================
SUBROUTINE GetRadialDistance2D(GlobalSideID,dir,origin,rmin,rmax)
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces       ,ONLY: GetSideBoundingBox
USE MOD_Particle_Mesh_Tools     ,ONLY: GetSideBoundingBoxTria
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: GlobalSideID, dir(3)
REAL, INTENT(IN)              :: origin(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: rmin,rmax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode, iNext
REAL                          :: BoundingBox(1:3,1:8)
REAL                          :: corners(2,4)               !> Bounding box corners in origin-shifted 2D
REAL                          :: point(2), vec(2), dist
LOGICAL                       :: r0inside
!===================================================================================================================================
! Get bounding box
IF (TrackingMethod.EQ.TRIATRACKING) THEN
  CALL GetSideBoundingBoxTria(GlobalSideID,BoundingBox)
ELSE
  CALL GetSideBoundingBox(GlobalSideID,BoundingBox)
END IF
IF(Symmetry%Axisymmetric) THEN
  ! Store the y-coordinate (=2) of bounding box only: the first two nodes have yMin (= 1), and the third has yMax (= 3)
  rmin = BoundingBox(2,1)
  rmax = BoundingBox(2,3)
ELSE
  ! Extract 4 bounding box corners in origin-shifted 2D surface coordinates
  ! Corner ordering: 1 (min,min), 2 (max,min), 3 (max,max), 4 (min,max)
  corners(1,1) = MINVAL(BoundingBox(dir(2),:)) - origin(1)
  corners(2,1) = MINVAL(BoundingBox(dir(3),:)) - origin(2)
  corners(1,2) = MAXVAL(BoundingBox(dir(2),:)) - origin(1)
  corners(2,2) = corners(2,1)
  corners(1,3) = corners(1,2)
  corners(2,3) = MAXVAL(BoundingBox(dir(3),:)) - origin(2)
  corners(1,4) = corners(1,1)
  corners(2,4) = corners(2,3)

  !-- rmax: maximum distance to the origin, which will always be at the corners
  rmax = 0.
  DO iNode = 1, 4
    dist = VECNORM2D(corners(:,iNode))
    IF (dist .GT. rmax) rmax = dist
  END DO

  !-- rmin: minimum distance to the origin, considering the distance to the edges, which might be minimal between corners
  rmin = HUGE(1.)
  DO iNode = 1, 4
    ! MOD(iNode,4)+1 connects 1 -> 2 (bottom), 2 -> 3 (right), 3 -> 4 (top), 4 -> 1 (left)
    iNext = MOD(iNode, 4) + 1
    point = corners(:, iNode)
    vec   = corners(:, iNext) - corners(:, iNode)
    ! Determine the closest point on the edge to the origin
    vec  = point + MIN(MAX(-DOT_PRODUCT(point, vec) / DOT_PRODUCT(vec, vec), 0.), 1.) * vec
    dist = VECNORM2D(vec)
    IF (dist .LT. rmin) rmin = dist
  END DO

  !-- Determine if the origin is inside of bounding box
  r0inside = .FALSE.
  IF ( (0. .GE. corners(1,1)) .AND. (0. .LE. corners(1,2)) .AND. (0. .GE. corners(2,1)) .AND. (0. .LE. corners(2,3)) ) THEN
    r0inside = .TRUE.
  END IF

  !-- Set rmin to zero to force the side to be classified as partially "inside", otherwise keep the smallest distance
  IF (r0inside) rmin = 0.
END IF

END SUBROUTINE GetRadialDistance2D


!===================================================================================================================================
!> Determines the minimum and maximum distance from a side's bounding box to a stadium (race track) shape.
!>
!>    corner 4 (min,max) -------- corner 3 (max,max)
!>           |                           |
!>           |                           |
!>           |                           |
!>    corner 1 (min,min) -------- corner 2 (max,min)
!>
!> The stadium is defined by a central line segment (origin +/- halfLength*dirVec) with a sweep radius.
!> The distances rmin/rmax are the closest/farthest distances from the bounding box to the central segment.
!===================================================================================================================================
SUBROUTINE GetRacetrackDistance2D(GlobalSideID, dir, origin, dirVec, halfLength, rmin, rmax)
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces       ,ONLY: GetSideBoundingBox
USE MOD_Particle_Mesh_Tools     ,ONLY: GetSideBoundingBoxTria
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: GlobalSideID, dir(3)
REAL, INTENT(IN)              :: origin(2)                  !< Center of the racetrack/stadium in surface plane coordinates
REAL, INTENT(IN)              :: dirVec(2)                  !< Normalized direction vector of the straight section
REAL, INTENT(IN)              :: halfLength                 !< Half-length of the straight section (> 0)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: rmin, rmax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode, iNext
REAL                          :: BoundingBox(1:3,1:8)
REAL                          :: corners(2,4)               !< Bounding box corners in origin-shifted 2D
REAL                          :: segA(2), segB(2)           !< Central segment endpoints in origin-shifted 2D
REAL                          :: d1, d2, d3, d4, edgeDist
LOGICAL                       :: segInside
!===================================================================================================================================
! Get bounding box
IF (TrackingMethod.EQ.TRIATRACKING) THEN
  CALL GetSideBoundingBoxTria(GlobalSideID, BoundingBox)
ELSE
  CALL GetSideBoundingBox(GlobalSideID, BoundingBox)
END IF

! Central segment endpoints in origin-shifted 2D surface coordinates
segA = -halfLength * dirVec
segB =  halfLength * dirVec

! Extract 4 bounding box corners in origin-shifted 2D surface coordinates
! Corner ordering: 1 (min,min), 2 (max,min), 3 (max,max), 4 (min,max)
corners(1,1) = MINVAL(BoundingBox(dir(2),:)) - origin(1)
corners(2,1) = MINVAL(BoundingBox(dir(3),:)) - origin(2)
corners(1,2) = MAXVAL(BoundingBox(dir(2),:)) - origin(1)
corners(2,2) = corners(2,1)
corners(1,3) = corners(1,2)
corners(2,3) = MAXVAL(BoundingBox(dir(3),:)) - origin(2)
corners(1,4) = corners(1,1)
corners(2,4) = corners(2,3)

!-- rmax: maximum distance from any corner to the central segment
rmax = 0.
DO iNode = 1, 4
  d1 = PointToSegmentDist2D(corners(:,iNode), segA, segB)
  IF (d1 .GT. rmax) rmax = d1
END DO

!-- rmin: minimum distance from any bounding box edge to the central segment
!   For non-intersecting segments, the minimum distance is attained at an endpoint.
!   We check 4 candidate distances per edge: 2 edge endpoints to central segment,
!   and 2 central segment endpoints to the edge.
rmin = HUGE(1.)
DO iNode = 1, 4
  ! MOD(iNode,4)+1 connects 1 -> 2 (bottom), 2 -> 3 (right), 3 -> 4 (top), 4 -> 1 (left)
  iNext = MOD(iNode, 4) + 1
  d1 = PointToSegmentDist2D(corners(:,iNode),  segA, segB)
  d2 = PointToSegmentDist2D(corners(:,iNext),  segA, segB)
  d3 = PointToSegmentDist2D(segA, corners(:,iNode), corners(:,iNext))
  d4 = PointToSegmentDist2D(segB, corners(:,iNode), corners(:,iNext))
  edgeDist = MIN(d1, d2, d3, d4)
  IF (edgeDist .LT. rmin) rmin = edgeDist
END DO

!-- Check if any part of the central segment lies inside the bounding box
segInside = .FALSE.

! Check if either central segment endpoint is inside the bounding box
IF ( (segA(1) .GE. corners(1,1)) .AND. (segA(1) .LE. corners(1,2)) .AND. &
     (segA(2) .GE. corners(2,1)) .AND. (segA(2) .LE. corners(2,3)) ) THEN
  segInside = .TRUE.
END IF
IF (.NOT. segInside) THEN
  IF ( (segB(1) .GE. corners(1,1)) .AND. (segB(1) .LE. corners(1,2)) .AND. &
       (segB(2) .GE. corners(2,1)) .AND. (segB(2) .LE. corners(2,3)) ) THEN
    segInside = .TRUE.
  END IF
END IF

! Check if the central segment intersects any bounding box edge
IF (.NOT. segInside) THEN
  DO iNode = 1, 4
    iNext = MOD(iNode, 4) + 1
    IF (SegmentsIntersect2D(segA, segB, corners(:,iNode), corners(:,iNext))) THEN
      segInside = .TRUE.
      EXIT
    END IF
  END DO
END IF

IF (segInside) rmin = 0.

END SUBROUTINE GetRacetrackDistance2D


!===================================================================================================================================
!> Computes the minimum distance from a 2D point to a line segment defined by endpoints segA and segB.
!> If the segment is degenerate (segA = segB), returns the distance to that point but should not happen, as this is treated with a
!> separate case with GetRadialDistance2D
!===================================================================================================================================
FUNCTION PointToSegmentDist2D(point, segA, segB) RESULT(dist)
! MODULES
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: point(2)   !< Query point
REAL, INTENT(IN)  :: segA(2)    !< Segment start point
REAL, INTENT(IN)  :: segB(2)    !< Segment end point
!-----------------------------------------------------------------------------------------------------------------------------------
! RESULT
REAL              :: dist
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: seg(2), rel(2), diff(2), relDist, segLenSq
!===================================================================================================================================
! Distance and length (squared) between half circles
seg      = segB - segA
segLenSq = DOT_PRODUCT(seg,seg)
! Vector between point and segment start
rel      = point - segA

! Degenerate segment: both endpoints coincide
IF (segLenSq .EQ. 0.) THEN
  dist = VECNORM2D(rel)
  RETURN
END IF

! Project point onto segment line, normalize to [0,1]
relDist    = DOT_PRODUCT(rel, seg) / segLenSq
! If distance is greater than 1 or smaller than 0, move the point of evaluation to the end of the respective segment
relDist    = MAX(0., MIN(1., relDist))
! Move point along the line segment to get the distance perpendicular to it
diff = point - (segA + relDist * seg)
dist = VECNORM2D(diff)

END FUNCTION PointToSegmentDist2D


!===================================================================================================================================
!> Tests whether two 2D line segments (a1,a2) and (b1,b2) have a proper (non-collinear) intersection.
!> Uses the cross-product orientation test.
!===================================================================================================================================
FUNCTION SegmentsIntersect2D(a1, a2, b1, b2) RESULT(intersect)
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: a1(2)   !< Segment A start
REAL, INTENT(IN)  :: a2(2)   !< Segment A end
REAL, INTENT(IN)  :: b1(2)   !< Segment B start
REAL, INTENT(IN)  :: b2(2)   !< Segment B end
!-----------------------------------------------------------------------------------------------------------------------------------
! RESULT
LOGICAL           :: intersect
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: da(2), db(2), d1, d2, d3, d4
!===================================================================================================================================
da = a2 - a1
db = b2 - b1

! Orientation of b1, b2 w.r.t. segment A
d1 = da(1) * (b1(2) - a1(2)) - da(2) * (b1(1) - a1(1))
d2 = da(1) * (b2(2) - a1(2)) - da(2) * (b2(1) - a1(1))

! Orientation of a1, a2 w.r.t. segment B
d3 = db(1) * (a1(2) - b1(2)) - db(2) * (a1(1) - b1(1))
d4 = db(1) * (a2(2) - b1(2)) - db(2) * (a2(1) - b1(1))

intersect = (d1 * d2 .LT. 0.) .AND. (d3 * d4 .LT. 0.)

END FUNCTION SegmentsIntersect2D


SUBROUTINE SampleSurfaceGroupProperties(SurfSideID,PartID,SpecID,SampleType,TorqueArray,ETrans,PartWeight)
!===================================================================================================================================
!> Sampling of torque and energy for surface group output
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Vars
USE MOD_Globals                   ,ONLY: abort
USE MOD_DSMC_Vars                 ,ONLY: useDSMC,PartIntEn, SpecDSMC
USE MOD_DSMC_Vars                 ,ONLY: CollisMode,DSMC
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfaceGroup
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: SurfSideID          !< Surface ID
INTEGER,INTENT(IN)       :: PartID              !< Particle ID
INTEGER,INTENT(IN)       :: SpecID              !< Particle species ID
CHARACTER(*),INTENT(IN)  :: SampleType
REAL,INTENT(IN)          :: TorqueArray(3)      !< Torque Array of impacting particle
REAL,INTENT(IN)          :: ETrans              !< Translational energy of impacting particle
REAL,INTENT(IN)          :: PartWeight          !< Particle weight (not including time step factor)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGroup
!-----------------------------------------------------------------------------------------------------------------------------------
iGroup = SurfaceGroup%SurfSide2GroupID(SurfSideID)
IF(iGroup.NE.0) THEN
  SurfaceGroup%SampState(1,iGroup) = SurfaceGroup%SampState(1,iGroup) + TorqueArray(1) * SurfaceGroup%SymmetryFactor(SurfSideID)
  SurfaceGroup%SampState(2,iGroup) = SurfaceGroup%SampState(2,iGroup) + TorqueArray(2) * SurfaceGroup%SymmetryFactor(SurfSideID)
  SurfaceGroup%SampState(3,iGroup) = SurfaceGroup%SampState(3,iGroup) + TorqueArray(3) * SurfaceGroup%SymmetryFactor(SurfSideID)
  SELECT CASE (TRIM(SampleType))
  CASE ('old')
    SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) + ETrans * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
    IF (useDSMC) THEN
      IF (CollisMode.GT.1) THEN
        IF ((Species(SpecID)%InterID.EQ.2).OR.Species(SpecID)%InterID.EQ.20) THEN
          !----  Sampling the internal (rotational) energy accommodation at walls
          SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) + PartIntEn(PartID)%ERot(1)* PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
          !----  Sampling for internal (vibrational) energy accommodation at walls
          SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) + PartIntEn(PartID)%EVib(1) * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
        END IF
        IF(DSMC%ElectronicModel.GT.0) THEN
          IF((Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized)) THEN
          !----  Sampling for internal (electronic) energy accommodation at walls
            SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) + PartIntEn(PartID)%EElec(1) * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
          END IF
        END IF
      END IF
    END IF
    ! Count the impacts to detect groups without any wall interaction during the sampling interval (only on impact)
    SurfaceGroup%Counter(iGroup) = SurfaceGroup%Counter(iGroup) + 1
  CASE ('new')
    SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) - ETrans * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
    IF (useDSMC) THEN
      IF (CollisMode.GT.1) THEN
        IF ((Species(SpecID)%InterID.EQ.2).OR.Species(SpecID)%InterID.EQ.20) THEN
          !----  Sampling the internal (rotational) energy accommodation at walls
          SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) - PartIntEn(PartID)%ERot(1) * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
          !----  Sampling for internal (vibrational) energy accommodation at walls
          SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) - PartIntEn(PartID)%EVib(1) * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
        END IF
        IF(DSMC%ElectronicModel.GT.0) THEN
          IF((Species(SpecID)%InterID.NE.4).AND.(.NOT.SpecDSMC(SpecID)%FullyIonized)) THEN
          !----  Sampling for internal (electronic) energy accommodation at walls
            SurfaceGroup%SampState(4,iGroup) = SurfaceGroup%SampState(4,iGroup) - PartIntEn(PartID)%EElec(1) * PartWeight * SurfaceGroup%SymmetryFactor(SurfSideID)
          END IF
        END IF
      END IF
    END IF
  CASE DEFAULT
    CALL abort(__STAMP__,'ERROR in CalcWallSample: wrong SampleType specified. Possible types -> ( old , new )')
  END SELECT
END IF

END SUBROUTINE SampleSurfaceGroupProperties

END MODULE MOD_Particle_Boundary_Tools
