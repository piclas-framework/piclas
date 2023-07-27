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
PUBLIC :: GetRadialDistance2D
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcWallSample(PartID,SurfSideID,SampleType,SurfaceNormal_opt)
!===================================================================================================================================
!> Sample the energy of particles before and after a wall interaction for the determination of macroscopic properties such as heat
!> flux and force per area
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars
USE MOD_Globals                   ,ONLY: abort,DOTPRODUCT
USE MOD_DSMC_Vars                 ,ONLY: SpecDSMC,useDSMC,PartStateIntEn,RadialWeighting
USE MOD_DSMC_Vars                 ,ONLY: CollisMode,DSMC,AmbipolElecVelo
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallState,CalcSurfaceImpact,SWIVarTimeStep
USE MOD_part_tools                ,ONLY: GetParticleWeight
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: PartID,SurfSideID
CHARACTER(*),INTENT(IN)            :: SampleType
REAL,INTENT(IN),OPTIONAL           :: SurfaceNormal_opt(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: ETrans, ETransAmbi, ERot, EVib, EElec, MomArray(1:3), MassIC, MPF
INTEGER         :: ETransID, ERotID, EVibID, EElecID, SpecID, SubP, SubQ
!===================================================================================================================================
MomArray(:)=0.
EVib = 0.
ERot = 0.
EElec = 0.

SubP = TrackInfo%p
SubQ = TrackInfo%q

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
  SampWallState(SAMPWALL_NVARS+SpecID,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_NVARS+SpecID,SubP,SubQ,SurfSideID) + 1
  ! Sampling of species-specific impact energies and angles
  IF(CalcSurfaceImpact) THEN
    IF (useDSMC) THEN
      IF (CollisMode.GT.1) THEN
        EVib = PartStateIntEn(1,PartID)
        ERot = PartStateIntEn(2,PartID)
        IF(DSMC%ElectronicModel.GT.0) THEN
          EElec = PartStateIntEn(3,PartID)
        END IF
      END IF
    END IF
    CALL SampleImpactProperties(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,EElec,TrackInfo%PartTrajectory,SurfaceNormal_opt)
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
        CALL SampleImpactProperties(SurfSideID,DSMC%AmbiDiffElecSpec,MPF,ETransAmbi,0.,0.,0.,TrackInfo%PartTrajectory,SurfaceNormal_opt)
      END IF
    END IF
  END IF
  ! Sample the time step for the correct determination of the heat flux
  IF (UseVarTimeStep) THEN
    SampWallState(SWIVarTimeStep,SubP,SubQ,SurfSideID) = SampWallState(SWIVarTimeStep,SubP,SubQ,SurfSideID) &
                                                              + PartTimeStep(PartID)
  ELSE IF(VarTimeStep%UseSpeciesSpecific) THEN
    SampWallState(SWIVarTimeStep,SubP,SubQ,SurfSideID) = SampWallState(SWIVarTimeStep,SubP,SubQ,SurfSideID) &
                                                              + Species(SpecID)%TimeStepFactor
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
  CALL abort(__STAMP__,'ERROR in CalcWallSample: wrong SampleType specified. Possible types -> ( old , new )')
END SELECT
!----  Sampling force at walls (correct sign is set above)
SampWallState(SAMPWALL_DELTA_MOMENTUMX,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMX,SubP,SubQ,SurfSideID) + MomArray(1)
SampWallState(SAMPWALL_DELTA_MOMENTUMY,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMY,SubP,SubQ,SurfSideID) + MomArray(2)
SampWallState(SAMPWALL_DELTA_MOMENTUMZ,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_DELTA_MOMENTUMZ,SubP,SubQ,SurfSideID) + MomArray(3)
!----  Sampling the energy (translation) accommodation at walls
SampWallState(ETransID ,SubP,SubQ,SurfSideID) = SampWallState(ETransID ,SubP,SubQ,SurfSideID) + ETrans * MPF
IF (useDSMC) THEN
  IF (CollisMode.GT.1) THEN
    IF ((SpecDSMC(SpecID)%InterID.EQ.2).OR.SpecDSMC(SpecID)%InterID.EQ.20) THEN
      !----  Sampling the internal (rotational) energy accommodation at walls
      SampWallState(ERotID ,SubP,SubQ,SurfSideID) = SampWallState(ERotID ,SubP,SubQ,SurfSideID) + PartStateIntEn(2,PartID) * MPF
      !----  Sampling for internal (vibrational) energy accommodation at walls
      SampWallState(EVibID ,SubP,SubQ,SurfSideID) = SampWallState(EVibID ,SubP,SubQ,SurfSideID) + PartStateIntEn(1,PartID) * MPF
    END IF
    IF(DSMC%ElectronicModel.GT.0) THEN
      !----  Sampling for internal (electronic) energy accommodation at walls
      SampWallState(EElecID ,SubP,SubQ,SurfSideID) = SampWallState(EElecID ,SubP,SubQ,SurfSideID) + PartStateIntEn(3,PartID) * MPF
    END IF
  END IF
END IF

END SUBROUTINE CalcWallSample


SUBROUTINE SampleImpactProperties(SurfSideID,SpecID,MPF,ETrans,EVib,ERot,EElec,PartTrajectory,SurfaceNormal)
!===================================================================================================================================
!> Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
!>
!===================================================================================================================================
!USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
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


SUBROUTINE StoreBoundaryParticleProperties(iPart,SpecID,PartPos,PartTrajectory,SurfaceNormal,iPartBound,mode,MPF_optIN,Velo_optIN)
!----------------------------------------------------------------------------------------------------------------------------------!
! Save particle position, velocity and species to PartDataBoundary container for writing to .h5 later
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                ,ONLY: abort
USE MOD_Particle_Vars          ,ONLY: usevMPF,PartMPF,Species,PartState
USE MOD_Particle_Boundary_Vars ,ONLY: PartStateBoundary,PartStateBoundaryVecLength,nVarPartStateBoundary
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
INTEGER,INTENT(IN) :: iPartBound  ! Part-BoundaryX on which the impact occurs
INTEGER,INTENT(IN) :: mode ! 1: particle impacts on BC (species is stored as positive value)
                           ! 2: particles is emitted from the BC into the simulation domain (species is stored as negative value)
REAL,INTENT(IN),OPTIONAL :: MPF_optIN ! Supply the MPF in special cases
REAL,INTENT(IN),OPTIONAL :: Velo_optIN(1:3) ! Supply the velocity in special cases
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF
INTEGER              :: dims(2)
! Temporary arrays
REAL, ALLOCATABLE    :: PartStateBoundary_tmp(:,:) ! (1:11,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle,iPartBound
!                                                  !                 2nd index: 1 to number of boundary-crossed particles
INTEGER              :: ALLOCSTAT
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

  ! Check if array maximum is reached.
  ! If this happens, re-allocate the arrays and increase their size (every time this barrier is reached, double the size)
  IF(iMax.GT.dims(2))THEN

    ! --- PartStateBoundary ---
    ALLOCATE(PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2)) = PartStateBoundary(1:nVarPartStateBoundary,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:2*dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:nVarPartStateBoundary,        1:  dims(2)) = PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2))
    PartStateBoundary(1:nVarPartStateBoundary,dims(2)+1:2*dims(2)) = 0.

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
  PartStateBoundary(8  ,iMax) = MPF
  PartStateBoundary(9  ,iMax) = time
  PartStateBoundary(10 ,iMax) = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(PartTrajectory,SurfaceNormal))))
  PartStateBoundary(11 ,iMax) = REAL(iPartBound)
END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties


SUBROUTINE GetRadialDistance2D(GlobalSideID,dir,origin,rmin,rmax)
!===================================================================================================================================
! Determines the radial distance to a given origin on a surface
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces       ,ONLY: GetSideBoundingBox
USE MOD_Particle_Mesh_Tools     ,ONLY: GetSideBoundingBoxTria
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Vars           ,ONLY: Symmetry
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: GlobalSideID, dir(3)
REAL, INTENT(IN)              :: origin(2)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: rmin,rmax
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode
REAL                          :: BoundingBox(1:3,1:8), point(2), vec(2)
REAL                          :: Vector1(3),Vector2(3),Vector3(3),xyzNod(3),corner(3),VecBoundingBox(3),radiusCorner(2,4)
LOGICAL                       :: r0inside
!===================================================================================================================================
! Determine which cells are inside/outside/partially inside the defined region
IF (TrackingMethod.EQ.TRIATRACKING) THEN
  CALL GetSideBoundingBoxTria(GlobalSideID,BoundingBox)
ELSE
  CALL GetSideBoundingBox(GlobalSideID,BoundingBox)
END IF
IF(Symmetry%Axisymmetric) THEN
  rmin = BoundingBox(2,1)
  rmax = BoundingBox(2,3)
ELSE
  r0inside=.FALSE.
  Vector1(:)=0.
  Vector2(:)=0.
  Vector3(:)=0.
  xyzNod(1)=MINVAL(BoundingBox(1,:))
  xyzNod(2)=MINVAL(BoundingBox(2,:))
  xyzNod(3)=MINVAL(BoundingBox(3,:))
  VecBoundingBox(1) = MAXVAL(BoundingBox(1,:)) -MINVAL(BoundingBox(1,:))
  VecBoundingBox(2) = MAXVAL(BoundingBox(2,:)) -MINVAL(BoundingBox(2,:))
  VecBoundingBox(3) = MAXVAL(BoundingBox(3,:)) -MINVAL(BoundingBox(3,:))
  Vector1(dir(2)) = VecBoundingBox(dir(2))
  Vector2(dir(2)) = VecBoundingBox(dir(2))
  Vector2(dir(3)) = VecBoundingBox(dir(3))
  Vector3(dir(3)) = VecBoundingBox(dir(3))
  !-- determine rmax (and corners)
  DO iNode=1,4
    SELECT CASE(iNode)
    CASE(1)
      corner = xyzNod
    CASE(2)
      corner = xyzNod + Vector1
    CASE(3)
      corner = xyzNod + Vector2
    CASE(4)
      corner = xyzNod + Vector3
    END SELECT
    corner(dir(2)) = corner(dir(2)) - origin(1)
    corner(dir(3)) = corner(dir(3)) - origin(2)
    radiusCorner(1,iNode)=SQRT(corner(dir(2))**2+corner(dir(3))**2)
  END DO !iNode
  rmax=MAXVAL(radiusCorner(1,1:4))
  !-- determine rmin
  DO iNode=1,4
    SELECT CASE(iNode)
    CASE(1)
      point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
      vec=(/Vector1(dir(2)),Vector1(dir(3))/)
    CASE(2)
      point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
      vec=(/Vector3(dir(2)),Vector3(dir(3))/)
    CASE(3)
      point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
      vec=(/-Vector1(dir(2)),-Vector1(dir(3))/)
    CASE(4)
      point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
      vec=(/-Vector3(dir(2)),-Vector3(dir(3))/)
    END SELECT
    vec=point + MIN(MAX(-DOT_PRODUCT(point,vec)/DOT_PRODUCT(vec,vec),0.),1.)*vec
    radiusCorner(2,iNode)=SQRT(DOT_PRODUCT(vec,vec)) !rmin
  END DO !iNode
  !-- determine if r0 is inside of bounding box
  IF ((origin(1) .GE. MINVAL(BoundingBox(dir(2),:))) .AND. &
      (origin(1) .LE. MAXVAL(BoundingBox(dir(2),:))) .AND. &
      (origin(2) .GE. MINVAL(BoundingBox(dir(3),:))) .AND. &
      (origin(2) .LE. MAXVAL(BoundingBox(dir(3),:))) ) THEN
      r0inside = .TRUE.
  END IF
  IF (r0inside) THEN
    rmin = 0.
  ELSE
    rmin=MINVAL(radiusCorner(2,1:4))
  END IF
END IF

END SUBROUTINE GetRadialDistance2D

END MODULE MOD_Particle_Boundary_Tools
