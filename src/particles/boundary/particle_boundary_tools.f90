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
PUBLIC :: MarkAuxBCElems
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
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallState,CalcSurfaceImpact
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
  IF (VarTimeStep%UseVariableTimeStep) THEN
    SampWallState(SAMPWALL_NVARS+nSpecies+1,SubP,SubQ,SurfSideID) = SampWallState(SAMPWALL_NVARS+nSpecies+1,SubP,SubQ,SurfSideID) &
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


SUBROUTINE StoreBoundaryParticleProperties(iPart,SpecID,PartPos,PartTrajectory,SurfaceNormal,iBC,mode,usevMPF_optIN)
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
INTEGER,INTENT(IN) :: iBC  ! Part-BoundaryX on which the impact occurs
INTEGER,INTENT(IN) :: mode ! 1: particle impacts on BC (species is stored as positive value)
                            ! 2: particles is emitted from the BC into the simulation domain (species is stored as negative value)
LOGICAL,INTENT(IN),OPTIONAL :: usevMPF_optIN ! For setting MPF for cases when PartMPF(iPart) might not yet be set during emission
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: MPF
INTEGER              :: dims(2)
! Temporary arrays
REAL, ALLOCATABLE    :: PartStateBoundary_tmp(:,:) ! (1:11,1:NParts) 1st index: x,y,z,vx,vy,vz,SpecID,Ekin,MPF,time,impact angle,iBC
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
    ALLOCATE(PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary_tmp temporary array!')
    ! Save old data
    PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2)) = PartStateBoundary(1:nVarPartStateBoundary,1:dims(2))

    ! Re-allocate PartStateBoundary to twice the size
    DEALLOCATE(PartStateBoundary)
    ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:2*dims(2)), STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
          __STAMP__&
          ,'ERROR in particle_boundary_tools.f90: Cannot allocate PartStateBoundary array!')
    PartStateBoundary(1:nVarPartStateBoundary,        1:  dims(2)) = PartStateBoundary_tmp(1:nVarPartStateBoundary,1:dims(2))
    PartStateBoundary(1:nVarPartStateBoundary,dims(2)+1:2*dims(2)) = 0.

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
  PartStateBoundary(11 ,iMax) = REAL(iBC)
END ASSOCIATE

END SUBROUTINE StoreBoundaryParticleProperties


SUBROUTINE MarkAuxBCElems()
!===================================================================================================================================
! check if auxBCs are inside BoundingBox of Elems
! -- plane: use plane equation f=a1*x+a2*y+a3*z+a4=0 and insert corresponding intervals of box -> fmin and fmax
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemHasAuxBCs
USE MOD_Particle_Mesh_Vars     ,ONLY: BoundsOfElem_Shared
USE MOD_Particle_Boundary_Vars ,ONLY: nAuxBCs,AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone
#if USE_MPI
USE MOD_MPI_Shared_Vars
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,iAuxBC,icoord,dir(3),positiontype,positiontype_tmp
REAL                     :: r_vec(3),n_vec(3),fmin,fmax,radius,BoundsBC(1:2,1:3)
REAL                     :: lmin,lmax,deltamin,deltamax,origin(2),halfangle
LOGICAL                  :: cartesian, backwards
!===================================================================================================================================

ALLOCATE(ElemHasAuxBCs(1:PP_nElems , 1:nAuxBCs))
ElemHasAuxBCs=.FALSE.

DO iAuxBC=1,nAuxBCs
  SELECT CASE (TRIM(AuxBCType(iAuxBC)))
  CASE ('plane')
    r_vec=AuxBC_plane(AuxBCMap(iAuxBC))%r_vec
    n_vec=AuxBC_plane(AuxBCMap(iAuxBC))%n_vec
    radius=AuxBC_plane(AuxBCMap(iAuxBC))%radius
    ! loop over all  elements
    DO iElem=1,PP_nElems
      ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
        fmin=-DOT_PRODUCT(r_vec,n_vec)
        fmax=fmin
        DO icoord=1,3
          IF (n_vec(icoord).GE.0) THEN
            fmin = fmin + n_vec(icoord)*Bounds(1,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(2,icoord)
          ELSE
            fmin = fmin + n_vec(icoord)*Bounds(2,icoord)
            fmax = fmax + n_vec(icoord)*Bounds(1,icoord)
          END IF
        END DO
        IF ((fmin.LE.0 .AND. fmax.GT.0).OR.(fmin.LT.0 .AND. fmax.GE.0)) THEN !plane intersects the box!
          !radius check needs to be implemented (compute intersection polygon and minimum radii): would sort out further elements!!!
          !quick, conservative solution: calculate bounding box of disc in space and compare with bb of element
          ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
          IF (radius .LT. 0.5*HUGE(radius)) THEN !huge was default
            BoundsBC(1,1:3) = r_vec - radius * SQRT(1.-(n_vec*n_vec))
            BoundsBC(2,1:3) = r_vec + radius * SQRT(1.-(n_vec*n_vec))
            DO icoord=1,3
              IF ( BoundsBC(2,icoord).LT.Bounds(1,icoord) .OR. BoundsBC(1,icoord).GT.Bounds(2,icoord) ) THEN
                ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
                EXIT
              END IF
            END DO
          END IF
        ELSE IF ((fmin.LT.0 .AND. fmax.LT.0).OR.(fmin.GT.0 .AND. fmax.GT.0)) THEN !plane does not intersect the box!
          ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
        ELSE !e.g. if elem has zero volume...
          CALL abort(&
            __STAMP__&
            ,'Error in MarkAuxBCElems for AuxBC:',iAuxBC)
        END IF
      END ASSOCIATE
    END DO
  CASE ('cylinder','cone')
    IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
      r_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cylinder(AuxBCMap(iAuxBC))%axis
      radius=AuxBC_cylinder(AuxBCMap(iAuxBC))%radius
      lmin=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cylinder(AuxBCMap(iAuxBC))%lmax
    ELSE !cone
      r_vec=AuxBC_cone(AuxBCMap(iAuxBC))%r_vec
      n_vec=AuxBC_cone(AuxBCMap(iAuxBC))%axis
      halfangle=AuxBC_cone(AuxBCMap(iAuxBC))%halfangle
      lmin=AuxBC_cone(AuxBCMap(iAuxBC))%lmin
      lmax=AuxBC_cone(AuxBCMap(iAuxBC))%lmax
    END IF
    cartesian=.TRUE.
    backwards=.FALSE.
    IF (ABS(n_vec(1)).EQ.1.) THEN
      dir(1)=1
      dir(2)=2
      dir(3)=3
      IF (n_vec(1).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(2)).EQ.1.) THEN
      dir(1)=2
      dir(2)=3
      dir(3)=1
      IF (n_vec(2).LT.0.) backwards=.TRUE.
    ELSE IF (ABS(n_vec(3)).EQ.1.) THEN
      dir(1)=3
      dir(2)=1
      dir(3)=2
      IF (n_vec(3).LT.0.) backwards=.TRUE.
    ELSE
      cartesian=.FALSE.
      SWRITE(*,*) 'WARNING in MarkAuxBCElems: all Elems are set to ElemHasAuxBCs=.TRUE. for AuxBC:',iAuxBC
      ElemHasAuxBCs(:,iAuxBC)=.TRUE. !actual intersection with box check to-be implemented!!!
    END IF
    IF (cartesian) THEN
      IF (backwards) THEN
        deltamin = -lmax
        deltamax = -lmin
      ELSE
        deltamin = lmin
        deltamax = lmax
      END IF
      origin(1) = r_vec(dir(2))
      origin(2) = r_vec(dir(3))
      ! loop over all  elements
      DO iElem=1,PP_nElems
        ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
          ! check for lmin and lmax
          IF ( r_vec(dir(1))+deltamax.LT.Bounds(1,dir(1)) .OR. r_vec(dir(1))+deltamin.GT.Bounds(2,dir(1)) ) THEN
            ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
          ELSE !between lmin and lmax
            IF (TRIM(AuxBCType(iAuxBC)).EQ.'cylinder') THEN
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
            ELSE !cone
              !local minimum radius
              IF (backwards) THEN
                radius = MAX(-Bounds(2,dir(1))+r_vec(dir(1)),lmin)*TAN(halfangle)
              ELSE
                radius = MAX(Bounds(1,dir(1))-r_vec(dir(1)),lmin)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype_tmp)
              !local maximum radius
              IF (backwards) THEN
                radius = MIN(-Bounds(1,dir(1))+r_vec(dir(1)),lmax)*TAN(halfangle)
              ELSE
                radius = MIN(Bounds(2,dir(1))-r_vec(dir(1)),lmax)*TAN(halfangle)
              END IF
              CALL CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
              !if both are type 0 or both are type 1 than the "total" type is not 2:
              IF ( .NOT.(positiontype_tmp.EQ.0 .AND. positiontype.EQ.0) &
                .AND. .NOT.(positiontype_tmp.EQ.1 .AND. positiontype.EQ.1) ) THEN
                positiontype=2
              END IF
            END IF
            IF (positiontype.EQ.2) THEN
              ElemHasAuxBCs(iElem,iAuxBC)=.TRUE.
            ELSE
              ElemHasAuxBCs(iElem,iAuxBC)=.FALSE.
            END IF
          END IF !check for lmin and lmax
        END ASSOCIATE
      END DO !iElem
    END IF !cartesian
  CASE('parabol')
    ElemHasAuxBCs(:,iAuxBC)=.TRUE. ! to be implemented!!!
  CASE DEFAULT
    SWRITE(*,*) ' AuxBC does not exist: ', TRIM(AuxBCType(iAuxBC))
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END SELECT
END DO

END SUBROUTINE MarkAuxBCElems


SUBROUTINE CheckBoundsWithCartRadius(Bounds,dir,origin,radius,positiontype)
!===================================================================================================================================
! checks how a cartesian bb is located with regard to a radius with cartesian axis (dir is cartesian axis and origin in orth. dirs)
!- positiontype=0 : complete bb is inside of radius
!- positiontype=1 : complete bb is outside of radius
!- positiontype=2 : bb is partly inside of radius
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN)           :: Bounds(1:2,1:3), origin(2), radius
INTEGER,INTENT(IN)        :: dir(3)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)       :: positiontype
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iDir1, iDir2, iDir3, iPoint
REAL                      :: BoundingBox(1:3,1:8), point(2), pointRadius
LOGICAL                   :: done, insideBound
!===================================================================================================================================
!-- convert minmax-values to bb-points
DO iDir1=0,1
  DO iDir2=0,1
      DO iDir3=0,1
        BoundingBox(1,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir1+1,1)
        BoundingBox(2,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir2+1,2)
        BoundingBox(3,iDir1*4 + iDir2*2 + iDir3+1) = Bounds(iDir3+1,3)
      END DO
  END DO
END DO

!-- check where the points are located relative to radius
done=.FALSE.
insideBound=.FALSE. ! Initialize
DO iDir1=0,1
  IF(done) EXIT
  DO iDir2=0,1
    IF(done) EXIT
    DO iDir3=0,1
      !-- coords orth. to axis of point:
      iPoint=iDir1*4 + iDir2*2 + iDir3+1
      point(1) = BoundingBox(dir(2),iPoint)-origin(1)
      point(2) = BoundingBox(dir(3),iPoint)-origin(2)
      pointRadius = SQRT( (point(1))**2+(point(2))**2 )
      IF (iPoint.EQ.1) THEN
        IF (pointRadius.LE.radius) THEN
          insideBound=.TRUE.
        ELSE !outside
          insideBound=.FALSE.
        END IF !in-/outside?
      ELSE !iPoint.GT.1: type must be 2 if state of point if different from last point
        IF (pointRadius.LE.radius) THEN
          IF (.NOT.insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        ELSE !outside
          IF (insideBound) THEN !different from last point
            positiontype=2
            done=.TRUE.
            EXIT
          END IF
        END IF !in-/outside?
      END IF !iPoint.EQ.1
    END DO !iDir3
  END DO !iDir2
END DO !iDir1
IF (.NOT.done) THEN
  IF (insideBound) THEN
    positiontype=0
  ELSE
    ! all points are outside of radius, but when radius is smaller than box, it can intersect it:
    IF ( origin(1) + radius .GE. Bounds(1,dir(2)) .AND. &
         origin(1) - radius .LE. Bounds(2,dir(2)) .AND. &
         origin(2) + radius .GE. Bounds(1,dir(3)) .AND. &
         origin(2) - radius .LE. Bounds(2,dir(3)) ) THEN !circle completely or partly inside box
      positiontype=2
    ELSE !points are really outside
      positiontype=1
    END IF
  END IF
END IF

END SUBROUTINE CheckBoundsWithCartRadius


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

END SUBROUTINE GetRadialDistance2D

END MODULE MOD_Particle_Boundary_Tools
