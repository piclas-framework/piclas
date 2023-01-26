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

MODULE MOD_Particle_Boundary_Condition
!===================================================================================================================================
!> Determines how particles interact with a given boundary condition. This routine is used by MOD_Part_Tools, hence, it cannot be
!> used here due to circular definitions!
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
PUBLIC :: GetBoundaryInteraction
PUBLIC :: GetBoundaryInteractionAuxBC
!PUBLIC :: PartSwitchElement
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(iPart,SideID,flip,ElemID,crossedBC,TriNum)
!===================================================================================================================================
!> Determines the post boundary state of a particle that interacts with a boundary condition
!> * Open: Particle is removed
!> * Reflective: Further treatment as part of the surface modelling (`src/particles/surfacemodel/`)
!> * Periodic (+ Rotational periodic): Further treatment in the specific routines
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals                  ,ONLY: abort
USE MOD_Mesh_Tools               ,ONLY: GetCNSideID
USE MOD_Part_Operations          ,ONLY: RemoveParticle
USE MOD_Particle_Surfaces        ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars            ,ONLY: PartSpecies,PDM,PEM
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod, TrackInfo, CountNbrOfLostParts, NbrOfLostParticles
USE MOD_Part_Tools               ,ONLY: StoreLostParticleProperties
USE MOD_Dielectric_vars          ,ONLY: DoDielectric,isDielectricElem
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Boundary_Vars   ,ONLY: PartBound,DoBoundaryParticleOutputHDF5
USE MOD_Particle_Surfaces_vars   ,ONLY: SideNormVec,SideType
USE MOD_SurfaceModel             ,ONLY: SurfaceModel,PerfectReflection
USE MOD_Particle_Vars            ,ONLY: LastPartPos
USE MOD_Particle_Boundary_Tools  ,ONLY: StoreBoundaryParticleProperties
#ifdef CODE_ANALYZE
USE MOD_Globals                  ,ONLY: myRank,UNIT_stdout
USE MOD_Mesh_Vars                ,ONLY: NGeo
USE MOD_Mesh_Tools               ,ONLY: GetCNElemID
USE MOD_Particle_Surfaces_Vars   ,ONLY: BezierControlPoints3D
USE MOD_Particle_Mesh_Vars       ,ONLY: ElemBaryNGeo_Shared
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip
INTEGER,INTENT(IN),OPTIONAL          :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID ! Global element index
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: CNSideID
REAL                                 :: n_loc(1:3)
#ifdef CODE_ANALYZE
REAL                                 :: v1(3),v2(3)
#endif /* CODE_ANALYZE */
!===================================================================================================================================

crossedBC    =.FALSE.

! Calculate normal vector
SELECT CASE(TrackingMethod)
CASE(REFMAPPING,TRACING)
  ! set BCSideID for normal vector calculation call with (curvi-)linear side description
  CNSideID = GetCNSideID(SideID)

  SELECT CASE(SideType(CNSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc = SideNormVec(1:3,CNSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=TrackInfo%xi,eta=TrackInfo%eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(  nVec=n_loc,xi=TrackInfo%xi,eta=TrackInfo%eta,SideID=SideID)
  END SELECT

  IF(flip.NE.0) n_loc=-n_loc

#ifdef CODE_ANALYZE
  ! check if normal vector points outwards
  v1 = 0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)  &
           + BezierControlPoints3D(:,NGeo,0   ,SideID)  &
           + BezierControlPoints3D(:,0   ,NGeo,SideID)  &
           + BezierControlPoints3D(:,NGeo,NGeo,SideID))
  v2 = v1  - ElemBaryNGeo(:,GetCNElemID(ElemID))

  IF (DOT_PRODUCT(v2,n_loc).LT.0) THEN
    IPWRITE(UNIT_stdout,*) 'Obtained wrong side orientation from flip. flip:',flip,'PartID:',iPart
    IPWRITE(UNIT_stdout,*) 'n_loc (flip)', n_loc,'n_loc (estimated):',v2
    CALL ABORT(__STAMP__,'SideID',SideID)
  END IF
#endif /* CODE_ANALYZE */

  ! Inserted particles are "pushed" inside the domain and registered as passing through the BC side. If they are very close to the
  ! boundary (first if) than the normal vector is compared with the trajectory. If the particle is entering the domain from outside
  ! it was inserted during surface flux and this routine shall not performed.
  ! Comparing the normal vector with the particle trajectory, if the particle trajectory is pointing inside the domain
  IF(DOT_PRODUCT(n_loc,TrackInfo%PartTrajectory).LE.0.) RETURN

CASE(TRIATRACKING)
  CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
END SELECT
! required for refmapping and tracing, optional for triatracking
crossedBC=.TRUE.

IF(.NOT.ALLOCATED(PartBound%MapToPartBC)) CALL abort(__STAMP__,' ERROR: PartBound not allocated!.')

ASSOCIATE( iBC => PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)) )
  ! Surface particle output to .h5
  IF(DoBoundaryParticleOutputHDF5.AND.PartBound%BoundaryParticleOutputHDF5(iBC))THEN
    CALL StoreBoundaryParticleProperties(iPart,PartSpecies(iPart), &
            LastPartPos(1:3,iPart)+TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha,TrackInfo%PartTrajectory(1:3),n_loc,iBC=iBC,mode=1)
  END IF

  ! Select the corresponding boundary condition and calculate particle treatment
  SELECT CASE(PartBound%TargetBoundCond(iBC))
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(1) ! PartBound%OpenBC
  !----------------------------------------------------------------------------------------------------------------------------------
    CALL RemoveParticle(iPart,BCID=iBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(2) ! PartBound%ReflectiveBC
  !-----------------------------------------------------------------------------------------------------------------------------------
  ! Decide which interaction (specular/diffuse reflection, species swap, SEE)
    CALL SurfaceModel(iPart,SideID,ElemID,n_loc)
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(3) ! PartBound%PeriodicBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL PeriodicBC(iPart,SideID,ElemID)
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(4) ! PartBound%SimpleAnodeBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)')
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(5) ! PartBound%SimpleCathodeBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL abort(__STAMP__,' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)')
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(6) ! PartBound%rot_periodic
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL RotPeriodicBC(iPart,SideID,ElemID)
    ! Sanity check: During rotational periodic movement, particles may enter a dielectric. Unfortunately, they must be deleted
    IF(DoDielectric)THEN
      IF(PDM%ParticleInside(iPart).AND.(isDielectricElem(PEM%LocalElemID(iPart))))THEN
        IF(CountNbrOfLostParts)THEN
          CALL StoreLostParticleProperties(iPart,ElemID)
          NbrOfLostParticles=NbrOfLostParticles+1
        END IF ! CountNbrOfLostParts
        CALL RemoveParticle(iPart,BCID=iBC)
      END IF ! PDM%ParticleInside(iPart).AND.(isDielectricElem(PEM%LocalElemID(iPart)))
    END IF ! DoDdielectric
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(10,11) ! PartBound%SymmetryBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL PerfectReflection(iPart,SideID,n_loc,opt_Symmetry=.TRUE.)
  CASE DEFAULT
    CALL abort(__STAMP__,' ERROR: PartBound not associated!. (unknown case)')
END SELECT !PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)
END ASSOCIATE

END SUBROUTINE GetBoundaryInteraction


SUBROUTINE GetBoundaryInteractionAuxBC(iPart,AuxBCIdx,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with an auxBC
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals                ,ONLY: abort,UNITVECTOR
USE MOD_SurfaceModel           ,ONLY: PerfectReflection, DiffuseReflection, SpeciesSwap
USE MOD_Particle_Vars          ,ONLY: PDM
USE MOD_Particle_Boundary_Vars ,ONLY: PartAuxBC
USE MOD_Particle_Boundary_Vars ,ONLY: AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_part_operations        ,ONLY: RemoveParticle
USE MOD_Particle_Tracking_Vars ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,AuxBCIdx
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: RanNum
REAL                                 :: n_loc(1:3)
REAL                                 :: intersec(3),r_vec(3),axis(3),cos2inv!,v_2(1:3),v_aux(1:3)
!===================================================================================================================================

crossedBC    =.FALSE.
SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
CASE ('plane')
  n_loc = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
CASE ('cylinder')
  intersec = LastPartPos(1:3,iPart) + TrackInfo%alpha*TrackInfo%PartTrajectory
  r_vec = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
  IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('cone')
  intersec = LastPartPos(1:3,iPart) + TrackInfo%alpha*TrackInfo%PartTrajectory
  r_vec = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
  cos2inv = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
  IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('parabol')
  intersec = LastPartPos(1:3,iPart) + TrackInfo%alpha*TrackInfo%PartTrajectory
  r_vec = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
  IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE DEFAULT
  CALL abort(__STAMP__,'AuxBC does not exist')
END SELECT
IF(DOT_PRODUCT(n_loc,TrackInfo%PartTrajectory).LT.0.)  THEN
  crossedBC=.FALSE.
  !RETURN
  CALL abort(__STAMP__,'Error in GetBoundaryInteractionAuxBC: Particle coming from outside!')
ELSE IF(DOT_PRODUCT(n_loc,TrackInfo%PartTrajectory).GT.0.)  THEN
  crossedBC=.TRUE.
ELSE
  CALL abort(__STAMP__,'Error in GetBoundaryInteractionAuxBC: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
END IF
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartAuxBC%TargetBoundCond(AuxBCIdx))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartAuxBC%OpenBC
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL RemoveParticle(iPart)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartAuxBC%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  !---- swap species?
!print*,'*********************'
!print*,AuxBCIdx
!print*,iPart,alpha,PartState(4:6,iPart)
!print*,iPart,alpha,LastPartPos(1:3,iPart),PartState(4:6,iPart)
  IF (PartAuxBC%NbrOfSpeciesSwaps(AuxBCIdx).gt.0) THEN
#ifndef IMPA
    CALL SpeciesSwap(PartID=iPart,SideID=-1,AuxBCIdx=AuxBCIdx)
#else
    CALL SpeciesSwap(PartID=iPart,SideID=-1)
#endif /*NOT IMPA*/
  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! particle did not Swap to species 0 !deleted particle -> particle swaped to species 0
    ! simple reflection (previously used wall interaction model, maxwellian scattering)
      CALL RANDOM_NUMBER(RanNum)
      IF(RanNum.GE.PartAuxBC%MomentumACC(AuxBCIdx)) THEN
        ! perfectly reflecting, specular re-emission
        CALL PerfectReflection(PartID=iPart,SideID=-1,n_loc=n_loc,AuxBCIdx=AuxBCIdx)
      ELSE
        CALL DiffuseReflection(PartID=iPart,SideID=-1,n_loc=n_loc,AuxBCIdx=AuxBCIdx)
      END IF
  END IF
!print*,iPart,alpha,LastPartPos(1:3,iPart),PartState(1:3,iPart)
!print*,iPart,alpha,PartState(4:6,iPart)
!print*,'*********************'
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
  CALL abort(__STAMP__,' ERROR: AuxBC bound not associated!. (unknown case)',999,999.)
END SELECT

! compiler warnings
RETURN
WRITE(*,*) 'AuxBCIdx', AuxBCIdx

END SUBROUTINE GetBoundaryInteractionAuxBC


SUBROUTINE PeriodicBC(PartID,SideID,ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY: BoundaryType
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
#if defined(ROS) || defined(IMPA)
USE MOD_Particle_Vars          ,ONLY: PEM
#endif
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
#if defined(IMPA)
USE MOD_TimeDisc_Vars          ,ONLY: ESDIRK_a,ERK_a
#endif /*IMPA */
#if defined(ROS)
USE MOD_TimeDisc_Vars          ,ONLY: RK_A
#endif /*ROS */
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: PartID, SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: PVID
!INTEGER                              :: moved(2),locSideID
!===================================================================================================================================

PVID = BoundaryType(SideInfo_Shared(SIDE_BCID,SideID),BC_ALPHA)

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' ParticlePosition: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' LastPartPos:      ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! set last particle position on face
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha
! perform the periodic movement
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + SIGN( GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
! update particle positon after periodic BC
PartState(1:3,PartID) = LastPartPos(1:3,PartID) + (TrackInfo%lengthPartTrajectory-TrackInfo%alpha)*TrackInfo%PartTrajectory
TrackInfo%lengthPartTrajectory  = TrackInfo%lengthPartTrajectory - TrackInfo%alpha

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' ParticlePosition-pp: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' LastPartPo-pp:       ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

#if defined(ROS) || defined(IMPA)
PEM%PeriodicMoved(PartID)=.TRUE.
#endif

! refmapping and tracing
! move particle from old element to new element
ElemID = SideInfo_Shared(SIDE_NBELEMID,SideID)

END SUBROUTINE PeriodicBC


SUBROUTINE RotPeriodicBC(PartID,SideID,ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Execution of the rotation periodic boundary condition:
! (1) rotate POI and velocity vector
! (2) calc particle position after rotating according new velo and remaining time after POI
! (3) move particle from old element to new element
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos,Species,PartSpecies
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_TImeDisc_Vars          ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars          ,ONLY: VarTimeStep
USE MOD_Particle_Boundary_Vars ,ONLY: RotPeriodicSideMapping, NumRotPeriodicNeigh, SurfSide2RotPeriodicSide, GlobalSide2SurfSide
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_part_operations        ,ONLY: RemoveParticle
USE MOD_part_tools             ,ONLY: StoreLostParticleProperties
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticles, TrackInfo, CountNbrOfLostParts,DisplayLostParticles
USE MOD_DSMC_Vars              ,ONLY: DSMC, AmbipolElecVelo
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: PartID, SideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iNeigh, RotSideID, GlobalElemID
REAL                                 :: adaptTimeStep
LOGICAL                              :: FoundInElem
REAL                                 :: LastPartPos_old(1:3),Velo_old(1:3), Velo_oldAmbi(1:3)
!===================================================================================================================================

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     RotPeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' ParticlePosition: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' LastPartPos:      ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

IF (VarTimeStep%UseVariableTimeStep) THEN
  adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
ELSE
  adaptTimeStep = 1.0
END IF

! set last particle position on face (POI)
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha
LastPartPos_old(1:3) = LastPartPos(1:3,PartID)
Velo_old(1:3) = PartState(4:6,PartID)
IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) Velo_oldAmbi(1:3) = AmbipolElecVelo(PartID)%ElecVelo(1:3)
END IF
! (1) perform the rotational periodic movement and adjust velocity vector
ASSOCIATE(rot_alpha => PartBound%RotPeriodicAngle(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))), &
          rot_alpha_delta => PartBound%RotPeriodicAngle(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))) * 0.9999999)
  SELECT CASE(GEO%RotPeriodicAxi)
    CASE(1) ! x-rotation axis: LastPartPos(1,PartID) = LastPartPos_old(1,PartID)
      LastPartPos(2,PartID) = COS(rot_alpha_delta)*LastPartPos_old(2) - SIN(rot_alpha_delta)*LastPartPos_old(3)
      LastPartPos(3,PartID) = SIN(rot_alpha_delta)*LastPartPos_old(2) + COS(rot_alpha_delta)*LastPartPos_old(3)
      PartState(5,PartID) = COS(rot_alpha)*Velo_old(2) - SIN(rot_alpha)*Velo_old(3)
      PartState(6,PartID) = SIN(rot_alpha)*Velo_old(2) + COS(rot_alpha)*Velo_old(3)
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
          AmbipolElecVelo(PartID)%ElecVelo(5) = COS(rot_alpha)*Velo_oldAmbi(2) - SIN(rot_alpha)*Velo_oldAmbi(3)
          AmbipolElecVelo(PartID)%ElecVelo(6) = SIN(rot_alpha)*Velo_oldAmbi(2) + COS(rot_alpha)*Velo_oldAmbi(3)
        END IF
      END IF
    CASE(2) ! y-rotation axis: LastPartPos(2,PartID) = LastPartPos_old(2,PartID)
      LastPartPos(1,PartID) = COS(rot_alpha_delta)*LastPartPos_old(1) + SIN(rot_alpha_delta)*LastPartPos_old(3)
      LastPartPos(3,PartID) =-SIN(rot_alpha_delta)*LastPartPos_old(1) + COS(rot_alpha_delta)*LastPartPos_old(3)
      PartState(4,PartID) = COS(rot_alpha)*Velo_old(1) + SIN(rot_alpha)*Velo_old(3)
      PartState(6,PartID) =-SIN(rot_alpha)*Velo_old(1) + COS(rot_alpha)*Velo_old(3)
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
          AmbipolElecVelo(PartID)%ElecVelo(4) = COS(rot_alpha)*Velo_oldAmbi(1) + SIN(rot_alpha)*Velo_oldAmbi(3)
          AmbipolElecVelo(PartID)%ElecVelo(6) = -SIN(rot_alpha)*Velo_oldAmbi(1) + COS(rot_alpha)*Velo_oldAmbi(3)
        END IF
      END IF
    CASE(3) ! z-rotation axis: LastPartPos(3,PartID) = LastPartPos_old(3,PartID)
      LastPartPos(1,PartID) = COS(rot_alpha_delta)*LastPartPos_old(1) - SIN(rot_alpha_delta)*LastPartPos_old(2)
      LastPartPos(2,PartID) = SIN(rot_alpha_delta)*LastPartPos_old(1) + COS(rot_alpha_delta)*LastPartPos_old(2)
      PartState(4,PartID) = COS(rot_alpha)*Velo_old(1) - SIN(rot_alpha)*Velo_old(2)
      PartState(5,PartID) = SIN(rot_alpha)*Velo_old(1) + COS(rot_alpha)*Velo_old(2)
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
          AmbipolElecVelo(PartID)%ElecVelo(4) = COS(rot_alpha)*Velo_oldAmbi(1) - SIN(rot_alpha)*Velo_oldAmbi(2)
          AmbipolElecVelo(PartID)%ElecVelo(5) = SIN(rot_alpha)*Velo_oldAmbi(1) + COS(rot_alpha)*Velo_oldAmbi(2)
        END IF
      END IF
  END SELECT
END ASSOCIATE
! (2) update particle position after periodic BC
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - TrackInfo%alpha/TrackInfo%lengthPartTrajectory) * dt*RKdtFrac &
                        * PartState(4:6,PartID) * adaptTimeStep
! compute moved particle || rest of movement
TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
TrackInfo%lengthPartTrajectory= VECNORM(TrackInfo%PartTrajectory)
IF(ALMOSTZERO(TrackInfo%lengthPartTrajectory))THEN
  TrackInfo%lengthPartTrajectory= 0.0
ELSE
  TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory
END IF
#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     RotPeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' ParticlePosition-pp: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(1X,G0))') ' LastPartPosition-pp: ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/
! (3) move particle from old element to new element
RotSideID = SurfSide2RotPeriodicSide(GlobalSide2SurfSide(SURF_SIDEID,SideID))
DO iNeigh=1,NumRotPeriodicNeigh(RotSideID)
  ! find rotational periodic elem through localization in all potential rotational periodic sides
  GlobalElemID = RotPeriodicSideMapping(RotSideID,iNeigh)
  IF(GlobalElemID.EQ.-1) CALL abort(__STAMP__,' ERROR: Halo-rot-periodic side has no corresponding element.')
  CALL ParticleInsideQuad3D(LastPartPos(1:3,PartID),GlobalElemID,FoundInElem)
  IF(FoundInElem) THEN
    ElemID = GlobalElemID
    EXIT
  END IF
END DO
IF(.NOT.FoundInElem) THEN
  ! Particle appears to have not crossed any of the checked sides. Deleted!
  IF(DisplayLostParticles)THEN
    IPWRITE(*,*) 'Error in RotPeriodicBC! Particle Number',PartID,'lost. Element:', ElemID,'(species:',PartSpecies(PartID),')'
    IPWRITE(*,*) 'LastPos: ', LastPartPos(1:3,PartID)
    IPWRITE(*,*) 'Pos:     ', PartState(1:3,PartID)
    IPWRITE(*,*) 'Velo:    ', PartState(4:6,PartID)
    IPWRITE(*,*) 'Particle deleted!'
  END IF ! DisplayLostParticles
  IF(CountNbrOfLostParts)THEN
    CALL StoreLostParticleProperties(PartID,ElemID)
    NbrOfLostParticles=NbrOfLostParticles+1
  END IF ! CountNbrOfLostParts
  CALL RemoveParticle(PartID,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
END IF

END SUBROUTINE RotPeriodicBC

END MODULE MOD_Particle_Boundary_Condition
