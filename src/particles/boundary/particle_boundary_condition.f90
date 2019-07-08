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

MODULE MOD_Particle_Boundary_Condition
!===================================================================================================================================
!! Determines how particles interact with a given boundary condition
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
INTERFACE GetBoundaryInteraction
  MODULE PROCEDURE GetBoundaryInteraction
END INTERFACE

INTERFACE GetBoundaryInteractionRef
  MODULE PROCEDURE GetBoundaryInteractionRef
END INTERFACE

INTERFACE GetBoundaryInteractionAuxBC
  MODULE PROCEDURE GetBoundaryInteractionAuxBC
END INTERFACE

INTERFACE PartSwitchElement
  MODULE PROCEDURE PartSwitchElement
END INTERFACE

PUBLIC::GetBoundaryInteraction,GetBoundaryInteractionRef,GetBoundaryInteractionAuxBC,PartSwitchElement
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,locSideID,ElemID,crossedBC&
                                  ,TriNum)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!  PeriodicBC              = 3
!  SimpleAnodeBC           = 4
!  SimpleCathodeBC         = 5
!  MPINeighborhoodBC       = 6
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,KeepWallParticles, UseCircularInflow, UseAdaptive, Species
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,nPorousBC
USE MOD_Particle_Boundary_Porous, ONLY: PorousBoundaryTreatment
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Mesh_Vars,              ONLY:BC
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
#if defined(LSERK)
USE MOD_TimeDisc_Vars,          ONLY:RK_a!,iStage
#endif
#if defined(IMPA)
USE MOD_Particle_Vars,          ONLY:PartIsImplicit
USE MOD_Particle_Vars,          ONLY:DoPartInNewton
#endif /*IMPA*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip,locSideID
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN)                   :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3),RanNum
INTEGER                              :: adsorbindex, iSpec, iSF
LOGICAL                              :: isSpeciesSwap, PorousReflection
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) THEN
CALL abort(&
__STAMP__&
,' ERROR: PartBound not allocated!.',999,999.)
END IF
IsSpeciesSwap=.FALSE.
crossedBC    =.FALSE.
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!----------------------------------------------------------------------------------------------------------------------------------
  IF (.NOT.TriaTracking) THEN
    IF(alpha/lengthPartTrajectory.LE.epsilontol)THEN !if particle is close to BC, it encounters the BC only if it leaves element/grid
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,SideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      END SELECT
      IF(flip.NE.0) n_loc=-n_loc
      IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
    END IF
  END IF
  IF(UseAdaptive) THEN
    iSpec = PartSpecies(iPart)
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.PartBound%MapToPartBC(BC(SideID))) THEN
          IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)  Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut = &
                                                                  Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut + 1
      END IF
    END DO
  END IF
  IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
#ifdef IMPA
  DoPartInNewton(iPart) = .FALSE.
  PartIsImplicit(iPart) = .FALSE.
#endif /*IMPA*/

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  !---- Treatment of adaptive and porous boundary conditions (deletion of particles in case of circular inflow or porous BC)
  PorousReflection = .FALSE.
  IF(UseCircularInflow) CALL SurfaceFluxBasedBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,lengthPartTrajectory,flip,xi,eta)
  IF(nPorousBC.GT.0) CALL PorousBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,PorousReflection)
  !---- swap species?
  IF (PartBound%NbrOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID))).gt.0) THEN
#ifndef IMPA
    CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,flip=flip,TriNum=TriNum)
#else
    CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap)
#endif /*NOT IMPA*/
  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! particle did not Swap to species 0 !deleted particle -> particle swaped to species 0
    ! decide if reactive or simple sattering
    IF (.NOT.PartBound%Reactive(PartBound%MapToPartBC(BC(SideID)))) THEN
    ! simple reflection (previously used wall interaction model, maxwellian scattering)
      CALL RANDOM_NUMBER(RanNum)
      IF((RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))).OR.PorousReflection) THEN
        ! perfectly reflecting, specular re-emission
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
          IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
      ELSE
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
          IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
      END IF
    ELSE
      ! chemical surface interaction (e.g. adsorption)
      adsorbindex = 0
      ! Decide which interaction (reflection, reaction, adsorption)            
      CALL ReactiveSurfaceTreatment(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap,adsorbindex&
                                            ,opt_Reflected=crossedBC,TriNum=TriNum)
      ! assign right treatment
      SELECT CASE (adsorbindex)
      CASE(1)
        ! 1: adsorption (is either removed or set to be on surface)
        IF (KeepWallParticles.AND.(adsorbindex.EQ.1)) THEN
          PDM%ParticleAtWall(iPart) = .TRUE.
        ELSE
          IF(CalcPartBalance) THEN
            nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
            PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
          END IF ! CalcPartBalance
          PDM%ParticleInside(iPart) = .FALSE.
#ifdef IMPA
          PartIsImplicit(iPart) = .FALSE.
          DoPartInNewton(iPart) = .FALSE.
#endif /*IMPA*/
          alpha=-1.
        END IF
      CASE(2)
        ! 2: Eley-Rideal reaction (particle is reflected in catalytic treatment routine)
      CASE(0) ! inelastic reflection
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
          IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
      CASE(-1) ! elastic reflection
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
          IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
      CASE DEFAULT ! should not happen
        WRITE(*,*)'boundary_condition: Adsorption error. wrong interactionindex chosen'
        CALL Abort(&
__STAMP__,&
'Boundary_Error: Adsorptionindex switched to unknown value.')
      END SELECT
    END IF
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID &
                        ,ElemID,opt_perimoved=crossedBC,TriNum=TriNum) ! opt_reflected is peri-moved

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) !PartBound%SimpleAnodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) !PartBound%SimpleCathodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(10) !PartBound%SymmetryBC
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap &
                            ,opt_Symmetry=.TRUE.,opt_Reflected=crossedBC,TriNum=TriNum)
CASE(100) !PartBound%AnalyzeBC
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL  SideAnalysis(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,locSideID,ElemID &
                    ,IsSpeciesSwap,opt_crossed=crossedBC)

CASE DEFAULT
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (unknown case)',999,999.)
END SELECT !PartBound%MapToPartBC(BC(SideID)

! compiler warnings
IF(1.EQ.2)THEN
  WRITE(*,*) 'ElemID', ElemID
END IF

END SUBROUTINE GetBoundaryInteraction


SUBROUTINE GetBoundaryInteractionRef(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with a boundary condition
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!  PeriodicBC              = 3
!  SimpleAnodeBC           = 4
!  SimpleCathodeBC         = 5
!  MPINeighborhoodBC       = 6
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals!,                ONLY:Abort
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,KeepWallParticles
USE MOD_Particle_Boundary_Vars, ONLY:PartBound
USE MOD_Particle_Surfaces_Vars, ONLY:SideType,SideNormVec,epsilontol
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Mesh_Vars,              ONLY:BC,nSides
USE MOD_Particle_Tracking_Vars, ONLY:CartesianPeriodic
USE MOD_Particle_Mesh_Vars,     ONLY:PartBCSideList
#if defined(IMPA)
USE MOD_Particle_Vars,          ONLY:PartIsImplicit
USE MOD_Particle_Vars,          ONLY:DoPartInNewton
#endif /*IMPA*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip
REAL,INTENT(IN)                      :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
INTEGER,INTENT(INOUT)                :: ElemID
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: RanNum,n_loc(1:3)
INTEGER                              :: BCSideID, adsorbindex
LOGICAL                              :: IsSpeciesSwap
!===================================================================================================================================

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) THEN
CALL abort(&
__STAMP__&
,' ERROR: PartBound not allocated!.',999,999.)
END IF
IsSpeciesSwap=.FALSE.
crossedBC    =.FALSE.
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartBound%OpenBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(alpha/lengthPartTrajectory.LE.epsilontol)THEN !if particle is close to BC, it encounters the BC only if it leaves element/grid
    BCSideID=PartBCSideList(SideID)
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,BCSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
  END IF

  IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.
#ifdef IMPA
  DoPartInNewton(iPart) = .FALSE.
  PartIsImplicit(iPart) = .FALSE.
#endif /*IMPA*/
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartBound%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  !---- swap species?
  BCSideID=PartBCSideList(SideID)
  IF (PartBound%NbrOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID))).gt.0) THEN
#ifndef IMPA
    CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,flip,BCSideID=BCSideID)
#else
    CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap)
#endif /*NOT IMPA*/
  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! particle did not Swap to species 0 !deleted particle -> particle swaped to species 0
    BCSideID=PartBCSideList(SideID)
    IF (.NOT.PartBound%Reactive(PartBound%MapToPartBC(BC(SideID)))) THEN
    ! simple reflection (previously used wall interaction model, maxwellian scattering)
      CALL RANDOM_NUMBER(RanNum)
      IF(RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))) THEN
        ! perfectly reflecting, specular re-emission
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap &
                              ,BCSideID=BCSideID,opt_reflected=crossedBC)
      ELSE
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap&
                              ,BCSideID=BCSideID,opt_reflected=crossedBC)
      END IF
    ELSE
      ! chemical surface interaction (such as adsorption)
      adsorbindex = 0
      ! Decide which interaction (reflection, reaction, adsorption)
      CALL ReactiveSurfaceTreatment(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap,adsorbindex&
                              ,BCSideID=BCSideID,opt_reflected=crossedBC)
      ! assign right treatment
      SELECT CASE (adsorbindex)
      CASE(1)
        ! 1: adsorption (is either removed or set to be on surface)
        IF (KeepWallParticles.AND.(adsorbindex.EQ.1)) THEN
          PDM%ParticleAtWall(iPart) = .TRUE.
        ELSE
          IF(CalcPartBalance) THEN
            nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
            PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
          END IF ! CalcPartBalance
          PDM%ParticleInside(iPart) = .FALSE.
#ifdef IMPA
          PartIsImplicit(iPart) = .FALSE.
          DoPartInNewton(iPart) = .FALSE.
#endif /*IMPA*/
          alpha=-1.
        END IF
      CASE(2)
        ! 2: Eley-Rideal reaction (particle is reflected in catalytic treatment routine)
        !CALL Particle_ER_Reflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap)
      CASE(0) ! inelastic reflection
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap&
                              ,BCSideID=BCSideID,opt_reflected=crossedBC)
      CASE(-1) !elastic reflection
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap&
                              ,BCSideID=BCSideID,opt_reflected=crossedBC)
      CASE DEFAULT ! should not happen
        WRITE(*,*)'Boundary_PIC: Adsorption error.'
        CALL Abort(&
__STAMP__,&
'Boundary_Error: Adsorptionindex switched to unknown value.')
      END SELECT
    END IF
  ELSE
    ! not inside any-more, removed in last step
    crossedBC=.TRUE.
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) !PartBound%PeriodicBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! sanity check
  IF(CartesianPeriodic) CALL abort(&
__STAMP__&
,' No periodic BCs for CartesianPeriodic!')
  ! move particle periodic distance
  BCSideID=PartBCSideList(SideID)
  CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID &
                        ,BCSideID=BCSideID,opt_perimoved=crossedBC) ! opt_reflected is peri-moved

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) !PartBound%SimpleAnodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5) !PartBound%SimpleCathodeBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6) !PartBound%MPINeighborhoodBC)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)',999,999.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(10) !PartBound%SymmetryBC
!-----------------------------------------------------------------------------------------------------------------------------------
  BCSideID=PartBCSideList(SideID)
  CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap &
                        ,BCSideID=BCSideID,opt_Symmetry=.TRUE.,opt_reflected=crossedBC)

CASE DEFAULT
CALL abort(&
__STAMP__&
,' ERROR: PartBound not associated!. BC(SideID)',BC(SideID),REAL(SideID/nSides))
END SELECT !PartBound%MapToPartBC(BC(SideID)

END SUBROUTINE GetBoundaryInteractionRef


SUBROUTINE GetBoundaryInteractionAuxBC(PartTrajectory,lengthPartTrajectory,alpha,iPart,AuxBCIdx,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with an auxBC
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,                ONLY:Abort
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies
USE MOD_Particle_Boundary_Vars, ONLY:PartAuxBC
!USE MOD_Particle_Surfaces_vars, ONLY:epsilontol
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
#if defined(LSERK)
USE MOD_TimeDisc_Vars,          ONLY:RK_a
#endif
!USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,AuxBCIdx
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: RanNum
LOGICAL                              :: isSpeciesSwap
!===================================================================================================================================

IsSpeciesSwap=.FALSE.
crossedBC    =.FALSE.
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartAuxBC%TargetBoundCond(AuxBCIdx))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartAuxBC%OpenBC
!-----------------------------------------------------------------------------------------------------------------------------------
!  IF(alpha/lengthPartTrajectory.LE.epsilontol)THEN !if particle is close to BC, it encounters the BC only if it leaves element/grid
!    IF (.NOT.TriaTracking) THEN
!      SELECT CASE(SideType(SideID))
!      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
!        n_loc=SideNormVec(1:3,SideID)
!      CASE(BILINEAR)
!        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
!      CASE(CURVED)
!        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
!      END SELECT
!      IF(flip.NE.0) n_loc=-n_loc
!      IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
!    END IF
!  END IF
  IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
      PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
  END IF ! CalcPartBalance
  PDM%ParticleInside(iPart) = .FALSE.
  alpha=-1.

!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartAuxBC%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  !---- swap species?
!print*,'*********************'
!print*,AuxBCIdx
!print*,iPart,alpha,PartState(iPart,4:6)
!print*,iPart,alpha,LastPartPos(iPart,1:3),PartState(iPart,1:3)
  IF (PartAuxBC%NbrOfSpeciesSwaps(AuxBCIdx).gt.0) THEN
! CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1, &
!      IsSpeciesSwap=IsSpeciesSwap,flip=-1,AuxBCIdx=AuxBCIdx)
#ifndef IMPA
    CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,IsSpeciesSwap=IsSpeciesSwap, &
      flip=-1,AuxBCIdx=AuxBCIdx)
#else
    CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,IsSpeciesSwap=IsSpeciesSwap)
#endif /*NOT IMPA*/

  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! particle did not Swap to species 0 !deleted particle -> particle swaped to species 0
      ! simple reflection (previously used wall interaction model, maxwellian scattering)
        CALL RANDOM_NUMBER(RanNum)
        IF(RanNum.GE.PartAuxBC%MomentumACC(AuxBCIdx)) THEN
          ! perfectly reflecting, specular re-emission
          CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,flip=-1, &
            IsSpeciesSwap=IsSpeciesSwap,opt_Reflected=crossedBC,AuxBCIdx=AuxBCIdx)
        ELSE
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,flip=-1, &
            IsSpeciesSwap=IsSpeciesSwap,opt_Reflected=crossedBC,AuxBCIdx=AuxBCIdx)
        END IF
  END IF
!print*,iPart,alpha,LastPartPos(iPart,1:3),PartState(iPart,1:3)
!print*,iPart,alpha,PartState(iPart,4:6)
!print*,'*********************'
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
CALL abort(&
__STAMP__&
,' ERROR: AuxBC bound not associated!. (unknown case)',999,999.)
END SELECT

! compiler warnings
IF(1.EQ.2)THEN
  WRITE(*,*) 'AuxBCIdx', AuxBCIdx
END IF

END SUBROUTINE GetBoundaryInteractionAuxBC


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,IsSpeciesSwap,BCSideID, &
  opt_Symmetry,opt_Reflected,TriNum,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,SurfMesh,SampWall,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars, ONLY:dXiEQ_SurfSample
USE MOD_Particle_Mesh_Vars,     ONLY:epsInCell
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,nSpecies,PartSpecies,Species,WriteMacroSurfaceValues,PartLorentzType
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_DSMC_Vars,              ONLY:DSMC
USE MOD_LD_Vars,                ONLY:useLD
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
USE MOD_TImeDisc_Vars,          ONLY:tend,time
USE MOD_Particle_Boundary_Vars, ONLY:AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
USE MOD_Equation_Vars,          ONLY:c2_inv
#if defined(LSERK)
USE MOD_Particle_Vars,          ONLY:Pt_temp,PDM
#elif (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars,          ONLY:PDM
#endif
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars,          ONLY:PEM
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID, flip!,ElemID
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_old(1:3),n_loc(1:3),WallVelo(3),intersec(3),r_vec(3),axis(3),cos2inv!,v_2(1:3),v_aux(1:3)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
!#if defined(LSERK)
!REAL                                 :: absPt_temp
!#endif
!REAL,PARAMETER                       :: oneMinus=0.99999999
!REAL                                 :: oneMinus!=0.99999999
REAL                                 :: LorentzFac, LorentzFacInv
REAL                                 :: epsLength
REAL                                 :: Xitild,EtaTild
INTEGER                              :: p,q, SurfSideID, locBCID
LOGICAL                              :: Symmetry, IsAuxBC
!===================================================================================================================================

IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF
IF (IsAuxBC) THEN
  SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
  CASE ('plane')
    n_loc = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
  CASE ('cylinder')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
    IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE ('cone')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
    cos2inv = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
    IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE ('parabol')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
    IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END SELECT
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    !RETURN
    CALL abort(&
      __STAMP__&
      ,'Error in PerfectReflection: Particle coming from outside!')
  ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  ELSE
    CALL abort(&
      __STAMP__&
      ,'Error in PerfectReflection: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
  END IF
  WallVelo=PartAuxBC%WallVelo(1:3,AuxBCIdx)
ELSE
  !OneMinus=1.0-MAX(epsInCell,epsilontol)
  epsLength=MAX(epsInCell,epsilontol)*lengthPartTrajectory
  WallVelo=PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
  locBCID=PartBound%MapToPartBC(BC(SideID))

  IF(PRESENT(BCSideID))THEN
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,BCSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
    END SELECT
  ELSE
    IF (TriaTracking) THEN
      CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
    ELSE
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,SideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
      END SELECT
      IF(flip.NE.0) n_loc=-n_loc
    END IF
  END IF
  IF(PRESENT(opt_Symmetry)) THEN
    Symmetry = opt_Symmetry
  ELSE
    Symmetry = .FALSE.
  END IF
  IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    RETURN
  ELSE
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  END IF
END IF !IsAuxBC

! SHOULD NOT BE NEEDED
!#if defined(ROS)
!IF(iStage.GT.0)THEN
!  !IF(iStage.GT.2)THEN
!  IF(RK_inflow(iStage).EQ.0)THEN
!    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
!    ! hence, we perform an evil hack and beam the particle on the intersection point
!    ! this ensures, that the particle is located within the mesh....
!    LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*0.999*alpha
!    PartState(PartID,1:3)   = LastPartPos(PartID,1:3)
!    PartTrajectory          = 0.
!    lengthPartTrajectory    = 0.
!    RETURN
!  END IF
!END IF
!#endif /*ROS*/

IF(SUM(ABS(WallVelo)).GT.0.)THEN
  SELECT CASE(PartLorentzType)
  CASE(3)
    v_old = PartState(PartID,4:6)
    PartState(PartID,4:6) = PartState(PartID,4:6) &
                          - 2.*DOT_PRODUCT(PartState(PartID,4:6),n_loc)*n_loc + WallVelo
    ! sanity check of new particle velocity
    LorentzFac=1.0-DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))*c2_inv
    IF(LorentzFac.LT.0.) CALL Abort(&
__STAMP__&
,'Particle exceeds speed of light! PartID ',PartID)
  CASE(5)
    ! map relativistic momentum to velocity
    LorentzFacInv         = 1.0+DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))*c2_inv
    LorentzFacInv         = 1.0/SQRT(LorentzFacInv)
    PartState(PartID,4)   = LorentzFacInv*PartState(PartID,4)
    PartState(PartID,5)   = LorentzFacInv*PartState(PartID,5)
    PartState(PartID,6)   = LorentzFacInv*PartState(PartID,6)
    v_old                 = PartState(PartID,4:6)
    ! update velocity
    PartState(PartID,4:6) = PartState(PartID,4:6) &
                          - 2.*DOT_PRODUCT(PartState(PartID,4:6),n_loc)*n_loc + WallVelo
    ! map back from velocity to relativistic momentum
    LorentzFac=1.0-DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6))*c2_inv
    IF(LorentzFac.LT.0.)THEN
CALL Abort(&
__STAMP__&
,'Particle exceeds speed of light! PartID ',PartID)
    END IF
    LorentzFac=1.0/SQRT(LorentzFac)
    PartState(PartID,4)   = LorentzFac*PartState(PartID,4)
    PartState(PartID,5)   = LorentzFac*PartState(PartID,5)
    PartState(PartID,6)   = LorentzFac*PartState(PartID,6)
  CASE DEFAULT
    v_old = PartState(PartID,4:6)
    PartState(PartID,4:6) = PartState(PartID,4:6) &
                          - 2.*DOT_PRODUCT(PartState(PartID,4:6),n_loc)*n_loc + WallVelo
  END SELECT
ELSE
  v_old = PartState(PartID,4:6)
  PartState(PartID,4:6) = PartState(PartID,4:6) &
                        - 2.*DOT_PRODUCT(PartState(PartID,4:6),n_loc)*n_loc
END IF

IF (.NOT.IsAuxBC) THEN
! Wall sampling Macrovalues
IF((.NOT.Symmetry).AND.(.NOT.UseLD)) THEN !surface mesh is not build for the symmetry BC!?!
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
    ! compute p and q
    ! correction of xi and eta, can only be applied if xi & eta are not used later!
    IF (TriaTracking) THEN
      p=1 ; q=1
    ELSE
      Xitild =MIN(MAX(-1.,xi ),0.99)
      Etatild=MIN(MAX(-1.,eta),0.99)
      p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
      q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
    END IF

  !----  Sampling Forces at walls
!       SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
!                                          * (v_old(1:3) - PartState(PartID,4:6)) * Species(PartSpecies(PartID))%MacroParticleFactor
    SampWall(SurfSideID)%State(10,p,q)= SampWall(SurfSideID)%State(10,p,q) + Species(PartSpecies(PartID))%MassIC &
                                        * (v_old(1) - PartState(PartID,4)) * Species(PartSpecies(PartID))%MacroParticleFactor
    SampWall(SurfSideID)%State(11,p,q)= SampWall(SurfSideID)%State(11,p,q) + Species(PartSpecies(PartID))%MassIC &
                                        * (v_old(2) - PartState(PartID,5)) * Species(PartSpecies(PartID))%MacroParticleFactor
    SampWall(SurfSideID)%State(12,p,q)= SampWall(SurfSideID)%State(12,p,q) + Species(PartSpecies(PartID))%MassIC &
                                        * (v_old(3) - PartState(PartID,6)) * Species(PartSpecies(PartID))%MacroParticleFactor
  !---- Counter for collisions (normal wall collisions - not to count if only Swaps to be counted, IsSpeciesSwap: already counted)
!       IF (.NOT.CalcSurfCollis%OnlySwaps) THEN
    IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
      SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q) = SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q) + 1
      IF (CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))) THEN
        AnalyzeSurfCollis%Number(PartSpecies(PartID)) = AnalyzeSurfCollis%Number(PartSpecies(PartID)) + 1
        AnalyzeSurfCollis%Number(nSpecies+1) = AnalyzeSurfCollis%Number(nSpecies+1) + 1
        IF (AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber) THEN
CALL Abort(&
__STAMP__&
,'maxSurfCollisNumber reached!')
        END IF
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) &
          = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
        !-- caution: for consistency with diffuse refl. v_old is used!
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) &
          = v_old(1)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) &
          = v_old(2)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) &
          = v_old(3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) &
          = LastPartPos(PartID,1)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) &
          = LastPartPos(PartID,2)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) &
          = LastPartPos(PartID,3)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = locBCID
      END IF
    END IF
  END IF
END IF
END IF !.NOT.IsAuxBC

! set particle position on face
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha

PartTrajectory(1:3)=PartTrajectory(1:3)-2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc
PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)

! #if !defined(IMPA) &&  !defined(ROS)
! compute moved particle || rest of movement
PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory
! #endif

#if defined(LSERK) || (PP_TimeDiscMethod==509)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
   ! correction for Runge-Kutta (correct position!!)
!---------- old ----------
!  absPt_temp=SQRT(Pt_temp(PartID,1)*Pt_temp(PartID,1)+Pt_temp(PartID,2)*Pt_temp(PartID,2)+Pt_temp(PartID,3)*Pt_temp(PartID,3))
!  ! scale PartTrajectory to new Pt_temp
!  Pt_temp(PartID,1:3)=absPt_temp*PartTrajectory(1:3)
!  ! deleate force history
!  Pt_temp(PartID,4:6)=0.
!  ! what happens with force term || acceleration?
!-------------------------
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#if defined(LSERK)
ELSE
  Pt_temp(PartID,1:3)=Pt_temp(PartID,1:3)-2.*DOT_PRODUCT(Pt_temp(PartID,1:3),n_loc)*n_loc
  IF (Symmetry) THEN !reflect also force history for symmetry
    Pt_temp(PartID,4:6)=Pt_temp(PartID,4:6)-2.*DOT_PRODUCT(Pt_temp(PartID,4:6),n_loc)*n_loc
  ELSE
    Pt_temp(PartID,4:6)=0. !produces best result compared to analytical solution in plate capacitor...
  END IF
#endif  /*LSERK*/
END IF
#endif  /*LSERK || (PP_TimeDiscMethod==509)*/

! rotation for IMEX and Rosenbrock Method (requires the rotation of the previous rk-stages... simplification of boundary condition)
! results in an order reduction
#ifdef IMPA
!IF(SUM(ABS(PEM%NormVec(PartID,1:3))).GT.0)THEN
!   IPWRITE(*,*) ' Caution: Field rotation for several reflection is not implemented!', iStage,PartIsImplicit(PartID), PartID
! END IF
PEM%NormVec(PartID,1:3)=n_loc
#endif /*IMPA*/
#ifdef ROS
! IF(SUM(ABS(PEM%NormVec(PartID,1:3))).GT.0)THEN
!   !IPWRITE(*,*) ' Caution: Field rotation for several reflection is not implemented!'
! END IF
PEM%NormVec(PartID,1:3)=n_loc
#endif /*ROS*/

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,IsSpeciesSwap,BCSideID &
  ,opt_Reflected,TriNum,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
! only implemented for DoRefMapping tracking
! PartBCs are reduced!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals,                ONLY:CROSSNORM,abort,UNITVECTOR
USE MOD_Globals_Vars,           ONLY:PI, BoltzmannConst
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,SurfMesh,SampWall,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars, ONLY:dXiEQ_SurfSample
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,Species,PartSpecies,nSpecies,WriteMacroSurfaceValues
#if defined(LSERK) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars,          ONLY:PDM
#endif
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,BezierControlPoints3D
USE MOD_Mesh_Vars,              ONLY:BC,NGEO
USE MOD_DSMC_Vars,              ONLY:SpecDSMC,CollisMode
USE MOD_DSMC_Vars,              ONLY:PartStateIntEn,DSMC, useDSMC
USE MOD_DSMC_Vars,              ONLY:PolyatomMolDSMC, VibQuantsPar
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues
USE MOD_TimeDisc_Vars,          ONLY:dt,tend,time,RKdtFrac
USE MOD_Particle_Boundary_Vars, ONLY:AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
#elif (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Vars,            ONLY: FPDoVibRelaxation
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID, flip
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: Opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: locBCID, vibQuant, vibQuantNew, VibQuantWall
REAL                                 :: VibQuantNewR                                                !
!REAL,PARAMETER                       :: oneMinus=0.99999999
REAL                                 :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
REAL                                 :: EtraWall, EtraNew
REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
REAL                                 :: n_loc(1:3), tang1(1:3),tang2(1:3), NewVelo(3)
REAL                                 :: ErotNew, ErotWall, EVibNew, Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                                 :: Xitild,EtaTild
!REAL                                 :: WallTransACC
INTEGER                              :: p,q, SurfSideID
REAL                                 :: POI_fak,TildTrajectory(3)
! Polyatomic Molecules
REAL, ALLOCATABLE                    :: RanNumPoly(:), VibQuantNewRPoly(:)
REAL                                 :: NormProb
INTEGER                              :: iPolyatMole, iDOF
INTEGER, ALLOCATABLE                 :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
! REAL, ALLOCATABLE                    :: VecXVibPolyFP(:), VecYVibPolyFP(:), CmrVibPolyFP(:)
! REAL, ALLOCATABLE                    :: EVPolyNewFP(:), EVPolyWallFP(:)
!REAL                                 :: ErotOldPoly(3), ErotNewPoly(3), ErotWallPoly(3), CmrRotPoly(3)
LOGICAL                              :: IsAuxBC
REAL                                 :: intersec(3), r_vec(3), axis(3), cos2inv
!===================================================================================================================================
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF
IF (IsAuxBC) THEN
  SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
  CASE ('plane')
    n_loc = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
  CASE ('cylinder')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
    IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE ('cone')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
    cos2inv = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
    IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE ('parabol')
    intersec = LastPartPos(PartID,1:3) + alpha*PartTrajectory
    r_vec = AuxBC_parabol(AuxBCMap(AuxBCIdx))%r_vec
    axis  = AuxBC_parabol(AuxBCMap(AuxBCIdx))%axis
    n_loc = UNITVECTOR( intersec - ( r_vec + axis*(DOT_PRODUCT(intersec-r_vec,axis)+0.5*AuxBC_parabol(AuxBCMap(AuxBCIdx))%zfac) ) )
    IF (.NOT.AuxBC_parabol(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
  CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,'AuxBC does not exist')
  END SELECT
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    !RETURN
    CALL abort(&
      __STAMP__&
      ,'Error in DiffuseReflection: Particle coming from outside!')
  ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  ELSE
    CALL abort(&
      __STAMP__&
      ,'Error in DiffuseReflection: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
  END IF
  WallVelo   = PartAuxBC%WallVelo(1:3,AuxBCIdx)
  WallTemp   = PartAuxBC%WallTemp(AuxBCIdx)
  TransACC   = PartAuxBC%TransACC(AuxBCIdx)
  VibACC     = PartAuxBC%VibACC(AuxBCIdx)
  RotACC     = PartAuxBC%RotACC(AuxBCIdx)
  IF (n_loc(3).NE.0.) THEN
    tang1(1) = 1.0
    tang1(2) = 1.0
    tang1(3) = -(n_loc(1)+n_loc(2))/n_loc(3)
  ELSE
    IF (n_loc(2).NE.0.) THEN
      tang1(1) = 1.0
      tang1(3) = 1.0
      tang1(2) = -(n_loc(1)+n_loc(3))/n_loc(2)
    ELSE
      IF (n_loc(1).NE.0.) THEN
        tang1(2) = 1.0
        tang1(3) = 1.0
        tang1(1) = -(n_loc(2)+n_loc(3))/n_loc(1)
      ELSE
        CALL abort(&
__STAMP__&
,'Error in DiffuseReflection, n_vec is zero for AuxBC',AuxBCIdx)
      END IF
    END IF
  END IF
  tang1=UNITVECTOR(tang1)
  tang2=CROSSNORM(n_loc,tang1)
ELSE
  !OneMinus=1.0-epsInCell

  ! additional states
  locBCID=PartBound%MapToPartBC(BC(SideID))
  ! get BC values
  WallVelo   = PartBound%WallVelo(1:3,locBCID)
  WallTemp   = PartBound%WallTemp(locBCID)
  TransACC   = PartBound%TransACC(locBCID)
  VibACC     = PartBound%VibACC(locBCID)
  RotACC     = PartBound%RotACC(locBCID)

  IF(PRESENT(BCSideID))THEN
    SELECT CASE(SideType(BCSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,BCSideID)
      tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID))
      tang2=CROSSNORM(n_loc,tang1)
      !tang2=BezierControlPoints3D(:,0,NGeo,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
      !   CALL abort(__STAMP__'nvec for bezier not implemented!',999,999.)
    END SELECT
  ELSE
    IF (TriaTracking) THEN
      CALL CalcNormAndTangTriangle(nVec=n_loc,tang1=tang1,tang2=tang2, &
          TriNum=TriNum,SideID=SideID)
    ELSE
      SELECT CASE(SideType(SideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
        n_loc=SideNormVec(1:3,SideID)
        tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,SideID)-BezierControlPoints3D(:,0,0,SideID))
        tang2=CROSSNORM(n_loc,tang1)
        !tang2=BezierControlPoints3D(:,0,NGeo,SideID)-BezierControlPoints3D(:,0,0,SideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,SideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,SideID)
        !   CALL abort(__STAMP__'nvec for bezier not implemented!',999,999.)
      END SELECT
      IF(flip.NE.0) n_loc=-n_loc
    END IF
  END IF

  IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
    IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
    RETURN
  ELSE
    IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
  END IF
END IF !IsAuxBC

! calculate new velocity vector (Extended Maxwellian Model)
VeloReal = SQRT(PartState(PartID,4) * PartState(PartID,4) + &
                PartState(PartID,5) * PartState(PartID,5) + &
                PartState(PartID,6) * PartState(PartID,6))

EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
CALL RANDOM_NUMBER(RanNum)
VeloCrad    = SQRT(-LOG(RanNum))
CALL RANDOM_NUMBER(RanNum)
VeloCz      = SQRT(-LOG(RanNum))
Fak_D       = VeloCrad**2 + VeloCz**2

EtraWall    = BoltzmannConst * WallTemp * Fak_D
EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
Cmr         = SQRT(2.0 * EtraNew / (Species(PartSpecies(PartID))%MassIC * Fak_D))
CALL RANDOM_NUMBER(RanNum)
Phi     = 2.0 * PI * RanNum
VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
VeloCz  = Cmr * VeloCz

IF (.NOT.IsAuxBC) THEN
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    !----  Sampling for energy (translation) accommodation at walls
    ! has to be corrected to new scheme
    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
    ! compute p and q
    ! correction of xi and eta, can only be applied if xi & eta are not used later!
    IF (TriaTracking) THEN
      p=1 ; q=1
    ELSE
      Xitild =MIN(MAX(-1.,xi ),0.99)
      Etatild=MIN(MAX(-1.,eta),0.99)
      p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
      q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
    END IF

    SampWall(SurfSideID)%State(1,p,q)= SampWall(SurfSideID)%State(1,p,q)+EtraOld       &
      *Species(PartSpecies(PartID))%MacroParticleFactor
    SampWall(SurfSideID)%State(2,p,q)= SampWall(SurfSideID)%State(2,p,q)+EtraWall      &
      *Species(PartSpecies(PartID))%MacroParticleFactor
    SampWall(SurfSideID)%State(3,p,q)= SampWall(SurfSideID)%State(3,p,q)+EtraNew       &
      *Species(PartSpecies(PartID))%MacroParticleFactor
  END IF
END IF !.NOT.IsAuxBC

!   Transformation local distribution -> global coordinates
! from flux comutaion
! v = nv*u+t1*v+t2*f3

!NewVelo = VeloCx*tang1+CROSS(-n_loc,tang1)*VeloCy-VeloCz*n_loc
NewVelo = VeloCx*tang1-tang2*VeloCy-VeloCz*n_loc

IF (.NOT.IsAuxBC) THEN !so far no internal DOF stuff for AuxBC!!!
!---- Internal energy accommodation
IF (useDSMC) THEN
IF (CollisMode.GT.1) THEN
IF ((SpecDSMC(PartSpecies(PartID))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartID))%InterID.EQ.20)) THEN

  !---- Rotational energy accommodation
    IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.2) THEN
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
    ELSE IF (SpecDSMC(PartSpecies(PartID))%Xi_Rot.EQ.3) THEN
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = RanNum*10. !the distribution function has only non-negligible  values betwenn 0 and 10
      NormProb = SQRT(ErotWall)*EXP(-ErotWall)/(SQRT(0.5)*EXP(-0.5))
      CALL RANDOM_NUMBER(RanNum)
      DO WHILE (RanNum.GE.NormProb)
        CALL RANDOM_NUMBER(RanNum)
        ErotWall = RanNum*10. !the distribution function has only non-negligible  values betwenn 0 and 10
        NormProb = SQRT(ErotWall)*EXP(-ErotWall)/(SQRT(0.5)*EXP(-0.5))
        CALL RANDOM_NUMBER(RanNum)
      END DO
      ErotWall = ErotWall*BoltzmannConst*WallTemp
    END IF
    ErotNew  = PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))

    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    !----  Sampling for internal energy accommodation at walls
      SampWall(SurfSideID)%State(4,p,q)=SampWall(SurfSideID)%State(4,p,q)+PartStateIntEn(PartID,2) &
                                                                          *Species(PartSpecies(PartID))%MacroParticleFactor
      SampWall(SurfSideID)%State(5,p,q)=SampWall(SurfSideID)%State(5,p,q)+ErotWall &
                                                                          * Species(PartSpecies(PartID))%MacroParticleFactor
      SampWall(SurfSideID)%State(6,p,q)=SampWall(SurfSideID)%State(6,p,q)+ErotNew &
                                                                          * Species(PartSpecies(PartID))%MacroParticleFactor
    END IF

    PartStateIntEn(PartID,2) = ErotNew

#if (PP_TimeDiscMethod==400)
    IF (BGKDoVibRelaxation) THEN
#elif (PP_TimeDiscMethod==300)
    IF (FPDoVibRelaxation) THEN
#endif
   !---- Vibrational energy accommodation
      IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) THEN
        EvibNew = 0.0
        iPolyatMole = SpecDSMC(PartSpecies(PartID))%SpecToPolyArray
        ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
                 VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
                 VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        END DO
        VibQuantNewRPoly(:) = VibQuantsPar(PartID)%Quants(:) + VibACC*(VibQuantWallPoly(:) - VibQuantsPar(PartID)%Quants(:))
        VibQuantNewPoly = INT(VibQuantNewRPoly)
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          CALL RANDOM_NUMBER(RanNum)
          IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
            EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
                      * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
          ELSE
            EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
                      * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
          END IF
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib) &
                    - DSMC%GammaQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
        DO WHILE (VibQuantWall.GE.SpecDSMC(PartSpecies(PartID))%MaxVibQuant)
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
        END DO
        VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
        VibQuantNew = INT(VibQuantNewR)
        CALL RANDOM_NUMBER(RanNum)
        IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
          EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
        ELSE
          EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
        END IF
      END IF

      IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
   !----  Sampling for internal energy accommodation at walls
        IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(PartID))%SpecToPolyArray
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            SampWall(SurfSideID)%State(7,p,q)= SampWall(SurfSideID)%State(7,p,q) &
                    + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(PartSpecies(PartID))%MacroParticleFactor
            SampWall(SurfSideID)%State(8,p,q)= SampWall(SurfSideID)%State(8,p,q) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(PartSpecies(PartID))%MacroParticleFactor
          END DO
        ELSE
          SampWall(SurfSideID)%State(7,p,q)= SampWall(SurfSideID)%State(7,p,q) + (VibQuant + DSMC%GammaQuant) &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * Species(PartSpecies(PartID))%MacroParticleFactor
          SampWall(SurfSideID)%State(8,p,q)= SampWall(SurfSideID)%State(8,p,q) + VibQuantWall &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * Species(PartSpecies(PartID))%MacroParticleFactor
        END IF
        SampWall(SurfSideID)%State(9,p,q)= SampWall(SurfSideID)%State(9,p,q) + &
                                            EvibNew * Species(PartSpecies(PartID))%MacroParticleFactor
      END IF
      IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) VibQuantsPar(PartID)%Quants(:) = VibQuantTemp(:)
      PartStateIntEn(PartID,1) = EvibNew
#if ((PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==300))
    END IF
#endif
    END IF
    END IF ! CollisMode > 1
    END IF ! useDSMC

IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
!----  Sampling force at walls
  SampWall(SurfSideID)%State(10,p,q)= SampWall(SurfSideID)%State(10,p,q) &
      + Species(PartSpecies(PartID))%MassIC * (PartState(PartID,4) - NewVelo(1)) * Species(PartSpecies(PartID))%MacroParticleFactor
  SampWall(SurfSideID)%State(11,p,q)= SampWall(SurfSideID)%State(11,p,q) &
      + Species(PartSpecies(PartID))%MassIC * (PartState(PartID,5) - NewVelo(2)) * Species(PartSpecies(PartID))%MacroParticleFactor
  SampWall(SurfSideID)%State(12,p,q)= SampWall(SurfSideID)%State(12,p,q) &
      + Species(PartSpecies(PartID))%MassIC * (PartState(PartID,6) - NewVelo(3)) * Species(PartSpecies(PartID))%MacroParticleFactor
 !---- Counter for collisions (normal wall collisions - not to count if only SpeciesSwaps to be counted)
  IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
    SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q)= SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q) +1
    IF (CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))) THEN
      AnalyzeSurfCollis%Number(PartSpecies(PartID)) = AnalyzeSurfCollis%Number(PartSpecies(PartID)) + 1
      AnalyzeSurfCollis%Number(nSpecies+1) = AnalyzeSurfCollis%Number(nSpecies+1) + 1
      IF (AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber) THEN
CALL Abort(&
__STAMP__&
,'maxSurfCollisNumber reached!')
      END IF
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) &
        = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) &
        = PartState(PartID,4)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) &
        = PartState(PartID,5)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) &
        = PartState(PartID,6)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) &
        = LastPartPos(PartID,1)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) &
        = LastPartPos(PartID,2)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) &
        = LastPartPos(PartID,3)
      AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) &
        = PartSpecies(PartID)
      AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) &
        = locBCID
    END IF
  END IF
END IF
END IF !.NOT.IsAuxBC

! intersection point with surface
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha

! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
!TildPos       =PartState(PartID,1:3)-dt*RKdtFrac*PartState(PartID,4:6)
TildTrajectory=dt*RKdtFrac*PartState(PartID,4:6)
POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
! travel rest of particle vector
!PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac * NewVelo(1:3)
IF (IsAuxBC) THEN
  IF (PartAuxBC%Resample(AuxBCIdx)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
ELSE
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
END IF
PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3)

!----  saving new particle velocity
PartState(PartID,4:6)   = NewVelo(1:3) + WallVelo(1:3)

! recompute trajectory etc
PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory
!lengthPartTrajectory=lengthPartTrajectory!+epsilontol

#if defined(LSERK) || (PP_TimeDiscMethod==509)
PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#endif

END SUBROUTINE DiffuseReflection

#ifndef IMPA
SUBROUTINE SpeciesSwap(PartTrajectory,alpha,xi,eta,PartID,SideID,IsSpeciesSwap,flip,BCSideID,TriNum,AuxBCIdx)
#else
SUBROUTINE SpeciesSwap(PartTrajectory,alpha,xi,eta,PartID,SideID,IsSpeciesSwap,AuxBCIdx)
#endif
    !----------------------------------------------------------------------------------------------------------------------------------!
! Computes the Species Swap on ReflectiveBC
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals,                ONLY:abort
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,SampWall,dXiEQ_SurfSample,SurfMesh,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,PartSpecies,PDM
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues,nSpecies,CollectCharges,nCollectChargesBCs,Species
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Particle_Analyze,       ONLY: CalcEkinPart
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_DSMC_Vars,              ONLY:DSMC
USE MOD_TimeDisc_Vars,          ONLY:TEnd,Time
#if defined(IMPA)
USE MOD_Particle_Vars,          ONLY:PartIsImplicit,DoPartInNewton
#endif /*IMPA*/
#ifndef IMPA
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType
#endif /*NOT IMPA*/

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID
#ifndef IMPA
INTEGER,INTENT(IN),OPTIONAL       :: flip,BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
#endif /*NOT IMPA*/
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: IsSpeciesSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: targetSpecies, iSwaps
REAL                              :: RanNum
REAL                              :: Xitild,EtaTild
INTEGER                           :: p,q,SurfSideID,locBCID
INTEGER                           :: iCC
#ifndef IMPA
REAL                              :: n_loc(1:3)
#endif /*NOT IMPA*/
LOGICAL                           :: IsAuxBC
!===================================================================================================================================
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF
IF (IsAuxBC) THEN
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.LE.PartAuxBC%ProbOfSpeciesSwaps(AuxBCIdx)) THEN
    targetSpecies=-1 !dummy init value
    DO iSwaps=1,PartAuxBC%NbrOfSpeciesSwaps(AuxBCIdx)
      IF (PartSpecies(PartID).eq.PartAuxBC%SpeciesSwaps(1,iSwaps,AuxBCIdx)) &
        targetSpecies = PartAuxBC%SpeciesSwaps(2,iSwaps,AuxBCIdx)
    END DO
    !swap species
    IF (targetSpecies.ge.0) IsSpeciesSwap=.TRUE.
    IF (targetSpecies.eq.0) THEN !delete particle -> same as PartAuxBC%OpenBC
      IF(CalcPartBalance) THEN
        nPartOut(PartSpecies(PartID))=nPartOut(PartSpecies(PartID)) + 1
        PartEkinOut(PartSpecies(PartID))=PartEkinOut(PartSpecies(PartID))+CalcEkinPart(PartID)
      END IF ! CalcPartBalance
      !---- Counter for collisions (normal wall collisions - not to count if only Swaps to be counted, IsSpeciesSwap: already counted)
      PDM%ParticleInside(PartID) = .FALSE.
      alpha=-1.
    ELSEIF (targetSpecies.gt.0) THEN !swap species
      PartSpecies(PartID)=targetSpecies
    END IF
  END IF !RanNum.LE.PartAuxBC%ProbOfSpeciesSwaps
ELSE
#ifndef IMPA
IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
  ELSE
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
  END IF
END IF

IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
  RETURN
!  CALL Abort(&
!    __STAMP__&
!    ,'SpeciesSwap was called for an particle coming from outside of the mesh!')
END IF
#endif /*NOT IMPA*/

locBCID = PartBound%MapToPartBC(BC(SideID))
CALL RANDOM_NUMBER(RanNum)
IF(RanNum.LE.PartBound%ProbOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID)))) THEN
  targetSpecies=-1 !dummy init value
  DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID)))
  IF (PartSpecies(PartID).eq.PartBound%SpeciesSwaps(1,iSwaps,PartBound%MapToPartBC(BC(SideID)))) &
    targetSpecies = PartBound%SpeciesSwaps(2,iSwaps,PartBound%MapToPartBC(BC(SideID)))
  END DO
  !swap species
  IF (targetSpecies.ge.0) IsSpeciesSwap=.TRUE.
  IF ( (targetSpecies.eq.0) .OR. (.NOT.CalcSurfCollis%Only0Swaps) ) THEN
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      !---- Counter for swap species collisions
      SurfSideID=SurfMesh%SideIDToSurfID(SideID)
      ! compute p and q
      ! correction of xi and eta, can only be applied if xi & eta are not used later!
      IF (TriaTracking) THEN
        p=1 ; q=1
      ELSE
        Xitild =MIN(MAX(-1.,xi ),0.99)
        Etatild=MIN(MAX(-1.,eta),0.99)
        p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
        q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
      END IF

      SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q) = SampWall(SurfSideID)%State(12+PartSpecies(PartID),p,q) + 1
      IF (CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))) THEN
        AnalyzeSurfCollis%Number(PartSpecies(PartID)) = AnalyzeSurfCollis%Number(PartSpecies(PartID)) + 1
        AnalyzeSurfCollis%Number(nSpecies+1) = AnalyzeSurfCollis%Number(nSpecies+1) + 1
        IF (AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber) THEN
CALL Abort(&
__STAMP__&
,'maxSurfCollisNumber reached!')
        END IF
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) &
          = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) &
          = PartState(PartID,4)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) &
          = PartState(PartID,5)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) &
          = PartState(PartID,6)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) &
          = LastPartPos(PartID,1)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) &
          = LastPartPos(PartID,2)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) &
          = LastPartPos(PartID,3)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = locBCID
      END IF
    END IF
  END IF
  IF (targetSpecies.eq.0) THEN !delete particle -> same as PartBound%OpenBC
    DO iCC=1,nCollectChargesBCs !-chargeCollect
      IF (CollectCharges(iCC)%BC .EQ. PartBound%MapToPartBC(BC(SideID))) THEN
        CollectCharges(iCC)%NumOfNewRealCharges = CollectCharges(iCC)%NumOfNewRealCharges &
          + Species(PartSpecies(PartID))%ChargeIC &
          * Species(PartSpecies(PartID))%MacroParticleFactor
        EXIT
      END IF
    END DO
    IF(CalcPartBalance) THEN
      nPartOut(PartSpecies(PartID))=nPartOut(PartSpecies(PartID)) + 1
      PartEkinOut(PartSpecies(PartID))=PartEkinOut(PartSpecies(PartID))+CalcEkinPart(PartID)
    END IF ! CalcPartBalance
    ! sample values of deleted species
    IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      SurfSideID=SurfMesh%SideIDToSurfID(SideID)
      IF (TriaTracking) THEN
        p=1 ; q=1
      ELSE
        Xitild =MIN(MAX(-1.,xi ),0.99)
        Etatild=MIN(MAX(-1.,eta),0.99)
        p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
        q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
      END IF
      !----  Sampling Forces at walls
  !       SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
  !                                          * (v_old(1:3) - PartState(PartID,4:6)) * Species(PartSpecies(PartID))%MacroParticleFactor
      SampWall(SurfSideID)%State(10,p,q)= SampWall(SurfSideID)%State(10,p,q) + Species(PartSpecies(PartID))%MassIC &
                                            * ( PartState(PartID,4)) * Species(PartSpecies(PartID))%MacroParticleFactor
      SampWall(SurfSideID)%State(11,p,q)= SampWall(SurfSideID)%State(11,p,q) + Species(PartSpecies(PartID))%MassIC &
                                            * ( PartState(PartID,5)) * Species(PartSpecies(PartID))%MacroParticleFactor
      SampWall(SurfSideID)%State(12,p,q)= SampWall(SurfSideID)%State(12,p,q) + Species(PartSpecies(PartID))%MassIC &
                                            * ( PartState(PartID,6)) * Species(PartSpecies(PartID))%MacroParticleFactor
    END IF
    !---- Counter for collisions (normal wall collisions - not to count if only Swaps to be counted, IsSpeciesSwap: already counted)
    PDM%ParticleInside(PartID) = .FALSE.
    alpha=-1.
#ifdef IMPA
    DoPartInNewton(PartID) = .FALSE.
    PartIsImplicit(PartID) = .FALSE.
#endif /*IMPA*/
  ELSEIF (targetSpecies.gt.0) THEN !swap species
    DO iCC=1,nCollectChargesBCs !-chargeCollect
      IF (CollectCharges(iCC)%BC .EQ. PartBound%MapToPartBC(BC(SideID))) THEN
        CollectCharges(iCC)%NumOfNewRealCharges = CollectCharges(iCC)%NumOfNewRealCharges &
          + (Species(PartSpecies(PartID))%ChargeIC-Species(targetSpecies)%ChargeIC) &
          * Species(PartSpecies(PartID))%MacroParticleFactor
        EXIT
      END IF
    END DO
    PartSpecies(PartID)=targetSpecies
  END IF
END IF !RanNum.LE.PartBound%ProbOfSpeciesSwaps
END IF !IsAuxBC
END SUBROUTINE SpeciesSwap


SUBROUTINE PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,ElemID,BCSideID,opt_perimoved,TriNum)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking,DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:epsInCell,GEO,SidePeriodicType
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,PEM
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Particle_Mesh_Vars,     ONLY:PartSideToElem
#if defined(IMPA)
USE MOD_TimeDisc_Vars,          ONLY:ESDIRK_a,ERK_a
#endif /*IMPA */
#if defined(ROS)
USE MOD_TimeDisc_Vars,          ONLY:RK_A
#endif /*ROS */
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars, ONLY:PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID!,ElemID
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_perimoved
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3)
REAL                                 :: epsLength
INTEGER                              :: PVID,moved(2),locSideID
!===================================================================================================================================

!OneMinus=1.0-MAX(epsInCell,epsilontol)
epsLength=MAX(epsInCell,epsilontol)*lengthPartTrajectory

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
  ELSE
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
  END IF
END IF

IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
  IF(PRESENT(opt_perimoved)) opt_perimoved=.FALSE.
  RETURN
ELSE
  IF(PRESENT(opt_perimoved)) opt_perimoved=.TRUE.
END IF

PVID = SidePeriodicType(SideID)

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition: ',PartState(PartID,1:3)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPos:      ',LastPartPos(PartID,1:3)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! set last particle position on face
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
! perform the periodic movement
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
! update particle positon after periodic BC
!PartState(PartID,1:3)   = PartState(PartID,1:3) + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
PartState(PartID,1:3) = LastPartPos(PartID,1:3) + (lengthPartTrajectory-alpha)*PartTrajectory
lengthPartTrajectory  = lengthPartTrajectory - alpha


#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition-pp: ',PartState(PartID,1:3)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPo-pp:       ',LastPartPos(PartID,1:3)
  END IF
END IF
#endif /*CODE_ANALYZE*/

#if defined(ROS) || defined(IMPA)
PEM%PeriodicMoved(PartID)=.TRUE.
#endif

! refmapping and tracing
! move particle from old element to new element
locSideID = PartSideToElem(S2E_LOC_SIDE_ID,SideID)
Moved     = PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
ElemID    = Moved(1)
#ifdef MPI
IF(ElemID.EQ.-1)THEN
  CALL abort(&
__STAMP__&
,' Halo region to small. Neighbor element is missing!')
END IF
#endif /*MPI*/
!ElemID   =PEM%Element(PartID)
IF (DoRefMapping) PEM%LastElement(PartID) = 0

IF (DoRefMapping) PEM%LastElement(PartID) = 0

IF(1.EQ.2)THEN
  alpha=0.2
END IF

END SUBROUTINE PeriodicBC


SUBROUTINE SideAnalysis(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,flip,locSideID,ElemID &
  , IsSpeciesSwap,BCSideID &
  , opt_crossed&
  , TriNum)
!----------------------------------------------------------------------------------------------------------------------------------!
! Analyze particle crossing (inner) side
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,CalcSurfCollis,AnalyzeSurfCollis
USE MOD_Particle_Mesh_Vars,     ONLY:epsInCell
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,nSpecies,PartSpecies,WriteMacroSurfaceValues
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_DSMC_Vars,              ONLY:DSMC
!USE MOD_LD_Vars,                ONLY:useLD
USE MOD_TImeDisc_Vars,          ONLY:tend,time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID, flip,locSideID
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL       :: BCSideID
INTEGER,INTENT(IN),OPTIONAL       :: TriNum
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL      :: opt_crossed
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3), WallVelo(3)
REAL                                 :: epsLength
INTEGER                              :: locBCID
INTEGER                              :: moved(2)
!===================================================================================================================================

epsLength=MAX(epsInCell,epsilontol)*lengthPartTrajectory
WallVelo=PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
locBCID=PartBound%MapToPartBC(BC(SideID))

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=SideID)
  ELSE
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
  END IF
END IF

IF(DOT_PRODUCT(PartTrajectory,n_loc).LE.0.) THEN
  IF(PRESENT(opt_crossed)) opt_crossed=.FALSE.
  RETURN
ELSE
  IF(PRESENT(opt_crossed)) opt_crossed=.TRUE.
END IF


! Wall sampling Macrovalues
!IF((.NOT.Symmetry).AND.(.NOT.UseLD)) THEN !surface mesh is not build for the symmetry BC!?!
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    !---- Counter for collisions (normal wall collisions - not to count if only Swaps to be counted, IsSpeciesSwap: already counted)
    IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
      IF (CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))) THEN
        AnalyzeSurfCollis%Number(PartSpecies(PartID)) = AnalyzeSurfCollis%Number(PartSpecies(PartID)) + 1
        AnalyzeSurfCollis%Number(nSpecies+1) = AnalyzeSurfCollis%Number(nSpecies+1) + 1
        IF (AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber) THEN
CALL Abort(&
__STAMP__&
,'maxSurfCollisNumber reached!')
        END IF
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) &
          = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4) &
          = PartState(PartID,4)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5) &
          = PartState(PartID,5)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6) &
          = PartState(PartID,6)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7) &
          = LastPartPos(PartID,1)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8) &
          = LastPartPos(PartID,2)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9) &
          = LastPartPos(PartID,3)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) &
          = locBCID
      END IF
    END IF
  END IF
!END IF

! refmapping and tracing
! move particle from old element to new element
Moved     = PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
ElemID    = Moved(1)
#ifdef MPI
IF(ElemID.EQ.-1)THEN
  CALL abort(&
__STAMP__&
,' Mesh-connectivity broken or halo region to small. Neighbor element is missing!')
END IF
#endif /*MPI*/

END SUBROUTINE SideAnalysis


FUNCTION PARTSWITCHELEMENT(xi,eta,locSideID,SideID,ElemID)
!===================================================================================================================================
! particle moves through face and switches element
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToElemAndSide
USE MOD_Mesh_Vars,              ONLY:MortarType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: locSideID, SideID,ElemID
REAL,INTENT(IN)     :: xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: PARTSWITCHELEMENT
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! move particle to new element
!     Type 1               Type 2              Type3
!      eta                  eta                 eta
!       ^                    ^                   ^
!       |                    |                   |
!   +---+---+            +---+---+           +---+---+
!   | 3 | 4 |            |   2   |           |   |   |
!   +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!   | 1 | 2 |            |   1   |           |   |   |
!   +---+---+            +---+---+           +---+---+

SELECT CASE(MortarType(1,SideID))
CASE(1)
  IF(Xi.GT.0.)THEN
    IF(Eta.GT.0.)THEN
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(4  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(4+4,locSideID,ElemID)
    ELSE
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
    END IF
  ELSE
    IF(Eta.GT.0.)THEN
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(3  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(3+4,locSideID,ElemID)
    ELSE
      PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
      PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
    END IF
  END IF
CASE(2)
  IF(Eta.GT.0.)THEN
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
  ELSE
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
  END IF
CASE(3)
  IF(Xi.LE.0.)THEN
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
  ELSE
    PARTSWITCHELEMENT(1)=PartElemToElemAndSide(2  ,locSideID,ElemID)
    PARTSWITCHELEMENT(2)=PartElemToElemAndSide(2+4,locSideID,ElemID)
  END IF
CASE DEFAULT ! normal side OR small mortar side
  PARTSWITCHELEMENT(1)=PartElemToElemAndSide(1  ,locSideID,ElemID)
  PARTSWITCHELEMENT(2)=PartElemToElemAndSide(1+4,locSideID,ElemID)
END SELECT

END FUNCTION PARTSWITCHELEMENT


SUBROUTINE ReactiveSurfaceTreatment(PartTrajectory,LengthPartTrajectory,alpha,xi,eta,PartID,sideID_IN,flip,IsSpeciesSwap,&
                              adsindex,BCSideID,Opt_Reflected,TriNum)
!===================================================================================================================================
!> Routine for Selection of Surface interaction
!===================================================================================================================================
USE MOD_Globals                ,ONLY: CROSSNORM,UNITVECTOR
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_Part_Tools             ,ONLY: VELOFROMDISTRIBUTION
USE MOD_DSMC_Analyze           ,ONLY: CalcWallSample
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues
USE MOD_Particle_Vars          ,ONLY: PartState,Species,PartSpecies
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_Mesh_Vars              ,ONLY: BC,NGeo
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Tools,ONLY: PartEnergyToSurface,SurfaceToPartEnergy
USE MOD_Particle_Boundary_Tools,ONLY: TSURUTACONDENSCOEFF
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, dXiEQ_SurfSample, Partbound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time, dt, RKdtFrac
USE MOD_Particle_Surfaces_vars ,ONLY: SideNormVec,SideType,BezierControlPoints3D
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, ModelERSpecular, SurfModel
USE MOD_SMCR                   ,ONLY: SMCR_PartAdsorb
USE MOD_SEE                    ,ONLY: SEE_PartDesorb
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)       :: adsindex
REAL,INTENT(INOUT)          :: PartTrajectory(1:3), LengthPartTrajectory, alpha
REAL,INTENT(IN)             :: xi, eta
INTEGER,INTENT(IN)          :: PartID
INTEGER,INTENT(IN)          :: sideID_IN
INTEGER,INTENT(IN)          :: flip
LOGICAL,INTENT(IN)          :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL :: BCSideID
INTEGER,INTENT(IN),OPTIONAL :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL :: Opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: n_loc(1:3), tang1(1:3),tang2(1:3)
REAL                             :: Adsorption_prob, Recombination_prob
INTEGER                          :: interactionCase
INTEGER                          :: SurfSideID, SpecID
REAL                             :: Norm_velo!, Norm_Ec
INTEGER                          :: outSpec(2)
! variables for Energy sampling
REAL                             :: TransArray(1:6),IntArray(1:6), reactionEnthalpie
REAL                             :: oldVelo(1:3)
INTEGER                          :: locBCID
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
REAL                             :: TransACC!, VibACC, RotACC
! Polyatomic Molecules
INTEGER                          :: iReact, RecombReactID
REAL                             :: VeloCrad, Fak_D, NewVelo(3)
REAL                             :: Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                             :: POI_fak, TildTrajectory(3)
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
!===================================================================================================================================

! find normal vector two perpendicular tangential vectors (normal_vector points outwards !!!)
IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
    tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,BCSideID)-BezierControlPoints3D(:,0,0,BCSideID))
    tang2=CROSSNORM(n_loc,tang1)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,tang1=tang1,tang2=tang2,TriNum=TriNum,SideID=sideID_IN)
  ELSE
    SELECT CASE(SideType(sideID_IN))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,sideID_IN)
        tang1=UNITVECTOR(BezierControlPoints3D(:,NGeo,0,sideID_IN)-BezierControlPoints3D(:,0,0,sideID_IN))
        tang2=CROSSNORM(n_loc,tang1)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,sideID_IN)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,sideID_IN)
    END SELECT
    IF(flip.NE.0) n_loc=-n_loc
  END IF
END IF

! check if BC was already crossed
IF(DOT_PRODUCT(n_loc,PartTrajectory).LT.0.)  THEN
  IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
  RETURN
ELSE
  IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
END IF

! additional states
locBCID=PartBound%MapToPartBC(BC(sideID_IN))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)

! initialize sampling arrays
TransArray(:) = 0.0
IntArray(:) = 0.0

! compute p and q
! correction of xi and eta, can only be applied if xi & eta are not used later!
IF (TriaTracking) THEN
  p=1 ; q=1
ELSE
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
END IF

SurfSideID = SurfMesh%SideIDToSurfID(sideID_IN)
SpecID = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)
! Update wallcollision counter
SurfModel%Info(SpecID)%WallCollCount = SurfModel%Info(SpecID)%WallCollCount + 1
IF (PartBound%SurfaceModel(locBCID).EQ.1) THEN
  SurfModel%Info(SpecID)%Accomodation = SurfModel%Info(SpecID)%Accomodation &
      + (PartBound%TransACC(locBCID) + PartBound%VibACC(locBCID)+ PartBound%RotACC(locBCID))/3.
END IF
#endif

interactionCase = 0
reactionEnthalpie = 0. ! negative at evaporation and positive at condensation
SELECT CASE(PartBound%SurfaceModel(locBCID))
CASE (1)
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption_prob.GE.RanNum) .AND. &
     (Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
    outSpec(2) = SpecID
    interactionCase = 1
  END IF
CASE (2)
  ! Set probabilities
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  Recombination_prob = Adsorption%ProbDes(p,q,SurfSideID,SpecID)
  ! check if still enough saved particles on surface
  IF (Adsorption%Coverage(p,q,SurfSideID,SpecID).LE.(-SurfModel%SumAdsorbPart(p,q,SurfSideID,SpecID))) THEN
    Adsorption_prob = Adsorption_prob + Recombination_prob
    Recombination_prob = 0.
  END IF
  ! Decide what happens to colliding particle
  CALL RANDOM_NUMBER(RanNum)
  IF ((Adsorption_prob+Recombination_prob).GE.RanNum) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF ((Adsorption_prob/(Adsorption_prob+Recombination_prob)).GE.RanNum) THEN
      interactionCase = 1
      outSpec(1) = 0
      outSpec(2) = SpecID
    ELSE
      interactionCase = 3
      DO iReact = Adsorption%DissNum+1,(Adsorption%ReactNum)
        RecombReactID = iReact-Adsorption%DissNum
        IF (Adsorption%RecombReact(2,RecombReactID,SpecID).EQ.Adsorption%ResultSpec(locBCID,SpecID)) THEN
          EXIT
        END IF
      END DO
      outSpec(1) = Adsorption%RecombReact(1,RecombReactID,SpecID)
      outSpec(2) = Adsorption%RecombReact(2,RecombReactID,SpecID)
      reactionEnthalpie = - Adsorption%EDissBond(iReact,SpecID) * Adsorption%ReactAccomodation(locBCID,SpecID) * BoltzmannConst
    END IF
  END IF
CASE (3)
  Norm_velo = DOT_PRODUCT(PartState(PartID,4:6),n_loc(1:3))
  !Norm_Ec = 0.5 * Species(SpecID)%MassIC * Norm_velo**2 + PartStateIntEn(PartID,1) + PartStateIntEn(PartID,2)
  CALL SMCR_PartAdsorb(p,q,SurfSideID,PartID,Norm_velo,interactionCase,outSpec,reactionEnthalpie)
CASE (4)
  ! TODO
CASE (5,6) ! Copied from CASE(1) and adjusted for secondary e- emission (SEE)
           ! 5: SEE by Levko2015
           ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! Get electron emission probability
  CALL SEE_PartDesorb(PartBound%SurfaceModel(locBCID),PartID,Adsorption_prob,interactionCase,outSpec)
  !Adsorption_prob = 1. !Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  ! CALL RANDOM_NUMBER(RanNum)
  ! IF(Adsorption_prob.GE.RanNum)THEN
  !    !  .AND. &
  !    !(Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
  !   outSpec(1) = SpecID
  !   outSpec(2) = 4!SpecID ! electron
  !   interactionCase = -2 ! perfect elastic scattering + particle creation
  !   WRITE (*,*) "SpecID,outSpec(1),outSpec(2),interactionCase =", SpecID,outSpec(1),outSpec(2),interactionCase
  ! END IF
CASE (101) ! constant condensation coefficient
  outSpec(2) = SpecID
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption_prob.GE.RanNum) ) THEN
    interactionCase = 1
  END IF
CASE (102) ! calculate condensation probability by tsuruta2005 and reflection distribution function
  outSpec(2) = SpecID
  interactionCase = 4
  velocityDistribution='liquid_refl'
  Norm_velo = DOT_PRODUCT(PartState(PartID,4:6),n_loc(1:3))
  CALL RANDOM_NUMBER(RanNum)
  IF ( (TSURUTACONDENSCOEFF(SpecID,Norm_velo,WallTemp).GE.RanNum) ) THEN
    interactionCase = 1
  END IF
END SELECT

SELECT CASE(interactionCase)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-4) ! Remove bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-3) ! Remove bombarding particle + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = 1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-2) ! Perfect elastic scattering + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-1) ! Perfect elastic scattering
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) ! Molecular adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
#endif
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) ! dissociative adsorption (particle dissociates on adsorption)
!-----------------------------------------------------------------------------------------------------------------------------------
  SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) + 1
  SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    !----  Sampling of reactionEnthalpie
    SampWall(SurfSideID)%SurfModelState(4,p,q) = SampWall(SurfSideID)%SurfModelState(4,p,q) &
                                           + reactionEnthalpie * Species(SpecID)%MacroParticleFactor
    ! Sample reaction counter
    DO iReact = 1,Adsorption%DissNum
      IF (Adsorption%DissocReact(1,iReact,SpecID).EQ.outSpec(1) .AND. Adsorption%DissocReact(2,iReact,SpecID).EQ.outSpec(2))THEN
        SampWall(SurfSideID)%SurfModelReactCount(iReact,SpecID,p,q)=SampWall(SurfSideID)%SurfModelReactCount(iReact,SpecID,p,q) + 1
      END IF
    END DO
  END IF
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! Eley-Rideal reaction (reflecting particle and changes species at contact and reaction partner removed from surface)
!-----------------------------------------------------------------------------------------------------------------------------------
  !SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = SurfModel%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) = SurfModel%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) - 1
  adsindex = 2
  ! --------
  ! sampling and analyze stuff for heat flux and reaction rates
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(outSpec(1))%NumOfDes = SurfModel%Info(outSpec(1))%NumOfDes + 1
  SurfModel%Info(outSpec(2))%NumOfDes = SurfModel%Info(outSpec(2))%NumOfDes + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    ! Sample recombination reaction counter
    DO iReact = 1,Adsorption%RecombNum
      IF (Adsorption%RecombReact(2,iReact,SpecID).EQ.outSpec(2))THEN
        SampWall(SurfSideID)%SurfModelReactCount(Adsorption%DissNum+iReact,SpecID,p,q) = &
            SampWall(SurfSideID)%SurfModelReactCount(Adsorption%DissNum+iReact,SpecID,p,q) + 1
      END IF
    END DO
    !----  Sampling of reactionEnthalpie
    reactionEnthalpie = reactionEnthalpie * Adsorption%ReactAccomodation(locBCID,SpecID)
    SampWall(SurfSideID)%SurfModelState(3,p,q) = SampWall(SurfSideID)%SurfModelState(3,p,q) &
                                           + reactionEnthalpie * Species(SpecID)%MacroParticleFactor
  END IF
  ! --------
  ! reflect particle and change its species
  oldVelo(1:3) = PartState(PartID,4:6)
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  IF (ModelERSpecular) THEN
    ! perfect velocity reflection
    NewVelo(1:3) = oldVelo(1:3) - 2.*DOT_PRODUCT(oldVelo(1:3),n_loc)*n_loc
    ! mass changes, therefore velocity is scaled because impuls remains the same
    NewVelo(1:3) = NewVelo(1:3) * (Species(outSpec(2))%MassIC/Species(SpecID)%MassIC)
  ELSE
    ! diffuse reflection
    TransACC   = PartBound%TransACC(locBCID)
    !VibACC     = PartBound%VibACC(locBCID)
    !RotACC     = PartBound%RotACC(locBCID)
    CALL RANDOM_NUMBER(RanNum)
    VeloCrad    = SQRT(-LOG(RanNum))
    CALL RANDOM_NUMBER(RanNum)
    VeloCz      = SQRT(-LOG(RanNum))
    Fak_D       = VeloCrad**2 + VeloCz**2
    EtraWall    = BoltzmannConst * WallTemp * Fak_D
    VeloReal    = SQRT(DOT_PRODUCT(oldVelo,oldVelo))
    EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
    EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
    Cmr         = SQRT(2.0 * EtraNew / (Species(outSpec(2))%MassIC * Fak_D))
    CALL RANDOM_NUMBER(RanNum)
    Phi     = 2.0 * PI * RanNum
    VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
    VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
    VeloCz  = Cmr * VeloCz
    NewVelo = VeloCx*tang1-tang2*VeloCy-VeloCz*n_loc
  END IF
  ! intersection point with surface
  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3)
  !----  saving new particle velocity
  PartState(PartID,4:6)   = NewVelo(1:3) + WallVelo(1:3)

  ! recompute trajectory etc
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory
  ! set new species
  PartSpecies(PartID) = OutSpec(2)

  CALL SurfaceToPartEnergy(PartID,OutSpec(2),WallTemp,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID,emission_opt=.TRUE.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4) ! Distribution function reflection  (reflecting particle according to defined distribution function)
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 2
  ! sampling and analyze stuff for heat flux and reaction rates
#if (PP_TimeDiscMethod==42)
  SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  SurfModel%ProperInfo(SpecID)%HeatFlux(1) = SurfModel%ProperInfo(SpecID)%HeatFlux(1) &
     + reactionEnthalpie/BoltzmannConst
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%SurfModelState(5,p,q) = SampWall(SurfSideID)%SurfModelState(5,p,q) &
                                           + ReactionEnthalpie * Species(SpecID)%MacroParticleFactor
  END IF
  oldVelo(1:3) = PartState(PartID,4:6)
  CALL PartEnergyToSurface(PartID,SpecID,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  ! sample new velocity
  NewVelo(1:3) = VELOFROMDISTRIBUTION(velocityDistribution,SpecID,WallTemp)
  ! important: n_loc points outwards
  PartState(PartID,4:6) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

  ! intersection point with surface
  LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution

  ! recompute trajectory etc
  PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * PartState(PartID,4:6)
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
  PartTrajectory=PartTrajectory/lengthPartTrajectory

  CALL SurfaceToPartEnergy(PartID,OutSpec(2),WallTemp,Transarray,IntArray)
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID,emission_opt=.TRUE.)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT ! diffuse reflection
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 0
END SELECT

END SUBROUTINE ReactiveSurfaceTreatment


SUBROUTINE SurfaceFluxBasedBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,lengthPartTrajectory,flip,xi,eta)
!===================================================================================================================================
! Treatment of particles at the boundary if adaptive surface BCs or circular inflows based on the surface flux are present
! Circular Inflow: Particles are deleted if within (allows multiple surface flux inflows defined by circles on a single boundary)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PDM, Species, LastPartPos, PartSpecies
USE MOD_Particle_Boundary_Vars, ONLY:PartBound
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_Particle_Analyze,       ONLY:CalcEkinPart
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: iPart, SideID, flip
REAL,INTENT(IN)                     :: PartTrajectory(1:3),lengthPartTrajectory,xi,eta
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                  :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: point(1:2), intersectionPoint(1:3), radius, n_loc(1:3)
INTEGER                             :: iSpec, iSF
!===================================================================================================================================
IF (.NOT.TriaTracking) THEN
  ! Inserted particles are "pushed" inside the domain and registered as passing through the BC side. If they are very close to the
  ! boundary (first if) than the normal vector is compared with the trajectory. If the particle is entering the domain from outside
  ! it was inserted during surface flux and this routine shall not performed.
  IF(alpha/lengthPartTrajectory.LE.epsilontol) THEN
    ! Determining the normal vector of the side, always pointing outside the domain
    SELECT CASE(SideType(SideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,SideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=SideID)
    END SELECT
    ! If flip is not zero, the normal vector of the side is pointing in the opposite direction
    IF(flip.NE.0) n_loc=-n_loc
    ! Comparing the normal vector with the particle trajectory, if the dot product is less/equal zero, the particle trajectory is
    ! pointing inside the domain
    IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
  END IF
END IF

iSpec = PartSpecies(iPart)
DO iSF=1,Species(iSpec)%nSurfacefluxBCs
  IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.PartBound%MapToPartBC(BC(SideID))) THEN
    IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      intersectionPoint(1:3) = LastPartPos(iPart,1:3) + alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(2))-Species(iSpec)%Surfaceflux(iSF)%origin(1)
      point(2)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(3))-Species(iSpec)%Surfaceflux(iSF)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
        PDM%ParticleInside(iPart)=.FALSE.
        alpha=-1.
        IF(CalcPartBalance) THEN
          nPartOut(iSpec)=nPartOut(iSpec) + 1
          PartEkinOut(iSpec)=PartEkinOut(iSpec)+CalcEkinPart(iPart)
        END IF ! CalcPartBalance
        ! Counting the particles leaving the domain through the constant mass flow boundary (circular inflow on reflective)
        IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut = &
                                                                Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut + 1
      END IF
    END IF
  END IF
END DO

END SUBROUTINE SurfaceFluxBasedBoundaryTreatment

END MODULE MOD_Particle_Boundary_Condition
