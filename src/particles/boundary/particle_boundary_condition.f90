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
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,KeepWallParticles, PartSurfaceModel, UseCircularInflow, UseAdaptive, Species
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

  ! Particle has left the element over an open BC
  crossedBC    =.TRUE.

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
    ! Decide if liquid or solid
    IF (PartBound%SolidState(PartBound%MapToPartBC(BC(SideID)))) THEN
      IF ((PartSurfaceModel.EQ.0) .OR. (.NOT.PartBound%SolidReactive(PartBound%MapToPartBC(BC(SideID))))) THEN
      ! simple reflection (previously used wall interaction model, maxwellian scattering)
        CALL RANDOM_NUMBER(RanNum)
        IF((RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))).OR.PorousReflection) THEN
          ! perfectly reflecting, specular re-emission
          CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum, opt_LocalSide= locSideID,opt_ElemID=ElemID)
        ELSE
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum, opt_LocalSide= locSideID,opt_ElemID=ElemID)
        END IF
      ELSE IF ((PartSurfaceModel.GT.0) .AND. (PartBound%SolidReactive(PartBound%MapToPartBC(BC(SideID))))) THEN
      ! chemical surface interaction (adsorption)
        adsorbindex = 0
        ! Decide which interaction (reflection, reaction, adsorption)
        CALL ReactiveSurfaceTreatment(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,adsorbindex,opt_Reflected=crossedBC&
                                              ,TriNum=TriNum)
        ! assign right treatment
        SELECT CASE (adsorbindex)
        CASE(1,2)
          ! 1: adsorption (is either removed or set to be on surface)
          ! 2: Eley-Rideal reaction (particle is removed and product inserted in surface flux)
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
  !         CALL Particle_ER_Reflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap)
        CASE(0) ! inelastic reflection
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
        CASE(-1) ! elastic reflection
          CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
        CASE DEFAULT ! should not happen
          WRITE(*,*)'Boundary_PIC: Adsorption error.'
          CALL Abort(&
__STAMP__,&
'Boundary_Error: Adsorptionindex switched to unknown value.')
        END SELECT
      END IF
    ELSE
      IF (PartSpecies(iPart).EQ.PartBound%LiquidSpec(PartBound%MapToPartBC(BC(SideID)))) THEN
        CALL ParticleCondensationCase(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,adsorbindex,TriNum=TriNum)
          IF (adsorbindex.EQ.1) THEN ! condensation (particle is removed)
            IF(CalcPartBalance) THEN
              nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
              PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
            END IF ! CalcPartBalance
            PDM%ParticleInside(iPart) = .FALSE.
            alpha=-1.
#ifdef IMPA
            PartIsImplicit(iPart) = .FALSE.
#endif /*IMPA*/
          ELSE IF (adsorbindex.EQ.0) THEN ! inelastic reflection
            CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
              IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
          END IF
      ELSE
        CALL RANDOM_NUMBER(RanNum)
        IF(RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))) THEN
          ! perfectly reflecting, specular re-emission
          CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
        ELSE
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip, &
            IsSpeciesSwap,opt_Reflected=crossedBC,TriNum=TriNum)
        END IF
      END IF
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
                            ,opt_Symmetry=.TRUE.,opt_Reflected=crossedBC,TriNum=TriNum,opt_LocalSide= locSideID,opt_ElemID=ElemID)
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
USE MOD_Particle_Vars,          ONLY:PDM,PartSpecies,KeepWallParticles, PartSurfaceModel
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

  ! Particle has left the element over an open BC
  crossedBC    =.TRUE.

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
    ! Decide if liquid or solid
    IF (PartBound%SolidState(PartBound%MapToPartBC(BC(SideID)))) THEN
      BCSideID=PartBCSideList(SideID)
      IF ((PartSurfaceModel.EQ.0) .OR. (.NOT.PartBound%SolidReactive(PartBound%MapToPartBC(BC(SideID))))) THEN
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
      ELSE IF ((PartSurfaceModel.GT.0) .AND. (PartBound%SolidReactive(PartBound%MapToPartBC(BC(SideID))))) THEN
      ! chemical surface interaction (adsorption)
        adsorbindex = 0
        ! Decide which interaction (reflection, reaction, adsorption)
        CALL ReactiveSurfaceTreatment(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,adsorbindex&
                                ,BCSideID=BCSideID,opt_reflected=crossedBC)
        ! assign right treatment
        SELECT CASE (adsorbindex)
        CASE(1,2)
          ! 1: adsorption (is either removed or set to be on surface)
          ! 2: Eley-Rideal reaction (particle is removed and inserted product inserted in surface flux)
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
      IF (PartSpecies(iPart).EQ.PartBound%LiquidSpec(PartBound%MapToPartBC(BC(SideID)))) THEN
        CALL ParticleCondensationCase(PartTrajectory,alpha,xi,eta,iPart,SideID,IsSpeciesSwap,adsorbindex,BCSideID=BCSideID)
          IF (adsorbindex.EQ.1) THEN ! condensation (particle is removed)
            IF(CalcPartBalance) THEN
              nPartOut(PartSpecies(iPart))=nPartOut(PartSpecies(iPart)) + 1
              PartEkinOut(PartSpecies(iPart))=PartEkinOut(PartSpecies(iPart))+CalcEkinPart(iPart)
            END IF ! CalcPartBalance
            PDM%ParticleInside(iPart) = .FALSE.
#ifdef IMPA
            PartIsImplicit(iPart) = .FALSE.
#endif /*IMPA*/
            alpha=-1.
          ELSE IF (adsorbindex.EQ.0) THEN ! inelastic reflection
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap&
                                ,BCSideID=BCSideID,opt_reflected=crossedBC)
          END IF
      ELSE
        CALL RANDOM_NUMBER(RanNum)
        IF(RanNum.GE.PartBound%MomentumACC(PartBound%MapToPartBC(BC(SideID)))) THEN
          ! perfectly reflecting, specular re-emission
          CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap &
                                ,BCSideID=BCSideID,opt_reflected=crossedBC)
        ELSE
          CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,IsSpeciesSwap&
                                ,BCSideID=BCSideID,opt_reflected=crossedBC)
        END IF
      END IF
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
  opt_Symmetry,opt_Reflected,TriNum,AuxBCIdx, opt_LocalSide, opt_ElemID)
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
USE MOD_Particle_Vars,          ONLY:VarTimeStep, usevMPF, PartMPF
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,epsilontol
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_DSMC_Vars,              ONLY:DSMC, RadialWeighting
USE MOD_DSMC_Symmetry2D,        ONLY:CalcRadWeightMPF
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
INTEGER,INTENT(IN),OPTIONAL       :: opt_LocalSide, opt_ElemID
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
REAL                                 :: MacroParticleFactor, POI_Y
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
      CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,ElemID_opt=opt_ElemID,LocSideID_opt=opt_LocalSide)
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

    IF (VarTimeStep%UseVariableTimeStep) THEN
      ! Sampling of the time step at the wall to get the correct time sample duration for the force per area calculation
      SampWall(SurfSideID)%State(12+nSpecies+1,p,q) = SampWall(SurfSideID)%State(12+nSpecies+1,p,q) &
                                                      + VarTimeStep%ParticleTimeStep(PartID)
    END IF

    IF(RadialWeighting%DoRadialWeighting) THEN
      POI_Y = LastPartPos(PartID,2) + PartTrajectory(2)*alpha
      IF (VarTimeStep%UseVariableTimeStep) THEN
        MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))*VarTimeStep%ParticleTimeStep(PartID)
      ELSE
        MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))
      END IF
    ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
      MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor*VarTimeStep%ParticleTimeStep(PartID)
    ELSE IF (usevMPF) THEN
      MacroParticleFactor = PartMPF(PartID)
    ELSE
      MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor
    END IF
  !----  Sampling Forces at walls
    SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
                                        * (v_old(1:3) - PartState(PartID,4:6)) * MacroParticleFactor
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
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
        !-- caution: for consistency with diffuse refl. v_old is used!
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = v_old(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(PartID,1:3)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
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
  ,opt_Reflected,TriNum,AuxBCIdx, opt_LocalSide, opt_ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
! only implemented for DoRefMapping tracking
! PartBCs are reduced!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals,                ONLY:CROSSNORM,abort,UNITVECTOR
USE MOD_Globals_Vars,           ONLY:PI, BoltzmannConst
USE MOD_Particle_Tracking_Vars, ONLY:TriaTracking, TrackInfo
USE MOD_Particle_Boundary_Vars, ONLY:PartBound,SurfMesh,SampWall,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars, ONLY:dXiEQ_SurfSample
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,Species,PartSpecies,nSpecies,WriteMacroSurfaceValues,Symmetry2D
USE MOD_Particle_Vars,          ONLY:Symmetry2DAxisymmetric, VarTimeStep, usevMPF, PartMPF
#if defined(LSERK) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars,          ONLY:PDM
#endif
USE MOD_Particle_Surfaces_vars, ONLY:SideNormVec,SideType,BezierControlPoints3D
USE MOD_Mesh_Vars,              ONLY:BC,NGEO
USE MOD_DSMC_Vars,              ONLY:SpecDSMC,CollisMode
USE MOD_DSMC_Vars,              ONLY:PartStateIntEn,DSMC, useDSMC, RadialWeighting
USE MOD_DSMC_Vars,              ONLY:PolyatomMolDSMC, VibQuantsPar
USE MOD_TimeDisc_Vars,          ONLY:dt,tend,time,RKdtFrac
USE MOD_Particle_Boundary_Vars, ONLY:AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
USE MOD_DSMC_Symmetry2D,        ONLY: CalcRadWeightMPF
USE MOD_Particle_Mesh_Vars,     ONLY: GEO, PartSideToElem
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
INTEGER,INTENT(IN),OPTIONAL       :: opt_LocalSide, opt_ElemID
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
! Symmetry
REAL                                :: rotVelY, rotVelZ, rotPosY, MacroParticleFactor, adaptTimeStep
REAL                                :: VelX, VelY, VelZ,VecX, VecY, VecZ
REAL                                :: Vector1(1:3), Vector2(1:3)
REAL                                :: POI_Y
REAL                                :: nx, ny, nz, nVal
INTEGER                             :: LocSideID, ElemID
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
          TriNum=TriNum,ElemID_opt=opt_ElemID,LocSideID_opt=opt_LocalSide)
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

IF(Symmetry2DAxisymmetric) THEN
  ! Storing the old and the new particle position (which is outside the domain), at this point the position is only in the xy-plane
  VelX = PartState(PartID,1) - LastPartPos(PartID,1)
  VelY = PartState(PartID,2) - LastPartPos(PartID,2)
  VelZ = PartState(PartID,3) - LastPartPos(PartID,3)

  ElemID = PartSideToElem(S2E_ELEM_ID,SideID)
  IF (ElemID .EQ. TrackInfo%CurrElem) THEN
    LocSideID = PartSideToElem(S2E_LOC_SIDE_ID,SideID)
  ELSE
    ElemID = PartSideToElem(S2E_NB_ELEM_ID,SideID)
    LocSideID = PartSideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  END IF

  ! Getting the vectors, which span the cell (1-2 and 1-4)
  Vector1(1:3)=GEO%NodeCoords(1:3,GEO%ElemSideNodeID(2,LocSideID,ElemID))-GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1,LocSideID,ElemID))
  Vector2(1:3)=GEO%NodeCoords(1:3,GEO%ElemSideNodeID(4,LocSideID,ElemID))-GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1,LocSideID,ElemID))

  ! Get the vector, which does NOT have the z-component
  IF (ABS(Vector1(3)).GT.ABS(Vector2(3))) THEN
    Vector1 = Vector2
  END IF
  ! Cross product of the two vectors is simplified as Vector1(3) is zero
  nx = Vector1(2)
  ny = -Vector1(1)
  nz = 0.
  ! Check for the correct orientation of the normal vectors (should be inwards)
  IF ((VelX*nx+VelY*ny).GT.0) THEN
    nx = -Vector1(2)
    ny = Vector1(1)
  END IF

  nVal = SQRT(nx*nx + ny*ny)
  nx = nx/nVal
  ny = ny/nVal
END IF

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

IF(RadialWeighting%DoRadialWeighting) THEN
  POI_Y = LastPartPos(PartID,2) + PartTrajectory(2)*alpha
  IF (VarTimeStep%UseVariableTimeStep) THEN
    MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))*VarTimeStep%ParticleTimeStep(PartID)
  ELSE
    MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))
  END IF
ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
  MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor*VarTimeStep%ParticleTimeStep(PartID)
ELSE IF (usevMPF) THEN
  MacroParticleFactor = PartMPF(PartID)
ELSE
  MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor
END IF

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

    SampWall(SurfSideID)%State(1,p,q) = SampWall(SurfSideID)%State(1,p,q) + EtraOld * MacroParticleFactor
    SampWall(SurfSideID)%State(2,p,q) = SampWall(SurfSideID)%State(2,p,q) + EtraWall * MacroParticleFactor
    SampWall(SurfSideID)%State(3,p,q) = SampWall(SurfSideID)%State(3,p,q) + EtraNew * MacroParticleFactor
  END IF
END IF !.NOT.IsAuxBC

!   Transformation local distribution -> global coordinates
! v = nv*u+t1*v+t2*f3

!NewVelo = VeloCx*tang1+CROSS(-n_loc,tang1)*VeloCy-VeloCz*n_loc
!---- Transformation local distribution -> global coordinates
IF(Symmetry2DAxisymmetric) THEN
  VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  VecZ = 0.
  NewVelo(1) = VecX*VeloCx + nx*VeloCz
  NewVelo(2) = VecY*VeloCx + ny*VeloCz
  NewVelo(3) = VeloCy
ELSE
  NewVelo(1:3) = VeloCx*tang1(1:3)-tang2(1:3)*VeloCy-VeloCz*n_loc(1:3)
END IF

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
      SampWall(SurfSideID)%State(4,p,q)=SampWall(SurfSideID)%State(4,p,q)+PartStateIntEn(PartID,2) * MacroParticleFactor
      SampWall(SurfSideID)%State(5,p,q)=SampWall(SurfSideID)%State(5,p,q)+ErotWall * MacroParticleFactor
      SampWall(SurfSideID)%State(6,p,q)=SampWall(SurfSideID)%State(6,p,q)+ErotNew * MacroParticleFactor
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
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * MacroParticleFactor
            SampWall(SurfSideID)%State(8,p,q)= SampWall(SurfSideID)%State(8,p,q) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * MacroParticleFactor
          END DO
        ELSE
          SampWall(SurfSideID)%State(7,p,q)= SampWall(SurfSideID)%State(7,p,q) + (VibQuant + DSMC%GammaQuant) &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * MacroParticleFactor
          SampWall(SurfSideID)%State(8,p,q)= SampWall(SurfSideID)%State(8,p,q) + VibQuantWall &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * MacroParticleFactor
        END IF
        SampWall(SurfSideID)%State(9,p,q)= SampWall(SurfSideID)%State(9,p,q) + EvibNew * MacroParticleFactor
! #endif
      END IF
      IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) VibQuantsPar(PartID)%Quants(:) = VibQuantTemp(:)
      PartStateIntEn(PartID,1) = EvibNew
#if ((PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==300))
    END IF
#endif
    END IF
    END IF ! CollisMode > 1
    END IF ! useDSMC



! Sampling of the time step at the wall to get the correct time sample duration for the surface values
IF (VarTimeStep%UseVariableTimeStep) THEN
  adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    SampWall(SurfSideID)%State(12+nSpecies+1,p,q) = SampWall(SurfSideID)%State(12+nSpecies+1,p,q) + adaptTimeStep
  END IF
ELSE
  adaptTimeStep = 1.
END IF

! intersection point with surface
LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha

! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
!TildPos       =PartState(PartID,1:3)-dt*RKdtFrac*PartState(PartID,4:6)
TildTrajectory=dt*RKdtFrac*PartState(PartID,4:6)*adaptTimeStep
POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
! travel rest of particle vector
!PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac * NewVelo(1:3)
IF (IsAuxBC) THEN
  IF (PartAuxBC%Resample(AuxBCIdx)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
ELSE
  IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
END IF
PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3) * adaptTimeStep

IF(Symmetry2DAxisymmetric) THEN
  ! Symmetry considerations --------------------------------------------------------
  rotPosY = SQRT(PartState(PartID,2)**2 + (PartState(PartID,3))**2)
  ! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
  !           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
  ! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
  ! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
  rotVelY = (NewVelo(2)*(PartState(PartID,2))+NewVelo(3)*PartState(PartID,3))/rotPosY
  rotVelZ = (-NewVelo(2)*PartState(PartID,3)+NewVelo(3)*(PartState(PartID,2)))/rotPosY

  PartState(PartID,2) = rotPosY
  PartState(PartID,3) = 0.0
  NewVelo(2) = rotVelY
  NewVelo(3) = rotVelZ
END IF

IF(Symmetry2D) THEN
  ! z-Variable is set to zero (should be for the axisymmetric case anyway after rotation)
  lastPartPos(PartID,3) = 0.0
  PartState(PartID,3)   = 0.0
END IF

IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
!----  Sampling force at walls
  SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) &
      + Species(PartSpecies(PartID))%MassIC * (PartState(PartID,4:6) - NewVelo(1:3)) * MacroParticleFactor
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
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = PartState(PartID,4:6)
      AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(PartID,1:3)
      AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
      AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
    END IF
  END IF
END IF
END IF !.NOT.IsAuxBC

!----  saving new particle velocity
PartState(PartID,4:6)   = NewVelo(1:3) + WallVelo(1:3)

! recompute trajectory etc
IF(Symmetry2DAxisymmetric) THEN
  PartTrajectory(1:2)=PartState(PartID,1:2) - LastPartPos(PartID,1:2)
  PartTrajectory(3) = 0.
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                          +PartTrajectory(2)*PartTrajectory(2))
  PartTrajectory=PartTrajectory/lengthPartTrajectory
ELSE
  PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                          +PartTrajectory(2)*PartTrajectory(2) &
                          +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
END IF

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
USE MOD_Particle_Vars,          ONLY:PartState,LastPartPos,PartSpecies,PDM,VarTimeStep, usevMPF, PartMPF
USE MOD_Particle_Vars,          ONLY:WriteMacroSurfaceValues,nSpecies,CollectCharges,nCollectChargesBCs,Species
USE MOD_Particle_Surfaces,      ONLY:CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Analyze_Vars,  ONLY:CalcPartBalance,nPartOut,PartEkinOut
USE MOD_Particle_Analyze,       ONLY: CalcEkinPart
USE MOD_Mesh_Vars,              ONLY:BC
USE MOD_DSMC_Vars,              ONLY:DSMC, RadialWeighting
USE MOD_DSMC_Symmetry2D,        ONLY:CalcRadWeightMPF
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
REAL                              :: MacroParticleFactor, POI_Y
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
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(PartID,1:3) + alpha * PartTrajectory(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = PartState(PartID,4:6)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(PartID,1:3)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
      END IF
    END IF
  END IF
  IF (targetSpecies.eq.0) THEN !delete particle -> same as PartBound%OpenBC
    IF(RadialWeighting%DoRadialWeighting) THEN
      POI_Y = LastPartPos(PartID,2) + PartTrajectory(2)*alpha
      IF (VarTimeStep%UseVariableTimeStep) THEN
        MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))*VarTimeStep%ParticleTimeStep(PartID)
      ELSE
        MacroParticleFactor = CalcRadWeightMPF(POI_Y,PartSpecies(PartID))
      END IF
    ELSE IF (VarTimeStep%UseVariableTimeStep) THEN
      MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor*VarTimeStep%ParticleTimeStep(PartID)
    ELSE IF (usevMPF) THEN
      MacroParticleFactor = PartMPF(PartID)
    ELSE
      MacroParticleFactor = Species(PartSpecies(PartID))%MacroParticleFactor
    END IF
    DO iCC=1,nCollectChargesBCs !-chargeCollect
      IF (CollectCharges(iCC)%BC .EQ. PartBound%MapToPartBC(BC(SideID))) THEN
        CollectCharges(iCC)%NumOfNewRealCharges = CollectCharges(iCC)%NumOfNewRealCharges &
          + Species(PartSpecies(PartID))%ChargeIC * MacroParticleFactor
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
      SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
                                            * PartState(PartID,4:6) * MacroParticleFactor
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
          + (Species(PartSpecies(PartID))%ChargeIC-Species(targetSpecies)%ChargeIC) * MacroParticleFactor
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
#if USE_MPI
IF(ElemID.EQ.-1)THEN
  CALL abort(&
__STAMP__&
,' Halo region to small. Neighbor element is missing!')
END IF
#endif /*USE_MPI*/
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
#if USE_MPI
IF(ElemID.EQ.-1)THEN
  CALL abort(&
__STAMP__&
,' Mesh-connectivity broken or halo region to small. Neighbor element is missing!')
END IF
#endif /*USE_MPI*/

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


SUBROUTINE ReactiveSurfaceTreatment(PartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,adsindex,BCSideID,Opt_Reflected,TriNum)
!===================================================================================================================================
! Routine for Selection of Surface interaction
!===================================================================================================================================
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_DSMC_Analyze           ,ONLY: CalcWallSample
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues, KeepWallParticles, PartSurfaceModel
USE MOD_Particle_Vars          ,ONLY: PartState,Species,PartSpecies
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
!USE MOD_Particle_Vars          ,ONLY : PDM, LastPartPos
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: CollisMode, PolyatomMolDSMC
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, dXiEQ_SurfSample, Partbound, SampWall
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
USE MOD_Particle_Surfaces_vars ,ONLY: SideNormVec,SideType
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_SMCR                   ,ONLY: SMCR_PartAdsorb
USE MOD_SEE                    ,ONLY: SEE_PartDesorb
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)            :: adsindex
REAL,INTENT(INOUT)               :: PartTrajectory(1:3), alpha
REAL,INTENT(IN)                  :: xi, eta
INTEGER,INTENT(IN)               :: PartID, GlobSideID
LOGICAL,INTENT(IN)               :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL      :: BCSideID
INTEGER,INTENT(IN),OPTIONAL      :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT),OPTIONAL     :: Opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: n_loc(1:3), tang1(1:3),tang2(1:3)
REAL                             :: Adsorption_prob, Recombination_prob
INTEGER                          :: adsorption_case
INTEGER                          :: SurfSideID, SpecID
REAL                             :: Norm_velo!, Norm_Ec
INTEGER                          :: outSpec(2)
! variables for Energy sampling
!   REAL                             :: IntersectionPos(1:3)
REAL                             :: TransArray(1:6),IntArray(1:6), AdsorptionEnthalpie
REAL                             :: VelXold, VelYold, VelZold
INTEGER                          :: locBCID, VibQuant, VibQuantWall
!   INTEGER                          :: VibQuantNew
!   REAL                             :: VibQuantNewR
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
!   REAL                             :: TransACC, VibACC, RotACC
REAL                             :: ErotNew, ErotWall, EVibNew
! Polyatomic Molecules
REAL, ALLOCATABLE                :: RanNumPoly(:)
INTEGER                          :: iPolyatMole, iDOF
INTEGER, ALLOCATABLE             :: VibQuantWallPoly(:)
!   REAL, ALLOCATABLE                :: VibQuantNewRPoly(:)
!   INTEGER, ALLOCATABLE             :: VibQuantNewPoly(:), VibQuantTemp(:)
INTEGER                          :: iReact
!===================================================================================================================================

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=GlobSideID)
  ELSE
    SELECT CASE(SideType(GlobSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,GlobSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,GlobSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,GlobSideID)
    END SELECT
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
locBCID=PartBound%MapToPartBC(BC(GlobSideID))
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

SurfSideID = SurfMesh%SideIDToSurfID(GlobSideID)
SpecID = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)
! Update wallcollision counter
Adsorption%AdsorpInfo(SpecID)%WallCollCount = Adsorption%AdsorpInfo(SpecID)%WallCollCount + 1
IF (PartSurfaceModel.EQ.1) THEN
  Adsorption%AdsorpInfo(SpecID)%Accomodation = Adsorption%AdsorpInfo(SpecID)%Accomodation &
      + (PartBound%TransACC(locBCID) + PartBound%VibACC(locBCID)+ PartBound%RotACC(locBCID))/3.
END IF
#endif

adsorption_case = 0
AdsorptionEnthalpie = 0.
SELECT CASE(PartSurfaceModel)
CASE (1)
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption_prob.GE.RanNum) .AND. &
     (Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
    outSpec(1) = SpecID
    adsorption_case = 1
  END IF
CASE (2)
  ! Set probabilities
  Adsorption_prob = Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  Recombination_prob = Adsorption%ProbDes(p,q,SurfSideID,SpecID)
  ! check if still enough saved particles on surface
  IF (Adsorption%Coverage(p,q,SurfSideID,SpecID).LE.(-Adsorption%SumAdsorbPart(p,q,SurfSideID,SpecID))) THEN
    Adsorption_prob = Adsorption_prob + Recombination_prob
    Recombination_prob = 0.
  END IF
  ! Decide what happens to colliding particle
  CALL RANDOM_NUMBER(RanNum)
  IF ((Adsorption_prob+Recombination_prob).GE.RanNum) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF ((Adsorption_prob/(Adsorption_prob+Recombination_prob)).GE.RanNum) THEN
      adsorption_case = 1
      outSpec(1) = SpecID
      outSpec(2) = 0
    ELSE
      adsorption_case = 3
      outSpec(1) = Adsorption%RecombData(1,SpecID)
      outSpec(2) = Adsorption%RecombData(2,SpecID)
      AdsorptionEnthalpie = - Adsorption%RecombEnergy(locBCID,SpecID) * Adsorption%RecombAccomodation(locBCID,SpecID) &
                          * BoltzmannConst
    END IF
  END IF
CASE (3)
  Norm_velo = PartState(PartID,4)*n_loc(1) + PartState(PartID,5)*n_loc(2) + PartState(PartID,6)*n_loc(3)
  !Norm_Ec = 0.5 * Species(SpecID)%MassIC * Norm_velo**2 + PartStateIntEn(PartID,1) + PartStateIntEn(PartID,2)
  CALL SMCR_PartAdsorb(p,q,SurfSideID,PartID,Norm_velo,adsorption_case,outSpec,AdsorptionEnthalpie)
CASE (4)
  ! TODO
CASE (5,6) ! Copied from CASE(1) and adjusted for secondary e- emission (SEE)
           ! 5: SEE by Levko2015
           ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! Get electron emission probability
  CALL SEE_PartDesorb(PartSurfaceModel,PartID,Adsorption_prob,adsorption_case,outSpec)
  !Adsorption_prob = 1. !Adsorption%ProbAds(p,q,SurfSideID,SpecID)
  ! CALL RANDOM_NUMBER(RanNum)
  ! IF(Adsorption_prob.GE.RanNum)THEN
  !    !  .AND. &
  !    !(Adsorption%Coverage(p,q,SurfSideID,SpecID).LT.Adsorption%MaxCoverage(SurfSideID,SpecID)) ) THEN
  !   outSpec(1) = SpecID
  !   outSpec(2) = 4!SpecID ! electron
  !   adsorption_case = -2 ! perfect elastic scattering + particle creation
  !   WRITE (*,*) "SpecID,outSpec(1),outSpec(2),adsorption_case =", SpecID,outSpec(1),outSpec(2),adsorption_case
  ! END IF
END SELECT

SELECT CASE(adsorption_case)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-4) ! Remove bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 2
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-3) ! Remove bombarding particle + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = 2
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-2) ! Perfect elastic scattering + electron creation
!-----------------------------------------------------------------------------------------------------------------------------------
  Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-1) ! Perfect elastic scattering
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) ! Molecular adsorption
!-----------------------------------------------------------------------------------------------------------------------------------
  Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) = Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(SpecID)%NumOfAds = Adsorption%AdsorpInfo(SpecID)%NumOfAds + 1
#endif
  IF ((KeepWallparticles).OR.&
      (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
!     ! Allocate particle belonging adsorbing side-index and side-subsurface-indexes
!     IF (KeepWallParticles) THEN
!       PDM%PartAdsorbSideIndx(1,PartID) = GlobSideID
!       PDM%PartAdsorbSideIndx(2,PartID) = p
!       PDM%PartAdsorbSideIndx(3,PartID) = q
!         LastPartPos(PartID,1) = PartState(PartID,1)
!         LastPartPos(PartID,2) = PartState(PartID,2)
!         LastPartPos(PartID,3) = PartState(PartID,3)
!     END IF
    VelXold = PartState(PartID,4)
    VelYold = PartState(PartID,5)
    VelZold = PartState(PartID,6)
    PartState(PartID,4)  = WallVelo(1)
    PartState(PartID,5)  = WallVelo(2)
    PartState(PartID,6)  = WallVelo(3)

    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(outSpec(1))%MassIC * VeloReal**2
    EtraWall = 0.0
    EtraNew = EtraWall

    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    ! must be old_velocity-new_velocity
    TransArray(4) = VelXold-PartState(PartID,4)
    TransArray(5) = VelYold-PartState(PartID,5)
    TransArray(6) = VelZold-PartState(PartID,6)

    !---- Internal energy accommodation
    IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(outSpec(1))%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = 0
      ErotNew  = 0
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      IF(SpecDSMC(outSpec(1))%PolyatomicMol) THEN
        EvibNew = 0.0
        iPolyatMole = SpecDSMC(outSpec(1))%SpecToPolyArray
        ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))!, &
!                   VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                   VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        END DO
!           VibQuantNewRPoly(:) = VibQuantWallPoly(:)
!           VibQuantNewPoly = INT(VibQuantNewRPoly)
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!             CALL RANDOM_NUMBER(RanNum)
!             IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!             ELSE
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!             END IF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(outSpec(1))%MacroParticleFactor
          IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
          IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib) &
                    - DSMC%GammaQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(1))%MaxVibQuant)
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        END DO
!           VibQuantNewR = VibQuantWall
!           VibQuantNew = INT(VibQuantNewR)
!           CALL RANDOM_NUMBER(RanNum)
!           IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!             EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!           ELSE
!             EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!           END IF
        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
      END IF
      SDEALLOCATE(RanNumPoly)
      SDEALLOCATE(VibQuantWallPoly)
!         SDEALLOCATE(VibQuantNewRPoly)
!         SDEALLOCATE(VibQuantNewPoly)
!         SDEALLOCATE(VibQuantTemp)
!         IntArray(6) = EvibNew
    END IF
    END IF
    !End internal energy accomodation

    !----  Sampling of energies
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,AdsorptionEnthalpie&
        ,locBCID)
    END IF
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) ! dissociative adsorption (particle dissociates on adsorption)
!-----------------------------------------------------------------------------------------------------------------------------------
  Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) = Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) + 1
  Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) = Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(2)) + 1
  adsindex = 1
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(SpecID)%NumOfAds = Adsorption%AdsorpInfo(SpecID)%NumOfAds + 1
#endif
  IF ((KeepWallparticles).OR.&
      (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    VelXold = PartState(PartID,4)
    VelYold = PartState(PartID,5)
    VelZold = PartState(PartID,6)
    PartState(PartID,4)  = WallVelo(1)
    PartState(PartID,5)  = WallVelo(2)
    PartState(PartID,6)  = WallVelo(3)

    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(SpecID)%MassIC * VeloReal**2
    EtraWall = 0.0
    EtraNew = EtraWall

    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    ! must be old_velocity-new_velocity
    TransArray(4) = VelXold-PartState(PartID,4)
    TransArray(5) = VelYold-PartState(PartID,5)
    TransArray(6) = VelZold-PartState(PartID,6)

    !---- Internal energy accommodation
    IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(SpecID)%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = 0 !- BoltzmannConst * WallTemp * LOG(RanNum)
      ErotNew  = 0 !PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      ! first add zero point energy of reactant to intarrray
      IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                        * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(SpecID)%MacroParticleFactor
        END DO
      ELSE
        VibQuant = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) &
                    - DSMC%GammaQuant)
        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(SpecID)%CharaTVib
      END IF
      ! then do first product
      EvibNew = 0.0
      IF (SpecDSMC(outSpec(1))%InterID.EQ.2) THEN
        IF(SpecDSMC(outSpec(1))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(outSpec(1))%SpecToPolyArray
!             ALLOCATE(VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))!, &
!                     VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                     VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
            CALL RANDOM_NUMBER(RanNumPoly)
            VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          END DO
!             VibQuantNewRPoly(:) = VibQuantWallPoly(:)
!             VibQuantNewPoly = INT(VibQuantNewRPoly)
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!               CALL RANDOM_NUMBER(RanNum)
!               IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!                 EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!                           * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!                 VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!               ELSE
!                 EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!                           * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!                 VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!               END IF
            IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                        * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
            IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
          END DO
        ELSE
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
          DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(1))%MaxVibQuant)
            CALL RANDOM_NUMBER(RanNum)
            VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
          END DO
!             VibQuantNewR = VibQuantWall
!             VibQuantNew = INT(VibQuantNewR)
!             CALL RANDOM_NUMBER(RanNum)
!             IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!               EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!             ELSE
!               EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!             END IF
          IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
          IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        END IF
      END IF
      ! then do second product
      IF (SpecDSMC(outSpec(2))%InterID.EQ.2) THEN
        IF(SpecDSMC(outSpec(2))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(outSpec(2))%SpecToPolyArray
!           ALLOCATE(VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))!, &
!                   VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                   VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
            CALL RANDOM_NUMBER(RanNumPoly)
            VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
          END DO
!           VibQuantNewRPoly(:) = VibQuantWallPoly(:)
!           VibQuantNewPoly = INT(VibQuantNewRPoly)
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!             CALL RANDOM_NUMBER(RanNum)
!             IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!             ELSE
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!             END IF
            IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                        * SpecDSMC(outSpec(2))%CharaTVib * Species(outSpec(2))%MacroParticleFactor
            IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(outSpec(2))%CharaTVib * Species(outSpec(2))%MacroParticleFactor
          END DO
        ELSE
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(2))%CharaTVib)
          DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(2))%MaxVibQuant)
            CALL RANDOM_NUMBER(RanNum)
            VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(2))%CharaTVib)
          END DO
!           VibQuantNewR = VibQuantWall
!           VibQuantNew = INT(VibQuantNewR)
!           CALL RANDOM_NUMBER(RanNum)
!           IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!             EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(outSpec(2))%CharaTVib
!           ELSE
!             EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(outSpec(2))%CharaTVib
!           END IF
          IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(2))%CharaTVib
          IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(2))%CharaTVib
        END IF
      END IF
      SDEALLOCATE(RanNumPoly)
      SDEALLOCATE(VibQuantWallPoly)
!         SDEALLOCATE(VibQuantNewRPoly)
!         SDEALLOCATE(VibQuantNewPoly)
!         SDEALLOCATE(VibQuantTemp)
      ! assign new vibrational energy
!         IntArray(6) = EvibNew
    END IF
    END IF
    !End internal energy accomodation

    !----  Sampling of energies
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,AdsorptionEnthalpie&
        ,locBCID)
    END IF
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! Eley-Rideal reaction (reflecting particle species change at contact and reaction partner removed from surface)
!-----------------------------------------------------------------------------------------------------------------------------------
  Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) = Adsorption%SumERDesorbed(p,q,SurfSideID,outSpec(2)) + 1
  Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) = Adsorption%SumAdsorbPart(p,q,SurfSideID,outSpec(1)) - 1
  adsindex = 2
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(outSpec(1))%NumOfDes = Adsorption%AdsorpInfo(outSpec(1))%NumOfDes + 1
  Adsorption%AdsorpInfo(outSpec(2))%NumOfDes = Adsorption%AdsorpInfo(outSpec(2))%NumOfDes + 1
#endif
  ! Sample recombination coefficient
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    IF (PartSurfaceModel.EQ.2) THEN
      DO iReact = 1,Adsorption%RecombNum
        IF (Adsorption%RecombData(2,SpecID).EQ.outSpec(2))THEN
          SampWall(SurfSideID)%Reaction(1,SpecID,p,q) = SampWall(SurfSideID)%Reaction(1,SpecID,p,q) + 1
        END IF
      END DO
    ELSE IF ( PartSurfaceModel.EQ.3) THEN
      DO iReact = 1,Adsorption%RecombNum
        IF (Adsorption%AssocReact(2,iReact,SpecID).EQ.outSpec(2))THEN
          SampWall(SurfSideID)%Reaction(iReact,SpecID,p,q) = SampWall(SurfSideID)%Reaction(iReact,SpecID,p,q) + 1
        END IF
      END DO
    END IF
  END IF

  IF ((KeepWallparticles).OR.&
      (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    VelXold = PartState(PartID,4)
    VelYold = PartState(PartID,5)
    VelZold = PartState(PartID,6)
    PartState(PartID,4)  = WallVelo(1)
    PartState(PartID,5)  = WallVelo(2)
    PartState(PartID,6)  = WallVelo(3)

    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(outSpec(1))%MassIC * VeloReal**2
    EtraWall = 0.0
    EtraNew = EtraWall

    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    ! must be old_velocity-new_velocity
    TransArray(4) = VelXold-PartState(PartID,4)
    TransArray(5) = VelYold-PartState(PartID,5)
    TransArray(6) = VelZold-PartState(PartID,6)

    !---- Internal energy accommodation
    IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(outSpec(1))%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = 0
      ErotNew  = 0
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      IF(SpecDSMC(outSpec(1))%PolyatomicMol) THEN
        EvibNew = 0.0
        iPolyatMole = SpecDSMC(outSpec(1))%SpecToPolyArray
        ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        END DO
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(outSpec(1))%MacroParticleFactor
          IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
          IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib) &
                    - DSMC%GammaQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(1))%MaxVibQuant)
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        END DO
        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
      END IF
      SDEALLOCATE(RanNumPoly)
      SDEALLOCATE(VibQuantWallPoly)
    END IF
    END IF
    !End internal energy accomodation
    !----  Sampling of energies
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,AdsorptionEnthalpie&
        ,locBCID)
    END IF
  END IF









!    ! perform perfect reflection
!    v_aux                   = -2.0*((LengthPartTrajectory-alpha)*DOT_PRODUCT(PartTrajectory(1:3),n_loc))*n_loc
!    LastPartPos(PartID,1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
!    PartState(PartID,1:3)   = PartState(PartID,1:3)+v_aux
!    ! new velocity vector (without chemistry)
!    v_2                  =(LengthPartTrajectory-alpha)*PartTrajectory(1:3)+v_aux
!    PartState(PartID,4:6)= SQRT(DOT_PRODUCT(PartState(PartID,4:6),PartState(PartID,4:6)))*&
!                           (1/(SQRT(DOT_PRODUCT(v_2,v_2))))*v_2 + WallVelo(1:3)
!    ! calculate new velocity vector (Extended Maxwellian Model)
!    VeloReal = SQRT(PartState(PartID,4) * PartState(PartID,4) + &
!                    PartState(PartID,5) * PartState(PartID,5) + &
!                    PartState(PartID,6) * PartState(PartID,6))
!    EtraOld = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
!    ! set internal energies and do accomodation
!    IF (CollisMode.GT.1) THEN
!    IF (SpecDSMC(SpecID)%InterID.EQ.2) THEN ! recombination particle is molecule
!      !---- Rotational energy accommodation
!      CALL RANDOM_NUMBER(RanNum)
!      ErotWall = 0 !- BoltzmannConst * WallTemp * LOG(RanNum)
!      ErotNew  = 0 !PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
!      IntArray(1) = PartStateIntEn(PartID,2)
!      IntArray(2) = ErotWall
!      IntArray(3) = ErotNew
!      PartStateIntEn(PartID,2) = ErotNew
!      !---- Vibrational energy accommodation
!      ! first do reactant
!      IF(SpecDSMC(SpecID)%PolyatomicMol) THEN
!        iPolyatMole = SpecDSMC(SpecID)%SpecToPolyArray
!        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
!                        * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(SpecID)%MacroParticleFactor
!        END DO
!      ELSE
!        VibQuant = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(SpecID)%CharaTVib) &
!                    - DSMC%GammaQuant)
!        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(SpecID)%CharaTVib
!      END IF
!      ! then do first product
!      EvibNew = 0.0
!      IF (SpecDSMC(outSpec(1))%InterID.EQ.2) THEN
!        IF(SpecDSMC(outSpec(1))%PolyatomicMol) THEN
!          iPolyatMole = SpecDSMC(outSpec(1))%SpecToPolyArray
!!             ALLOCATE(VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
!          ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))!, &
!!                     VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!!                     VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
!          CALL RANDOM_NUMBER(RanNumPoly)
!          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
!          DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
!            CALL RANDOM_NUMBER(RanNumPoly)
!            VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
!          END DO
!!             VibQuantNewRPoly(:) = VibQuantWallPoly(:)
!!             VibQuantNewPoly = INT(VibQuantNewRPoly)
!          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!!               CALL RANDOM_NUMBER(RanNum)
!!               IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!!                 EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!!                           * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!!                 VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!!               ELSE
!!                 EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!!                           * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!!                 VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!!               END IF
!            IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
!                        * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
!            IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
!                        * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
!          END DO
!        ELSE
!          CALL RANDOM_NUMBER(RanNum)
!          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
!          DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(1))%MaxVibQuant)
!            CALL RANDOM_NUMBER(RanNum)
!            VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
!          END DO
!!             VibQuantNewR = VibQuantWall
!!             VibQuantNew = INT(VibQuantNewR)
!!             CALL RANDOM_NUMBER(RanNum)
!!             IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!!               EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!!             ELSE
!!               EvibNew = EvibNew + (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!!             END IF
!          IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
!          IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
!        END IF
!      END IF
!    ELSE ! recombination particle is atomic
!
!    END IF
!    END IF
!    ! set new species and velocity (after chemistry)
!    PartSpecies(PartID) = outSpec(2)
!    VeloRealNew = SQRT(2 * EtraOld / Species(PartSpecies(PartID))%MassIC)
!    PartState(PartID,4:6)= VeloRealNew/VeloReal * PartState(PartID,4:6)
!    ! correction for particle position
!    PartTrajectory       = VeloRealNew/VeloReal * (PartState(PartID,1:3) - LastPartPos(PartID,1:3))
!    PartState(PartID,1:3)= PartTrajectory + LastPartPos(PartID,1:3)
!    ! determine length of part trajectory
!    lengthPartTrajectory = SQRT(PartTrajectory(1)*PartTrajectory(1) &
!                           +PartTrajectory(2)*PartTrajectory(2) &
!                           +PartTrajectory(3)*PartTrajectory(3) )
!    ! normalized trajectory
!    PartTrajectory       = PartTrajectory/lengthPartTrajectory
!    lengthPartTrajectory = lengthPartTrajectory

!    !----  Sampling of energies
!    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
!      CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,AdsorptionEnthalpie&
!        ,locBCID)
!    END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT ! diffuse reflection
!-----------------------------------------------------------------------------------------------------------------------------------
  adsindex = 0
END SELECT

END SUBROUTINE ReactiveSurfaceTreatment

SUBROUTINE ParticleCondensationCase(PartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,condensindex,BCSideID,TriNum)
!===================================================================================================================================
! Routine for Selection of Liquid interaction (Condensation or Reflection)
!===================================================================================================================================
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_DSMC_Analyze           ,ONLY: CalcWallSample
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: PartState,Species,PartSpecies
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: Liquid
USE MOD_DSMC_Vars              ,ONLY: CollisMode, PolyatomMolDSMC
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, SpecDSMC, DSMC, VibQuantsPar
USE MOD_Particle_Boundary_Vars ,ONLY: SurfMesh, dXiEQ_SurfSample, Partbound
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, time
USE MOD_Particle_Surfaces_vars ,ONLY: SideNormVec,SideType
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)            :: condensindex
REAL,INTENT(INOUT)               :: PartTrajectory(1:3), alpha
REAL,INTENT(IN)                  :: xi, eta
INTEGER,INTENT(IN)               :: PartID, GlobSideID
LOGICAL,INTENT(IN)               :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL      :: BCSideID
INTEGER,INTENT(IN),OPTIONAL      :: TriNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: n_loc(1:3), tang1(1:3),tang2(1:3)
REAL                             :: Condensation_prob
INTEGER                          :: adsorption_case
INTEGER                          :: SurfSideID, SpecID
INTEGER                          :: outSpec(2)
! variables for Energy sampling
REAL                             :: TransArray(1:6),IntArray(1:6), EvaporationEnthalpie
REAL                             :: VelXold, VelYold, VelZold
INTEGER                          :: locBCID, VibQuant, VibQuantWall
!   INTEGER                          :: VibQuantNew
!   REAL                             :: VibQuantNewR
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
!   REAL                             :: TransACC, VibACC, RotACC
REAL                             :: ErotNew, ErotWall, EVibNew
! Polyatomic Molecules
REAL, ALLOCATABLE                :: RanNumPoly(:)
INTEGER                          :: iPolyatMole, iDOF
INTEGER, ALLOCATABLE             :: VibQuantWallPoly(:)
!   REAL, ALLOCATABLE                :: VibQuantNewRPoly(:)
!   INTEGER, ALLOCATABLE             :: VibQuantNewPoly(:), VibQuantTemp(:)
!===================================================================================================================================

! additional states
locBCID=PartBound%MapToPartBC(BC(GlobSideID))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)

TransArray(:) = 0.0
IntArray(:) = 0.0

SurfSideID = SurfMesh%SideIDToSurfID(GlobSideID)
SpecID = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)
! Update wallcollision counter
Liquid%Info(SpecID)%WallCollCount = Liquid%Info(SpecID)%WallCollCount + 1
#endif
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

IF(PRESENT(BCSideID))THEN
  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
  END SELECT
ELSE
  IF (TriaTracking) THEN
    CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=GlobSideID)
  ELSE
    SELECT CASE(SideType(GlobSideID))
    CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
      n_loc=SideNormVec(1:3,GlobSideID)
    CASE(BILINEAR)
      CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,GlobSideID)
    CASE(CURVED)
      CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,GlobSideID)
    END SELECT
  END IF
END IF

EvaporationEnthalpie = 0. ! negative at evaporation and positive at condensation

Condensation_prob = Liquid%ProbCondens(p,q,SurfSideID,SpecID)
CALL RANDOM_NUMBER(RanNum)
IF ( (Condensation_prob.GE.RanNum) ) THEN
  outSpec(1) = SpecID
  adsorption_case = 1
END IF

SELECT CASE(adsorption_case)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(-1) ! perfect elastic scattering
!-----------------------------------------------------------------------------------------------------------------------------------
  condensindex = -1
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) ! molecular condensation
!-----------------------------------------------------------------------------------------------------------------------------------
  Liquid%SumCondensPart(p,q,SurfSideID,outSpec(1)) = Liquid%SumCondensPart(p,q,SurfSideID,outSpec(1)) + 1
  condensindex = 1
#if (PP_TimeDiscMethod==42)
  Liquid%Info(SpecID)%NumOfAds = Liquid%Info(SpecID)%NumOfAds + 1
#endif
  IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    VelXold = PartState(PartID,4)
    VelYold = PartState(PartID,5)
    VelZold = PartState(PartID,6)
    PartState(PartID,4)  = WallVelo(1)
    PartState(PartID,5)  = WallVelo(2)
    PartState(PartID,6)  = WallVelo(3)

    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(outSpec(1))%MassIC * VeloReal**2
    EtraWall = 0.0
    EtraNew = EtraWall

    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    ! must be old_velocity-new_velocity
    TransArray(4) = VelXold-PartState(PartID,4)
    TransArray(5) = VelYold-PartState(PartID,5)
    TransArray(6) = VelZold-PartState(PartID,6)

    !---- Internal energy accommodation
    IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(outSpec(1))%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = 0
      ErotNew  = 0
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      IF(SpecDSMC(outSpec(1))%PolyatomicMol) THEN
        EvibNew = 0.0
        iPolyatMole = SpecDSMC(outSpec(1))%SpecToPolyArray
        ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF))!, &
!                   VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                   VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        END DO
!           VibQuantNewRPoly(:) = VibQuantWallPoly(:)
!           VibQuantNewPoly = INT(VibQuantNewRPoly)
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!             CALL RANDOM_NUMBER(RanNum)
!             IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!             ELSE
!               EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!                         * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!               VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!             END IF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(outSpec(1))%MacroParticleFactor
          IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
          IntArray(6) = IntArray(6) + (VibQuantWallPoly(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(outSpec(1))%CharaTVib * Species(outSpec(1))%MacroParticleFactor
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib) &
                    - DSMC%GammaQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        DO WHILE (VibQuantWall.GE.SpecDSMC(outSpec(1))%MaxVibQuant)
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(outSpec(1))%CharaTVib)
        END DO
!           VibQuantNewR = VibQuantWall
!           VibQuantNew = INT(VibQuantNewR)
!           CALL RANDOM_NUMBER(RanNum)
!           IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!             EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!           ELSE
!             EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(outSpec(1))%CharaTVib
!           END IF
        IntArray(4) = (VibQuant + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
        IntArray(6) = (VibQuantWall + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(outSpec(1))%CharaTVib
      END IF
      SDEALLOCATE(RanNumPoly)
      SDEALLOCATE(VibQuantWallPoly)
!         SDEALLOCATE(VibQuantNewRPoly)
!         SDEALLOCATE(VibQuantNewPoly)
!         SDEALLOCATE(VibQuantTemp)
!         IntArray(6) = EvibNew
    END IF
    END IF
    !End internal energy accomodation

    !----  Sampling of energies
    IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
      CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,EvaporationEnthalpie&
        ,locBCID)
    END IF
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT ! diffuse reflection
!-----------------------------------------------------------------------------------------------------------------------------------
  condensindex = 0
END SELECT

END SUBROUTINE ParticleCondensationCase


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
