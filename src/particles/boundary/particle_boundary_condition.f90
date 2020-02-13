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
! Determines how particles interact with a given boundary condition. This routine is used by MOD_Part_Tools, hence, it cannot be 
! used here due to circular definitions! 
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

INTERFACE GetBoundaryInteractionAuxBC
  MODULE PROCEDURE GetBoundaryInteractionAuxBC
END INTERFACE

INTERFACE PartSwitchElement
  MODULE PROCEDURE PartSwitchElement
END INTERFACE

PUBLIC::GetBoundaryInteraction,GetBoundaryInteractionAuxBC,PartSwitchElement
!===================================================================================================================================

CONTAINS

SUBROUTINE GetBoundaryInteraction(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,flip,ElemID,crossedBC&
                                  ,TriNum,locSideID)
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
USE MOD_Globals                  ,ONLY: abort
USE MOD_Particle_Surfaces        ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars            ,ONLY: PDM, UseCircularInflow
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackingMethod
USE MOD_Particle_Mesh_Vars       ,ONLY: PartBCSideList
USE MOD_Particle_Boundary_Vars   ,ONLY: PartBound,nPorousBC,DoBoundaryParticleOutput
USE MOD_Particle_Boundary_Porous ,ONLY: PorousBoundaryTreatment
USE MOD_Particle_Surfaces_vars   ,ONLY: SideNormVec,SideType
USE MOD_SurfaceModel             ,ONLY: ReactiveSurfaceTreatment
USE MOD_part_operations          ,ONLY: RemoveParticle
USE MOD_Mesh_Vars                ,ONLY: BC
#if defined(IMPA)
USE MOD_Particle_Vars            ,ONLY: PartIsImplicit
USE MOD_Particle_Vars            ,ONLY: DoPartInNewton
#endif /*IMPA*/
USE MOD_Dielectric_Vars          ,ONLY: DoDielectricSurfaceCharge
USE MOD_Particle_Vars            ,ONLY: LastPartPos
USE MOD_Particle_Boundary_Tools  ,ONLY: BoundaryParticleOutput,DielectricSurfaceCharge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: iPart,SideID,flip
REAL,INTENT(IN)                      :: xi,eta
INTEGER,INTENT(IN),OPTIONAL          :: TriNum
INTEGER,INTENT(IN),OPTIONAL          :: locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                :: ElemID
REAL,INTENT(INOUT)                   :: alpha,PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT)                  :: crossedBC
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: n_loc(1:3),RanNum
INTEGER                              :: ReflectionIndex, BCSideID
LOGICAL                              :: isSpeciesSwap, ElasticReflectionAtPorousBC
!===================================================================================================================================

IsSpeciesSwap=.FALSE.
crossedBC    =.FALSE.

! Calculate normal vector
BCSideID=SideID
SELECT CASE(TrackingMethod)
CASE(REFMAPPING,TRACING)
  ! set BCSideID for normal vector calculation call with (curvi-)linear side description
  IF (TrackingMethod.EQ.REFMAPPING) BCSideID=PartBCSideList(SideID)

  SELECT CASE(SideType(BCSideID))
  CASE(PLANAR_RECT,PLANAR_NONRECT,PLANAR_CURVED)
    n_loc=SideNormVec(1:3,BCSideID)
  CASE(BILINEAR)
    CALL CalcNormAndTangBilinear(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  CASE(CURVED)
    CALL CalcNormAndTangBezier(nVec=n_loc,xi=xi,eta=eta,SideID=BCSideID)
  END SELECT

  IF(flip.NE.0) n_loc=-n_loc

  ! Inserted particles are "pushed" inside the domain and registered as passing through the BC side. If they are very close to the
  ! boundary (first if) than the normal vector is compared with the trajectory. If the particle is entering the domain from outside
  ! it was inserted during surface flux and this routine shall not performed.
  ! Comparing the normal vector with the particle trajectory, if the particle trajectory is pointing inside the domain
  IF(DOT_PRODUCT(n_loc,PartTrajectory).LE.0.) RETURN
CASE(TRIATRACKING)
  CALL CalcNormAndTangTriangle(nVec=n_loc,TriNum=TriNum,SideID=BCSideID)
END SELECT
! required for refmapping and tracing, optional for triatracking
crossedBC=.TRUE.

IF (.NOT. ALLOCATED(PartBound%MapToPartBC)) THEN
  CALL abort(&
  __STAMP__&
  ,' ERROR: PartBound not allocated!.')
END IF

ASSOCIATE( iBC => PartBound%MapToPartBC(BC(SideID)) )
  ! Surface particle output to .h5
  IF(DoBoundaryParticleOutput.AND.PartBound%BoundaryParticleOutput(iBC))THEN
    CALL BoundaryParticleOutput(iPart,LastPartPos(1:3,iPart)+PartTrajectory(1:3)*alpha,PartTrajectory(1:3),n_loc)
  END IF

  ! Select the corresponding boundary condition and calculate particle treatment
  SELECT CASE(PartBound%TargetBoundCond(iBC))
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(1) !PartBound%OpenBC)
  !----------------------------------------------------------------------------------------------------------------------------------
  CALL RemoveParticle(iPart,BCID=iBC,alpha=alpha)
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(2) !PartBound%ReflectiveBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
  !---- Treatment of adaptive and porous boundary conditions (deletion of particles in case of circular inflow or porous BC)
  ElasticReflectionAtPorousBC = .FALSE.
  IF(UseCircularInflow) CALL SurfaceFluxBasedBoundaryTreatment(iPart,SideID,alpha,PartTrajectory)
  IF(nPorousBC.GT.0) CALL PorousBoundaryTreatment(iPart,SideID,alpha,PartTrajectory,ElasticReflectionAtPorousBC)

  !---- Dielectric particle-surface interaction
  IF(DoDielectricSurfaceCharge.AND.PartBound%Dielectric(iBC)) CALL DielectricSurfaceCharge(iPart,ElemID,PartTrajectory,alpha)

  !---- swap species?
  IF (PartBound%NbrOfSpeciesSwaps(iBC).gt.0) THEN
    CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,n_loc,iPart,SideID,IsSpeciesSwap)
  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! Particle did not Swap to species 0 (deleted particle -> particle is swapped to species 0)
    IF (PartBound%Reactive(iBC)) THEN
      ! Decide which interaction (reflection, reaction, adsorption)
      CALL ReactiveSurfaceTreatment(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,IsSpeciesSwap &
          ,ReflectionIndex)
    ELSE
      ! simple reflection (Maxwellian scattering)
      ReflectionIndex=2 ! diffuse reflection
      CALL RANDOM_NUMBER(RanNum)
      IF(RanNum.GE.PartBound%MomentumACC(iBC).OR.ElasticReflectionAtPorousBC) ReflectionIndex=1 ! perfect reflection
    END IF
    ! assign right treatment
    SELECT CASE (ReflectionIndex)
    CASE(-2) ! special case for double check that needs to be performed because particle moves away from surface
      ! can happen if particle was reflected or inserted on the surface and consequently alpha is almost 0
    CASE(1) !elastic reflection
      CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,IsSpeciesSwap)
    CASE(2) ! inelastic reflection
      CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,IsSpeciesSwap)
    CASE(3) ! reflection performed in reactive treatment routine
    CASE DEFAULT
      CALL abort(&
          __STAMP__&
          ,'ERROR: wrong interaction case in Boundary Condition! ReflectionIndex=',IntInfoOpt=ReflectionIndex)
    END SELECT
  END IF
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(3) !PartBound%PeriodicBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,ElemID) ! opt_reflected is peri-moved
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(4) !PartBound%SimpleAnodeBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (PartBound%SimpleAnodeBC)')
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(5) !PartBound%SimpleCathodeBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (PartBound%SimpleCathodeBC)')
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(6) !PartBound%MPINeighborhoodBC)
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (PartBound%MPINeighborhoodBC)')
  !-----------------------------------------------------------------------------------------------------------------------------------
  CASE(10,11) !PartBound%SymmetryBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL  PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,iPart,SideID,n_loc,IsSpeciesSwap,opt_Symmetry=.TRUE.)
  CASE(100) !PartBound%AnalyzeBC
  !-----------------------------------------------------------------------------------------------------------------------------------
    CALL  SideAnalysis(PartTrajectory,alpha,xi,eta,iPart,SideID,locSideID,ElemID,IsSpeciesSwap)
  CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,' ERROR: PartBound not associated!. (unknown case)')
END SELECT !PartBound%MapToPartBC(BC(SideID)
END ASSOCIATE

END SUBROUTINE GetBoundaryInteraction


SUBROUTINE GetBoundaryInteractionAuxBC(PartTrajectory,lengthPartTrajectory,alpha,iPart,AuxBCIdx,crossedBC)
!===================================================================================================================================
! Computes the post boundary state of a particle that interacts with an auxBC
!  OpenBC                  = 1
!  ReflectiveBC            = 2
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals                ,ONLY: abort,UNITVECTOR
USE MOD_Particle_Vars          ,ONLY: PDM
USE MOD_Particle_Boundary_Vars ,ONLY: PartAuxBC
USE MOD_Particle_Boundary_Vars ,ONLY: AuxBCType,AuxBCMap,AuxBC_plane,AuxBC_cylinder,AuxBC_cone,AuxBC_parabol
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_part_operations        ,ONLY: RemoveParticle
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
REAL                                 :: n_loc(1:3)
REAL                                 :: intersec(3),r_vec(3),axis(3),cos2inv!,v_2(1:3),v_aux(1:3)
!===================================================================================================================================

IsSpeciesSwap=.FALSE.
crossedBC    =.FALSE.
SELECT CASE (TRIM(AuxBCType(AuxBCIdx)))
CASE ('plane')
  n_loc = AuxBC_plane(AuxBCMap(AuxBCIdx))%n_vec
CASE ('cylinder')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
  r_vec = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cylinder(AuxBCMap(AuxBCIdx))%axis
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis) ) )
  IF (.NOT.AuxBC_cylinder(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('cone')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
  r_vec = AuxBC_cone(AuxBCMap(AuxBCIdx))%r_vec
  axis  = AuxBC_cone(AuxBCMap(AuxBCIdx))%axis
  cos2inv = 1./COS(AuxBC_cone(AuxBCMap(AuxBCIdx))%halfangle)**2
  n_loc = UNITVECTOR( intersec - ( r_vec + axis*DOT_PRODUCT(intersec-r_vec,axis)*cos2inv ) )
  IF (.NOT.AuxBC_cone(AuxBCMap(AuxBCIdx))%inwards) n_loc=-n_loc
CASE ('parabol')
  intersec = LastPartPos(1:3,iPart) + alpha*PartTrajectory
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
  crossedBC=.FALSE.
  !RETURN
  CALL abort(&
    __STAMP__&
    ,'Error in GetBoundaryInteractionAuxBC: Particle coming from outside!')
ELSE IF(DOT_PRODUCT(n_loc,PartTrajectory).GT.0.)  THEN
  crossedBC=.TRUE.
ELSE
  CALL abort(&
    __STAMP__&
    ,'Error in GetBoundaryInteractionAuxBC: n_vec is perpendicular to PartTrajectory for AuxBC',AuxBCIdx)
END IF
! Select the corresponding boundary condition and calculate particle treatment
SELECT CASE(PartAuxBC%TargetBoundCond(AuxBCIdx))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1) !PartAuxBC%OpenBC
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL RemoveParticle(iPart,alpha=alpha)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2) !PartAuxBC%ReflectiveBC)
!-----------------------------------------------------------------------------------------------------------------------------------
  !---- swap species?
!print*,'*********************'
!print*,AuxBCIdx
!print*,iPart,alpha,PartState(4:6,iPart)
!print*,iPart,alpha,LastPartPos(1:3,iPart),PartState(4:6,iPart)
  IF (PartAuxBC%NbrOfSpeciesSwaps(AuxBCIdx).gt.0) THEN
! CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1, &
!      IsSpeciesSwap=IsSpeciesSwap,flip=-1,AuxBCIdx=AuxBCIdx)
#ifndef IMPA
    CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,n_loc=n_loc,PartID=iPart,SideID=-1,IsSpeciesSwap=IsSpeciesSwap, &
      AuxBCIdx=AuxBCIdx)
#else
    CALL SpeciesSwap(PartTrajectory,alpha,xi=-1.,eta=-1.,n_loc=n_loc,PartID=iPart,SideID=-1,IsSpeciesSwap=IsSpeciesSwap)
#endif /*NOT IMPA*/
  END IF
  IF (PDM%ParticleInside(iPart)) THEN ! particle did not Swap to species 0 !deleted particle -> particle swaped to species 0
    ! simple reflection (previously used wall interaction model, maxwellian scattering)
      CALL RANDOM_NUMBER(RanNum)
      IF(RanNum.GE.PartAuxBC%MomentumACC(AuxBCIdx)) THEN
        ! perfectly reflecting, specular re-emission
        CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,n_loc=n_loc, &
          IsSpeciesSwap=IsSpeciesSwap,AuxBCIdx=AuxBCIdx)
      ELSE
        CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi=-1.,eta=-1.,PartID=iPart,SideID=-1,n_loc=n_loc, &
          IsSpeciesSwap=IsSpeciesSwap,AuxBCIdx=AuxBCIdx)
      END IF
  END IF
!print*,iPart,alpha,LastPartPos(1:3,iPart),PartState(1:3,iPart)
!print*,iPart,alpha,PartState(4:6,iPart)
!print*,'*********************'
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
CALL abort(&
__STAMP__&
,' ERROR: AuxBC bound not associated!. (unknown case)',999,999.)
END SELECT

! compiler warnings
RETURN
WRITE(*,*) 'AuxBCIdx', AuxBCIdx

END SUBROUTINE GetBoundaryInteractionAuxBC


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_Loc,IsSpeciesSwap, &
                             opt_Symmetry,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,SurfMesh,SampWall,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,nSpecies,PartSpecies,Species,WriteMacroSurfaceValues,PartLorentzType
USE MOD_Particle_Vars           ,ONLY: VarTimeStep
USE MOD_Mesh_Vars               ,ONLY: BC
USE MOD_DSMC_Vars               ,ONLY: DSMC,RadialWeighting,PartStateIntEn
USE MOD_Particle_Vars           ,ONLY: WriteMacroSurfaceValues, usevMPF
USE MOD_TImeDisc_Vars           ,ONLY: tend,time
USE MOD_Equation_Vars           ,ONLY: c2_inv
#if defined(LSERK)
USE MOD_Particle_Vars           ,ONLY: Pt_temp,PDM
#elif (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars           ,ONLY: PEM
#endif
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact
USE MOD_part_tools              ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID !,ElemID
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_old(1:3),WallVelo(3)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
!#if defined(LSERK)
!REAL                                 :: absPt_temp
!#endif
!REAL,PARAMETER                       :: oneMinus=0.99999999
!REAL                                 :: oneMinus!=0.99999999
REAL                                 :: LorentzFac, LorentzFacInv
REAL                                 :: Xitild,EtaTild
INTEGER                              :: p,q, SurfSideID, locBCID
LOGICAL                              :: Symmetry, IsAuxBC
REAL                                 :: MacroParticleFactor
LOGICAL                              :: DoSample
REAL                                 :: EtraOld,EvibOld,ErotOld
!===================================================================================================================================

IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF
IF (IsAuxBC) THEN
  WallVelo=PartAuxBC%WallVelo(1:3,AuxBCIdx)
ELSE
  WallVelo=PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
  locBCID=PartBound%MapToPartBC(BC(SideID))

  IF(PRESENT(opt_Symmetry)) THEN
    Symmetry = opt_Symmetry
  ELSE
    Symmetry = .FALSE.
  END IF
END IF !IsAuxBC

IF(SUM(ABS(WallVelo)).GT.0.)THEN
  SELECT CASE(PartLorentzType)
  CASE(3)
    v_old = PartState(4:6,PartID)
    PartState(4:6,PartID) = PartState(4:6,PartID) &
                          - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
    ! sanity check of new particle velocity
    LorentzFac=1.0-DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    IF(LorentzFac.LT.0.) CALL Abort(&
__STAMP__&
,'Particle exceeds speed of light! PartID ',PartID)
  CASE(5)
    ! map relativistic momentum to velocity
    LorentzFacInv         = 1.0+DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    LorentzFacInv         = 1.0/SQRT(LorentzFacInv)
    PartState(4:6,PartID) = LorentzFacInv*PartState(4:6,PartID)
    v_old                 = PartState(4:6,PartID)
    ! update velocity
    PartState(4:6,PartID) = PartState(4:6,PartID) &
                          - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
    ! map back from velocity to relativistic momentum
    LorentzFac=1.0-DOT_PRODUCT(PartState(4:6,PartID),PartState(4:6,PartID))*c2_inv
    IF(LorentzFac.LT.0.)THEN
CALL Abort(&
__STAMP__&
,'Particle exceeds speed of light! PartID ',PartID)
    END IF
    LorentzFac=1.0/SQRT(LorentzFac)
    PartState(4:6,PartID) = LorentzFac*PartState(4:6,PartID)
  CASE DEFAULT
    v_old = PartState(4:6,PartID)
    PartState(4:6,PartID) = PartState(4:6,PartID) &
                          - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc + WallVelo
  END SELECT
ELSE
  v_old = PartState(4:6,PartID)
  PartState(4:6,PartID) = PartState(4:6,PartID) &
                        - 2.*DOT_PRODUCT(PartState(4:6,PartID),n_loc)*n_loc
END IF

DoSample = (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)

IF (.NOT.IsAuxBC) THEN
  ! Wall sampling Macrovalues
  IF(.NOT.Symmetry) THEN !surface mesh is not built for the symmetry BC!?!
    IF (DoSample) THEN ! DoSample
      SurfSideID=SurfMesh%SideIDToSurfID(SideID)
      ! compute p and q
      ! correction of xi and eta, can only be applied if xi & eta are not used later!
      IF (TrackingMethod.EQ.TRIATRACKING) THEN
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
      IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
        MacroParticleFactor = GetParticleWeight(PartID)
      ELSE
        MacroParticleFactor = GetParticleWeight(PartID)*Species(PartSpecies(PartID))%MacroParticleFactor
      END IF
      !----  Sampling Forces at walls
      SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
          * (v_old(1:3) - PartState(4:6,PartID)) * MacroParticleFactor
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
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(1:3,PartID) + alpha * PartTrajectory(1:3)
          !-- caution: for consistency with diffuse refl. v_old is used!
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = v_old(1:3)
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(1:3,PartID)
          AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
          AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
        END IF
      END IF

      ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
      IF(CalcSurfaceImpact) THEN
        EtraOld = 0.5*Species(PartSpecies(PartID))%MassIC*VECNORM(v_old)**2
        IF(ALLOCATED(PartStateIntEn))THEN
          EvibOld=PartStateIntEn(1,PartID)
          ErotOld=PartStateIntEn(2,PartID)
        ELSE
          EvibOld=0.
          ErotOld=0.
        END IF ! ALLOCATED(PartStateIntEn)
        CALL CountSurfaceImpact(SurfSideID,PartSpecies(PartID),MacroParticleFactor,EtraOld,EvibOld,ErotOld,PartTrajectory,n_loc,p,q)
      END IF ! CalcSurfaceImpact
    END IF ! DoSample
  END IF ! .NOT.Symmetry
END IF !.NOT.IsAuxBC


! set particle position on face
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

PartTrajectory(1:3)=PartTrajectory(1:3)-2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc
PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)

! #if !defined(IMPA) &&  !defined(ROS)
! compute moved particle || rest of movement
PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                         +PartTrajectory(2)*PartTrajectory(2) &
                         +PartTrajectory(3)*PartTrajectory(3) )
PartTrajectory=PartTrajectory/lengthPartTrajectory
! #endif

#if defined(LSERK) || (PP_TimeDiscMethod==509)
!#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=506)
   ! correction for Runge-Kutta (correct position!!)
!---------- old ----------
!  absPt_temp=SQRT(Pt_temp(1,PartID)*Pt_temp(1,PartID)+Pt_temp(2,PartID)*Pt_temp(2,PartID)+Pt_temp(3,PartID)*Pt_temp(3,PartID))
!  ! scale PartTrajectory to new Pt_temp
!  Pt_temp(1:3,PartID)=absPt_temp*PartTrajectory(1:3)
!  ! deleate force history
!  Pt_temp(4:6,PartID)=0.
!  ! what happens with force term || acceleration?
!-------------------------
IF (.NOT.ALMOSTZERO(DOT_PRODUCT(WallVelo,WallVelo))) THEN
  PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#if defined(LSERK)
ELSE
  Pt_temp(1:3,PartID)=Pt_temp(1:3,PartID)-2.*DOT_PRODUCT(Pt_temp(1:3,PartID),n_loc)*n_loc
  IF (Symmetry) THEN !reflect also force history for symmetry
    Pt_temp(4:6,PartID)=Pt_temp(4:6,PartID)-2.*DOT_PRODUCT(Pt_temp(4:6,PartID),n_loc)*n_loc
  ELSE
    Pt_temp(4:6,PartID)=0. !produces best result compared to analytical solution in plate capacitor...
  END IF
#endif  /*LSERK*/
END IF
#endif  /*LSERK || (PP_TimeDiscMethod==509)*/

! rotation for IMEX and Rosenbrock Method (requires the rotation of the previous rk-stages... simplification of boundary condition)
! results in an order reduction
#ifdef IMPA
!IF(SUM(ABS(PEM%NormVec(1:3,PartID))).GT.0)THEN
!   IPWRITE(*,*) ' Caution: Field rotation for several reflection is not implemented!', iStage,PartIsImplicit(PartID), PartID
! END IF
PEM%NormVec(1:3,PartID)=n_loc
#endif /*IMPA*/
#ifdef ROS
! IF(SUM(ABS(PEM%NormVec(1:3,PartID))).GT.0)THEN
!   !IPWRITE(*,*) ' Caution: Field rotation for several reflection is not implemented!'
! END IF
PEM%NormVec(1:3,PartID)=n_loc
#endif /*ROS*/

END SUBROUTINE PerfectReflection


SUBROUTINE DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_loc,IsSpeciesSwap,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the diffuse reflection in 3D
! only implemented for RefMapping tracking
! PartBCs are reduced!
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                 ,ONLY: abort, OrthoNormVec,VECNORM
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod, TrackInfo
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,SurfMesh,SampWall,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact, GetWallTemperature
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,Species,PartSpecies,nSpecies,WriteMacroSurfaceValues,Symmetry
USE MOD_Particle_Vars           ,ONLY: VarTimeStep, usevMPF
#if defined(LSERK) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
USE MOD_Mesh_Vars               ,ONLY: BC
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC,CollisMode
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn,DSMC, useDSMC, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: PolyatomMolDSMC, VibQuantsPar
USE MOD_TimeDisc_Vars           ,ONLY: dt,tend,time,RKdtFrac
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO, PartSideToElem
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Vars                ,ONLY: BGKDoVibRelaxation
#elif (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Vars             ,ONLY: FPDoVibRelaxation
#endif
USE MOD_part_tools              ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: locBCID, vibQuant, vibQuantNew, VibQuantWall
REAL                                 :: VibQuantNewR                                                !
!REAL,PARAMETER                       :: oneMinus=0.99999999
REAL                                 :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
REAL                                 :: EtraWall, EtraNew
REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
REAL                                 :: tang1(1:3),tang2(1:3), NewVelo(3)
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
! Symmetry
REAL                                :: rotVelY, rotVelZ, rotPosY, MacroParticleFactor, adaptTimeStep
REAL                                :: VelX, VelY, VelZ,VecX, VecY, VecZ
REAL                                :: Vector1(1:3), Vector2(1:3)
REAL                                :: nx, ny, nz, nVal
INTEGER                             :: LocSideID, ElemID, iLoop
LOGICAL                             :: DoSample
REAL                                :: EvibOld,ErotOld
!===================================================================================================================================
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF
IF (IsAuxBC) THEN
  WallVelo   = PartAuxBC%WallVelo(1:3,AuxBCIdx)
  WallTemp   = PartAuxBC%WallTemp(AuxBCIdx)
  TransACC   = PartAuxBC%TransACC(AuxBCIdx)
  VibACC     = PartAuxBC%VibACC(AuxBCIdx)
  RotACC     = PartAuxBC%RotACC(AuxBCIdx)
ELSE
  ! additional states
  locBCID=PartBound%MapToPartBC(BC(SideID))
  ! get BC values
  WallVelo   = PartBound%WallVelo(1:3,locBCID)
  WallTemp   = GetWallTemperature(PartID,locBCID,PartTrajectory,alpha)
  TransACC   = PartBound%TransACC(locBCID)
  VibACC     = PartBound%VibACC(locBCID)
  RotACC     = PartBound%RotACC(locBCID)

END IF !IsAuxBC
CALL OrthoNormVec(n_loc,tang1,tang2)

IF(Symmetry%Axisymmetric) THEN
  ! Storing the old and the new particle position (which is outside the domain), at this point the position is only in the xy-plane
  VelX = PartState(1,PartID) - LastPartPos(1,PartID)
  VelY = PartState(2,PartID) - LastPartPos(2,PartID)
  VelZ = PartState(3,PartID) - LastPartPos(3,PartID)

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
VeloReal = VECNORM(PartState(4:6,PartID))

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

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  MacroParticleFactor = GetParticleWeight(PartID)
ELSE
  MacroParticleFactor = GetParticleWeight(PartID)*Species(PartSpecies(PartID))%MacroParticleFactor
END IF

DoSample=DSMC%CalcSurfaceVal.AND.((Time.GE.(1.-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroSurfaceValues)

IF (.NOT.IsAuxBC) THEN
  IF (DoSample) THEN
    !----  Sampling for energy (translation) accommodation at walls
    ! has to be corrected to new scheme
    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
    ! compute p and q
    ! correction of xi and eta, can only be applied if xi & eta are not used later!
    IF (TrackingMethod.EQ.TRIATRACKING) THEN
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

    ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
    IF(ALLOCATED(PartStateIntEn))THEN
      EvibOld=PartStateIntEn(1,PartID)
      ErotOld=PartStateIntEn(2,PartID)
    ELSE
      EvibOld=0.
      ErotOld=0.
    END IF ! ALLOCATED(PartStateIntEn)
    IF(CalcSurfaceImpact) CALL CountSurfaceImpact(SurfSideID,PartSpecies(PartID),MacroParticleFactor,EtraOld,EvibOld,ErotOld,&
                                                  PartTrajectory,n_loc,p,q)
  END IF
END IF !.NOT.IsAuxBC


!   Transformation local distribution -> global coordinates
! v = nv*u+t1*v+t2*f3

!NewVelo = VeloCx*tang1+CROSS(-n_loc,tang1)*VeloCy-VeloCz*n_loc
!---- Transformation local distribution -> global coordinates
IF(Symmetry%Axisymmetric) THEN
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
        ErotNew  = PartStateIntEn(2,PartID) + RotACC *(ErotWall - PartStateIntEn(2,PartID))

        IF (DoSample) THEN
          !----  Sampling for internal energy accommodation at walls
          SampWall(SurfSideID)%State(4,p,q)=SampWall(SurfSideID)%State(4,p,q)+PartStateIntEn(2,PartID) * MacroParticleFactor
          SampWall(SurfSideID)%State(5,p,q)=SampWall(SurfSideID)%State(5,p,q)+ErotWall * MacroParticleFactor
          SampWall(SurfSideID)%State(6,p,q)=SampWall(SurfSideID)%State(6,p,q)+ErotNew * MacroParticleFactor
        END IF

        PartStateIntEn(2,PartID) = ErotNew

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
            VibQuant     = NINT(PartStateIntEn(1,PartID)/(BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib) &
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

          IF (DoSample) THEN
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
          PartStateIntEn(1,PartID) = EvibNew
#if (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==300)
        END IF ! FPDoVibRelaxation || BGKDoVibRelaxation
#endif
      END IF
    END IF ! CollisMode > 1
  END IF ! useDSMC

  ! Sampling of the time step at the wall to get the correct time sample duration for the surface values
  IF (VarTimeStep%UseVariableTimeStep) THEN
    adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
    IF (DoSample) THEN
      SampWall(SurfSideID)%State(12+nSpecies+1,p,q) = SampWall(SurfSideID)%State(12+nSpecies+1,p,q) + adaptTimeStep
    END IF
  ELSE
    adaptTimeStep = 1.
  END IF ! VarTimeStep%UseVariableTimeStep

  ! intersection point with surface
  LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  !TildPos       =PartState(1:3,PartID)-dt*RKdtFrac*PartState(4:6,PartID)
  TildTrajectory=dt*RKdtFrac*PartState(4:6,PartID)*adaptTimeStep
  POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  !PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac * NewVelo(1:3)
  IF (IsAuxBC) THEN
    IF (PartAuxBC%Resample(AuxBCIdx)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  ELSE
    IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  END IF ! IsAuxBC
  PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3) * adaptTimeStep

  IF(Symmetry%Axisymmetric) THEN
    ! Symmetry considerations --------------------------------------------------------
    rotPosY = SQRT(PartState(2,PartID)**2 + (PartState(3,PartID))**2)
    ! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
    !           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
    ! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
    ! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
    rotVelY = (NewVelo(2)*(PartState(2,PartID))+NewVelo(3)*PartState(3,PartID))/rotPosY
    rotVelZ = (-NewVelo(2)*PartState(3,PartID)+NewVelo(3)*(PartState(2,PartID)))/rotPosY

    PartState(2,PartID) = rotPosY
    PartState(3,PartID) = 0.0
    NewVelo(2) = rotVelY
    NewVelo(3) = rotVelZ
  END IF ! Symmetry%Axisymmetric

  ! IF(Symmetry%Order.EQ.2) THEN
  !   ! z-Variable is set to zero (should be for the axisymmetric case anyway after rotation)
  !   lastPartPos(3,PartID) = 0.0
  !   PartState(3,PartID)   = 0.0
  ! END IF ! Symmetry%Order.EQ.2
  DO iLoop=Symmetry%Order+1,3
    LastPartPos(iLoop,PartID) = 0.0
    PartState(iLoop,PartID) = 0.0
  END DO

  IF (DoSample) THEN
    !----  Sampling force at walls
    SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) &
        + Species(PartSpecies(PartID))%MassIC * (PartState(4:6,PartID) - NewVelo(1:3)) * MacroParticleFactor
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
        END IF ! AnalyzeSurfCollis%Number(nSpecies+1) .GT. AnalyzeSurfCollis%maxPartNumber
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(1:3,PartID) + alpha * PartTrajectory(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = PartState(4:6,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(1:3,PartID)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
      END IF ! CalcSurfCollis%AnalyzeSurfCollis .AND. (ANY(AnalyzeSurfCollis%BCs.EQ.0) .OR. ANY(AnalyzeSurfCollis%BCs.EQ.locBCID))
    END IF ! .NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap
  END IF ! DoSample
END IF !.NOT.IsAuxBC

!----  saving new particle velocity
PartState(4:6,PartID)   = NewVelo(1:3) + WallVelo(1:3)

! recompute trajectory etc
IF(Symmetry%Axisymmetric) THEN
  PartTrajectory(1:2)=PartState(1:2,PartID) - LastPartPos(1:2,PartID)
  PartTrajectory(3) = 0.
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                          +PartTrajectory(2)*PartTrajectory(2))
ELSE
  PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                          +PartTrajectory(2)*PartTrajectory(2) &
                          +PartTrajectory(3)*PartTrajectory(3) )
END IF
IF(ABS(lengthPartTrajectory).GT.0.) PartTrajectory=PartTrajectory/lengthPartTrajectory

#if defined(LSERK) || (PP_TimeDiscMethod==509)
PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#endif

END SUBROUTINE DiffuseReflection


SUBROUTINE SpeciesSwap(PartTrajectory,alpha,xi,eta,n_Loc,PartID,SideID,IsSpeciesSwap,AuxBCIdx,targetSpecies_IN)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the Species Swap on ReflectiveBC
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                 ,ONLY: abort,VECNORM
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,SampWall,dXiEQ_SurfSample,SurfMesh,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,PartSpecies,usevMPF
USE MOD_Particle_Vars           ,ONLY: WriteMacroSurfaceValues,nSpecies,CollectCharges,nCollectChargesBCs,Species
USE MOD_Mesh_Vars               ,ONLY: BC
USE MOD_DSMC_Vars               ,ONLY: DSMC, RadialWeighting
USE MOD_TimeDisc_Vars           ,ONLY: TEnd,Time
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn
#if defined(IMPA)
USE MOD_Particle_Vars           ,ONLY: PartIsImplicit,DoPartInNewton
#endif /*IMPA*/
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_part_operations         ,ONLY: RemoveParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
INTEGER,INTENT(IN),OPTIONAL       :: targetSpecies_IN
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
LOGICAL                           :: IsAuxBC
REAL                              :: MacroParticleFactor
LOGICAL                           :: DoSample
REAL                              :: EtraOld,EvibOld,ErotOld
!===================================================================================================================================
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF

IF (IsAuxBC) THEN
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.LE.PartAuxBC%ProbOfSpeciesSwaps(AuxBCIdx)) THEN
    targetSpecies=-1 ! Dummy initialization value
    IF(PRESENT(targetSpecies_IN))THEN
      targetSpecies = targetSpecies_IN
    ELSE ! Normal swap routine
      DO iSwaps=1,PartAuxBC%NbrOfSpeciesSwaps(AuxBCIdx)
        IF (PartSpecies(PartID).eq.PartAuxBC%SpeciesSwaps(1,iSwaps,AuxBCIdx)) &
            targetSpecies = PartAuxBC%SpeciesSwaps(2,iSwaps,AuxBCIdx)
      END DO
    END IF ! PRESENT(targetSpecies_IN)
    !swap species
    IF (targetSpecies.ge.0) IsSpeciesSwap=.TRUE.
    IF (targetSpecies.eq.0) THEN !delete particle -> same as PartAuxBC%OpenBC
      CALL RemoveParticle(PartID,alpha=alpha)
    ELSEIF (targetSpecies.gt.0) THEN !swap species
      PartSpecies(PartID)=targetSpecies
    END IF
  END IF !RanNum.LE.PartAuxBC%ProbOfSpeciesSwaps
ELSE
  DoSample = (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)
  locBCID = PartBound%MapToPartBC(BC(SideID))
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.LE.PartBound%ProbOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID)))) THEN
    targetSpecies=-1 ! Dummy initialization value
    IF(PRESENT(targetSpecies_IN))THEN
      targetSpecies = targetSpecies_IN
    ELSE ! Normal swap routine
      DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(PartBound%MapToPartBC(BC(SideID)))
        IF (PartSpecies(PartID).eq.PartBound%SpeciesSwaps(1,iSwaps,PartBound%MapToPartBC(BC(SideID)))) &
            targetSpecies = PartBound%SpeciesSwaps(2,iSwaps,PartBound%MapToPartBC(BC(SideID)))
      END DO
    END IF ! PRESENT(targetSpecies_IN)
    !swap species
    IF (targetSpecies.ge.0) IsSpeciesSwap=.TRUE.
    IF ( (targetSpecies.eq.0) .OR. (.NOT.CalcSurfCollis%Only0Swaps) ) THEN
      IF (DoSample) THEN
        !---- Counter for swap species collisions
        SurfSideID=SurfMesh%SideIDToSurfID(SideID)
        ! compute p and q
        ! correction of xi and eta, can only be applied if xi & eta are not used later!
        IF (TrackingMethod.EQ.TRIATRACKING) THEN
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
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(1:3,PartID) + alpha * PartTrajectory(1:3)
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4:6) = PartState(4:6,PartID)
          AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7:9) = LastPartPos(1:3,PartID)
          AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1)) = PartSpecies(PartID)
          AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1)) = locBCID
        END IF
      END IF ! DoSample
    END IF
    IF (targetSpecies.eq.0) THEN !delete particle -> same as PartBound%OpenBC
      IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
        MacroParticleFactor = GetParticleWeight(PartID)
      ELSE
        MacroParticleFactor = GetParticleWeight(PartID)*Species(PartSpecies(PartID))%MacroParticleFactor
      END IF
      DO iCC=1,nCollectChargesBCs !-chargeCollect
        IF (CollectCharges(iCC)%BC .EQ. PartBound%MapToPartBC(BC(SideID))) THEN
          CollectCharges(iCC)%NumOfNewRealCharges = CollectCharges(iCC)%NumOfNewRealCharges &
              + Species(PartSpecies(PartID))%ChargeIC * MacroParticleFactor
          EXIT
        END IF
      END DO
      ! sample values of deleted species
      IF (DoSample) THEN
        SurfSideID=SurfMesh%SideIDToSurfID(SideID)
        IF (TrackingMethod.EQ.TRIATRACKING) THEN
          p=1 ; q=1
        ELSE
          Xitild =MIN(MAX(-1.,xi ),0.99)
          Etatild=MIN(MAX(-1.,eta),0.99)
          p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
          q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
        END IF
        !----  Sampling Forces at walls
        SampWall(SurfSideID)%State(10:12,p,q)= SampWall(SurfSideID)%State(10:12,p,q) + Species(PartSpecies(PartID))%MassIC &
            * PartState(4:6,PartID) * MacroParticleFactor
      END IF

      ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z), angle and number of impacts
      IF(CalcSurfaceImpact) THEN
        EtraOld = 0.5*Species(PartSpecies(PartID))%MassIC*VECNORM(PartState(4:6,PartID))**2
        IF(ALLOCATED(PartStateIntEn))THEN
          EvibOld=PartStateIntEn(1,PartID)
          ErotOld=PartStateIntEn(2,PartID)
        ELSE
          EvibOld=0.
          ErotOld=0.
        END IF ! ALLOCATED(PartStateIntEn)
        CALL CountSurfaceImpact(SurfSideID,PartSpecies(PartID),MacroParticleFactor,EtraOld,EvibOld,ErotOld,PartTrajectory,n_loc,p,q)
      END IF ! CalcSurfaceImpact
      CALL RemoveParticle(PartID,BCID=PartBound%MapToPartBC(BC(SideID)),alpha=alpha)
    ELSE IF (targetSpecies.gt.0) THEN !swap species
      DO iCC=1,nCollectChargesBCs !-chargeCollect
        IF (CollectCharges(iCC)%BC .EQ. PartBound%MapToPartBC(BC(SideID))) THEN
          CollectCharges(iCC)%NumOfNewRealCharges = CollectCharges(iCC)%NumOfNewRealCharges &
              + (Species(PartSpecies(PartID))%ChargeIC-Species(targetSpecies)%ChargeIC) * MacroParticleFactor
          EXIT
        END IF
      END DO
      PartSpecies(PartID)=targetSpecies
    END IF ! targetSpecies.eq.0
  END IF ! RanNum.LE.PartBound%ProbOfSpeciesSwaps
END IF ! IsAuxBC

END SUBROUTINE SpeciesSwap


SUBROUTINE PeriodicBC(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,ElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO,SidePeriodicType
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos,PEM
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
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
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), lengthPartTrajectory, alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID!,ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: PVID,moved(2),locSideID
!===================================================================================================================================

PVID = SidePeriodicType(SideID)

#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPos:      ',LastPartPos(1:3,PartID)
  END IF
END IF
#endif /*CODE_ANALYZE*/

! set last particle position on face
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha
! perform the periodic movement
LastPartPos(1:3,PartID) = LastPartPos(1:3,PartID) + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
! update particle positon after periodic BC
!PartState(1:3,PartID)   = PartState(1:3,PartID) + SIGN(GEO%PeriodicVectors(1:3,ABS(PVID)),REAL(PVID))
PartState(1:3,PartID) = LastPartPos(1:3,PartID) + (lengthPartTrajectory-alpha)*PartTrajectory
lengthPartTrajectory  = lengthPartTrajectory - alpha


#ifdef CODE_ANALYZE
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    IPWRITE(UNIT_stdout,'(I0,A)') '     PeriodicBC: '
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' ParticlePosition-pp: ',PartState(1:3,PartID)
    IPWRITE(UNIT_stdout,'(I0,A,3(X,G0))') ' LastPartPo-pp:       ',LastPartPos(1:3,PartID)
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
IF (TrackingMethod.EQ.REFMAPPING) PEM%LastElement(PartID) = 0

IF(1.EQ.2)THEN
  alpha=0.2
END IF

END SUBROUTINE PeriodicBC


SUBROUTINE SideAnalysis(PartTrajectory,alpha,xi,eta,PartID,SideID,locSideID,ElemID,IsSpeciesSwap)
!----------------------------------------------------------------------------------------------------------------------------------!
! Analyze particle crossing (inner) side
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,CalcSurfCollis,AnalyzeSurfCollis
USE MOD_Particle_Vars          ,ONLY: PartState,LastPartPos,nSpecies,PartSpecies,WriteMacroSurfaceValues
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_TImeDisc_Vars          ,ONLY: tend,time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(INOUT)                :: PartTrajectory(1:3), alpha
REAL,INTENT(IN)                   :: xi, eta
INTEGER,INTENT(IN)                :: PartID, SideID,locSideID
LOGICAL,INTENT(IN)                :: IsSpeciesSwap
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT),OPTIONAL    :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: WallVelo(3)
INTEGER                              :: locBCID
INTEGER                              :: moved(2)
!===================================================================================================================================

WallVelo=PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(SideID)))
locBCID=PartBound%MapToPartBC(BC(SideID))

! Wall sampling Macrovalues
!IF(.NOT.Symmetry) THEN !surface mesh is not build for the symmetry BC!?!
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
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),1:3) = LastPartPos(1:3,PartID) + alpha * PartTrajectory(1:3)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),4)   = PartState(4,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),5)   = PartState(5,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),6)   = PartState(6,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),7)   = LastPartPos(1,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),8)   = LastPartPos(2,PartID)
        AnalyzeSurfCollis%Data(AnalyzeSurfCollis%Number(nSpecies+1),9)   = LastPartPos(3,PartID)
        AnalyzeSurfCollis%Spec(AnalyzeSurfCollis%Number(nSpecies+1))     = PartSpecies(PartID)
        AnalyzeSurfCollis%BCid(AnalyzeSurfCollis%Number(nSpecies+1))     = locBCID
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
USE MOD_Particle_Mesh_Vars ,ONLY: PartElemToElemAndSide
USE MOD_Mesh_Vars          ,ONLY: MortarType
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


SUBROUTINE SurfaceFluxBasedBoundaryTreatment(iPart,SideID,alpha,PartTrajectory)
!===================================================================================================================================
! Treatment of particles at the boundary if adaptive surface BCs or circular inflows based on the surface flux are present
! Circular Inflow: Particles are deleted if within (allows multiple surface flux inflows defined by circles on a single boundary)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species, LastPartPos, PartSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_part_operations        ,ONLY: RemoveParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: iPart, SideID
REAL,INTENT(IN)                     :: PartTrajectory(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                  :: alpha
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: point(1:2), intersectionPoint(1:3), radius
INTEGER                             :: iSpec, iSF
!===================================================================================================================================

iSpec = PartSpecies(iPart)
DO iSF=1,Species(iSpec)%nSurfacefluxBCs
  IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.PartBound%MapToPartBC(BC(SideID))) THEN
    IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      intersectionPoint(1:3) = LastPartPos(1:3,iPart)+ alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(2))-Species(iSpec)%Surfaceflux(iSF)%origin(1)
      point(2)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(3))-Species(iSpec)%Surfaceflux(iSF)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
        CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(BC(SideID)),alpha=alpha)
      END IF
    END IF
  END IF
END DO

END SUBROUTINE SurfaceFluxBasedBoundaryTreatment

END MODULE MOD_Particle_Boundary_Condition
