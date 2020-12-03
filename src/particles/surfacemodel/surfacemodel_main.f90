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

MODULE MOD_SurfaceModel
!===================================================================================================================================
!> Main Module for surface model
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
PUBLIC :: SurfaceTreatment
PUBLIC :: PerfectReflection, DiffuseReflection, SpeciesSwap
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceTreatment(PartTrajectory,LengthPartTrajectory,alpha,xi,eta,PartID,SideID,ElemID,n_Loc,IsSpeciesSwap)
!===================================================================================================================================
!> Routine for Selection of Surface interaction
!===================================================================================================================================
USE MOD_Globals                 ,ONLY: abort,UNITVECTOR,OrthoNormVec
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
USE MOD_Part_Tools              ,ONLY: VeloFromDistribution
USE MOD_part_operations         ,ONLY: CreateParticle, RemoveParticle
USE MOD_Particle_Vars           ,ONLY: PartState,Species,PartSpecies
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_Particle_Vars           ,ONLY: LastPartPos, PEM
USE MOD_Particle_Boundary_Tools ,ONLY: SurfaceToPartEnergyInternal, CalcWallSample, AnalyzeSurfaceCollisions
USE MOD_Particle_Boundary_Tools ,ONLY: AddPartInfoToSample,CalcRotWallVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample, Partbound, CalcSurfaceImpact, GlobalSide2SurfSide
USE MOD_TimeDisc_Vars           ,ONLY: dt, RKdtFrac
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars       ,ONLY: SurfModel
USE MOD_SEE                     ,ONLY: SecondaryElectronEmission
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared

USE MOD_Particle_Boundary_Porous,ONLY: PorousBoundaryTreatment
USE MOD_Particle_Boundary_Tools ,ONLY: DielectricSurfaceCharge

USE MOD_Particle_Vars           ,ONLY: PDM
USE MOD_Particle_Vars           ,ONLY: UseCircularInflow
USE MOD_Particle_Boundary_Vars  ,ONLY: nPorousBC
USE MOD_Dielectric_Vars         ,ONLY: DoDielectricSurfaceCharge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)          :: PartTrajectory(1:3), LengthPartTrajectory, alpha
REAL,INTENT(IN)             :: xi, eta
REAL,INTENT(IN)             :: n_loc(1:3)
INTEGER,INTENT(IN)          :: PartID, SideID, ElemID
LOGICAL,INTENT(INOUT)       :: IsSpeciesSwap
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: ProductSpec(2)   !< 1: product species of incident particle (also used for simple reflection)
                                                     !< 2: additional species added or removed from surface
                                                          !< If productSpec is negative, then the respective particles are adsorbed
                                                          !< If productSpec is positive the particle is reflected/emitted
                                                          !< with respective species
INTEGER                          :: ProductSpecNbr   !< number of emitted particles for ProductSpec(1)
CHARACTER(30)                    :: velocityDistribution(2)   !< specifying keyword for velocity distribution
REAL                             :: TempErgy(2)               !< temperature, energy or velocity used for VeloFromDistribution
REAL                             :: PartTrajectory2(1:3)
INTEGER                          :: NewPartID
REAL                             :: RanNum
REAL                             :: Xitild,EtaTild
INTEGER                          :: p,q
REAL                             :: tang1(1:3),tang2(1:3)
INTEGER                          :: SurfSideID, SpecID
! variables for Energy sampling
REAL                             :: TransArray(1:6),IntArray(1:6)
REAL                             :: oldVelo(1:3)
INTEGER                          :: locBCID
REAL                             :: VeloReal, EtraOld
REAL                             :: EtraWall, EtraNew
REAL                             :: WallVelo(1:3), WallTemp
REAL                             :: TransACC!, VibACC, RotACC
! Polyatomic Molecules
REAL                             :: VeloCrad, Fak_D, NewVelo(3)
REAL                             :: Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                             :: POI_fak, TildTrajectory(3),POI_vec(3)
INTEGER                          :: iNewPart ! particle counter for newly created particles

INTEGER                         :: iBC, ReflectionIndex
LOGICAL                         :: ElasticReflectionAtPorousBC
!===================================================================================================================================
iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
! =============================
! Workflow:
!
!  1.  Initial surface checks:  Check incident velocity vector and surface normal
!  2.  Select surface model:    Determine what happens at the surface
!  3.  (New) Particle handling: Perform reflection/removal of incident particle, create (multiple) possible new particles
!==============================

!---- Treatment of adaptive and porous boundary conditions (deletion of particles in case of circular inflow or porous BC)
ElasticReflectionAtPorousBC = .FALSE.
IF(UseCircularInflow) CALL SurfaceFluxBasedBoundaryTreatment(PartID,SideID,alpha,PartTrajectory)
IF(nPorousBC.GT.0) CALL PorousBoundaryTreatment(PartID,SideID,alpha,PartTrajectory,ElasticReflectionAtPorousBC)

!---- Dielectric particle-surface interaction
IF(DoDielectricSurfaceCharge.AND.PartBound%Dielectric(iBC)) CALL DielectricSurfaceCharge(PartID,ElemID,PartTrajectory,alpha)

!---- swap species?
IF (PartBound%NbrOfSpeciesSwaps(iBC).gt.0) THEN
  CALL SpeciesSwap(PartTrajectory,alpha,xi,eta,n_loc,PartID,SideID,IsSpeciesSwap)
  ! Particle was deleted during the species swap, leave the routine 
  IF(.NOT.PDM%ParticleInside(PartID)) RETURN
END IF

!===================================================================================================================================
! 1.) Initial surface checks
! find normal vector two perpendicular tangential vectors (normal_vector points outwards !!!)
!===================================================================================================================================
CALL OrthoNormVec(n_loc,tang1,tang2)

! additional states
locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
! get BC values
WallVelo     = PartBound%WallVelo(1:3,locBCID)
WallTemp     = PartBound%WallTemp(locBCID)

POI_vec(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

IF(PartBound%RotVelo(locBCID)) THEN
  CALL CalcRotWallVelo(locBCID,PartID,POI_vec,WallVelo)
END IF

! initialize sampling arrays
TransArray(:) = 0.0
IntArray(:) = 0.0

SpecID = PartSpecies(PartID)

ReflectionIndex = -1 ! has to be reset in SurfaceModel, otherwise abort() will be called
ProductSpec(1) = SpecID
ProductSpec(2) = 0
ProductSpecNbr = 0
velocityDistribution(1:2)=''
TempErgy(1:2)=WallTemp

!===================================================================================================================================
! 2.) Select surface model
! Here, the surfacemodel decides how the particle is treated on the surface
!===================================================================================================================================
SELECT CASE(PartBound%SurfaceModel(locBCID))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (0) ! Maxwellian scattering (diffuse/specular reflection)
!-----------------------------------------------------------------------------------------------------------------------------------
  ReflectionIndex=2 ! diffuse reflection
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.GE.PartBound%MomentumACC(iBC).OR.ElasticReflectionAtPorousBC) ReflectionIndex=1 ! perfect reflection
  ! assign right treatment
  SELECT CASE (ReflectionIndex)
  CASE(1) ! elastic reflection
    CALL PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_loc,IsSpeciesSwap)
  CASE(2) ! inelastic reflection
    CALL DiffuseReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_loc,IsSpeciesSwap)
  CASE DEFAULT
    CALL abort(&
        __STAMP__&
        ,'ERROR: wrong interaction case in surface treatment! ReflectionIndex=',IntInfoOpt=ReflectionIndex)
  END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (5,6,7) ! Copied from CASE(1) and adjusted for secondary e- emission (SEE)
             ! 5: SEE by Levko2015
             ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
             ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Update wallcollision counter (currently here as SurfModel is only allocated for SurfaceModel GT 0)
  SurfModel%Info(SpecID)%WallCollCount = SurfModel%Info(SpecID)%WallCollCount + 1

  ! Get electron emission probability
  CALL SecondaryElectronEmission(PartBound%SurfaceModel(locBCID),PartID,locBCID,ReflectionIndex,ProductSpec,&
  ProductSpecNbr,TempErgy(2),velocityDistribution)
CASE DEFAULT
    CALL abort(&
      __STAMP__&
      ,'Not implemented anymore!')
END SELECT

!===================================================================================================================================
! 3.) (New) Particle handling
! Here, the incident particle is reflected/adsorbed and an additional product is emitted/adsorbed
!===================================================================================================================================
SELECT CASE(ReflectionIndex)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1,2) ! (particle is treated in boundary condition)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! 1: Diffuse reflection
  ! 2: Perfect elastic scattering
  RETURN
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3) ! reactive interaction case
!-----------------------------------------------------------------------------------------------------------------------------------
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
  SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
  ! Old particle
  IF (ProductSpec(1).LT.0) THEN
    SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  END IF

  ! New particle
  IF (ProductSpec(2).LT.0) THEN
    SurfModel%Info(SpecID)%NumOfAds = SurfModel%Info(SpecID)%NumOfAds + 1
  END IF
  !-----------------------------------------------------------
  ! Treat incident particle
  CALL AddPartInfoToSample(PartID,TransArray,IntArray,'old')
  ! Sample momentum, heatflux and collision counter on surface
  CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,&
                      impact_opt=CalcSurfaceImpact,PartTrajectory_opt=PartTrajectory,SurfaceNormal_opt=n_loc)
  CALL AnalyzeSurfaceCollisions(PartID,PartTrajectory,alpha,IsSpeciesSwap,locBCID)

  IF (ProductSpec(1).LE.0) THEN
    CALL RemoveParticle(PartID,alpha=alpha)
  ELSE
    oldVelo(1:3) = PartState(4:6,PartID)
    IF(TRIM(velocityDistribution(1)).NE.'') THEN
      ! sample new velocity for reflected particle
      NewVelo(1:3) = VeloFromDistribution(velocityDistribution(1),ProductSpec(1),TempErgy(1))
      ! important: n_loc points outwards
      PartState(4:6,PartID) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

      ! intersection point with surface
      LastPartPos(1:3,PartID) = POI_vec(1:3)
      ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
      TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
      POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
      ! travel rest of particle vector
      IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution

      ! recompute trajectory etc
      PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * PartState(4:6,PartID)
    ELSE
      IF (PartBound%MomentumACC(locBCID).GT.0.0) THEN
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
        Cmr         = SQRT(2.0 * EtraNew / (Species(ProductSpec(1))%MassIC * Fak_D))
        CALL RANDOM_NUMBER(RanNum)
        Phi     = 2.0 * PI * RanNum
        VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
        VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
        VeloCz  = Cmr * VeloCz
        NewVelo = VeloCx*tang1-tang2*VeloCy-VeloCz*n_loc
      ELSE
        ! perfect velocity reflection
        NewVelo(1:3) = oldVelo(1:3) - 2.*DOT_PRODUCT(oldVelo(1:3),n_loc)*n_loc
        ! mass changes, therefore velocity is scaled because momentum remains the same
        NewVelo(1:3) = NewVelo(1:3) * (Species(ProductSpec(1))%MassIC/Species(SpecID)%MassIC)
      END IF
      ! intersection point with surface
      LastPartPos(1:3,PartID) = POI_vec(1:3)
      ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
      TildTrajectory=dt*RKdtFrac*oldVelo(1:3)
      POI_fak=1.- (lengthPartTrajectory-alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
      ! travel rest of particle vector
      IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
      PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3)
      !----  saving new particle velocity
      PartState(4:6,PartID)   = NewVelo(1:3) + WallVelo(1:3)
    END IF
    PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
    lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
    PartTrajectory=PartTrajectory/lengthPartTrajectory

    ! set new species of reflected particle
    PartSpecies(PartID) = ProductSpec(1)
    ! Adding the energy that is transferred from the surface onto the internal energies of the particle
    CALL SurfaceToPartEnergyInternal(PartID,WallTemp)
    CALL AddPartInfoToSample(PartID,TransArray,IntArray,'new')
    ! Sample momentum, heatflux and collision counter on surface
    CALL CalcWallSample(PartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,emission_opt=.TRUE.)
  END IF

  !-----------------------------------------------------------
  ! Create new particles
  IF (ProductSpec(2).GT.0) THEN
    DO iNewPart = 1, ProductSpecNbr
      SurfModel%Info(ProductSpec(2))%NumOfDes = SurfModel%Info(ProductSpec(2))%NumOfDes + 1
      ! create new particle and assign correct energies
      ! sample newly created velocity
      NewVelo(1:3) = VeloFromDistribution(velocityDistribution(2),ProductSpec(2),TempErgy(2))
      ! Rotate velocity vector from global coordinate system into the surface local coordinates (important: n_loc points outwards)
      NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_Loc(1:3)*NewVelo(3) + WallVelo(1:3)

      PartTrajectory2=UNITVECTOR(NewVelo(1:3))

      CALL CreateParticle(ProductSpec(2),LastPartPos(1:3,PartID),PEM%GlobalElemID(PartID),NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
      ! Adding the energy that is transferred from the surface onto the internal energies of the particle
      CALL SurfaceToPartEnergyInternal(NewPartID,WallTemp)

      CALL AddPartInfoToSample(NewPartID,TransArray,IntArray,'new')
      CALL CalcWallSample(NewPartID,SurfSideID,p,q,Transarray,IntArray,IsSpeciesSwap,emission_opt=.TRUE.)
    END DO ! iNewPart = 1, ProductSpecNbr
  END IF
END SELECT

END SUBROUTINE SurfaceTreatment


SUBROUTINE PerfectReflection(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,SideID,n_Loc,IsSpeciesSwap, &
                             opt_Symmetry,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState,GlobalSide2SurfSide
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,nSpecies,PartSpecies,Species,WriteMacroSurfaceValues,PartLorentzType
USE MOD_Particle_Vars           ,ONLY: VarTimeStep
USE MOD_DSMC_Vars               ,ONLY: DSMC,RadialWeighting,PartStateIntEn, AmbipolElecVelo
USE MOD_Particle_Vars           ,ONLY: WriteMacroSurfaceValues,usevMPF
USE MOD_TImeDisc_Vars           ,ONLY: tend,time
USE MOD_Globals_Vars            ,ONLY: c2_inv
#if defined(LSERK)
USE MOD_Particle_Vars           ,ONLY: Pt_temp,PDM
#elif (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars           ,ONLY: PEM
#endif
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact,CalcRotWallVelo
USE MOD_part_tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
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
REAL                                 :: v_old(1:3),WallVelo(3), v_old_Ambi(1:3)
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
REAL                                 :: POI_vec(1:3)
REAL                                 :: adaptTimeStep
!===================================================================================================================================
! Initialize
adaptTimeStep = 1.0
IsAuxBC=.FALSE.
Symmetry = .FALSE.

IF (PRESENT(AuxBCIdx)) IsAuxBC=.TRUE.

IF (IsAuxBC) THEN
  WallVelo = PartAuxBC%WallVelo(1:3,AuxBCIdx)
ELSE
  WallVelo = PartBound%WallVelo(1:3,PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
  locBCID  = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  IF(PRESENT(opt_Symmetry)) Symmetry = opt_Symmetry
END IF !IsAuxBC

! Get Point Of Intersection
POI_vec(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

!IF(PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
!  CALL CalcRotWallVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)),PartID,POI_vec,WallVelo)
!END IF

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
  IF (DSMC%DoAmbipolarDiff) THEN
    IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
      v_old_Ambi = AmbipolElecVelo(PartID)%ElecVelo(1:3)
      AmbipolElecVelo(PartID)%ElecVelo(1:3) = AmbipolElecVelo(PartID)%ElecVelo(1:3) &
                     - 2.*DOT_PRODUCT(AmbipolElecVelo(PartID)%ElecVelo(1:3),n_loc)*n_loc
    END IF
  END IF
END IF

DoSample = (DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)


IF (.NOT.IsAuxBC) THEN
  ! Wall sampling Macrovalues
  IF(.NOT.Symmetry) THEN !surface mesh is not built for the symmetry BC!?!
    IF (DoSample) THEN ! DoSample
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
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
        adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
        ! Sampling of the time step at the wall to get the correct time sample duration for the force per area calculation
        SampWallState(12+nSpecies+1,p,q,SurfSideID) = SampWallState(12+nSpecies+1,p,q,SurfSideID) &
            + adaptTimeStep
      END IF
      IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
        MacroParticleFactor = GetParticleWeight(PartID)
      ELSE
        MacroParticleFactor = GetParticleWeight(PartID)*Species(PartSpecies(PartID))%MacroParticleFactor
      END IF
      !----  Sampling Forces at walls
      SampWallState(10:12,p,q,SurfSideID) = SampWallState(10:12,p,q,SurfSideID) + Species(PartSpecies(PartID))%MassIC &
          * (v_old(1:3) - PartState(4:6,PartID)) * MacroParticleFactor
      !---- Counter for collisions (normal wall collisions - not to count if only Swaps to be counted, IsSpeciesSwap: already counted)
      !       IF (.NOT.CalcSurfCollis%OnlySwaps) THEN
      IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
        SampWallState(12+PartSpecies(PartID),p,q,SurfSideID) = SampWallState(12+PartSpecies(PartID),p,q,SurfSideID) + 1
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
LastPartPos(1:3,PartID) = POI_vec(1:3)

PartTrajectory(1:3)     = PartTrajectory(1:3)-2.*DOT_PRODUCT(PartTrajectory(1:3),n_loc)*n_loc

! Check if rotated velo is used, otherwise mirror the LastPartPos
!IF(PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))))THEN
!  PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac &
!                          * PartState(4:6,PartID) * adaptTimeStep
!ELSE
  PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*(lengthPartTrajectory - alpha)
!END IF ! PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))

! #if !defined(IMPA) &&  !defined(ROS)
! compute moved particle || rest of movement
PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)

IF(ALMOSTZERO(VECNORM(PartTrajectory)))THEN
  lengthPartTrajectory= 0.0
ELSE
  lengthPartTrajectory=SQRT(PartTrajectory(1)*PartTrajectory(1) &
                           +PartTrajectory(2)*PartTrajectory(2) &
                           +PartTrajectory(3)*PartTrajectory(3) )
  PartTrajectory=PartTrajectory/lengthPartTrajectory
END IF
! #endif

#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
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
#endif  /*LSERK || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)*/

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
USE MOD_Globals                 ,ONLY: ABORT, OrthoNormVec,VECNORM
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC,CollisMode
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn,DSMC, useDSMC, RadialWeighting, AmbipolElecVelo
USE MOD_DSMC_Vars               ,ONLY: PolyatomMolDSMC, VibQuantsPar
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Part_Tools              ,ONLY: GetParticleWeight
USE MOD_Particle_Boundary_Vars  ,ONLY: dXiEQ_SurfSample,CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact, GetWallTemperature,CalcRotWallVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState,GlobalSide2SurfSide
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Surfaces       ,ONLY: CalcNormAndTangTriangle,CalcNormAndTangBilinear,CalcNormAndTangBezier
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,Species,PartSpecies,nSpecies,WriteMacroSurfaceValues,Symmetry
USE MOD_Particle_Vars           ,ONLY: VarTimeStep, usevMPF
USE MOD_TimeDisc_Vars           ,ONLY: dt,tend,time,RKdtFrac
USE MOD_DSMC_ElectronicModel    ,ONLY: RelaxElectronicShellWall
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
#if (PP_TimeDiscMethod==400)
USE MOD_BGK_Vars                ,ONLY: BGKDoVibRelaxation
#elif (PP_TimeDiscMethod==300)
USE MOD_FPFlow_Vars             ,ONLY: FPDoVibRelaxation
#endif
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
REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC, ElecACC
REAL                                 :: tang1(1:3),tang2(1:3), NewVelo(3), POI_vec(1:3)
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
INTEGER                             :: LocSideID, CNElemID
LOGICAL                             :: DoSample
REAL                                :: EvibOld,ErotOld
REAL                                :: EtraOldAmbi, EtraNewAmbi, EtraWallAmbi, NewVeloAmbi(3), VeloCxAmbi, VeloCyAmbi, VeloCzAmbi
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
  ElecACC    = PartAuxBC%ElecACC(AuxBCIdx)
ELSE
  ! additional states
  locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  ! get BC values
  WallVelo   = PartBound%WallVelo(1:3,locBCID)
  WallTemp   = GetWallTemperature(PartID,locBCID,PartTrajectory,alpha)
  TransACC   = PartBound%TransACC(locBCID)
  VibACC     = PartBound%VibACC(locBCID)
  RotACC     = PartBound%RotACC(locBCID)
  ElecACC    = PartBound%ElecACC(locBCID)
END IF !IsAuxBC

POI_vec(1:3) = LastPartPos(1:3,PartID) + PartTrajectory(1:3)*alpha

IF(PartBound%RotVelo(locBCID)) THEN
  CALL CalcRotWallVelo(locBCID,PartID,POI_vec,WallVelo)
END IF


CALL OrthoNormVec(n_loc,tang1,tang2)

IF(Symmetry%Axisymmetric) THEN
  ! Storing the old and the new particle position (which is outside the domain), at this point the position is only in the xy-plane
  VelX = PartState(1,PartID) - LastPartPos(1,PartID)
  VelY = PartState(2,PartID) - LastPartPos(2,PartID)
  VelZ = PartState(3,PartID) - LastPartPos(3,PartID)

  CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,SideID))
  LocSideID = SideInfo_Shared(SIDE_LOCALID,SideID)

  ! Getting the vectors, which span the cell (1-2 and 1-4)
  Vector1(1:3)=NodeCoords_Shared(1:3,ElemSideNodeID_Shared(2,LocSideID,CNElemID)+1)-NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)
  Vector2(1:3)=NodeCoords_Shared(1:3,ElemSideNodeID_Shared(4,LocSideID,CNElemID)+1)-NodeCoords_Shared(1:3,ElemSideNodeID_Shared(1,LocSideID,CNElemID)+1)

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

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
    VeloReal = VECNORM(AmbipolElecVelo(PartID)%ElecVelo(1:3))
    EtraOldAmbi = 0.5 * Species(DSMC%AmbiDiffElecSpec)%MassIC * VeloReal**2
    CALL RANDOM_NUMBER(RanNum)
    VeloCrad    = SQRT(-LOG(RanNum))
    CALL RANDOM_NUMBER(RanNum)
    VeloCzAmbi      = SQRT(-LOG(RanNum))
    Fak_D       = VeloCrad**2 + VeloCzAmbi**2
    EtraWallAmbi= BoltzmannConst * WallTemp * Fak_D
    EtraNewAmbi = EtraOldAmbi + TransACC * (EtraWallAmbi - EtraOldAmbi)
    Cmr         = SQRT(2.0 * EtraNewAmbi / (Species(DSMC%AmbiDiffElecSpec)%MassIC * Fak_D))
    CALL RANDOM_NUMBER(RanNum)
    Phi     = 2.0 * PI * RanNum
    VeloCxAmbi  = Cmr * VeloCrad * COS(Phi) ! tang1
    VeloCyAmbi  = Cmr * VeloCrad * SIN(Phi) ! tang2
    VeloCzAmbi  = Cmr * VeloCzAmbi
  END IF
END IF

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
    SurfSideID=GlobalSide2SurfSide(SURF_SIDEID,SideID)
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

    SampWallState(1,p,q,SurfSideID) = SampWallState(1,p,q,SurfSideID) + EtraOld  * MacroParticleFactor
    SampWallState(2,p,q,SurfSideID) = SampWallState(2,p,q,SurfSideID) + EtraWall * MacroParticleFactor
    SampWallState(3,p,q,SurfSideID) = SampWallState(3,p,q,SurfSideID) + EtraNew  * MacroParticleFactor
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
        SampWallState(1,p,q,SurfSideID) = SampWallState(1,p,q,SurfSideID) + EtraOldAmbi  * MacroParticleFactor
        SampWallState(2,p,q,SurfSideID) = SampWallState(2,p,q,SurfSideID) + EtraWallAmbi * MacroParticleFactor
        SampWallState(3,p,q,SurfSideID) = SampWallState(3,p,q,SurfSideID) + EtraNewAmbi  * MacroParticleFactor
      END IF
    END IF

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
! add wall velocity
NewVelo(1:3) = NewVelo(1:3) + WallVelo(1:3)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
    IF(Symmetry%Axisymmetric) THEN
      NewVeloAmbi(1) = VecX*VeloCxAmbi + nx*VeloCzAmbi
      NewVeloAmbi(2) = VecY*VeloCxAmbi + ny*VeloCzAmbi
      NewVeloAmbi(3) = VeloCyAmbi
    ELSE
      NewVeloAmbi(1:3) = VeloCxAmbi*tang1(1:3)-tang2(1:3)*VeloCyAmbi-VeloCzAmbi*n_loc(1:3)
    END IF
    NewVeloAmbi(1:3) = NewVeloAmbi(1:3) + WallVelo(1:3)
  END IF
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
          SampWallState(4,p,q,SurfSideID)=SampWallState(4,p,q,SurfSideID)+PartStateIntEn(2,PartID) * MacroParticleFactor
          SampWallState(5,p,q,SurfSideID)=SampWallState(5,p,q,SurfSideID)+ErotWall * MacroParticleFactor
          SampWallState(6,p,q,SurfSideID)=SampWallState(6,p,q,SurfSideID)+ErotNew * MacroParticleFactor
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
                SampWallState(7,p,q,SurfSideID)= SampWallState(7,p,q,SurfSideID) &
                    + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * MacroParticleFactor
                SampWallState(8,p,q,SurfSideID)= SampWallState(8,p,q,SurfSideID) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                    * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * MacroParticleFactor
              END DO
            ELSE
              SampWallState(7,p,q,SurfSideID)= SampWallState(7,p,q,SurfSideID) + (VibQuant + DSMC%GammaQuant) &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * MacroParticleFactor
              SampWallState(8,p,q,SurfSideID)= SampWallState(8,p,q,SurfSideID) + VibQuantWall &
                  * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib * MacroParticleFactor
            END IF
            SampWallState(9,p,q,SurfSideID)= SampWallState(9,p,q,SurfSideID) + EvibNew * MacroParticleFactor
            ! #endif
          END IF
          IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) VibQuantsPar(PartID)%Quants(:) = VibQuantTemp(:)
          PartStateIntEn(1,PartID) = EvibNew
#if (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==300)
        END IF ! FPDoVibRelaxation || BGKDoVibRelaxation
#endif
      END IF

      IF ( DSMC%ElectronicModel ) THEN
        IF((SpecDSMC(PartSpecies(PartID))%InterID.NE.4).AND.(.NOT.SpecDSMC(PartSpecies(PartID))%FullyIonized)) THEN
          CALL RANDOM_NUMBER(RanNum)
          IF (RanNum.LT.ElecACC) THEN
            IF (DoSample) SampWallState(4,p,q,SurfSideID)=SampWallState(4,p,q,SurfSideID)+PartStateIntEn(3,PartID) * MacroParticleFactor
            PartStateIntEn(3,PartID) = RelaxElectronicShellWall(PartID, WallTemp)
            IF (DoSample) THEN
              SampWallState(5,p,q,SurfSideID)=SampWallState(5,p,q,SurfSideID)+PartStateIntEn(3,PartID) * MacroParticleFactor
              SampWallState(6,p,q,SurfSideID)=SampWallState(6,p,q,SurfSideID)+PartStateIntEn(3,PartID) * MacroParticleFactor
            END IF
          END IF
        END IF
      END IF
    END IF ! CollisMode > 1
  END IF ! useDSMC

  ! Sampling of the time step at the wall to get the correct time sample duration for the surface values
  IF (VarTimeStep%UseVariableTimeStep) THEN
    adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
    IF (DoSample) THEN
      SampWallState(12+nSpecies+1,p,q,SurfSideID) = SampWallState(12+nSpecies+1,p,q,SurfSideID) + adaptTimeStep
    END IF
  ELSE
    adaptTimeStep = 1.
  END IF ! VarTimeStep%UseVariableTimeStep

  ! intersection point with surface
  LastPartPos(1:3,PartID) = POI_vec(1:3)

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
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) THEN
        rotVelY = (NewVeloAmbi(2)*(PartState(2,PartID))+NewVeloAmbi(3)*PartState(3,PartID))/rotPosY
        rotVelZ = (-NewVeloAmbi(2)*PartState(3,PartID)+NewVeloAmbi(3)*(PartState(2,PartID)))/rotPosY

        NewVeloAmbi(2) = rotVelY
        NewVeloAmbi(3) = rotVelZ
      END IF
    END IF
    rotVelY = (NewVelo(2)*(PartState(2,PartID))+NewVelo(3)*PartState(3,PartID))/rotPosY
    rotVelZ = (-NewVelo(2)*PartState(3,PartID)+NewVelo(3)*(PartState(2,PartID)))/rotPosY

    PartState(2,PartID) = rotPosY
    PartState(3,PartID) = 0.0
    NewVelo(2) = rotVelY
    NewVelo(3) = rotVelZ
  END IF ! Symmetry%Axisymmetric

  IF(Symmetry%Order.LT.3) THEN
    ! y/z-variable is set to zero for the different symmetry cases
    LastPartPos(Symmetry%Order+1:3,PartID) = 0.0
    PartState(Symmetry%Order+1:3,PartID) = 0.0
  END IF

  IF (DoSample) THEN
    !----  Sampling force at walls
    SampWallState(10:12,p,q,SurfSideID)= SampWallState(10:12,p,q,SurfSideID) &
        + Species(PartSpecies(PartID))%MassIC * (PartState(4:6,PartID) - NewVelo(1:3)) * MacroParticleFactor
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0)  SampWallState(10:12,p,q,SurfSideID)= SampWallState(10:12,p,q,SurfSideID) &
        + Species(DSMC%AmbiDiffElecSpec)%MassIC * (AmbipolElecVelo(PartID)%ElecVelo(1:3) - NewVeloAmbi(1:3)) * MacroParticleFactor
    END IF
    !---- Counter for collisions (normal wall collisions - not to count if only SpeciesSwaps to be counted)
    IF (.NOT.CalcSurfCollis%OnlySwaps .AND. .NOT.IsSpeciesSwap) THEN
      SampWallState(12+PartSpecies(PartID),p,q,SurfSideID)= SampWallState(12+PartSpecies(PartID),p,q,SurfSideID) +1
      IF (DSMC%DoAmbipolarDiff) THEN
        IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0)  SampWallState(12+DSMC%AmbiDiffElecSpec,p,q,SurfSideID)= &
          SampWallState(12+DSMC%AmbiDiffElecSpec,p,q,SurfSideID) +1
      END IF
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
PartState(4:6,PartID)   = NewVelo(1:3)
IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(PartSpecies(PartID))%ChargeIC.GT.0.0) AmbipolElecVelo(PartID)%ElecVelo(1:3) = NewVeloAmbi(1:3) 
END IF

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

#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
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
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,dXiEQ_SurfSample,CalcSurfCollis,AnalyzeSurfCollis,PartAuxBC
USE MOD_Particle_Boundary_Vars  ,ONLY: SampWallState,GlobalSide2SurfSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,PartSpecies,usevMPF
USE MOD_Particle_Vars           ,ONLY: WriteMacroSurfaceValues,nSpecies,Species
USE MOD_DSMC_Vars               ,ONLY: DSMC, RadialWeighting
USE MOD_TimeDisc_Vars           ,ONLY: TEnd,Time
USE MOD_Particle_Boundary_Vars  ,ONLY: CalcSurfaceImpact
USE MOD_Particle_Boundary_Tools ,ONLY: CountSurfaceImpact
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn
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
  locBCID = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  CALL RANDOM_NUMBER(RanNum)
  IF(RanNum.LE.PartBound%ProbOfSpeciesSwaps(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
    targetSpecies=-1 ! Dummy initialization value
    IF(PRESENT(targetSpecies_IN))THEN
      targetSpecies = targetSpecies_IN
    ELSE ! Normal swap routine
      DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
        IF (PartSpecies(PartID).eq.PartBound%SpeciesSwaps(1,iSwaps,PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) &
            targetSpecies = PartBound%SpeciesSwaps(2,iSwaps,PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
      END DO
    END IF ! PRESENT(targetSpecies_IN)
    !swap species
    IF (targetSpecies.ge.0) IsSpeciesSwap=.TRUE.
    IF ( (targetSpecies.eq.0) .OR. (.NOT.CalcSurfCollis%Only0Swaps) ) THEN
      IF (DoSample) THEN
        !---- Counter for swap species collisions
        SurfSideID=GlobalSide2SurfSide(SURF_SIDEID,SideID)
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

        SampWallState(12+PartSpecies(PartID),p,q,SurfSideID) = SampWallState(12+PartSpecies(PartID),p,q,SurfSideID) + 1
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
      ! sample values of deleted species
      IF (DoSample) THEN
        SurfSideID=GlobalSide2SurfSide(SURF_SIDEID,SideID)
        IF (TrackingMethod.EQ.TRIATRACKING) THEN
          p=1 ; q=1
        ELSE
          Xitild =MIN(MAX(-1.,xi ),0.99)
          Etatild=MIN(MAX(-1.,eta),0.99)
          p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
          q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1
        END IF
        !----  Sampling Forces at walls
        SampWallState(10:12,p,q,SurfSideID)= SampWallState(10:12,p,q,SurfSideID) + Species(PartSpecies(PartID))%MassIC &
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
      CALL RemoveParticle(PartID,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)),alpha=alpha)
    ELSE IF (targetSpecies.gt.0) THEN !swap species
      PartSpecies(PartID)=targetSpecies
    END IF ! targetSpecies.eq.0
  END IF ! RanNum.LE.PartBound%ProbOfSpeciesSwaps
END IF ! IsAuxBC

END SUBROUTINE SpeciesSwap


SUBROUTINE SurfaceFluxBasedBoundaryTreatment(iPart,SideID,alpha,PartTrajectory)
!===================================================================================================================================
! Treatment of particles at the boundary if adaptive surface BCs or circular inflows based on the surface flux are present
! Circular Inflow: Particles are deleted if within (allows multiple surface flux inflows defined by circles on a single boundary)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species, LastPartPos, PartSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
!USE MOD_Mesh_Vars              ,ONLY: BC
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
  IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))) THEN
    IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      intersectionPoint(1:3) = LastPartPos(1:3,iPart)+ alpha*PartTrajectory(1:3)
      point(1)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(2))-Species(iSpec)%Surfaceflux(iSF)%origin(1)
      point(2)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(3))-Species(iSpec)%Surfaceflux(iSF)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
        CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)),alpha=alpha)
      END IF
    END IF
  END IF
END DO

END SUBROUTINE SurfaceFluxBasedBoundaryTreatment


END MODULE MOD_SurfaceModel
