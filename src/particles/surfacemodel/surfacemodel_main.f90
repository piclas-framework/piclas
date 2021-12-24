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
PUBLIC :: SurfaceModel, MaxwellScattering, PerfectReflection, DiffuseReflection, SpeciesSwap
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceModel(PartID,SideID,GlobElemID,n_Loc)
!===================================================================================================================================
!> Selection and execution of a gas-surface interaction model
!> 1.) Initial surface pre-treatment: Porous BC, species swap and charge deposition on dielectrics
!> 2.) Count and sample the properties BEFORE the surface interaction
!> 3.) Perform the selected gas-surface interaction, currently implemented models:
!           0: Maxwell Scattering
!       5/6/7: Secondary Electron Emission
!> 4.) Count and sample the properties AFTER the surface interaction
!===================================================================================================================================
!USE MOD_Globals           ,ONLY: myrank,ierror,mpi_comm_world
USE MOD_Globals                   ,ONLY: abort,UNITVECTOR,OrthoNormVec
USE MOD_Particle_Vars             ,ONLY: PartSpecies, WriteMacroSurfaceValues
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackingMethod, TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: Partbound, GlobalSide2SurfSide, dXiEQ_SurfSample, PartBound
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars             ,ONLY: PDM, LastPartPos
USE MOD_Particle_Vars             ,ONLY: UseCircularInflow
USE MOD_Dielectric_Vars           ,ONLY: DoDielectricSurfaceCharge
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcSurfCollCounter, SurfAnalyzeCount, SurfAnalyzeNumOfAds, SurfAnalyzeNumOfDes
USE MOD_SurfaceModel_Tools        ,ONLY: SurfaceModel_ParticleEmission
USE MOD_SEE                       ,ONLY: SecondaryElectronEmission
USE MOD_SurfaceModel_Porous       ,ONLY: PorousBoundaryTreatment
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_PICDepo_Tools             ,ONLY: DepositParticleOnNodes
USE MOD_part_operations           ,ONLY: RemoveParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)
INTEGER,INTENT(IN) :: PartID, SideID
INTEGER,INTENT(IN) :: GlobElemID  !< Global element ID of the particle impacting the surface
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ProductSpec(1:2) !< 1: product species of incident particle (also used for simple reflection)
                                       !< 2: additional species added or removed from surface
                                       !< If productSpec is negative, then the respective particles are adsorbed
                                       !< If productSpec is positive the particle is reflected/emitted
                                       !< with respective species
INTEGER            :: ProductSpecNbr   !< number of emitted particles for ProductSpec(2)
REAL               :: TempErgy(2)      !< temperature, energy or velocity used for VeloFromDistribution
REAL               :: Xitild,Etatild
INTEGER            :: SpecID, locBCID
INTEGER            :: iBC, SurfSideID
LOGICAL            :: SpecularReflectionOnly,DoSample
!===================================================================================================================================
iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))

!===================================================================================================================================
! 1.) Initial surface pre-treatment
!===================================================================================================================================
!---- Treatment of adaptive and porous boundary conditions (deletion of particles in case of circular inflow or porous BC)
SpecularReflectionOnly = .FALSE.
IF(UseCircularInflow) CALL SurfaceFluxBasedBoundaryTreatment(PartID,SideID)
IF(nPorousBC.GT.0) CALL PorousBoundaryTreatment(PartID,SideID,SpecularReflectionOnly)

!---- Dielectric particle-surface interaction
IF(DoDielectricSurfaceCharge.AND.PartBound%Dielectric(iBC)) THEN
  CALL DepositParticleOnNodes(PartID,LastPartPos(1:3,PartID)+TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha,GlobElemID)
END IF

!---- swap species?
IF (PartBound%NbrOfSpeciesSwaps(iBC).gt.0) THEN
  CALL SpeciesSwap(PartID,SideID)
END IF

locBCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SpecID = PartSpecies(PartID)
ProductSpec(1) = SpecID
ProductSpec(2) = 0
ProductSpecNbr = 0
TempErgy(1:2)=PartBound%WallTemp(locBCID)
!===================================================================================================================================
! 2.) Count and sample the properties BEFORE the surface interaction
!===================================================================================================================================
! Counter for surface analyze
IF(CalcSurfCollCounter) SurfAnalyzeCount(SpecID) = SurfAnalyzeCount(SpecID) + 1
! Sampling
DoSample = (DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)
IF(DoSample) THEN
  IF (TrackingMethod.EQ.TRIATRACKING) THEN
    TrackInfo%p = 1 ; TrackInfo%q = 1
  ELSE
    Xitild  = MIN(MAX(-1.,TrackInfo%xi ),0.99)
    Etatild = MIN(MAX(-1.,TrackInfo%eta),0.99)
    TrackInfo%p = INT((Xitild +1.0)/dXiEQ_SurfSample)+1
    TrackInfo%q = INT((Etatild+1.0)/dXiEQ_SurfSample)+1
  END IF
  SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
  ! Sample momentum, heatflux and collision counter on surface
  CALL CalcWallSample(PartID,SurfSideID,'old',SurfaceNormal_opt=n_loc)
END IF
!===================================================================================================================================
! Particle was deleted during the species swap/circular flow, leave the routine
!===================================================================================================================================
IF(.NOT.PDM%ParticleInside(PartID)) THEN
  ! Increase the counter for deleted/absorbed/adsorbed particles
  IF(CalcSurfCollCounter) SurfAnalyzeNumOfAds(SpecID) = SurfAnalyzeNumOfAds(SpecID) + 1
  RETURN
END IF
!===================================================================================================================================
! 3.) Perform the selected gas-surface interaction
!===================================================================================================================================
SELECT CASE(PartBound%SurfaceModel(locBCID))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (0) ! Maxwellian scattering (diffuse/specular reflection)
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL MaxwellScattering(PartID,SideID,n_Loc,SpecularReflectionOnly)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (SEE_MODELS_ID)
  ! 5: SEE by Levko2015
  ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
  ! 8: SEE-E (bombarding electrons are reflected, e- on dielectric materials is considered for SEE and three different outcomes)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Get electron emission probability
  CALL SecondaryElectronEmission(PartID,locBCID,ProductSpec,ProductSpecNbr,TempErgy(2))
  !IF(myrank.eq.0) read*; CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
  ! Decide the fate of the impacting particle
  IF (ProductSpec(1).LE.0) THEN
    CALL RemoveParticle(PartID)
  ELSE
    CALL MaxwellScattering(PartID,SideID,n_Loc)
  END IF
  ! Emit the secondary electrons
  IF (ProductSpec(2).GT.0) THEN
    CALL SurfaceModel_ParticleEmission(n_loc, PartID, SideID, ProductSpec, ProductSpecNbr, TempErgy, GlobElemID)
  END IF
CASE DEFAULT
  CALL abort(__STAMP__,'Unknown surface model. PartBound%SurfaceModel(locBCID) = ',IntInfoOpt=PartBound%SurfaceModel(locBCID))
END SELECT

!===================================================================================================================================
! 4.) Count and sample the properties AFTER the surface interaction
!===================================================================================================================================
! Counter for surface analyze
IF(CalcSurfCollCounter) THEN
  ! Old particle was deleted/absorbed/adsorbed at the wall (in case it happend during a surface model and not during the 
  ! SpeciesSwap routine)
  IF (ProductSpec(1).LE.0) THEN
    SurfAnalyzeNumOfAds(SpecID) = SurfAnalyzeNumOfAds(SpecID) + 1
  END IF
  ! New particle was created/desorbed at the wall
  IF (ProductSpec(2).GT.0) THEN
    SurfAnalyzeNumOfDes(ProductSpec(2)) = SurfAnalyzeNumOfDes(ProductSpec(2)) + ProductSpecNbr
  END IF
END IF
! Sampling
IF(DoSample) THEN
  ! Sample momentum, heatflux and collision counter on surface (only the original impacting particle, not the newly created parts
  ! through SurfaceModel_ParticleEmission)
  CALL CalcWallSample(PartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
END IF

END SUBROUTINE SurfaceModel


SUBROUTINE MaxwellScattering(PartID,SideID,n_loc,SpecularReflectionOnly_opt)
!===================================================================================================================================
!> SurfaceModel = 0, classic DSMC gas-surface interaction model choosing between a perfect specular and a complete diffuse
!> reflection by comparing the given momentum accommodation coefficient (MomentumACC) with a random number
!===================================================================================================================================
USE MOD_Particle_Boundary_Vars  ,ONLY: Partbound
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: n_loc(1:3)
INTEGER,INTENT(IN)          :: PartID, SideID
LOGICAL,INTENT(IN),OPTIONAL :: SpecularReflectionOnly_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: RanNum,ACC
INTEGER :: iBC
!===================================================================================================================================
iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
ACC = PartBound%MomentumACC(iBC)

! Check if optional parameter was supplied
ASSOCIATE( SpecularReflectionOnly => MERGE(SpecularReflectionOnly_opt, .FALSE., PRESENT(SpecularReflectionOnly_opt)) ) 
  IF (SpecularReflectionOnly.OR.ACC.EQ.0.0) THEN
    CALL PerfectReflection(PartID,SideID,n_loc)
  ELSE IF(ACC.LT.1.0) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF(RanNum.GE.ACC) THEN
      CALL PerfectReflection(PartID,SideID,n_loc)
    ELSE
      CALL DiffuseReflection(PartID,SideID,n_loc)
    END IF
  ELSE
    CALL DiffuseReflection(PartID,SideID,n_loc)
  END IF
END ASSOCIATE

END SUBROUTINE MaxwellScattering


SUBROUTINE PerfectReflection(PartID,SideID,n_Loc,opt_Symmetry,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the perfect reflection in 3D
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,PartAuxBC
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,PartSpecies,Species,PartLorentzType
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
USE MOD_Globals_Vars            ,ONLY: c2_inv
#if defined(LSERK)
USE MOD_Particle_Vars           ,ONLY: Pt_temp,PDM
#elif (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Vars           ,ONLY: PEM
#endif
USE MOD_SurfaceModel_Tools      ,ONLY: CalcRotWallVelo
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID !,ElemID
LOGICAL,INTENT(IN),OPTIONAL       :: opt_Symmetry
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: v_old(1:3),WallVelo(3), v_old_Ambi(1:3)
REAL                                 :: LorentzFac, LorentzFacInv
INTEGER                              :: locBCID
LOGICAL                              :: Symmetry, IsAuxBC
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
POI_vec(1:3) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha

!IF(PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))) THEN
!  CALL CalcRotWallVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)),POI_vec,WallVelo)
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

! set particle position on face
LastPartPos(1:3,PartID) = POI_vec(1:3)

TrackInfo%PartTrajectory(1:3)     = TrackInfo%PartTrajectory(1:3)-2.*DOT_PRODUCT(TrackInfo%PartTrajectory(1:3),n_loc)*n_loc

! Check if rotated velo is used, otherwise mirror the LastPartPos
!IF(PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))))THEN
!  PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac &
!                          * PartState(4:6,PartID) * adaptTimeStep
!ELSE
  PartState(1:3,PartID) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*(TrackInfo%lengthPartTrajectory - TrackInfo%alpha)
!END IF ! PartBound%RotVelo(PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))

! #if !defined(IMPA) &&  !defined(ROS)
! compute moved particle || rest of movement
TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)

TrackInfo%lengthPartTrajectory = VECNORM(TrackInfo%PartTrajectory)
IF(ALMOSTZERO(TrackInfo%lengthPartTrajectory)) THEN
  TrackInfo%lengthPartTrajectory= 0.0
ELSE
  TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory
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


SUBROUTINE DiffuseReflection(PartID,SideID,n_loc,AuxBCIdx)
!----------------------------------------------------------------------------------------------------------------------------------!
!> Computes the new particle state (position, velocity, and energy) after a diffuse reflection
!> 1.) Get the wall velocity, temperature and accommodation coefficients
!> 2.) Get the tangential vectors
!> 3.) Calculate new velocity vector (Extended Maxwellian Model)+
!> 4.) Perform vector transformation from the local to the global coordinate system and add wall velocity
!> 5.) Perform internal energy accommodation at the wall
!> 6.) Determine the new particle position after the reflection
!> 7.) Axisymmetric simulation: Rotate the vector back into the symmetry plane
!> 8.) Saving new particle velocity and recompute the trajectory based on new and old particle position
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_Globals                 ,ONLY: ABORT, OrthoNormVec, VECNORM, DOTPRODUCT
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
USE MOD_SurfaceModel_Tools      ,ONLY: GetWallTemperature, CalcRotWallVelo
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,PartAuxBC
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,Species,PartSpecies,Symmetry
USE MOD_Particle_Vars           ,ONLY: VarTimeStep
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars           ,ONLY: PDM
#endif
USE MOD_SurfaceModel_Tools      ,ONLY: CalcPostWallCollVelo, SurfaceModel_EnergyAccommodation
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
REAL,INTENT(IN)                   :: n_loc(1:3)
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: LocSideID, CNElemID, locBCID, SpecID
REAL                              :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC, ElecACC
REAL                              :: tang1(1:3), tang2(1:3), NewVelo(3), POI_vec(1:3), NewVeloAmbi(3), VeloC(1:3), VeloCAmbi(1:3)
REAL                              :: POI_fak, TildTrajectory(3), adaptTimeStep
LOGICAL                           :: IsAuxBC
! Symmetry
REAL                              :: rotVelY, rotVelZ, rotPosY
REAL                              :: nx, ny, nVal, VelX, VelY, VecX, VecY, Vector1(1:3), Vector2(1:3)
!===================================================================================================================================
IF (PRESENT(AuxBCIdx)) THEN
  IsAuxBC=.TRUE.
ELSE
  IsAuxBC=.FALSE.
END IF

! 1.) Get the wall velocity, temperature and accommodation coefficients
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
  WallTemp   = GetWallTemperature(PartID,locBCID,SideID)
  TransACC   = PartBound%TransACC(locBCID)
  VibACC     = PartBound%VibACC(locBCID)
  RotACC     = PartBound%RotACC(locBCID)
  ElecACC    = PartBound%ElecACC(locBCID)
END IF !IsAuxBC

SpecID = PartSpecies(PartID)

POI_vec(1:3) = LastPartPos(1:3,PartID) + TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha

IF(PartBound%RotVelo(locBCID)) THEN
  CALL CalcRotWallVelo(locBCID,POI_vec,WallVelo)
END IF

! 2.) Get the tangential vectors
IF(Symmetry%Axisymmetric) THEN
  ! Storing the old and the new particle position (which is outside the domain), at this point the position is only in the xy-plane
  VelX = PartState(1,PartID) - LastPartPos(1,PartID)
  VelY = PartState(2,PartID) - LastPartPos(2,PartID)

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
  ! Check for the correct orientation of the normal vectors (should be inwards)
  IF ((VelX*nx+VelY*ny).GT.0) THEN
    nx = -Vector1(2)
    ny = Vector1(1)
  END IF

  nVal = SQRT(nx*nx + ny*ny)
  nx = nx/nVal
  ny = ny/nVal
ELSE
  CALL OrthoNormVec(n_loc,tang1,tang2)
END IF

! 3.) Calculate new velocity vector (Extended Maxwellian Model)
VeloC(1:3) = CalcPostWallCollVelo(SpecID,DOTPRODUCT(PartState(4:6,PartID)),WallTemp,TransACC)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    VeloCAmbi(1:3) = CalcPostWallCollVelo(DSMC%AmbiDiffElecSpec,DOTPRODUCT(AmbipolElecVelo(PartID)%ElecVelo(1:3)),WallTemp,TransACC)
  END IF
END IF

! 4.) Perform vector transformation from the local to the global coordinate system and add wall velocity
!     NewVelo = VeloCx*tang1+CROSS(-n_loc,tang1)*VeloCy-VeloCz*n_loc
IF(Symmetry%Axisymmetric) THEN
  VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2)
  NewVelo(1) = VecX*VeloC(1) + nx*VeloC(3)
  NewVelo(2) = VecY*VeloC(1) + ny*VeloC(3)
  NewVelo(3) = VeloC(2)
ELSE
  NewVelo(1:3) = VeloC(1)*tang1(1:3)-tang2(1:3)*VeloC(2)-VeloC(3)*n_loc(1:3)
END IF

NewVelo(1:3) = NewVelo(1:3) + WallVelo(1:3)

IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
    IF(Symmetry%Axisymmetric) THEN
      NewVeloAmbi(1) = VecX*VeloCAmbi(1) + nx*VeloCAmbi(3)
      NewVeloAmbi(2) = VecY*VeloCAmbi(1) + ny*VeloCAmbi(3)
      NewVeloAmbi(3) = VeloCAmbi(2)
    ELSE
      NewVeloAmbi(1:3) = VeloCAmbi(1)*tang1(1:3)-tang2(1:3)*VeloCAmbi(2)-VeloCAmbi(3)*n_loc(1:3)
    END IF
    NewVeloAmbi(1:3) = NewVeloAmbi(1:3) + WallVelo(1:3)
  END IF
END IF

! 5.) Perform internal energy accommodation at the wall
IF (.NOT.IsAuxBC) THEN !so far no internal DOF stuff for AuxBC!!!
  CALL SurfaceModel_EnergyAccommodation(PartID,locBCID,WallTemp)


! 6.) Determine the new particle position after the reflection
  LastPartPos(1:3,PartID) = POI_vec(1:3)

  IF (VarTimeStep%UseVariableTimeStep) THEN
    adaptTimeStep = VarTimeStep%ParticleTimeStep(PartID)
  ELSE
    adaptTimeStep = 1.
  END IF ! VarTimeStep%UseVariableTimeStep
  ! recompute initial position and ignoring preceding reflections and trajectory between current position and recomputed position
  !TildPos       =PartState(1:3,PartID)-dt*RKdtFrac*PartState(4:6,PartID)
  TildTrajectory=dt*RKdtFrac*PartState(4:6,PartID)*adaptTimeStep
  POI_fak=1.- (TrackInfo%lengthPartTrajectory-TrackInfo%alpha)/SQRT(DOT_PRODUCT(TildTrajectory,TildTrajectory))
  ! travel rest of particle vector
  !PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - alpha/lengthPartTrajectory) * dt*RKdtFrac * NewVelo(1:3)
  IF (IsAuxBC) THEN
    IF (PartAuxBC%Resample(AuxBCIdx)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  ELSE
    IF (PartBound%Resample(locBCID)) CALL RANDOM_NUMBER(POI_fak) !Resample Equilibirum Distribution
  END IF ! IsAuxBC
  PartState(1:3,PartID)   = LastPartPos(1:3,PartID) + (1.0 - POI_fak) * dt*RKdtFrac * NewVelo(1:3) * adaptTimeStep

! 7.) Axisymmetric simulation: Rotate the vector back into the symmetry plane
  IF(Symmetry%Axisymmetric) THEN
    ! Symmetry considerations --------------------------------------------------------
    rotPosY = SQRT(PartState(2,PartID)**2 + (PartState(3,PartID))**2)
    ! Rotation: Vy' =   Vy * cos(alpha) + Vz * sin(alpha) =   Vy * y/y' + Vz * z/y'
    !           Vz' = - Vy * sin(alpha) + Vz * cos(alpha) = - Vy * z/y' + Vz * y/y'
    ! Right-hand system, using new y and z positions after tracking, position vector and velocity vector DO NOT have to
    ! coincide (as opposed to Bird 1994, p. 391, where new positions are calculated with the velocity vector)
    IF (DSMC%DoAmbipolarDiff) THEN
      IF(Species(SpecID)%ChargeIC.GT.0.0) THEN
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
END IF !.NOT.IsAuxBC

! 8.) Saving new particle velocity and recompute the trajectory based on new and old particle position
PartState(4:6,PartID)   = NewVelo(1:3)
IF (DSMC%DoAmbipolarDiff) THEN
  IF(Species(SpecID)%ChargeIC.GT.0.0) AmbipolElecVelo(PartID)%ElecVelo(1:3) = NewVeloAmbi(1:3)
END IF

! Recompute trajectory etc
IF(Symmetry%Axisymmetric) THEN
  TrackInfo%PartTrajectory(1:2)=PartState(1:2,PartID) - LastPartPos(1:2,PartID)
  TrackInfo%PartTrajectory(3) = 0.
  TrackInfo%lengthPartTrajectory=SQRT(TrackInfo%PartTrajectory(1)**2 + TrackInfo%PartTrajectory(2)**2)
ELSE
  TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)
  TrackInfo%lengthPartTrajectory=VECNORM(TrackInfo%PartTrajectory(1:3))
END IF
IF(ABS(TrackInfo%lengthPartTrajectory).GT.0.) TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory

#if defined(LSERK) || (PP_TimeDiscMethod==508) || (PP_TimeDiscMethod==509)
PDM%IsNewPart(PartID)=.TRUE. !reconstruction in timedisc during push
#endif

END SUBROUTINE DiffuseReflection


SUBROUTINE SpeciesSwap(PartID,SideID,AuxBCIdx,targetSpecies_IN)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the Species Swap on ReflectiveBC
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                 ,ONLY: abort,VECNORM
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound,PartAuxBC
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars           ,ONLY: PartSpecies
USE MOD_part_operations         ,ONLY: RemoveParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: AuxBCIdx
INTEGER,INTENT(IN),OPTIONAL       :: targetSpecies_IN
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: targetSpecies, iSwaps
REAL                              :: RanNum
INTEGER                           :: locBCID
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
    IF (targetSpecies.eq.0) THEN !delete particle -> same as PartAuxBC%OpenBC
      CALL RemoveParticle(PartID)
    ELSEIF (targetSpecies.gt.0) THEN !swap species
      PartSpecies(PartID)=targetSpecies
    END IF
  END IF !RanNum.LE.PartAuxBC%ProbOfSpeciesSwaps
ELSE
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
    IF (targetSpecies.eq.0) THEN
      ! Delete particle -> same as PartBound%OpenBC
      CALL RemoveParticle(PartID,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
    ELSE IF (targetSpecies.gt.0) THEN
      ! Swap species
      PartSpecies(PartID)=targetSpecies
    END IF ! targetSpecies.eq.0
  END IF ! RanNum.LE.PartBound%ProbOfSpeciesSwaps
END IF ! IsAuxBC

END SUBROUTINE SpeciesSwap


SUBROUTINE SurfaceFluxBasedBoundaryTreatment(iPart,SideID)
!===================================================================================================================================
! Treatment of particles at the boundary if adaptive surface BCs or circular inflows based on the surface flux are present
! Circular Inflow: Particles are deleted if within (allows multiple surface flux inflows defined by circles on a single boundary)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species, LastPartPos, PartSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_part_operations        ,ONLY: RemoveParticle
USE MOD_Particle_Tracking_Vars ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: iPart, SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: point(1:2), intersectionPoint(1:3), radius
INTEGER                             :: iSpec, iSF
!===================================================================================================================================

iSpec = PartSpecies(iPart)
DO iSF=1,Species(iSpec)%nSurfacefluxBCs
  IF(Species(iSpec)%Surfaceflux(iSF)%BC.EQ.PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))) THEN
    IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      intersectionPoint(1:3) = LastPartPos(1:3,iPart) + TrackInfo%alpha*TrackInfo%PartTrajectory(1:3)
      point(1)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(2))-Species(iSpec)%Surfaceflux(iSF)%origin(1)
      point(2)=intersectionPoint(Species(iSpec)%Surfaceflux(iSF)%dir(3))-Species(iSpec)%Surfaceflux(iSF)%origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
        CALL RemoveParticle(iPart,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
      END IF
    END IF
  END IF
END DO

END SUBROUTINE SurfaceFluxBasedBoundaryTreatment


END MODULE MOD_SurfaceModel
