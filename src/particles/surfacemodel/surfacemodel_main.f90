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
PUBLIC :: SurfaceModelling
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Selection and execution of a gas-surface interaction model
!> 0.) Initial surface pre-treatment: Porous BC and initialization of charge deposition on dielectrics
!> 1.) Count and sample the properties BEFORE the surface interaction
!> 2.) Perform the species swap
!> 3.) Perform the selected gas-surface interaction, currently implemented models:
!           0: Maxwell Scattering
!          20: Adsorption or Eley-Rideal reaction
!       5/6/7/8/9/10/11: Secondary Electron Emission
!> 4.) PIC ONLY: Deposit charges on dielectric surface (when activated), if these were removed/changed in SpeciesSwap or SurfaceModel
!> 5.) Count and sample the properties AFTER the surface interaction
!===================================================================================================================================
SUBROUTINE SurfaceModelling(PartID,SideID,GlobalElemID,n_Loc)
! MODULES
USE MOD_Globals                   ,ONLY: abort,UNITVECTOR,OrthoNormVec
#if USE_MPI
USE MOD_Globals                   ,ONLY: myrank
#endif /*USE_MPI*/
USE MOD_Globals_Vars              ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars             ,ONLY: PartSpecies,WriteMacroSurfaceValues,Species,usevMPF,PartMPF
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackingMethod, TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, GlobalSide2SurfSide, dXiEQ_SurfSample
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC, SurfModEnergyDistribution
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars             ,ONLY: PDM, LastPartPos, PEM
USE MOD_Particle_Vars             ,ONLY: UseCircularInflow
USE MOD_Dielectric_Vars           ,ONLY: DoDielectricSurfaceCharge
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive, DoRadialWeighting, DoLinearWeighting, DoCellLocalWeighting
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcSurfCollCounter, SurfAnalyzeCount, SurfAnalyzeNumOfAds, SurfAnalyzeNumOfDes
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: CalcBoundaryParticleOutput
USE MOD_SurfaceModel_Tools        ,ONLY: MaxwellScattering, SurfaceModelParticleEmission
USE MOD_SurfaceModel_Chemistry    ,ONLY: SurfaceModelChemistry, SurfaceModelEventProbability
USE MOD_SEE                       ,ONLY: SecondaryElectronEmissionYield
USE MOD_SurfaceModel_Porous       ,ONLY: PorousBoundaryTreatment
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_PICDepo_Tools             ,ONLY: DepositParticleOnNodes
USE MOD_part_operations           ,ONLY: RemoveParticle, CreateParticle
USE MOD_part_tools                ,ONLY: CalcRadWeightMPF, CalcVarWeightMPF, VeloFromDistribution, GetParticleWeight
USE MOD_PICDepo_Vars              ,ONLY: DoDeposition
USE MOD_Part_Operations           ,ONLY: UpdateBPO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)
INTEGER,INTENT(IN) :: PartID, SideID
INTEGER,INTENT(IN) :: GlobalElemID  !< Global element ID of the particle impacting the surface
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
REAL               :: TempErgy         !< temperature, energy or velocity used for VeloFromDistribution
REAL               :: Xitild,Etatild
INTEGER            :: PartSpecImpact, locBCID, SurfSideID
LOGICAL            :: SpecularReflectionOnly,DoSample
REAL               :: ChargeImpact,PartPosImpact(1:3) !< Charge and position of impact of bombarding particle
REAL               :: ChargeRefl                      !< Charge of reflected particle
REAL               :: MPF                             !< macro-particle factor
REAL               :: ChargeHole                      !< Charge of SEE electrons holes
INTEGER            :: iProd,iNewPart,NewPartID,iElem
REAL               :: NewVelo(3), NewPos(1:3)
!===================================================================================================================================
!===================================================================================================================================
! 0.) Initial surface pre-treatment
!===================================================================================================================================
!---- Treatment of adaptive and porous boundary conditions (deletion of particles in case of circular inflow or porous BC)
SpecularReflectionOnly = .FALSE.
IF(UseCircularInflow) CALL SurfaceFluxBasedBoundaryTreatment(PartID,SideID)
IF(nPorousBC.GT.0) CALL PorousBoundaryTreatment(PartID,SideID,SpecularReflectionOnly)

locBCID        = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
PartSpecImpact = PartSpecies(PartID)
ProductSpec(1) = PartSpecImpact
ProductSpec(2) = 0
ProductSpecNbr = 0

! Store info of impacting particle for possible surface charging
PartPosImpact(1:3) = LastPartPos(1:3,PartID)+TrackInfo%PartTrajectory(1:3)*TrackInfo%alpha
IF(DoDielectricSurfaceCharge.AND.PartBound%Dielectric(locBCID)) THEN ! Surface charging active + dielectric surface contact
  IF(usevMPF)THEN
    MPF = PartMPF(PartID)
  ELSE
    MPF = Species(PartSpecies(PartID))%MacroParticleFactor
  END IF ! usevMPF
  ChargeImpact = Species(PartSpecies(PartID))%ChargeIC*MPF
END IF
!===================================================================================================================================
! 1.) Count and sample the properties BEFORE the surface interaction
!===================================================================================================================================
! Counter for surface analyze
IF(CalcSurfCollCounter) SurfAnalyzeCount(PartSpecImpact) = SurfAnalyzeCount(PartSpecImpact) + 1
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
  ! Sample momentum, heatflux and collision counter on surface (Check if particle is still inside is required, since particles can
  ! be removed in the case of UseCircularInflow and nPorousBC. These particles shall not be sampled.)
  IF(PDM%ParticleInside(PartID)) CALL CalcWallSample(PartID,SurfSideID,'old',SurfaceNormal_opt=n_loc)
END IF
!===================================================================================================================================
! 2.) Species Swap
!===================================================================================================================================
IF (PartBound%NbrOfSpeciesSwaps(locBCID).GT.0) CALL SpeciesSwap(PartID,SideID)
!===================================================================================================================================
! Particle was deleted during the species swap/circular flow, leave the routine
!===================================================================================================================================
IF(.NOT.PDM%ParticleInside(PartID)) THEN
  ! Increase the counter for deleted/absorbed/adsorbed particles
  IF(CalcSurfCollCounter) SurfAnalyzeNumOfAds(PartSpecImpact) = SurfAnalyzeNumOfAds(PartSpecImpact) + 1
  IF(DoDeposition.AND.DoDielectricSurfaceCharge.AND.PartBound%Dielectric(locBCID)) &
      CALL DepositParticleOnNodes(ChargeImpact, PartPosImpact, GlobalElemID)
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
CASE (1)  ! Sticking coefficient model using tabulated, empirical values
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL StickingCoefficientModel(PartID,SideID,n_Loc)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (2)  ! Event probability model
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL SurfaceModelEventProbability(PartID,SideID,GlobalElemID,n_loc,PartPosImpact(1:3))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(20)  ! Catalytic gas-surface interaction: Adsorption or Eley-Rideal reaction
!-----------------------------------------------------------------------------------------------------------------------------------
  CALL SurfaceModelChemistry(PartID,SideID,GlobalElemID,n_Loc,PartPosImpact(1:3))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (SEE_MODELS_ID)
  ! 3: SEE by square-fit
  ! 4: SEE by power-fit
  ! 5: SEE by Levko2015
  ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
  ! 8: SEE-E (bombarding electrons are reflected, e- on dielectric materials is considered for SEE and three different outcomes)
  ! 9: SEE-I when Ar+ ion bombards surface with 0.01 probability and fixed SEE electron energy of 6.8 eV
  !10: SEE-I (bombarding electrons are removed, Ar+ on copper is considered for SEE)
  !11: SEE-E by e- on quartz (SiO2) is considered
  !12: SEE-E Seiler, H. (1983). Secondary electron emission in the scanning electron microscope.
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Determine the yield and consequently the number of secondaries to be emitted
  CALL SecondaryElectronEmissionYield(PartID,locBCID,ProductSpec,ProductSpecNbr,TempErgy)

  ! Emit the secondary electrons and deposit the opposite charge on the boundary (in case of a dielectric/VDL)
  IF (ProductSpec(2).GT.0) THEN
    CALL SurfaceModelParticleEmission(n_loc, PartID, SideID, ProductSpec(2), ProductSpecNbr, TempErgy, GlobalElemID, &
                                      PartPosImpact(1:3),EnergyDistribution=SurfModEnergyDistribution(locBCID))
    ! Deposit opposite charge of SEE on node
    IF(DoDeposition.AND.DoDielectricSurfaceCharge)THEN

      ! Method 1: PartBound%Dielectric = T (dielectric boundary)
      IF(PartBound%Dielectric(locBCID))THEN
        ! Get MPF
        IF (usevMPF) THEN
          IF (DoRadialWeighting) THEN
            MPF = CalcRadWeightMPF(PartPosImpact(2),ProductSpec(2))
          ELSE IF (DoLinearWeighting.OR.DoCellLocalWeighting) THEN
            iElem = PEM%LocalElemID(PartID)
            MPF = CalcVarWeightMPF(PartPosImpact(:),iElem)
          ELSE
            MPF = Species(ProductSpec(2))%MacroParticleFactor
          END IF
        ELSE
            MPF = Species(ProductSpec(2))%MacroParticleFactor
        END IF ! usevMPF
        ! Calculate the opposite charge
        ChargeHole = -Species(ProductSpec(2))%ChargeIC*MPF
        ! Deposit the charge(s)
        DO iProd = 1, ProductSpecNbr
          CALL DepositParticleOnNodes(ChargeHole, PartPosImpact, GlobalElemID)
        END DO ! iProd = 1, ProductSpecNbr
      END IF ! PartBound%Dielectric(locBCID)

#if USE_HDG
      ! Method 2: Virtual dielectric layer (VDL)
      IF(ABS(PartBound%PermittivityVDL(locBCID)).GT.0.0)THEN
        ! Create new particles
        DO iNewPart = 1, ProductSpecNbr
          ! Set velocity to zero
          NewVelo(1:3) = 0.0
          ! Create new position by using POI
          NewPos(1:3) = PartPosImpact(1:3)
          ! Create new particle: in case of vMPF or VarTimeStep, new particle inherits the values of the old particle
          ! Routine sets PartState and LastPartPos to the given position (in this case the POI)
          CALL CreateParticle(ProductSpec(2),NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0.,OldPartID=PartID,NewPartID=NewPartID)
          ! This routines moves the particle away from the boundary using LastPartPos as the starting point and changes the species index
          ! and before the particle is deposited and removed, the species index cannot be used anymore
          CALL VirtualDielectricLayerDisplacement(NewPartID,SideID,n_Loc)
          ! Invert species index to invert the charge later in the deposition step
          PartSpecies(NewPartID) = -PartSpecies(NewPartID)
        END DO ! iNewPart = 1, ProductSpecNbr
      END IF ! ABS(PartBound%PermittivityVDL(locBCID)).GT.0.0
#endif /*USE_HDG*/
    END IF ! DoDeposition.AND.DoDielectricSurfaceCharge
  END IF ! ProductSpec(2).GT.0

  ! Decide the fate of the impacting particle
  IF (ProductSpec(1).LE.0) THEN
    ! This routine also calls UpdateBPO()
    CALL RemoveParticle(PartID,BCID=PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID)))
  ELSE
    CALL MaxwellScattering(PartID,SideID,n_Loc)
  END IF
#if USE_HDG
!-----------------------------------------------------------------------------------------------------------------------------------
CASE (VDL_MODEL_ID)  ! Virtual dielectric layer (VDL)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Check if BPO boundary is encountered: Do this here because the BC info is lost during MPI communication of these particles
  ! Normally, this happens in RemoveParticle()
  IF(CalcBoundaryParticleOutput) CALL UpdateBPO(PartID,SideID)
  ! Set the LastPartPos to the POI
  LastPartPos(1:3,PartID) = PartPosImpact(1:3)
  ! This routines moves the particle away from the boundary using LastPartPos as the starting point and changes the species index
  ! and before the particle is deposited and removed, the species index cannot be used anymore
  CALL VirtualDielectricLayerDisplacement(PartID,SideID,n_Loc)
#endif /*USE_HDG*/
CASE DEFAULT
  CALL abort(__STAMP__,'Unknown surface model. PartBound%SurfaceModel(locBCID) = ',IntInfoOpt=PartBound%SurfaceModel(locBCID))
END SELECT

!===================================================================================================================================
! 4.) PIC ONLY: Deposit charges on dielectric surface (when activated), if these were removed/changed in SurfaceModel
!===================================================================================================================================
IF(DoDeposition.AND.DoDielectricSurfaceCharge) THEN ! Surface charging active

  ! Method 1: PartBound%Dielectric = T (dielectric boundary)
  IF(PartBound%Dielectric(locBCID))THEN
    ! Check what happened to the impacting particle
    IF(.NOT.PDM%ParticleInside(PartID))THEN
      ! Particle was deleted on surface contact: deposit impacting charge
      CALL DepositParticleOnNodes(ChargeImpact, PartPosImpact, GlobalElemID)
    ELSEIF(PDM%ParticleInside(PartID))THEN
      ! Sanity check
      IF(PartSpecies(PartID).LT.0)THEN
        IPWRITE (*,*) "PartID        :", PartID
        IPWRITE (*,*) "global ElemID :", GlobalElemID
        CALL abort(__STAMP__,'SurfaceModel() -> DepositParticleOnNodes(): Negative PartSpecies')
      END IF
      ! Particle may have been swapped: check difference in charge
      IF(usevMPF)THEN
        MPF = PartMPF(PartID)
      ELSE
        MPF = Species(PartSpecies(PartID))%MacroParticleFactor
      END IF ! usevMPF
      ChargeRefl = Species(PartSpecies(PartID))%ChargeIC*MPF
      ! Calculate the charge difference between the impacting and reflecting particle
      CALL DepositParticleOnNodes(ChargeImpact-ChargeRefl, PartPosImpact, GlobalElemID)
    END IF ! .NOT.PDM%ParticleInside(PartID)
  END IF ! PartBound%Dielectric(locBCID)

#if USE_HDG
  ! Method 2: Virtual dielectric layer (VDL)
  IF(ABS(PartBound%PermittivityVDL(locBCID)).GT.0.0)THEN
    ! Check what happened to the impacting particle
    IF(.NOT.PDM%ParticleInside(PartID))THEN
      ! Particle was deleted on surface contact: deposit impacting charge but shifted by the VDL offset away from the surface
      ! Set velocity to zero
      NewVelo(1:3) = 0.0
      ! Create new position by using POI
      NewPos(1:3) = PartPosImpact(1:3)
      ! Create new particle: in case of vMPF or VarTimeStep, new particle inherits the values of the old particle
      ! Routine sets PartState and LastPartPos to the given position (in this case the POI)
      CALL CreateParticle(ABS(ProductSpec(1)),NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0.,OldPartID=PartID,NewPartID=NewPartID)
      ! This routines moves the particle away from the boundary using LastPartPos as the starting point and changes the species index
      ! and before the particle is deposited and removed, the species index cannot be used anymore
      CALL VirtualDielectricLayerDisplacement(NewPartID,SideID,n_Loc)
    ELSEIF(PDM%ParticleInside(PartID))THEN
      ! Particle is still inside and was not deleted
      SELECT CASE(PartBound%SurfaceModel(locBCID))
      CASE (VDL_MODEL_ID)  ! Virtual dielectric layer (VDL)
        ! Particle is still inside because it possibly needs to be communicated via MPI to a different process where it is killed
      CASE DEFAULT
        CALL abort(__STAMP__,'Reflected particles not implemented for VDL in combination with, e.g., SEE')
      END SELECT
    END IF ! .NOT.PDM%ParticleInside(PartID)
  END IF ! ABS(PartBound%PermittivityVDL(locBCID)).GT.0.0
#endif /*USE_HDG*/
END IF ! DoDeposition.AND.DoDielectricSurfaceCharge

!===================================================================================================================================
! 5.) Count and sample the properties AFTER the surface interaction
!===================================================================================================================================
! Counter for surface analyze
IF(CalcSurfCollCounter) THEN
  ! Old particle was deleted/absorbed/adsorbed at the wall (in case it happend during a surface model and not during the
  ! SpeciesSwap routine)
  IF (ProductSpec(1).LE.0) THEN
    SurfAnalyzeNumOfAds(PartSpecImpact) = SurfAnalyzeNumOfAds(PartSpecImpact) + 1
  END IF
  ! New particle was created/desorbed at the wall
  IF (ProductSpec(2).GT.0) THEN
    SurfAnalyzeNumOfDes(ProductSpec(2)) = SurfAnalyzeNumOfDes(ProductSpec(2)) + ProductSpecNbr
  END IF
END IF
! Sampling
IF(DoSample) THEN
  ! Sample momentum, heatflux and collision counter on surface (only the original impacting particle, not the newly created parts
  ! through SurfaceModelParticleEmission), checking if the particle was deleted/absorbed/adsorbed at the wall
  IF(PDM%ParticleInside(PartID)) CALL CalcWallSample(PartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
END IF

END SUBROUTINE SurfaceModelling

SUBROUTINE SpeciesSwap(PartID,SideID,targetSpecies_IN)
!----------------------------------------------------------------------------------------------------------------------------------!
! Computes the Species Swap on ReflectiveBC
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals                 ,ONLY: abort,VECNORM
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Vars           ,ONLY: PartSpecies
USE MOD_part_operations         ,ONLY: RemoveParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: PartID, SideID
INTEGER,INTENT(IN),OPTIONAL       :: targetSpecies_IN
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: targetSpecies, iSwaps, locBCID
REAL                              :: RanNum
LOGICAL                           :: PerformSwap
!===================================================================================================================================

PerformSwap = .TRUE.
locBCID = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
IF(PartBound%ProbOfSpeciesSwaps(locBCID).LT.1.) THEN
  CALL RANDOM_NUMBER(RanNum)
  PerformSwap = RanNum.LE.PartBound%ProbOfSpeciesSwaps(locBCID)
END IF

IF(PerformSwap) THEN
  targetSpecies=-1 ! Dummy initialization value
  IF(PRESENT(targetSpecies_IN))THEN
    targetSpecies = targetSpecies_IN
  ELSE ! Normal swap routine
    DO iSwaps=1,PartBound%NbrOfSpeciesSwaps(locBCID)
      IF (PartSpecies(PartID).EQ.PartBound%SpeciesSwaps(1,iSwaps,locBCID)) targetSpecies = PartBound%SpeciesSwaps(2,iSwaps,locBCID)
    END DO
  END IF ! PRESENT(targetSpecies_IN)
  ! Swap species
  IF (targetSpecies.EQ.0) THEN
    ! Delete particle -> same as PartBound%OpenBC
    CALL RemoveParticle(PartID,BCID=locBCID)
  ELSE IF (targetSpecies.GT.0) THEN
    ! Swap species
    PartSpecies(PartID)=targetSpecies
  END IF ! targetSpecies.EQ.0
END IF ! RanNum.LE.PartBound%ProbOfSpeciesSwaps

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

#if USE_HDG
!===================================================================================================================================
!> Shift impacting particles by a specific displacement away from the boundary in the normal direction. Change the species ID to
!> flag such particles so that later, after MPI communication, they can be deleted and deposited at the target position to form a
!> surface charge on a (virtual) dielectric layer.
!> Note: LastPartPos must be set to the impact position on the surface
!===================================================================================================================================
SUBROUTINE VirtualDielectricLayerDisplacement(PartID,SideID,n_Loc)
! MODULES
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Particle_Vars          ,ONLY: SpeciesOffsetVDL
USE MOD_Particle_Vars          ,ONLY: LastPartPos, PartSpecies, PartState
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_Particle_Tracking_Vars ,ONLY: TrackInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)
INTEGER,INTENT(IN) :: PartID, SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPartBound
!===================================================================================================================================
! Get particle boundary index
iPartBound = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))

! Set new particle position: Shift away from face by ratio of real VDL thickness and permittivity.
! Note the negative sign that is required as the normal vector points outwards
PartState(1:3,PartID) = LastPartPos(1:3,PartID) - (PartBound%ThicknessVDL(iPartBound)/PartBound%PermittivityVDL(iPartBound))*n_loc

! Set tracking variables
TrackInfo%PartTrajectory=PartState(1:3,PartID) - LastPartPos(1:3,PartID)

TrackInfo%lengthPartTrajectory = VECNORM(TrackInfo%PartTrajectory)
IF(ALMOSTZERO(TrackInfo%lengthPartTrajectory)) THEN
  TrackInfo%lengthPartTrajectory= 0.0
ELSE
  TrackInfo%PartTrajectory=TrackInfo%PartTrajectory/TrackInfo%lengthPartTrajectory
END IF

! Encode species index: Set a temporarily invalid number, which holds the information that the particle has interacted with a VDL.
! The Particle is removed after MPI communication because the new position might be on a different process due to the displacement
! Check if the particle has an encoded species index
IF(PartSpecies(PartID).GT.SpeciesOffsetVDL)THEN ! Check if already encoded
  ! Do not encode twice (re-bouncing lost particles)
ELSE
  PartSpecies(PartID) = PartSpecies(PartID) + SpeciesOffsetVDL
END IF ! PartSpecies(iPart).GE.SpeciesOffsetVDL

END SUBROUTINE VirtualDielectricLayerDisplacement
#endif /*USE_HDG*/


SUBROUTINE StickingCoefficientModel(PartID,SideID,n_Loc)
!===================================================================================================================================
!> Empirical sticking coefficient model using the product of a non-bounce probability (angle dependence with a cut-off angle) and
!> condensation probability (linear temperature dependence, using different temperature limits). Particle sticking to the wall
!> will be simply deleted, transferring the complete energy to the wall heat flux
!> Assumed data input is [upper impact angle limit (deg), alpha_B (deg), T_1 (K), T_2 (K)]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, SampWallState, GlobalSide2SurfSide, SWIStickingCoefficient
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_part_operations         ,ONLY: RemoveParticle
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
USE MOD_SurfaceModel_Vars       ,ONLY: StickingCoefficientData
USE MOD_DSMC_Vars               ,ONLY: DSMC, SamplingActive
USE MOD_Particle_Vars           ,ONLY: WriteMacroSurfaceValues
USE MOD_SurfaceModel_Tools      ,ONLY: MaxwellScattering, GetWallTemperature
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: n_loc(1:3)
INTEGER,INTENT(IN)              :: PartID, SideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iVal, SubP, SubQ, locBCID, SurfSideID
REAL                            :: ImpactAngle, NonBounceProb, Prob, alpha_B, LowerTemp, UpperTemp, RanNum, WallTemp
LOGICAL                         :: CalcStickCoeff, ParticleSticks
!===================================================================================================================================

CalcStickCoeff = .TRUE.
ParticleSticks = .FALSE.
Prob = 0.
locBCID = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
! Get the wall temperature
WallTemp = GetWallTemperature(PartID,locBCID,SideID)
! Determine impact angle of particle
ImpactAngle = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(TrackInfo%PartTrajectory,n_loc))))
! Select the appropriate temperature coefficients and bounce angle
IF(UBOUND(StickingCoefficientData,dim=2).GE.2) THEN
  ! More than 1 row of coefficients
  DO iVal = 1, UBOUND(StickingCoefficientData,dim=2)
    IF(ImpactAngle.LE.StickingCoefficientData(1,iVal)) THEN
      alpha_B   = StickingCoefficientData(2,iVal)
      LowerTemp = StickingCoefficientData(3,iVal)
      UpperTemp = StickingCoefficientData(4,iVal)
      CalcStickCoeff = .TRUE.
      EXIT
    ELSE
      ! Impact angle is outside coefficient data range
      CalcStickCoeff = .FALSE.
    END IF
  END DO
ELSE
  ! 1 row of coefficients
  IF(ImpactAngle.LE.StickingCoefficientData(1,1)) THEN
    alpha_B   = StickingCoefficientData(2,1)
    LowerTemp = StickingCoefficientData(3,1)
    UpperTemp = StickingCoefficientData(4,1)
  ELSE
    ! Impact angle is outside coefficient data range
    CalcStickCoeff = .FALSE.
  END IF
END IF

! Calculate the sticking coefficient and compare with random number
IF(CalcStickCoeff) THEN
  ! Calculate the non-bounce probability
  IF(ImpactAngle.GE.alpha_B) THEN
    NonBounceProb = (90.-ImpactAngle) / (90.-alpha_B)
  ELSE
    NonBounceProb = 1.
  END IF
  ! Calculate the sticking probability, using the wall temperature
  IF(WallTemp.LT.LowerTemp) THEN
    Prob = NonBounceProb
  ELSE IF(WallTemp.GT.UpperTemp) THEN
    Prob = 0.
  ELSE
    Prob = NonBounceProb * (UpperTemp - WallTemp) / (UpperTemp - LowerTemp)
  END IF
  ! Determine whether the particles sticks to the boundary
  IF(Prob.GT.0.) THEN
    CALL RANDOM_NUMBER(RanNum)
    IF(Prob.GT.RanNum) ParticleSticks = .TRUE.
  END IF
END IF

IF(ParticleSticks) THEN
  ! Remove the particle from the simulation (total energy was added to the sampled heat flux before the interaction)
  CALL RemoveParticle(PartID,BCID=locBCID)
ELSE
  ! Perform regular Maxwell scattering
  CALL MaxwellScattering(PartID,SideID,n_Loc)
END IF

! Sampling the sticking coefficient
IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
  SubP = TrackInfo%p
  SubQ = TrackInfo%q
  SampWallState(SWIStickingCoefficient,SubP,SubQ,SurfSideID) = SampWallState(SWIStickingCoefficient,SubP,SubQ,SurfSideID) + Prob
END IF

END SUBROUTINE StickingCoefficientModel

END MODULE MOD_SurfaceModel
