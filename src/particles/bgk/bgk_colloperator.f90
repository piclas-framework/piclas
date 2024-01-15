!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_BGK_CollOperator
!===================================================================================================================================
! Module approximating the collision operator using the Bhatnagar-Gross-Krook model
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE BGK_CollisionOperator
  MODULE PROCEDURE BGK_CollisionOperator
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: BGK_CollisionOperator, ARShakhov
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_CollisionOperator(iPartIndx_Node, nPart, NodeVolume, AveragingValues)
!===================================================================================================================================
!> Subroutine for the cell-local BGK collision operator:
!> 1.) Moment calculation: Summing up the relative velocities and their squares
!> 2.) Calculation of the relaxation frequency of the distribution function towards the target distribution function
!> 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency
!> 4.) Determine the number of particles undergoing a relaxation (including vibration and rotation)
!> 5.) Determine the new rotational and vibrational state of molecules undergoing a relaxation
!> 6.) Sample new particle velocities from the target distribution function, depending on the chosen model
!> 7.) Determine the new bulk velocity and the new relative velocity of the particles
!> 8.) Treatment of the vibrational energy of molecules
!> 9.) Determine the new PartState (for molecules, including rotational energy)
!> 9.) Scaling of the rotational energy of molecules
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: DOTPRODUCT, CROSS
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies, usevMPF, UseVarTimeStep
USE MOD_Particle_Vars         ,ONLY: UseRotRefFrame, RotRefFrameOmega, PartVeloRotRef, PDM
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, DSMC, PartStateIntEn, PolyatomMolDSMC, RadialWeighting, CollInf
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_BGK_Vars              ,ONLY: SpecBGK, BGKDoVibRelaxation, BGKMovingAverage
USE MOD_BGK_Vars              ,ONLY: BGK_MeanRelaxFactor, BGK_MeanRelaxFactorCounter, BGK_MaxRelaxFactor, BGK_MaxRotRelaxFactor
USE MOD_BGK_Vars              ,ONLY: BGK_PrandtlNumber, BGK_Viscosity, BGK_ThermalConductivity
USE MOD_part_tools            ,ONLY: GetParticleWeight
#ifdef CODE_ANALYZE
USE MOD_Globals               ,ONLY: abort,unit_stdout,myrank
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                        :: NodeVolume
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(INOUT), OPTIONAL           :: AveragingValues(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: vBulk(3), u0ij(3,3), u2, V_rel(3), dtCell
REAL                  :: alpha, alphaRot(nSpecies), CellTemp, dens, InnerDOF, NewEn, OldEn, Prandtl, relaxfreq
REAL                  :: dynamicvis, thermalcond
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelax(:),iPartIndx_NodeRelaxTemp(:),iPartIndx_NodeRelaxRot(:),iPartIndx_NodeRelaxVib(:)
INTEGER               :: iLoop, iPart, nRelax, nXiVibDOF
REAL, ALLOCATABLE     :: Xi_vib_DOF(:,:), VibEnergyDOF(:,:)
INTEGER               :: iSpec, nSpec(nSpecies), jSpec, nRotRelax, nVibRelax
REAL                  :: OldEnRot, NewEnRot(nSpecies), NewEnVib(nSpecies)
REAL                  :: TotalMass, u2Spec(nSpecies), u2i(3), vBulkAll(3)
REAL                  :: SpecTemp(nSpecies)
#ifdef CODE_ANALYZE
REAL                  :: Energy_old,Energy_new,Momentum_old(3),Momentum_new(3)
INTEGER               :: iMom
REAL,PARAMETER        :: RelMomTol=1e-6  ! Relative tolerance applied to conservation of momentum before/after reaction
REAL,PARAMETER        :: RelEneTol=1e-12 ! Relative tolerance applied to conservation of energy before/after reaction
#endif /* CODE_ANALYZE */
REAL                  :: totalWeightSpec(nSpecies), totalWeight, totalWeight2, partWeight, CellTemptmp, MassIC_Mixture
REAL                  :: EVibSpec(nSpecies), Xi_VibSpec(nSpecies), Xi_VibSpecNew(nSpecies)
REAL                  :: ERotSpec(nSpecies), Xi_RotSpec(nSpecies), Xi_RotTotal
REAL                  :: CellTempRel, TEqui
REAL                  :: TVibSpec(nSpecies), TRotSpec(nSpecies), VibRelaxWeightSpec(nSpecies), RotRelaxWeightSpec(nSpecies)
REAL                  :: collisionfreqSpec(nSpecies), rotrelaxfreqSpec(nSpecies), vibrelaxfreqSpec(nSpecies)
REAL                  :: betaR(nSpecies), betaV(nSpecies)
!===================================================================================================================================
#ifdef CODE_ANALYZE
! Momentum and energy conservation check: summing up old values
Momentum_new = 0.0; Momentum_old = 0.0; Energy_new = 0.0; Energy_old = 0.0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  Momentum_old(1:3) = Momentum_old(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  Energy_old = Energy_old + DOTPRODUCT(PartState(4:6,iPart))*0.5*Species(iSpec)%MassIC*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    ! Add internal energies (vibration, rotation) for molecules and molecular ions
    Energy_old = Energy_old + (PartStateIntEn(1,iPart) + PartStateIntEn(2,iPart))*partWeight
  END IF
END DO
#endif

IF(nPart.LE.2) RETURN

! 1.) Moment calculation: Summing up the relative velocities and their squares
CALL CalcMoments(nPart, iPartIndx_Node, nSpec, vBulkAll, totalWeight, totalWeight2, totalWeightSpec, TotalMass,  u2, u2Spec, u0ij, &
                 u2i, OldEn, EVibSpec, ERotSpec, CellTemp, SpecTemp, dtCell)

IF((CellTemp.LE.0.0).OR.(MAXVAL(nSpec(:)).EQ.1).OR.(totalWeight.LE.0.0)) RETURN

IF(UseVarTimeStep) THEN
  dtCell = dt * dtCell / totalWeight
ELSE
  dtCell = dt
END IF

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  dens = totalWeight / NodeVolume
ELSE
  ! MPF is the same for all species
  dens = totalWeight * Species(1)%MacroParticleFactor / NodeVolume
END IF

IF (BGKMovingAverage) THEN
  CALL DoAveraging(dens, u2, u0ij, u2i, CellTemp, AveragingValues)
END IF

CALL CalcInnerDOFs(nSpec, EVibSpec, ERotSpec, totalWeightSpec, TVibSpec, TRotSpec, InnerDOF, Xi_VibSpec, Xi_RotSpec)

! 2.) Calculation of the relaxation frequency of the distribution function towards the target distribution function
CALL CalcGasProperties(nSpec, dens, InnerDOF, totalWeightSpec, totalWeight, TotalMass, u2Spec, SpecTemp, CellTemp, Xi_VibSpec, &
    Xi_RotSpec, Prandtl, relaxfreq, dynamicvis, thermalcond, MassIC_Mixture)

IF(DSMC%CalcQualityFactors) THEN
  BGK_MeanRelaxFactor         = BGK_MeanRelaxFactor + relaxfreq * dtCell
  BGK_MeanRelaxFactorCounter  = BGK_MeanRelaxFactorCounter + 1
  BGK_MaxRelaxFactor          = MAX(BGK_MaxRelaxFactor,relaxfreq*dtCell)
  BGK_PrandtlNumber           = BGK_PrandtlNumber + Prandtl
  BGK_Viscosity               = BGK_Viscosity + dynamicvis
  BGK_ThermalConductivity     = BGK_ThermalConductivity + thermalcond
END IF

! 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency using the collision frequency,
!     which is not the same as the relaxation frequency of distribution function, calculated above
IF(ANY(SpecDSMC(:)%InterID.EQ.2).OR.ANY(SpecDSMC(:)%InterID.EQ.20)) THEN
  collisionfreqSpec = 0.0
  DO iSpec = 1, nSpecies
    DO jSpec = 1, nSpecies
      IF (iSpec.EQ.jSpec) THEN
        CellTemptmp = CellTemp !SpecTemp(iSpec)
      ELSE
        CellTemptmp = CellTemp
      END IF
      ! Sum up collision frequencies of species i with itself and the other species
      ! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 87f
      ! For SpecBGK(iSpec)%CollFreqPreFactor(jSpec) see bgk_init.f90
      ! VHS according to M. Pfeiffer, "Extending the particle ellipsoidal statistical Bhatnagar-Gross-Krook method to diatomic
      ! molecules including quantized vibrational energies", Phys. Fluids 30, 116103 (2018), Eq. (18)
      collisionfreqSpec(iSpec) = collisionfreqSpec(iSpec) + SpecBGK(iSpec)%CollFreqPreFactor(jSpec) * totalWeightSpec(jSpec) &
              * (Dens / totalWeight) *CellTemptmp**(-CollInf%omega(iSpec,jSpec) +0.5)
    END DO
  END DO
  ! Calculate relaxation frequencies of rotation and vibration with relaxation properties
  ! M. Pfeiffer, "Extending the particle ellipsoidal statistical Bhatnagar-Gross-Krook method to diatomic molecules including
  ! quantized vibrational energies", Phys. Fluids 30, 116103 (2018)
  ! N.E. Gimelshein et. al, "Vibrational relaxation rates in the direct simulation Monte Carlo method", Phys. Fluids 14, 4452 (2018)
  ! relaxfreqSpec = collisionfreqSpec / collision number Z with RelaxProb = 1/Z
  rotrelaxfreqSpec(:) = collisionfreqSpec(:) * DSMC%RotRelaxProb
  vibrelaxfreqSpec(:) = collisionfreqSpec(:) * DSMC%VibRelaxProb
  
  IF(DSMC%CalcQualityFactors) THEN
    BGK_MaxRotRelaxFactor          = MAX(BGK_MaxRotRelaxFactor,MAXVAL(rotrelaxfreqSpec(:))*dtCell)
  END IF
END IF

! 4.) Determine the relaxation temperatures and the number of particles undergoing a relaxation (including vibration + rotation)
ALLOCATE(iPartIndx_NodeRelax(nPart), iPartIndx_NodeRelaxTemp(nPart))
iPartIndx_NodeRelax = 0; iPartIndx_NodeRelaxTemp = 0
ALLOCATE(iPartIndx_NodeRelaxRot(nPart),iPartIndx_NodeRelaxVib(nPart))
iPartIndx_NodeRelaxRot = 0; iPartIndx_NodeRelaxVib = 0

! Allocate Xi_vib_DOF
IF(BGKDoVibRelaxation) THEN
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    nXiVibDOF = MAXVAL(PolyatomMolDSMC(:)%VibDOF)
    ALLOCATE(Xi_vib_DOF(DSMC%NumPolyatomMolecs,nXiVibDOF))
  END IF
END IF

CALL CalcTRelax(ERotSpec, Xi_RotSpec, Xi_VibSpec, EVibSpec, totalWeightSpec, totalWeight, totalWeight2, dtCell, CellTemp, &
    TRotSpec, TVibSpec, relaxfreq, rotrelaxfreqSpec, vibrelaxfreqSpec, Xi_VibSpecNew, Xi_vib_DOF, nXiVibDOF, CellTempRel, TEqui, &
    betaR, betaV)

CALL DetermineRelaxPart(nPart, iPartIndx_Node, relaxfreq, dtCell, nRelax, nRotRelax, nVibRelax, &
    RotRelaxWeightSpec, VibRelaxWeightSpec, iPartIndx_NodeRelax, iPartIndx_NodeRelaxTemp, iPartIndx_NodeRelaxRot, &
    iPartIndx_NodeRelaxVib, vBulk, OldEnRot, OldEn, rotrelaxfreqSpec, vibrelaxfreqSpec, betaR, betaV)

! Allocate VibEnergyDOF
IF(BGKDoVibRelaxation) THEN
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    ALLOCATE(VibEnergyDOF(nVibRelax,nXiVibDOF))
    VibEnergyDOF = 0.0
  END IF
END IF

! Return if no particles are undergoing a relaxation
IF ((nRelax.EQ.0).AND.(nRotRelax.EQ.0).AND.(nVibRelax.EQ.0)) RETURN

! 5.) Determine the new rotational and vibrational states of molecules undergoing a relaxation
IF(ANY(SpecDSMC(:)%InterID.EQ.2).OR.ANY(SpecDSMC(:)%InterID.EQ.20)) THEN
  CALL RelaxInnerEnergy(nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib, iPartIndx_NodeRelaxRot, nXiVibDOF, Xi_vib_DOF, &
    Xi_VibSpecNew, Xi_RotSpec, VibEnergyDOF, TEqui, NewEnVib, NewEnRot)
END IF

! 6.) Sample new particle velocities from the target distribution function, depending on the chosen model
CALL SampleFromTargetDistr(nRelax, iPartIndx_NodeRelax, Prandtl, u2, u0ij, u2i, vBulkAll, CellTempRel, CellTemp, vBulk, &
  MassIC_Mixture)

NewEn = 0.

! Calculation of the new bulk velocity
vBulk = vBulk/TotalMass

! Loop over all relaxing particles for calculation of the new energy
DO iLoop = 1, nRelax
  iPart = iPartIndx_NodeRelax(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  ! Thermal velocity of all relaxing particles
  V_rel(1:3) = PartState(4:6,iPart) - vBulk(1:3)
  ! Sum up kinetic energies
  NewEn = NewEn + (V_rel(1)**(2.) + V_rel(2)**(2.) + V_rel(3)**(2.))*0.5*Species(iSpec)%MassIC*partWeight
END DO
! Loop over all non-relaxing particles for calculation of the new energy
DO iLoop = 1, nPart-nRelax
  iPart = iPartIndx_NodeRelaxTemp(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  ! Thermal velocity of all non-relaxing particles
  V_rel(1:3) = PartState(4:6,iPart) - vBulk(1:3)
  ! Sum up kinetic energies
  NewEn = NewEn + (V_rel(1)**(2.) + V_rel(2)**(2.) + V_rel(3)**(2.))*0.5*Species(iSpec)%MassIC*partWeight
END DO

! 7.) Vibrational energy of the molecules: Ensure energy conservation by scaling the new vibrational states with the factor alpha
IF(ANY(SpecDSMC(:)%InterID.EQ.2).OR.ANY(SpecDSMC(:)%InterID.EQ.20)) THEN
  CALL EnergyConsVib(nPart, totalWeight, nVibRelax, VibRelaxWeightSpec, iPartIndx_NodeRelaxVib, NewEnVib, OldEn, nXiVibDOF, &
    VibEnergyDOF, Xi_VibSpecNew, TEqui)
END IF

! Remaining vibrational + translational energy + old rotational energy for translation and rotation
OldEn = OldEn + OldEnRot
Xi_RotTotal = 0.0
DO iSpec = 1, nSpecies
  ! Sum of relaxing rotational degrees of freedom
  Xi_RotTotal = Xi_RotTotal + Xi_RotSpec(iSpec)*RotRelaxWeightSpec(iSpec)
END DO

! 8.) Determine the new particle state and ensure energy conservation by scaling the new velocities with the factor alpha
! Calculation of scaling factor alpha
alpha = SQRT(OldEn/NewEn*(3.*(totalWeight-1.))/(Xi_RotTotal+3.*(totalWeight-1.)))
! Calculation of the final particle velocities with vBulkAll (average flow velocity before relaxation), scaling factor alpha,
! the particle velocity PartState(4:6,iPart) after the relaxation but before the energy conservation and vBulk (average value of
! the latter)
DO iLoop = 1, nRelax
  iPart = iPartIndx_NodeRelax(iLoop)
  PartState(4:6,iPart) = vBulkAll(1:3) + alpha*(PartState(4:6,iPart)-vBulk(1:3))
END DO
DO iLoop = 1, nPart-nRelax
  iPart = iPartIndx_NodeRelaxTemp(iLoop)
  PartState(4:6,iPart) = vBulkAll(1:3) + alpha*(PartState(4:6,iPart)-vBulk(1:3))
END DO

IF(UseRotRefFrame) THEN
  ! Resetting the velocity in the rotational frame of reference for particles that underwent a relaxation
  DO iLoop = 1, nRelax
    iPart = iPartIndx_NodeRelax(iLoop)
    IF(PDM%InRotRefFrame(iPart)) THEN
      PartVeloRotRef(1:3,iPart) = PartState(4:6,iPart) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,iPart))
    END IF
  END DO
END IF

! 9.) Rotation: Scale the new rotational state of the molecules to ensure energy conservation
DO iSpec = 1, nSpecies
  ! Calculate scaling factor alpha per species, see F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal
  ! statistical Bhatnagar-Gross-Krook method including internal degrees of freedom", subitted to Phys. Fluids, August 2023
  IF (NewEnRot(iSpec).GT.0.0) THEN
    alphaRot(iSpec) = OldEn/NewEnRot(iSpec)*(Xi_RotSpec(iSpec)*RotRelaxWeightSpec(iSpec)/(Xi_RotTotal+3.*(totalWeight-1.)))
  ELSE
    alphaRot(iSpec) = 0.0
  END IF
END DO
DO iLoop = 1, nRotRelax
  iPart = iPartIndx_NodeRelaxRot(iLoop)
  iSpec = PartSpecies(iPart)
  ! Scaling of rotational energy with factor alpha
  PartStateIntEn( 2,iPart) = alphaRot(iSpec)*PartStateIntEn( 2,iPart)
END DO

! CODE ANALYZE: Compare the old momentum and energy of the cell with the new, abort if relative difference is above the limits
#ifdef CODE_ANALYZE
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  Momentum_new(1:3) = Momentum_new(1:3) + (PartState(4:6,iPart)) * Species(iSpec)%MassIC*partWeight
  Energy_new = Energy_new + DOTPRODUCT((PartState(4:6,iPart)))*0.5*Species(iSpec)%MassIC*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    Energy_new = Energy_new + (PartStateIntEn(1,iPart) + PartStateIntEn(2,iPart))*partWeight
  END IF
END DO
! Check for energy difference
IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,RelEneTol)) THEN
  WRITE(UNIT_StdOut,*) '\n'
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_old-Energy_new
  ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
    IF(energy.GT.0.0)THEN
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_old-Energy_new)/energy
    END IF
  END ASSOCIATE
  IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelEneTol
  IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha           : ", OldEn, alpha
  IPWRITE(UNIT_StdOut,*)                     " nPart, nRelax, nRotRelax, nVibRelax: ", nPart, nRelax
  CALL abort(&
      __STAMP__&
      ,'CODE_ANALYZE: BGK_CollisionOperator is not energy conserving!')
END IF
! Check for momentum difference
DO iMom=1,3
  IF (.NOT.ALMOSTEQUALRELATIVE(Momentum_old(iMom),Momentum_new(iMom),RelMomTol)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Direction (x,y,z)        : ",iMom
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_old             : ",Momentum_old(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_new             : ",Momentum_new(iMom)
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Momentum difference : ",Momentum_old(iMom)-Momentum_new(iMom)
    ASSOCIATE( Momentum => MAX(ABS(Momentum_old(iMom)),ABS(Momentum_new(iMom))) )
      IF(Momentum.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Momentum difference : ",(Momentum_old(iMom)-Momentum_new(iMom))/Momentum
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance        : ",RelMomTol
    IPWRITE(UNIT_StdOut,*)                     " OldEn, alpha             : ", OldEn, alpha
    IPWRITE(UNIT_StdOut,*)                     " nPart, nRelax, nRotRelax, nVibRelax: ", nPart, nRelax
    CALL abort(&
        __STAMP__&
        ,'CODE_ANALYZE: BGK_CollisionOperator is not momentum conserving!')
  END IF
END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE BGK_CollisionOperator


SUBROUTINE CalcMoments(nPart, iPartIndx_Node, nSpec, vBulkAll, totalWeight, totalWeight2, totalWeightSpec, TotalMass, u2, u2Spec, &
    u0ij, u2i, OldEn, EVibSpec, ERotSpec, CellTemp, SpecTemp, dtCell)
!===================================================================================================================================
!> Moment calculation: Summing up the relative velocities and their squares
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies, UseVarTimeStep, PartTimeStep
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation, BGKCollModel
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Globals               ,ONLY: DOTPRODUCT
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, iPartIndx_Node(nPart)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)          :: nSpec(nSpecies)
REAL, INTENT(OUT)             :: u2Spec(nSpecies),u0ij(3,3), OldEn, EVibSpec(nSpecies), ERotSpec(nSpecies), u2i(3), u2
REAL, INTENT(OUT)             :: CellTemp, SpecTemp(nSpecies), totalWeightSpec(nSpecies)
REAL, INTENT(OUT)             :: vBulkAll(3), totalWeight, totalWeight2, TotalMass, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iPart, iSpec, fillMa1, fillMa2
REAL                          :: V_rel(1:3), vmag2, EnerTotal, ThermEner, totalWeightSpec2(nSpecies), vBulkSpec(3,nSpecies)
REAL                          :: partWeight, tempweight, tempweight2, tempmass, vBulkTemp(3), totalWeight3
LOGICAL                       :: validSpec(nSpecies)
!===================================================================================================================================
totalWeightSpec = 0.0; totalWeightSpec2=0.0; vBulkAll=0.0; TotalMass=0.0; vBulkSpec=0.0; nSpec=0; dtCell=0.0
! Loop over all simulation particles to sum up particle velocities to calculate bulk velocities
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  totalWeightSpec(iSpec) = totalWeightSpec(iSpec) + partWeight
  totalWeightSpec2(iSpec) =   totalWeightSpec2(iSpec) + partWeight*partWeight
  vBulkAll(1:3)  =  vBulkAll(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  TotalMass = TotalMass + Species(iSpec)%MassIC*partWeight
  vBulkSpec(1:3,iSpec) = vBulkSpec(1:3,iSpec) + PartState(4:6,iPart)*partWeight
  nSpec(iSpec) = nSpec(iSpec) + 1 ! Count number of simulation particles per species
  IF(UseVarTimeStep) THEN
    dtCell = dtCell + PartTimeStep(iPart)*partWeight
  END IF
END DO
totalWeight = SUM(totalWeightSpec)
totalWeight2 = SUM(totalWeightSpec2)

! Calculate total bulk velocity and bulk velocities per species
vBulkAll(1:3) = vBulkAll(1:3) / TotalMass
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).GT.0) vBulkSpec(:,iSpec) = vBulkSpec(:,iSpec) / totalWeightSpec(iSpec)
END DO

totalWeight3 = 0.; u2Spec=0.0; u0ij=0.0; u2i=0.0; OldEn=0.0; EVibSpec=0.0; ERotSpec=0.0
! Loop over all simulation particles to sum up relative velocities
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  iSpec = PartSpecies(iPart)
  ! Calculate thermal velocity with bulk velocity of the species
  V_rel(1:3)=PartState(4:6,iPart)-vBulkSpec(1:3,iSpec)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  ! Summing up thermal velocities (squared) of all particles per species
  u2Spec(iSpec) = u2Spec(iSpec) + vmag2*partWeight

  ! Calculate thermal velocity with bulk velocity of the gas
  V_rel(1:3)=PartState(4:6,iPart)-vBulkAll(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  IF (BGKCollModel.EQ.1) THEN ! ESBGK
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ij(fillMa1, fillMa2)= u0ij(fillMa1, fillMa2) &
            + V_rel(fillMa1)*V_rel(fillMa2)*Species(iSpec)%MassIC*partWeight
      END DO
    END DO
  END IF
  IF (BGKCollModel.EQ.2) THEN ! Shakhov
    u2i(1:3) = u2i(1:3) + V_rel(1:3)*vmag2 * partWeight*Species(iSpec)%MassIC
    totalWeight3 = totalWeight3 + partWeight*partWeight*partWeight
  END IF

  ! Sum up old energy of thermal velocities and sum up internal energies --> E_T
  OldEn = OldEn + 0.5*Species(iSpec)%MassIC * vmag2*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(BGKDoVibRelaxation) THEN
      ! EVib without zero-point energy
      EVibSpec(iSpec) = EVibSpec(iSpec) + (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
    END IF
    ERotSpec(iSpec) = ERotSpec(iSpec) + PartStateIntEn(2,iPart) * partWeight
  END IF
END DO

IF (BGKCollModel.EQ.1)  THEN ! ESBGK
  ! Pressure tensor
  u0ij = u0ij* totalWeight / (TotalMass*(totalWeight - totalWeight2/totalWeight))
END IF
IF (BGKCollModel.EQ.2)  THEN ! Shakhov
  ! Heatflux
  u2i = u2i*totalWeight**3/(TotalMass*(totalWeight**3-3.*totalWeight*totalWeight2+2.*totalWeight3))
END IF

! Calculation of cell temperature and bulk velocity (squared)
IF (nSpecies.GT.1) THEN ! mixture
  validSpec = .FALSE.
  SpecTemp = 0.0
  EnerTotal = 0.0
  tempweight = 0.0; tempweight2 = 0.0; tempmass = 0.0; vBulkTemp = 0.0
  DO iSpec = 1, nSpecies
    ! At least two particles and non-zero squared thermal velocity needed for a valid species
    IF ((nSpec(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(u2Spec(iSpec)))) THEN
      validSpec = .TRUE.
      ! Calculation of the species temperature --> translational temperatures of the different species
      SpecTemp(iSpec) = Species(iSpec)%MassIC * u2Spec(iSpec) &
          /(3.0*BoltzmannConst*(totalWeightSpec(iSpec) - totalWeightSpec2(iSpec)/totalWeightSpec(iSpec)))
      ! Thermal energy
      EnerTotal =  EnerTotal + 3./2.*BoltzmannConst*SpecTemp(iSpec) * totalWeightSpec(iSpec)
      vmag2 = DOTPRODUCT(vBulkSpec(1:3,iSpec))
      ! Add kinetic energy
      EnerTotal = EnerTotal + totalWeightSpec(iSpec) * Species(iSpec)%MassIC / 2. * vmag2
      tempweight = tempweight + totalWeightSpec(iSpec)
      tempweight2 = tempweight2 + totalWeightSpec2(iSpec)
      tempmass = tempmass +  totalWeightSpec(iSpec) * Species(iSpec)%MassIC
      vBulkTemp(1:3) = vBulkTemp(1:3) + vBulkSpec(1:3,iSpec)*totalWeightSpec(iSpec) * Species(iSpec)%MassIC
    END IF
  END DO
  IF (ANY(validSpec)) THEN
    vBulkTemp(1:3) = vBulkTemp(1:3) / tempmass
    ! Squared bulk velocity of the mixture
    vmag2 = DOTPRODUCT(vBulkTemp(1:3))
    ! EnerTotal = kinetic energy (tempmass / 2. * vmag2) + thermal energy (3. * tempweight * BoltzmannConst * CellTemp / 2)
    ThermEner = EnerTotal - tempmass / 2. * vmag2
    ! Calculation of the cell temperature from the thermal energy --> translational temperature of the mixture
    CellTemp = 2. * ThermEner / (3.*tempweight*BoltzmannConst)
    ! Mean squared thermal velocity c^2 of a particle, calculated with the cell temperature and the density-averaged mass
    u2 = 3. * CellTemp * BoltzmannConst * (tempweight - tempweight2/tempweight) / tempmass
  ELSE ! only one part per species or cloned species with u2spec = 0 because PartState(4:6) = vBulkAll
    u2 = OldEn / (TotalMass*(1. - totalWeight2/totalWeight**2)) * 2. ! variance-free
    CellTemp = TotalMass/totalWeight * u2 / (3.0*BoltzmannConst)
  END IF
ELSE ! single species gas
  u2 = u2Spec(1) / (totalWeight - totalWeight2/totalWeight)
  CellTemp = Species(1)%MassIC * u2 / (3.0*BoltzmannConst)
END IF

END SUBROUTINE CalcMoments


SUBROUTINE DoAveraging(dens, u2, u0ij, u2i, CellTemp, AverageValues)
!===================================================================================================================================
!> BGK moving average
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: Species
USE MOD_BGK_Vars              ,ONLY: BGKCollModel, BGKMovingAverageFac
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: dens, u0ij(3,3), u2i(3), u2, CellTemp, AverageValues(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!===================================================================================================================================
IF (BGKCollModel.EQ.1) THEN
  IF (AverageValues(1).EQ.0.0) THEN
    AverageValues(1) = dens
    AverageValues(2) = u2
    AverageValues(3) = u0ij(1,1); AverageValues(4) = u0ij(2,2); AverageValues(5) = u0ij(3,3);
    AverageValues(6) = u0ij(1,2); AverageValues(7) = u0ij(1,3); AverageValues(8) = u0ij(2,3);
  ELSE
    AverageValues(1) = BGKMovingAverageFac*dens + (1.-BGKMovingAverageFac)*AverageValues(1)
    AverageValues(2) = BGKMovingAverageFac*u2 + (1.-BGKMovingAverageFac)*AverageValues(2)
    AverageValues(3) = BGKMovingAverageFac*u0ij(1,1) + (1.-BGKMovingAverageFac)*AverageValues(3)
    AverageValues(4) = BGKMovingAverageFac*u0ij(2,2) + (1.-BGKMovingAverageFac)*AverageValues(4)
    AverageValues(5) = BGKMovingAverageFac*u0ij(3,3) + (1.-BGKMovingAverageFac)*AverageValues(5)
    AverageValues(6) = BGKMovingAverageFac*u0ij(1,2) + (1.-BGKMovingAverageFac)*AverageValues(6)
    AverageValues(7) = BGKMovingAverageFac*u0ij(1,3) + (1.-BGKMovingAverageFac)*AverageValues(7)
    AverageValues(8) = BGKMovingAverageFac*u0ij(2,3) + (1.-BGKMovingAverageFac)*AverageValues(8)
  END IF
  dens = AverageValues(1)  
  u2 = AverageValues(2)
  u0ij(1,1)=AverageValues(3); u0ij(2,2)=AverageValues(4); u0ij(3,3)=AverageValues(5);
  u0ij(1,2)=AverageValues(6); u0ij(1,3)=AverageValues(7); u0ij(2,3)=AverageValues(8);
ELSE IF (BGKCollModel.EQ.2) THEN
  IF (AverageValues(1).EQ.0.0) THEN
    AverageValues(1) = dens
    AverageValues(2) = u2
    AverageValues(3:5) = u2i(1:3)
  ELSE
    AverageValues(1) = BGKMovingAverageFac*dens + (1.-BGKMovingAverageFac)*AverageValues(1)
    AverageValues(2) = BGKMovingAverageFac*u2 + (1.-BGKMovingAverageFac)*AverageValues(2)
    AverageValues(3:5) = BGKMovingAverageFac*u2i(1:3) + (1.-BGKMovingAverageFac)*AverageValues(3:5)
  END IF
  dens = AverageValues(1)
  u2 = AverageValues(2)
  u2i(1:3) = AverageValues(3:5)
ELSE IF (BGKCollModel.EQ.3) THEN
  IF (AverageValues(1).EQ.0.0) THEN
    AverageValues(1) = dens
    AverageValues(2) = u2
  ELSE
    AverageValues(1) = BGKMovingAverageFac*dens + (1.-BGKMovingAverageFac)*AverageValues(1)
    AverageValues(2) = BGKMovingAverageFac*u2 + (1.-BGKMovingAverageFac)*AverageValues(2)
  END IF
  dens = AverageValues(1)
  u2 = AverageValues(2)
END IF
CellTemp = Species(1)%MassIC * u2 / (3.0*BoltzmannConst)

END SUBROUTINE DoAveraging


SUBROUTINE CalcInnerDOFs(nSpec, EVibSpec, ERotSpec, totalWeightSpec, TVibSpec, TRotSpec, InnerDOF, Xi_VibSpec, Xi_RotSpec)
!===================================================================================================================================
!> Determine the internal degrees of freedom and the respective temperature (rotation/vibration) for diatomic/polyatomic species
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: SpecDSMC
USE MOD_BGK_Vars               ,ONLY: BGKDoVibRelaxation
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Analyze_Tools ,ONLY: CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nSpec(nSpecies)
REAL, INTENT(IN)              :: EVibSpec(nSpecies), ERotSpec(nSpecies), totalWeightSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: TVibSpec(nSpecies), TRotSpec(nSpecies), InnerDOF, Xi_VibSpec(nSpecies)
REAL, INTENT(OUT)             :: Xi_RotSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec
!===================================================================================================================================
Xi_VibSpec=0.; InnerDOF=0.; Xi_RotSpec=0.; TVibSpec=0.; TRotSpec=0.
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).EQ.0) CYCLE
  ! Only for molecules
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(BGKDoVibRelaxation) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN ! polyatomic
        ! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
        TVibSpec(iSpec) = CalcTVibPoly(EVibSpec(iSpec)/totalWeightSpec(iSpec) + SpecDSMC(iSpec)%EZeroPoint, iSpec)
      ELSE ! diatomic
        ! Calculation of vibrational temperature and DOFs from Pfeiffer, Physics of Fluids 30, 116103 (2018), "Extending the
        ! particle ellipsoidal statistical Bhatnagar-Gross-Krook method to diatomic molecules including quantized vibrational
        ! energies"
        ! TVibSpec = vibrational energy without zero-point energy
        TVibSpec(iSpec) = EVibSpec(iSpec) / (totalWeightSpec(iSpec)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
        IF (TVibSpec(iSpec).GT.0.0) TVibSpec(iSpec) = SpecDSMC(iSpec)%CharaTVib/LOG(1. + 1./(TVibSpec(iSpec)))
      END IF
      IF (TVibSpec(iSpec).GT.0.0) Xi_VibSpec(iSpec) = 2.* EVibSpec(iSpec) / (totalWeightSpec(iSpec)*BoltzmannConst*TVibSpec(iSpec))
    END IF
    Xi_RotSpec(iSpec) = SpecDSMC(iSpec)%Xi_Rot
    ! Calculation of rotational temperature from Pfeiffer, Physics of Fluids 30, 116103 (2018), "Extending the particle ellipsoidal
    ! statistical Bhatnagar-Gross-Krook method to diatomic molecules including quantized vibrational energies"
    TRotSpec(iSpec) = 2.*ERotSpec(iSpec)/(Xi_RotSpec(iSpec)*totalWeightSpec(iSpec)*BoltzmannConst)
  END IF
  InnerDOF = InnerDOF + Xi_RotSpec(iSpec)  + Xi_VibSpec(iSpec)
END DO

END SUBROUTINE CalcInnerDOFs


SUBROUTINE CalcGasProperties(nSpec, dens, InnerDOF, totalWeightSpec, totalWeight, TotalMass, u2Spec, SpecTemp, CellTemp, &
    Xi_VibSpec, Xi_RotSpec, Prandtl, relaxfreq, dynamicvis, thermalcond, MassIC_Mixture)
!===================================================================================================================================
!> Calculate the reference dynamic viscosity, Prandtl number and the resulting relaxation frequency of the distribution function
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: Species, nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, CollInf, DSMC
USE MOD_BGK_Vars              ,ONLY: BGK_ExpectedPrandtlNumber, BGKMixtureModel, BGKCollModel
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst, Pi
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nSpec(nSpecies)
REAL, INTENT(IN)              :: totalWeightSpec(nSpecies), totalWeight, TotalMass, u2Spec(nSpecies), SpecTemp(nSpecies), CellTemp
REAL, INTENT(IN)              :: Xi_VibSpec(nSpecies), Xi_RotSpec(nSpecies), dens, InnerDOF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: Prandtl, relaxfreq, dynamicvis, thermalcond, MassIC_Mixture
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec, jSpec
REAL                          :: MolarFraction(1:nSpecies), DOFFraction(1:nSpecies), MassFraction(1:nSpecies)
REAL                          :: PrandtlCorrection, dynamicvisSpec(nSpecies), thermalcondSpec(nSpecies), Phi(nSpecies)
REAL                          :: TotalDOFWeight, C_P, CellTempSpec(nSpecies+1)
!===================================================================================================================================
MassIC_Mixture = TotalMass / totalWeight
thermalcond=0.
IF (nSpecies.GT.1) THEN ! gas mixture
  MolarFraction(1:nSpecies) = totalWeightSpec(1:nSpecies) / totalWeight
  MassFraction(1:nSpecies) = MolarFraction(1:nSpecies) * Species(1:nSpecies)%MassIC / MassIC_Mixture
  DOFFraction(1:nSpecies) = totalWeightSpec(1:nSpecies) * (5.+Xi_RotSpec(1:nSpecies)+Xi_VibSpec(1:nSpecies))
  TotalDOFWeight = SUM(DOFFraction)

  PrandtlCorrection = 0.
  C_P = 0.0
  DO iSpec = 1, nSpecies
    IF (nSpec(iSpec).EQ.0) CYCLE
    ! Correction of Pr for calculation of relaxation frequency, see alpha - F. Hild, M. Pfeiffer, "Multi-species modeling in the
    ! particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method including internal degrees of freedom", subitted to Phys.
    ! Fluids, August 2023
    PrandtlCorrection = PrandtlCorrection + DOFFraction(iSpec)*MassIC_Mixture/Species(iSpec)%MassIC/TotalDOFWeight
    C_P = C_P + ((5. + (Xi_VibSpec(iSpec)+Xi_RotSpec(iSpec)))/2.) * BoltzmannConst / Species(iSpec)%MassIC * MassFraction(iSpec) 
  END DO

  SELECT CASE(BGKMixtureModel)
  ! Both cases are described in Pfeiffer et. al., Physics of Fluids 33, 036106 (2021),
  ! "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method for monatomic gas species"
  ! Extension to poyatomic mixtures according to F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal
  ! statistical Bhatnagar-Gross-Krook method including internal degrees of freedom", subitted to Phys. Fluids, August 2023

  CASE (1)  ! Wilke's mixing rules
    DO iSpec = 1, nSpecies
      ! Dynamic viscosity per species
      ! Omega = OmegaVHS + 0.5
      IF ((nSpec(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(u2Spec(iSpec)))) THEN
        ! Species temperature: Sufficient number of particles per species are available
        dynamicvisSpec(iSpec) = 30.*SQRT(Species(iSpec)%MassIC* BoltzmannConst*CollInf%Tref(iSpec,iSpec)/Pi) &
              /(4.*(4.- 2.*CollInf%omega(iSpec,iSpec)) * (6. - 2.*CollInf%omega(iSpec,iSpec))* CollInf%dref(iSpec,iSpec)**(2.) &
              *CollInf%Tref(iSpec,iSpec)**(CollInf%omega(iSpec,iSpec) + 0.5)*SpecTemp(iSpec)**(-CollInf%omega(iSpec,iSpec) - 0.5))
      ELSE
        ! Cell temperature: Low particle number case
        dynamicvisSpec(iSpec) = 30.*SQRT(Species(iSpec)%MassIC* BoltzmannConst*CollInf%Tref(iSpec,iSpec)/Pi) &
              /(4.*(4.- 2.*CollInf%omega(iSpec,iSpec)) * (6. - 2.*CollInf%omega(iSpec,iSpec))* CollInf%dref(iSpec,iSpec)**(2.) &
              *CollInf%Tref(iSpec,iSpec)**(CollInf%omega(iSpec,iSpec) + 0.5)*CellTemp**(-CollInf%omega(iSpec,iSpec) - 0.5))
      END IF
      ! Thermal conductivity per species (Eucken's formula with a correction by Hirschfelder for the internal degrees of freedom)
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN ! inner DOF
        ! Istomin et. al., "Eucken correction in high-temperature gases with electronic excitation", J. Chem. Phys. 140,
        ! 184311 (2014)
        thermalcondspec(iSpec) = 0.25 * (15. + 2. * (Xi_VibSpec(iSpec)+Xi_RotSpec(iSpec))  * 1.328) &
                                            * dynamicvisSpec(iSpec) * BoltzmannConst / Species(iSpec)%MassIC
      ELSE ! atoms
        thermalcondspec(iSpec) = 0.25 * 15. * dynamicvisSpec(iSpec) * BoltzmannConst / Species(iSpec)%MassIC
      END IF
    END DO
    Phi= 0.0
    ! Calculation of factor phi, depending on mass ratios and ratios of dynamic viscosities
    DO iSpec = 1, nSpecies
      DO jSpec = 1, nSpecies
        Phi(iSpec) =  Phi(iSpec) &
                   + REAL(totalWeightSpec(jSpec)) &
                   * ( 1.0+SQRT(dynamicvisSpec(iSpec)/dynamicvisSpec(jSpec)) &
                   * (Species(jSpec)%MassIC/Species(iSpec)%MassIC)**(0.25) )**(2.0) &
                   / ( SQRT(8.0 * (1.0 + Species(iSpec)%MassIC/Species(jSpec)%MassIC)) )
      END DO
    END DO
    dynamicvis = 0.0
    thermalcond = 0.0
    ! Sum up dynamic viscosities and thermal conductivities of species
    DO iSpec = 1, nSpecies
      IF (nSpec(iSpec).EQ.0) CYCLE
      dynamicvis = dynamicvis + REAL(totalWeightSpec(iSpec)) * dynamicvisSpec(iSpec) / Phi(iSpec)
      thermalcond = thermalcond + REAL(totalWeightSpec(iSpec)) * thermalcondspec(iSpec) / Phi(iSpec)
    END DO
  
  CASE(2) ! Collision integrals (VHS)
    DO iSpec = 1, nSpecies
      IF ((nSpec(iSpec).LT.2).OR.ALMOSTZERO(u2Spec(iSpec))) THEN
        CellTempSpec(iSpec) = CellTemp
      ELSE
        CellTempSpec(iSpec) = SpecTemp(iSpec)
      END IF
    END DO
    CellTempSpec(nSpecies+1) = CellTemp
    CALL CalcViscosityThermalCondColIntVHS(CellTempSpec(1:nSpecies+1), MolarFraction(1:nSpecies),dens, Xi_RotSpec, Xi_VibSpec, &
      dynamicvis, thermalcond)
  END SELECT

  ! Calculation of Prandtl number
  Prandtl = C_P*dynamicvis/thermalcond*PrandtlCorrection

  IF(DSMC%CalcQualityFactors) BGK_ExpectedPrandtlNumber = BGK_ExpectedPrandtlNumber + Prandtl

  ! Calculation of relaxation frequency
  relaxfreq = Prandtl*dens*BoltzmannConst*CellTemp/dynamicvis

ELSE ! single species gas
  ! Calculation of reference dynamic viscosity
  dynamicvis = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*CollInf%Tref(1,1)/Pi) &
             / (4.*(4.- 2.*CollInf%omega(1,1)) * (6. - 2.*CollInf%omega(1,1))* CollInf%dref(1,1)**2.)
  ! Calculation of Prandtl number: Pr = cp*mu/K with inner DOF, atoms: Pr = 2/3
  Prandtl = 2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
  ! Calculation of relaxation frequency using the exponential ansatz of the viscosity mu
  IF (BGKCollModel.EQ.1) THEN ! ESBGK: relaxfreq nu = Pr*n*kB*T/mu
    relaxfreq = Prandtl*dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omega(1,1) + 0.5) &
        /dynamicvis*CellTemp**(-CollInf%omega(1,1) +0.5)
  ELSE ! relaxfreq nu = n*kB*T/mu
    relaxfreq = dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omega(1,1) + 0.5) &
        /dynamicvis*CellTemp**(-CollInf%omega(1,1) +0.5)
  END IF
END IF

END SUBROUTINE CalcGasProperties


SUBROUTINE CalcTRelax(ERotSpec, Xi_RotSpec, Xi_VibSpec, EVibSpec, totalWeightSpec, totalWeight, totalWeight2, dtCell, CellTemp, &
    TRotSpec, TVibSpec, relaxfreq, rotrelaxfreqSpec, vibrelaxfreqSpec, Xi_VibSpecNew, Xi_vib_DOF, nXiVibDOF, CellTempRel, TEqui, &
    betaR, betaV)
!===================================================================================================================================
!> Calculate the relaxation energies and temperatures
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: PolyatomMolDSMC, SpecDSMC, DSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Globals               ,ONLY: abort
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nXiVibDOF
REAL, INTENT(IN)              :: TRotSpec(nSpecies), ERotSpec(nSpecies), Xi_RotSpec(nSpecies)
REAL, INTENT(IN)              :: TVibSpec(nSpecies), EVibSpec(nSpecies), Xi_VibSpec(nSpecies)
REAL, INTENT(IN)              :: totalWeightSpec(nSpecies), totalWeight, totalWeight2, CellTemp, dtCell
REAL, INTENT(IN)              :: relaxfreq, rotrelaxfreqSpec(nSpecies), vibrelaxfreqSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: Xi_vib_DOF(DSMC%NumPolyatomMolecs,nXiVibDOF), Xi_VibSpecNew(nSpecies), CellTempRel, TEqui
REAL, INTENT(OUT)             :: betaR(nSpecies), betaV(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec, iDOF, iPolyatMole, i, j
REAL                          :: RotFracSpec(nSpecies), VibFracSpec(nSpecies)
REAL                          :: ERotSpecMean(nSpecies), EVibSpecMean(nSpecies), ETransRelMean
REAL                          :: ERotTtransSpecMean(nSpecies), EVibTtransSpecMean(nSpecies)
REAL                          :: TEqui_Old, TEqui_Old2, TEquiNum, TEquiDenom, exparg
REAL                          :: eps_prec=1.0E-0
!===================================================================================================================================
! According to J. Mathiaud et. al., "An ES-BGK model for diatomic gases with correct relaxation rates for internal energies",
! European Journal of Mechanics - B/Fluids, 96, pp. 65-77, 2022
! For implentation, see F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-
! Krook method including internal degrees of freedom", subitted to Phys. Fluids, August 2023

RotFracSpec=0.0; VibFracSpec=0.0; Xi_vib_DOF=0.0; Xi_VibSpecNew=0.0; betaR=1.0; betaV=1.0
ERotSpecMean=0.0; ERotTtransSpecMean=0.0; EVibSpecMean=0.0; EVibTtransSpecMean=0.0; ETransRelMean=0.0; CellTempRel=0.0

DO iSpec = 1, nSpecies
  IF(totalWeightSpec(iSpec).GT.0.) THEN
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN ! molecules
      ! Mean rotational energy per particle of a species
      ERotSpecMean(iSpec) = ERotSpec(iSpec)/totalWeightSpec(iSpec)
      ! Mean rotational energy per particle of a species for the mixture translational temperature, ERot(Ttrans)
      ERotTtransSpecMean(iSpec) = CellTemp * Xi_RotSpec(iSpec) * BoltzmannConst / 2.
      ! Calculate number of rotational relaxing molecules
      RotFracSpec(iSpec) = totalWeightSpec(iSpec)*(rotrelaxfreqSpec(iSpec)/relaxfreq)*(1.-EXP(-relaxfreq*dtCell))

      IF(BGKDoVibRelaxation) THEN
        ! Mean vibrational energy per particle of a species
        EVibSpecMean(iSpec) = EVibSpec(iSpec)/totalWeightSpec(iSpec)
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN ! polyatomic
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          ! Loop over all vibrational DOF
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            ! Mean vibrational energy per particle of a species for the mixture translational temperature, EVib(Ttrans)
            EVibTtransSpecMean(iSpec) = EVibTtransSpecMean(iSpec) + BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
              / (EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/CellTemp) - 1.)
          END DO
        ELSE ! diatomic
          ! Mean vibrational energy per particle of a species for the mixture translational temperature, EVib(Ttrans)
          EVibTtransSpecMean(iSpec) = BoltzmannConst * SpecDSMC(iSpec)%CharaTVib / (EXP(SpecDSMC(iSpec)%CharaTVib/CellTemp) - 1.)
        END IF
        ! Calculate number of vibrational relaxing molecules
        VibFracSpec(iSpec) = totalWeightSpec(iSpec)*(vibrelaxfreqSpec(iSpec)/relaxfreq)*(1.-EXP(-relaxfreq*dtCell))
      END IF

      ! Mean translational energy per particle to satisfy the Landau-Teller equation
      ETransRelMean = ETransRelMean + (3./2. * BoltzmannConst * CellTemp - (rotrelaxfreqSpec(iSpec)/relaxfreq) * &
        (ERotTtransSpecMean(iSpec)-ERotSpecMean(iSpec))) * totalWeightSpec(iSpec)/totalWeight
      IF (BGKDoVibRelaxation) THEN
        ETransRelMean = ETransRelMean - (vibrelaxfreqSpec(iSpec)/relaxfreq)*(EVibTtransSpecMean(iSpec)-EVibSpecMean(iSpec)) * &
          totalWeightSpec(iSpec)/totalWeight
      END IF
    ELSE ! atomic
      ! Mean translational energy per particle to satisfy the Landau-Teller equation
      ETransRelMean = ETransRelMean + 3./2. * BoltzmannConst * CellTemp * totalWeightSpec(iSpec)/totalWeight
    END IF
  END IF
END DO

! Calculation of the cell temperature with ETransRelMean to satisfy the Landau-Teller equation
IF (ETransRelMean.GT.0.0) THEN
  CellTempRel = 2. * ETransRelMean / (3. * BoltzmannConst)
ELSE
  CellTempRel = CellTemp
END IF

! Calculation of equilibrium temperature for relaxation and energy conservation
TEqui_Old = 0.0
TEquiNum = 3.*(totalWeight - totalWeight2/totalWeight)*CellTemp
TEquiDenom = 3.*(totalWeight - totalWeight2/totalWeight)
! Sum up over all species
DO iSpec = 1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    TEquiNum = TEquiNum + Xi_RotSpec(iSpec)*RotFracSpec(iSpec)*TRotSpec(iSpec)
    TEquiDenom = TEquiDenom + Xi_RotSpec(iSpec)*RotFracSpec(iSpec)
    IF(BGKDoVibRelaxation) THEN
      TEquiNum = TEquiNum + Xi_VibSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
      TEquiDenom = TEquiDenom + Xi_VibSpec(iSpec)*VibFracSpec(iSpec)
    END IF
  END IF
END DO
TEqui = TEquiNum/TEquiDenom

i=0
! Solving of equation system until accuracy eps_prec is reached
outerLoop: DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  i = i + 1
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      ! if difference small: equilibrium, no beta
      IF (ABS(TRotSpec(iSpec)-TEqui).GT.1E-3) THEN
        betaR(iSpec) = (TRotSpec(iSpec)-CellTemp)/(TRotSpec(iSpec)-TEqui)
        IF (betaR(iSpec).LT.0.0) THEN
          betaR(iSpec) = 1.
        END IF
        ! new calculation of number of rotational relaxing molecules with betaR
        RotFracSpec(iSpec) = totalWeightSpec(iSpec)*(rotrelaxfreqSpec(iSpec)/relaxfreq)*(1.-EXP(-relaxfreq*dtCell))*betaR(iSpec)
      END IF
      IF(BGKDoVibRelaxation) THEN
        ! if difference small: equilibrium, no beta
        IF (ABS(TVibSpec(iSpec)-TEqui).GT.1E-3) THEN
          betaV(iSpec) = (TVibSpec(iSpec)-CellTemp)/(TVibSpec(iSpec)-TEqui)
          IF (betaV(iSpec).LT.0.0) THEN
            betaV(iSpec) = 1.
          END IF
          ! new calculation of number of vibrational relaxing molecules
          VibFracSpec(iSpec) = totalWeightSpec(iSpec)*(vibrelaxfreqSpec(iSpec)/relaxfreq)*(1.-EXP(-relaxfreq*dtCell))*betaV(iSpec)
        END IF

        ! new calculation of the vibrational degrees of freedom per species
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          ! Loop over all vibrational degrees of freedom to calculate them using TEqui
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            exparg = PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui
            ! Check if the exponent is within the range of machine precision for calculation of vibrational degrees of freedom
            IF(CHECKEXP(exparg)) THEN
              IF(exparg.gt.0.) THEN ! positive overflow: exp -> inf
                Xi_vib_DOF(iSpec,iDOF) = 2.*exparg/(EXP(exparg)-1.)
              ELSE ! negative overflow: exp -> 0
                Xi_vib_DOF(iSpec,iDOF) = 2.*exparg/(-1.)
              END IF ! exparg.gt.0.
            ELSE
              Xi_vib_DOF(iSpec,iDOF) = 0.0
            END IF ! CHECKEXP(exparg)
          END DO
          Xi_VibSpecNew(iSpec) = SUM(Xi_vib_DOF(iSpec,1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
        ELSE ! diatomic
          exparg = SpecDSMC(iSpec)%CharaTVib/TEqui
          ! Check if the exponent is within the range of machine precision for calculation of vibrational degrees of freedom
          IF(CHECKEXP(exparg)) THEN
            Xi_VibSpecNew(iSpec) = 2.*SpecDSMC(iSpec)%CharaTVib/TEqui/(EXP(exparg)-1.)
          ELSE
            Xi_VibSpecNew(iSpec) = 0.0
          END IF ! CHECKEXP(exparg)
        END IF
      END IF
    END IF
  END DO
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui
  ! new calculation of equilibrium temperature with new RotFracSpec, new VibFracSpec, new VibDOF(TEqui) in denominator
  TEquiNum = 3.*(totalWeight - totalWeight2/totalWeight)*CellTemp
  TEquiDenom = 3.*(totalWeight - totalWeight2/totalWeight)
  ! Sum up over all species
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      TEquiNum = TEquiNum + Xi_RotSpec(iSpec)*RotFracSpec(iSpec)*TRotSpec(iSpec)
      TEquiDenom = TEquiDenom + Xi_RotSpec(iSpec)*RotFracSpec(iSpec)
      IF(BGKDoVibRelaxation) THEN
        TEquiNum = TEquiNum + Xi_VibSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
        TEquiDenom = TEquiDenom + Xi_VibSpecNew(iSpec)*VibFracSpec(iSpec)
      END IF
    END IF
  END DO
  TEqui = TEquiNum/TEquiDenom
  IF(BGKDoVibRelaxation) THEN
    j=0
    ! accuracy eps_prec not reached yet
    innerLoop: DO WHILE ( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      j=j+1
      ! mean value of old and new equilibrium temperature
      TEqui = (TEqui + TEqui_Old2) * 0.5
      DO iSpec = 1, nSpecies
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          ! new calculation of the vibrational degrees of freedom per species
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            ! Loop over all vibrational degrees of freedom to calculate them using TEqui
            DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
              exparg = PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui
              ! Check if the exponent is within the range of machine precision for calculation of vibrational degrees of freedom
              IF(CHECKEXP(exparg)) THEN
                IF(exparg.gt.0.) THEN ! positive overflow: exp -> inf
                  Xi_vib_DOF(iSpec,iDOF) = 2.*exparg/(EXP(exparg)-1.)
                ELSE ! negative overflow: exp -> 0
                  Xi_vib_DOF(iSpec,iDOF) = 2.*exparg/(-1.)
                END IF ! exparg.gt.0.
              ELSE
                Xi_vib_DOF(iSpec,iDOF) = 0.0
              END IF ! CHECKEXP(exparg)
            END DO
            Xi_VibSpecNew(iSpec) = SUM(Xi_vib_DOF(iSpec,1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
          ELSE ! diatomic
            exparg = SpecDSMC(iSpec)%CharaTVib/TEqui
            ! Check if the exponent is within the range of machine precision for calculation of vibrational degrees of freedom
            IF(CHECKEXP(exparg)) THEN
              Xi_VibSpecNew(iSpec) = 2.*SpecDSMC(iSpec)%CharaTVib/TEqui/(EXP(exparg)-1.)
            ELSE
              Xi_VibSpecNew(iSpec) = 0.0
            END IF ! CHECKEXP(exparg)
          END IF
        END IF
      END DO
      TEqui_Old2 = TEqui
      ! new calculation of equilibrium temperature with new RotFracSpec, new VibFracSpec, new VibDOF(TEqui) in denominator
      TEquiNum = 3.*(totalWeight - totalWeight2/totalWeight)*CellTemp
      TEquiDenom = 3.*(totalWeight - totalWeight2/totalWeight)
      ! Sum up over all species
      DO iSpec = 1, nSpecies
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          TEquiNum = TEquiNum + Xi_RotSpec(iSpec)*RotFracSpec(iSpec)*TRotSpec(iSpec) + &
            Xi_VibSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
          TEquiDenom = TEquiDenom + Xi_RotSpec(iSpec)*RotFracSpec(iSpec) + Xi_VibSpecNew(iSpec)*VibFracSpec(iSpec)
        END IF
      END DO
      TEqui = TEquiNum/TEquiDenom
      IF (j.EQ.30) EXIT innerLoop
    END DO innerLoop
  END IF
  IF (i.EQ.30) EXIT outerLoop
END DO outerLoop

END SUBROUTINE CalcTRelax


SUBROUTINE DetermineRelaxPart(nPart, iPartIndx_Node, relaxfreq, dtCell, nRelax, nRotRelax, nVibRelax, &
    RotRelaxWeightSpec, VibRelaxWeightSpec, iPartIndx_NodeRelax, iPartIndx_NodeRelaxTemp, iPartIndx_NodeRelaxRot, &
    iPartIndx_NodeRelaxVib, vBulk, OldEnRot, OldEn, rotrelaxfreqSpec, vibrelaxfreqSpec, betaR, betaV)
!===================================================================================================================================
!> Determine the number of particles undergoing a relaxation (including vibration and rotation)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: Species, PartSpecies, PartState, nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, PartStateIntEn
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, iPartIndx_Node(nPart)
REAL, INTENT(IN)              :: relaxfreq, dtCell, rotrelaxfreqSpec(nSpecies), vibrelaxfreqSpec(nSpecies)
REAL, INTENT(IN)              :: betaR(nSpecies), betaV(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)          :: iPartIndx_NodeRelax(:), iPartIndx_NodeRelaxTemp(:)
INTEGER, INTENT(OUT)          :: iPartIndx_NodeRelaxRot(:), iPartIndx_NodeRelaxVib(:)
INTEGER, INTENT(OUT)          :: nRelax, nRotRelax, nVibRelax
REAL, INTENT(OUT)             :: vBulk(3), OldEnRot, RotRelaxWeightSpec(nSpecies), VibRelaxWeightSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT-OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: OldEn
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iSpec, iLoop, nNotRelax
REAL                          :: ProbAddPartTrans, iRan, partWeight, ProbAddPartRot(nSpecies), ProbAddPartVib(nSpecies)
!===================================================================================================================================
VibRelaxWeightSpec=0.0; RotRelaxWeightSpec=0.0; nRelax=0; nNotRelax=0; vBulk=0.0; nRotRelax=0; nVibRelax=0; OldEnRot=0.0

! Calculate probability of relaxation of a particle towards the target distribution function
ProbAddPartTrans = 1.-EXP(-relaxfreq*dtCell)
! Calculate probabilities of relaxation of a particle in the rotation and vibration
! See F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method
! including internal degrees of freedom", subitted to Phys. Fluids, August 2023
ProbAddPartRot(:) = ProbAddPartTrans * rotrelaxfreqSpec(:)/relaxfreq*betaR(:)
ProbAddPartVib(:) = ProbAddPartTrans * vibrelaxfreqSpec(:)/relaxfreq*betaV(:)

! Loop over all simulation particles
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  CALL RANDOM_NUMBER(iRan)
  ! Count particles that are undergoing a relaxation
  IF (ProbAddPartTrans.GT.iRan) THEN
    nRelax = nRelax + 1
    iPartIndx_NodeRelax(nRelax) = iPart
  ! Count particles that are not undergoing a relaxation
  ELSE
    nNotRelax = nNotRelax + 1
    iPartIndx_NodeRelaxTemp(nNotRelax) = iPart
    ! Sum up velocities of non-relaxing particles for bulk velocity
    vBulk(1:3) = vBulk(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  END IF

  ! For molecules: relaxation of inner DOF
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    ! Rotation
    CALL RANDOM_NUMBER(iRan)
    ! Count particles that are undergoing a relaxation, in total and per species
    IF (ProbAddPartRot(iSpec).GT.iRan) THEN
      nRotRelax = nRotRelax + 1
      RotRelaxWeightSpec(iSpec) = RotRelaxWeightSpec(iSpec) + partWeight
      iPartIndx_NodeRelaxRot(nRotRelax) = iPart
      ! Sum up total rotational energy
      OldEnRot = OldEnRot + PartStateIntEn(2,iPart) * partWeight
    END IF
    ! Vibration
    IF(BGKDoVibRelaxation) THEN
      CALL RANDOM_NUMBER(iRan)
      ! Count particles that are undergoing a relaxation, in total and per species
      IF (ProbAddPartVib(iSpec).GT.iRan) THEN
        nVibRelax = nVibRelax + 1
        VibRelaxWeightSpec(iSpec) = VibRelaxWeightSpec(iSpec) + partWeight
        iPartIndx_NodeRelaxVib(nVibRelax) = iPart
        ! Sum up total vibrational energy of all relaxing particles, considering zero-point energy, and add to translational energy
        OldEn = OldEn + (PartStateIntEn(1,iPartIndx_NodeRelaxVib(nVibRelax)) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
      END IF
    END IF
  END IF
END DO

END SUBROUTINE DetermineRelaxPart


SUBROUTINE RelaxInnerEnergy(nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib, iPartIndx_NodeRelaxRot, nXiVibDOF, Xi_vib_DOF, &
    Xi_VibSpec, Xi_RotSpec, VibEnergyDOF, TEqui, NewEnVib, NewEnRot)
!===================================================================================================================================
!> Determine the new rotational and vibrational energy of relaxing particles
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, PartStateIntEn, PolyatomMolDSMC, DSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nXiVibDOF
INTEGER, INTENT(IN)           :: nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib(nVibRelax), iPartIndx_NodeRelaxRot(nRotRelax)
REAL, INTENT(IN)              :: Xi_vib_DOF(DSMC%NumPolyatomMolecs,nXiVibDOF), Xi_VibSpec(nSpecies), Xi_RotSpec(nSpecies)
REAL, INTENT(IN)              :: TEqui
REAL, INTENT(INOUT)           :: NewEnVib(nSpecies), NewEnRot(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VibEnergyDOF(nVibRelax,nXiVibDOF)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iPart, iDOF, iPolyatMole, iSpec
REAL                          :: partWeight, iRan
!===================================================================================================================================
NewEnVib = 0.0; NewEnRot=0.0
IF(BGKDoVibRelaxation) THEN
  ! Loop over all particles undergoing a relaxation in the vibration
  DO iLoop = 1, nVibRelax
    iPart = iPartIndx_NodeRelaxVib(iLoop)
    iSpec = PartSpecies(iPart)
    partWeight = GetParticleWeight(iPart)
    ! polyatomic, more than one vibrational DOF
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      PartStateIntEn(1,iPart) = 0.0
      ! Sum up the new vibrational energy over all DOFs, see M. Pfeiffer et. al., "Extension of Particle-based BGK Models to
      ! Polyatomic Species in Hypersonic Flow around a Flat-faced Cylinder", AIP Conference Proceedings 2132, 100001 (2019)
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        CALL RANDOM_NUMBER(iRan)
        VibEnergyDOF(iLoop,iDOF) = - LOG(iRan)*Xi_vib_DOF(iPolyatMole,iDOF)/2.*TEqui*BoltzmannConst
        PartStateIntEn(1,iPart) = PartStateIntEn(1,iPart)+VibEnergyDOF(iLoop,iDOF)
      END DO
    ! ELSE: diatomic, only one vibrational DOF, calculate new vibrational energy according to M. Pfeiffer, "Extending the particle
    ! ellipsoidal statistical Bhatnagar-Gross-Krook method to diatomic molecules including quantized vibrational energies",
    ! Phys. Fluids 30, 116103 (2018)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn( 1,iPart) = -LOG(iRan)*Xi_VibSpec(iSpec)/2.*TEqui*BoltzmannConst
    END IF
    ! Sum up new vibrational energy per species
    NewEnVib(iSpec) = NewEnVib(iSpec) + PartStateIntEn(1,iPart) * partWeight
  END DO
END IF
! Loop over all particles undergoing a relaxation in the rotation
DO iLoop = 1, nRotRelax
  iPart = iPartIndx_NodeRelaxRot(iLoop)
  iSpec = PartSpecies(iPart)
  partWeight = GetParticleWeight(iPart)
  CALL RANDOM_NUMBER(iRan)
  ! Calculate new rotational energy according to M. Pfeiffer et. al., "Extension of Particle-based BGK Models to Polyatomic Species
  ! in Hypersonic Flow around a Flat-faced Cylinder", AIP Conference Proceedings 2132, 100001 (2019)
  PartStateIntEn( 2,iPart) = -Xi_RotSpec(iSpec) / 2. * BoltzmannConst*TEqui*LOG(iRan)
  ! Sum up new rotational energy per species
  NewEnRot(iSpec) = NewEnRot(iSpec) + PartStateIntEn( 2,iPart) * partWeight
END DO

END SUBROUTINE RelaxInnerEnergy


SUBROUTINE SampleFromTargetDistr(nRelax, iPartIndx_NodeRelax, Prandtl, u2, u0ij, u2i, vBulkAll, CellTempRel, CellTemp, vBulk, &
    MassIC_Mixture)
!===================================================================================================================================
!> Sample new particle velocities from the target distribution function, depending on the chosen model
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, PartState
USE MOD_BGK_Vars              ,ONLY: BGKCollModel, ESBGKModel
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Globals               ,ONLY: abort
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nRelax, iPartIndx_NodeRelax(:)
REAL, INTENT(IN)              :: Prandtl, u2, u0ij(3,3), u2i(3), vBulkAll(3), CellTempRel, CellTemp, MassIC_Mixture
REAL, INTENT(INOUT)           :: vBulk(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, fillMa1, fillMa2, INFO, iLoop, iSpec
REAL                          :: iRanPart(3, nRelax), A(3,3), KronDelta, SMat(3,3), W(3), Work(100), tempVelo(3), partWeight
!===================================================================================================================================
! According to M. Pfeiffer, "Particle-based fluid dynamics: Comparison of different Bhatnagar-Gross-Krook models and the direct
! simulation Monte Carlo method for hypersonic flows", Phys. Fluids 30, 106106 (2018) and F. Hild, M. Pfeiffer, "Multi-species
! modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method including internal degrees of freedom",
! subitted to Phys. Fluids, August 2023
IF (nRelax.GT.0) THEN
  SELECT CASE(BGKCollModel)
  CASE (1)  ! Ellipsoidal Statistical BGK
    IF (ESBGKModel.EQ.1) THEN
      ! Approximated solution
      DO fillMa1 = 1, 3
        DO fillMa2 = fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          ! Fill symmetric transformation matrix SMat with anisotopic matrix A = SS
          SMat(fillMa1, fillMa2) = KronDelta - (1.-Prandtl)/(2.*Prandtl) &
            *(3./u2*u0ij(fillMa1, fillMa2)-KronDelta)
        END DO
      END DO
      SMat(2,1)=SMat(1,2)
      SMat(3,1)=SMat(1,3)
      SMat(3,2)=SMat(2,3)
      ! Generate random normals for the sampling of new velocities of all relaxing particles
      CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
    ELSE
      DO fillMa1 = 1, 3
        DO fillMa2 = fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          ! Fill anisotopic matrix A
          A(fillMa1, fillMa2) = KronDelta*CellTempRel*BoltzmannConst/MassIC_Mixture - (1.-Prandtl)/Prandtl * &
            (CellTempRel/CellTemp*u0ij(fillMa1, fillMa2) - KronDelta*CellTempRel*BoltzmannConst/MassIC_Mixture)
        END DO
      END DO
      IF (ESBGKModel.EQ.2) THEN
        ! Exact solution
        ! Compute eigenvalues and eigenvectors of matrix A --> output: W is the array that contains the eigenvalues, A then contains
        ! the orthonormal eigenvectors of anisotropic matrix A
        CALL DSYEV('V','U',3,A,3,W,Work,100,INFO)
        SMat = 0.0
        IF (W(3).LT.0.0) THEN
          ! Due to ascending order of eigenvalues, all three eigenvalues are lower than zero here
          ! Fallback to Maxwell BGK
          SMat(1,1) = SQRT(BoltzmannConst*CellTempRel/MassIC_Mixture)
          SMat(2,2) = SQRT(BoltzmannConst*CellTempRel/MassIC_Mixture)
          SMat(3,3) = SQRT(BoltzmannConst*CellTempRel/MassIC_Mixture)
        ELSE
          ! At least W(3) is not negative
          ! Set negative eigenvalues to zero to get positive semidefinite matrix
          IF (W(1).LT.0.0) THEN
            W(1) = 0.0
            IF (W(2).LT.0.0) W(2) = 0.0
          END IF
          ! SMat with square roots of the eigenvalues as diagonal elements
          SMat(1,1) = SQRT(W(1))
          SMat(2,2) = SQRT(W(2))
          SMat(3,3) = SQRT(W(3))
          ! Diagonalisation of anisotropic matrix, SMat is square root of anisotropic matrix
          SMat = MATMUL(A, SMat)
          SMat = MATMUL(SMat, TRANSPOSE(A))
        END IF
        ! Generate random normals for the sampling of new velocities of all relaxing particles
        CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
      ELSE IF (ESBGKModel.EQ.3) THEN
        ! Metropolis-Hastings
        A(2,1)=A(1,2)
        A(3,1)=A(1,3)
        A(3,2)=A(2,3)
        CALL MetropolisES(nRelax, iRanPart, A*MassIC_Mixture, CellTempRel)
      END IF
    END IF

  CASE (2)  ! Shakov BGK
    ! Acceptance-rejection method
    CALL ARShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)

  CASE (3)  ! Standard BGK (Maxwell target distribution)
    ! Generate random normals for the sampling of new velocities of all relaxing particles
    CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
  END SELECT
  
  ! Loop over all particles undergoing a relaxation towards the target distribution function
  DO iLoop = 1, nRelax
    iPart = iPartIndx_NodeRelax(iLoop)
    iSpec = PartSpecies(iPart)
    ! Calculation of new velocities of all particles
    IF ((BGKCollModel.EQ.1).AND.(ESBGKModel.EQ.1)) THEN
      ! Transformation of normalized thermal velocity vector tempVelo (sampled from a Maxwellian distribution) to a thermal velocity
      ! vector sampled from the ESBGK target distribution function (anisotropic Gaussian distribution)
      tempVelo(1:3) = SQRT(BoltzmannConst*CellTempRel/Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
      PartState(4:6,iPart) = vBulkAll(1:3) + MATMUL(SMat,tempVelo)
    ELSE IF ((BGKCollModel.EQ.1).AND.(ESBGKModel.EQ.2)) THEN
      tempVelo(1:3) = SQRT(MassIC_Mixture/Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
      PartState(4:6,iPart) = vBulkAll(1:3) + MATMUL(SMat,tempVelo)
    ELSE IF ((BGKCollModel.EQ.1).AND.(ESBGKModel.EQ.3)) THEN
      ! New thermal velocity (in x,y,z) of particle with mass scaling multiplied by normal distributed random vector
      PartState(4:6,iPart) = vBulkAll(1:3) + SQRT(1./Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
    ELSE
      ! New thermal velocity (in x,y,z) of particle is sqrt(k_B*T/m) multiplied by normal distributed random vector
      PartState(4:6,iPart) = vBulkAll(1:3) + SQRT(BoltzmannConst*CellTempRel/Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
    END IF
    partWeight = GetParticleWeight(iPart)
    ! Sum up new velocities of relaxing particles for bulk velocity, velocities of non-relaxing particles already calculated in
    ! subroutine DetermineRelaxPart
    vBulk(1:3) = vBulk(1:3) + PartState(4:6,iPart)*Species(iSpec)%MassIC*partWeight
  END DO
END IF ! nRelax.GT.0

END SUBROUTINE SampleFromTargetDistr


SUBROUTINE EnergyConsVib(nPart, totalWeight, nVibRelax, VibRelaxWeightSpec, iPartIndx_NodeRelaxVib, NewEnVib, OldEn, nXiVibDOF, &
    VibEnergyDOF, Xi_VibSpec, TEqui)
!===================================================================================================================================
!> Routine to ensure energy conservation when including vibrational degrees of freedom (continuous and quantized)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, nSpecies
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, PolyatomMolDSMC, VibQuantsPar, DSMC, SpecDSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation, BGKUseQuantVibEn
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, nXiVibDOF
INTEGER, INTENT(IN)           :: nVibRelax, iPartIndx_NodeRelaxVib(nPart)
REAL, INTENT(IN)              :: VibRelaxWeightSpec(nSpecies), Xi_VibSpec(nSpecies), totalWeight
REAL, INTENT(IN)              :: NewEnVib(nSpecies), VibEnergyDOF(nVibRelax,nXiVibDOF), TEqui!, EVibTtransSpecMean(nSpecies)
REAL, INTENT(INOUT)           :: OldEn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iLoop, iDOF, iSpec, iQuant, iQuaMax, iPolyatMole
REAL                          :: Xi_VibTotal, alpha(nSpecies), partWeight, betaV, iRan, MaxColQua
!===================================================================================================================================
! According to F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook
! method including internal degrees of freedom", subitted to Phys. Fluids, August 2023
IF(BGKDoVibRelaxation) THEN
  ! Vibrational energy is positive for at least one species + there are vibrational relaxations
  IF (ANY(NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)) THEN
    Xi_VibTotal = 0.0
    DO iSpec = 1, nSpecies
      ! Sum of relaxing vibrational degrees of freedom
      Xi_VibTotal = Xi_VibTotal + Xi_VibSpec(iSpec)*VibRelaxWeightSpec(iSpec)
    END DO
    ! Calculate scaling factor alpha per species
    ! EVibTtransSpecMean(iSpec)*VibRelaxWeightSpec(iSpec) is energy that should be in vibration
    DO iSpec = 1, nSpecies
      IF (NewEnVib(iSpec).GT.0.0) THEN
        alpha(iSpec) = OldEn/NewEnVib(iSpec)*(Xi_VibSpec(iSpec)*VibRelaxWeightSpec(iSpec)/(3.*(totalWeight-1.)+Xi_VibTotal))
      ELSE
        alpha(iSpec) = 0.
      END IF
    END DO
    ! Quantized vibrational energy
    IF (BGKUseQuantVibEn) THEN
      DO iLoop = 1, nVibRelax
        iPart = iPartIndx_NodeRelaxVib(iLoop)
        partWeight = GetParticleWeight(iPart)
        iSpec = PartSpecies(iPart)
        ! Polyatomic molecules
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          PartStateIntEn(1,iPart) = 0.0
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          ! Loop over all vibrational DOF
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            ! Energy per vibrational mode alpha*VibEnergyDOF is reformulated to a quantum number iQuant
            betaV = alpha(iSpec)*VibEnergyDOF(iLoop,iDOF)/(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst)
            CALL RANDOM_NUMBER(iRan)
            iQuant = INT(betaV+iRan)
            ! Check maximum vibrational quantum number
            IF(iQuant.GT.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)) iQuant=PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)
            ! Remaining energy negative, new quantum number needs to be calculated
            IF ((OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight).LT.0.0) THEN
              ! Maximum quantum number
              MaxColQua = OldEn/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*partWeight)
              ! OldEn < k_B*CharaTVibDOF --> iQuant < 1
              IF (INT(MaxColQua).EQ.0) THEN
                iQuant = 0
              ELSE
                CALL RANDOM_NUMBER(iRan)
                ! Calculation of new iQuant
                iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                ! Determine maximum quantum number
                iQuaMax = MIN(INT(MaxColQua)+1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
                ! Calculation of new iQuant as long as iQuant > maximum quantum number
                DO WHILE (iQuant.GE.iQuaMax)
                  CALL RANDOM_NUMBER(iRan)
                  iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                END DO
              END IF
            END IF
            ! Sum up the vibrational energy over all vibrational DOF
            PartStateIntEn( 1,iPart)  = PartStateIntEn( 1,iPart) &
               + iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            VibQuantsPar(iPart)%Quants(iDOF) = iQuant
            ! Remaining OldEn for remaining particles
            OldEn = OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight
          END DO
          ! Add zero-point energy
          PartStateIntEn( 1,iPart)  = PartStateIntEn( 1,iPart) + SpecDSMC(iSpec)%EZeroPoint
        ELSE  ! Diatomic molecules
          ! Vibrational energy is reformulated to a quantum number iQuant
          betaV = alpha(iSpec)*PartStateIntEn( 1,iPart)/(SpecDSMC(iSpec)%CharaTVib*BoltzmannConst)
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT(betaV+iRan)
          ! Check maximum vibrational quantum number
          IF (iQuant.GT.SpecDSMC(iSpec)%MaxVibQuant) iQuant = SpecDSMC(iSpec)%MaxVibQuant
          PartStateIntEn( 1,iPart)  = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
          ! Remaining energy negative, new quantum number needs to be calculated
          IF ((OldEn - (PartStateIntEn( 1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight).LT.0.0) THEN
            ! Maximum quantum number
            MaxColQua = OldEn/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib*partWeight)
            ! OldEn < k_B*CharaTVib --> iQuant < 1
            IF (INT(MaxColQua).EQ.0) THEN
              iQuant = 0
            ELSE
              CALL RANDOM_NUMBER(iRan)
              ! Calculation of new iQuant
              iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(iSpec)%CharaTVib)
              ! Determine maximum quantum number
              iQuaMax = MIN(INT(MaxColQua)+1, SpecDSMC(iSpec)%MaxVibQuant)
              ! Calculation of new iQuant as long as iQuant > maximum quantum number
              DO WHILE (iQuant.GE.iQuaMax)
                CALL RANDOM_NUMBER(iRan)
                iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(iSpec)%CharaTVib)
              END DO
            END IF
            ! Calculate vibrational energy including zero-point energy
            PartStateIntEn( 1,iPart)  = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
          END IF
          ! Remaining OldEn for remaining particles
          OldEn = OldEn - (PartStateIntEn( 1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
        END IF
      END DO
    ELSE ! Continuous treatment of vibrational energy
      DO iLoop = 1, nVibRelax
        iPart = iPartIndx_NodeRelaxVib(iLoop)
        iSpec = PartSpecies(iPart)
        partWeight = GetParticleWeight(iPart)
        ! Scaling of vibrational energy with factor alpha + zero-point energy
        PartStateIntEn( 1,iPart) = alpha(iSpec)*PartStateIntEn( 1,iPart) + SpecDSMC(iSpec)%EZeroPoint
        ! Remaining OldEn for remaining particles
        OldEn = OldEn - (PartStateIntEn( 1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
      END DO
    END IF ! BGKUseQuantVibEn
  ! NewEnVib = 0 for all species, relaxation towards the vibrational ground-state (new state is simply the zero-point energy)
  ELSE IF (nVibRelax.GT.0) THEN
    ! Set zero-point energy as vibrational energy for all particles with vibrational relaxation
    DO iLoop = 1, nVibRelax
      iPart = iPartIndx_NodeRelaxVib(iLoop)
      iSpec = PartSpecies(iPart)
      PartStateIntEn( 1,iPart) = SpecDSMC(iSpec)%EZeroPoint
    END DO
  END IF ! (NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)
END IF ! BGKDoVibRelaxation

END SUBROUTINE EnergyConsVib


SUBROUTINE MetropolisES(nPart, iRanPart, A, CellTempRel)
!===================================================================================================================================
!> Sampling from ESBGK target distribution function by using a Metropolis-Hastings method
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Basis                 ,ONLY: INV33
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: A(3,3), CellTempRel
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: iRanPartTemp(3), V2, iRan, NewProb, OldProb, NormProb, prefacor
INTEGER                        :: iLoop, iPart, iRun
LOGICAL                        :: Changed
REAL                           :: AC(3), AInvers(3,3), detA
!===================================================================================================================================
! Generate normal distributed random vector as start vector for the thermal velocity
prefacor = SQRT(BoltzmannConst*CellTempRel)
iRanPart(1,1) = rnor()*prefacor
iRanPart(2,1) = rnor()*prefacor
iRanPart(3,1) = rnor()*prefacor
! Inverse matrix of A
CALL INV33(A, AInvers, detA)
AC(1:3) = MATMUL(AInvers, iRanPart(1:3,1))
V2 = iRanPart(1,1)*AC(1) + iRanPart(2,1)*AC(2) + iRanPart(3,1)*AC(3)
OldProb = EXP(-0.5*V2)
! Burn-in phase, 35 initial steps
DO iLoop = 1, 35
  ! Generate normal distributed random vector for the thermal velocity
  iRanPartTemp(1) = rnor()*prefacor
  iRanPartTemp(2) = rnor()*prefacor
  iRanPartTemp(3) = rnor()*prefacor
  AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
  V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
  NewProb = EXP(-0.5*V2)
  NormProb = MIN(1.,NewProb/OldProb)
  CALL RANDOM_NUMBER(iRan)
  ! Acceptance of new sample with probability NormProb
  IF (NormProb.GT.iRan) THEN
    iRanPart(1:3,1) = iRanPartTemp(1:3)
    OldProb = NewProb
  END IF
END DO
! Main phase, for all following particles
DO iPart = 2, nPart
  ! Normal distributed random vector from previous particle
  iRanPart(1,iPart) = iRanPart(1,iPart-1)
  iRanPart(2,iPart) = iRanPart(2,iPart-1)
  iRanPart(3,iPart) = iRanPart(3,iPart-1)
  iRun = 0
  Changed = .FALSE.
  ! For acception: velocity should be changed at least once and at least ten steps in the Markov chain should be taken
  DO WHILE ((iRun.LT.10).OR.(.NOT.Changed))
    iRun = iRun + 1
    ! Generate normal distributed random vector for the thermal velocity
    iRanPartTemp(1) = rnor()*prefacor
    iRanPartTemp(2) = rnor()*prefacor
    iRanPartTemp(3) = rnor()*prefacor
    AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
    V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
    NewProb = EXP(-0.5*V2)
    NormProb = MIN(1.,NewProb/OldProb)
    CALL RANDOM_NUMBER(iRan)
    ! Acceptance of new sample with probability NormProb, velocity is changed
    IF (NormProb.GT.iRan) THEN
      Changed = .TRUE.
      iRanPart(1:3,iPart) = iRanPartTemp(1:3)
      OldProb = NewProb
    END IF
  END DO
END DO

END SUBROUTINE MetropolisES


SUBROUTINE ARShakhov(nPart, iRanPart, Vtherm, HeatVec, Prandtl)
!===================================================================================================================================
!> Acceptance-rejection method for sampling from the Shakhov distribution function
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Prandtl, Vtherm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vheat, V2, iRan, OldProb, Envelope
INTEGER                        :: iPart
!===================================================================================================================================
! Calculate envelope function
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope =  1.+4.*Envelope

! Loop over all relaxing particles
DO iPart = 1, nPart
  ! Generate random normals
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + (1.-Prandtl)*Vheat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
  CALL RANDOM_NUMBER(iRan)
  ! Acception if Envelope*iRan < OldProb
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb = (1. + (1.-Prandtl)*Vheat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARShakhov


SUBROUTINE BGK_BuildTransGaussNums(nPart, iRanPart)
!===================================================================================================================================
!> Generate normal distributed random vector for sampling of new velocities of all relaxing particles relaxing
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iLoop
!===================================================================================================================================
! Generate three normal distributed random values for all relaxing simulation particles
DO iLoop = 1, nPart
  iRanPart(1,iLoop) = rnor()
  iRanPart(2,iLoop) = rnor()
  iRanPart(3,iLoop) = rnor()
END DO

END SUBROUTINE BGK_BuildTransGaussNums


SUBROUTINE CalcViscosityThermalCondColIntVHS(CellTemp, Xi, dens, Xi_RotSpec, Xi_VibSpec, Visc, ThermalCond)
!===================================================================================================================================
!> Determination of the mixture viscosity and thermal conductivity using collision integrals (derived for the Variable Hard
!> Sphere model). Solving an equation system depending on the number of species.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : CollInf, SpecDSMC
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species, nSpecies
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp(nSpecies+1), Xi(nSpecies), dens, Xi_RotSpec(nSpecies), Xi_VibSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Visc,ThermalCond
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL        :: Sigma_11, Sigma_22, B_12(nSpecies,nSpecies), A_12(nSpecies,nSpecies), InteractDiam, cv, DiffCoef(nSpecies, nSpecies)
REAL        :: Mass, ViscSpec(nSpecies), ThermalCondSpec(nSpecies), TVHS, omegaVHS, E_12, CellTemptmp
REAL        :: ThermalCondSpec_Vib(nSpecies), ThermalCondSpec_Rot(nSpecies), cv_rot, cv_vib, rhoSpec
REAL        :: Xj_Dij(nSpecies,nSpecies), Xi_Dij_tot
REAL        :: ViscMat(nSpecies, nSpecies), RHSSolve(nSpecies), m0, pressure
INTEGER     :: iSpec, jSpec, kSpec, IPIV(nSpecies), info_dgesv
!===================================================================================================================================
ViscSpec = 0.; ThermalCondSpec = 0.; ThermalCondSpec_Vib = 0.; ThermalCondSpec_Rot = 0.; DiffCoef =0.; A_12 = 0.; B_12 = 0.
Xj_Dij = 0.; cv_rot = 0.; cv_vib = 0.; E_12 = 0.
! Loop over all species combinations
DO iSpec = 1, nSpecies
  IF (Xi(iSpec).LE.0.0) CYCLE
  ! Calculate cv with rotational and vibrational degrees of freedom
  IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    cv_rot = (Xi_RotSpec(iSpec)*BoltzmannConst)/(2.*Species(iSpec)%MassIC)
    cv_vib = (Xi_VibSpec(iSpec)*BoltzmannConst)/(2.*Species(iSpec)%MassIC)
  END IF
  DO jSpec = iSpec, nSpecies
    IF (Xi(jSpec).LE.0.0) CYCLE
    ! Interaction parameters
    InteractDiam = CollInf%dref(iSpec,jSpec)
    Mass = Species(iSpec)%MassIC*Species(jSpec)%MassIC/(Species(iSpec)%MassIC + Species(jSpec)%MassIC) ! reduced mass
    TVHS = CollInf%Tref(iSpec,jSpec)
    omegaVHS = CollInf%omega(iSpec,jSpec)
    IF (iSpec.EQ.jSpec) THEN
      CellTemptmp = CellTemp(iSpec) ! Species temperature or cell temperature for nSpec<2 or u2spec=0
    ELSE
      CellTemptmp = CellTemp(nSpecies+1) ! Cell temperature
    END IF
    ! Calculation of collision integral Sigma_22
    CALL CalcSigma_22VHS(CellTemptmp, InteractDiam, Mass, TVHS, omegaVHS, Sigma_22)
    IF (iSpec.EQ.jSpec) THEN
      cv= 3./2.*BoltzmannConst/(2.*Mass) ! DOF = 3, translational part
      ! Calculation of the viscosity and thermal conductivity
      ! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 160
      ViscSpec(iSpec) = (5./8.)*(BoltzmannConst*CellTemp(iSpec))/Sigma_22
      ThermalCondSpec(iSpec) = (25./16.)*(cv*BoltzmannConst*CellTemp(iSpec))/Sigma_22
      ! Results in the same as ThermalCondSpec(iSpec) = (15./4.)*BoltzmannConst/(2.*Mass)*ViscSpec(iSpec)
      ! Additional calculation of Sigma_11VHS and the diffusion coefficient for molecular species
      IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        CALL CalcSigma_11VHS(CellTemp(nSpecies+1), InteractDiam, Mass, TVHS, omegaVHS, Sigma_11)
        E_12 = BoltzmannConst*CellTemp(nSpecies+1)/(8.*Species(iSpec)%MassIC*Species(jSpec)%MassIC &
          /(Species(iSpec)%MassIC+Species(jSpec)%MassIC)**2.*Sigma_11)
        DiffCoef(iSpec,jSpec) = 3.*E_12/(2.*(Species(iSpec)%MassIC+Species(jSpec)%MassIC)*dens)
      END IF
    ELSE
      ! Calculation of collision integral Sigma_11
      CALL CalcSigma_11VHS(CellTemp(nSpecies+1), InteractDiam, Mass, TVHS, omegaVHS, Sigma_11)
      ! Parameters for calculation of contribution of species to mixture transport coefficients
      ! Pfeiffer et. al., Physics of Fluids 33, 036106 (2021), "Multi-species modeling in the particle-based ellipsoidal
      ! statistical Bhatnagar-Gross-Krook method for monatomic gas species"
      B_12(iSpec,jSpec) = (5.*GAMMA(4.-omegaVHS)-GAMMA(5.-omegaVHS))/(5.*GAMMA(3.-omegaVHS))
      B_12(jSpec,iSpec) = B_12(iSpec,jSpec)
      A_12(iSpec,jSpec) = Sigma_22 / (5.*Sigma_11)
      A_12(jSpec,iSpec) = A_12(iSpec,jSpec)
      E_12 = BoltzmannConst*CellTemp(nSpecies+1)/(8.*Species(iSpec)%MassIC*Species(jSpec)%MassIC &
          /(Species(iSpec)%MassIC+Species(jSpec)%MassIC)**2.*Sigma_11)
      DiffCoef(iSpec,jSpec) = 3.*E_12/(2.*(Species(iSpec)%MassIC+Species(jSpec)%MassIC)*dens)
      DiffCoef(jSpec,iSpec) = DiffCoef(iSpec,jSpec)
    END IF
    Xj_Dij(iSpec,jSpec) = Xi(jSpec)/DiffCoef(iSpec,jSpec)
    Xj_Dij(jSpec,iSpec) = Xj_Dij(iSpec,jSpec)
  END DO
  IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    ! Calculation of thermal conductivity of rotation and vibration for each molecular species
    ! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 254f
    ! F. Hild, M. Pfeiffer, "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method
    ! including internal degrees of freedom", subitted to Phys. Fluids, August 2023
    Xi_Dij_tot = SUM(Xj_Dij(iSpec,:))
    rhoSpec = dens * Species(iSpec)%MassIC * Xi(iSpec)
    ThermalCondSpec_Rot(iSpec) = (rhoSpec*cv_rot/Xi_Dij_tot)
    ThermalCondSpec_Vib(iSpec) = (rhoSpec*cv_vib/Xi_Dij_tot)
  END IF
END DO

! Calculate mixture viscosity by solving a system of linear equations with matrices
! Pfeiffer et. al., Physics of Fluids 33, 036106 (2021),
! "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method for monatomic gas species"
! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 352
ViscMat = 0.0
DO iSpec = 1, nSpecies
  IF (Xi(iSpec).LE.0.0) THEN
    ViscMat(iSpec,iSpec) = 1. ! Ensure invertibility of ViscMat
    CYCLE
  END IF
  DO jSpec = 1, nSpecies
    IF (Xi(jSpec).LE.0.0) CYCLE
    IF (iSpec.EQ.jSpec) THEN
      ViscMat(iSpec, jSpec) = Xi(iSpec)/ViscSpec(iSpec)
      DO kSpec = 1, nSpecies
        IF (Xi(kSpec).LE.0.0) CYCLE
        IF (kSpec.EQ.iSpec) CYCLE
        ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + 3.*Xi(kSpec) / ((Species(iSpec)%MassIC*dens &
          + Species(kSpec)%MassIC*dens)*DiffCoef(iSpec,kSpec))*(2./3.+Species(kSpec)%MassIC/Species(iSpec)%MassIC*A_12(iSpec,kSpec))
      END DO
    ELSE
      ViscMat(iSpec, jSpec) = -Xi(iSpec)*3. / ((Species(iSpec)%MassIC*dens &
        + Species(jSpec)%MassIC*dens)*DiffCoef(iSpec,jSpec))*(2./3.-A_12(iSpec,jSpec))
    END IF
  END DO
END DO
RHSSolve(:) = Xi(:)
CALL DGESV(nSpecies, 1, ViscMat, nSpecies, IPIV, RHSSolve, nSpecies, info_dgesv)
Visc = SUM(RHSSolve)

! Calculate mixture thermal conductivity by solving a system of linear equations with matrices
! Pfeiffer et. al., Physics of Fluids 33, 036106 (2021),
! "Multi-species modeling in the particle-based ellipsoidal statistical Bhatnagar-Gross-Krook method for monatomic gas species"
! S. Chapman and T.G. Cowling, "The mathematical Theory of Non-Uniform Gases", Cambridge University Press, 1970, S. 350f
pressure = BoltzmannConst*dens*CellTemp(nSpecies+1)
ViscMat = 0.0
DO iSpec = 1, nSpecies
  IF (Xi(iSpec).LE.0.0) THEN
    ViscMat(iSpec,iSpec) = 1. ! Ensure invertibility of ViscMat
    CYCLE
  END IF
  DO jSpec = 1, nSpecies
    IF (Xi(jSpec).LE.0.0) CYCLE
    IF (iSpec.EQ.jSpec) THEN
      ViscMat(iSpec, jSpec) = Xi(iSpec)/ThermalCondSpec(iSpec)
      DO kSpec = 1, nSpecies
        IF (Xi(kSpec).LE.0.0) CYCLE
        IF (kSpec.EQ.iSpec) CYCLE
        m0 = Species(iSpec)%MassIC+Species(kSpec)%MassIC
        ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + CellTemp(nSpecies+1)*Xi(kSpec)/(5.*pressure*DiffCoef(iSpec,kSpec)) &
          * (6.*Species(iSpec)%MassIC**2./m0**2.+(5.-4.*B_12(iSpec,kSpec))*Species(kSpec)%MassIC**2./m0**2. &
          + 8.*Species(iSpec)%MassIC*Species(kSpec)%MassIC/m0**2.*A_12(iSpec, kSpec))
      END DO
    ELSE
      m0 = Species(iSpec)%MassIC+Species(jSpec)%MassIC
      ViscMat(iSpec, jSpec) = -Xi(iSpec)*CellTemp(nSpecies+1) * (Species(iSpec)%MassIC*Species(jSpec)%MassIC/m0**2.) &
        /(5.*pressure*DiffCoef(iSpec,jSpec)) *(11.-4.*B_12(iSpec,jSpec)-8.*A_12(iSpec,jSpec))
    END IF
  END DO
END DO
RHSSolve(:) = Xi(:)
CALL DGESV(nSpecies, 1, ViscMat, nSpecies, IPIV, RHSSolve, nSpecies, info_dgesv)
! Thermal conductivity from translation, rotation and vibration
ThermalCond = SUM(RHSSolve) + SUM(ThermalCondSpec_Rot) + SUM(ThermalCondSpec_Vib)

END SUBROUTINE CalcViscosityThermalCondColIntVHS


SUBROUTINE CalcSigma_11VHS(CellTemp, Dref, Mass, Tref, omegaVHS, Sigma_11)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars            ,ONLY: Pi, BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, Dref, Mass, Tref, omegaVHS
REAL, INTENT(OUT)               :: Sigma_11
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Prefactor
!===================================================================================================================================
  ! See Stephani et. al., Physics of Fluids 24, 077101 (2012),
  ! “Consistent treatment of transport properties for five-species air direct simulation Monte Carlo/Navier-Stokes applications”
  Prefactor = Pi/2.*Dref*Dref*SQRT(BoltzmannConst/(2.*Pi*Mass))*Tref**omegaVHS*GAMMA(3.-omegaVHS)/GAMMA(2.-omegaVHS)
  Sigma_11 = Prefactor*CellTemp**(0.5-omegaVHS)

END SUBROUTINE CalcSigma_11VHS


SUBROUTINE CalcSigma_22VHS(CellTemp, Dref, Mass, Tref, omegaVHS, Sigma_22)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars            ,ONLY: Pi, BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, Dref, Mass, Tref, omegaVHS
REAL, INTENT(OUT)               :: Sigma_22
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: Prefactor
!===================================================================================================================================
  ! See Stephani et. al., Physics of Fluids 24, 077101 (2012),
  ! “Consistent treatment of transport properties for five-species air direct simulation Monte Carlo/Navier-Stokes applications”
  Prefactor = Pi/3.*Dref*Dref*SQRT(BoltzmannConst/(2.*Pi*Mass))*Tref**omegaVHS*GAMMA(4.-omegaVHS)/GAMMA(2.-omegaVHS)
  Sigma_22 = Prefactor*CellTemp**(0.5-omegaVHS)

END SUBROUTINE CalcSigma_22VHS

END MODULE MOD_BGK_CollOperator
