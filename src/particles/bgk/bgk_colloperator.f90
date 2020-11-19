!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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
PUBLIC :: BGK_CollisionOperator, ARShakhov, CalcTEquiPoly, CalcTEqui
!===================================================================================================================================

CONTAINS

SUBROUTINE BGK_CollisionOperator(iPartIndx_Node, nPart, NodeVolume, AveragingPara, CorrectStep)
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
!> 9.) Determine the new DSMC_RHS (for molecules, including rotational energy)
!> 9.) Scaling of the rotational energy of molecules
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies, usevMPF, VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: DSMC_RHS, SpecDSMC, DSMC, PartStateIntEn, PolyatomMolDSMC, RadialWeighting, CollInf
USE MOD_TimeDisc_Vars         ,ONLY: dt
USE MOD_BGK_Vars              ,ONLY: SpecBGK, BGKMovingAverageLength, BGKDoVibRelaxation
USE MOD_BGK_Vars              ,ONLY: BGK_MeanRelaxFactor, BGK_MeanRelaxFactorCounter, BGK_MaxRelaxFactor, BGK_MaxRotRelaxFactor
USE MOD_BGK_Vars              ,ONLY: BGK_PrandtlNumber
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
REAL, INTENT(INOUT), OPTIONAL           :: AveragingPara(5,BGKMovingAverageLength)
INTEGER, INTENT(INOUT), OPTIONAL        :: CorrectStep
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: vBulk(3), u0ij(3,3), u2, V_rel(3), dtCell
REAL                  :: alpha, CellTemp, dens, InnerDOF, NewEn, OldEn, Prandtl, relaxfreq, TEqui
INTEGER, ALLOCATABLE  :: iPartIndx_NodeRelax(:),iPartIndx_NodeRelaxTemp(:),iPartIndx_NodeRelaxRot(:),iPartIndx_NodeRelaxVib(:)
INTEGER               :: iLoop, nRelax, iPolyatMole
REAL, ALLOCATABLE     :: Xi_vib_DOF(:), VibEnergyDOF(:,:)
INTEGER               :: iSpec, nSpec(nSpecies), jSpec, nRotRelax, nVibRelax
REAL                  :: OldEnRot, NewEnRot, NewEnVib
REAL                  :: vBulkSpec(1:3, nSpecies), TotalMass, u2Spec(nSpecies), u2i(3), vBulkAll(3)
REAL                  :: SpecTemp(nSpecies)
#ifdef CODE_ANALYZE
REAL                  :: Energy_old,Energy_new,Momentum_old(3),Momentum_new(3)
INTEGER               :: iMom
REAL,PARAMETER        :: RelMomTol=1e-6  ! Relative tolerance applied to conservation of momentum before/after reaction
REAL,PARAMETER        :: RelEneTol=1e-12 ! Relative tolerance applied to conservation of energy before/after reaction
#endif /* CODE_ANALYZE */
REAL                  :: totalWeightSpec(nSpecies), totalWeight, partWeight, totalWeightSpec2(nSpecies), CellTemptmp
REAL                  :: EVibSpec(nSpecies), ERotSpec(nSpecies), Xi_VibSpec(nSpecies), Xi_RotSpec(nSpecies),Xi_Vib_oldSpec(nSpecies)
REAL                  :: TVibSpec(nSpecies), TRotSpec(nSpecies), RotExpSpec(nSpecies), VibExpSpec(nSpecies)
REAL                  :: collisionfreqSpec(nSpecies),rotrelaxfreqSpec(nSpecies), vibrelaxfreqSpec(nSpecies), Xi_RotTotal
INTEGER               :: nVibRelaxSpec(nSpecies), nRotRelaxSpec(nSpecies)
!===================================================================================================================================
#ifdef CODE_ANALYZE
! Momentum and energy conservation check: summing up old values
Momentum_new = 0.0; Momentum_old = 0.0; Energy_new = 0.0; Energy_old = 0.0
DO iLoop = 1, nPart
  iSpec = PartSpecies(iPartIndx_Node(iLoop))
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  Momentum_old(1:3) = Momentum_old(1:3) + PartState(4:6,iPartIndx_Node(iLoop))*Species(iSpec)%MassIC*partWeight
  Energy_old = Energy_old + (PartState(4,iPartIndx_Node(iLoop))**(2.) + PartState(5,iPartIndx_Node(iLoop))**(2.) &
           + PartState(6,iPartIndx_Node(iLoop))**(2.))*0.5*Species(iSpec)%MassIC*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    Energy_old = Energy_old + (PartStateIntEn(1,iPartIndx_Node(iLoop)) + PartStateIntEn(2,iPartIndx_Node(iLoop)))*partWeight
  END IF
END DO
#endif

IF(nPart.LE.2) RETURN

! 1.) Moment calculation: Summing up the relative velocities and their squares
CALL CalcMoments(nPart, iPartIndx_Node, nSpec, vBulkSpec, vBulkAll, totalWeight, totalWeightSpec, & 
    totalWeightSpec2, TotalMass,  u2, u2Spec, u0ij, u2i, OldEn, EVibSpec, ERotSpec, CellTemp, SpecTemp, dtCell)
IF((CellTemp.LE.0).OR.(MAXVAL(nSpec(:)).EQ.1).OR.(totalWeight.LE.0.0)) RETURN

IF(VarTimeStep%UseVariableTimeStep) THEN
  dtCell = dt * dtCell / totalWeight
ELSE
  dtCell = dt
END IF

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  ! totalWeight contains the weighted particle number
  dens = totalWeight / NodeVolume
ELSE
  dens = totalWeight * Species(1)%MacroParticleFactor / NodeVolume
END IF

! Calculation of the rotational and vibrational degrees of freedom for molecules
IF (nSpecies.EQ.1) THEN
  IF((SpecDSMC(1)%InterID.EQ.2).OR.(SpecDSMC(1)%InterID.EQ.20)) THEN
    IF(SpecDSMC(1)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(1)%SpecToPolyArray
      ALLOCATE(Xi_vib_DOF(PolyatomMolDSMC(iPolyatMole)%VibDOF))
      Xi_vib_DOF(:) = 0.
    END IF
  END IF
END IF

CALL CalcInnerDOFs(nSpec, EVibSpec, ERotSpec, totalWeightSpec, TVibSpec, TRotSpec, InnerDOF, Xi_VibSpec, Xi_Vib_oldSpec & 
    ,Xi_RotSpec)

! 2.) Calculation of the relaxation frequency of the distribution function towards the target distribution function
CALL CalcGasProperties(nSpec, dens, InnerDOF, totalWeightSpec, totalWeight, TotalMass, u2Spec, u0ij, u2, SpecTemp, CellTemp, &
    Xi_VibSpec, Xi_RotSpec, Prandtl, relaxfreq)

IF(DSMC%CalcQualityFactors) THEN
  BGK_MeanRelaxFactor         = BGK_MeanRelaxFactor + relaxfreq * dtCell
  BGK_MeanRelaxFactorCounter  = BGK_MeanRelaxFactorCounter + 1
  BGK_MaxRelaxFactor          = MAX(BGK_MaxRelaxFactor,relaxfreq*dtCell)
  BGK_PrandtlNumber           = BGK_PrandtlNumber + Prandtl
END IF

! 3.) Treatment of molecules: determination of the rotational and vibrational relaxation frequency using the collision frequency,
!     which is not the same as the relaxation frequency of distribution function, calculated above.
IF(ANY(SpecDSMC(:)%InterID.EQ.2).OR.ANY(SpecDSMC(:)%InterID.EQ.20)) THEN
  collisionfreqSpec = 0.0
  DO iSpec = 1, nSpecies
    DO jSpec = 1, nSpecies
      IF (iSpec.EQ.jSpec) THEN
        CellTemptmp = CellTemp !SpecTemp(iSpec)
      ELSE
        CellTemptmp = CellTemp
      END IF
      collisionfreqSpec(iSpec) = collisionfreqSpec(iSpec) + SpecBGK(iSpec)%CollFreqPreFactor(jSpec) * totalWeightSpec(iSpec)*totalWeightSpec(jSpec) &
              *Dens *CellTemptmp**(-CollInf%omega(iSpec,jSpec) +0.5) /(totalWeight*totalWeight)
    END DO
  END DO
  rotrelaxfreqSpec(:) = collisionfreqSpec(:) * DSMC%RotRelaxProb
  vibrelaxfreqSpec(:) = collisionfreqSpec(:) * DSMC%VibRelaxProb
  RotExpSpec=0.; VibExpSpec=0.

  IF(SpecDSMC(1)%PolyatomicMol) THEN
    CALL CalcTEquiPoly(nPart, CellTemp, TRotSpec(1), TVibSpec(1), Xi_vib_DOF, Xi_Vib_oldSpec(1), RotExpSpec(1), VibExpSpec(1), &
                        TEqui, rotrelaxfreqSpec(1), vibrelaxfreqSpec(1), dtCell)
    Xi_VibSpec(1) = SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
  ELSE
    CALL CalcTEquiMulti(nPart, nSpec, CellTemp, TRotSpec, TVibSpec, Xi_VibSpec, Xi_Vib_oldSpec, RotExpSpec, VibExpSpec,  &
      TEqui, rotrelaxfreqSpec, vibrelaxfreqSpec, dtCell)
  END IF
  IF(DSMC%CalcQualityFactors) THEN
    BGK_MaxRotRelaxFactor          = MAX(BGK_MaxRotRelaxFactor,MAXVAL(rotrelaxfreqSpec(:))*dtCell)
  END IF
END IF

! 4.) Determine the number of particles undergoing a relaxation (including vibration and rotation)
ALLOCATE(iPartIndx_NodeRelax(nPart), iPartIndx_NodeRelaxTemp(nPart))
iPartIndx_NodeRelax = 0; iPartIndx_NodeRelaxTemp = 0
ALLOCATE(iPartIndx_NodeRelaxRot(nPart),iPartIndx_NodeRelaxVib(nPart))
iPartIndx_NodeRelaxRot = 0; iPartIndx_NodeRelaxVib = 0

CALL DetermineRelaxPart(nPart, iPartIndx_Node, relaxfreq, dtCell, RotExpSpec, VibExpSpec, nRelax, nRotRelax, nVibRelax, &
    nRotRelaxSpec, nVibRelaxSpec, iPartIndx_NodeRelax, iPartIndx_NodeRelaxTemp, iPartIndx_NodeRelaxRot, &
    iPartIndx_NodeRelaxVib, vBulk, OldEnRot, OldEn)
IF ((nRelax.EQ.0).AND.(nRotRelax.EQ.0).AND.(nVibRelax.EQ.0)) RETURN

IF(BGKDoVibRelaxation) THEN
   IF(SpecDSMC(1)%PolyatomicMol) THEN
     ALLOCATE(VibEnergyDOF(nVibRelax,PolyatomMolDSMC(iPolyatMole)%VibDOF))
   END IF
END IF
! 5.) Determine the new rotational and vibrational state of molecules undergoing a relaxation
CALL RelaxInnerEnergy(nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib, iPartIndx_NodeRelaxRot, Xi_vib_DOF, Xi_VibSpec, &
    Xi_RotSpec , TEqui, VibEnergyDOF, NewEnVib, NewEnRot)

! 6.) Sample new particle velocities from the target distribution function, depending on the chosen model
CALL SampleFromTargetDistr(nRelax, ipartindx_noderelax, Prandtl, u2, u0ij, u2i, vBulkAll, CellTemp, vBulk)

NewEn = 0.
vBulk = vBulk/TotalMass
DO iLoop = 1, nRelax 
  iSpec = PartSpecies(iPartIndx_NodeRelax(iLoop))
  partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
  V_rel(1:3) = DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop)) - vBulk(1:3)
  NewEn = NewEn + (V_rel(1)**(2.) + V_rel(2)**(2.) + V_rel(3)**(2.))*0.5*Species(iSpec)%MassIC*partWeight
END DO
DO iLoop = 1, nPart-nRelax 
  iSpec = PartSpecies(iPartIndx_NodeRelaxTemp(iLoop))
  partWeight = GetParticleWeight(iPartIndx_NodeRelaxTemp(iLoop))
  V_rel(1:3) = PartState(4:6,iPartIndx_NodeRelaxTemp(iLoop)) - vBulk(1:3)
  NewEn = NewEn + (V_rel(1)**(2.) + V_rel(2)**(2.) + V_rel(3)**(2.))*0.5*Species(iSpec)%MassIC*partWeight
END DO

! 7.) Vibrational energy of the molecules: Ensure energy conservation by scaling the new vibrational states with the factor alpha
IF(ANY(SpecDSMC(:)%InterID.EQ.2).OR.ANY(SpecDSMC(:)%InterID.EQ.20)) THEN
  CALL EnergyConsVib(nPart, nVibRelax, nVibRelaxSpec, iPartIndx_NodeRelaxVib, NewEnVib, OldEn, Xi_VibSpec, VibEnergyDOF, TEqui)
END IF

OldEn = OldEn + OldEnRot
! 8.) Determine the new particle state and ensure energy conservation by scaling the new velocities with the factor alpha.
!     The actual update of particle velocity happens in the TimeDisc through the change in the velocity (DSMC_RHS)
Xi_RotTotal = 0.0
DO iSpec = 1, nSpecies
  Xi_RotTotal = Xi_RotTotal + Xi_RotSpec(iSpec)*nRotRelaxSpec(iSpec)
END DO
alpha = SQRT(OldEn/NewEn*(3.*(nPart-1.))/(Xi_RotTotal+3.*(nPart-1.)))
DO iLoop = 1, nRelax
  DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop)) = vBulkAll(1:3) + alpha*(DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop))-vBulk(1:3)) &
                      - PartState(4:6,iPartIndx_NodeRelax(iLoop))
END DO
DO iLoop = 1, nPart-nRelax
  DSMC_RHS(1:3,iPartIndx_NodeRelaxTemp(iLoop)) = vBulkAll(1:3) &
                      + alpha*(PartState(4:6,iPartIndx_NodeRelaxTemp(iLoop))-vBulk(1:3)) &
                      - PartState(4:6,iPartIndx_NodeRelaxTemp(iLoop))
END DO

! 9.) Rotation: Scale the new rotational state of the molecules to ensure energy conservation
IF ( (nRotRelax.GT.0)) alpha = OldEn/NewEnRot*(Xi_RotTotal/(Xi_RotTotal+3.*(nPart-1.)))
DO iLoop = 1, nRotRelax
  PartStateIntEn( 2,iPartIndx_NodeRelaxRot(iLoop)) = alpha*PartStateIntEn( 2,iPartIndx_NodeRelaxRot(iLoop))
END DO

! CODE ANALYZE: Compare the old momentum and energy of the cell with the new, abort if relative difference is above the limits
#ifdef CODE_ANALYZE
DO iLoop = 1, nPart
  iSpec = PartSpecies(iPartIndx_Node(iLoop))
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  Momentum_new(1:3) = Momentum_new(1:3) + (DSMC_RHS(1:3,iPartIndx_Node(iLoop)) + PartState(4:6,iPartIndx_Node(iLoop))) & 
          * Species(iSpec)%MassIC*partWeight
  Energy_new = Energy_new &
          + ((DSMC_RHS(1,iPartIndx_Node(iLoop)) + PartState(4,iPartIndx_Node(iLoop)))**(2.) &
          +  (DSMC_RHS(2,iPartIndx_Node(iLoop)) + PartState(5,iPartIndx_Node(iLoop)))**(2.) &
          +  (DSMC_RHS(3,iPartIndx_Node(iLoop)) + PartState(6,iPartIndx_Node(iLoop)))**(2.))*0.5*Species(iSpec)%MassIC*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    Energy_new = Energy_new + (PartStateIntEn(1,iPartIndx_Node(iLoop)) + PartStateIntEn(2,iPartIndx_Node(iLoop)))*partWeight
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

SUBROUTINE CalcMoments(nPart, iPartIndx_Node, nSpec, vBulkSpec, vBulkAll, totalWeight, totalWeightSpec, & 
    totalWeightSpec2, TotalMass,  u2, u2Spec, u0ij, u2i, OldEn, EVibSpec, ERotSpec, CellTemp, SpecTemp, dtCell)
!===================================================================================================================================
!> Moment calculation: Summing up the relative velocities and their squares
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, Species, PartSpecies, nSpecies, VarTimeStep
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, SpecDSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation, BGKCollModel
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart, iPartIndx_Node(nPart)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)          :: nSpec(nSpecies)
REAL, INTENT(OUT)             :: u2Spec(nSpecies),u0ij(3,3), OldEn, EVibSpec(nSpecies), ERotSpec(nSpecies), u2i(3), u2
REAL, INTENT(OUT)             :: CellTemp, SpecTemp(nSpecies), totalWeightSpec(nSpecies), totalWeightSpec2(nSpecies)
REAL, INTENT(OUT)             :: vBulkAll(3), vBulkSpec(3,nSpecies), totalWeight, TotalMass, dtCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iSpec, fillMa1, fillMa2
REAL                          :: V_rel(1:3), vmag2, partWeight, EnerTotal
REAL                          :: tempweight, tempweight2, tempmass, vBulkTemp(3), totalWeight2, totalWeight3
!===================================================================================================================================
totalWeightSpec = 0.0; totalWeightSpec2=0.0; vBulkAll=0.0; TotalMass=0.0; vBulkSpec=0.0; nSpec=0; dtCell=0.0
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  iSpec = PartSpecies(iPartIndx_Node(iLoop))
  totalWeightSpec(iSpec) = totalWeightSpec(iSpec) + partWeight
  totalWeightSpec2(iSpec) =   totalWeightSpec2(iSpec) + partWeight*partWeight
  vBulkAll(1:3)  =  vBulkAll(1:3) + PartState(4:6,iPartIndx_Node(iLoop))*Species(iSpec)%MassIC*partWeight
  TotalMass = TotalMass + Species(iSpec)%MassIC*partWeight
  vBulkSpec(1:3,iSpec) = vBulkSpec(1:3,iSpec) + PartState(4:6,iPartIndx_Node(iLoop))*partWeight
  nSpec(iSpec) = nSpec(iSpec) + 1
  IF(VarTimeStep%UseVariableTimeStep) THEN
    dtCell = dtCell + VarTimeStep%ParticleTimeStep(iPartIndx_Node(iLoop))*partWeight
  END IF
END DO
totalWeight = SUM(totalWeightSpec)
totalWeight2 = SUM(totalWeightSpec2)
IF ((MAXVAL(nSpec(:)).EQ.1).OR.(totalWeight.LE.0.0)) RETURN
vBulkAll(1:3) = vBulkAll(1:3) / TotalMass
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).GT.0) vBulkSpec(:,iSpec) = vBulkSpec(:,iSpec) /totalWeightSpec(iSpec)
END DO

totalWeight3 = 0.; u2Spec=0.0; u0ij=0.0; u2i=0.0; OldEn=0.0; EVibSpec=0.0; ERotSpec=0.0
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  iSpec = PartSpecies(iPartIndx_Node(iLoop))  
  V_rel(1:3)=PartState(4:6,iPartIndx_Node(iLoop))-vBulkSpec(1:3,iSpec)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2  
  u2Spec(iSpec) = u2Spec(iSpec) + vmag2*partWeight

  V_rel(1:3)=PartState(4:6,iPartIndx_Node(iLoop))-vBulkAll(1:3)  
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2 
  IF (BGKCollModel.EQ.1) THEN 
    DO fillMa1 =1, 3
      DO fillMa2 =fillMa1, 3
        u0ij(fillMa1, fillMa2)= u0ij(fillMa1, fillMa2) & 
            + V_rel(fillMa1)*V_rel(fillMa2)*Species(iSpec)%MassIC*partWeight
      END DO
    END DO
  END IF
  IF (BGKCollModel.EQ.2) THEN
    u2i(1:3) = u2i(1:3) + V_rel(1:3)*vmag2 * partWeight*Species(iSpec)%MassIC
    totalWeight3 = totalWeight3 + partWeight*partWeight*partWeight
  END IF
  OldEn = OldEn + 0.5*Species(iSpec)%MassIC * vmag2*partWeight
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(BGKDoVibRelaxation) THEN
      EVibSpec(iSpec) = EVibSpec(iSpec) + (PartStateIntEn(1,iPartIndx_Node(iLoop)) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
    END IF 
    ERotSpec(iSpec) = ERotSpec(iSpec) + PartStateIntEn(2,iPartIndx_Node(iLoop)) * partWeight
  END IF
END DO

u0ij = u0ij* totalWeight / (TotalMass*(totalWeight - totalWeight2/totalWeight))
IF (BGKCollModel.EQ.2)  THEN
  u2i = u2i*totalWeight**3/(TotalMass*(totalWeight**3-3.*totalWeight*totalWeight2+2.*totalWeight3))
END IF

IF (nSpecies.GT.1) THEN
  SpecTemp = 0.0
  EnerTotal = 0.0 
  tempweight = 0.0; tempweight2 = 0.0; tempmass = 0.0; vBulkTemp = 0.0
  DO iSpec = 1, nSpecies
    IF ((nSpec(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(u2Spec(iSpec)))) THEN
      SpecTemp(iSpec) = Species(iSpec)%MassIC * u2Spec(iSpec) &
          /(3.0*BoltzmannConst*(totalWeightSpec(iSpec) - totalWeightSpec2(iSpec)/totalWeightSpec(iSpec)))
      EnerTotal =  EnerTotal + 3./2.*BoltzmannConst*SpecTemp(iSpec) * totalWeightSpec(iSpec)
      vmag2 = vBulkSpec(1,iSpec)**(2.) + vBulkSpec(2,iSpec)**(2.) + vBulkSpec(3,iSpec)**(2.)
      EnerTotal = EnerTotal + totalWeightSpec(iSpec) * Species(iSpec)%MassIC / 2. * vmag2
      tempweight = tempweight + totalWeightSpec(iSpec)
      tempweight2 = tempweight2 + totalWeightSpec2(iSpec)
      tempmass = tempmass +  totalWeightSpec(iSpec) * Species(iSpec)%MassIC 
      vBulkTemp(1:3) = vBulkTemp(1:3) + vBulkSpec(1:3,iSpec)*totalWeightSpec(iSpec) * Species(iSpec)%MassIC 
    END IF   
  END DO

  vBulkTemp(1:3) = vBulkTemp(1:3) / tempmass
  vmag2 = vBulkTemp(1)*vBulkTemp(1) + vBulkTemp(2)*vBulkTemp(2) + vBulkTemp(3)*vBulkTemp(3)
  EnerTotal = EnerTotal -  tempmass / 2. * vmag2
  CellTemp = 2. * EnerTotal / (3.*tempweight*BoltzmannConst)
  u2 = 3. * CellTemp * BoltzmannConst * (tempweight - tempweight2/tempweight) / tempmass
ELSE
  u2 = u2Spec(1) / (totalWeight - totalWeight2/totalWeight)
  CellTemp = Species(1)%MassIC * u2 / (3.0*BoltzmannConst)
END IF

END SUBROUTINE CalcMoments

SUBROUTINE CalcInnerDOFs(nSpec, EVibSpec, ERotSpec, totalWeightSpec, TVibSpec, TRotSpec, InnerDOF, Xi_VibSpec, Xi_Vib_oldSpec & 
    ,Xi_RotSpec)
!===================================================================================================================================
!> Determine the internal degrees of freedom and the respective temperature (rotation/vibration) for diatomic/polyatomic species
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_DSMC_Analyze          ,ONLY: CalcTVibPoly
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nSpec(nSpecies)
REAL, INTENT(IN)              :: EVibSpec(nSpecies), ERotSpec(nSpecies), totalWeightSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: TVibSpec(nSpecies), TRotSpec(nSpecies), InnerDOF, Xi_VibSpec(nSpecies), Xi_Vib_oldSpec(nSpecies)
REAL, INTENT(OUT)             :: Xi_RotSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPolyatMole, iSpec, iDOF
!===================================================================================================================================
Xi_VibSpec=0.; InnerDOF=0.; Xi_RotSpec=0.; Xi_Vib_oldSpec=0.; TVibSpec=0.; TRotSpec=0.
DO iSpec = 1, nSpecies
  IF (nSpec(iSpec).EQ.0) CYCLE
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(BGKDoVibRelaxation) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN        
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        TVibSpec(iSpec) = CalcTVibPoly(EVibSpec(iSpec)/totalWeightSpec(iSpec), 1)
        IF (TVibSpec(iSpec).GT.0.0) THEN
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            Xi_VibSpec(iSpec) = Xi_VibSpec(iSpec) + 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVibSpec(iSpec) &
                                /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TVibSpec(iSpec)) - 1.)
          END DO
        END IF
      ELSE
        TVibSpec(iSpec)=EVibSpec(iSpec) / (totalWeightSpec(iSpec)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
        IF (TVibSpec(iSpec).GT.0.0) THEN
          TVibSpec(iSpec)= SpecDSMC(iSpec)%CharaTVib/LOG(1. + 1./(TVibSpec(iSpec)))
          Xi_VibSpec(iSpec) = 2.* EVibSpec(iSpec) / (totalWeightSpec(iSpec)*BoltzmannConst*TVibSpec(iSpec))
        END IF
      END IF
      Xi_Vib_oldSpec(iSpec) = Xi_VibSpec(iSpec)
    END IF
    Xi_RotSpec(iSpec) = SpecDSMC(iSpec)%Xi_Rot
    TRotSpec(iSpec) = 2.*ERotSpec(iSpec)/(Xi_RotSpec(iSpec)*totalWeightSpec(iSpec)*BoltzmannConst)    
  END IF
  InnerDOF = InnerDOF + Xi_RotSpec(iSpec)  + Xi_VibSpec(iSpec)
END DO

END SUBROUTINE CalcInnerDOFs

SUBROUTINE CalcGasProperties(nSpec, dens, InnerDOF, totalWeightSpec, totalWeight, TotalMass, u2Spec, u0ij, u2, SpecTemp, CellTemp, &
    Xi_VibSpec, Xi_RotSpec, Prandtl, relaxfreq)
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
REAL, INTENT(IN)              :: u0ij(3,3), u2, Xi_VibSpec(nSpecies), Xi_RotSpec(nSpecies), dens, InnerDOF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: Prandtl, relaxfreq
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSpec, jSpec, INFO
REAL                          :: MolarFraction(1:nSpecies), MassCoef, MassIC_Mixture, MassFraction(nSpecies)
REAL                          :: PrandtlCorrection, dynamicvisSpec(nSpecies), thermalcondSpec(nSpecies), Phi(nSpecies)
REAL                          :: dynamicvis, thermalcond, C_P, nu, A(3,3), W(3), Work(100), Theta, CellTempSpec(nSpecies+1)
!===================================================================================================================================
IF (nSpecies.GT.1) THEN
  MolarFraction(1:nSpecies) = totalWeightSpec(1:nSpecies) / totalWeight
  MassIC_Mixture = TotalMass / totalWeight
  MassCoef = 0.0
  DO iSpec = 1, nSpecies
    MassCoef=MassCoef + REAL(totalWeightSpec(iSpec))/REAL(totalWeight)*Species(iSpec)%MassIC
    MassFraction(iSpec) = MolarFraction(iSpec) * Species(iSpec)%MassIC / MassIC_Mixture
  END DO

  PrandtlCorrection = 0.
  DO iSpec = 1, nSpecies
    PrandtlCorrection = PrandtlCorrection + REAL(totalWeightSpec(iSpec))/REAL(totalWeight)*MassCoef/Species(iSpec)%MassIC
  END DO

  IF (BGKMixtureModel.EQ.1) THEN
    DO iSpec = 1, nSpecies
      IF ((nSpec(iSpec).GE.2).AND.(.NOT.ALMOSTZERO(u2Spec(iSpec)))) THEN
        ! Species temperature: Sufficient number of particles per species
        dynamicvisSpec(iSpec) = 30.*SQRT(Species(iSpec)%MassIC* BoltzmannConst*CollInf%Tref(iSpec,iSpec)/Pi) &
              /(4.*(4.- 2.*CollInf%omega(iSpec,iSpec)) * (6. - 2.*CollInf%omega(iSpec,iSpec))* CollInf%dref(iSpec,iSpec)**(2.) &
              *CollInf%Tref(iSpec,iSpec)**(CollInf%omega(iSpec,iSpec) + 0.5)*SpecTemp(iSpec)**(-CollInf%omega(iSpec,iSpec) - 0.5))
      ELSE
        ! Cell temperature: Low particle number case
        dynamicvisSpec(iSpec) = 30.*SQRT(Species(iSpec)%MassIC* BoltzmannConst*CollInf%Tref(iSpec,iSpec)/Pi) &
              /(4.*(4.- 2.*CollInf%omega(iSpec,iSpec)) * (6. - 2.*CollInf%omega(iSpec,iSpec))* CollInf%dref(iSpec,iSpec)**(2.) &
              *CollInf%Tref(iSpec,iSpec)**(CollInf%omega(iSpec,iSpec) + 0.5)*CellTemp**(-CollInf%omega(iSpec,iSpec) - 0.5))
      END IF
      ! innerdof pro spec !
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        thermalcondspec(iSpec) = 0.25 * (15. + 2. * (Xi_VibSpec(iSpec)+Xi_RotSpec(iSpec))  * 1.328) &
                                          * dynamicvisSpec(iSpec) &
                                          * BoltzmannConst / Species(iSpec)%MassIC
      ELSE
        thermalcondspec(iSpec) = 0.25 * 15.* dynamicvisSpec(iSpec) &
                                          * BoltzmannConst / Species(iSpec)%MassIC
      END IF
    END DO
    Phi= 0.0
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
    C_P = 0.0
    DO iSpec = 1, nSpecies
      IF (nSpec(iSpec).EQ.0) CYCLE
      dynamicvis = dynamicvis + REAL(totalWeightSpec(iSpec)) * dynamicvisSpec(iSpec) / Phi(iSpec)
      thermalcond = thermalcond + REAL(totalWeightSpec(iSpec)) * thermalcondspec(iSpec) / Phi(iSpec)
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        C_P = C_P + ((5. + (Xi_VibSpec(iSpec)+Xi_RotSpec(iSpec)))/2.) * BoltzmannConst / Species(iSpec)%MassIC * MassFraction(iSpec)
      ELSE
        C_P = C_P + (5./2.) * BoltzmannConst / Species(iSpec)%MassIC * MassFraction(iSpec)
      END IF
    END DO
  ELSE IF (BGKMixtureModel.EQ.2) THEN
    C_P = 0.0
    DO iSpec = 1, nSpecies
      IF ((nSpec(iSpec).LT.2).OR.ALMOSTZERO(u2Spec(iSpec))) THEN
        CellTempSpec(iSpec) = CellTemp
      ELSE
        C_P = C_P + REAL(totalWeightSpec(iSpec))/REAL(totalWeight)*Species(iSpec)%MassIC
        CellTempSpec(iSpec) = SpecTemp(iSpec)
      END IF
    END DO
    C_P = 5./2.*BoltzmannConst/C_P
    CellTempSpec(nSpecies+1) = CellTemp
    CALL CalcViscosityThermalCondColIntVHS(CellTempSpec(1:nSpecies+1), totalWeightSpec(1:nSpecies)/totalWeight,dens, dynamicvis, thermalcond)
  END IF

  Prandtl = C_P*dynamicvis/thermalcond*PrandtlCorrection
  IF(DSMC%CalcQualityFactors) BGK_ExpectedPrandtlNumber = BGK_ExpectedPrandtlNumber + Prandtl
  A = u0ij
  CALL DSYEV('N','U',3,A,3,W,Work,100,INFO)
  Theta = u2 / 3.

  nu = 1.-1./Prandtl
  nu= MAX(nu,-Theta/(W(3)-Theta))
  Prandtl = 1./(1.-nu)
  relaxfreq = Prandtl*dens*BoltzmannConst*CellTemp/dynamicvis
ELSE
  dynamicvis = 30.*SQRT(Species(1)%MassIC* BoltzmannConst*CollInf%Tref(1,1)/Pi) &
             / (4.*(4.- 2.*CollInf%omega(1,1)) * (6. - 2.*CollInf%omega(1,1))* CollInf%dref(1,1)**2.)
  Prandtl =2.*(InnerDOF + 5.)/(2.*InnerDOF + 15.)
  IF (BGKCollModel.EQ.1) THEN
    relaxfreq = Prandtl*dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omega(1,1) + 0.5) &
        /dynamicvis*CellTemp**(-CollInf%omega(1,1) +0.5)
  ELSE
    relaxfreq = dens*BoltzmannConst*CollInf%Tref(1,1)**(CollInf%omega(1,1) + 0.5) &
        /dynamicvis*CellTemp**(-CollInf%omega(1,1) +0.5)
  END IF
END IF

END SUBROUTINE CalcGasProperties

SUBROUTINE DetermineRelaxPart(nPart, iPartIndx_Node, relaxfreq, dtCell, RotExpSpec, VibExpSpec, nRelax, nRotRelax, nVibRelax, &
    nRotRelaxSpec, nVibRelaxSpec, iPartIndx_NodeRelax, iPartIndx_NodeRelaxTemp, iPartIndx_NodeRelaxRot, &
    iPartIndx_NodeRelaxVib, vBulk, OldEnRot, OldEn)
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
REAL, INTENT(IN)              :: relaxfreq, dtCell, RotExpSpec(nSpecies), VibExpSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)          :: nRelax, iPartIndx_NodeRelax(nPart), iPartIndx_NodeRelaxTemp(nPart)
INTEGER, INTENT(OUT)          :: iPartIndx_NodeRelaxRot(nPart), iPartIndx_NodeRelaxVib(nPart)
INTEGER, INTENT(OUT)          :: nRotRelax, nVibRelax, nRotRelaxSpec(nSpecies), nVibRelaxSpec(nSpecies)
REAL, INTENT(OUT)             :: vBulk(3), OldEnRot
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT-OUTPUT VARIABLES
REAL, INTENT(INOUT)           :: OldEn
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: nNotRelax, iSpec, iLoop
REAL                          :: ProbAddPartTrans, iRan, partWeight
!===================================================================================================================================
nVibRelaxSpec =0; nRotRelaxSpec =0; nRelax=0; nNotRelax=0; vBulk=0.0; nRotRelax=0; nVibRelax=0; OldEnRot=0.0
ProbAddPartTrans = 1.-EXP(-relaxfreq*dtCell)
DO iLoop = 1, nPart  
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  CALL RANDOM_NUMBER(iRan)  
  IF (ProbAddPartTrans.GT.iRan) THEN
    nRelax = nRelax + 1
    iPartIndx_NodeRelax(nRelax) = iPartIndx_Node(iLoop)
  ELSE
    iSpec = PartSpecies(iPartIndx_Node(iLoop))  
    nNotRelax = nNotRelax + 1
    iPartIndx_NodeRelaxTemp(nNotRelax) = iPartIndx_Node(iLoop)
    vBulk(1:3) = vBulk(1:3) + PartState(4:6,iPartIndx_Node(iLoop))*Species(iSpec)%MassIC*partWeight
  END IF
  iSpec = PartSpecies(iPartIndx_Node(iLoop)) 
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    !Rotation
    CALL RANDOM_NUMBER(iRan)
    IF ((1.-RotExpSpec(iSpec)).GT.iRan) THEN
      nRotRelax = nRotRelax + 1
      nRotRelaxSpec(iSpec) = nRotRelaxSpec(iSpec) + 1
      iPartIndx_NodeRelaxRot(nRotRelax) = iPartIndx_Node(iLoop)
      OldEnRot = OldEnRot + PartStateIntEn(2,iPartIndx_Node(iLoop)) * partWeight
    END IF
    ! Vibration
    IF(BGKDoVibRelaxation) THEN
      CALL RANDOM_NUMBER(iRan)
      IF ((1.-VibExpSpec(iSpec)).GT.iRan) THEN
        nVibRelax = nVibRelax + 1
        nVibRelaxSpec(iSpec) = nVibRelaxSpec(iSpec) + 1
        iPartIndx_NodeRelaxVib(nVibRelax) = iPartIndx_Node(iLoop)
        OldEn = OldEn + (PartStateIntEn(1,iPartIndx_NodeRelaxVib(nVibRelax)) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
      END IF
    END IF
  END IF
END DO

END SUBROUTINE DetermineRelaxPart

SUBROUTINE RelaxInnerEnergy(nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib, iPartIndx_NodeRelaxRot, Xi_vib_DOF, Xi_VibSpec, &
    Xi_RotSpec , TEqui, VibEnergyDOF, NewEnVib, NewEnRot)
!===================================================================================================================================
!> Determine the new rotational and vibrational energy of relaxing particles
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, nSpecies
USE MOD_DSMC_Vars             ,ONLY: SpecDSMC, PartStateIntEn, PolyatomMolDSMC
USE MOD_BGK_Vars              ,ONLY: BGKDoVibRelaxation
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nVibRelax, nRotRelax, iPartIndx_NodeRelaxVib(nVibRelax), iPartIndx_NodeRelaxRot(nRotRelax)
REAL, INTENT(IN)              :: Xi_vib_DOF(:), TEqui, Xi_VibSpec(nSpecies), Xi_RotSpec(nSpecies)
REAL, INTENT(INOUT)           :: NewEnVib, NewEnRot
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VibEnergyDOF(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iDOF, iPolyatMole, iSpec
REAL                          :: partWeight, iRan
!===================================================================================================================================
! VIB Relaxation
NewEnVib = 0.0; NewEnRot=0.0
IF(BGKDoVibRelaxation) THEN
  DO iLoop = 1, nVibRelax
    iSpec = PartSpecies(iPartIndx_NodeRelaxVib(iLoop))
    partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
       iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
       PartStateIntEn(1,iPartIndx_NodeRelaxVib(iLoop)) = 0.0
       DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
         CALL RANDOM_NUMBER(iRan)
         VibEnergyDOF(iLoop,iDOF) = - LOG(iRan)*Xi_vib_DOF(iDOF)/2.*TEqui*BoltzmannConst
         PartStateIntEn(1,iPartIndx_NodeRelaxVib(iLoop)) = PartStateIntEn(1,iPartIndx_NodeRelaxVib(iLoop))+VibEnergyDOF(iLoop,iDOF)
       END DO
    ELSE      
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) = -LOG(iRan)*Xi_VibSpec(iSpec)/2.*TEqui*BoltzmannConst
    END IF
    NewEnVib = NewEnVib + PartStateIntEn(1,iPartIndx_NodeRelaxVib(iLoop)) * partWeight
  END DO
END IF
! ROT Relaxation
DO iLoop = 1, nRotRelax
  iSpec = PartSpecies(iPartIndx_NodeRelaxRot(iLoop)) 
  partWeight = GetParticleWeight(iPartIndx_NodeRelaxRot(iLoop))
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn( 2,iPartIndx_NodeRelaxRot(iLoop)) = -Xi_RotSpec(iSpec) / 2. * BoltzmannConst*TEqui*LOG(iRan)
  NewEnRot = NewEnRot + PartStateIntEn( 2,iPartIndx_NodeRelaxRot(iLoop)) * partWeight
END DO

END SUBROUTINE RelaxInnerEnergy

SUBROUTINE SampleFromTargetDistr(nRelax, ipartindx_noderelax, Prandtl, u2, u0ij, u2i, vBulkAll, CellTemp, vBulk)
!===================================================================================================================================
!> Sample new particle velocities from the target distribution function, depending on the chosen model
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species
USE MOD_DSMC_Vars             ,ONLY: DSMC_RHS
USE MOD_BGK_Vars              ,ONLY: BGKCollModel, ESBGKModel
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nRelax, ipartindx_noderelax(:)
REAL, INTENT(IN)              :: Prandtl, u2, u0ij(3,3), u2i(3), vBulkAll(3), CellTemp
REAL, INTENT(INOUT)           :: vBulk(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: fillMa1, fillMa2, INFO, iLoop, iSpec
REAL                          :: iRanPart(3, nRelax), A(3,3), KronDelta, SMat(3,3), W(3), Work(100), tempVelo(3), partWeight
!===================================================================================================================================
IF (nRelax.GT.0) THEN
  SELECT CASE(BGKCollModel)
  CASE (1)  ! Ellipsoidal Statistical
    IF (ESBGKModel.EQ.1) THEN
      !! Approximated Solution
      DO fillMa1 =1, 3
        DO fillMa2 =fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandtl)/(2.*Prandtl) &
            *(3./u2*u0ij(fillMa1, fillMa2)-KronDelta)
        END DO
      END DO
      SMat(2,1)=SMat(1,2)
      SMat(3,1)=SMat(1,3)
      SMat(3,2)=SMat(2,3)
      CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
    ELSE
      !! Exact Solution
      DO fillMa1 =1, 3
        DO fillMa2 =fillMa1, 3
          IF (fillMa1.EQ.fillMa2) THEN
            KronDelta = 1.0
          ELSE
            KronDelta = 0.0
          END IF
          A(fillMa1, fillMa2) = KronDelta - (1.-Prandtl)/Prandtl*(3.*u0ij(fillMa1, fillMa2)/u2 - KronDelta)
        END DO
      END DO
      IF (ESBGKModel.EQ.2) THEN
        CALL DSYEV('V','U',3,A,3,W,Work,100,INFO)
        SMat = 0.0
        IF (W(1).LT.0.0) THEN
          W(1) = 0.0
          IF (W(2).LT.0) W(2) = 0.0
        END IF
        IF (W(3).LT.0) THEN
          W(3) = 0.0
          DO fillMa1 =1, 3
            DO fillMa2 =fillMa1, 3
              IF (fillMa1.EQ.fillMa2) THEN
                KronDelta = 1.0
              ELSE
                KronDelta = 0.0
              END IF
              SMat(fillMa1, fillMa2)= KronDelta - (1.-Prandtl)/(2.*Prandtl) &
                *(3./u2*u0ij(fillMa1, fillMa2)-KronDelta)
            END DO
          END DO
          SMat(2,1)=SMat(1,2)
          SMat(3,1)=SMat(1,3)
          SMat(3,2)=SMat(2,3)
        ELSE
          SMat(1,1) = SQRT(W(1))
          SMat(2,2) = SQRT(W(2))
          SMat(3,3) = SQRT(W(3))
          SMat = MATMUL(A, SMat)
          SMat = MATMUL(SMat, TRANSPOSE(A))
        END IF
        CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
      ELSE IF (ESBGKModel.EQ.3) THEN
        A(2,1)=A(1,2)
        A(3,1)=A(1,3)
        A(3,2)=A(2,3)
        CALL MetropolisES(nRelax, iRanPart, A)
      END IF
    END IF
  CASE (2)  ! Shakov
!    CALL MetropolisShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
    CALL ARShakhov(nRelax, iRanPart, u2/3., u2i, Prandtl)
  CASE (3)  ! Standard BGK (Maxwell target distribution)
    CALL BGK_BuildTransGaussNums(nRelax, iRanPart)
  END SELECT
  DO iLoop = 1, nRelax
    iSpec = PartSpecies(iPartIndx_NodeRelax(iLoop))
    IF ((BGKCollModel.EQ.1).AND.(ESBGKModel.NE.3)) THEN
      tempVelo(1:3) = SQRT(BoltzmannConst*CellTemp/Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
      DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop)) = vBulkAll(1:3) + MATMUL(SMat,tempVelo)
    ELSE
      DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop)) = vBulkAll(1:3) &
        + SQRT(BoltzmannConst*CellTemp/Species(iSpec)%MassIC)*iRanPart(1:3,iLoop)
    END IF
    partWeight = GetParticleWeight(iPartIndx_NodeRelax(iLoop))
    vBulk(1:3) = vBulk(1:3) + DSMC_RHS(1:3,iPartIndx_NodeRelax(iLoop))*Species(iSpec)%MassIC*partWeight
  END DO
END IF ! nRelax.GT.0

END SUBROUTINE SampleFromTargetDistr


SUBROUTINE EnergyConsVib(nPart, nVibRelax, nVibRelaxSpec, iPartIndx_NodeRelaxVib, NewEnVib, OldEn, Xi_VibSpec, VibEnergyDOF, TEqui)
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
INTEGER, INTENT(IN)           :: nPart, nVibRelax, iPartIndx_NodeRelaxVib(:), nVibRelaxSpec(nSpecies)
REAL, INTENT(IN)              :: NewEnVib, VibEnergyDOF(:,:), Xi_VibSpec(nSpecies), TEqui
REAL, INTENT(INOUT)           :: OldEn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iLoop, iDOF, iSpec, iQuant, iQuaMax, iPolyatMole
REAL                          :: alpha, partWeight, betaV, iRan, MaxColQua, Xi_VibTotal
!===================================================================================================================================
IF(BGKDoVibRelaxation) THEN
  IF ((NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)) THEN
    IF (BGKUseQuantVibEn) THEN
      alpha = OldEn/NewEnVib*(Xi_VibSpec(1)*nVibRelax/(3.*(nPart-1.)+Xi_VibSpec(1)*nVibRelax))
      DO iLoop = 1, nVibRelax
        partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
        iSpec = PartSpecies(iPartIndx_NodeRelaxVib(iLoop))
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN      
          PartStateIntEn(1,iPartIndx_NodeRelaxVib(iLoop)) = 0.0
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            betaV = alpha*VibEnergyDOF(iLoop,iDOF)/(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst)
            CALL RANDOM_NUMBER(iRan)
            iQuant = INT(betaV+iRan)
            IF(iQuant.GT.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)) iQuant=PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)
            IF ((OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight).LT.0.0) THEN
              MaxColQua = OldEn/(BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*partWeight)
              IF (INT(MaxColQua).EQ.0) THEN
                iQuant = 0
              ELSE
                CALL RANDOM_NUMBER(iRan)
                iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                iQuaMax = MIN(INT(MaxColQua)+1, PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF))
                DO WHILE (iQuant.GE.iQuaMax)
                  CALL RANDOM_NUMBER(iRan)
                  iQuant = INT(-LOG(iRan)*TEqui/PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF))
                END DO
              END IF
            END IF
            PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop))  = PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) &
               + iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            VibQuantsPar(iPartIndx_NodeRelaxVib(iLoop))%Quants(iDOF) = iQuant
            OldEn = OldEn - iQuant*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight
          END DO
          PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop))  = PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) &
               + SpecDSMC(iSpec)%EZeroPoint
        ELSE  ! Diatomic molecules
          betaV = alpha*PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop))/(SpecDSMC(iSpec)%CharaTVib*BoltzmannConst)
          CALL RANDOM_NUMBER(iRan)
          iQuant = INT(betaV+iRan)
          IF (iQuant.GT.SpecDSMC(iSpec)%MaxVibQuant) iQuant = SpecDSMC(iSpec)%MaxVibQuant
          PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop))  = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
          IF ((OldEn - (PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) - SpecDSMC(iSpec)%EZeroPoint)*partWeight).LT.0.0) THEN
            MaxColQua = OldEn/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib*partWeight)
            IF (INT(MaxColQua).EQ.0) THEN
              iQuant = 0
            ELSE
              CALL RANDOM_NUMBER(iRan)
              iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(iSpec)%CharaTVib)
              iQuaMax = MIN(INT(MaxColQua)+1, SpecDSMC(iSpec)%MaxVibQuant)
              DO WHILE (iQuant.GE.iQuaMax)
                CALL RANDOM_NUMBER(iRan)
                iQuant = INT(-LOG(iRan)*TEqui/SpecDSMC(iSpec)%CharaTVib)
              END DO
            END IF
            PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop))  = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
          END IF
          OldEn = OldEn - (PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
        END IF ! SpecDSMC(1)%PolyatomicMol
      END DO
    ELSE ! Continuous treatment of vibrational energy
      Xi_VibTotal = 0.0
      DO iSpec = 1, nSpecies
        Xi_VibTotal = Xi_VibTotal + Xi_VibSpec(iSpec)*nVibRelaxSpec(iSpec)
      END DO
      alpha = OldEn/NewEnVib*(Xi_VibTotal/(3.*(nPart-1.)+Xi_VibTotal))
      DO iLoop = 1, nVibRelax
        iSpec = PartSpecies(iPartIndx_NodeRelaxVib(iLoop))
        partWeight = GetParticleWeight(iPartIndx_NodeRelaxVib(iLoop))
        PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) = alpha*PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) &
          + SpecDSMC(iSpec)%EZeroPoint
        OldEn = OldEn - (PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
      END DO
    END IF ! BGKUseQuantVibEn
  ELSE IF (nVibRelax.GT.0) THEN ! Relaxation towards the vibrational ground-state (new state is simply the zero-point energy)
    DO iLoop = 1, nVibRelax
      iSpec = PartSpecies(iPartIndx_NodeRelaxVib(iLoop))
      PartStateIntEn( 1,iPartIndx_NodeRelaxVib(iLoop)) = SpecDSMC(iSpec)%EZeroPoint
    END DO
  END IF ! (NewEnVib.GT.0.0).AND.(nVibRelax.GT.0)
END IF ! BGKDoVibRelaxation

END SUBROUTINE EnergyConsVib

SUBROUTINE ARGrads13(nPart, iRanPart, Vtherm, HeatVec, PressTens)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Vtherm, PressTens(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vheat, V2, iRan, OldProb, Envelope, Envelope2, cMat, KronDelta
INTEGER                        :: iPart, fillMa1, fillMa2
!===================================================================================================================================
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope2 = MAX(ABS(PressTens(1,2)),ABS(PressTens(1,3)),ABS(PressTens(2,3)))/Vtherm
Envelope =  1.+3.*MAX(Envelope, Envelope2)

DO iPart = 1, nPart
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  cMat = 0.0
  DO fillMa1 =1, 3
    DO fillMa2 =1, 3
      IF (fillMa1.EQ.fillMa2) THEN
        KronDelta = 1.0
      ELSE
        KronDelta = 0.0
      END IF
      cMat = cMat + iRanPart(fillMa1,iPart)*iRanPart(fillMa2,iPart)*(PressTens(fillMa1,fillMa2)-KronDelta*Vtherm)
    END DO
  END DO
!  cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
!  cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
!  cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + cMat/(2.*Vtherm) + VHeat/(Vtherm**(3./2.))*(V2/5.-1.))
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    cMat = 0.0
    DO fillMa1 =1, 3
      DO fillMa2 =1, 3
        IF (fillMa1.EQ.fillMa2) THEN
          KronDelta = 1.0
        ELSE
          KronDelta = 0.0
        END IF
        cMat = cMat + iRanPart(fillMa1,iPart)*iRanPart(fillMa2,iPart)*(PressTens(fillMa1,fillMa2)-KronDelta*Vtherm)
      END DO
    END DO
!    cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
!    cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
!    cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb =  (1. + cMat/(2.*Vtherm) + VHeat/(Vtherm**(3./2.))*(V2/5.-1.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARGrads13

SUBROUTINE ARChapEnsk(nPart, iRanPart, Vtherm, HeatVec, PressTens)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Vtherm, PressTens(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vheat, V2, iRan, OldProb, Envelope, Envelope2, cMat, cPress
INTEGER                        :: iPart
!===================================================================================================================================
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope2 = MAX(ABS(PressTens(1,2)),ABS(PressTens(1,3)),ABS(PressTens(2,3)))/Vtherm
Envelope =  1.+4.*MAX(Envelope, Envelope2)

DO iPart = 1, nPart
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  cMat = 0.0
  cPress = 0.0
  cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
  cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
  cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
  cPress=cPress + (PressTens(1,1)-Vtherm)*(iRanPart(1,iPart)*iRanPart(1,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
  cPress=cPress + (PressTens(2,2)-Vtherm)*(iRanPart(2,iPart)*iRanPart(2,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + cMat/Vtherm + cPress/(2.*Vtherm) + VHeat/(2.*Vtherm**(3./2.))*(V2/5.-1.))
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    cMat = 0.0
    cPress = 0.0
    cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
    cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
    cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
    cPress=cPress + (PressTens(1,1)-Vtherm)*(iRanPart(1,iPart)*iRanPart(1,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
    cPress=cPress + (PressTens(2,2)-Vtherm)*(iRanPart(2,iPart)*iRanPart(2,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb =  (1. + cMat/Vtherm + cPress/(2.*Vtherm) + VHeat/(2.*Vtherm**(3./2.))*(V2/5.-1.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARChapEnsk

SUBROUTINE MetropolisES(nPart, iRanPart, A)
!===================================================================================================================================
!> description
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Basis ,ONLY: INV33
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: A(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: iRanPartTemp(3), V2, iRan, NewProb, OldProb, NormProb
INTEGER                        :: iLoop, iPart, iRun
LOGICAL                        :: Changed
REAL                           :: AC(3), AInvers(3,3), detA
!===================================================================================================================================
iRanPart(1,1) = rnor()
iRanPart(2,1) = rnor()
iRanPart(3,1) = rnor()
CALL INV33(A,AInvers, detA)
AC(1:3) = MATMUL(AInvers, iRanPart(1:3,1))
V2 = iRanPart(1,1)*AC(1) + iRanPart(2,1)*AC(2) + iRanPart(3,1)*AC(3)
OldProb = EXP(-0.5*V2)
!Burn in
DO iLoop = 1, 35 !50
  iRanPartTemp(1) = rnor()
  iRanPartTemp(2) = rnor()
  iRanPartTemp(3) = rnor()
  AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
  V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
  NewProb = EXP(-0.5*V2)
  NormProb = MIN(1.,NewProb/OldProb)
  CALL RANDOM_NUMBER(iRan)
  IF (NormProb.GT.iRan) THEN
    iRanPart(1:3,1) = iRanPartTemp(1:3)
    OldProb = NewProb
  END IF
END DO
! All the others
DO iPart = 2, nPart
  iRanPart(1,iPart) = iRanPart(1,iPart-1)
  iRanPart(2,iPart) = iRanPart(2,iPart-1)
  iRanPart(3,iPart) = iRanPart(3,iPart-1)
  iRun = 0
  Changed = .FALSE.
  DO WHILE ((iRun.LT.10).OR.(.NOT.Changed))
    iRun = iRun + 1
    iRanPartTemp(1) = rnor()
    iRanPartTemp(2) = rnor()
    iRanPartTemp(3) = rnor()
    AC(1:3) = MATMUL(AInvers, iRanPartTemp(1:3))
    V2 = iRanPartTemp(1)*AC(1) + iRanPartTemp(2)*AC(2) + iRanPartTemp(3)*AC(3)
    NewProb = EXP(-0.5*V2)
    NormProb = MIN(1.,NewProb/OldProb)
    CALL RANDOM_NUMBER(iRan)
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
!> description
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
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope =  1.+4.*Envelope

DO iPart = 1, nPart
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb = (1. + (1.-Prandtl)*VHeat/(5.*Vtherm**(3./2.))*(V2/2.-5./2.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARShakhov

SUBROUTINE BGK_BuildTransGaussNums(nPart, iRanPart)
!===================================================================================================================================
!> description
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
DO iLoop = 1, nPart
  iRanPart(1,iLoop) = rnor()
  iRanPart(2,iLoop) = rnor()
  iRanPart(3,iLoop) = rnor()
END DO

END SUBROUTINE BGK_BuildTransGaussNums

SUBROUTINE CalcTEqui(nPart, CellTemp, TRot, TVib, Xi_Vib, Xi_Vib_old, RotExp, VibExp,  &
      TEqui, rotrelaxfreq, vibrelaxfreq, dtCell, DoVibRelaxIn)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY: SpecDSMC
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq, dtCell
INTEGER, INTENT(IN)             :: nPart
LOGICAL, OPTIONAL, INTENT(IN)   :: DoVibRelaxIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_vib, TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaR, betaV, RotFrac, VibFrac, TEqui_Old2
REAL                            :: eps_prec=1.0E-0
REAL                            :: correctFac, correctFacRot, maxexp   !, Xi_rel
LOGICAL                         :: DoVibRelax
!===================================================================================================================================
IF (PRESENT(DoVibRelaxIn)) THEN
  DoVibRelax = DoVibRelaxIn
ELSE
  DoVibRelax = BGKDoVibRelaxation
END IF
maxexp = LOG(HUGE(maxexp))
!  Xi_rel = 2.*(2. - CollInf%omega(1,1))
!  correctFac = 1. + (2.*SpecDSMC(1)%CharaTVib / (CellTemp*(EXP(SpecDSMC(1)%CharaTVib / CellTemp)-1.)))**(2.) &
!        * EXP(SpecDSMC(1)%CharaTVib /CellTemp) / (2.*Xi_rel)
!  correctFacRot = 1. + 2./Xi_rel

correctFac = 1.
correctFacRot = 1.
RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
RotFrac = nPart*(1.-RotExp)
IF(DoVibRelax) THEN
  VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
  VibFrac = nPart*(1.-VibExp)
ELSE
  VibExp = 0.0
  VibFrac = 0.0
  Xi_vib = 0.0
END IF
TEqui_Old = 0.0
TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib_old*VibFrac)
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  IF (ABS(TRot-TEqui).LT.1E-3) THEN
    RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
  ELSE
    betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dtCell/correctFacRot
    IF (-betaR.GT.0.0) THEN
      RotExp = 0.
    ELSE IF (betaR.GT.maxexp) THEN
      RotExp = 0.
    ELSE
      RotExp = exp(-betaR)
    END IF
  END IF
  RotFrac = nPart*(1.-RotExp)
  IF(DoVibRelax) THEN
    IF (ABS(TVib-TEqui).LT.1E-3) THEN
      VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
    ELSE
      betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dtCell/correctFac
      IF (-betaV.GT.0.0) THEN
        VibExp = 0.
      ELSE IF (betaV.GT.maxexp) THEN
        VibExp = 0.
      ELSE
        VibExp = exp(-betaV)
      END IF
    END IF
    IF ((SpecDSMC(1)%CharaTVib/TEqui).GT.maxexp) THEN
      Xi_Vib = 0.0
    ELSE
      Xi_vib = 2.*SpecDSMC(1)%CharaTVib/TEqui/(EXP(SpecDSMC(1)%CharaTVib/TEqui)-1.)
    END IF
    VibFrac = nPart*(1.-VibExp)
  END IF
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui
  TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib*VibFrac)
  IF(DoVibRelax) THEN
    DO WHILE( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      TEqui =(TEqui + TEqui_Old2)*0.5
      IF ((SpecDSMC(1)%CharaTVib/TEqui).GT.maxexp) THEN
        Xi_Vib = 0.0
      ELSE
        Xi_vib = 2.*SpecDSMC(1)%CharaTVib/TEqui/(EXP(SpecDSMC(1)%CharaTVib/TEqui)-1.)
      END IF
      TEqui_Old2 = TEqui
      TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib) / (3.*(nPart-1.)+2.*RotFrac+Xi_vib*VibFrac)
    END DO
  END IF
END DO

END SUBROUTINE CalcTEqui

SUBROUTINE CalcTEquiMulti(nPart, nSpec, CellTemp, TRotSpec, TVibSpec, Xi_VibSpec, Xi_Vib_oldSpec, RotExpSpec, VibExpSpec,  &
      TEqui, rotrelaxfreqSpec, vibrelaxfreqSpec, dtCell, DoVibRelaxIn)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY: SpecDSMC
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
USE MOD_Particle_Vars,          ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TRotSpec(nSpecies), TVibSpec(nSpecies), Xi_Vib_oldSpec(nSpecies)
REAL, INTENT(IN)                :: rotrelaxfreqSpec(nSpecies), vibrelaxfreqSpec(nSpecies), dtCell
INTEGER, INTENT(IN)             :: nPart, nSpec(nSpecies)
LOGICAL, OPTIONAL, INTENT(IN)   :: DoVibRelaxIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_VibSpec(nSpecies), TEqui, RotExpSpec(nSpecies), VibExpSpec(nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaR, betaV, RotFracSpec(nSpecies), VibFracSpec(nSpecies), TEqui_Old2
REAL                            :: eps_prec=1.0E-0
REAL                            :: correctFac, correctFacRot, maxexp, TEquiNumDof   !, Xi_rel, 
LOGICAL                         :: DoVibRelax
INTEGER                         :: iSpec
!===================================================================================================================================
IF (PRESENT(DoVibRelaxIn)) THEN
  DoVibRelax = DoVibRelaxIn
ELSE
  DoVibRelax = BGKDoVibRelaxation
END IF
maxexp = LOG(HUGE(maxexp))
!  Xi_rel = 2.*(2. - CollInf%omega(1,1))
!  correctFac = 1. + (2.*SpecDSMC(1)%CharaTVib / (CellTemp*(EXP(SpecDSMC(1)%CharaTVib / CellTemp)-1.)))**(2.) &
!        * EXP(SpecDSMC(1)%CharaTVib /CellTemp) / (2.*Xi_rel)
!  correctFacRot = 1. + 2./Xi_rel

correctFac = 1.
correctFacRot = 1.
RotFracSpec = 0.0
VibFracSpec = 0.0
DO iSpec=1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    RotExpSpec(iSpec) = exp(-rotrelaxfreqSpec(iSpec)*dtCell/correctFacRot)
    RotFracSpec(iSpec) = nSpec(iSpec)*(1.-RotExpSpec(iSpec))
    IF(DoVibRelax) THEN
      VibExpSpec(iSpec) = exp(-vibrelaxfreqSpec(iSpec)*dtCell/correctFac)
      VibFracSpec(iSpec) = nSpec(iSpec)*(1.-VibExpSpec(iSpec))
    ELSE
      VibExpSpec(iSpec) = 0.0
      VibFracSpec(iSpec) = 0.0
      Xi_VibSpec(iSpec) = 0.0
    END IF
  END IF
END DO
TEqui_Old = 0.0
TEqui = 3.*(nPart-1.)*CellTemp
TEquiNumDof = 3.*(nPart-1.)
DO iSpec=1, nSpecies
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    TEqui = TEqui + 2.*RotFracSpec(iSpec)*TRotSpec(iSpec)+Xi_Vib_oldSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
    TEquiNumDof = TEquiNumDof + 2.*RotFracSpec(iSpec) + Xi_Vib_oldSpec(iSpec)*VibFracSpec(iSpec)
  END IF
END DO
TEqui = TEqui / TEquiNumDof
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  DO iSpec = 1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      IF (ABS(TRotSpec(iSpec)-TEqui).LT.1E-3) THEN
        RotExpSpec(iSpec) = exp(-rotrelaxfreqSpec(iSpec)*dtCell/correctFacRot)
      ELSE
        betaR = ((TRotSpec(iSpec)-CellTemp)/(TRotSpec(iSpec)-TEqui))*rotrelaxfreqSpec(iSpec)*dtCell/correctFacRot
        IF (-betaR.GT.0.0) THEN
          RotExpSpec(iSpec) = 0.
        ELSE IF (betaR.GT.maxexp) THEN
          RotExpSpec(iSpec) = 0.
        ELSE
          RotExpSpec(iSpec) = exp(-betaR)
        END IF
      END IF
      RotFracSpec(iSpec) = nSpec(iSpec)*(1.-RotExpSpec(iSpec))
      IF(DoVibRelax) THEN
        IF (ABS(TVibSpec(iSpec)-TEqui).LT.1E-3) THEN
          VibExpSpec(iSpec) = exp(-vibrelaxfreqSpec(iSpec)*dtCell/correctFac)
        ELSE
          betaV = ((TVibSpec(iSpec)-CellTemp)/(TVibSpec(iSpec)-TEqui))*vibrelaxfreqSpec(iSpec)*dtCell/correctFac
          IF (-betaV.GT.0.0) THEN
            VibExpSpec(iSpec) = 0.
          ELSE IF (betaV.GT.maxexp) THEN
            VibExpSpec(iSpec) = 0.
          ELSE
            VibExpSpec(iSpec) = exp(-betaV)
          END IF
        END IF
        IF ((SpecDSMC(iSpec)%CharaTVib/TEqui).GT.maxexp) THEN
          Xi_VibSpec(iSpec) = 0.0
        ELSE
          Xi_VibSpec(iSpec) = 2.*SpecDSMC(iSpec)%CharaTVib/TEqui/(EXP(SpecDSMC(iSpec)%CharaTVib/TEqui)-1.)
        END IF
        VibFracSpec(iSpec) = nSpec(iSpec)*(1.-VibExpSpec(iSpec))
      END IF
    END IF
  END DO
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui

  TEqui = 3.*(nPart-1.)*CellTemp
  TEquiNumDof = 3.*(nPart-1.)
  DO iSpec=1, nSpecies
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      TEqui = TEqui + 2.*RotFracSpec(iSpec)*TRotSpec(iSpec)+Xi_Vib_oldSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
      TEquiNumDof = TEquiNumDof + 2.*RotFracSpec(iSpec) + Xi_VibSpec(iSpec)*VibFracSpec(iSpec)
    END IF
  END DO
  TEqui = TEqui / TEquiNumDof
  IF(DoVibRelax) THEN
    DO WHILE( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      TEqui =(TEqui + TEqui_Old2)*0.5
      DO iSpec=1, nSpecies
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          IF ((SpecDSMC(iSpec)%CharaTVib/TEqui).GT.maxexp) THEN
            Xi_VibSpec(iSpec) = 0.0
          ELSE
            Xi_VibSpec(iSpec) = 2.*SpecDSMC(iSpec)%CharaTVib/TEqui/(EXP(SpecDSMC(iSpec)%CharaTVib/TEqui)-1.)
          END IF
        END IF
      END DO
      TEqui_Old2 = TEqui
      TEqui = 3.*(nPart-1.)*CellTemp
      TEquiNumDof = 3.*(nPart-1.)
      DO iSpec=1, nSpecies
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          TEqui = TEqui + 2.*RotFracSpec(iSpec)*TRotSpec(iSpec)+Xi_Vib_oldSpec(iSpec)*VibFracSpec(iSpec)*TVibSpec(iSpec)
          TEquiNumDof = TEquiNumDof + 2.*RotFracSpec(iSpec) + Xi_VibSpec(iSpec)*VibFracSpec(iSpec)
        END IF
      END DO
      TEqui = TEqui / TEquiNumDof
    END DO
  END IF
END DO
END SUBROUTINE CalcTEquiMulti



SUBROUTINE CalcTEquiPoly(nPart, CellTemp, TRot, TVib, Xi_Vib_DOF, Xi_Vib_old, RotExp, VibExp, TEqui, rotrelaxfreq, vibrelaxfreq, &
      dtCell, DoVibRelaxIn)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY: SpecDSMC, PolyatomMolDSMC
USE MOD_BGK_Vars,               ONLY: BGKDoVibRelaxation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp, TRot, TVib, Xi_Vib_old, rotrelaxfreq, vibrelaxfreq
INTEGER, INTENT(IN)             :: nPart
REAL, INTENT(IN)                :: dtCell
LOGICAL, OPTIONAL, INTENT(IN)   :: DoVibRelaxIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Xi_vib_DOF(:), TEqui, RotExp, VibExp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                            :: TEqui_Old, betaR, betaV, RotFrac, VibFrac, Xi_Rot, TEqui_Old2
REAL                            :: eps_prec=1.0
REAL                            :: correctFac, correctFacRot, maxexp
INTEGER                         :: iDOF, iPolyatMole
LOGICAL                         :: DoVibRelax
!===================================================================================================================================
IF (PRESENT(DoVibRelaxIn)) THEN
  DoVibRelax = DoVibRelaxIn
ELSE
  DoVibRelax = BGKDoVibRelaxation
END IF

maxexp = LOG(HUGE(maxexp))
Xi_Rot =   SpecDSMC(1)%Xi_Rot
iPolyatMole = SpecDSMC(1)%SpecToPolyArray
!  Xi_rel = 2.*(2. - CollInf%omega(1,1))
!  correctFac = 0.0
!  DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!    correctFac = correctFac &
!        + (2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (CellTemp           &
!        *(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / CellTemp)-1.)))**(2.)  &
!        * EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / CellTemp) / 2.
!  END DO
!  correctFac = 1. + correctFac/Xi_rel
!  correctFacRot = 1. + Xi_Rot/Xi_rel

correctFac = 1.
correctFacRot = 1.
RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
RotFrac = nPart*(1.-RotExp)
IF(DoVibRelax) THEN
  VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
  VibFrac = nPart*(1.-VibExp)
ELSE
  VibExp = 0.0
  VibFrac = 0.0
  Xi_vib_DOF = 0.0
END IF
TEqui_Old = 0.0
TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)/(3.*(nPart-1.)+2.*RotFrac+Xi_Vib_old*VibFrac)
DO WHILE ( ABS( TEqui - TEqui_Old ) .GT. eps_prec )
  IF (ABS(TRot-TEqui).LT.1E-3) THEN
    RotExp = exp(-rotrelaxfreq*dtCell/correctFacRot)
  ELSE
    betaR = ((TRot-CellTemp)/(TRot-TEqui))*rotrelaxfreq*dtCell/correctFacRot
    IF (-betaR.GT.0.0) THEN
      RotExp = 0.
    ELSE IF (betaR.GT.maxexp) THEN
      RotExp = 0.
    ELSE
      RotExp = exp(-betaR)
    END IF
  END IF
  RotFrac = nPart*(1.-RotExp)
  IF(DoVibRelax) THEN
    IF (ABS(TVib-TEqui).LT.1E-3) THEN
      VibExp = exp(-vibrelaxfreq*dtCell/correctFac)
    ELSE
      betaV = ((TVib-CellTemp)/(TVib-TEqui))*vibrelaxfreq*dtCell/correctFac
      IF (-betaV.GT.0.0) THEN
        VibExp = 0.
      ELSE IF (betaV.GT.maxexp) THEN
        VibExp = 0.
      ELSE
        VibExp = exp(-betaV)
      END IF
    END IF
    DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
      IF ((PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui).LT.maxexp) THEN
        Xi_vib_DOF(iDOF) = 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui &
                                    /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui)-1.)
      ELSE
        Xi_vib_DOF(iDOF) = 0.0
      END IF
    END DO
    VibFrac = nPart*(1.-VibExp)
  END IF
  TEqui_Old = TEqui
  TEqui_Old2 = TEqui
  TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)  &
          / (3.*(nPart-1.)+2.*RotFrac+SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))*VibFrac)
  IF(DoVibRelax) THEN
    DO WHILE( ABS( TEqui - TEqui_Old2 ) .GT. eps_prec )
      TEqui =(TEqui + TEqui_Old2)*0.5
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF ((PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui).LT.maxexp) THEN
          Xi_vib_DOF(iDOF) = 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui &
                                      /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/TEqui)-1.)
        ELSE
          Xi_vib_DOF(iDOF) = 0.0
        END IF
      END DO
      TEqui_Old2 = TEqui
      TEqui = (3.*(nPart-1.)*CellTemp+2.*RotFrac*TRot+Xi_Vib_old*VibFrac*TVib)  &
          / (3.*(nPart-1.)+2.*RotFrac+SUM(Xi_vib_DOF(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))*VibFrac)
    END DO
  END IF
END DO

END SUBROUTINE CalcTEquiPoly

SUBROUTINE CalcViscosityThermalCondColIntVHS(CellTemp, Xi, dens, Visc, ThermalCond)
!===================================================================================================================================
!> Determination of the mixture viscosity and thermal conductivity using collision integrals (derived for the Variable Hard
!> Sphere model). Solving an equation system depending on the number of species.
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : CollInf
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species, nSpecies
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp(nSpecies+1), Xi(nSpecies), dens
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)               :: Visc,ThermalCond
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL        :: Sigma_11, Sigma_22, B_12(nSpecies,nSpecies), A_12(nSpecies,nSpecies), InteractDiam, cv, DiffCoef(nSpecies, nSpecies)
REAL        :: Mass, ViscSpec(nSpecies), ThermalCondSpec(nSpecies), TVHS, omegaVHS, E_12, CellTemptmp
REAL        :: ViscMat(nSpecies, nSpecies), RHSSolve(nSpecies), m0, pressure
INTEGER     :: iSpec, jSpec, kSpec, IPIV(nSpecies), info_dgesv
!===================================================================================================================================
ViscSpec = 0.; ThermalCondSpec = 0.; DiffCoef =0.; A_12 = 0.; B_12 = 0.
DO iSpec = 1, nSpecies
  DO jSpec = iSpec, nSpecies
    InteractDiam = CollInf%dref(iSpec,jSpec)
    Mass = Species(iSpec)%MassIC*Species(jSpec)%MassIC/(Species(iSpec)%MassIC + Species(jSpec)%MassIC)
    TVHS = CollInf%Tref(iSpec,jSpec)
    omegaVHS = CollInf%omega(iSpec,jSpec)
    IF (iSpec.EQ.jSpec) THEN
      CellTemptmp = CellTemp(iSpec)
    ELSE
      CellTemptmp = CellTemp(nSpecies+1)
    END IF
    Sigma_22 = CalcSigma_22VHS(CellTemptmp,InteractDiam,Mass,TVHS, omegaVHS)
    IF (iSpec.EQ.jSpec) THEN
      cv= 3./2.*BoltzmannConst/(2.*Mass)
      ViscSpec(iSpec) = (5./8.)*(BoltzmannConst*CellTemp(iSpec))/Sigma_22
      ThermalCondSpec(iSpec) = (25./16.)*(cv*BoltzmannConst*CellTemp(iSpec))/Sigma_22
      !ThermalCondSpec(iSpec) = (15./4.)*BoltzmannConst/(2.*Mass)*ViscSpec(iSpec) 
    ELSE     
      CALL CalcSigma_11VHS(CellTemp(nSpecies+1),InteractDiam,Mass,TVHS, omegaVHS, Sigma_11)
      B_12(iSpec,jSpec) = (5.*GAMMA(4.-omegaVHS)-GAMMA(5.-omegaVHS))/(5.*GAMMA(3.-omegaVHS))
      B_12(jSpec,iSpec) = B_12(iSpec,jSpec)
      A_12(iSpec,jSpec) = Sigma_22 / (5.*Sigma_11)
      A_12(jSpec,iSpec) = A_12(iSpec,jSpec)
      E_12 = BoltzmannConst*CellTemp(nSpecies+1)/(8.*Species(iSpec)%MassIC*Species(jSpec)%MassIC & 
          /(Species(iSpec)%MassIC+Species(jSpec)%MassIC)**2.*Sigma_11)
      DiffCoef(iSpec,jSpec) = 3.*E_12/(2.*(Species(iSpec)%MassIC+Species(jSpec)%MassIC)*dens)
      DiffCoef(jSpec,iSpec) = DiffCoef(iSpec,jSpec)
    END IF
  END DO
END DO

ViscMat = 0.0
DO iSpec = 1, nSpecies
  DO jSpec = 1, nSpecies
    IF (iSpec.EQ.jSpec) THEN
      ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + Xi(iSpec)/ViscSpec(iSpec)
      DO kSpec = 1, nSpecies
        IF(kSpec.EQ.iSpec) CYCLE
        ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + 3.*Xi(kSpec) / ((Species(iSpec)%MassIC*dens &
        + Species(kSpec)%MassIC*dens)*DiffCoef(iSpec,kSpec))*(2./3.+Species(kSpec)%MassIC/Species(iSpec)%MassIC*A_12(iSpec,kSpec))
      END DO
    ELSE
      ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) -  Xi(iSpec)*3. / ((Species(iSpec)%MassIC*dens &
        + Species(jSpec)%MassIC*dens)*DiffCoef(iSpec,jSpec))*(2./3.-A_12(iSpec,jSpec))
    END IF
  END DO
  RHSSolve(iSpec) = Xi(iSpec)
END DO
CALL DGESV(nSpecies, 1, ViscMat, nSpecies, IPIV, RHSSolve, nSpecies, info_dgesv)
Visc = SUM(RHSSolve)

pressure = BoltzmannConst*dens*CellTemp(nSpecies+1)
ViscMat = 0.0
DO iSpec = 1, nSpecies
  DO jSpec = 1, nSpecies
    IF (iSpec.EQ.jSpec) THEN
      ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + Xi(iSpec)/ThermalCondSpec(iSpec)
      DO kSpec = 1, nSpecies
        IF(kSpec.EQ.iSpec) CYCLE
        m0 = Species(iSpec)%MassIC+Species(kSpec)%MassIC
        ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) + CellTemp(nSpecies+1)*Xi(kSpec)/(5.*pressure*DiffCoef(iSpec,kSpec)) &
          * (6.*Species(iSpec)%MassIC**2./m0**2.+(5.-4.*B_12(iSpec,kSpec))*Species(kSpec)%MassIC**2./m0**2. &
          + 8.*Species(iSpec)%MassIC*Species(kSpec)%MassIC/m0**2.*A_12(iSpec, kSpec))
      END DO
    ELSE
      m0 = Species(iSpec)%MassIC+Species(jSpec)%MassIC
      ViscMat(iSpec, jSpec) = ViscMat(iSpec, jSpec) - Xi(iSpec)*CellTemp(nSpecies+1) &
        *(Species(iSpec)%MassIC*Species(jSpec)%MassIC/m0**2.)/(5.*pressure*DiffCoef(iSpec,jSpec)) &
        *(11.-4.*B_12(iSpec,jSpec)-8.*A_12(iSpec,jSpec))
    END IF
  END DO
  RHSSolve(iSpec) = Xi(iSpec)
END DO
CALL DGESV(nSpecies, 1, ViscMat, nSpecies, IPIV, RHSSolve, nSpecies, info_dgesv)
ThermalCond = SUM(RHSSolve)

END SUBROUTINE CalcViscosityThermalCondColIntVHS


SUBROUTINE CalcSigma_11VHS(CellTemp,Dref,Mass,Tref, omegaVHS, Sigma_11)
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars          ,ONLY: Pi, BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp,Dref,Mass,Tref, omegaVHS
REAL, INTENT(OUT)               :: Sigma_11
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL          :: Prefactor
!===================================================================================================================================
  Prefactor = Pi/2.*Dref*Dref*SQRT(BoltzmannConst/(2.*Pi*Mass))*Tref**omegaVHS*GAMMA(3.-omegaVHS)/GAMMA(2.-omegaVHS)
  Sigma_11 = Prefactor*CellTemp**(0.5-omegaVHS)

END SUBROUTINE CalcSigma_11VHS

REAL FUNCTION CalcSigma_22VHS(CellTemp,Dref,Mass,Tref, omegaVHS)
!===================================================================================================================================
!> 
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars          ,ONLY: Pi, BoltzmannConst
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: CellTemp,Dref,Mass,Tref, omegaVHS
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Prefactor
!===================================================================================================================================
  Prefactor = Pi/3.*Dref*Dref*SQRT(BoltzmannConst/(2.*Pi*Mass))*Tref**omegaVHS*GAMMA(4.-omegaVHS)/GAMMA(2.-omegaVHS)
  CalcSigma_22VHS = Prefactor*CellTemp**(0.5-omegaVHS)

END FUNCTION CalcSigma_22VHS

END MODULE MOD_BGK_CollOperator
