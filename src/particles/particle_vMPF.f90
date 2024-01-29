!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer and Asim Mirza
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

MODULE MOD_vMPF
!===================================================================================================================================
! Module controlling particle number by merge and split routines
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
PUBLIC :: SplitAndMerge
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Main routine for split and merge particles
!> Loop over all elements:
!> 1.) build partindx list for cell
!> 2.) build partindx list for species
!> 3.) Call split or merge routine
!===================================================================================================================================
SUBROUTINE SplitAndMerge()
! MODULES
USE MOD_PARTICLE_Vars         ,ONLY: vMPFMergeThreshold, vMPFSplitThreshold, PEM, nSpecies, PartSpecies,PDM
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_part_tools            ,ONLY: UpdateNextFreePosition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, iLoop, iPart, nPartCell, iSpec
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:,:), nPart(:)
#if USE_LOADBALANCE
REAL                  :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
ALLOCATE(nPart(nSpecies))
DO iElem = 1, nElems
  nPart(:) = 0
  nPartCell = PEM%pNumber(iElem)
  ALLOCATE(iPartIndx_Node(nSpecies,nPartCell))
  DO iSpec = 1, nSpecies
    iPartIndx_Node(iSpec,1:nPartCell) = 0
  END DO
  iPart = PEM%pStart(iElem)

  ! 1.) build partindx list for cell
  DO iLoop = 1, nPartCell
    IF (.NOT.PDM%ParticleInside(iPart)) THEN
      iPart = PEM%pNext(iPart)
      CYCLE
    END IF
    iSpec = PartSpecies(iPart)
    nPart(iSpec) = nPart(iSpec) + 1
    iPartIndx_Node(iSpec,nPart(iSpec)) = iPart
    iPart = PEM%pNext(iPart)
  END DO

  DO iSpec = 1, nSpecies
    ! Skip default values
    IF((vMPFMergeThreshold(iSpec).EQ.0).AND.(vMPFSplitThreshold(iSpec).EQ.0)) CYCLE
    ! Skip when no particles are present
    IF(nPart(iSpec).EQ.0) CYCLE

    ! 2.) Call split or merge routine
    IF(nPart(iSpec).GT.vMPFMergeThreshold(iSpec).AND.(vMPFMergeThreshold(iSpec).NE.0)) THEN   ! Merge
      CALL MergeParticles(iPartIndx_Node(iSpec,1:nPart(iSpec)), nPart(iSpec), vMPFMergeThreshold(iSpec),iElem)
    ELSE IF(nPart(iSpec).LT.vMPFSplitThreshold(iSpec)) THEN                                   ! Split
      CALL SplitParticles(iPartIndx_Node(iSpec,1:nPart(iSpec)), nPart(iSpec), vMPFSplitThreshold(iSpec))
    END IF

  END DO
  DEALLOCATE(iPartIndx_Node)

#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO
CALL UpdateNextFreePosition()
DEALLOCATE(nPart)

END SUBROUTINE SplitAndMerge


!===================================================================================================================================
!> Routine for merge particles
!> 1.) Calc bulk velocity v_bulk (for momentum conservation)
!> 2.) Calc temperature, energy and degree of fredoms (for energy conservation)
!> 2.1) Calc energies (E_trans, E_elec, E_vib, E_rot)
!> 2.2) Calc temperature and degree of fredoms
!> 3.) Delete particles randomly (until nPartNew is reached)
!> 4.) Calc bulkvelocity v_bulk_new after deleting
!> 5.) Calc energy after deleting
!> 5.1) E_trans_new
!> 5.2) E_elec_new
!> 5.3) E_vib_new
!> 5.4) E_rot_new
!> 6.) Ensuring momentum and energy conservation
!> 6.1) E_elec
!> 6.2) E_vib
!> 6.3) E_rot
!> 6.4) E_trans
!===================================================================================================================================
SUBROUTINE MergeParticles(iPartIndx_Node_in, nPart, nPartNew, iElem)
! MODULES
USE MOD_Globals               ,ONLY: ISFINITE
USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PartMPF, PartSpecies, Species, CellEelec_vMPF, CellEvib_vMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, CollisMode, SpecDSMC, DSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_Particle_Analyze_Tools,ONLY: CalcTelec, CalcTVibPoly
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Globals               ,ONLY: LOG_RAN
#ifdef CODE_ANALYZE
USE MOD_Globals               ,ONLY: unit_stdout,myrank,abort
USE MOD_Particle_Vars         ,ONLY: Symmetry
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
INTEGER, INTENT(IN)           :: iPartIndx_Node_in(nPart)
INTEGER, INTENT(IN)           :: nPartNew, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: iRan
INTEGER               :: iPartIndx_Node(nPart), iLoop, nDelete, nTemp, iPart, iSpec, iQua, iPolyatMole, iDOF, iQuaCount
REAL                  :: partWeight, totalWeight2, totalWeight, alpha
REAL                  :: V_rel(3), vmag2, vBulk(3)
REAL                  :: vBulk_new(3)
REAL                  :: T_elec, T_vib, T_rot, DOF_elec, DOF_vib, DOF_rot
REAL                  :: E_trans, E_trans_new
REAL                  :: E_elec, E_elec_new
REAL                  :: E_vib, E_vib_new
REAL                  :: E_rot, E_rot_new
REAL                  :: Energy_Sum, E_elec_upper, E_elec_lower, betaV
REAL, ALLOCATABLE     :: DOF_vib_poly(:), EnergyTemp_vibPoly(:,:)
#ifdef CODE_ANALYZE
REAL,PARAMETER        :: RelMomTol    = 5e-9  ! Relative tolerance applied to conservation of momentum before/after reaction
REAL,PARAMETER        :: RelEneTol    = 1e-12 ! Relative tolerance applied to conservation of energy before/after reaction
REAL,PARAMETER        :: RelWeightTol = 5e-12 ! Relative tolerance applied to conservation of mass before/after reaction
REAL                  :: Energy_old, Momentum_old(3),Energy_new, Momentum_new(3),totalWeightNew
INTEGER               :: iMomDim, iMom
#endif /* CODE_ANALYZE */
REAL                  :: lostWeight
!===================================================================================================================================
vBulk = 0.0; vBulk_new = 0.0;  totalWeight = 0.0; totalWeight2 = 0.0
E_trans = 0.0; E_trans_new = 0.0
E_elec = 0.0; E_elec_new = 0.0; DOF_elec = 0.0
E_vib = 0.0; E_vib_new = 0.0; DOF_vib = 0.0
E_rot = 0.0; E_rot_new = 0.0; DOF_rot = 0.0

iPartIndx_Node = iPartIndx_Node_in

iSpec = PartSpecies(iPartIndx_Node(1))  ! in iPartIndx_Node all particles are from same species

#ifdef CODE_ANALYZE
! Consider the remainder of the vibrational and electronic energy from the previous time step, which was not distributed due to
! quantized energy levels
Energy_old = CellEvib_vMPF(iSpec, iElem) + CellEelec_vMPF(iSpec, iElem)
Momentum_old = 0.0; Momentum_new = 0.0
totalWeightNew = 0.
#endif /* CODE_ANALYZE */

! 1.) calc bulk velocity (for momentum conservation)
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  totalWeight = totalWeight + partWeight
  totalWeight2 = totalWeight2 + partWeight*partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(4:6,iPart) * partWeight

#ifdef CODE_ANALYZE
  ! Energy conservation
  Energy_old = Energy_old + 0.5 * Species(iSpec)%MassIC &
  * DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart)) * partWeight
  IF(CollisMode.GT.1) THEN
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      Energy_old = Energy_old + (PartStateIntEn(1,iPart) +  PartStateIntEn(2,iPart)) * partWeight
    END IF
    IF(DSMC%ElectronicModel.GT.0) Energy_old = Energy_old + PartStateIntEn(3,iPart)*partWeight
  END IF
  ! Momentum conservation
  Momentum_old(1:3) = Momentum_old(1:3) + Species(iSpec)%MassIC * PartState(4:6,iPart) * partWeight
#endif /* CODE_ANALYZE */
END DO

vBulk(1:3) = vBulk(1:3) / totalWeight

! 2.) Calc energy, temperature and degree of freedom (for energy conservation)
! 2.1) Calc energy
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  V_rel(1:3)=PartState(4:6,iPart)-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  E_trans = E_trans + 0.5 * vmag2 * partWeight * Species(iSpec)%MassIC
  IF(CollisMode.GT.1) THEN
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      ! Rotational and vibrational energy
      E_vib = E_vib + (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
      E_rot = E_rot + partWeight * PartStateIntEn(2,iPart)
    END IF
    ! Electronic energy
    IF(DSMC%ElectronicModel.GT.0.AND.SpecDSMC(iSpec)%InterID.NE.4) THEN
      E_elec = E_elec + partWeight * PartStateIntEn(3,iPart)
    END IF
  END IF
END DO
IF((E_vib + CellEvib_vMPF(iSpec, iElem)).GT.0.0) THEN
  E_vib = E_vib + CellEvib_vMPF(iSpec, iElem)
  CellEvib_vMPF(iSpec, iElem) = 0.0
END IF
IF((E_elec + CellEelec_vMPF(iSpec, iElem)).GT.0.0) THEN
  E_elec = E_elec + CellEelec_vMPF(iSpec, iElem)
  CellEelec_vMPF(iSpec, iElem) = 0.0
END IF


! 2.2) Calc temperature and degree of freedoms
IF(CollisMode.GT.1) THEN
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      ALLOCATE(DOF_vib_poly(PolyatomMolDSMC(iPolyatMole)%VibDOF))
      DOF_vib_poly = 0.0
      ALLOCATE(EnergyTemp_vibPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF,nPartNew))
      EnergyTemp_vibPoly = 0.0
      T_vib = CalcTVibPoly(E_vib/totalWeight+SpecDSMC(iSpec)%EZeroPoint, iSpec)
      IF (T_vib.GT.0.0) THEN
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          DOF_vib= DOF_vib + 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/T_vib &
                              /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/T_vib) - 1.)
          DOF_vib_poly(iDOF) = 2.*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/T_vib &
                              /(EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/T_vib) - 1.)
        END DO
      END IF
    ELSE
      T_vib=E_vib / (totalWeight*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
      IF (T_vib.GT.0.0) THEN
        T_vib= SpecDSMC(iSpec)%CharaTVib/LOG(1. + 1./(T_vib))
        DOF_vib = 2.* E_vib / (totalWeight*BoltzmannConst*T_vib)
      END IF
    END IF
    DOF_rot = SpecDSMC(iSpec)%Xi_Rot
    T_rot = 2.*E_rot/(DOF_rot*totalWeight*BoltzmannConst)
  END IF
  IF(DSMC%ElectronicModel.GT.0.AND.SpecDSMC(iSpec)%InterID.NE.4) THEN
    T_elec = CalcTelec(E_elec/totalWeight, iSpec)
    IF (T_elec.GT.0.0) DOF_elec = 2.*E_elec/(totalWeight*BoltzmannConst*T_elec)
  END IF
END IF

! 3.) delete particles randomly (until nPartNew is reached)
lostWeight = 0.
nTemp = nPart
nDelete = nPart - nPartNew
DO iLoop = 1, nDelete
  CALL RANDOM_NUMBER(iRan)
  iPart = INT(iRan*nTemp) + 1
  partWeight = GetParticleWeight(iPartIndx_Node(iPart))
  lostWeight = lostWeight + partWeight
  PDM%ParticleInside(iPartIndx_Node(iPart)) = .FALSE.
  iPartIndx_Node(iPart) = iPartIndx_Node(nTemp)
  nTemp = nTemp - 1
END DO

! 4.) calc bulk velocity after deleting and set new MPF
DO iLoop = 1, nPartNew
  iPart = iPartIndx_Node(iLoop)
  ! Set new particle weight by distributing the total weight onto the remaining particles, proportional to their current weight
  PartMPF(iPart) = totalWeight * PartMPF(iPart) / (totalWeight - lostWeight)
  partWeight = GetParticleWeight(iPart)
  vBulk_new(1:3) = vBulk_new(1:3) + PartState(4:6,iPart) * partWeight
END DO
vBulk_new(1:3) = vBulk_new(1:3) / totalWeight

! 5.) calc energy after deleting
DO iLoop = 1, nPartNew
  iPart = iPartIndx_Node(iLoop)
  partWeight = GetParticleWeight(iPart)
  V_rel(1:3)=PartState(4:6,iPart)-vBulk_new(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  E_trans_new = E_trans_new + 0.5 * vmag2 * partWeight * Species(iSpec)%MassIC
  IF(CollisMode.GT.1) THEN
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      ! Rotational and vibrational energy
      E_vib_new = E_vib_new + (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint) * partWeight
      E_rot_new = E_rot_new + partWeight * PartStateIntEn(2,iPart)
    END IF
    ! Electronic energy
    IF(DSMC%ElectronicModel.GT.0.AND.SpecDSMC(iSpec)%InterID.NE.4) THEN
      E_elec_new = E_elec_new + partWeight * PartStateIntEn(3,iPart)
    END IF
  END IF
END DO

! 6.) ensuring momentum and energy conservation
! 6.1) ensuring electronic excitation
IF(CollisMode.GT.1) THEN
  IF(DSMC%ElectronicModel.GT.0.AND.SpecDSMC(iSpec)%InterID.NE.4) THEN
    Energy_Sum = E_elec
    IF (E_elec.GT.0.0) THEN
      IF (E_elec_new.EQ.0.0) THEN
!        E_elec_new = 0.0
        DO iLoop = 1, nPartNew  ! temporal continuous energy distribution
          iPart = iPartIndx_Node(iLoop)
          CALL RANDOM_NUMBER(iRan)
          PartStateIntEn(3,iPart) = -LOG_RAN()*DOF_elec*0.5*T_elec*BoltzmannConst
!          PartStateIntEn(3,iPart) = -LOG_RAN() * (E_elec/totalWeight)
          partWeight = GetParticleWeight(iPart)
          E_elec_new = E_elec_new + partWeight * PartStateIntEn(3,iPart)
        END DO
      END IF
      alpha = E_elec/E_elec_new
      DO iLoop = 1, nPartNew
!        alpha = E_elec/E_elec_new
!        Test_E = PartStateIntEn(3,iPart)
        iPart = iPartIndx_Node(iLoop)
        PartStateIntEn(3,iPart) = alpha * PartStateIntEn(3,iPart)
        iQuaCount = 0
        DO iQua = 1, SpecDSMC(iSpec)%MaxElecQuant - 1
          iQuaCount = iQuaCount + 1
          IF(BoltzmannConst*SpecDSMC(iSpec)%ElectronicState(2,iQua).GT.PartStateIntEn(3,iPart)) EXIT ! iQuant increased if EXIT?
        END DO
        E_elec_upper = SpecDSMC(iSpec)%ElectronicState(2,iQuaCount) * BoltzmannConst
        E_elec_lower = SpecDSMC(iSpec)%ElectronicState(2,iQuaCount-1) * BoltzmannConst
        CALL RANDOM_NUMBER(iRan)
        IF((PartStateIntEn(3,iPart)-E_elec_lower)/(E_elec_upper-E_elec_lower).LE.iRan) THEN
          PartStateIntEn(3,iPart) = E_elec_lower
        ELSE
          PartStateIntEn(3,iPart) = E_elec_upper
        END IF
        partWeight = GetParticleWeight(iPart)
        Energy_Sum = Energy_Sum - PartStateIntEn(3,iPart) * partWeight
!        IF(PartStateIntEn(3,iPart).EQ.0.0) then
!          E_elec_new = E_elec_new - Test_E * partWeight
!        END IF
      END DO
    END IF
    CellEelec_vMPF(iSpec, iElem) = CellEelec_vMPF(iSpec, iElem) + Energy_Sum
    Energy_Sum = 0.0
  END IF
END IF

! 6.2) ensuring vibrational excitation
IF(CollisMode.GT.1) THEN
  IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
    Energy_Sum = Energy_Sum + E_vib
    IF (E_vib.GT.0.0) THEN
      IF (E_vib_new.EQ.0.0) THEN
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DO iLoop = 1, nPartNew  ! temporal continuous energy distribution
            iPart = iPartIndx_Node(iLoop)
            PartStateIntEn(1,iPart) = 0.0
            DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
              CALL RANDOM_NUMBER(iRan)
              EnergyTemp_vibPoly(iDOF,iLoop) = -LOG_RAN()*DOF_vib_poly(iDOF)*0.5*T_vib*BoltzmannConst
              PartStateIntEn(1,iPart) = PartStateIntEn(1,iPart) + EnergyTemp_vibPoly(iDOF,iLoop)
            END DO
            partWeight = GetParticleWeight(iPart)
            E_vib_new = E_vib_new + partWeight * PartStateIntEn(1,iPart)
          END DO
        ELSE
          DO iLoop = 1, nPartNew  ! temporal continuous energy distribution
            iPart = iPartIndx_Node(iLoop)
            CALL RANDOM_NUMBER(iRan)
            PartStateIntEn(1,iPart) = -LOG_RAN()*DOF_vib*0.5*T_vib*BoltzmannConst
            partWeight = GetParticleWeight(iPart)
            E_vib_new = E_vib_new + partWeight * PartStateIntEn(1,iPart)
          END DO
        END IF
      ELSE
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DO iLoop = 1, nPartNew  ! temporal continuous energy distribution
            iPart = iPartIndx_Node(iLoop)
            DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
              EnergyTemp_vibPoly(iDOF,iLoop) = VibQuantsPar(iPart)%Quants(iDOF) * BoltzmannConst &
                                             * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            END DO
          END DO
        END IF
      END IF
      alpha = E_vib/E_vib_new
      DO iLoop = 1, nPartNew
        iPart = iPartIndx_Node(iLoop)
        partWeight = GetParticleWeight(iPart)
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          PartStateIntEn(1,iPart) = 0.0
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
            betaV = alpha*EnergyTemp_vibPoly(iDOF,iLoop)/(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! CALL RANDOM_NUMBER(iRan)
            ! iQua = INT(betaV+iRan)
            iQua = INT(betaV)
            IF(iQua.GT.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)) iQua=PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF(iDOF)
            PartStateIntEn( 1,iPart)  = PartStateIntEn( 1,iPart) &
               + iQua*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst
            VibQuantsPar(iPart)%Quants(iDOF) = iQua
            Energy_Sum = Energy_Sum - iQua*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)*BoltzmannConst*partWeight
          END DO
          PartStateIntEn( 1,iPart)  = PartStateIntEn( 1,iPart) &
               + SpecDSMC(iSpec)%EZeroPoint
        ELSE  ! Diatomic molecules
          betaV = alpha*PartStateIntEn(1,iPart)/(SpecDSMC(iSpec)%CharaTVib*BoltzmannConst)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! CALL RANDOM_NUMBER(iRan)
          ! iQua = INT(betaV+iRan)
          ! IF((betaV-INT(betaV)).LE.iRan) THEN
          !   iQua = INT(betaV)
          ! ELSE
          !   iQua = INT(betaV) + 1
          ! END IF
          iQua = INT(betaV)
          IF (iQua.GT.SpecDSMC(iSpec)%MaxVibQuant) iQua = SpecDSMC(iSpec)%MaxVibQuant
          PartStateIntEn(1,iPart)  = (iQua + DSMC%GammaQuant)*SpecDSMC(iSpec)%CharaTVib*BoltzmannConst
          Energy_Sum = Energy_Sum - (PartStateIntEn(1,iPart) - SpecDSMC(iSpec)%EZeroPoint)*partWeight
        END IF ! SpecDSMC(1)%PolyatomicMol
      END DO
    END IF
    CellEvib_vMPF(iSpec, iElem) = CellEvib_vMPF(iSpec, iElem) + Energy_Sum
    Energy_Sum = 0.0

! 6.3) ensuring rotational excitation
    alpha = E_rot/E_rot_new
    DO iLoop = 1, nPartNew
      iPart = iPartIndx_Node(iLoop)
      partWeight = GetParticleWeight(iPart)
      PartStateIntEn(2,iPart)  = alpha * PartStateIntEn(2,iPart)
    END DO
  END IF
END IF

! 6.4) new translation energy
! Sanity check: catch problem when bulk of particles consists solely of clones (all have the same velocity vector)
alpha = 0.! Initialize
IF (E_trans.GT.0. .AND. E_trans_new.GT.0. ) THEN
  IF (ISFINITE(SQRT(E_trans/E_trans_new))) THEN; alpha = SQRT(E_trans/E_trans_new)
  ELSE                                         ; alpha = 0.; END IF
END IF
#ifdef CODE_ANALYZE
Energy_new = CellEvib_vMPF(iSpec, iElem) + CellEelec_vMPF(iSpec, iElem)
#endif /* CODE_ANALYZE */

DO iLoop = 1, nPartNew
  iPart = iPartIndx_Node(iLoop)
  PartState(4:6,iPart) = vBulk(1:3) + alpha*(PartState(4:6,iPart)-vBulk_new(1:3))
#ifdef CODE_ANALYZE
  partWeight = GetParticleWeight(iPart)
  ! Energy conservation
  Energy_new = Energy_new + 0.5*Species(iSpec)%MassIC * DOT_PRODUCT(PartState(4:6,iPart),PartState(4:6,iPart)) * partWeight
  IF(CollisMode.GT.1) THEN
    IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
      Energy_new = Energy_new + (PartStateIntEn(1,iPart) + PartStateIntEn(2,iPart)) * partWeight
    END IF
    IF(DSMC%ElectronicModel.GT.0) Energy_new = Energy_new + PartStateIntEn(3,iPart)*partWeight
  END IF
  ! Momentum conservation
  Momentum_new(1:3) = Momentum_new(1:3) + Species(iSpec)%MassIC * PartState(4:6,iPart) * partWeight
  totalWeightNew = totalWeightNew + partWeight
#endif /* CODE_ANALYZE */
END DO

#ifdef CODE_ANALYZE
  ! Check weights
  IF (.NOT.ALMOSTEQUALRELATIVE(totalWeight,totalWeightNew,RelWeightTol)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " totalWeight            : ",totalWeight
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " totalWeightNew         : ",totalWeightNew
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. weight difference : ",totalWeightNew-totalWeight
    ASSOCIATE( weight => MAX(ABS(totalWeight),ABS(totalWeightNew)) )
      IF(weight.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. weight difference : ",(totalWeightNew-totalWeight)/weight
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelWeightTol
    IPWRITE(UNIT_StdOut,*)                     " Old/new particle number: ",nPart, nPartNew
    IPWRITE(UNIT_StdOut,*)                     " Species                : ",iSpec
    IPWRITE(UNIT_StdOut,*)                     " iElem                  : ",iElem
    CALL abort(__STAMP__,'CODE_ANALYZE: part merge is not mass conserving!')
  END IF

  ! Check for energy difference
  IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,RelEneTol)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_new-Energy_old
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " alpha                  : ",alpha
    ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
      IF(energy.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_new-Energy_old)/energy
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelEneTol
    IPWRITE(UNIT_StdOut,*)                     " Old/new particle number: ",nPart, nPartNew
    IPWRITE(UNIT_StdOut,*)                     " Species                : ",iSpec
    CALL abort(__STAMP__,'CODE_ANALYZE: part merge is not energy conserving!')
  END IF
  ! Check for momentum difference
  IF(Symmetry%Order.EQ.3) THEN
    ! Do not check the momentum in z as it can be very small (close to machine precision), leading to greater relative errors
    iMomDim = 3
  ELSE IF(Symmetry%Order.EQ.2) THEN
    iMomDim = 2
  ELSE
    iMomDim = 1
  END IF
  DO iMom=1,iMomDim
    IF (.NOT.ALMOSTEQUALRELATIVE(Momentum_old(iMom),Momentum_new(iMom),RelMomTol)) THEN
      WRITE(UNIT_StdOut,*) '\n'
      IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Direction (x,y,z)        : ",iMom
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_old             : ",Momentum_old(iMom)
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Momentum_new             : ",Momentum_new(iMom)
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Momentum difference : ",Momentum_new(iMom)-Momentum_old(iMom)
      ASSOCIATE( Momentum => MAX(ABS(Momentum_old(iMom)),ABS(Momentum_new(iMom))) )
        IF(Momentum.GT.0.0)THEN
          IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Momentum difference : ",(Momentum_new(iMom)-Momentum_old(iMom))/Momentum
        END IF
      END ASSOCIATE
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelMomTol
      CALL abort(__STAMP__,'CODE_ANALYZE: part merge is not momentum conserving!')
    END IF
  END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE MergeParticles


!===================================================================================================================================
!> Routine for split particles
!> Split particle == clone particle randomly until new particle number is reached and adjust MPF accordingly
!===================================================================================================================================
SUBROUTINE SplitParticles(iPartIndx_Node, nPartIn, nPartNew)
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PartMPF, PartSpecies, PEM, PartPosRef, vMPFSplitLimit
USE MOD_Particle_Vars         ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_DSMC_Vars             ,ONLY: PartStateIntEn, CollisMode, SpecDSMC, DSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_Particle_Tracking_Vars,ONLY: TrackingMethod
USE MOD_Part_Tools            ,ONLY: GetNextFreePosition
!#ifdef CODE_ANALYZE
!USE MOD_Globals               ,ONLY: unit_stdout,myrank,abort
!USE MOD_Particle_Vars         ,ONLY: Symmetry
!#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: nPartIn, nPartNew
INTEGER, INTENT(INOUT)               :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: iRan
INTEGER               :: nSplit, iPart, iNewPart, PartIndx, PositionNbr, LocalElemID, nPart
!===================================================================================================================================
! split particles randomly (until nPartNew is reached)
iNewPart = 0
nPart = nPartIn
nSplit = nPartNew - nPart
DO iNewPart=1,nSplit
  ! Get a random particle (only from the initially available)
  CALL RANDOM_NUMBER(iRan)
  iPart = INT(iRan*nPart) + 1
  PartIndx = iPartIndx_Node(iPart)
  ! Check whether the weighting factor would drop below vMPFSplitLimit (can be below one)
  ! If the resulting MPF is less than 1 you go sub-atomic where the concept of time and space become irrelevant...
  IF((PartMPF(PartIndx) / 2.).LT.vMPFSplitLimit) THEN
    ! Skip this particle
    iPartIndx_Node(iPart) = iPartIndx_Node(nPart)
    nPart = nPart - 1
    ! Leave the DO WHILE loop, if every original particle has been split up to MPF = vMPFSplitLimit (can be below one)
    IF(nPart.EQ.0) EXIT
    ! Limit is reached, therefore go to next particle in reduced list
    CYCLE
  END IF
  IF((PartMPF(PartIndx) / 2.).LT.vMPFSplitLimit) CALL abort(__STAMP__,'Particle split below limit: PartMPF(PartIndx) / 2.=',&
      RealInfoOpt=PartMPF(PartIndx) / 2.)
  PartMPF(PartIndx) = PartMPF(PartIndx) / 2.   ! split particle
  PositionNbr = GetNextFreePosition()
  PartState(1:6,PositionNbr) = PartState(1:6,PartIndx)
  IF(TrackingMethod.EQ.REFMAPPING) PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,PartIndx)
  PartSpecies(PositionNbr) = PartSpecies(PartIndx)
  IF(CollisMode.GT.1) THEN
    PartStateIntEn(1:2,PositionNbr) = PartStateIntEn(1:2,PartIndx)
    IF(SpecDSMC(PartSpecies(PositionNbr))%PolyatomicMol) THEN
      IF(ALLOCATED(VibQuantsPar(PositionNbr)%Quants)) DEALLOCATE(VibQuantsPar(PositionNbr)%Quants)
      ALLOCATE(VibQuantsPar(PositionNbr)%Quants(PolyatomMolDSMC(SpecDSMC(PartSpecies(PositionNbr))%SpecToPolyArray)%VibDOF))
      VibQuantsPar(PositionNbr)%Quants(:) = VibQuantsPar(PartIndx)%Quants(:)
    END IF
    IF(DSMC%ElectronicModel.GT.0) PartStateIntEn(3,PositionNbr) = PartStateIntEn(3,PartIndx)
  END IF
  PartMPF(PositionNbr) = PartMPF(PartIndx)
  IF(UseVarTimeStep) PartTimeStep(PositionNbr) = PartTimeStep(PartIndx)
  PEM%GlobalElemID(PositionNbr) = PEM%GlobalElemID(PartIndx)
  PEM%LastGlobalElemID(PositionNbr) = PEM%GlobalElemID(PartIndx)
  LocalElemID = PEM%LocalElemID(PositionNbr)
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PDM%IsNewPart(PositionNbr)       = .FALSE.
  PDM%dtFracPush(PositionNbr)      = .FALSE.
  PEM%pNext(PEM%pEnd(LocalElemID)) = PositionNbr     ! Next Particle of same Elem (Linked List)
  PEM%pEnd(LocalElemID) = PositionNbr
  PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
END DO


END SUBROUTINE SplitParticles


#ifdef WIP
!===================================================================================================================================
!> Calculation of distribution moments
!> 1.) Calc bulk velocity
!> 2.) Summing up the relative velocities and their squares to calculate the moments (PressTens, HeatVec)
!> 3.) Fill missing entries in PressTens
!===================================================================================================================================
SUBROUTINE CalculateDistMoments(iPartIndx_Node, nPart, vBulk, Vtherm2, PressTens, HeatVec, Energy)
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(INOUT)                       :: vBulk(3), Energy, Vtherm2, PressTens(3,3), HeatVec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: V_rel(3), vmag2
INTEGER               :: iLoop,fillMa1, fillMa2
REAL                  :: partWeight, totalWeight
!===================================================================================================================================
Vtherm2 = 0.0; PressTens = 0.0; HeatVec = 0.0
vBulk = 0.0; totalWeight = 0.0; Energy = 0.

! 1.) calc bulkvelocity
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(4:6,iPartIndx_Node(iLoop)) * partWeight
END DO
vBulk(1:3) = vBulk(1:3)/ totalWeight

! 2.) Summing up the relative velocities and their square to calculate the moments (PressTens, HeatVec)
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  V_rel(1:3)=PartState(4:6,iPartIndx_Node(iLoop))-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  Vtherm2 = Vtherm2 + vmag2 * partWeight
  DO fillMa1 =1, 3
    DO fillMa2 =fillMa1, 3
      PressTens(fillMa1, fillMa2)= PressTens(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2) * partWeight
    END DO
  END DO
  HeatVec(1:3) = HeatVec(1:3) + V_rel(1:3)*vmag2 * partWeight
  Energy = Energy + 0.5 * vmag2 * partWeight
  ! sample inner energies here!
END DO
IF(nPart.GT.2) THEN
  HeatVec = HeatVec*nPart*nPart/((nPart-1.)*(nPart-2.)*totalWeight)
ELSE
  HeatVec = 0.0
END IF
Vtherm2 = Vtherm2*nPart/((nPart-1.)*totalWeight)
! 3.) Fill missing entries in PressTens
PressTens(2,1)=PressTens(1,2)
PressTens(3,1)=PressTens(1,3)
PressTens(3,2)=PressTens(2,3)
PressTens = PressTens/totalWeight

END SUBROUTINE CalculateDistMoments
#endif /*WIP*/


END MODULE MOD_vMPF
