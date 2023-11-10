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

MODULE MOD_DSMC_Collis
!===================================================================================================================================
! Module including collisions, relaxation and reaction decision
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_perform_collision
  MODULE PROCEDURE DSMC_perform_collision
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_perform_collision
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_Elastic_Col(iPair)
!===================================================================================================================================
! Performs simple elastic collision (CollisMode = 1)
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, RadialWeighting
USE MOD_Particle_Vars           ,ONLY: PartSpecies, PartState, UseVarTimeStep, Species, usevMPF
USE MOD_DSMC_CollisVec          ,ONLY: PostCollVec
USE MOD_part_tools              ,ONLY: GetParticleWeight
#ifdef CODE_ANALYZE
USE MOD_Globals                 ,ONLY: Abort
USE MOD_Globals                 ,ONLY: unit_stdout,myrank
USE MOD_Particle_Vars           ,ONLY: Symmetry
#endif /* CODE_ANALYZE */
USE MOD_DSMC_Vars               ,ONLY: DSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: cRelaNew(3)                      ! post-collision relative velocities
  INTEGER                       :: iPart1, iPart2, iSpec1, iSpec2   ! Colliding particles 1 and 2, their species
#ifdef CODE_ANALYZE
REAL,PARAMETER                :: RelMomTol=5e-9  ! Relative tolerance applied to conservation of momentum before/after reaction
REAL                          :: Momentum_old(3),Momentum_new(3)
INTEGER                       :: iMom, iMomDim
#endif /* CODE_ANALYZE */
!===================================================================================================================================
! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
IF (DSMC%ReservoirSimu.AND.DSMC%ReservoirSimuRate) RETURN

iPart1 = Coll_pData(iPair)%iPart_p1
iPart2 = Coll_pData(iPair)%iPart_p2

iSpec1 = PartSpecies(iPart1)
iSpec2 = PartSpecies(iPart2)

#ifdef CODE_ANALYZE
  ! Momentum conservation
  Momentum_old(1:3) = Species(iSpec1)%MassIC * PartState(4:6,iPart1) * GetParticleWeight(iPart1) &
                    + Species(iSpec2)%MassIC * PartState(4:6,iPart2) * GetParticleWeight(iPart2)
#endif /* CODE_ANALYZE */

  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    FracMassCent1 = Species(iSpec1)%MassIC * GetParticleWeight(iPart1) / (Species(iSpec1)%MassIC &
                  * GetParticleWeight(iPart1) + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
    FracMassCent2 = Species(iSpec2)%MassIC *GetParticleWeight(iPart2) / (Species(iSpec1)%MassIC  &
                  * GetParticleWeight(iPart1) + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(PartSpecies(iPart1), Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(PartSpecies(iPart2), Coll_pData(iPair)%PairType)
  END IF
  !Calculation of center of mass velocity
  VeloMx = FracMassCent1 * PartState(4,iPart1) + FracMassCent2 * PartState(4,iPart2)
  VeloMy = FracMassCent1 * PartState(5,iPart1) + FracMassCent2 * PartState(5,iPart2)
  VeloMz = FracMassCent1 * PartState(6,iPart1) + FracMassCent2 * PartState(6,iPart2)

  cRelaNew(1:3) = PostCollVec(iPair)

 ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  PartState(4,iPart1) = VeloMx + FracMassCent2 * cRelaNew(1)
  PartState(5,iPart1) = VeloMy + FracMassCent2 * cRelaNew(2)
  PartState(6,iPart1) = VeloMz + FracMassCent2 * cRelaNew(3)
 ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  PartState(4,iPart2) = VeloMx - FracMassCent1 * cRelaNew(1)
  PartState(5,iPart2) = VeloMy - FracMassCent1 * cRelaNew(2)
  PartState(6,iPart2) = VeloMz - FracMassCent1 * cRelaNew(3)
#ifdef CODE_ANALYZE
  Momentum_new(1:3) = Species(iSpec2)%MassIC* (/VeloMx - FracMassCent1*cRelaNew(1),&
                                                VeloMy - FracMassCent1*cRelaNew(2),&
                                                VeloMz - FracMassCent1*cRelaNew(3)/) * GetParticleWeight(iPart2) &
                    + Species(iSpec1)%MassIC* (/VeloMx + FracMassCent2*cRelaNew(1),&
                                                VeloMy + FracMassCent2*cRelaNew(2),&
                                                VeloMz + FracMassCent2*cRelaNew(3)/) * GetParticleWeight(iPart1)
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
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Momentum difference : ",Momentum_old(iMom)-Momentum_new(iMom)
      ASSOCIATE( Momentum => MAX(ABS(Momentum_old(iMom)),ABS(Momentum_new(iMom))) )
        IF(Momentum.GT.0.0)THEN
          IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Momentum difference : ",(Momentum_old(iMom)-Momentum_new(iMom))/Momentum
        END IF
      END ASSOCIATE
      IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",RelMomTol
      IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 1              : ",iSpec1
      IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 2              : ",iSpec2
      CALL abort(&
          __STAMP__&
          ,'CODE_ANALYZE: DSMC_Elastic_Col is not momentum conserving')
    END IF
  END DO
#endif /* CODE_ANALYZE */

END SUBROUTINE DSMC_Elastic_Col


SUBROUTINE DSMC_Relax_Col_LauxTSHO(iPair)
!===================================================================================================================================
! Performs inelastic collisions with energy exchange (CollisMode = 2/3), allows the relaxation of both collision partners
! Vibrational (of the relaxing molecule), rotational and relative translational energy (of both molecules) is redistributed (V-R-T)
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, RadialWeighting
USE MOD_Particle_Vars         ,ONLY: PartSpecies, PartState, Species, UseVarTimeStep, PEM, usevMPF
USE MOD_DSMC_ElectronicModel  ,ONLY: ElectronicEnergyExchange, TVEEnergyExchange
USE MOD_DSMC_PolyAtomicModel  ,ONLY: DSMC_RotRelaxPoly, DSMC_VibRelaxPoly
USE MOD_DSMC_Relaxation       ,ONLY: DSMC_VibRelaxDiatomic, DSMC_calc_P_rot, DSMC_calc_P_vib, DSMC_calc_P_elec
USE MOD_DSMC_CollisVec        ,ONLY: PostCollVec
USE MOD_part_tools            ,ONLY: GetParticleWeight
USE MOD_MCC_Vars              ,ONLY: UseMCC, SpecXSec
USE MOD_MCC_XSec              ,ONLY: XSec_CalcElecRelaxProb, XSec_ElectronicRelaxation
USE MOD_MCC_Vars              ,ONLY: XSec_Relaxation
USE MOD_Particle_Analyze_Vars ,ONLY: CalcRelaxProb
#ifdef CODE_ANALYZE
USE MOD_Globals               ,ONLY: Abort
USE MOD_Globals               ,ONLY: unit_stdout,myrank
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
REAL (KIND=8)                 :: iRan                             ! Random number
LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2   ! Check whether rot or vib relax is performed
REAL (KIND=8)                 :: Xi_rel, Xi, FakXi                ! Factors of DOF
REAL                          :: cRelaNew(3)                       ! random relative velocity
REAL                          :: ReducedMass
REAL                          :: ProbRot1, ProbRotMax1, ProbRot2, ProbRotMax2, ProbVib1, ProbVib2, ProbElec1, ProbElec2
INTEGER                       :: iCase, iSpec1, iSpec2, iPart1, iPart2, iElem ! Colliding particles 1 and 2 and their species
! variables for electronic level relaxation and transition
INTEGER                       :: ElecLevelRelax
LOGICAL                       :: DoElec1, DoElec2
#ifdef CODE_ANALYZE
REAL                          :: Energy_old,Energy_new
REAL                          :: Weight1, Weight2
#endif /* CODE_ANALYZE */
!===================================================================================================================================
  ElecLevelRelax = 0 ! initialize

  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2
  iSpec1 = PartSpecies(iPart1)
  iSpec2 = PartSpecies(iPart2)
  iElem  = PEM%LocalElemID(iPart1)
  iCase = CollInf%Coll_Case(iSpec1,iSpec2)

  DoRot1  = .FALSE.
  DoRot2  = .FALSE.
  DoVib1  = .FALSE.
  DoVib2  = .FALSE.
  DoElec1 = .FALSE.
  DoElec2 = .FALSE.

  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    ReducedMass = (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) * Species(iSpec2)%MassIC*GetParticleWeight(iPart2)) &
                / (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) + Species(iSpec2)%MassIC*GetParticleWeight(iPart2))
  ELSE
    ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF

#ifdef CODE_ANALYZE
  Weight1 = GetParticleWeight(iPart1)
  Weight2 = GetParticleWeight(iPart2)
  ! Energy conservation
  Energy_old=0.5*Species(iSpec1)%MassIC*DOT_PRODUCT(PartState(4:6,iPart1),PartState(4:6,iPart1)) * Weight1 &
            +0.5*Species(iSpec2)%MassIC*DOT_PRODUCT(PartState(4:6,iPart2),PartState(4:6,iPart2)) * Weight2 &
            + (PartStateIntEn(1,iPart1) + PartStateIntEn(2,iPart1)) * Weight1 &
            + (PartStateIntEn(1,iPart2) + PartStateIntEn(2,iPart2)) * Weight2
  IF(DSMC%ElectronicModel.GT.0) Energy_old=Energy_old + PartStateIntEn(3,iPart1)*Weight1 + PartStateIntEn(3,iPart2) * Weight2
#endif /* CODE_ANALYZE */
  Xi_rel = 2*(2. - CollInf%omega(iSpec1,iSpec2))
  ! DOF of relative motion in VHS model

  Coll_pData(iPair)%Ec = 0.5 * ReducedMass* Coll_pData(iPair)%cRela2

  Xi = Xi_rel !Xi are all DOF in the collision

!--------------------------------------------------------------------------------------------------!
! Decision if Rotation, Vibration and Electronic Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! ELECTRONIC
!--------------------------------------------------------------------------------------------------!
  IF ((DSMC%ElectronicModel.EQ.1).OR.(DSMC%ElectronicModel.EQ.2).OR.(DSMC%ElectronicModel.EQ.3)) THEN
    ! Model 1/2
    IF((SpecDSMC(iSpec1)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec1)%FullyIonized)) THEN
      SELECT CASE(DSMC%ElectronicModel)
      CASE(1)
        CALL RANDOM_NUMBER(iRan)
        CALL DSMC_calc_P_elec(iSpec1, iSpec2, ProbElec1)
        IF (ProbElec1.GT.iRan) DoElec1 = .TRUE.
      CASE(2)
        DoElec1 = .TRUE.
      END SELECT
    END IF
    IF((SpecDSMC(iSpec2)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec2)%FullyIonized)) THEN
      SELECT CASE(DSMC%ElectronicModel)
      CASE(1)
        CALL RANDOM_NUMBER(iRan)
        CALL DSMC_calc_P_elec(iSpec2, iSpec1, ProbElec2)
        IF (ProbElec2.GT.iRan) DoElec2 = .TRUE.
      CASE(2)
        DoElec2 = .TRUE.
      END SELECT
    END IF
    ! Model 3: Cross-section based relaxation
    IF(UseMCC) THEN
      IF(SpecXSec(iCase)%UseElecXSec) THEN
        CALL XSec_ElectronicRelaxation(iPair,iCase,iPart1,iPart2,DoElec1,DoElec2,ElecLevelRelax)
      END IF      ! SpecXSec(iCase)%UseElecXSec
    END IF        ! UseMCC
  END IF
!--------------------------------------------------------------------------------------------------!
! ROTATIONAL + VIBRATIONAL
!--------------------------------------------------------------------------------------------------!
  IF((SpecDSMC(iSpec1)%InterID.EQ.2).OR.(SpecDSMC(iSpec1)%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_rot(iSpec1, iSpec2, iPair, Coll_pData(iPair)%iPart_p1, Xi_rel, ProbRot1, ProbRotMax1)
    IF(ProbRot1.GT.iRan) THEN
      DoRot1 = .TRUE.
      IF ((DSMC%DoTEVRRelaxation).OR.(.NOT.(DoElec1.OR.DoElec2))) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart1) * GetParticleWeight(iPart1)
        Xi = Xi + SpecDSMC(iSpec1)%Xi_Rot
      END IF
    END IF
    IF(DSMC%CalcQualityFactors.AND.(DSMC%RotRelaxProb.GE.2)) THEN
      DSMC%CalcRotProb(iSpec1,2) = MAX(DSMC%CalcRotProb(iSpec1,2),ProbRot1)
      DSMC%CalcRotProb(iSpec1,1) = DSMC%CalcRotProb(iSpec1,1) + ProbRot1
      DSMC%CalcRotProb(iSpec1,3) = DSMC%CalcRotProb(iSpec1,3) + 1
    END IF

    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_vib(iPair, iSpec1, iSpec2, Xi_rel, iElem, ProbVib1)
    IF(ProbVib1.GT.iRan) DoVib1 = .TRUE.
  END IF

  IF (DSMC%ReservoirSimu) THEN
    IF(CalcRelaxProb) THEN
      IF(XSec_Relaxation) THEN
        IF(DoVib1) THEN
          SpecXSec(iCase)%VibCount = SpecXSec(iCase)%VibCount + 1.0
        END IF
      END IF
    END IF
  END IF

  IF((SpecDSMC(iSpec2)%InterID.EQ.2).OR.(SpecDSMC(iSpec2)%InterID.EQ.20)) THEN
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_rot(iSpec2, iSpec1, iPair, Coll_pData(iPair)%iPart_p2, Xi_rel, ProbRot2, ProbRotMax2)
    IF(ProbRot2.GT.iRan) THEN
      DoRot2 = .TRUE.
      IF ((DSMC%DoTEVRRelaxation).OR.(.NOT.(DoElec1.OR.DoElec2))) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart2) * GetParticleWeight(iPart2)
        Xi = Xi + SpecDSMC(iSpec2)%Xi_Rot
      END IF
    END IF
    IF(DSMC%CalcQualityFactors.AND.(DSMC%RotRelaxProb.GE.2)) THEN
      DSMC%CalcRotProb(iSpec2,2) = MAX(DSMC%CalcRotProb(iSpec2,2),ProbRot2)
      DSMC%CalcRotProb(iSpec2,1) = DSMC%CalcRotProb(iSpec2,1) + ProbRot2
      DSMC%CalcRotProb(iSpec2,3) = DSMC%CalcRotProb(iSpec2,3) + 1
    END IF
    CALL RANDOM_NUMBER(iRan)
    CALL DSMC_calc_P_vib(iPair, iSpec2, iSpec1, Xi_rel, iElem, ProbVib2)
    IF(ProbVib2.GT.iRan) DoVib2 = .TRUE.
  END IF

! Cross-section model only allows one relaxation process during a single collision
IF (DSMC%ElectronicModel.EQ.3) THEN
  IF(DoElec1.OR.DoElec2) THEN
    IF(SpecXSec(iCase)%UseVibXSec) THEN
      DoVib1 = .FALSE.
      DoVib2 = .FALSE.
    END IF
  END IF
END IF

  FakXi = 0.5*Xi  - 1.  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

IF (DSMC%ReservoirSimu) THEN
  IF(CalcRelaxProb) THEN
    IF(XSec_Relaxation) THEN
      IF(DoVib2) THEN
        SpecXSec(iCase)%VibCount = SpecXSec(iCase)%VibCount + 1.0
      END IF
    END IF
  END IF
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
  IF (DSMC%ReservoirSimuRate) RETURN
END IF

!--------------------------------------------------------------------------------------------------!
! Electronic Relaxation / Transition
!--------------------------------------------------------------------------------------------------!
  IF (DSMC%DoTEVRRelaxation) THEN
    IF(.NOT.SpecDSMC(iSpec1)%PolyatomicMol) THEN
      IF(DoElec1.AND.DoVib1) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart1)  + PartStateIntEn(1,iPart1)
        CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,iPart1,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart1)  - PartStateIntEn(1,iPart1)
        DoElec1=.false.
        DoVib1=.false.
      END IF !DoElec1.AND.DoVib1
    END IF ! .NOT.SpecDSMC(iSpec1)%PolyatomicMol

    IF(.NOT.SpecDSMC(iSpec2)%PolyatomicMol) THEN
      IF(DoElec2.AND.DoVib2) THEN
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart2) + PartStateIntEn(1,iPart2)
        CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,iPart2,FakXi)
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart2) - PartStateIntEn(1,iPart2)
        DoElec2=.false.
        DoVib2=.false.
      END IF ! .NOT.SpecDSMC(iSpec2)%PolyatomicMol
    END IF !DoElec2.AND.DoVib2
  END IF ! DSMC%DoTEVRRelaxation

  ! Relaxation of first particle
  IF ( DoElec1 ) THEN
    ! calculate energy for electronic relaxation of particle 1
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart1) * GetParticleWeight(iPart1)
    CALL ElectronicEnergyExchange(iPair, Coll_pData(iPair)%iPart_p1, FakXi, XSec_Level=ElecLevelRelax)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart1) * GetParticleWeight(iPart1)
  END IF

  ! Electronic relaxation of second particle
  IF ( DoElec2 ) THEN
    ! calculate energy for electronic relaxation of particle 2
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart2) * GetParticleWeight(iPart2)
    CALL ElectronicEnergyExchange(iPair, Coll_pData(iPair)%iPart_p2, FakXi, XSec_Level=ElecLevelRelax)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart2) * GetParticleWeight(iPart2)
  END IF

  IF ((DoElec1.OR.DoElec2).AND.DoRot1.AND.(.NOT.DSMC%DoTEVRRelaxation)) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart1) * GetParticleWeight(iPart1)
    FakXi = FakXi + 0.5*SpecDSMC(iSpec1)%Xi_Rot
  END IF
  IF ((DoElec1.OR.DoElec2).AND.DoRot2.AND.(.NOT.DSMC%DoTEVRRelaxation)) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart2) * GetParticleWeight(iPart2)
    FakXi = FakXi + 0.5*SpecDSMC(iSpec2)%Xi_Rot
  END IF

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoVib1) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart1) * GetParticleWeight(iPart1)
    IF(SpecDSMC(iSpec1)%PolyatomicMol) THEN
      CALL DSMC_VibRelaxPoly(iPair, Coll_pData(iPair)%iPart_p1,FakXi)
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair, iPart1,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(1,iPart1) * GetParticleWeight(iPart1)
  END IF

  IF(DoVib2) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart2) * GetParticleWeight(iPart2)
    IF(SpecDSMC(iSpec2)%PolyatomicMol) THEN
      CALL DSMC_VibRelaxPoly(iPair, Coll_pData(iPair)%iPart_p2,FakXi)
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair, iPart2,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(1,iPart2) * GetParticleWeight(iPart2)
  END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!
  IF(DoRot1) THEN
    IF(SpecDSMC(iSpec1)%PolyatomicMol.AND.(SpecDSMC(iSpec1)%Xi_Rot.EQ.3)) THEN
      FakXi = FakXi - 0.5*SpecDSMC(iSpec1)%Xi_Rot
      CALL DSMC_RotRelaxPoly(iPair, Coll_pData(iPair)%iPart_p1, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,Coll_pData(iPair)%iPart_p1)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(2,Coll_pData(iPair)%iPart_p1) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,Coll_pData(iPair)%iPart_p1)
      FakXi = FakXi - 0.5*SpecDSMC(iSpec1)%Xi_Rot
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
      PartStateIntEn(2,iPart1) = PartStateIntEn(2,iPart1)/GetParticleWeight(iPart1)
    END IF
  END IF

  IF(DoRot2) THEN
    IF(SpecDSMC(iSpec2)%PolyatomicMol.AND. &
        (SpecDSMC(iSpec2)%Xi_Rot.EQ.3)) THEN
      FakXi = FakXi - 0.5*SpecDSMC(iSpec2)%Xi_Rot
      CALL DSMC_RotRelaxPoly(iPair, Coll_pData(iPair)%iPart_p2, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,Coll_pData(iPair)%iPart_p2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(2,Coll_pData(iPair)%iPart_p2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,Coll_pData(iPair)%iPart_p2)
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
      PartStateIntEn(2,iPart2) = PartStateIntEn(2,iPart2)/GetParticleWeight(iPart2)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    FracMassCent1 = Species(iSpec1)%MassIC *GetParticleWeight(iPart1)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
    FracMassCent2 = Species(iSpec2)%MassIC *GetParticleWeight(iPart2)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(iSpec1, Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(iSpec2, Coll_pData(iPair)%PairType)
  END IF

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(4,iPart1) + FracMassCent2 * PartState(4,iPart2)
  VeloMy = FracMassCent1 * PartState(5,iPart1) + FracMassCent2 * PartState(5,iPart2)
  VeloMz = FracMassCent1 * PartState(6,iPart1) + FracMassCent2 * PartState(6,iPart2)

  Coll_pData(iPair)%cRela2 = 2. * Coll_pData(iPair)%Ec/ReducedMass
  cRelaNew(1:3) = PostCollVec(iPair)

  ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  PartState(4,iPart1) = VeloMx + FracMassCent2*cRelaNew(1)
  PartState(5,iPart1) = VeloMy + FracMassCent2*cRelaNew(2)
  PartState(6,iPart1) = VeloMz + FracMassCent2*cRelaNew(3)
  ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  PartState(4,iPart2) = VeloMx - FracMassCent1*cRelaNew(1)
  PartState(5,iPart2) = VeloMy - FracMassCent1*cRelaNew(2)
  PartState(6,iPart2) = VeloMz - FracMassCent1*cRelaNew(3)

#ifdef CODE_ANALYZE
  Energy_new= 0.5*Species(iSpec2)%MassIC*((VeloMx - FracMassCent1*cRelaNew(1))**2 &
                                        + (VeloMy - FracMassCent1*cRelaNew(2))**2 &
                                        + (VeloMz - FracMassCent1*cRelaNew(3))**2) * Weight2 &
             +0.5*Species(iSpec1)%MassIC*((VeloMx + FracMassCent2*cRelaNew(1))**2 &
                                        + (VeloMy + FracMassCent2*cRelaNew(2))**2 &
                                        + (VeloMz + FracMassCent2*cRelaNew(3))**2) * Weight1 &
                        + (PartStateIntEn(1,iPart1) + PartStateIntEn(2,iPart1)) * Weight1 &
                        + (PartStateIntEn(1,iPart2) + PartStateIntEn(2,iPart2)) * Weight2
  IF(DSMC%ElectronicModel.GT.0) Energy_new = Energy_new + PartStateIntEn(3,iPart1) * Weight1 &
                                                   + PartStateIntEn(3,iPart2) * Weight2
  ! Check for energy difference
  IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,1.0e-12)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_old-Energy_new
    ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
      IF(energy.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_old-Energy_new)/energy
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",1.0e-12
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 1              : ",iSpec1
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 2              : ",iSpec2
    CALL abort(__STAMP__,'CODE_ANALYZE: DSMC_Relaxation with SelectionProcedure = 1 is not energy conserving!')
  END IF
#endif /* CODE_ANALYZE */

END SUBROUTINE DSMC_Relax_Col_LauxTSHO


SUBROUTINE DSMC_Relax_Col_Gimelshein(iPair)
!===================================================================================================================================
! Performs inelastic collisions with energy exchange (CollisMode = 2/3)
! Selection procedure according to Gimelshein et al. (Physics of Fluids, V 14, No 12, 2002: 'Vibrational Relaxation rates in the
! DSMC Method') For further understanding see Zhang, Schwarzentruber, Physics of Fluids, V25, 2013: 'inelastic collision selection
! procedures for DSMC calculation of gas mixtures')
!===================================================================================================================================
! MODULES
  USE MOD_Globals,                ONLY : Abort
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC, PolyatomMolDSMC, SpecDSMC, PartStateIntEn
  USE MOD_DSMC_Vars,              ONLY : RadialWeighting
  USE MOD_Particle_Vars,          ONLY : PartSpecies, PartState, PEM, usevMPF, UseVarTimeStep, Species
  USE MOD_DSMC_PolyAtomicModel,   ONLY : DSMC_RotRelaxPoly, DSMC_VibRelaxPoly, DSMC_VibRelaxPolySingle
  USE MOD_DSMC_Relaxation,        ONLY : DSMC_VibRelaxDiatomic, DSMC_calc_P_rot, DSMC_calc_P_vib, DSMC_calc_P_elec
  USE MOD_DSMC_CollisVec,         ONLY : PostCollVec
  USE MOD_DSMC_ElectronicModel,   ONLY: ElectronicEnergyExchange
  USE MOD_part_tools            ,ONLY: GetParticleWeight
#ifdef CODE_ANALYZE
  USE MOD_Globals                ,ONLY : unit_stdout,myrank
  USE MOD_Particle_Vars          ,ONLY : Species
  USE MOD_part_tools             ,ONLY : GetParticleWeight
#endif /* CODE_ANALYZE */
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iPair
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2                 ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz                       ! center of mass velo
  INTEGER                       :: iDOF, iPolyatMole, DOFRelax, iElem
  REAL (KIND=8)                 :: iRan
  LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2               ! Check whether rot or vib relax is performed
  LOGICAL                       :: DoElec1, DoElec2
  REAL (KIND=8)                 :: FakXi, Xi_rel                                ! Factors of DOF
  REAL                          :: cRelaNew(3),ReducedMass                      ! post collision relative velocity
  REAL                          :: ProbFrac1, ProbFrac2, ProbFrac3, ProbFrac4   ! probability-fractions according to Zhang
  REAL                          :: ProbFrac5, ProbFrac6                         ! probability-fractions according to Zhang
  REAL                          :: ProbRot1, ProbRot2, ProbVib1, ProbVib2       ! probabilities for rot-/vib-relax for part 1/2
  REAL                          :: ProbElec1, ProbElec2
  REAL                          :: BLCorrFact, ProbRotMax1, ProbRotMax2         ! Correction factor for BL-redistribution of energy
  INTEGER                       :: iPart1, iPart2, iSpec1, iSpec2               ! Colliding particles 1 and 2 and their species
#ifdef CODE_ANALYZE
  REAL                          :: Energy_old,Energy_new
  REAL                          :: Weight1, Weight2
#endif /* CODE_ANALYZE */
!===================================================================================================================================
! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
IF (DSMC%ReservoirSimu.AND.DSMC%ReservoirSimuRate) RETURN

  iPart1 = Coll_pData(iPair)%iPart_p1
  iPart2 = Coll_pData(iPair)%iPart_p2
  iSpec1 = PartSpecies(iPart1)
  iSpec2 = PartSpecies(iPart2)
  iElem  = PEM%LocalElemID(Coll_pData(iPair)%iPart_p1)

  ! set some initial values
  DoRot1  = .FALSE.
  DoRot2  = .FALSE.
  DoVib1  = .FALSE.
  DoVib2  = .FALSE.
  DoElec1 = .FALSE.
  DoElec2 = .FALSE.
  ProbVib1 = 0.
  ProbRot1 = 0.
  ProbVib2 = 0.
  ProbRot2 = 0.
  ProbElec1 = 0. ; ProbElec2 = 0.

  Xi_rel = 2.*(2. - CollInf%omega(iSpec1,iSpec2)) ! DOF of relative motion in VHS model
  FakXi  = 0.5*Xi_rel - 1.

  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    ReducedMass = (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) * Species(iSpec2)%MassIC*GetParticleWeight(iPart2)) &
                / (Species(iSpec1)%MassIC*GetParticleWeight(iPart1) + Species(iSpec2)%MassIC*GetParticleWeight(iPart2))
  ELSE
    ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
  END IF

  Coll_pData(iPair)%Ec = 0.5 * ReducedMass* Coll_pData(iPair)%cRela2

#ifdef CODE_ANALYZE
  Weight1 = GetParticleWeight(iPart1)
  Weight2 = GetParticleWeight(iPart2)
  ! Energy conservation
  Energy_old=0.5*Species(iSpec1)%MassIC*DOT_PRODUCT(PartState(4:6,iPart1),PartState(4:6,iPart1)) * Weight1 &
            +0.5*Species(iSpec2)%MassIC*DOT_PRODUCT(PartState(4:6,iPart2),PartState(4:6,iPart2)) * Weight2 &
            + (PartStateIntEn(1,iPart1) + PartStateIntEn(2,iPart1)) * Weight1 &
            + (PartStateIntEn(1,iPart2) + PartStateIntEn(2,iPart2)) * Weight2
  IF(DSMC%ElectronicModel.GT.0) Energy_old=Energy_old + PartStateIntEn(3,iPart1)*Weight1 + PartStateIntEn(3,iPart2) * Weight2
#endif /* CODE_ANALYZE */
!--------------------------------------------------------------------------------------------------!
! Decision if Rotation, Vibration and Electronic Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!

  ! calculate probability for rotational/vibrational relaxation for both particles
  IF ((SpecDSMC(iSpec1)%InterID.EQ.2).OR.(SpecDSMC(iSpec1)%InterID.EQ.20)) THEN
    CALL DSMC_calc_P_vib(iPair, iSpec1, iSpec2, Xi_rel, iElem, ProbVib1)
    CALL DSMC_calc_P_rot(iSpec1, iSpec2, iPair, iPart1, Xi_rel, ProbRot1, ProbRotMax1)
  ELSE
    ProbVib1 = 0.
    ProbRot1 = 0.
  END IF
  IF ((SpecDSMC(iSpec2)%InterID.EQ.2).OR.(SpecDSMC(iSpec2)%InterID.EQ.20)) THEN
    CALL DSMC_calc_P_vib(iPair, iSpec2, iSpec1, Xi_rel, iElem, ProbVib2)
    CALL DSMC_calc_P_rot(iSpec2, iSpec1, iPair, iPart2, Xi_rel, ProbRot2, ProbRotMax2)
  ELSE
    ProbVib2 = 0.
    ProbRot2 = 0.
  END IF
  IF (DSMC%ElectronicModel.EQ.1) THEN
    ! Model 1/2
    IF((SpecDSMC(iSpec1)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec1)%FullyIonized)) THEN
      CALL DSMC_calc_P_elec(iSpec1, iSpec2, ProbElec1)
    END IF
    IF((SpecDSMC(iSpec2)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec2)%FullyIonized)) THEN
      CALL DSMC_calc_P_elec(iSpec2, iSpec1, ProbElec2)
    END IF
  END IF

  ! Calculate probability fractions
  IF(SpecDSMC(iSpec1)%PolyatomicMol) THEN
    IF(DSMC%PolySingleMode) THEN
      ! If single-mode relaxation is considered, every mode has its own probability while it accumulates all the previous ones
      ! Here, the last mode of the molecule has the highest probability, later it is found which exact mode is going to be relaxed
      iPolyatMole = SpecDSMC(iSpec1)%SpecToPolyArray
      ProbFrac1 = PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(PolyatomMolDSMC(iPolyatMole)%VibDOF)
    ELSE
      ProbFrac1 = ProbVib1
    END IF
  ELSE
    ProbFrac1 = ProbVib1
  END IF
  IF(SpecDSMC(iSpec2)%PolyatomicMol) THEN
    IF(DSMC%PolySingleMode) THEN
      iPolyatMole = SpecDSMC(iSpec2)%SpecToPolyArray
      ProbFrac2   = ProbFrac1 + PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(PolyatomMolDSMC(iPolyatMole)%VibDOF)
    ELSE
      ProbFrac2   = ProbFrac1 + ProbVib2
    END IF
  ELSE
    ProbFrac2     = ProbFrac1 + ProbVib2
  END IF
  ProbFrac3       = ProbFrac2 + ProbRot1
  ProbFrac4       = ProbFrac3 + ProbRot2
  IF (DSMC%ElectronicModel.EQ.1) THEN
    ProbFrac5 = ProbFrac4 + ProbElec1
    ProbFrac6 = ProbFrac5 + ProbElec2
  END IF

  ! Check if sum of probabilities is less than 1.
  IF (ProbFrac4.GT. 1.0) THEN
    CALL Abort(__STAMP__,'Error! Sum of internal relaxation probabilities > 1.0 for iPair ',iPair)
  END IF
  IF ((DSMC%ElectronicModel.EQ.1).OR.(DSMC%ElectronicModel.EQ.2)) THEN
    ! Check if sum of probabilities is less than 1.
    IF (ProbFrac6.GT. 1.0) THEN
      CALL Abort(__STAMP__,'Error! Sum of internal relaxation probabilities > 1.0 for iPair ',iPair)
    END IF
  END IF

  ! Select relaxation procedure (vibration, rotation)
  CALL RANDOM_NUMBER(iRan)
  IF(iRan .LT. ProbFrac1) THEN                    !            R1 < A1
    IF (SpecDSMC(iSpec1)%PolyatomicMol.AND.DSMC%PolySingleMode) THEN
      ! Determination which vibrational mode should be relaxed for single-mode polyatomic relaxation
      iPolyatMole = SpecDSMC(iSpec1)%SpecToPolyArray
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF(PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF).GT.iRan) THEN
          DoVib1 = .TRUE.
          DOFRelax = iDOF
          EXIT
        END IF
      END DO
    ELSE
      DoVib1 = .TRUE.
    END IF
  ELSEIF(iRan .LT. ProbFrac2) THEN                !      A1 <= R1 < A2
    IF (SpecDSMC(iSpec2)%PolyatomicMol.AND.DSMC%PolySingleMode) THEN
      ! Determination which vibrational mode should be relaxed for single-mode polyatomic relaxation
      iPolyatMole = SpecDSMC(iSpec2)%SpecToPolyArray
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        IF(ProbFrac1 + PolyatomMolDSMC(iPolyatMole)%VibRelaxProb(iDOF).GT.iRan) THEN
          DoVib2 = .TRUE.
          DOFRelax = iDOF
          EXIT
        END IF
      END DO
    ELSE
      DoVib2 = .TRUE.
    END IF
  ELSEIF(iRan .LT. ProbFrac3) THEN                !      A2 <= R1 < A3
    DoRot1 = .TRUE.
  ELSEIF(iRan .LT. ProbFrac4) THEN                !      A3 <= R1 < A4
    DoRot2 = .TRUE.
  ELSEIF (DSMC%ElectronicModel.EQ.1) THEN
    IF (iRan .LT. ProbFrac5) THEN                !      A3 <= R1 < A4
      DoElec1 = .TRUE.
    ELSEIF (iRan .LT. ProbFrac6) THEN                !      A3 <= R1 < A4
      DoElec2 = .TRUE.
    END IF
  END IF

  IF (DSMC%ElectronicModel.EQ.2) THEN
    IF((SpecDSMC(iSpec1)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec1)%FullyIonized)) DoElec1 = .TRUE.
    IF((SpecDSMC(iSpec2)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec2)%FullyIonized)) DoElec2 = .TRUE.
  END IF

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

  IF(DoVib1) THEN
    ! check if correction term for BL redistribution (depending on relaxation model) is needed
    BLCorrFact = 1.
    ! Adding the interal energy of the particle to be redistributed
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart1)*GetParticleWeight(iPart1)
    IF(SpecDSMC(PartSpecies(iPart1))%PolyatomicMol) THEN
      IF (.NOT.DSMC%PolySingleMode) THEN
        ! --------------------------------------------------------------------------------------------------!
        !  Multi-mode relaxation with the Metropolis-Hastings method
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPoly(iPair,iPart1,FakXi)
      ELSE
        ! --------------------------------------------------------------------------------------------------!
        !  Single-mode relaxation of a previously selected mode
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPolySingle(iPair,iPart1,FakXi,DOFRelax)
      END IF
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair,iPart1,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(1,iPart1)*GetParticleWeight(iPart1)
  END IF

  IF(DoVib2) THEN
    ! check if correction term for BL redistribution (depending on relaxation model) is needed
    BLCorrFact = 1.
    ! Adding the interal energy of the particle to be redistributed (not if single-mode polyatomic relaxation is enabled)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart2)*GetParticleWeight(iPart2)
    IF(SpecDSMC(PartSpecies(iPart2))%PolyatomicMol) THEN
      IF (.NOT.DSMC%PolySingleMode) THEN
        ! --------------------------------------------------------------------------------------------------!
        !  Multi-mode relaxation with the Metropolis-Hastings method
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPoly(iPair,iPart2,FakXi)
      ELSE
        ! --------------------------------------------------------------------------------------------------!
        !  Single-mode relaxation of a previously selected mode
        ! --------------------------------------------------------------------------------------------------!
        CALL DSMC_VibRelaxPolySingle(iPair,iPart2,FakXi,DOFRelax)
      END IF
    ELSE
      CALL DSMC_VibRelaxDiatomic(iPair,iPart2,FakXi)
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(1,iPart2)*GetParticleWeight(iPart2)
  END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!
  IF(DoRot1) THEN
    !check if correction term in distribution (depending on relaxation model) is needed
    IF(DSMC%RotRelaxProb.EQ.3.0) THEN
      BLCorrFact = ProbRot1 / ProbRotMax1
    ELSE
      BLCorrFact = 1.
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart1)*GetParticleWeight(iPart1)
    IF(SpecDSMC(iSpec1)%PolyatomicMol.AND.(SpecDSMC(iSpec1)%Xi_Rot.EQ.3)) THEN
      CALL DSMC_RotRelaxPoly(iPair, iPart1, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,iPart1)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      FakXi = FakXi + 0.5*SpecDSMC(iSpec1)%Xi_Rot
      PartStateIntEn(2,iPart1) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi)*BLCorrFact)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,iPart1)
      FakXi = FakXi - 0.5*SpecDSMC(iSpec1)%Xi_Rot
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
      PartStateIntEn(2,iPart1) = PartStateIntEn(2,iPart1)/GetParticleWeight(iPart1)
    END IF
  END IF

  IF(DoRot2) THEN
    !check if correction term in distribution (depending on relaxation model) is needed
    IF(DSMC%RotRelaxProb.EQ.3.0) THEN
      BLCorrFact = ProbRot2 / ProbRotMax2
    ELSE
      BLCorrFact = 1.
    END IF
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(2,iPart2)*GetParticleWeight(iPart2)
    IF(SpecDSMC(iSpec2)%PolyatomicMol.AND.(SpecDSMC(iSpec2)%Xi_Rot.EQ.3)) THEN
      CALL DSMC_RotRelaxPoly(iPair, iPart2, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,iPart2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      FakXi = FakXi + 0.5*SpecDSMC(iSpec2)%Xi_Rot
      PartStateIntEn(2,iPart2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi)*BLCorrFact)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(2,iPart2)
      FakXi = FakXi - 0.5*SpecDSMC(iSpec2)%Xi_Rot
    END IF
    IF(RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
      PartStateIntEn(2,iPart2) = PartStateIntEn(2,iPart2)/GetParticleWeight(iPart2)
    END IF
  END IF

!--------------------------------------------------------------------------------------------------!
! Electronical Relaxation
!--------------------------------------------------------------------------------------------------!
  ! Relaxation of first particle
  IF ( DoElec1 ) THEN
    ! calculate energy for electronic relaxation of particle 1
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)
    CALL ElectronicEnergyExchange(iPair,iPart1,FakXi)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart1)*GetParticleWeight(iPart1)
  END IF

  ! Electronic relaxation of second particle
  IF ( DoElec2 ) THEN
    ! calculate energy for electronic relaxation of particle 2
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,iPart2)*GetParticleWeight(iPart2)
    CALL ElectronicEnergyExchange(iPair,iPart2,FakXi)
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(3,iPart2)*GetParticleWeight(iPart2)
  END IF
!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  IF (RadialWeighting%DoRadialWeighting.OR.UseVarTimeStep.OR.usevMPF) THEN
    FracMassCent1 = Species(iSpec1)%MassIC *GetParticleWeight(iPart1)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
    FracMassCent2 = Species(iSpec2)%MassIC *GetParticleWeight(iPart2)/(Species(iSpec1)%MassIC *GetParticleWeight(iPart1) &
          + Species(iSpec2)%MassIC *GetParticleWeight(iPart2))
  ELSE
    FracMassCent1 = CollInf%FracMassCent(iSpec1, Coll_pData(iPair)%PairType)
    FracMassCent2 = CollInf%FracMassCent(iSpec2, Coll_pData(iPair)%PairType)
  END IF

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(4,iPart1) + FracMassCent2 * PartState(4,iPart2)
  VeloMy = FracMassCent1 * PartState(5,iPart1) + FracMassCent2 * PartState(5,iPart2)
  VeloMz = FracMassCent1 * PartState(6,iPart1) + FracMassCent2 * PartState(6,iPart2)

  Coll_pData(iPair)%cRela2 = 2. * Coll_pData(iPair)%Ec/ReducedMass
  cRelaNew(1:3) = PostCollVec(iPair)

  ! deltaV particle 1 (post collision particle 1 velocity in laboratory frame)
  PartState(4,iPart1) = VeloMx + FracMassCent2*cRelaNew(1)
  PartState(5,iPart1) = VeloMy + FracMassCent2*cRelaNew(2)
  PartState(6,iPart1) = VeloMz + FracMassCent2*cRelaNew(3)
  ! deltaV particle 2 (post collision particle 2 velocity in laboratory frame)
  PartState(4,iPart2) = VeloMx - FracMassCent1*cRelaNew(1)
  PartState(5,iPart2) = VeloMy - FracMassCent1*cRelaNew(2)
  PartState(6,iPart2) = VeloMz - FracMassCent1*cRelaNew(3)

#ifdef CODE_ANALYZE
  Energy_new= 0.5*Species(PartSpecies(iPart2))%MassIC*((VeloMx - FracMassCent1*cRelaNew(1))**2 &
                                                     + (VeloMy - FracMassCent1*cRelaNew(2))**2 &
                                                     + (VeloMz - FracMassCent1*cRelaNew(3))**2) * Weight2 &
             +0.5*Species(PartSpecies(iPart1))%MassIC*((VeloMx + FracMassCent2*cRelaNew(1))**2 &
                                                     + (VeloMy + FracMassCent2*cRelaNew(2))**2 &
                                                     + (VeloMz + FracMassCent2*cRelaNew(3))**2) * Weight1 &
                        + (PartStateIntEn(1,iPart1) + PartStateIntEn(2,iPart1)) * Weight1 &
                        + (PartStateIntEn(1,iPart2) + PartStateIntEn(2,iPart2)) * Weight2
  IF(DSMC%ElectronicModel.GT.0) Energy_new = Energy_new + PartStateIntEn(3,iPart1) * Weight1 &
                                                   + PartStateIntEn(3,iPart2) * Weight2
  ! Check for energy difference
  IF (.NOT.ALMOSTEQUALRELATIVE(Energy_old,Energy_new,1.0e-12)) THEN
    WRITE(UNIT_StdOut,*) '\n'
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_old             : ",Energy_old
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Energy_new             : ",Energy_new
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " abs. Energy difference : ",Energy_old-Energy_new
    ASSOCIATE( energy => MAX(ABS(Energy_old),ABS(Energy_new)) )
      IF(energy.GT.0.0)THEN
        IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')" rel. Energy difference : ",(Energy_old-Energy_new)/energy
      END IF
    END ASSOCIATE
    IPWRITE(UNIT_StdOut,'(I0,A,ES25.14E3)')    " Applied tolerance      : ",1.0e-12
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 1              : ",iSpec1
    IPWRITE(UNIT_StdOut,'(I0,A,I0)')           " Species 2              : ",iSpec2
    CALL abort(&
        __STAMP__&
        ,'CODE_ANALYZE: gimelshein DSMC_relaxation is not energy conserving')
  END IF
#endif /* CODE_ANALYZE */

END SUBROUTINE DSMC_Relax_Col_Gimelshein


SUBROUTINE DSMC_perform_collision(iPair, iElem, NodeVolume, NodePartNum)
!===================================================================================================================================
! Collision mode is selected (1: Elastic, 2: Non-elastic, 3: Non-elastic with chemical reactions)
!===================================================================================================================================
! MODULES
USE MOD_Globals               ,ONLY: Abort, CROSS
USE MOD_DSMC_Vars             ,ONLY: CollisMode, Coll_pData, SelectionProc
USE MOD_DSMC_Vars             ,ONLY: DSMC
USE MOD_Particle_Vars         ,ONLY: PartState, WriteMacroVolumeValues, Symmetry
USE MOD_Particle_Vars         ,ONLY: UseRotRefFrame, PDM, PartVeloRotRef, RotRefFrameOmega
USE MOD_TimeDisc_Vars         ,ONLY: TEnd, Time
USE MOD_DSMC_Vars             ,ONLY: RadialWeighting
USE MOD_Particle_Vars         ,ONLY: usevMPF, Species, PartSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: CalcCollRates
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(IN)           :: iElem
REAL, INTENT(IN), OPTIONAL    :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: RelaxToDo
INTEGER                       :: iPart1, iPart2                         ! Colliding particles 1 and 2
REAL                          :: Distance
REAL                          :: MacroParticleFactor, PairWeight
!===================================================================================================================================
IF (DSMC%ReservoirSimu) THEN
  IF(CalcCollRates) THEN
    PairWeight = (GetParticleWeight(Coll_pData(iPair)%iPart_p1) + GetParticleWeight(Coll_pData(iPair)%iPart_p2))/2.
    IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
      ! Weighting factor already included in the PairWeight
      MacroParticleFactor = 1.
    ELSE
      ! Weighting factor should be the same for all species anyway (BGG: first species is the non-BGG particle species)
      MacroParticleFactor = Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MacroParticleFactor
    END IF
    DSMC%NumColl(Coll_pData(iPair)%PairType) = DSMC%NumColl(Coll_pData(iPair)%PairType) + PairWeight*MacroParticleFactor
  END IF
END IF

iPart1 = Coll_pData(iPair)%iPart_p1
iPart2 = Coll_pData(iPair)%iPart_p2

IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    IF(Symmetry%Order.EQ.3) THEN
      Distance = SQRT((PartState(1,iPart1) - PartState(1,iPart2))**2 + (PartState(2,iPart1) - PartState(2,iPart2))**2 &
                    + (PartState(3,iPart1) - PartState(3,iPart2))**2)
    ELSE IF(Symmetry%Order.EQ.2) THEN
      Distance = SQRT((PartState(1,iPart1) - PartState(1,iPart2))**2 + (PartState(2,iPart1) - PartState(2,iPart2))**2)
    ELSE
      Distance = ABS(PartState(1,iPart1) - PartState(1,iPart2))
    END IF
    DSMC%CollSepDist = DSMC%CollSepDist + Distance
    DSMC%CollSepCount = DSMC%CollSepCount + 1
  END IF
END IF

SELECT CASE(CollisMode)
  CASE(1) ! elastic collision
    CALL DSMC_Elastic_Col(iPair)
  CASE(2) ! collision with relaxation
    SELECT CASE(SelectionProc)
      CASE(1)
        CALL DSMC_Relax_Col_LauxTSHO(iPair)
      CASE(2)
        CALL DSMC_Relax_Col_Gimelshein(iPair)
      CASE DEFAULT
        CALL Abort(__STAMP__,'ERROR in DSMC_perform_collision: Wrong Selection Procedure:',SelectionProc)
    END SELECT
  CASE(3) ! chemical reactions
    RelaxToDo = .TRUE.
    IF (PRESENT(NodeVolume).AND.PRESENT(NodePartNum)) THEN
      CALL ReactionDecision(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
    ELSE
      CALL ReactionDecision(iPair, RelaxToDo, iElem)
    END IF
    IF (RelaxToDo) THEN
      SELECT CASE(SelectionProc)
        CASE(1)
          CALL DSMC_Relax_Col_LauxTSHO(iPair)
        CASE(2)
          CALL DSMC_Relax_Col_Gimelshein(iPair)
        CASE DEFAULT
          CALL Abort(__STAMP__,'ERROR in DSMC_perform_collision: Wrong Selection Procedure:',SelectionProc)
      END SELECT
    END IF
  CASE DEFAULT
    CALL Abort(__STAMP__,'ERROR in DSMC_perform_collision: Wrong Collision Mode:',CollisMode)
END SELECT

! Rotational frame of reference
IF(UseRotRefFrame) THEN
  ! Transform the new velocity in the inertial frame and reset the velocity in the rotational frame
  IF(PDM%InRotRefFrame(iPart1)) THEN
    PartVeloRotRef(1:3,iPart1) = PartState(4:6,iPart1) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,iPart1))
    PartVeloRotRef(1:3,iPart2) = PartState(4:6,iPart2) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,iPart2))
  END IF
END IF

END SUBROUTINE DSMC_perform_collision


SUBROUTINE ReactionDecision(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
!===================================================================================================================================
!> Decision of reaction path to perform
!> 1.) Calculate the TCE reaction probabilities/test whether any QK reactions are possible (XSec probabilities are treated in the
!>     DSMC_prob_calc subroutine)
!> 2.) Determine which TCE reaction is most likely to occur
!> 3.) Treat QK and TCE reaction paths
!>    a. Logical array which indicates which reactions are possible (QK: reactions above the threshold energy, TCE: reaction
!>       chosen based on probability)
!>    b. If any reaction is to occur, a XSec-based reaction path will not be considered.
!> 4.) Treat XSec reaction paths
!===================================================================================================================================
! MODULES
USE MOD_Globals                 ,ONLY: Abort
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, CollInf, ChemReac, RadialWeighting
USE MOD_MCC_Vars                ,ONLY: SpecXSec
USE MOD_Particle_Vars           ,ONLY: Species, PartSpecies, PEM, UseVarTimeStep, usevMPF, VarTimeStep
USE MOD_DSMC_ChemReact          ,ONLY: CalcReactionProb, DSMC_Chemistry
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars               ,ONLY: offsetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_DSMC_QK_Chemistry       ,ONLY: QK_TestReaction
USE MOD_MCC_XSec                ,ONLY: XSec_CalcReactionProb
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(IN)           :: iElem
LOGICAL, INTENT(INOUT)        :: RelaxToDo
REAL, INTENT(IN), OPTIONAL    :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart1, iPart2, nPartNode, nPair, iCase, ReacTest, iPath, ReacCounter, iSpec
REAL                          :: Volume, NumDens, ReactionProb, iRan, ReactionProbSum
REAL, ALLOCATABLE             :: ReactionProbArray(:)
LOGICAL,ALLOCATABLE           :: PerformReaction(:)
!===================================================================================================================================
iPart1 = Coll_pData(iPair)%iPart_p1
iPart2 = Coll_pData(iPair)%iPart_p2
iCase = CollInf%Coll_Case(PartSpecies(iPart1), PartSpecies(iPart2))
IF (ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths.EQ.0) RETURN
ALLOCATE(PerformReaction(ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths))
! 0.) Get the correct volume and particle number (in the case of an octree the volume and particle number are given as optionals)
IF (PRESENT(NodeVolume)) THEN
  Volume = NodeVolume
ELSE
  Volume = ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
END IF
IF (PRESENT(NodePartNum)) THEN
  nPartNode = NodePartNum
ELSE
  nPartNode = PEM%pNumber(iElem)
END IF
nPair = INT(nPartNode/2)
IF(RadialWeighting%DoRadialWeighting.OR.usevMPF) THEN
  NumDens = SUM(CollInf%Coll_SpecPartNum(:)) / Volume
ELSE IF (UseVarTimeStep) THEN
  NumDens = SUM(CollInf%Coll_SpecPartNum(:)) / Volume * Species(1)%MacroParticleFactor
ELSE
  NumDens = nPartNode / Volume * Species(1)%MacroParticleFactor
END IF
! 1.) Calculate the reaction probabilities/test whether any QK reactions are possible
ReactionProbSum = 0.
ALLOCATE(ReactionProbArray(ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths))
ReactionProbArray = 0.
! Reset the complete array (only populated for the specific collision case)
PerformReaction = .FALSE.
DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
  ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
  IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'QK') THEN
    CALL QK_TestReaction(iPair,ReacTest,PerformReaction(iPath))
  ELSE IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'TCE') THEN
    CALL CalcReactionProb(iPair,ReacTest,ReactionProbArray(iPath),nPair,NumDens)
    ReactionProbSum = ReactionProbSum + ReactionProbArray(iPath)
  END IF
END DO

! 2.) Determine which TCE reaction is most likely to occur
IF(ReactionProbSum.GT.0.) THEN
  ReactionProb = 0.
  CALL RANDOM_NUMBER(iRan)
  ! Check if the reaction probability is greater than a random number
  IF (ReactionProbSum.GT.iRan) THEN
    ! Decide which reaction should occur
    CALL RANDOM_NUMBER(iRan)
    DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
      ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
      IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'TCE') THEN
        ReactionProb = ReactionProb + ReactionProbArray(iPath)
        IF((ReactionProb/ReactionProbSum).GT.iRan) THEN
          PerformReaction(iPath) = .TRUE.
          EXIT
        END IF
      END IF
    END DO
  END IF
END IF

! 3.) Decide which reaction to perform: TCE- and QK-based chemistry
ReactionProb = 0.; ReacCounter = 0
ReacCounter = COUNT(PerformReaction(:))
IF(ReacCounter.GT.0) THEN
  IF(ReacCounter.GT.1) CALL RANDOM_NUMBER(iRan)
  DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
    IF(PerformReaction(iPath)) THEN
      ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
      ! Do not perform a relaxation after the chemical reaction
      RelaxToDo = .FALSE.
      IF(ReacCounter.GT.1) THEN
        ! Determine which reaction will occur, perform it and leave the loop
        ReactionProb = ReactionProb + 1./REAL(ReacCounter)
        IF(ReactionProb.GT.iRan) THEN
          CALL DSMC_Chemistry(iPair, ReacTest)
          ! Exit the routine
          RETURN
        END IF
      ELSE
        CALL DSMC_Chemistry(iPair, ReacTest)
        ! Exit the routine
        RETURN
      END IF
    END IF
  END DO
END IF

! 4.) Cross-section based chemistry (XSec)
IF(ChemReac%CollCaseInfo(iCase)%HasXSecReaction) THEN
  IF(SpecXSec(iCase)%UseCollXSec) THEN
    ! Interpolate the reaction cross-section at the current collision energy
    CALL XSec_CalcReactionProb(iPair,iCase,iElem)
    ReactionProbSum = SpecXSec(iCase)%CrossSection
  ELSE
    ! Reaction probabilities were saved and added to the total collision probability
    ReactionProbSum = Coll_pData(iPair)%Prob
  END IF
  ReactionProb = 0.
  ! Decide which reaction should occur
  CALL RANDOM_NUMBER(iRan)
  DO iPath = 1, ChemReac%CollCaseInfo(iCase)%NumOfReactionPaths
    ReacTest = ChemReac%CollCaseInfo(iCase)%ReactionIndex(iPath)
    IF(TRIM(ChemReac%ReactModel(ReacTest)).EQ.'XSec') THEN
      ReactionProb = ReactionProb + ChemReac%CollCaseInfo(iCase)%ReactionProb(iPath)
      IF((ReactionProb/ReactionProbSum).GT.iRan) THEN
        CALL DSMC_Chemistry(iPair, ReacTest)
        ! Do not perform a relaxation after the chemical reaction
        RelaxToDo = .FALSE.
        ! Exit the routine
        RETURN
      END IF
    END IF
  END DO
  ! Reducing the collision probability (might be used for the determination of the relaxation probability) by the reaction
  ! probability sum if no reaction occurred
  IF(SpecXSec(iCase)%UseCollXSec) THEN
    SpecXSec(iCase)%CrossSection = SpecXSec(iCase)%CrossSection - SUM(ChemReac%CollCaseInfo(iCase)%ReactionProb(:))
  ELSE
    Coll_pData(iPair)%Prob = Coll_pData(iPair)%Prob - SUM(ChemReac%CollCaseInfo(iCase)%ReactionProb(:))
  END IF
END IF

! Species-specific time step: considering the probability that the particle with the lower time step might not have undergone the
! collision, disabling the following collision process. This works with a constant background gas as the change of the background
! gas particle is not tracked anyway.
IF(VarTimeStep%UseSpeciesSpecific) THEN
  IF(VarTimeStep%DisableForMCC) THEN
      iSpec = PartSpecies(iPart1)
      IF(Species(iSpec)%TimeStepFactor.LT.1.) THEN
        CALL RANDOM_NUMBER(iRan)
        IF(iRan.GT.Species(iSpec)%TimeStepFactor) RelaxToDo = .FALSE.
      END IF
  END IF
END IF

END SUBROUTINE ReactionDecision


!--------------------------------------------------------------------------------------------------!
END MODULE MOD_DSMC_Collis
