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

MODULE MOD_DSMC_QK_PROCEDURES
!===================================================================================================================================
! module including qk procedures
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE QK_dissociation
  MODULE PROCEDURE QK_dissociation
END INTERFACE

INTERFACE QK_recombination
  MODULE PROCEDURE QK_recombination
END INTERFACE

INTERFACE QK_exchange
  MODULE PROCEDURE QK_exchange
END INTERFACE

INTERFACE QK_ImpactIonization
  MODULE PROCEDURE QK_ImpactIonization
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! Public Part
PUBLIC :: QK_dissociation, QK_recombination, QK_exchange, QK_ImpactIonization, QK_IonRecombination
!===================================================================================================================================
CONTAINS

SUBROUTINE QK_dissociation(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! Decide whether a dissociation occurs: Check if the relative translation kinetic energy plus vibrational energy of the molecule
! under consideration is larger than the dissociation energy
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY: PartSpecies
USE MOD_DSMC_ChemReact,         ONLY: DSMC_Chemistry
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE                                                                                    !
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
LOGICAL, INTENT(INOUT)        :: RelaxToDo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, PartToExec, PartReac2
!===================================================================================================================================
! Determining, which collision partner is the dissociating particle (always the first molecule in the DefinedReact array)
IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  PartToExec = Coll_pData(iPair)%iPart_p1
  PartReac2  = Coll_pData(iPair)%iPart_p2
ELSE
  PartToExec = Coll_pData(iPair)%iPart_p2
  PartReac2  = Coll_pData(iPair)%iPart_p1
END IF
! Determine the collision energy (relative translational + vibrational energy of dissociating molecule)
Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 + PartStateIntEn(1,PartToExec)
! Correction for second collision partner
IF ((SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.20)) THEN
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - SpecDSMC(PartSpecies(PartReac2))%EZeroPoint
END IF
! Determination of the quantum number corresponding to the collision energy
iQuaMax   = INT(Coll_pData(iPair)%Ec / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - DSMC%GammaQuant)
! Comparing the collision quantum number with the dissociation quantum number
IF ( iQuaMax .GT. SpecDSMC(PartSpecies(PartToExec))%DissQuant ) THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
  IF (.NOT. DSMC%ReservoirSimuRate  ) THEN
# endif
    CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
  END IF
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coefficient
  END IF
# endif
  RelaxToDo = .FALSE.
END IF

END SUBROUTINE QK_dissociation


SUBROUTINE QK_recombination(iPair,iReac,iPart_p3,RelaxToDo,NodeVolume,NodePartNum)
!===================================================================================================================================
! tests for molecular recombination of two colliding atoms by the use of Birds QK theory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY: PartSpecies, Species, PEM, PartState
USE MOD_DSMC_ChemReact,         ONLY: DSMC_Chemistry
USE MOD_Particle_Mesh_Vars,     ONLY: ElemVolume_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARiABLES
INTEGER, INTENT(IN)           :: iPair, iReac,iPart_p3
LOGICAL, INTENT(INOUT)        :: RelaxToDo
REAL, INTENT(IN), OPTIONAL    :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARiABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxColQua, iQua
REAL                          :: Volume, nPartNode, ReactionProb, iRan, Xi, FakXi
LOGICAL                       :: recomb
!===================================================================================================================================

IF (PRESENT(NodeVolume)) THEN
  Volume = NodeVolume
ELSE
  Volume = ElemVolume_Shared(PEM%CNElemID(iPart_p3))
END IF
IF (PRESENT(NodePartNum)) THEN
  nPartNode = NodePartNum
ELSE
  nPartNode = PEM%pNumber(PEM%LocalElemID(iPart_p3))
END IF
! select Q-K Model // do not use Gallis
SELECT CASE (ChemReac%QKMethod(iReac))
  CASE(1) ! Bird and Propability
  ! density of cell
  ! calculate collision temperature
    ReactionProb = ChemReac%QKCoeff(iReac,1) * ( (2 -CollInf%omega(ChemReac%DefinedReact(iReac,1,1)                         &
                 , ChemReac%DefinedReact(iReac,1,2)))**ChemReac%QKCoeff(iReac,2) )                                              &
                 * GAMMA(2-CollInf%omega(ChemReac%DefinedReact(iReac,1,1)                                                   &
                 , ChemReac%DefinedReact(iReac,1,2))) / GAMMA(2-CollInf%omega(ChemReac%DefinedReact(iReac,1,1)              &
                 , ChemReac%DefinedReact(iReac,1,2))+ChemReac%QKCoeff(iReac,2) ) * nPartNode / Volume                           &
                 *  Species(PartSpecies(iPart_p3))%MacroParticleFactor                                                           &
                 *  ( ( ( CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                                   &
                 /              ( 2 * BoltzmannConst * ( 2 - CollInf%omega(ChemReac%DefinedReact(iReac,1,1)                  &
                 , ChemReac%DefinedReact(iReac,1,2)) ) )  )                                                                      &
                 /  SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib)**ChemReac%QKCoeff(iReac,2) )                           &
                 * 1.0/6.0 * PI * ( CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p1))&
                 +           CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p2),PartSpecies(Coll_pData(iPair)%iPart_p2))       &
                 +           CollInf%dref(PartSpecies(iPart_p3),PartSpecies(iPart_p3))       )**3
    IF ( ReactionProb .ge. 1 ) THEN
      IPWRITE(UNIT_stdOut,*) 'ERROR: Recombination probability  >1'
      IPWRITE(UNIT_stdOut,*) 'iReac: ',iReac
      IPWRITE(UNIT_stdOut,*) 'Probability: ', ReactionProb
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT. DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
    END IF
# endif
    CALL RANDOM_NUMBER(iRan)
    IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
    ! calculate collision energy as required to performe the chemical reaction (non-qk)
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + 0.5 * Species(PartSpecies(iPart_p3))%MassIC * &
                              DOTPRODUCT(PartState(4:6,iPart_p3)) &
                             + PartStateIntEn(1,iPart_p3) + PartStateIntEn(2,iPart_p3)
        CALL DSMC_Chemistry(iPair, iReac, iPart_p3)
#if (PP_TimeDiscMethod==42)
      END IF
      IF ( DSMC%ReservoirRateStatistic ) THEN
        ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
      END IF
# endif
      RelaxToDo = .FALSE.
    END IF
!-----------------------------------------------------------------------------------------------------------------------------------
  CASE(2) ! Gallis and trial general LB redistribution
    recomb = .FALSE.
    Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType) * Coll_pData(iPair)%CRela2 + &
                                         ChemReac%EForm(iReac)
    Xi = 2.* ( 2. - CollInf%omega(ChemReac%DefinedReact(iReac,2,1),ChemReac%DefinedReact(iReac,2,1))) &
       + SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
            ! vibrational relaxation
    FakXi = 0.5*Xi -1. ! exponent factor of DOF, substitute of Xi_c - Xi_vib, lax diss page 40
    MaxColQua = INT(Coll_pData(iPair)%Ec/(BoltzmannConst *  SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib) &
                      - DSMC%GammaQuant)
    iQuaMax = MIN( MaxColQua + 1, SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%MaxVibQuant )
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
    END DO
    IF ( iQua .EQ. 0 ) THEN
      ReactionProb = nPartNode * Species(PartSpecies(iPart_p3))%MacroParticleFactor / Volume * 1.0/6.0 * PI         &
                   * (CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p1)) &
                   + CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p2),PartSpecies(Coll_pData(iPair)%iPart_p2))  &
                   + CollInf%dref(PartSpecies(iPart_p3),PartSpecies(iPart_p3)))**3
      CALL RANDOM_NUMBER(iRan)
      IF ( ReactionProb .gt. iRan) THEN
        IF ( MaxColQua .gt. 0 ) THEN
          recomb = .TRUE.
        END IF
      END IF
    END IF
    IF (recomb .eqv. .TRUE. ) THEN
#if (PP_TimeDiscMethod==42)
        ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF (.NOT.DSMC%ReservoirSimuRate ) THEN
# endif
        ! calculate collision energy as required to performe the chemical reaction (non-qk)
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                              + 0.5 * Species(PartSpecies(iPart_p3))%MassIC * &
                              DOTPRODUCT(PartState(4:6,iPart_p3)) &
                              + PartStateIntEn(1,iPart_p3) + PartStateIntEn(2,iPart_p3)
        CALL DSMC_Chemistry(iPair, iReac, iPart_p3)
        RelaxToDo = .FALSE.
#if (PP_TimeDiscMethod==42)
      END IF
      IF ( DSMC%ReservoirRateStatistic ) THEN
        ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
      END IF
# endif
    END IF

    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'ERROR in DSMC_collis: Only two recombination methods of the Q-K method available. ')
    END SELECT

END SUBROUTINE QK_recombination


SUBROUTINE QK_exchange(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! QK_exchange reactions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY : PartSpecies
USE MOD_DSMC_ChemReact,         ONLY : DSMC_Chemistry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
LOGICAL, INTENT(INOUT)        :: RelaxToDo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iQuaMax, MaxColQua, iQua, PartToExec, PartReac2, iQuaDiss, iQua1, iQua2
REAL                          :: ReactionProb, iRan, Xi, FakXi, denominator, coeffT, Tcoll
!===================================================================================================================================

IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  PartToExec = Coll_pData(iPair)%iPart_p1
  PartReac2  = Coll_pData(iPair)%iPart_p2
ELSE
  PartToExec = Coll_pData(iPair)%iPart_p2
  PartReac2  = Coll_pData(iPair)%iPart_p1
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
SELECT CASE (ChemReac%QKMethod(iReac))
!-----------------------------------------------------------------------------------------------------------------------------------
  CASE(1) ! Bird Paper 2011, endothermic exchange reaction
  ! collision energy with respect to the vibrational energy
  ! the energy of the formed vibrational ground state has to substracted from the collision energy
  !           and takes no part
    Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                         + PartStateIntEn(1,PartToExec) - DSMC%GammaQuant * BoltzmannConst            &
                                                          * SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
    iQuaMax = INT(Coll_pData(iPair)%Ec / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - &
                          DSMC%GammaQuant )
    IF ( Coll_pData(iPair)%Ec .gt. ChemReac%EActiv(iReac) ) THEN
      denominator = 0
      DO iQua = 0 , iQuaMax
        denominator = denominator + &
                      (1 - iQua*BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib  / &
                       Coll_pData(iPair)%Ec)**(1-CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec)))
      END DO
      ! normalized ReactionProbe
      ReactionProb =( (  1 - ChemReac%EActiv(iReac) /Coll_pData(iPair)%Ec)&
                        **(1-CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec))) ) / denominator
    ELSE
      ReactionProb = 0.
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT. DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
    END IF
# endif
    CALL RANDOM_NUMBER(iRan)
    IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
    ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
    IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
  ! recalculate collision energy as required for the performance of the exchange reaction
      Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                           + PartStateIntEn(1,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(1,Coll_pData(iPair)%iPart_p2) &
                           + PartStateIntEn(2,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(2,Coll_pData(iPair)%iPart_p2)
      CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
    END IF
    IF ( DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
    END IF
# endif
    RelaxToDo = .FALSE.
    END IF
!-----------------------------------------------------------------------------------------------------------------------------------
  CASE(2) ! Gallis and trial general LB redistribution
    ! only endothermic exchange reaction
    ! this procedure is not valid for exothermic exchange reactions
    Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                        + PartStateIntEn(1,PartToExec)
    iQuaDiss = INT(ChemReac%EActiv(iReac)/(BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib) )
    ! trial GLB redistribution
    Xi = 2.* ( 2. - CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec))) + SpecDSMC(PartSpecies(PartToExec))%Xi_Rot
    FakXi = 0.5*Xi -1. ! exponent factor of DOF, substitute of Xi_c - Xi_vib, Laux1996 diss page 40
    MaxColQua = INT(Coll_pData(iPair)%Ec/(BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib) &
                    - DSMC%GammaQuant)
    iQuaMax = MIN( MaxColQua + 1, SpecDSMC(PartSpecies(PartToExec))%MaxVibQuant )
    CALL RANDOM_NUMBER(iRan)
    iQua = INT(iRan * iQuaMax)
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
      !laux diss page 31
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
    END DO
    ! from thinking should be greather than
    IF ( iQua == iQuadiss) THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
    ! recalculate collision energy as required for the performance of the exchange reaction
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                                  + PartStateIntEn(1,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(1,Coll_pData(iPair)%iPart_p2) &
                                  + PartStateIntEn(2,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(2,Coll_pData(iPair)%iPart_p2)
        CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
      END IF
      IF ( DSMC%ReservoirRateStatistic ) THEN
        ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
      END IF
# endif
    END IF
!-----------------------------------------------------------------------------------------------------------------------------------
  CASE(3) ! Bird 2013, unpublished part of his new book
  ! exothermic exchange reaction
    SELECT CASE (SpecDSMC(PartSpecies(PartReac2))%InterID)
      CASE(2,20,200)
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + PartStateIntEn(1,PartToExec) + PartStateIntEn(2,PartToExec) &
                             + PartStateIntEn(2,PartReac2)  + PartStateIntEn(2,PartReac2)
        iQua1 = INT(PartStateIntEn(1,PartToExec) / ( BoltzmannConst &
              * CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec))) - DSMC%GammaQuant)
        iQua2 = INT(PartStateIntEn(1,PartReac2)  / ( BoltzmannConst &
              * CollInf%omega(PartSpecies(PartReac2),PartSpecies(PartReac2))) - DSMC%GammaQuant)
        coeffT = (2. - CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec)) &
               + 0.5*SpecDSMC(PartSpecies(PartToExec))%Xi_Rot + 0.5*SpecDSMC(PartSpecies(PartReac2 ))%Xi_Rot   + &
                                                iQua1 * log(REAL(1 + 1/iQua1)) + iQua2 * log(real(1 + 1/iQua2) ) )
      CASE DEFAULT
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + PartStateIntEn(1,PartToExec) + PartStateIntEn(2,PartToExec)
        iQua1 = INT(PartStateIntEn(1,PartToExec) / ( BoltzmannConst &
              * CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec))) - DSMC%GammaQuant)
        coeffT = (2. - CollInf%omega(PartSpecies(PartToExec),PartSpecies(PartToExec)) &
               + 0.5*SpecDSMC(PartSpecies(PartToExec))%Xi_Rot + iQua1 * log(real(1 + 1/iQua1) ) )
    END SELECT
    Tcoll = Coll_pData(iPair)%Ec / ( BoltzmannConst * coeffT )
    ! reaction probability
    ReactionProb = ( ChemReac%QKCoeff(iReac,1) * Tcoll ** ChemReac%QKCoeff(iReac,2) ) * &
                         exp( - ChemReac%EActiv(iReac) / ( BoltzmannConst * Tcoll ) )
!            ReactionProb = ChemReac%QKCoeff(iReac,1)                                                             * &
!                           ( ( Coll_pData(iPair)%Ec / ( BoltzmannConst * coeffT ) )**ChemReac%QKCoeff(iReac,2) ) * &
!                           exp( -  ChemReac%EActiv(iReac) * coeffT /  Coll_pData(iPair)%Ec  )
    IF ( ReactionProb .gt. 1 ) THEN
      IPWRITE(UNIT_stdOut,*) 'Error: ReactionProb in chemical reaction >1!'
      IPWRITE(UNIT_stdOut,*) 'ReactionProb:',ReactionProb
      IPWRITE(UNIT_stdOut,*) 'iReac:',iReac
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT. DSMC%ReservoirRateStatistic) THEN
      ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
    END IF
# endif
    CALL RANDOM_NUMBER(iRan)
      IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
        ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
        IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
        ! recalculate collision energy as required for the performance of the exchange reaction
          Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                               + PartStateIntEn(1,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(1,Coll_pData(iPair)%iPart_p2) &
                               + PartStateIntEn(2,Coll_pData(iPair)%iPart_p1) + PartStateIntEn(2,Coll_pData(iPair)%iPart_p2)
          CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
        END IF
        IF ( DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
        END IF
# endif
          RelaxToDo = .FALSE.
      END IF
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'ERROR in DSMC_collis: Only two exchange reaction methods for the Q-K method available. ')
  END SELECT

END SUBROUTINE QK_exchange


SUBROUTINE QK_ImpactIonization(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! Impact ionization by means of QK theory
! derived from the work of Liechty 2010-02
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars             ,ONLY: Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac, DSMC, RadialWeighting
USE MOD_DSMC_ChemReact        ,ONLY: DSMC_Chemistry
USE MOD_Particle_Vars         ,ONLY: PartSpecies, Species, VarTimeStep, usevMPF
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL, INTENT(INOUT)        :: RelaxToDo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: React1Inx, React2Inx, MaxElecQua
REAL                          :: IonizationEnergy, Weight1, Weight2, ReducedMass
!===================================================================================================================================


IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  React1Inx = Coll_pData(iPair)%iPart_p1
  React2Inx = Coll_pData(iPair)%iPart_p2
ELSE
  React1Inx = Coll_pData(iPair)%iPart_p2
  React2Inx = Coll_pData(iPair)%iPart_p1
END IF
! this is based on the idea of the QK method but used accordingly to the dissociation
! this time it is not possible to use quantizied levels as they are not equally spaced
! therefore we use the energy

Weight1 = GetParticleWeight(React1Inx)
Weight2 = GetParticleWeight(React2Inx)

IF (usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarTimeStep%UseVariableTimeStep) THEN
  ReducedMass = (Species(PartSpecies(React1Inx))%MassIC*Weight1 * Species(PartSpecies(React2Inx))%MassIC*Weight2) &
              / (Species(PartSpecies(React1Inx))%MassIC*Weight1 + Species(PartSpecies(React2Inx))%MassIC*Weight2)
ELSE
  ReducedMass = CollInf%MassRed(Coll_pData(iPair)%PairType)
END IF

Coll_pData(iPair)%Ec = 0.5*ReducedMass*Coll_pData(iPair)%CRela2

IF(DSMC%ElectronicModel) Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(3,React1Inx)*Weight1

! ionization level is last known energy level of species
MaxElecQua=SpecDSMC(PartSpecies(React1Inx))%MaxElecQuant - 1
IonizationEnergy=SpecDSMC(PartSpecies(React1Inx))%ElectronicState(2,MaxElecQua)*BoltzmannConst*(2*Weight1 + Weight2)/3.
! if you have electronic levels above the ionization limit, such limits should be used instead of
! the pure energy comparison

IF(Coll_pData(iPair)%Ec .GT. IonizationEnergy)THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate ) THEN
#endif
    ! neu machen
    CALL DSMC_Chemistry(iPair, iReac)
#if (PP_TimeDiscMethod==42)
  END IF
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
  IF ( DSMC%ReservoirRateStatistic) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
#endif
  RelaxToDo = .FALSE.
END IF
END SUBROUTINE


SUBROUTINE QK_IonRecombination(iPair,iReac,iPart_p3,RelaxToDo,NodeVolume,NodePartNum)
!===================================================================================================================================
! Check if colliding ion + electron recombines to neutral atom/ molecule
! requires a third collision partner, which has to take a part of the energy
! this approach is similar to Birds QK molecular recombination, but modified for ion-electron recombination
!===================================================================================================================================
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac!, DSMC_RHS
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars,              ONLY: DSMC
#endif
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_Particle_Vars,          ONLY: PartSpecies, Species, PEM
USE MOD_DSMC_ChemReact,         ONLY: DSMC_Chemistry
USE MOD_Globals_Vars,           ONLY: Pi
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY: ElemVolume_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iPair, iReac,iPart_p3
REAL, INTENT(IN), OPTIONAL          :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL       :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL, INTENT(INOUT)              :: RelaxToDo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iQuaMax1, iQuaMax2,MaxElecQuant, iQua
REAL                                :: Volume, nPartNode, ReactionProb, iRan, Vref, coeffT,Tcoll, acorrection
INTEGER                             :: PartReac1,PartReac2
!===================================================================================================================================


IF (PRESENT(NodeVolume)) THEN
  Volume = NodeVolume
ELSE
  Volume = ElemVolume_Shared(PEM%CNElemID(iPart_p3))
END IF
IF (PRESENT(NodePartNum)) THEN
  nPartNode = NodePartNum
ELSE
  nPartNode = PEM%pNumber(PEM%LocalElemID(iPart_p3))
END IF

IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  PartReac1 = Coll_pData(iPair)%iPart_p1
  PartReac2 = Coll_pData(iPair)%iPart_p2
ELSE
  PartReac1 = Coll_pData(iPair)%iPart_p2
  PartReac2 = Coll_pData(iPair)%iPart_p1
END IF

! Determine max electronic quant of first collision partner
MaxElecQuant = SpecDSMC(PartSpecies(PartReac1))%MaxElecQuant - 1
! determine old Quant
DO iQua = 0, MaxElecQuant
  IF ( PartStateIntEn(3,PartReac1) / BoltzmannConst .ge. &
    SpecDSMC(PartSpecies(PartReac1))%ElectronicState(2,iQua) ) THEN
    iQuaMax1 = iQua
  ELSE
  ! exit loop
    EXIT
  END IF
END DO

! energy of collision
!IF(SpecDSMC(PartSpecies(PartReac2))%InterID.NE.4)THEN
!  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
!                             + PartStateIntEn(3,PartReac1) + PartStateIntEn(3,PartReac2)

!  ! Determine max electronic quant of second collision partner
!  MaxElecQuant = SpecDSMC(PartSpecies(PartReac2))%MaxElecQuant - 1
!  ! determine old Quant
!  DO iQua = 0, MaxElecQuant
!    IF ( PartStateIntEn(3,PartReac2) / BoltzmannConst .ge. &
!      SpecDSMC(PartSpecies(PartReac2))%ElectronicState(2,iQua) ) THEN
!      iQuaMax2 = iQua
!    ELSE
!    ! exit loop
!      EXIT
!    END IF
!  END DO
!ELSE ! second collision partner is an electron
  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + PartStateIntEn(3,PartReac1)
  iQuamax2=0
!END IF

! Bird's recombination approach modified to ion - electron recombination
! reference radius
Vref = 1.0/6.0 * PI * ( CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p1)) + &
                        CollInf%dref(PartSpecies(Coll_pData(iPair)%iPart_p2),PartSpecies(Coll_pData(iPair)%iPart_p2)) + &
                        CollInf%dref(PartSpecies(iPart_p3),PartSpecies(iPart_p3))       )**3

! temperature coeff; therefore,shift iQuaMax, because iQua of first level is zero
coeffT = 2 - CollInf%omega(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)) + (iQuaMax1+1) &
       * LOG(1.0 + 1.0/REAL(iQuaMax1+1) ) + (iQuaMax2+1) * LOG(1.0 + 1.0/REAL(iQuaMax2+1) )

Tcoll = (Coll_pData(iPair)%Ec/(2*BoltzmannConst)) / coeffT

! correction of a
acorrection = (coeffT**ChemReac%QKCoeff(iReac,2))*GAMMA(coeffT)/GAMMA(coeffT+ChemReac%QKCoeff(iReac,2) )

! correct version
! IF(acorrection.LT.0)THEN
!   IPWRITE(UNIT_stdOut,*) 'ERROR: correction term below zero!'
!   IPWRITE(UNIT_stdOut,*) 'iReac: ',iReac
!   IPWRITE(UNIT_stdOut,*) 'acorr: ',acorrection
!   IPWRITE(UNIT_stdOut,*) 'This method can not be employed!'
!   CALL abort(__STAMP__&
!       ,' Error, correction term below zero! ')
! END IF
! simple fix
IF(acorrection.LT.0) acorrection=1.0

! reaction probablility
ReactionProb = nPartNode/Volume * Species(PartSpecies(iPart_p3))%MacroParticleFactor &
             * acorrection* ChemReac%QKCoeff(iReac,1) * Vref *Tcoll**ChemReac%QKCoeff(iReac,2)


!! debug
!print*,'acorrection',acorrection
!print*,'collision energy',Coll_pData(iPair)%Ec
!print*,'a',ChemReac%QKCoeff(iReac,1)
!print*,'b',ChemReac%QKCoeff(iReac,2)
!print*,'vref',vref
!print*,'nPartNode',nPartNode
!print*,'Tcoll',Tcoll
!print*,'Tpower',Tcoll**ChemReac%QKCoeff(iReac,2)
!print*,'Volume',Volume
!print*,'ReactionProb',ReactionProb
!read*

! IF((ReactionProb.GE.1).OR.(ReactionProb.LT.0))THEN
!  IPWRITE(UNIT_stdOut,*) ' ERROR: 1<Recombination probability <0'
!  IPWRITE(UNIT_stdOut,*) ' iReac: ',iReac
!  IPWRITE(UNIT_stdOut,*) ' Probability: ', ReactionProb
! !                STOP
! END IF

#if (PP_TimeDiscMethod==42)
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
  END IF
#endif

CALL RANDOM_NUMBER(iRan)
!ReactionProb=0.0
IF (ReactionProb.GT.iRan) THEN
!evor = 0.5* Species(PartSpecies(PartReac1))%MassIC &
!   * DOTPRODUCT(4:6,PartReac1) &
!    + 0.5* Species(PartSpecies(PartReac2))%MassIC* (PartState(4,PartReac2)**2+PartState(5,PartReac2)**2+PartState(6,PartReac2)**2)&
!    + 0.5* Species(PartSpecies(iPart_p3))%MassIC* (PartState(4,iPart_p3)**2+PartState(5,iPart_p3)**2+PartState(6,iPart_p3)**2) &
!    + PartStateIntEn(3,PartReac1)
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate  ) THEN
#endif

    ! Relative velocity square between mean velocity of pseudo molecule AB and X
!    CRela2X = ((PartState(4,Coll_pData(iPair)%iPart_p1) + PartState(4,Coll_pData(iPair)%iPart_p2))/2 - PartState(4,iPart_p3))**2&
!             +((PartState(5,Coll_pData(iPair)%iPart_p1) + PartState(5,Coll_pData(iPair)%iPart_p2))/2 - PartState(5,iPart_p3))**2&
!             +((PartState(6,Coll_pData(iPair)%iPart_p1) + PartState(6,Coll_pData(iPair)%iPart_p2))/2 - PartState(6,iPart_p3))**2

 ! calculate collision energy as required to performe the chemical reaction (non-qk)

    ! old
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
!      + 0.5 * Species(PartSpecies(iPart_p3))%MassIC * DOTPRODUCT(PartState(4:6,iPart_p3))

    ! new
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
!                         + 0.5*CRela2X * (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC*Species(PartSpecies(iPart_p3))%MassIC) &
!                                        /(Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(PartSpecies(iPart_p3))%MassIC)


!    IF(SpecDSMC(PartSpecies(PartReac1))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,PartReac1) + PartStateIntEn(2,PartReac1)
!    END IF
!    IF(SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,PartReac2) + PartStateIntEn(2,PartReac2)
!    END IF
!    IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(1,iPart_p3) + PartStateIntEn(2,iPart_p3)
!    END IF
!    IF (usevMPF) THEN
!      CALL IonRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
!    ELSE
     !CALL IonRecomb(iReac, iPair, iPart_p3, iElem)
     CALL DSMC_Chemistry(iPair, iReac, iPart_p3)
!    END IF
#if (PP_TimeDiscMethod==42)
  END IF
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
#endif
  RelaxToDo = .FALSE.
!enach = 0.5* Species(PartSpecies(PartReac1))%MassIC &
!  * ((PartState(4,PartReac1)+DSMC_RHS(1,PartReac1))**2+(PartState(5,PartReac1) &
!  +DSMC_RHS(2,PartReac1))**2+(PartState(6,PartReac1)+DSMC_RHS(3,PartReac1))**2) &
!  + 0.5* Species(PartSpecies(iPart_p3))%MassIC* ((PartState(4,iPart_p3)+DSMC_RHS(1,iPart_p3))**2 &
! +(PartState(5,iPart_p3)+DSMC_RHS(2,iPart_p3))**2+(PartState(6,iPart_p3)+DSMC_RHS(3,iPart_p3))**2) &
!  + PartStateIntEn(3,PartReac1)
!print*, evor, enach, evor-enach
!read*
END IF
END SUBROUTINE QK_IonRecombination

END MODULE MOD_DSMC_QK_PROCEDURES
