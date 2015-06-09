#include "boltzplatz.h"

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
PUBLIC :: QK_dissociation, QK_recombination, QK_exchange, QK_ImpactIonization, lacz_gamma
!===================================================================================================================================
CONTAINS

SUBROUTINE QK_dissociation(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! test for QK dissociation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY: PartSpecies, BoltzmannConst
USE MOD_DSMC_ChemReact,         ONLY: MolecDissoc
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
INTEGER                       :: iQuaMax, iQuaDiss, PartToExec, PartReac2
!===================================================================================================================================

IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  PartToExec = Coll_pData(iPair)%iPart_p1
  PartReac2  = Coll_pData(iPair)%iPart_p2
ELSE
  PartToExec = Coll_pData(iPair)%iPart_p2
  PartReac2  = Coll_pData(iPair)%iPart_p1
END IF
! Check if the relative translation kinetic energy plus vibrational energy of the molecule under consideration
! is larger than the activation energy
Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                     + PartStateIntEn(PartToExec,1)
! correcture for second collision partner
IF (SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2) THEN
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - DSMC%GammaQuant *BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib 
END IF

iQuaMax   = INT(Coll_pData(iPair)%Ec   / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - &
                                      DSMC%GammaQuant)
iQuaDiss  = INT(ChemReac%EActiv(iReac) / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) )
IF ( iQuaMax .GT. iQuaDiss ) THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate  ) THEN
# endif
      ! calculate the collision energy as required input for the performance of the dissociation reaction
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,2) + &
                          PartStateIntEn(PartReac2,1) + PartStateIntEn(PartReac2,2)
    CALL MolecDissoc(iReac, iPair)
#if (PP_TimeDiscMethod==42)
  END IF
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
# endif
  RelaxToDo = .FALSE.
END IF
END SUBROUTINE


SUBROUTINE QK_recombination(iPair,iReac,iPart_p3,RelaxToDo,iElem,NodeVolume,NodePartNum)
!===================================================================================================================================
! tests for molecular recombination of two colliding atoms by the use of Birds QK theory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY: PartSpecies, BoltzmannConst, Species, PEM, GEO, PartState,  usevMPF
USE MOD_DSMC_ChemReact,         ONLY: AtomRecomb
USE MOD_vmpf_collision,         ONLY: AtomRecomb_vMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARiABLES
INTEGER, INTENT(IN)           :: iPair, iReac,iPart_p3
LOGICAL, INTENT(INOUT)        :: RelaxToDo
INTEGER, INTENT(IN)           :: iElem
REAL, INTENT(IN), OPTIONAL    :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARiABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                       :: iQuaMax, MaxColQua, iQua
REAL                          :: Volume, nPartNode, omegaAB, ReactionProb, iRan, Xi, FakXi
LOGICAL                       :: recomb
!===================================================================================================================================

IF (PRESENT(NodeVolume)) THEN
  Volume = NodeVolume
ELSE
  Volume = GEO%Volume(PEM%Element(iPart_p3))
END IF
IF (PRESENT(NodePartNum)) THEN
  nPartNode = NodePartNum
ELSE
  nPartNode = PEM%pNumber(PEM%Element(iPart_p3))
END IF
! select Q-K Model // do not use Gallis
SELECT CASE (ChemReac%QKMethod(iReac))
  CASE(1) ! Bird and Propability
  ! density of cell
  ! calculate collision temperature
    omegaAB = 0.5 * ( SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS          &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%omegaVHS  )
    ReactionProb = ChemReac%QKCoeff(iReac,1) * ( (2 -omegaAB)**ChemReac%QKCoeff(iReac,2) )*lacz_gamma(2-omegaAB) / &
                   lacz_gamma(2-omegaAB+ChemReac%QKCoeff(iReac,2) ) * nPartNode / Volume                         * &
                   Species(PartSpecies(iPart_p3))%MacroParticleFactor                                            * &
                   ( ( ( CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                    / &
                               ( 2 * BoltzmannConst * ( 2 - omegaAB ) )  )                                       / &
                             SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib)**ChemReac%QKCoeff(iReac,2) )  * &
                            1.0/6.0 * PI * ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%DrefVHS           + &
                                      SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%DrefVHS                  + &
                                      SpecDSMC(PartSpecies(iPart_p3))%DrefVHS       )**3
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
                                 ( PartState(iPart_p3,4)**2 + PartState(iPart_p3,5)**2 + PartState(iPart_p3,6)**2 ) &
                             + PartStateIntEn(iPart_p3,1) + PartStateIntEn(iPart_p3,2)
        IF (usevMPF) THEN
          CALL AtomRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
        ELSE
          CALL AtomRecomb(iReac, iPair, iPart_p3)
        END IF
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
    Xi = 2.* ( 2. - SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%omegaVHS ) + SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%Xi_Rot
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
      ReactionProb = nPartNode * Species(PartSpecies(iPart_p3))%MacroParticleFactor / Volume     * &
                     1.0/6.0 * PI * ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%DrefVHS + &
                                      SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%DrefVHS + &
                                      SpecDSMC(PartSpecies(iPart_p3))%DrefVHS       )**3
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
                              ( PartState(iPart_p3,4)**2 + PartState(iPart_p3,5)**2 + PartState(iPart_p3,6)**2 ) &
                              + PartStateIntEn(iPart_p3,1) + PartStateIntEn(iPart_p3,2)
        IF (usevMPF) THEN
          CALL AtomRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
        ELSE
          CALL AtomRecomb(iReac, iPair, iPart_p3)
        END IF
        RelaxToDo = .FALSE.
#if (PP_TimeDiscMethod==42)
      END IF
      IF ( DSMC%ReservoirRateStatistic ) THEN
        ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
      END IF
# endif
    END IF

    CASE DEFAULT
      CALL abort(__STAMP__&
          ,'ERROR in DSMC_collis: Only two recombination methods of the Q-K method available. ')
    END SELECT

END SUBROUTINE QK_recombination


SUBROUTINE QK_exchange(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! QK_exchange reactions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst
USE MOD_DSMC_ChemReact,         ONLY : MolecExch
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
                         + PartStateIntEn(PartToExec,1) - DSMC%GammaQuant * BoltzmannConst            &
                                                          * SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
    iQuaMax = INT(Coll_pData(iPair)%Ec / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - &
                          DSMC%GammaQuant )
    IF ( Coll_pData(iPair)%Ec .gt. ChemReac%EActiv(iReac) ) THEN
      denominator = 0
      DO iQua = 0 , iQuaMax
        denominator = denominator + &
                      (1 - iQua*BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib  / &
                       Coll_pData(iPair)%Ec)**(1-SpecDSMC(PartSpecies(PartToExec))%omegaVHS)
      END DO
      ! normalized ReactionProbe
      ReactionProb =( (  1 - ChemReac%EActiv(iReac) /Coll_pData(iPair)%Ec)&
                        **(1-SpecDSMC(PartSpecies(PartToExec))%omegaVHS) ) / denominator
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
                           + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                           + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
      CALL MolecExch(iReac, iPair)
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
                        + PartStateIntEn(PartToExec,1)
    iQuaDiss = INT(ChemReac%EActiv(iReac)/(BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib) )
    ! trial GLB redistribution
    Xi = 2.* ( 2. - SpecDSMC(PartSpecies(PartToExec))%omegaVHS ) + SpecDSMC(PartSpecies(PartToExec))%Xi_Rot
    FakXi = 0.5*Xi -1. ! exponent factor of DOF, substitute of Xi_c - Xi_vib, lax diss page 40
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
                                  + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                                  + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
        CALL MolecExch(iReac, iPair)
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
                             + PartStateIntEn(PartToExec,1) + PartStateIntEn(PartToExec,2) &
                             + PartStateIntEn(PartReac2,2)  + PartStateIntEn(PartReac2,2)
        iQua1 = INT(PartStateIntEn(PartToExec,1) / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%omegaVHS ) - &
                      DSMC%GammaQuant)
        iQua2 = INT(PartStateIntEn(PartReac2,1)  / ( BoltzmannConst * SpecDSMC(PartSpecies(PartReac2 ))%omegaVHS ) - &
                      DSMC%GammaQuant)
        coeffT = (2. - SpecDSMC(PartSpecies(PartToExec))%omegaVHS + 0.5*SpecDSMC(PartSpecies(PartToExec))%Xi_Rot + &
                                                                  0.5*SpecDSMC(PartSpecies(PartReac2 ))%Xi_Rot   + &
                                                iQua1 * log(REAL(1 + 1/iQua1)) + iQua2 * log(real(1 + 1/iQua2) ) )
      CASE DEFAULT
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + PartStateIntEn(PartToExec,1) + PartStateIntEn(PartToExec,2)
        iQua1 = INT(PartStateIntEn(PartToExec,1) / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%omegaVHS ) - &
                       DSMC%GammaQuant)
        coeffT = (2. - SpecDSMC(PartSpecies(PartToExec))%omegaVHS + 0.5*SpecDSMC(PartSpecies(PartToExec))%Xi_Rot + &
                                      iQua1 * log(real(1 + 1/iQua1) ) )
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
                               + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                               + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
          CALL MolecExch(iReac, iPair)
#if (PP_TimeDiscMethod==42)
        END IF
        IF ( DSMC%ReservoirRateStatistic) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
        END IF
# endif
          RelaxToDo = .FALSE.
      END IF
    CASE DEFAULT
      CALL abort(__STAMP__&
          ,'ERROR in DSMC_collis: Only two exchange reaction methods for the Q-K method available. ')
  END SELECT

END SUBROUTINE QK_exchange


SUBROUTINE QK_ImpactIonization(iPair,iReac,RelaxToDo)
!===================================================================================================================================
! Impact ionization by means of QK theory
! derived from the work of Liechty 2010-02
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst
USE MOD_DSMC_ChemReact,         ONLY : ElecImpactIoni
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars,              ONLY : DSMC
#endif
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
INTEGER                       :: PartToExec, PartReac2
REAL                          :: JToEv
!===================================================================================================================================

JToEv = 1.602176565E-19 

IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
  PartToExec = Coll_pData(iPair)%iPart_p1
  PartReac2 = Coll_pData(iPair)%iPart_p2
ELSE
  PartToExec = Coll_pData(iPair)%iPart_p2
  PartReac2 = Coll_pData(iPair)%iPart_p1
END IF
! this is based on the idea of the QK method but used accordingly to the dissociation
! this time it is not possible to use quantizied levels as they are not equally spaced
! therefore we use the energy
Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                           + PartStateIntEn(PartToExec,3)
! if you have electronic levels above the ionization limit, such limits should be used instead of
! the pure energy comparison
IF ( Coll_pData(iPair)%Ec .gt. SpecDSMC(PartSpecies(PartToExec))%Eion_eV * JToEv ) THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
    ! neu machen
    CALL ElecImpactIoni(iReac, iPair)
#if (PP_TimeDiscMethod==42)
  END IF
  IF ( DSMC%ReservoirRateStatistic) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
# endif
  RelaxToDo = .FALSE.
END IF
END SUBROUTINE

!-----------------------------------------------------------------------------------------------------------------------------------

RECURSIVE FUNCTION lacz_gamma(a) result(g)
!===================================================================================================================================
! gamma function taken from
! http://rosettacode.org/wiki/Gamma_function#Fortran
! variefied against build in and compiled with double precision
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    !
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(KIND=8), INTENT(IN) :: a
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER, PARAMETER :: cg = 7
REAL(KIND=8) :: g
! these precomputed values are taken by the sample code in Wikipedia,
! and the sample itself takes them from the GNU Scientific Library
REAL(KIND=8), DIMENSION(0:8), PARAMETER :: p = &
     (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
     771.32342877765313, -176.61502916214059, 12.507343278686905, &
     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)

REAL(KIND=8) :: t, w, x
INTEGER :: icount
!===================================================================================================================================

x = a

IF ( x < 0.5 ) THEN
   g = PI / ( SIN(pi*x) * lacz_gamma(1.0-x) )
ELSE
   x = x - 1.0
   t = p(0)
   DO icount=1, cg+2
      t = t + p(icount)/(x+REAL(icount))
   END DO
   w = x + REAL(cg) + 0.5
   g = SQRT(2.0*PI) * w**(x+0.5) * EXP(-w) * t
END IF
END FUNCTION lacz_gamma

END MODULE MOD_DSMC_QK_PROCEDURES
