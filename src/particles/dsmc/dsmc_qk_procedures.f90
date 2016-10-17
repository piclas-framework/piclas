#include "boltzplatz.h"

MODULE MOD_DSMC_QK_PROCEDURES
!===================================================================================================================================
! module including qk procedures
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE QK_dissociation ! this routine has not been varified yet
  MODULE PROCEDURE QK_dissociation
END INTERFACE

INTERFACE QK_recombination ! this routine has not been varified yet
  MODULE PROCEDURE QK_recombination
END INTERFACE

INTERFACE QK_exchange ! this routine has not been varified yet
  MODULE PROCEDURE QK_exchange
END INTERFACE

INTERFACE QK_ImpactIonization
  MODULE PROCEDURE QK_ImpactIonization
END INTERFACE


INTERFACE QK_IonRecombination ! this routine has not been varified yet
  MODULE PROCEDURE QK_IonRecombination
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part
! Public Part 
PUBLIC :: QK_dissociation, QK_recombination, QK_exchange, QK_ImpactIonization, lacz_gamma, QK_IonRecombination
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
! correction for second collision partner
IF (SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2) THEN
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - DSMC%GammaQuant *BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib 
END IF

iQuaMax   = INT(Coll_pData(iPair)%Ec   / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ) - &
                                      DSMC%GammaQuant)
iQuaDiss  = INT(ChemReac%EActiv(iReac) / ( BoltzmannConst * SpecDSMC(PartSpecies(PartToExec))%CharaTVib ))

IF (SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2) THEN ! DEBUG
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + DSMC%GammaQuant *BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib 
END IF

IF ( iQuaMax .GT. iQuaDiss ) THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to perform the reaction
  IF (.NOT. DSMC%ReservoirSimuRate  ) THEN
# endif
      ! calculate the collision energy as required input for the performance of the dissociation reaction
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,2) + &
                          PartStateIntEn(PartReac2,1) + PartStateIntEn(PartReac2,2)
    CALL DSMC_Chemistry(iReac, iPair)
#if (PP_TimeDiscMethod==42)
  END IF
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1 ! DEBUG
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coefficient
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
USE MOD_Particle_Vars,          ONLY: PartSpecies, BoltzmannConst, Species, PEM, PartState,  usevMPF
USE MOD_Particle_Mesh_Vars,     ONLY: GEO
USE MOD_DSMC_ChemReact,         ONLY: DSMC_Chemistry
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
          CALL DSMC_Chemistry(iReac, iPair, iPart_p3)
        END IF
#if (PP_TimeDiscMethod==42)
      END IF
      ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1 ! DEBUG
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
          CALL DSMC_Chemistry(iReac, iPair, iPart_p3)
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
USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst
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
    ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1 ! DEBUG
    IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
  ! recalculate collision energy as required for the performance of the exchange reaction
      Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                           + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                           + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
      CALL DSMC_Chemistry(iReac, iPair)
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
        CALL DSMC_Chemistry(iReac, iPair)
#if (PP_TimeDiscMethod==42)
      END IF
      ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1 ! DEBUG
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
    ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1 ! DEBUG
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
          CALL DSMC_Chemistry(iReac, iPair)
#if (PP_TimeDiscMethod==42)
        END IF
        ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
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
USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac
USE MOD_Particle_Vars,          ONLY : PartSpecies, BoltzmannConst
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
INTEGER                       :: PartToExec, PartReac2, MaxElecQua
REAL                          :: IonizationEnergy
!===================================================================================================================================


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
! ionization level is last known energy level of species
MaxElecQua=SpecDSMC(PartSpecies(PartToExec))%MaxElecQuant - 1
IonizationEnergy=SpecDSMC(PartSpecies(PartToExec))%ElectronicState(2,MaxElecQua)*BoltzmannConst
! if you have electronic levels above the ionization limit, such limits should be used instead of
! the pure energy comparison
IF(Coll_pData(iPair)%Ec .GT. IonizationEnergy)THEN
#if (PP_TimeDiscMethod==42)
  ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate ) THEN
# endif
    ! neu machen
    CALL ElecImpactIoniQK(iReac, iPair)
#if (PP_TimeDiscMethod==42)
  END IF
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
  IF ( DSMC%ReservoirRateStatistic) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
# endif
  RelaxToDo = .FALSE.
END IF
END SUBROUTINE


SUBROUTINE ElecImpactIoniQK(iReac, iPair)
!===================================================================================================================================
! Perfoms the electron impact ionization with Q-K
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,             ONLY : Coll_pData, DSMC_RHS, CollInf, SpecDSMC, DSMCSumOfFormedParticles
USE MOD_DSMC_Vars,             ONLY : ChemReac, PartStateIntEn !, Debug_Energy
USE MOD_Particle_Vars,         ONLY : BoltzmannConst, PartSpecies, PartState, PDM, PEM, NumRanVec, RandomVec
USE MOD_Particle_Vars,         ONLY : usevMPF, PartPosRef
USE MOD_DSMC_ElectronicModel,  ONLY : ElectronicEnergyExchange, TVEEnergyExchange
USE MOD_DSMC_Relaxation,       ONLY : DSMC_VibRelaxDiatomic
USE MOD_DSMC_PolyAtomicModel,  ONLY : DSMC_VibRelaxPoly, DSMC_RotRelaxPoly, DSMC_RelaxVibPolyProduct
USE MOD_Particle_Tracking_vars, ONLY: DoRefMapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
  INTEGER, INTENT(IN)           :: iPair, iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  INTEGER                       :: iVec
  REAL                          :: iRan, FacEtraDistri
  REAL                          :: ERel_React1_React2, ERel_React2_Elec
  INTEGER                       :: PositionNbr, React1Inx, ElecInx
  REAL                          :: VxPseuAtom, VyPseuAtom, VzPseuAtom
  INTEGER                       :: MaxElecQua
  REAL                          :: IonizationEnergy
  REAL                          :: FakXi, Xi!, Debug_Energy(2)
!===================================================================================================================================


!..Get the index of react1 and the electron
IF (SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.eq.4) THEN
  ElecInx = Coll_pData(iPair)%iPart_p1
  React1Inx = Coll_pData(iPair)%iPart_p2 
ELSE
  React1Inx = Coll_pData(iPair)%iPart_p1
  ElecInx = Coll_pData(iPair)%iPart_p2
END IF

!Debug_Energy=0.0
!Debug_Energy(1)  = Debug_Energy(1) +&
!      0.5* Species(PartSpecies(React1Inx))%MassIC*(PartState(React1Inx,4)**2+PartState(React1Inx,5)**2+PartState(React1Inx,6)**2)&
!    + 0.5* Species(PartSpecies(ElecInx))%MassIC*(PartState(ElecInx,4)**2+PartState(ElecInx,5)**2+PartState(ElecInx,6)**2)&
!    + PartStateIntEn(React1Inx,3)


! ionization level is last known energy level of species
MaxElecQua=SpecDSMC(PartSpecies(React1Inx))%MaxElecQuant - 1
IonizationEnergy=SpecDSMC(PartSpecies(React1Inx))%ElectronicState(2,MaxElecQua)*BoltzmannConst
! nullify
!ElecTransfer = 0.

IF(usevMPF) CALL abort(&
       __STAMP__&
        ,' Reaction not implemented with vMPF ',iReac)

IF (SpecDSMC(PartSpecies(ElecInx))%InterID.NE.4) CALL abort(&
       __STAMP__&
        ,' Only electron impact ionization. Further collision partner not implemented!  ',iReac)


! remove ionization energy from collision
Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                           + PartStateIntEn(React1Inx,3)- IonizationEnergy

Xi  = 4.*(2. - SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS) 

! if first collison partner is a molecule add its internal energy to the collision energy
PartSpecies(React1Inx) = ChemReac%DefinedReact(iReac,2,1)
IF((SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(React1Inx))%InterID.EQ.20)) THEN  
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(React1Inx,2)
  Xi = Xi + SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
  FakXi = 0.5*Xi  - 1 
  CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi )
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,3) + PartStateIntEn(React1Inx,1) 
  IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN      
    CALL DSMC_VibRelaxPoly(iPair, React1Inx,FakXi)
  ELSE
    CALL DSMC_VibRelaxDiatomic(iPair,React1Inx,FakXi)
  END IF
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,1) 
  IF(SpecDSMC(PartSpecies(React1Inx))%PolyatomicMol) THEN 
    FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(React1Inx))%Xi_Rot
    CALL DSMC_RotRelaxPoly(iPair, React1Inx, FakXi)
  ELSE
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(React1Inx,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
  END IF
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,2) 
ELSE
  ! store the electronic energy of the ionized particle
  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40  
  CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,React1Inx,FakXi)
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(React1Inx,3) 
  !ElecTransfer = PartStateIntEn(React1Inx,3)
END IF

! distribute Etra of pseudo neutral particle (i+e) and old electron
CALL RANDOM_NUMBER(iRan)
FacEtraDistri = iRan
CALL RANDOM_NUMBER(iRan)
! laux diss page 40, omegaVHS only of one species
DO WHILE ((4 *FacEtraDistri*(1-FacEtraDistri))**(1-SpecDSMC(PartSpecies(React1Inx))%omegaVHS).LT.iRan)
  CALL RANDOM_NUMBER(iRan)
  FacEtraDistri = iRan
  CALL RANDOM_NUMBER(iRan)
END DO
ERel_React1_React2 = Coll_pData(iPair)%Ec * FacEtraDistri
ERel_React2_Elec = Coll_pData(iPair)%Ec - ERel_React1_React2

!.... Get free particle index for the 3rd particle produced
DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
IF (PositionNbr.EQ.0) THEN
   CALL Abort(&
    __STAMP__,&
    ' New Particle Number greater than max particle number!')
END IF

!Set new Species of electron
PDM%ParticleInside(PositionNbr) = .true.
PartSpecies(PositionNbr) = ChemReac%DefinedReact(iReac,2,3)
PartState(PositionNbr,1:3) = PartState(React1Inx,1:3)
IF(DoRefMapping)THEN ! copy particle position in reference state to new particle
  PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,React1Inx)
END IF

PartStateIntEn(PositionNbr, 1) = 0.
PartStateIntEn(PositionNbr, 2) = 0.
PartStateIntEn(PositionNbr, 3) = 0.
PEM%Element(PositionNbr) = PEM%Element(React1Inx)

!Scattering of pseudo atom (e-i) and collision partner e (scattering of e) 
FracMassCent1 = CollInf%FracMassCent(ChemReac%DefinedReact(iReac,1,1), &
    CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),PartSpecies(ElecInx)))
FracMassCent2 = CollInf%FracMassCent(PartSpecies(ElecInx), &
    CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),PartSpecies(ElecInx)))

!Calculation of velo from center of mass
VeloMx = FracMassCent1 * PartState(React1Inx, 4) &
       + FracMassCent2 * PartState(ElecInx, 4)
VeloMy = FracMassCent1 * PartState(React1Inx, 5) &
       + FracMassCent2 * PartState(ElecInx, 5)
VeloMz = FracMassCent1 * PartState(React1Inx, 6) &
       + FracMassCent2 * PartState(ElecInx, 6)

!calculate random vec and new squared velocities
Coll_pData(iPair)%CRela2 = 2 * ERel_React1_React2 /CollInf%MassRed(Coll_pData(iPair)%PairType)
CALL RANDOM_NUMBER(iRan)
iVec = INT(NumRanVec * iRan + 1)
RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

! deltaV particle 2
DSMC_RHS(ElecInx,1) = VeloMx - FracMassCent1*RanVelox - PartState(ElecInx, 4)
DSMC_RHS(ElecInx,2) = VeloMy - FracMassCent1*RanVeloy - PartState(ElecInx, 5)
DSMC_RHS(ElecInx,3) = VeloMz - FracMassCent1*RanVeloz - PartState(ElecInx, 6)

!Debug_Energy(2) = Debug_Energy(2)&
!      + 0.5* Species(PartSpecies(ElecInx))%MassIC * (&
!           (VeloMx - FracMassCent1*RanVelox)**2 &
!          +(VeloMy - FracMassCent1*RanVeloy)**2 &
!          +(VeloMz - FracMassCent1*RanVeloz)**2) 

!Set velocity of pseudo atom (i+e)
VxPseuAtom = (VeloMx + FracMassCent2*RanVelox)
VyPseuAtom = (VeloMy + FracMassCent2*RanVeloy)
VzPseuAtom = (VeloMz + FracMassCent2*RanVeloz)

!Scattering of i + e
FracMassCent1 = CollInf%FracMassCent(PartSpecies(React1Inx), &
              &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
FracMassCent2 = CollInf%FracMassCent(PartSpecies(PositionNbr), & 
              &  CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))

!calculate random vec and new squared velocities
Coll_pData(iPair)%CRela2 = 2 *  ERel_React2_Elec / & 
        CollInf%MassRed(CollInf%Coll_Case(PartSpecies(React1Inx),PartSpecies(PositionNbr)))
CALL RANDOM_NUMBER(iRan)
iVec = INT(NumRanVec * iRan + 1)
RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

!deltaV particle 1
DSMC_RHS(React1Inx,1) = VxPseuAtom + FracMassCent2*RanVelox - PartState(React1Inx, 4)
DSMC_RHS(React1Inx,2) = VyPseuAtom + FracMassCent2*RanVeloy - PartState(React1Inx, 5)
DSMC_RHS(React1Inx,3) = VzPseuAtom + FracMassCent2*RanVeloz - PartState(React1Inx, 6)

!deltaV new formed particle
PartState(PositionNbr,4:6) = 0.
DSMC_RHS(PositionNbr,1) = VxPseuAtom - FracMassCent1*RanVelox 
DSMC_RHS(PositionNbr,2) = VyPseuAtom - FracMassCent1*RanVeloy 
DSMC_RHS(PositionNbr,3) = VzPseuAtom - FracMassCent1*RanVeloz 

!Debug_Energy(2) = Debug_Energy(2)&
!      + 0.5* Species(PartSpecies(React1Inx))%MassIC  * (&
!           (VxPseuAtom + FracMassCent2*RanVelox)**2 &
!          +(VyPseuAtom + FracMassCent2*RanVeloy)**2 &
!          +(VzPseuAtom + FracMassCent2*RanVeloz)**2) &
!      + 0.5* Species(PartSpecies(PositionNbr))%MassIC  * (&
!           (VxPseuAtom - FracMassCent1*RanVelox)**2 &
!          +(VyPseuAtom - FracMassCent1*RanVeloy)**2 &
!          +(VzPseuAtom - FracMassCent1*RanVeloz)**2) &
!        + PartStateIntEn(React1Inx,3) + IonizationEnergy

END SUBROUTINE ElecImpactIoniQK


SUBROUTINE QK_IonRecombination(iPair,iReac,iPart_p3,RelaxToDo,iElem,NodeVolume,NodePartNum)
!===================================================================================================================================
! Check if colliding ion + electron recombines to neutral atom/ molecule
! requires a third collision partner, which has to take a part of the energy
! this approach is similar to Birds QK molecular recombination, but modified for ion-electron recombination
!===================================================================================================================================
USE MOD_DSMC_Vars,              ONLY: Coll_pData, CollInf, SpecDSMC, PartStateIntEn, ChemReac!, DSMC_RHS
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars,              ONLY: DSMC
#endif
USE MOD_Particle_Vars,          ONLY: PartSpecies, BoltzmannConst, Species, PEM
USE MOD_Particle_Mesh_Vars,     ONLY: GEO
USE MOD_DSMC_ChemReact,         ONLY: DSMC_Chemistry
USE MOD_Globals_Vars,           ONLY: Pi
USE MOD_Globals
!USE MOD_vmpf_collision,         ONLY: IonRecomb_vMPF
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
INTEGER, INTENT(IN)                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: iQuaMax1, iQuaMax2,MaxElecQuant, iQua
REAL                                :: Volume, nPartNode, omegaAB, ReactionProb, iRan, Vref, coeffT,Tcoll, acorrection!, CRela2X
INTEGER                             :: PartReac1,PartReac2
!REAL :: evor, enach
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
  IF ( PartStateIntEn(PartReac1,3) / BoltzmannConst .ge. &
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
!                             + PartStateIntEn(PartReac1,3) + PartStateIntEn(PartReac2,3)

!  ! Determine max electronic quant of second collision partner
!  MaxElecQuant = SpecDSMC(PartSpecies(PartReac2))%MaxElecQuant - 1
!  ! determine old Quant
!  DO iQua = 0, MaxElecQuant
!    IF ( PartStateIntEn(PartReac2,3) / BoltzmannConst .ge. &
!      SpecDSMC(PartSpecies(PartReac2))%ElectronicState(2,iQua) ) THEN
!      iQuaMax2 = iQua
!    ELSE
!    ! exit loop
!      EXIT
!    END IF
!  END DO
!ELSE ! second collision partner is an electron
  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                             + PartStateIntEn(PartReac1,3)
  iQuamax2=0
!END IF

! Bird's recombination approach modified to ion - electron recombination
! reference radius
Vref = 1.0/6.0 * PI * ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%DrefVHS + &
                        SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%DrefVHS + &
                        SpecDSMC(PartSpecies(iPart_p3))%DrefVHS       )**3

! omega VHS
omegaAB = 0.5 * ( SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS          &
                + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%omegaVHS  )

! temperature coeff; therefore,shift iQuaMax, because iQua of first level is zero
coeffT = 2 - omegaAB + (iQuaMax1+1) * LOG(1.0 + 1.0/REAL(iQuaMax1+1) ) + (iQuaMax2+1) * LOG(1.0 + 1.0/REAL(iQuaMax2+1) ) 

Tcoll = (Coll_pData(iPair)%Ec/(2*BoltzmannConst)) / coeffT

! correction of a
!acorrection = (coeffT**ChemReac%QKCoeff(iReac,2))*lacz_gamma(coeffT)/lacz_gamma(coeffT+ChemReac%QKCoeff(iReac,2) ) 
acorrection = (coeffT**ChemReac%QKCoeff(iReac,2))*gamma(coeffT)/gamma(coeffT+ChemReac%QKCoeff(iReac,2) ) 

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
# endif

CALL RANDOM_NUMBER(iRan)
!ReactionProb=0.0
IF (ReactionProb.GT.iRan) THEN
!evor = 0.5* Species(PartSpecies(PartReac1))%MassIC &
!   * (PartState(PartReac1,4)**2+PartState(PartReac1,5)**2+PartState(PartReac1,6)**2) &
!    + 0.5* Species(PartSpecies(PartReac2))%MassIC* (PartState(PartReac2,4)**2+PartState(PartReac2,5)**2+PartState(PartReac2,6)**2)&
!    + 0.5* Species(PartSpecies(iPart_p3))%MassIC* (PartState(iPart_p3,4)**2+PartState(iPart_p3,5)**2+PartState(iPart_p3,6)**2) &
!    + PartStateIntEn(PartReac1,3)
#if (PP_TimeDiscMethod==42)
! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
  IF (.NOT. DSMC%ReservoirSimuRate  ) THEN
# endif

    ! Relative velocity square between mean velocity of pseudo molecule AB and X
!    CRela2X = ((PartState(Coll_pData(iPair)%iPart_p1,4) + PartState(Coll_pData(iPair)%iPart_p2,4))/2 - PartState(iPart_p3,4))**2&
!             +((PartState(Coll_pData(iPair)%iPart_p1,5) + PartState(Coll_pData(iPair)%iPart_p2,5))/2 - PartState(iPart_p3,5))**2&
!             +((PartState(Coll_pData(iPair)%iPart_p1,6) + PartState(Coll_pData(iPair)%iPart_p2,6))/2 - PartState(iPart_p3,6))**2

 ! calculate collision energy as required to performe the chemical reaction (non-qk)

    ! old
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
!      + 0.5 * Species(PartSpecies(iPart_p3))%MassIC * ( PartState(iPart_p3,4)**2 + PartState(iPart_p3,5)**2 &
!                                                                                 + PartState(iPart_p3,6)**2 )

    ! new
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec &
!                         + 0.5*CRela2X * (Species(ChemReac%DefinedReact(iReac,2,1))%MassIC*Species(PartSpecies(iPart_p3))%MassIC) &
!                                        /(Species(ChemReac%DefinedReact(iReac,2,1))%MassIC+Species(PartSpecies(iPart_p3))%MassIC)


!    IF(SpecDSMC(PartSpecies(PartReac1))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(PartReac1,1) + PartStateIntEn(PartReac1,2)
!    END IF
!    IF(SpecDSMC(PartSpecies(PartReac2))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(PartReac2,1) + PartStateIntEn(PartReac2,2)
!    END IF
!    IF(SpecDSMC(PartSpecies(iPart_p3))%InterID.EQ.2)THEN
!      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(iPart_p3,1) + PartStateIntEn(iPart_p3,2)
!    END IF
!    IF (usevMPF) THEN
!      CALL IonRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
!    ELSE
     !CALL IonRecomb(iReac, iPair, iPart_p3, iElem)
     CALL DSMC_Chemistry(iReac, iPair, iPart_p3)
!    END IF
#if (PP_TimeDiscMethod==42)
  END IF
  ChemReac%ReacCount(iReac) = ChemReac%ReacCount(iReac) + 1
  IF ( DSMC%ReservoirRateStatistic ) THEN
    ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
  END IF
# endif
  RelaxToDo = .FALSE.
!enach = 0.5* Species(PartSpecies(PartReac1))%MassIC &
!  * ((PartState(PartReac1,4)+DSMC_RHS(PartReac1,1))**2+(PartState(PartReac1,5) &
!  +DSMC_RHS(PartReac1,2))**2+(PartState(PartReac1,6)+DSMC_RHS(PartReac1,3))**2) &
!  + 0.5* Species(PartSpecies(iPart_p3))%MassIC* ((PartState(iPart_p3,4)+DSMC_RHS(iPart_p3,1))**2 &
! +(PartState(iPart_p3,5)+DSMC_RHS(iPart_p3,2))**2+(PartState(iPart_p3,6)+DSMC_RHS(iPart_p3,3))**2) &
!  + PartStateIntEn(PartReac1,3)
!print*, evor, enach, evor-enach
!read*
END IF
END SUBROUTINE QK_IonRecombination


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
   !DO icount=1, cg+2
   DO icount=1, cg+1
      t = t + p(icount)/(x+REAL(icount))
   END DO
   w = x + REAL(cg) + 0.5
   g = SQRT(2.0*PI) * w**(x+0.5) * EXP(-w) * t
END IF
END FUNCTION lacz_gamma

END MODULE MOD_DSMC_QK_PROCEDURES
