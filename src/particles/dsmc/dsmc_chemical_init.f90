#include "boltzplatz.h"

MODULE MOD_DSMC_ChemInit
!===================================================================================================================================
! Initialization of chemical module
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_chemical_init
  MODULE PROCEDURE DSMC_chemical_init
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_chemical_init, InitPartitionFunction
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_chemical_init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,          ONLY : ChemReac,CollisMode, DSMC
  USE MOD_ReadInTools
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies, BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(LEN=3)         :: hilf 
  INTEGER               :: iReac, iReac2, iReac3, iReac4, iSpec, iChemDir
  INTEGER, ALLOCATABLE  :: PairCombID(:,:), DummyRecomb(:,:)
  LOGICAL, ALLOCATABLE  :: YetDefined_Help(:)
  INTEGER               :: Reactant1, Reactant2, Reactant3, MaxSpecies
!===================================================================================================================================
  
! reading reaction values   
  ChemReac%NumOfReact = GETINT('DSMC-NumOfReactions','0')
  IF(CollisMode.EQ.3)THEN
    IF(ChemReac%NumOfReact.EQ.0)THEN
      CALL Abort(&
          __STAMP__,&
          ' Collisions with chemical reactions require chemical reaction database! ',CollisMode,REAL(ChemReac%NumOfReact))
    END IF
  END IF
  ! Calculation of the backward reaction rates 
  IF(DSMC%BackwardReacRate) THEN
   ChemReac%NumOfReact = 2 * ChemReac%NumOfReact
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (ChemReac%NumOfReact.GT.0) THEN
    ALLOCATE(ChemReac%NumReac(ChemReac%NumOfReact))
    ChemReac%NumReac = 0
#if (PP_TimeDiscMethod==42)
    ALLOCATE(ChemReac%ReacCount(ChemReac%NumOfReact))
    ChemReac%ReacCount = 0
    ALLOCATE(ChemReac%ReacCollMean(ChemReac%NumOfReact))
    ChemReac%ReacCollMean = 0.0
    ALLOCATE(ChemReac%ReacCollMeanCount(ChemReac%NumOfReact))
    ChemReac%ReacCollMeanCount = 0
#endif
    ALLOCATE(ChemReac%QKProcedure(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%QKMethod(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%QKCoeff(ChemReac%NumOfReact,2))
    ALLOCATE(YetDefined_Help(ChemReac%NumOfReact))
    YetDefined_Help = .FALSE.
    ALLOCATE(ChemReac%ReactType(ChemReac%NumOfReact))
    ChemReac%ReactType = '0'
    ALLOCATE(ChemReac%DefinedReact(ChemReac%NumOfReact,2,3))
    ChemReac%DefinedReact = 0
    ALLOCATE(ChemReac%ReactCase(nSpecies,nSpecies))
    ChemReac%ReactCase = 0
    ALLOCATE(ChemReac%ReactNum(nSpecies, nSpecies, nSpecies))
    ChemReac%ReactNum = 0
    ALLOCATE(ChemReac%ReactNumRecomb(nSpecies, nSpecies, nSpecies))
    ChemReac%ReactNumRecomb = 0
    ALLOCATE(ChemReac%Arrhenius_Prefactor(ChemReac%NumOfReact),&
             ChemReac%Arrhenius_Powerfactor(ChemReac%NumOfReact),&
             ChemReac%EActiv(ChemReac%NumOfReact),&
             ChemReac%EForm(ChemReac%NumOfReact),&
             ChemReac%Hab(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MeanEVibQua_PerIter(nSpecies))
    ALLOCATE(ChemReac%MeanEVib_PerIter(nSpecies))
    ChemReac%Hab=0.0
    ALLOCATE(DummyRecomb(nSpecies,nSpecies))
    DummyRecomb = 0
    ALLOCATE(ChemReac%CEXa(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%CEXb(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MEXa(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MEXb(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%ReactInfo(ChemReac%NumOfReact))

    DO iReac = 1, ChemReac%NumOfReact
      WRITE(UNIT=hilf,FMT='(I3)') iReac
      ChemReac%ReactType(iReac)             = GETSTR('DSMC-Reaction'//TRIM(hilf)//'-ReactionType','0')
      ChemReac%QKProcedure(iReac)           = GETLOGICAL('DSMC-Reaction'//TRIM(hilf)//'-QKProcedure','.FALSE.')
      CHemReac%QKMethod(iReac)       = GETINT('DSMC-Reaction'//TRIM(hilf)//'-QK-Method','0') 
      ChemReac%QKCoeff(iReac,1)      = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-QK-Coeff1','0')
      ChemReac%QKCoeff(iReac,2)      = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-QK-Coeff2','0')
      ChemReac%DefinedReact(iReac,1,:)      = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Reactants',3,'0,0,0')
      ChemReac%DefinedReact(iReac,2,:)      = GETINTARRAY('DSMC-Reaction'//TRIM(hilf)//'-Products',3,'0,0,0')
      ChemReac%Arrhenius_Prefactor(iReac)   = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Prefactor','0')
      ChemReac%Arrhenius_Powerfactor(iReac) = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Arrhenius-Powerfactor','0')
      ChemReac%EActiv(iReac)                = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-Activation-Energy_K','0')*BoltzmannConst
      ChemReac%EForm(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-HeatOfFormation_K','0')*BoltzmannConst
      ChemReac%CEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXa','0')
      ChemReac%CEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-CEXb','0')
      ChemReac%MEXa(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXa','0')
      ChemReac%MEXb(iReac)                 = GETREAL('DSMC-Reaction'//TRIM(hilf)//'-MEXb','0')
      ! Filling up ChemReac-Array with forward rate coeff., switching reactants with products and setting new energies
      IF(DSMC%BackwardReacRate) THEN
        IF(ChemReac%QKProcedure(iReac)) THEN
          CALL abort(__STAMP__,&
          'Automatic calculation of backward reaction rate not supported with Q-K reaction:',iReac)
        END IF
        IF(iReac.GT.(ChemReac%NumOfReact/2)) THEN
          IF(TRIM(ChemReac%ReactType(iReac-ChemReac%NumOfReact/2)).EQ.'D') ChemReac%ReactType(iReac) = 'R'
          IF(TRIM(ChemReac%ReactType(iReac-ChemReac%NumOfReact/2)).EQ.'E') ChemReac%ReactType(iReac) = 'E'
          IF(TRIM(ChemReac%ReactType(iReac-ChemReac%NumOfReact/2)).EQ.'i') THEN
            CALL abort(__STAMP__,&
          'Automatic calculation of backward reaction rate not supported with ionization reactions:',iReac)
          END IF
          IF(TRIM(ChemReac%ReactType(iReac-ChemReac%NumOfReact/2)).EQ.'x') THEN
            CALL abort(__STAMP__,&
          'Automatic calculation of backward reaction rate not supported with CEX/MEX reactions:',iReac)
          END IF
          ChemReac%DefinedReact(iReac,1,:)      = ChemReac%DefinedReact(iReac-ChemReac%NumOfReact/2,2,:)
          ChemReac%DefinedReact(iReac,2,:)      = ChemReac%DefinedReact(iReac-ChemReac%NumOfReact/2,1,:)
          ChemReac%Arrhenius_Prefactor(iReac)   = ChemReac%Arrhenius_Prefactor(iReac-ChemReac%NumOfReact/2)
          ChemReac%Arrhenius_Powerfactor(iReac) = ChemReac%Arrhenius_Powerfactor(iReac-ChemReac%NumOfReact/2)
          ChemReac%EActiv(iReac) = ChemReac%EForm(iReac-ChemReac%NumOfReact/2) + ChemReac%EActiv(iReac-ChemReac%NumOfReact/2)
          ChemReac%EForm(iReac) = -ChemReac%EForm(iReac-ChemReac%NumOfReact/2)
        END IF
        ! Calculation of stoichiometric coefficients
        ALLOCATE(ChemReac%ReactInfo(iReac)%StoichCoeff(nSpecies,2))
        ChemReac%ReactInfo(iReac)%StoichCoeff(1:nSpecies,1:2) = 0
        DO iSpec=1, nSpecies
          DO iChemDir=1,2
            IF(ChemReac%DefinedReact(iReac,iChemDir,1).EQ.iSpec) THEN
              ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
            END IF
            IF(ChemReac%DefinedReact(iReac,iChemDir,2).EQ.iSpec) THEN
              ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
            END IF
            IF(ChemReac%DefinedReact(iReac,iChemDir,3).EQ.iSpec) THEN
              ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) = ChemReac%ReactInfo(iReac)%StoichCoeff(iSpec,iChemDir) + 1
            END IF
          END DO
        END DO
      END IF
      ! Proof of reactant definition
      IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)*ChemReac%DefinedReact(iReac,1,3)).EQ.0) THEN
          CALL abort(&
              __STAMP__,&
          'Error in Definition Reactants of Chemical Reaction: ',iReac)
        END IF
      ELSE 
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)).EQ.0) THEN
          CALL abort(&
              __STAMP__,&
          'Error in Definition Reactants of Chemical Reaction: ',iReac)
        END IF
      END IF
      ! Proof of product definition
      IF (TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)*ChemReac%DefinedReact(iReac,2,3)).EQ.0) THEN
          CALL abort(&
              __STAMP__,&
          'Error in Definition Reactants of Chemical Reaction: ',iReac)
        END IF
      ELSE 
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)).EQ.0) THEN
          CALL abort(&
              __STAMP__,&
          'Error in Definition Reactants of Chemical Reaction: ',iReac)
        END IF
      END IF
      MaxSpecies = MAXVAL(ChemReac%DefinedReact(iReac,1:2,1:3))
      IF(MaxSpecies.GT.nSpecies) THEN
        CALL abort(__STAMP__,&
          'Error in Definition Reactants of Chemical Reaction, wrong species index: ',iReac)
      END IF
    END DO
    
    ALLOCATE(PairCombID(nSpecies, nSpecies))
    PairCombID = 0
    CALL DSMC_BuildChem_IDArray(PairCombID)

    ! Case 19: only ion recombination
    DO iReac = 1, ChemReac%NumOfReact
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'rQK') THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        Reactant3 = ChemReac%DefinedReact(iReac,1,3)
        IF(.NOT.YetDefined_Help(iReac)) THEN
          ChemReac%ReactCase(Reactant1, Reactant2) = 19
          ChemReac%ReactCase(Reactant2, Reactant1) = 19
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
        END IF
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
        DummyRecomb(Reactant1,Reactant2) = iReac
      END IF
    END DO  
    ! Case 18: only electron impact ionization
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'iQK').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 18
          ChemReac%ReactCase(Reactant2, Reactant1) = 18
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO  
    ! Case 16: only simple CEX/MEX possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'x').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 16
          ChemReac%ReactCase(Reactant2, Reactant1) = 16
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO  
    ! Case 6/17: associative ionization (N + N -> N2(ion) + e) and recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'i').AND.(.NOT.YetDefined_Help(iReac))) THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        ChemReac%ReactCase(Reactant1, Reactant2) = 6
        ChemReac%ReactCase(Reactant2, Reactant1) = 6
        ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
        ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
        YetDefined_Help(iReac) = .TRUE.
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
              Reactant1 = ChemReac%DefinedReact(iReac,1,1)
              Reactant2 = ChemReac%DefinedReact(iReac,1,2)
              ChemReac%ReactCase(Reactant1, Reactant2) = 17
              ChemReac%ReactCase(Reactant2, Reactant1) = 17
              ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
              ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
              ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
              ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
              YetDefined_Help(iReac2) = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 5/7/8/9/10/11/12: At least two dissociations possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF (iReac.NE.iReac2) THEN
            IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac2))) THEN
              IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                  PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                    ! Case 5: Two dissociations only
                    Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                    Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                    ChemReac%ReactCase(Reactant1, Reactant2) = 5
                    ChemReac%ReactCase(Reactant2, Reactant1) = 5
                    ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                    ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                    ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                    ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                    DO iReac3 = 1, ChemReac%NumOfReact
                      IF ((iReac3.NE.iReac2).AND.(iReac3.NE.iReac)) THEN
                        IF(DSMC%NumPolyatomMolecs.GT.0) THEN ! 3 diss only possible with polyatomic molecules
                          IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                            IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                  ! Case 7: Three dissociations
                                  Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                  Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                  ChemReac%ReactCase(Reactant1, Reactant2) = 7
                                  ChemReac%ReactCase(Reactant2, Reactant1) = 7
                                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                  ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                  ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                  ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                  ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                  YetDefined_Help(iReac3) = .TRUE.
                              DO iReac4 = 1, ChemReac%NumOfReact
                                IF ((iReac4.NE.iReac3).AND.(iReac4.NE.iReac2).AND.(iReac4.NE.iReac)) THEN
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 8: Four dissociations
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 8
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 8
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 9: Three dissociations and one exchange
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 9
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 9
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
  !                                IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac4))) THEN
  !                                  IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
  !                                      PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
  !                                        ! Case 16: Three dissociations and one recombination
  !                                        Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
  !                                        Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
  !                                        ChemReac%ReactCase(Reactant1, Reactant2) = 16
  !                                        ChemReac%ReactCase(Reactant2, Reactant1) = 16
  !                                        ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
  !                                        ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
  !                                        ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
  !                                        ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
  !                                        ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
  !                                        ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
  !                                        YetDefined_Help(iReac4) = .TRUE.
  !                                  END IF
  !                                END IF
                                END IF
                              END DO
                            END IF
                          END IF
                        END IF
                        IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                          IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                              PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                ! Case 10: Two dissociations and one exchange
                                Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                ChemReac%ReactCase(Reactant1, Reactant2) = 10
                                ChemReac%ReactCase(Reactant2, Reactant1) = 10
                                ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                YetDefined_Help(iReac3) = .TRUE.
                                DO iReac4 = 1, ChemReac%NumOfReact
                                  IF ((TRIM(ChemReac%ReactType(iReac4)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac4))) THEN
                                    IF (PairCombID(ChemReac%DefinedReact(iReac4,1,1),ChemReac%DefinedReact(iReac4,1,2)).EQ.&
                                        PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2))) THEN
                                          ! Case 11: Two dissociations, one exchange and recombination
                                          Reactant1 = ChemReac%DefinedReact(iReac4,1,1)
                                          Reactant2 = ChemReac%DefinedReact(iReac4,1,2)
                                          ChemReac%ReactCase(Reactant1, Reactant2) = 11
                                          ChemReac%ReactCase(Reactant2, Reactant1) = 11
                                          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                          ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                          ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                          ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                          ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                          ChemReac%ReactNum(Reactant1, Reactant2, 4) = iReac4
                                          ChemReac%ReactNum(Reactant2, Reactant1, 4) = iReac4
                                          YetDefined_Help(iReac4) = .TRUE.
                                    END IF
                                  END IF
                                END DO
                          END IF
                        END IF
                        IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                          IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                              PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                                ! Case 12: Two dissociations and one recombination
                                Reactant1 = ChemReac%DefinedReact(iReac3,1,1)
                                Reactant2 = ChemReac%DefinedReact(iReac3,1,2)
                                ChemReac%ReactCase(Reactant1, Reactant2) = 12
                                ChemReac%ReactCase(Reactant2, Reactant1) = 12
                                ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                                ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                                ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                                ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                                ChemReac%ReactNum(Reactant1, Reactant2, 3) = iReac3
                                ChemReac%ReactNum(Reactant2, Reactant1, 3) = iReac3
                                YetDefined_Help(iReac3) = .TRUE.
                          END IF
                        END IF
                      END IF
                    END DO
                    YetDefined_Help(iReac) = .TRUE.
                    YetDefined_Help(iReac2) = .TRUE.
                    ChemReac%MeanEVib_Necc = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Cases 4/13/14: Dissociation, recombination and exchange reactions possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  ! Case 4: One dissociation and one exchange
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 4
                  ChemReac%ReactCase(Reactant2, Reactant1) = 4
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                  ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
                  ChemReac%MeanEVib_Necc = .TRUE.
                  DO iReac3 = 1, ChemReac%NumOfReact
                    IF ((TRIM(ChemReac%ReactType(iReac3)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac3))) THEN
                      IF (PairCombID(ChemReac%DefinedReact(iReac3,1,1),ChemReac%DefinedReact(iReac3,1,2)).EQ.&
                          PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                            ! Case 13: One dissociation, one recombination and one exchange
                            Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                            Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                            ChemReac%ReactCase(Reactant1, Reactant2) = 13
                            ChemReac%ReactCase(Reactant2, Reactant1) = 13
                            ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                            ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                            ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                            ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                            YetDefined_Help(iReac3) = .TRUE.
                      END IF
                    END IF
                  END DO
            END IF
          END IF
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  ! Case 14: One dissociation and one recombination
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 14
                  ChemReac%ReactCase(Reactant2, Reactant1) = 14
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
                  ChemReac%MeanEVib_Necc = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 15: One exchange and one recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((TRIM(ChemReac%ReactType(iReac2)).EQ.'R').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                  Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                  Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                  ChemReac%ReactCase(Reactant1, Reactant2) = 15
                  ChemReac%ReactCase(Reactant2, Reactant1) = 15
                  ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                  ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                  YetDefined_Help(iReac) = .TRUE.
                  YetDefined_Help(iReac2) = .TRUE.
                  ChemReac%MeanEVib_Necc = .TRUE.
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 3: only exchange reaction possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'E').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 3
          ChemReac%ReactCase(Reactant2, Reactant1) = 3
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
          ChemReac%MeanEVib_Necc = .TRUE.
      END IF
    END DO
    ! Case 2: only dissociation possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((TRIM(ChemReac%ReactType(iReac)).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 2
          ChemReac%ReactCase(Reactant2, Reactant1) = 2
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
          ChemReac%MeanEVib_Necc = .TRUE.
      END IF
    END DO
    ! Case 1: only molecular recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF(TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        Reactant3 = ChemReac%DefinedReact(iReac,1,3)
        IF(.NOT.YetDefined_Help(iReac)) THEN
          ChemReac%ReactCase(Reactant1, Reactant2) = 1
          ChemReac%ReactCase(Reactant2, Reactant1) = 1
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
          ChemReac%MeanEVib_Necc = .TRUE.
        END IF
        ChemReac%ReactNumRecomb(Reactant1, Reactant2, Reactant3) = iReac
        ChemReac%ReactNumRecomb(Reactant2, Reactant1, Reactant3) = iReac
        DummyRecomb(Reactant1,Reactant2) = iReac
      END IF
    END DO
    ! Filling empty values of ReactNumRecomb with the last recombination reaction for that collision pair
    DO iReac = 1, ChemReac%NumOfReact
      IF((TRIM(ChemReac%ReactType(iReac)).EQ.'R').OR.(TRIM(ChemReac%ReactType(iReac)).EQ.'rQK')) THEN
        Reactant1 = ChemReac%DefinedReact(iReac,1,1)
        Reactant2 = ChemReac%DefinedReact(iReac,1,2)
        DO iSpec = 1, nSpecies
          IF(ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec).EQ.0) THEN
            ChemReac%ReactNumRecomb(Reactant1, Reactant2, iSpec) = DummyRecomb(Reactant1,Reactant2)
            ChemReac%ReactNumRecomb(Reactant2, Reactant1, iSpec) = DummyRecomb(Reactant1,Reactant2)
          END IF
        END DO
      END IF
    END DO

    DEALLOCATE(PairCombID)
    DEALLOCATE(DummyRecomb)
    CALL Calc_Arrhenius_Factors()
  ELSE
    SWRITE(*,'(A)') 'NO REACTIONS DEFINED!'
  END IF

END SUBROUTINE DSMC_chemical_init


SUBROUTINE DSMC_BuildChem_IDArray(PairCombID)
!===================================================================================================================================
! Builds array for the identification of different species combinations
!===================================================================================================================================
! MODULES
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(INOUT)  :: PairCombID(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: iSpec, jSpec, iComb
!===================================================================================================================================
  iComb = 1
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec,  nSpecies
        PairCombID(iSpec, jSpec) = iComb
        PairCombID(jSpec, iSpec) = iComb
        iComb = iComb + 1
    END DO
  END DO 

END SUBROUTINE DSMC_BuildChem_IDArray


SUBROUTINE Calc_Arrhenius_Factors()
!===================================================================================================================================
! Calculates variables for Arrhenius calculation (H_ab, Beta) for diatomic chemical reactions 
!===================================================================================================================================
! MODULES
  USE MOD_Globals,            ONLY : Abort
  USE MOD_DSMC_Vars,          ONLY : ChemReac, DSMC, SpecDSMC, CollInf
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies, BoltzmannConst
  USE MOD_Globals_Vars,       ONLY : Pi
  USE MOD_DSMC_Analyze,       ONLY : CalcTVib
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: iQuaMax1, iQuaMax2, iQua1, iQua2, iReac, iQuaMax1_temp, iQuaMax2_temp
  INTEGER                 :: iQuaMax3, iQua3, iSpec
  REAL                    :: Xi_vib1, Xi_vib2, Xi_vib3
!  REAL                    :: EVib1, TVib1, EVib2, TVib2              ! needed for TSHO
!===================================================================================================================================

  DO iReac = 1, ChemReac%NumOfReact
    ! Calculate the Arrhenius arrays only if the reaction is not a QK reaction
    IF (.NOT.ChemReac%QKProcedure(iReac)) THEN
      ! Compute VHS Factor H_ab necessary for reaction probs, only defined for one omega for all species (laux diss page 24)
      ChemReac%Hab(iReac) = GAMMA(2.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) * 2.0 &
              * CollInf%Cab(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2))) &
              / ((1 + CollInf%KronDelta(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))) &
              * SQRT(Pi)) * (2.0 * BoltzmannConst &
              / CollInf%MassRed(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1), ChemReac%DefinedReact(iReac,1,2)))) &
              ** (0.5 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
      IF((.NOT.SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%PolyatomicMol) &
          .AND.(.NOT.SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%PolyatomicMol)) THEN
        !Recombination
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'R') THEN
          ChemReac%MeanEVib_Necc = .TRUE.
          IF(SpecDSMC(ChemReac%DefinedReact(iReac,1,3))%PolyatomicMol &
            .OR.SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%PolyatomicMol) CYCLE
          iQuaMax3 = MAXVAL(SpecDSMC(:)%MaxVibQuant)
          ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(nSpecies,0:iQuaMax3-1))
          ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(:,:) = 0.0
          DO iSpec=1,nSpecies
            IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
              DO iQua3 = 0 , iQuaMax3-1
                IF (iQua3.NE.0) THEN     
                  Xi_vib3 = 2 * iQua3 * LOG(1.0/ iQua3 + 1.0 )
                ELSE
                  Xi_vib3 = 0
                END IF
                ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iSpec,iQua3) = ChemReac%Arrhenius_Prefactor(iReac) &
                                         * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                         * GAMMA(( 3.0 + Xi_vib3 + SpecDSMC(iSpec)%Xi_Rot)/2 &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 ) &
                                         / (ChemReac%Hab(iReac) * GAMMA(( 3.0 + Xi_vib3 &
                                         + SpecDSMC(iSpec)%Xi_Rot)/2 &
                                         + 1.5 + ChemReac%Arrhenius_Powerfactor(iReac) ))
              END DO
            ELSE
              ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iSpec,0) = ChemReac%Arrhenius_Prefactor(iReac) &
                                         * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                         * GAMMA(1.5 &
                                         - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 ) &
                                         / ( ChemReac%Hab(iReac) * GAMMA(3.0 + ChemReac%Arrhenius_Powerfactor(iReac)))
            END IF
          END DO
        END IF
        !Dissociation
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'D') THEN
          ChemReac%MeanEVib_Necc = .TRUE.
          iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
          iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
          iQuaMax1_temp = 100
          iQuaMax2_temp = 100
          !!!!!!would be needed in TSHO Case!!!!!
          ! iQuaMax1_temp = iQuaMax1 - 1 
          ! iQuaMax2_temp = iQuaMax2 - 1 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(iQuaMax2.EQ.0) THEN
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(0:iQuaMax1-1,0:1))
          ELSE
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:iQuaMax2-1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
          END IF
          DO iQua1 = 0 , iQuaMax1_temp      
            !Calculation of the vibrational DOF in TSHO Model
            IF (iQua1.NE.0) THEN     
              Xi_vib1 = 2 * iQua1 * LOG(1.0/ iQua1 + 1.0 )
              !!!!!!would be needed in TSHO Case!!!!!
              !  EVib1 = (iQua1 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib
              !  TVib1 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib & 
              !          , EVib1, iQuaMax1) 
              !  Xi_vib1 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib / TVib1) - 1)) &
              !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1 / TVib1) - 1))
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
            ELSE
              Xi_vib1 = 0
            END IF
            
            IF(iQuaMax2.EQ.0) THEN
              Xi_vib2 = 0
              ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
              IF (Xi_vib1.GT.0) THEN
                IF (SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor.EQ.0) THEN
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2))
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(Xi_vib1/2) * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                      + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                      + Xi_vib1/2))
                END IF
              ELSE
                ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = 0
              END IF
            ELSE
              DO iQua2 = 0, iQuaMax2_temp
                !Calculation of the vibrational DOF in TSHO Model
                IF (iQua2.NE.0) THEN
                  Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                  !!!!!!would be needed in TSHO Case!!!!!
                  !  EVib2 = (iQua2 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib 
                  ! TVib2 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib & 
                  !          , EVib2, iQuaMax2) 
                  !  Xi_vib2 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib / TVib2) - 1)) &
                  !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2 / TVib2) - 1))
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                ELSE
                  Xi_vib2 = 0
                END IF
                !also only defined for one omega VHS for each species
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2) = Xi_vib1 + Xi_vib2  &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) 
                IF (Xi_vib1.GT.0) THEN
                  IF (SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor.EQ.0) THEN
                    ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2))
                  ELSE
                    ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(Xi_vib1/2) * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                        + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                        + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        * GAMMA(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                        + Xi_vib1/2))
                  END IF
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = 0
                END IF
              END DO
            END IF
          END DO
          !! would not be needed in TSHO Case
          IF (iQuaMax2.EQ.0) THEN
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
          ELSE
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(100,100)
          END IF
        END IF
      
       ! Exchange Reaction
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'E') THEN
          ChemReac%MeanEVib_Necc = .TRUE.
          iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
          iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
          IF(iQuaMax1.EQ.0) THEN
            CALL Abort(__STAMP__,&
              'ERROR in DSMCSpecies.ini: The first defined particle in exchange reaction has to be a molecule', iReac)
          END IF
          iQuaMax1_temp = 100
          iQuaMax2_temp = 100
          !!!!!!would be needed in TSHO Case!!!!!
          ! iQuaMax1_temp = iQuaMax1 - 1 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(iQuaMax2.EQ.0) THEN
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:iQuaMax1-1,0:1))
          ELSE
            ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:iQuaMax2-1))
            ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
          END IF

          DO iQua1 = 0 , iQuaMax1_temp      
            !Calculation of the vibrational DOF in TSHO Model
            IF (iQua1.NE.0) THEN     
              Xi_vib1 = 2 * iQua1 * LOG(1.0/ iQua1 + 1.0 )
              !!!!!!would be needed in TSHO Case!!!!!
              !  EVib1 = (iQua1 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib
              !  TVib1 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib & 
              !          , EVib1, iQuaMax1) 
              !  Xi_vib1 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib / TVib1) - 1)) &
              !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1  &
              !          / (TVib1 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%CharaTVib * iQuaMax1 / TVib1) - 1))
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
            ELSE
              Xi_vib1 = 0
            END IF

            IF(iQuaMax2.EQ.0) THEN
              Xi_vib2 = 0
              ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
              IF (Xi_vib1.GT.0) THEN
                ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                    / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2))
              ELSE
                ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,0) = 0
              END IF
            ELSE
              DO iQua2 = 0, iQuaMax2_temp
                !Calculation of the vibrational DOF in TSHO Model
                IF (iQua2.NE.0) THEN
                  Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                  !!!!!!would be needed in TSHO Case!!!!!
                  !  EVib2 = (iQua2 + DSMC%GammaQuant) * BoltzmannConst * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib 
                  ! TVib2 = CalcTVib(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib & 
                  !          , EVib2, iQuaMax2) 
                  !  Xi_vib2 = 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib / TVib2) - 1)) &
                  !          - 2 * SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2  &
                  !          / (TVib2 * (EXP(SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%CharaTVib * iQuaMax2 / TVib2) - 1))
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
                ELSE
                  Xi_vib2 = 0
                END IF
                !also only defined for one omega VHS for each species
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2) = Xi_vib1 + Xi_vib2  &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) 
                IF (Xi_vib1.GT.0) THEN
                  ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2))
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1,iQua2) = 0
                END IF
              END DO
            END IF
          END DO
          !! would not be needed in TSHO Case
          IF (iQuaMax2.EQ.0) THEN
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
          ELSE
            ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(100,100)
          END IF
        END IF

        ! Ionization
        IF (TRIM(ChemReac%ReactType(iReac)).EQ.'i') THEN
          ChemReac%MeanEVib_Necc = .TRUE.
          iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
          iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
          iQuaMax1_temp = 100
          iQuaMax2_temp = 100
          IF(iQuaMax1.EQ.0) THEN
            IF(iQuaMax2.EQ.0) THEN
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:1, 0:1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0:1,0:1))
            ELSE
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:1, 0:iQuaMax2-1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0:1,0:iQuaMax2-1))
            END IF
          ELSE
            IF(iQuaMax2.EQ.0) THEN
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0:iQuaMax1-1,0:1))
            ELSE
              ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:iQuaMax2-1))
              ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0:iQuaMax1-1,0:iQuaMax2-1))
            END IF
          END IF
          IF(iQuaMax1.EQ.0) THEN
            IF(iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(0,0) = 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
              ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(0,0)/2) &
                    / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(0,0)/2))
            ELSE
              DO iQua2 = 0, iQuaMax2_temp
                ! Calculation of the vibrational DOF in TSHO Model
                IF (iQua2.NE.0) THEN
                  Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                ELSE
                  Xi_vib2 = 0
                END IF
                ! also only defined for one omega VHS for each species
                ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2) = Xi_vib2  &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                        + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) 
                ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(0,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2)/2) &
                    / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(0,iQua2)/2))
              END DO
            END IF
          ELSE
            DO iQua1 = 0 , iQuaMax1_temp      
              ! Calculation of the vibrational DOF in TSHO Model
              IF (iQua1.NE.0) THEN     
                Xi_vib1 = 2 * iQua1 * LOG(1.0/ iQua1 + 1.0 )
              ELSE
                Xi_vib1 = 0
              END IF
              
              IF(iQuaMax2.EQ.0) THEN
                Xi_vib2 = 0
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                          + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
                IF (Xi_vib1.GT.0) THEN
                  ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                      *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                      - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                      * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2) &
                      / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                      + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2))
                ELSE
                  ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(iQua1,0) = 0
                END IF
              ELSE
                DO iQua2 = 0, iQuaMax2_temp
                  ! Calculation of the vibrational DOF in TSHO Model
                  IF (iQua2.NE.0) THEN
                    Xi_vib2 = 2 * iQua2 * LOG(1.0/ iQua2 + 1.0 )
                  ELSE
                    Xi_vib2 = 0
                  END IF
                  ! also only defined for one omega VHS for each species
                  ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2) = Xi_vib1 + Xi_vib2  &
                          + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                          + SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%Xi_Rot &
                          + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) 
                  IF (Xi_vib1.GT.0) THEN
                    ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                        *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                        - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                        * GAMMA(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2) &
                        / (ChemReac%Hab(iReac) * GAMMA(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                        + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2))
                  ELSE
                    ChemReac%ReactInfo(iReac)%Beta_Ion_Arrhenius(iQua1,iQua2) = 0
                  END IF
                END DO
              END IF
            END DO
          END IF
    !      !! would not be needed in TSHO Case
          IF (iQuaMax1.EQ.0) THEN
            IF (iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(0,0)
            ELSE
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(0,100)
            END IF
          ELSE
            IF (iQuaMax2.EQ.0) THEN
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
            ELSE
              ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,101:iQuaMax2-1) = ChemReac%ReactInfo(iReac)%Xi_Total(100,100)
            END IF
          END IF
        END IF
      END IF ! First and second educt not polyatomic
    END IF ! Not Q-K
  END DO

END SUBROUTINE Calc_Arrhenius_Factors


SUBROUTINE InitPartitionFunction(iSpec, PartitionArraySize)
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,       ONLY: Pi, PlanckConst
USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY: BoltzmannConst, Species
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iSpec, PartitionArraySize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: iInter, iPolyatMole, iDOF
  REAL                            :: Qtra, Qrot, Qvib, Qelec, Temp
!===================================================================================================================================

  DO iInter = 1, PartitionArraySize
    Temp = iInter * DSMC%PartitionInterval
    Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
    IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
          Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
        ELSE
          Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
                                                                          * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
                                                                          * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
        END IF
        Qvib = 1.
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                  / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
        END DO
      ELSE
        Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
        Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
      END IF
    ELSE
      Qrot = 1.
      Qvib = 1.
    END IF
    Qelec = 0.
    DO iDOF=1, SpecDSMC(iSpec)%NumElecLevels
      Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
    END DO
    SpecDSMC(iSpec)%PartitionFunction(iInter) = Qtra * Qrot * Qvib * Qelec
  END DO


END SUBROUTINE InitPartitionFunction


END MODULE MOD_DSMC_ChemInit
