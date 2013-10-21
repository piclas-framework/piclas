MODULE MOD_DSMC_ChemInit
!===================================================================================================================================
! Add comments please!
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

PUBLIC :: DSMC_chemical_init
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_chemical_init()

  USE MOD_DSMC_Vars,          ONLY : ChemReac
  USE MOD_ReadInTools
  USE MOD_Globals,            ONLY : UNIT_StdOut
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies, BoltzmannConst
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf 
  INTEGER               :: iReac, iReac2
  INTEGER, ALLOCATABLE  :: PairCombID(:,:)
  LOGICAL, ALLOCATABLE  :: YetDefined_Help(:)
  INTEGER               :: Reactant1, Reactant2, Reactant3
!#ifdef MPI
!#endif
!===================================================================================================================================
  
! reading reaction values   
  ChemReac%NumOfReact = GETINT('DSMC-NumOfReactions','0')
  IF (ChemReac%NumOfReact.GT.0) THEN    
    ALLOCATE(ChemReac%NumReac(ChemReac%NumOfReact))
    ChemReac%NumReac = 0
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
    ALLOCATE(ChemReac%Arrhenius_Prefactor(ChemReac%NumOfReact),&
             ChemReac%Arrhenius_Powerfactor(ChemReac%NumOfReact),&
             ChemReac%EActiv(ChemReac%NumOfReact),&
             ChemReac%EForm(ChemReac%NumOfReact))
    ALLOCATE(ChemReac%MeanEVibQua_PerIter(nSpecies))

    DO iReac = 1, ChemReac%NumOfReact
      WRITE(UNIT=hilf,FMT='(I2)') iReac
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
      !proof of reactant definition
      IF (ChemReac%ReactType(iReac).EQ.'R') THEN
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)*ChemReac%DefinedReact(iReac,1,3)).EQ.0) THEN
          WRITE(*,*) 'Error in Definition Reactants of Chemical Reaction: ', iReac
          STOP
        END IF
      ELSE 
        IF ((ChemReac%DefinedReact(iReac,1,1)*ChemReac%DefinedReact(iReac,1,2)).EQ.0) THEN
          WRITE(*,*) 'Error in Definition Reactants of Chemical Reaction: ', iReac
          STOP
        END IF
      END IF
      !proof of product definition
      IF ((ChemReac%ReactType(iReac).EQ.'i').OR.(ChemReac%ReactType(iReac).EQ.'D')) THEN
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)*ChemReac%DefinedReact(iReac,2,3)).EQ.0) THEN
          WRITE(*,*) 'Error in Definition Products of Chemical Reaction: ', iReac
          STOP
        END IF
      ELSE 
        IF ((ChemReac%DefinedReact(iReac,2,1)*ChemReac%DefinedReact(iReac,2,2)).EQ.0) THEN
          WRITE(*,*) 'Error in Definition Products of Chemical Reaction: ', iReac
          STOP
        END IF
      END IF
    END DO
    
    ALLOCATE(PairCombID(nSpecies, nSpecies))
    PairCombID = 0
    CALL DSMC_BuildChem_IDArray(PairCombID)
  
    ! Case 6: only ionization possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'I').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 6
          ChemReac%ReactCase(Reactant2, Reactant1) = 6
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 4: Dissociation and exchange reaction possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF ((ChemReac%ReactType(iReac2).EQ.'E').AND.(.NOT.YetDefined_Help(iReac2))) THEN
            IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
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
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 5: Two Dissociations possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
        DO iReac2 = 1, ChemReac%NumOfReact
          IF (iReac.NE.iReac2) THEN
            IF ((ChemReac%ReactType(iReac2).EQ.'D').AND.(.NOT.YetDefined_Help(iReac2))) THEN
              IF (PairCombID(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)).EQ.&
                  PairCombID(ChemReac%DefinedReact(iReac2,1,1),ChemReac%DefinedReact(iReac2,1,2))) THEN
                    Reactant1 = ChemReac%DefinedReact(iReac,1,1)
                    Reactant2 = ChemReac%DefinedReact(iReac,1,2)
                    ChemReac%ReactCase(Reactant1, Reactant2) = 5
                    ChemReac%ReactCase(Reactant2, Reactant1) = 5
                    ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
                    ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
                    ChemReac%ReactNum(Reactant1, Reactant2, 2) = iReac2
                    ChemReac%ReactNum(Reactant2, Reactant1, 2) = iReac2
                    YetDefined_Help(iReac) = .TRUE.
                    YetDefined_Help(iReac2) = .TRUE.
              END IF
            END IF
          END IF
        END DO
      END IF
    END DO
    ! Case 3: only exchange reaction possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'E').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 3
          ChemReac%ReactCase(Reactant2, Reactant1) = 3
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 2: only dissociation possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'D').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          ChemReac%ReactCase(Reactant1, Reactant2) = 2
          ChemReac%ReactCase(Reactant2, Reactant1) = 2
          ChemReac%ReactNum(Reactant1, Reactant2, 1) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, 1) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO
    ! Case 1: only molecular recombination possible
    DO iReac = 1, ChemReac%NumOfReact
      IF ((ChemReac%ReactType(iReac).EQ.'R').AND.(.NOT.YetDefined_Help(iReac))) THEN
          Reactant1 = ChemReac%DefinedReact(iReac,1,1)
          Reactant2 = ChemReac%DefinedReact(iReac,1,2)
          Reactant3 = ChemReac%DefinedReact(iReac,1,3)
          ChemReac%ReactCase(Reactant1, Reactant2) = 1
          ChemReac%ReactCase(Reactant2, Reactant1) = 1
          ChemReac%ReactNum(Reactant1, Reactant2, Reactant3) = iReac
          ChemReac%ReactNum(Reactant2, Reactant1, Reactant3) = iReac
          YetDefined_Help(iReac) = .TRUE.
      END IF
    END DO

    DEALLOCATE(PairCombID)
    CALL Calc_Arrhenius_Factors()
  ELSE
    WRITE(*,'(A)') 'NO REACTIONS DEFINED!'
  END IF

END SUBROUTINE DSMC_chemical_init

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_BuildChem_IDArray(PairCombID)

  USE MOD_DSMC_Vars,          ONLY : ChemReac
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER                 :: iSpec, jSpec, iComb
  INTEGER, INTENT(INOUT)  :: PairCombID(:,:)
!#ifdef MPI
!#endif
!===============================================================================================
  iComb = 1
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec,  nSpecies
        PairCombID(iSpec, jSpec) = iComb
        PairCombID(jSpec, iSpec) = iComb
        iComb = iComb + 1
    END DO
  END DO 

END SUBROUTINE DSMC_BuildChem_IDArray

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE Calc_Arrhenius_Factors()

  USE MOD_DSMC_Vars,          ONLY : ChemReac, DSMC, SpecDSMC, CollInf
  USE MOD_PARTICLE_Vars,      ONLY : nSpecies, BoltzmannConst
  USE MOD_Equation_Vars,      ONLY : Pi
  USE MOD_DSMC_Analyze,       ONLY : CalcTVib
  USE nr,                     ONLY : gammln
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER                 :: iQuaMax1, iQuaMax2, iQua1, iQua2, iReac, iQuaMax1_temp, iQuaMax2_temp
  INTEGER                 :: iQuaMax3, iQuaMax3_temp, iQua3
  REAL                    :: EVib1, EVib2, TVib1, TVib2, Xi_vib1, Xi_vib2, H_ab, Xi_vib3
!#ifdef MPI
!#endif
!===============================================================================================

ALLOCATE(ChemReac%ReactInfo(ChemReac%NumOfReact))

DO iReac = 1, ChemReac%NumOfReact
  ! calculate the Arrhenius Arrays only if the reaction is not a QK reaction
  IF ( ChemReac%QKProcedure(iReac) .EQV. .false. ) THEN

    ! Compute VHS Factor H_ab necessary for reaction probs (laux diss page 24)
    ! only defined for one omega for al species
    H_ab = exp(gammln(2-SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) * 2 &
          * CollInf%Cab(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2))) &
          / ((1 + CollInf%KronDelta(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1),ChemReac%DefinedReact(iReac,1,2)))) &
          * SQRT(Pi)) &
          * (2 * BoltzmannConst &
          / CollInf%MassRed(CollInf%Coll_Case(ChemReac%DefinedReact(iReac,1,1), ChemReac%DefinedReact(iReac,1,2)))) &
          ** (0.5 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
    !Recombination
    IF (ChemReac%ReactType(iReac).EQ.'R') THEN
      ChemReac%MeanEVib_Necc = .TRUE.
      iQuaMax3 = SpecDSMC(ChemReac%DefinedReact(iReac,1,3))%MaxVibQuant
      iQuaMax3_temp = 100
      !!!!!!would be needed in TSHO Case!!!!!
      ! iQuaMax3_temp = iQuaMax3 - 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(iQuaMax3.EQ.0) THEN
        ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(0:1))
        ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(0) = ChemReac%Arrhenius_Prefactor(iReac) &
                                   * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                   - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                   * exp(gammln(1.5 &
                                   - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 )) &
                                   / ( H_ab * exp(gammln(3.0 + ChemReac%Arrhenius_Powerfactor(iReac))))
      ELSE
        ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(0:iQuaMax3-1))
        DO iQua3 = 0 , iQuaMax3_temp
          IF (iQua3.NE.0) THEN     
            Xi_vib3 = 2 * iQua3 * LOG(1.0/ iQua3 + 1.0 )
            !!!!!!would be needed in TSHO Case!!!!!
            !  ...
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
          ELSE
            Xi_vib3 = 0
          END IF 
          ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(iQua3) = ChemReac%Arrhenius_Prefactor(iReac) &
                                   * BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                                   - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS) &
                                   * exp(gammln(( 3.0 + Xi_vib3 + SpecDSMC(ChemReac%DefinedReact(iReac,1,3))%Xi_Rot)/2 &
                                   - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + 2 )) &
                                   / ( H_ab * exp(gammln(( 3.0 + Xi_vib3 + SpecDSMC(ChemReac%DefinedReact(iReac,1,3))%Xi_Rot)/2 &
                                   + 1.5 + ChemReac%Arrhenius_Powerfactor(iReac) )))
        END DO
      END IF 
    END IF
    !Dissociation
    IF (ChemReac%ReactType(iReac).EQ.'D') THEN
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
                  * exp(gammln(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
                  / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                  + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)))
            ELSE
              ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,0) = ChemReac%Arrhenius_Prefactor(iReac) & 
                  *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                  - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                  * exp(gammln(Xi_vib1/2)) * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                  + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
                  / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                  + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
                  * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                  + Xi_vib1/2)))
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
                    * exp(gammln(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2)) &
                    / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2)))
              ELSE
                ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(iQua1,iQua2) = ChemReac%Arrhenius_Prefactor(iReac) & 
                    *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                    - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                    * exp(gammln(Xi_vib1/2)) * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                    + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2)) &
                    / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                    + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                    + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,iQua2)/2)) &
                    * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                    + Xi_vib1/2)))
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
    IF (ChemReac%ReactType(iReac).EQ.'E') THEN
      ChemReac%MeanEVib_Necc = .TRUE.
      iQuaMax1 = SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%MaxVibQuant
      iQuaMax2 = SpecDSMC(ChemReac%DefinedReact(iReac,1,2))%MaxVibQuant
      IF (iQuaMax2.NE.0) THEN
        WRITE(*,*) 'ERROR: The second Partner of an exchange reactions must be an atom! ReactNum:', iReac
        STOP
      END IF
      iQuaMax1_temp = 100
      !!!!!!would be needed in TSHO Case!!!!!
      ! iQuaMax1_temp = iQuaMax1 - 1 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(ChemReac%ReactInfo(iReac)%Xi_Total(0:iQuaMax1-1, 0:1))
      ALLOCATE(ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(0:iQuaMax1-1))
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
        Xi_vib2 = 0
        ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0) = Xi_vib1+ SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%Xi_Rot &
                  + 2*(2 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)
        IF (Xi_vib1.GT.0) THEN
  ! Das folgende IF ist, ist nur eine Übergangslösung, damit das PHI3 in austauschreaktionen immer 0 ist. Sollte
  ! sich dies mal ändern, muß man hier noch Phi3 pro reaktion einlesen.
  !        IF (SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor.EQ.0) THEN
            ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1) = ChemReac%Arrhenius_Prefactor(iReac) & 
                *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
                - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
                * exp(gammln(ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
                / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
                + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + &
                ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)))
  !        ELSE
  !          ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1) = ChemReac%Arrhenius_Prefactor(iReac) & 
  !              *(BoltzmannConst**(0.5 - ChemReac%Arrhenius_Powerfactor(iReac) &
  !              - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS)) &
  !              * exp(gammln(Xi_vib1/2)) * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
  !              + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
  !              / (H_ab * exp(gammln(ChemReac%Arrhenius_Powerfactor(iReac) - 0.5 &
  !              + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS + ChemReac%ReactInfo(iReac)%Xi_Total(iQua1,0)/2)) &
  !              * exp(gammln(SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
  !              + Xi_vib1/2)))
  !        END IF
        ELSE
          ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(iQua1) = 0
        END IF
      END DO
      !! would not be needed in TSHO Case
      ChemReac%ReactInfo(iReac)%Xi_Total(101:iQuaMax1-1,0) = ChemReac%ReactInfo(iReac)%Xi_Total(100,0)
    END IF
  END IF
END DO

END SUBROUTINE Calc_Arrhenius_Factors

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_DSMC_ChemInit
