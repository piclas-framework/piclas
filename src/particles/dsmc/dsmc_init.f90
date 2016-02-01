#include "boltzplatz.h"

MODULE MOD_DSMC_Init
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitDSMC
  MODULE PROCEDURE InitDSMC
END INTERFACE

INTERFACE FinalizeDSMC
  MODULE PROCEDURE FinalizeDSMC
END INTERFACE

INTERFACE DSMC_SetInternalEnr_LauxVFD
  MODULE PROCEDURE DSMC_SetInternalEnr_LauxVFD
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitDSMC, DSMC_SetInternalEnr_LauxVFD,FinalizeDSMC
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDSMC()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,             ONLY : nElems
USE MOD_Globals_Vars,          ONLY:PI
USE MOD_ReadInTools
USE MOD_DSMC_ElectronicModel,  ONLY: ReadSpeciesLevel
USE MOD_DSMC_Vars
USE MOD_PARTICLE_Vars,         ONLY: nSpecies, BoltzmannConst, Species, PDM, PartSpecies, useVTKFileBGG
USE MOD_DSMC_Analyze,          ONLY: InitHODSMC
USE MOD_TimeDisc_Vars,         ONLY: TEnd
USE MOD_DSMC_ChemInit,         ONLY: DSMC_chemical_init
USE MOD_DSMC_PolyAtomicModel,  ONLY: InitPolyAtomicMolecs, DSMC_SetInternalEnr_Poly, DSMC_SetInternalEnr_PolyFast 
USE MOD_DSMC_PolyAtomicModel,  ONLY: DSMC_SetInternalEnr_PolyFastPart2
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf , hilf2
  INTEGER               :: iCase, iSpec, jSpec, nCase, iPart , iINit
  REAL                  :: A1, A2     ! species constant for cross section (p. 24 Laux)
  REAL                  :: JToEv
  INTEGER,ALLOCATABLE   :: CalcSurfCollis_SpeciesRead(:) !help array for reading surface stuff
  REAL                  :: BGGasEVib
#if ( PP_TimeDiscMethod ==42 )
  CHARACTER(LEN=64)     :: DebugElectronicStateFilename
  INTEGER               :: ii
#endif
!===================================================================================================================================
  JToEv = 1.602176565E-19  
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' DSMC INIT ...'
  
! reading/writing OutputMesh stuff
  DSMC%OutputMeshInit = GETLOGICAL('Particles-DSMC-OutputMeshInit','.FALSE.')
  DSMC%OutputMeshSamp = GETLOGICAL('Particles-DSMC-OutputMeshSamp','.FALSE.')
! reading and reset general DSMC values
  CollisMode = GETINT('Particles-DSMC-CollisMode','1') !0: no collis, 1:elastic col, 2:elast+rela, 3:chem
  DSMC%GammaQuant   = GETREAL('Particles-DSMC-GammaQuant', '0.5')
  DSMC%TimeFracSamp = GETREAL('Part-TimeFracForSampling','0')
  DSMC%NumOutput = GETINT('Particles-NumberForDSMCOutputs','1')
  DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim','.FALSE.')
  DSMC%ReservoirSimuRate = GETLOGICAL('Particles-DSMCReservoirSimRate','.FALSE.')
  DSMC%ReservoirRateStatistic = GETLOGICAL('Particles-DSMCReservoirStatistic','.FALSE.')
  DSMC%VibEnergyModel = GETINT('Particles-ModelForVibrationEnergy','0')
  DSMC%ElectronicStateDatabase = GETSTR('Particles-DSMCElectronicDatabase','none')
  LD_MultiTemperaturMod=GETINT('LD-ModelForMultiTemp','0')
  DSMC%ElectronicState = .FALSE.
  IF ( (DSMC%ElectronicStateDatabase .ne. 'none') .AND. (CollisMode .GT. 1)  ) THEN 
    DSMC%ElectronicState = .TRUE.
  END IF
  IF ((DSMC%VibEnergyModel.EQ.1).AND.(CollisMode.EQ.3)) THEN
    CALL Abort(&
     __STAMP__,&
    'TSHO Model is not working with Chemical Reactions !!!')
  ELSE IF ((DSMC%VibEnergyModel.GT.1).OR.(DSMC%VibEnergyModel.LT.0)) THEN
    CALL Abort(&
     __STAMP__,&
    'ERROR in ModelForVibrationEnergy Flag!')
  END IF
  IF (DSMC%NumOutput.NE.0) THEN
    DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
  END IF
  DSMC%NumPolyatomMolecs = 0
  ALLOCATE(DSMC%CollProbSamp(nElems))
  DSMC%CollProbSamp(1:nElems) = 0.0
  ALLOCATE(DSMC%CollProbOut(nElems,2))
  DSMC%CollProbOut(1:nElems,1:2) = 0.0

! definition of DSMC particle values
  ALLOCATE(DSMC_RHS(PDM%maxParticleNumber,3))
  DSMC_RHS = 0

  IF (nSpecies.LE.0) THEN 
    CALL Abort(&
     __STAMP__,&
     "ERROR: nSpecies .LE. 0:", nSpecies)
  END IF

  IF (CollisMode.EQ.0) THEN
#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001) || (PP_TimeDiscMethod==42)
    CALL Abort(&
     __STAMP__,&
     "Free Molecular Flow (CollisMode=0) is not supported for LD or DEBUG!")
#endif
  ELSE !CollisMode.GT.0

! reading species data of ini_2
  ALLOCATE(SpecDSMC(nSpecies))
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    SpecDSMC(iSpec)%Name    = GETSTR('Part-Species'//TRIM(hilf)//'-SpeciesName','none')
    SpecDSMC(iSpec)%InterID = GETINT('Part-Species'//TRIM(hilf)//'-InteractionID','0')
    SpecDSMC(iSpec)%TrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceTemp','0')
    SpecDSMC(iSpec)%DrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceDiam','0')
    IF((SpecDSMC(iSpec)%InterID*SpecDSMC(iSpec)%TrefVHS*SpecDSMC(iSpec)%DrefVHS).eq.0) THEN
    CALL Abort(&
     __STAMP__,&
     "ERROR in species data ini_2")
    END IF
    SpecDSMC(iSpec)%omegaVHS = GETREAL('Part-Species'//TRIM(hilf)//'-omegaVHS','0') ! default case HS 
  ! reading electronic state informations from HDF5 file
    IF ( DSMC%ElectronicState .and. SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
      CALL ReadSpeciesLevel(SpecDSMC(iSpec)%Name,iSpec)
    END IF
  END DO

! species and case assignment arrays

  ALLOCATE(CollInf%Coll_Case(nSpecies,nSpecies))
  iCase = 0
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = iCase + 1
      CollInf%Coll_Case(iSpec,jSpec) = iCase
      CollInf%Coll_Case(jSpec,iSpec) = iCase
    END DO
  END DO
  nCase = iCase
  CollInf%NumCase = nCase 
  ALLOCATE(DSMC%NumColl(nCase +1))
  DSMC%NumColl = 0
  ALLOCATE(CollInf%Coll_CaseNum(nCase))
  CollInf%Coll_CaseNum = 0
  ALLOCATE(CollInf%Coll_SpecPartNum(nSpecies))
  CollInf%Coll_SpecPartNum = 0

  ALLOCATE(CollInf%FracMassCent(nSpecies, nCase)) ! Calculation of mx/(mx+my) and reduced mass
  CollInf%FracMassCent = 0
  ALLOCATE(CollInf%MassRed(nCase))
  CollInf%MassRed = 0
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      CollInf%FracMassCent(iSpec, iCase)  = Species(iSpec)%MassIC &
                                          / (Species(iSpec)%MassIC + Species(jSpec)%MassIC)
      CollInf%FracMassCent(jSpec, iCase)  = Species(jSpec)%MassIC &
                                          / (Species(iSpec)%MassIC + Species(jSpec)%MassIC)
      CollInf%MassRed(iCase) = (Species(iSpec)%MassIC*Species(jSpec)%MassIC) &
                             / (Species(iSpec)%MassIC+Species(jSpec)%MassIC)
    END DO
  END DO

! calculation of factors for pc

  ALLOCATE(CollInf%Cab(nCase))
  ALLOCATE(CollInf%KronDelta(nCase))
  CollInf%Cab = 0
  CollInf%KronDelta = 0
  
  DO iSpec = 1, nSpecies
    DO jSpec = iSpec, nSpecies
      iCase = CollInf%Coll_Case(iSpec,jSpec)
      IF (iSpec.eq.jSpec) THEN
        CollInf%KronDelta(iCase) = 1
      ELSE
        CollInf%KronDelta(iCase) = 0
      END IF
      A1 = 0.5 * SQRT(Pi) * SpecDSMC(iSpec)%DrefVHS*(2*(2-SpecDSMC(iSpec)%omegaVHS) &
            * BoltzmannConst * SpecDSMC(iSpec)%TrefVHS)**(SpecDSMC(iSpec)%omegaVHS*0.5)
      A2 = 0.5 * SQRT(Pi) * SpecDSMC(jSpec)%DrefVHS*(2*(2-SpecDSMC(jSpec)%omegaVHS) &
            * BoltzmannConst * SpecDSMC(jSpec)%TrefVHS)**(SpecDSMC(jSpec)%omegaVHS*0.5)
      CollInf%Cab(iCase) = (A1 + A2)**2 * ((Species(iSpec)%MassIC + Species(jSpec)%MassIC) &
            / (Species(iSpec)%MassIC * Species(jSpec)%MassIC))**SpecDSMC(iSpec)%omegaVHS 
            !the omega should be the same for both in vhs!!!
    END DO
  END DO

!-----------------------------------------------------------------------------------------------------------------------------------
! reading BG Gas stuff (required for the temperature definition in iInit=0)
!-----------------------------------------------------------------------------------------------------------------------------------
  BGGas%BGGasSpecies  = GETINT('Particles-DSMCBackgroundGas','0')
  BGGas%BGGasDensity  = GETREAL('Particles-DSMCBackgroundGasDensity','0')

  IF (useVTKFileBGG) THEN
    IF (Species(BGGas%BGGasSpecies)%Init(0)%velocityDistribution.NE.'maxwell_lpn') & !(use always Init 0 for BGG !!!)
    CALL abort(&
    __STAMP__&
    ,'only maxwell_lpn is implemened as velocity-distribution for BGG from VTK-File!')
    CALL ReadinVTKFileBGG()
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing molecular stuff
!-----------------------------------------------------------------------------------------------------------------------------------
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN ! perform relaxation (molecular) reactions
  ! allocate internal enery arrays
    IF ( DSMC%ElectronicState ) THEN
      ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,3))
      ALLOCATE(PartElecQua(PDM%maxParticleNumber))
      PartElecQua = 0
    ELSE
      ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,2))
    ENDIF
    ! reading molecular stuff
    SpecDSMC(1:nSpecies)%Xi_Rot = 0
    SpecDSMC(1:nSpecies)%MaxVibQuant = 0
    SpecDSMC(1:nSpecies)%CharaTVib = 0
    SpecDSMC(1:nSpecies)%PolyatomicMol=.false.
    SpecDSMC(1:nSpecies)%SpecToPolyArray = 0
    
    DO iSpec = 1, nSpecies
      IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
        WRITE(UNIT=hilf,FMT='(I2)') iSpec
        SpecDSMC(iSpec)%PolyatomicMol=GETLOGICAL('Part-Species'//TRIM(hilf)//'-PolyatomicMol','.FALSE.')
        IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%ElectronicState)  THEN
          CALL Abort(&
            __STAMP__,&
            '! Simulation of Polyatomic Molecules and Electronic States are not possible yet!!!')
        END IF
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
          SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
        ELSE
          SpecDSMC(iSpec)%Xi_Rot     = 2
          SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib','0.')
          SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')
          IF(SpecDSMC(iSpec)%Ediss_eV*SpecDSMC(iSpec)%CharaTVib.EQ.0) THEN
            CALL Abort(&
              __STAMP__,&
              'Error! Ediss_eV or CharaTVib is not set or equal to zero!')
          ELSE
            IF (DSMC%VibEnergyModel.EQ.0) THEN
              SpecDSMC(iSpec)%MaxVibQuant = 200
            ELSE
              SpecDSMC(iSpec)%MaxVibQuant = INT(SpecDSMC(iSpec)%Ediss_eV*JToEv/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)) + 1
            END IF
          END IF
        END IF
        SpecDSMC(iSpec)%VFD_Phi3_Factor = GETREAL('Part-Species'//TRIM(hilf)//'-VFDPhi3','0.')
        ! Setting the values of Rot-/Vib-RelaxProb to a fix value!
        ! This should be changed to a calculated value for every coll pair/situation!!!1
        SpecDSMC(iSpec)%RotRelaxProb  = 0.2     !0.2
        SpecDSMC(iSpec)%VibRelaxProb  = 0.02    !0.02
        ! multi init stuff
        ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
        DO iInit = 0, Species(iSpec)%NumberOfInits
          IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
            IF (iInit .EQ. 0) THEN !0. entry := old style parameter def. (default values if not def., some values might be needed)
              hilf2=TRIM(hilf)
            ELSE ! iInit >0
              WRITE(UNIT=hilf2,FMT='(I2)') iInit
              hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
            END IF ! iInit
            SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','0.')
            SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','0.')
            IF (SpecDSMC(iSpec)%Init(iInit)%TRot*SpecDSMC(iSpec)%Init(iInit)%TVib.EQ.0.) THEN
              IF (iInit.EQ.0)THEN
                IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                  CALL Abort(&
                    __STAMP__,&
                    'Error! TVib and TRot need to defined in Part-SpeciesXX-TempVib/TempRot for iSpec',iSpec)
                ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                  CALL Abort(&
                    __STAMP__,&
                    'Error! TVib and TRot need to defined in Part-SpeciesXX-TempVib/TempRot for BGGas')
                END IF
              ELSE ! iInit >0
                CALL Abort(&
                 __STAMP__,&
                 'Error! TVib and TRot need to defined in Part-SpeciesXX-InitXX-TempVib/TempRot for iSpec, iInit',iSpec,REAL(iInit))
              END IF
            END IF
          END IF
        END DO !Inits
      END IF !Molecule
    END DO !Species

    ! init of electronic state
    IF ( DSMC%ElectronicState ) THEN
      DO iSpec = 1, nSpecies
        ! electrons do not have an electron hull
        IF(SpecDSMC(iSpec)%InterID.EQ.4) CYCLE
        IF(.NOT.ALLOCATED(SpecDSMC(iSpec)%Init)) &
          ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
        WRITE(UNIT=hilf,FMT='(I2)') iSpec
        DO iInit = 0, Species(iSpec)%NumberOfInits
          IF (iInit .EQ. 0) THEN !0. entry := old style parameter def. (default values if not def., some values might be needed)
            hilf2=TRIM(hilf)
          ELSE ! iInit >0
            WRITE(UNIT=hilf2,FMT='(I2)') iInit
            hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
          END IF ! iInit
          SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-Tempelec','0.')
          IF (SpecDSMC(iSpec)%Init(iInit)%Telec.EQ.0.) THEN
            IF (iInit.EQ.0)THEN
              IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                CALL Abort(&
                  __STAMP__,&
                  ' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for Species',iSpec)
              ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                CALL Abort(&
                  __STAMP__,&
                  ' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for BGGas')
              END IF
            ELSE ! iInit >0
              CALL Abort(&
                __STAMP__,&
                ' Error! Telec needs to defined in Part-SpeciesXX-InitXX-Tempelc for iSpec, iInit',iSpec,REAL(iInit))
            END IF
          END IF
        END DO !Inits
        SpecDSMC(iSpec)%ElecRelaxProb = 0.01    !or 0.02 | Bird: somewhere in range 0.01 .. 0.02
      END DO ! iSpec
    END IF ! DSMC%ElectronicState

    
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      ALLOCATE(VibQuantsPar(PDM%maxParticleNumber))
      ALLOCATE(PolyatomMolDSMC(DSMC%NumPolyatomMolecs))
      DO iSpec = 1, nSpecies
        DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
          IF (SpecDSMC(iSpec)%PolyatomicMol) CALL InitPolyAtomicMolecs(iSpec,iInit)
        END DO
      END DO
    END IF
    
#if ( PP_TimeDiscMethod ==42 )
    IF ( DSMC%ElectronicState ) THEN
      DO iSpec = 1, nSpecies
        IF ( SpecDSMC(iSpec)%InterID .eq. 4) THEN
          SpecDSMC(iSpec)%MaxElecQuant = 0
        ELSE
          ALLOCATE( SpecDSMC(iSpec)%levelcounter         ( 0:size(SpecDSMC(iSpec)%ElectronicState,2)-1) , &
                    SpecDSMC(iSpec)%dtlevelcounter       ( 0:size(SpecDSMC(ispec)%ElectronicState,2)-1) , &
                    SpecDSMC(iSpec)%ElectronicTransition ( 1:nSpecies                                   , &
                                                          0:size(SpecDSMC(ispec)%ElectronicState,2)-1,   &
                                                          0:size(SpecDSMC(ispec)%ElectronicState,2)-1)   )
          SpecDSMC(iSpec)%levelcounter         = 0
          SpecDSMC(iSpec)%dtlevelcounter       = 0
          SpecDSMC(iSpec)%ElectronicTransition = 0
        END IF
      END DO
    END IF
#endif

    DO iPart = 1, PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (Species(PartSpecies(iPart))%NumberOfInits.EQ.0) THEN
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            iInit = PDM%PartInit(iPart)
            CALL DSMC_SetInternalEnr_PolyFastPart2(PartSpecies(iPart),iInit,iPart) ! wait for polyDSMC merge
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),0,iPart)
          END IF
        ELSE
          iInit = PDM%PartInit(iPart)
          IF(SpecDSMC(PartSpecies(iPart))%InterID.EQ.2) THEN
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iInit,iPart)
          END IF
        END IF
      END IF
    END DO
    
#if ( PP_TimeDiscMethod ==42 )
    ! Debug Output for initialized electronic state
    IF ( DSMC%ElectronicState ) THEN
      DO iSpec = 1, nSpecies
        print*,SpecDSMC(iSpec)%InterID
        IF ( SpecDSMC(iSpec)%InterID .ne. 4 ) THEN
          IF (  SpecDSMC(iSpec)%levelcounter(0) .ne. 0) THEN
            WRITE(DebugElectronicStateFilename,'(I2.2)') iSpec
            DebugElectronicStateFilename = 'Initial_Electronic_State_Species_'//trim(DebugElectronicStateFilename)//'.dat'
            open(unit=483,file=DebugElectronicStateFilename,form='formatted',status='unknown')
            DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant - 1
              WRITE(483,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / &
                                               REAL( Species(iSpec)%Init(0)%initialParticleNumber )
            END DO
            close(unit=483)
          END IF
        END IF
      END DO
    END IF
#endif

#if (PP_TimeDiscMethod!=1000) && (PP_TimeDiscMethod!=1001)
    DEALLOCATE(PDM%PartInit)
#endif
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! Define chemical Reactions
!-----------------------------------------------------------------------------------------------------------------------------------
  ChemReac%MeanEVib_Necc = .FALSE. ! this is only true, if a dis react is defined
  IF (CollisMode.EQ.3) THEN ! perform chemical reactions
    ! reading, molecular reaction stuff, ionization stuff,    
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I2)') iSpec
      SpecDSMC(iSpec)%Eion_eV               = GETREAL('Part-Species'//TRIM(hilf)//'-IonizationEn_eV','0')    
      SpecDSMC(iSpec)%RelPolarizability     = GETREAL('Part-Species'//TRIM(hilf)//'-RelPolarizability','0')
      SpecDSMC(iSpec)%NumEquivElecOutShell  = GETINT('Part-Species'//TRIM(hilf)//'-NumEquivElecOutShell','0')
      SpecDSMC(iSpec)%NumOfPro              = GETINT('Part-Species'//TRIM(hilf)//'-NumOfProtons','0')
      IF((SpecDSMC(iSpec)%Eion_eV*SpecDSMC(iSpec)%RelPolarizability*SpecDSMC(iSpec)%NumEquivElecOutShell &
              *SpecDSMC(iSpec)%NumOfPro).eq.0) THEN    
        SWRITE(*,*) "Ionization parameters are not defined for species:", iSpec     
      END IF
    END DO
    CALL DSMC_chemical_init()
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing Octree stuff
!-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','80')
  DSMC%UseOctree = GETLOGICAL('Particles-DSMC-UseOctree','.FALSE.')

!-----------------------------------------------------------------------------------------------------------------------------------
! Set mean VibQua of BGGas for dissoc reaction
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (BGGas%BGGasSpecies.NE.0) THEN
    IF(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.2) THEN
      BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib & 
                + BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib  &
                /  (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1) &
                - BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                / (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1)
      BGGasEVib = BGGasEVib/(BoltzmannConst*SpecDSMC(BGGas%BGGasSpecies)%CharaTVib) - DSMC%GammaQuant
      BGGas%BGMeanEVibQua = MIN(INT(BGGasEVib) + 1, SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant)    
    ELSE
      BGGas%BGMeanEVibQua = 0
    END IF
  END IF

  END IF !CollisMode.GT.0
!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing Surface stuff
!-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%CalcSurfaceVal = GETLOGICAL('Particles-DSMC-CalcSurfaceVal','.FALSE.')
  IF (DSMC%CalcSurfaceVal) THEN
  CALL DSMC_BuildSurfaceOutputMapping()
#ifdef MPI
  CALL DSMC_BuildHaloSurfaceOutputMapping()
#endif
  END IF
  ! Species-dependent calculations
  ALLOCATE(DSMC%CalcSurfCollis_SpeciesFlags(1:nSpecies))
  DSMC%CalcSurfCollis_NbrOfSpecies = GETINT('Particles-DSMC-CalcSurfCollis_NbrOfSpecies','0')
  IF ( (DSMC%CalcSurfCollis_NbrOfSpecies.GT.0) .AND. (DSMC%CalcSurfCollis_NbrOfSpecies.LE.nSpecies) ) THEN
    ALLOCATE(CalcSurfCollis_SpeciesRead(1:DSMC%CalcSurfCollis_NbrOfSpecies))
    hilf2=''
    DO iSpec=1,DSMC%CalcSurfCollis_NbrOfSpecies !build default string: 1 - CSC_NoS
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      hilf2=TRIM(hilf2)//TRIM(hilf)
      IF (ispec.NE.DSMC%CalcSurfCollis_NbrOfSpecies) hilf2=TRIM(hilf2)//','
    END DO
    CalcSurfCollis_SpeciesRead = GETINTARRAY('Particles-DSMC-CalcSurfCollis_Species',DSMC%CalcSurfCollis_NbrOfSpecies,hilf2)
    DSMC%CalcSurfCollis_SpeciesFlags(:)=.FALSE.
    DO iSpec=1,DSMC%CalcSurfCollis_NbrOfSpecies
      DSMC%CalcSurfCollis_SpeciesFlags(CalcSurfCollis_SpeciesRead(ispec))=.TRUE.
    END DO
    DEALLOCATE(CalcSurfCollis_SpeciesRead)
  ELSE IF (DSMC%CalcSurfCollis_NbrOfSpecies.EQ.0) THEN !default
    DSMC%CalcSurfCollis_SpeciesFlags(:)=.TRUE.
  ELSE
    CALL abort(__STAMP__,&
              'Error in Particles-DSMC-CalcSurfCollis_NbrOfSpecies!')
  END IF
  DSMC%CalcSurfCollis_OnlySwaps = GETLOGICAL('Particles-DSMC-CalcSurfCollis_OnlySwaps','.FALSE.')
  DSMC%CalcSurfCollis_Only0Swaps = GETLOGICAL('Particles-DSMC-CalcSurfCollis_Only0Swaps','.FALSE.')
  DSMC%CalcSurfCollis_Output = GETLOGICAL('Particles-DSMC-CalcSurfCollis_Output','.FALSE.')
  IF (DSMC%CalcSurfCollis_Only0Swaps) DSMC%CalcSurfCollis_OnlySwaps=.TRUE.

  SWRITE(UNIT_stdOut,'(A)')' INIT DSMC DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitDSMC


SUBROUTINE DSMC_SetInternalEnr_LauxVFD(iSpecies, iInit, iPart)
!===================================================================================================================================
! Energy distribution according to dissertation of Laux
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC
#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001)
  USE MOD_DSMC_Vars,            ONLY : LD_MultiTemperaturMod
#endif
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
  INTEGER, INTENT(IN)           :: iSpecies, iInit, iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan
  INTEGER                       :: iQuant
!===================================================================================================================================

  ! set vibrational energy
  IF (SpecDSMC(iSpecies)%InterID.EQ.2) THEN
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*SpecDSMC(iSpecies)%Init(iInit)%TVib/SpecDSMC(iSpecies)%CharaTVib)
    DO WHILE (iQuant.GE.SpecDSMC(iSpecies)%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*SpecDSMC(iSpecies)%Init(iInit)%TVib/SpecDSMC(iSpecies)%CharaTVib)
    END DO
    !evtl muss partstateinten nochmal geändert werden, mpi, resize etc..
    PartStateIntEn(iPart, 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpecies)%CharaTVib*BoltzmannConst
  ! set rotational energy
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpecies)%Init(iInit)%TRot*LOG(iRan)
  ELSE
    PartStateIntEn(iPart, 1) = 0
    PartStateIntEn(iPart, 2) = 0
  END IF
  ! set electronic energy
  IF ( DSMC%ElectronicState .AND. SpecDSMC(iSpecies)%InterID .NE. 4) THEN
    CALL InitElectronShell(iSpecies,iPart,iInit)
  ENDIF
#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001)
  IF (LD_MultiTemperaturMod .EQ. 3 ) THEN ! no discret vib levels for this LD method
    IF (SpecDSMC(iSpecies)%InterID.EQ.2) THEN
      PartStateIntEn(iPart, 1) = (BoltzmannConst*SpecDSMC(iSpecies)%CharaTVib) &
                               / (EXP(SpecDSMC(iSpecies)%CharaTVib/SpecDSMC(iSpecies)%Init(iInit)%TVib) - 1.0) &
                               + DSMC%GammaQuant * BoltzmannConst*SpecDSMC(iSpecies)%CharaTVib
    !set rotational energy
      PartStateIntEn(iPart, 2) = BoltzmannConst*SpecDSMC(iSpecies)%Init(iInit)%TRot
    ELSE
      PartStateIntEn(iPart, 1) = 0
      PartStateIntEn(iPart, 2) = 0
    END IF
  END IF
#endif

END SUBROUTINE DSMC_SetInternalEnr_LauxVFD


SUBROUTINE DSMC_BuildSurfaceOutputMapping()
!===================================================================================================================================
! Perform mapping for surface output
!===================================================================================================================================
! MODULES
 USE MOD_Globals
! USE MOD_Preproc
! USE MOD_Particle_Tracking_vars, ONLY:DoRefMapping
! USE MOD_Mesh_Vars,              ONLY:SurfElem
! USE MOD_DSMC_Vars,              ONLY:nSurfSample,XiEQ_Surf,deltaXiEQ_Surf
! USE MOD_Basis,                  ONLY:BarycentricWeights,InitializeVandermonde
! USE MOD_Interpolation_Vars,     ONLY:xGP,wBary
! USE MOD_ChangeBasis,            ONLY:ChangeBasis2D
! USE MOD_Particle_Mesh_Vars,     ONLY:PartBCSideList,nTotalBCSides,nTotalSides

!  USE MOD_Mesh_Vars,          ONLY:nBCSides, SideToElem, BC
!  USE MOD_Particle_Vars,      ONLY:nSpecies
!  USE MOD_Particle_Mesh_Vars, ONLY:GEO,PartBound
!  USE MOD_DSMC_Vars,          ONLY:SurfMesh, SampWall
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! INTEGER                   :: p,q, iSide
! REAL,ALLOCATABLE          :: VDM_NGP_NSurfEQ(:,:)!, wBary_NSurfSample(:)



!  INTEGER                 :: iElem, iLocSide, iSide, iNode, iNode2, iSampWallAlloc
!  INTEGER, ALLOCATABLE    :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:)
!  REAL,ALLOCATABLE        :: TempSurfaceArea(:)
!  LOGICAL                 :: IsSortedSurfNode
!===================================================================================================================================

!IF(.NOT.DoRefMapping) CALL abort(__STAMP__,&
!    ' Analyze Surfaces requires DoRefMapping=.TRUE.')
!
!nSurfSample = GETINT('DSMC-nSurfSample','0')
!ALLOCATE(XiEQ_Surf(0:nSurfSample))
!IF(nSurfSample.EQ.0)THEN
!  XiEQ_Surf=0.
!  DeltaXiEQ_Surf=2.
!ELSE
!  DO q=0,nSurfSample
!    XiEQ_Surf(q) = 2./REAL(nSurfSample) * REAL(q) - 1. 
!  END DO
!  DeltaXi_NGeo=2./NGeo_in
!ENDIF
!
!ALLOCATE(Vdm_NGP_NSurfEQ(0:NSurfSample,0:PP_N) & 
!        ,wBary_NSurfSample(0:NSurfSample)      &
!        ,SurfMesh%SideToSurfID(1:nTotalSides)  )
!
!!CALL BarycentricWeights(NSurfSample,XiEQ_Surf,wBary_NSurfSample)
!CALL InitializeVandermonde(PP_N,NSurfSample,wBary,xGP,XiEQ_Surf,Vdm_NGP_NSurfEQ)
!
!! first, get number of bc sides
!SurfMesh%nSurfaceBCSides = 0
!DO iSide=1,nTotalSides
!  IF(BC(iSide).LE.1) CYCLE
!  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
!    SurfMesh%nSurfaceBCSides = SurfMesh%nSurfaceBCSides + 1
!    SurfMesh%SideToSurfID(iSide) = SurfMesh%nSurfaceBCSides
!  END IF
!END DO ! iSide=1,nTotalSides
!
!!CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:),tmp)
!
!
!!
!!  ! ChangeBasis3D to lower or higher polynomial degree
!!  DO iElem=1,PP_nElems
!!    CALL ChangeBasis3D(BGDataSize,N_In,NBG,Vdm_BGFieldIn_BGField,BGField_tmp(:,:,:,:,iElem),BGField(:,:,:,:,iElem))
!!  END DO ! iElem
!
!
!
!DEALLOCATE(Vdm_NGP_NSurfEQ, wBary_NSurfSample)
!

      !  STOP

    CALL abort(&
    __STAMP__&
    ,' Subroutine not implemented!')

!  ALLOCATE(TempBCSurfNodes(4*nBCSides))
!  ALLOCATE(TempSideSurfNodeMap(1:4,1:nBCSides))
!  ALLOCATE(SurfMesh%GlobSideToSurfSideMap(nBCSides))
!  ALLOCATE(TempSurfaceArea(nBCSides))
!  SurfMesh%nSurfaceNode=0
!  SurfMesh%nSurfaceBCSides=0
!  SurfMesh%GlobSideToSurfSideMap(1:nBCSides)=0
!
!  DO iSide=1, nBCSides
!    IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
!      SurfMesh%nSurfaceBCSides = SurfMesh%nSurfaceBCSides + 1
!      SurfMesh%GlobSideToSurfSideMap(iSide) = SurfMesh%nSurfaceBCSides
!      iElem = SideToElem(1,iSide)
!      IF (iElem.LT.1) THEN
!        iElem = SideToElem(2,iSide)
!        iLocSide = SideToElem(4,iSide)
!      ELSE
!        iLocSide = SideToElem(3,iSide)
!      END IF
!      TempSurfaceArea(SurfMesh%nSurfaceBCSides) = CalcArea(iLocSide, iElem)
!      DO iNode2 = 1, 4
!      IsSortedSurfNode = .false.
!        DO iNode = 1, SurfMesh%nSurfaceNode 
!          IF (GEO%ElemSideNodeID(iNode2, iLocSide, iElem).EQ.TempBCSurfNodes(iNode)) THEN
!          TempSideSurfNodeMap(iNode2,SurfMesh%nSurfaceBCSides) = iNode
!          IsSortedSurfNode = .true.
!          EXIT
!          END IF
!        END DO
!        IF(.NOT.IsSortedSurfNode) THEN
!          SurfMesh%nSurfaceNode = SurfMesh%nSurfaceNode + 1
!          TempBCSurfNodes(SurfMesh%nSurfaceNode) = GEO%ElemSideNodeID(iNode2, iLocSide, iElem)
!          TempSideSurfNodeMap(iNode2,SurfMesh%nSurfaceBCSides) = SurfMesh%nSurfaceNode
!        END IF
!      END DO  
!    END IF
!  END DO

!  ALLOCATE(SurfMesh%BCSurfNodes(1:SurfMesh%nSurfaceNode))
!  SurfMesh%BCSurfNodes(1:SurfMesh%nSurfaceNode) = TempBCSurfNodes(1:SurfMesh%nSurfaceNode)
!  ALLOCATE(SurfMesh%SideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides))
!  SurfMesh%SideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides)
!  ! still required !!!
!  ALLOCATE(SurfMesh%SurfaceArea(1:SurfMesh%nSurfaceBCSides))
!  SurfMesh%SurfaceArea(1:SurfMesh%nSurfaceBCSides)=TempSurfaceArea(1:SurfMesh%nSurfaceBCSides)
!  ALLOCATE(SampWall(1:SurfMesh%nSurfaceBCSides))
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(1) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(2) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(3) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(4) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(5) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(6) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(7) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(8) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(9) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Force(1) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Force(2) = 0.0
!  SampWall(1:SurfMesh%nSurfaceBCSides)%Force(3) = 0.0
!  DO iSampWallAlloc=1,SurfMesh%nSurfaceBCSides
!    ALLOCATE(SampWall(iSampWallAlloc)%Counter(1:nSpecies))
!    SampWall(iSampWallAlloc)%Counter(1:nSpecies) = 0.0
!  END DO
!  DEALLOCATE(TempBCSurfNodes)
!  DEALLOCATE(TempSideSurfNodeMap)
!  DEALLOCATE(TempSurfaceArea)

END SUBROUTINE DSMC_BuildSurfaceOutputMapping


#ifdef MPI
SUBROUTINE DSMC_BuildHaloSurfaceOutputMapping()
!===================================================================================================================================
! Perform mapping for halo surface output of MPI case
!===================================================================================================================================
! MODULES
USE MOD_Globals
!  USE MOD_Particle_Vars,      ONLY : nSpecies
!  USE MOD_DSMC_Vars,          ONLY : SampWallHaloCell
!  USE MOD_DSMC_Vars,          ONLY : SurfMesh
!  USE MOD_Particle_Mesh_Vars, ONLY : PartBound
  !USE MOD_part_MPI_Vars,      ONLY : MPIGEO
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!  INTEGER                 :: iSide, iSampWallAlloc
!===================================================================================================================================

    CALL abort(&
    __STAMP__&
    ,' Subroutine not implemented!')
  !ALLOCATE(SurfMesh%HaloSideIDToSurfSideMap(SIZE(MPIGEO%BC,2)))
  !SurfMesh%HaloSideIDToSurfSideMap(1:SIZE(MPIGEO%BC,2))=0
  !SurfMesh%nHaloSurfaceBCSides = 0

  !DO iSide=1, SIZE(MPIGEO%BC,2)
  !  IF ((MPIGEO%BC(1,iSide).NE.0).AND.(MPIGEO%BC(1,iSide).NE.-1).AND.(MPIGEO%BC(1,iSide).NE.424242)) THEN
  !    IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(MPIGEO%BC(1,iSide))).EQ.PartBound%ReflectiveBC) THEN
  !      SurfMesh%nHaloSurfaceBCSides = SurfMesh%nHaloSurfaceBCSides + 1
  !      SurfMesh%HaloSideIDToSurfSideMap(iSide) = SurfMesh%nHaloSurfaceBCSides
  !    END IF
  !  END IF
  !END DO

  !ALLOCATE(SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides))
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(1) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(2) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(3) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(4) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(5) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(6) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(7) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(8) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(9) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(1) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(2) = 0.0
  !SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(3) = 0.0
  !DO iSampWallAlloc=1,SurfMesh%nHaloSurfaceBCSides
  !  ALLOCATE(SampWallHaloCell(iSampWallAlloc)%Counter(1:nSpecies))
  !  SampWallHaloCell(iSampWallAlloc)%Counter(1:nSpecies) = 0.0
  !END DO

END SUBROUTINE DSMC_BuildHaloSurfaceOutputMapping
#endif


SUBROUTINE ReadinVTKFileBGG
!===================================================================================================================================
! Readin of custom VTK file for non-constant background gas distribution
!===================================================================================================================================
! MODULES
  USE MOD_Globals
!  USE MOD_ReadInTools
#ifndef MPI
!  USE MOD_Mesh_Vars,              ONLY : nNodes
#endif
!  USE MOD_Mesh_Vars,              ONLY : nElems
!  USE MOD_Particle_Vars,          ONLY : BGGdataAtElem
!  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!  INTEGER                :: unit_in
!  INTEGER                :: os !openStatus
!  CHARACTER(255)         :: VTKfile
!  CHARACTER(LEN=255)     :: cdummy, varname
!  INTEGER                :: npoints, iNode, ncells, icell, VI, icoord, iElem
!  REAL, ALLOCATABLE      :: VTKNodes(:,:), VTK_BGGdata_Cells(:,:),VTKCellsSP(:,:)
!  LOGICAL, ALLOCATABLE   :: IsAssociated1(:),IsAssociated2(:)
!  REAL                   :: x(3), Elem_SP(3)
!  INTEGER                :: iVTKcell, CellX, CellY, CellZ
!  INTEGER, ALLOCATABLE  :: VTKCells(:,:)
!  LOGICAL                :: InElementCheck
!===================================================================================================================================

    CALL abort(&
    __STAMP__&
    ,' Subroutine not implemented!')


!  SWRITE(UNIT_stdOut,'(132("~"))')
!  SWRITE(UNIT_stdOut,'(A)')'Reading VTK file for BGG...'
!
!  VTKfile = GETSTR('BGG-VTK-File','blubb')
!  IF(TRIM(VTKfile).EQ.'blubb')THEN 
!    CALL abort(__STAMP__,&
!    'ERROR: No VTK-Filename for Background-Gas defined!')
!  END IF 
!
!  unit_in = 1123
!  OPEN(UNIT   = unit_in,              &
!       FILE   = VTKfile,              &
!       IOSTAT = os,                   &
!       STATUS = 'OLD',                &
!       ACTION = 'READ',               &
!       ACCESS = 'SEQUENTIAL'          )
!
!  IF(os.NE.0) THEN  ! File Error
!    CALL abort(__STAMP__,&
!    'ERROR: cannot open VTK file: '//trim(VTKfile))
!  END IF
!    
!  DO iNode=1,7
!    READ(unit_in, '(A)') cdummy
!  END DO
!  READ(unit_in,*) cdummy,npoints,cdummy  ! POINTS ???? float
!  ALLOCATE(VTKNodes(1:3,npoints))
!  DO iNode = 0,INT(npoints/3)-1
!    READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2),VTKNodes(:,3*iNode+3)
!  END DO
!  IF (MOD(npoints,3).EQ.1) THEN
!    READ(unit_in,*) VTKNodes(:,3*iNode+1)
!  ELSE IF (MOD(npoints,3).EQ.2) THEN
!    READ(unit_in,*) VTKNodes(:,3*iNode+1),VTKNodes(:,3*iNode+2)
!  END IF
!  READ(unit_in,*) cdummy,ncells,cdummy  ! CELLS ???? ????
!  ALLOCATE (VTKCells(1:8, ncells))
!  DO icell = 1,ncells
!    READ(unit_in,*), cdummy, VTKCells(:,icell)
!  END DO
!  VTKCells(:,:) = VTKCells(:,:) + 1
!  READ(unit_in, '(A)') cdummy  ! blank line
!  READ(unit_in, '(A)') cdummy  ! blank line
!  DO icell = 1,ncells
!    READ(unit_in,*) cdummy !skip cells
!  END DO
!  !var1
!  DO icell=1,2
!    READ(unit_in, '(A)') cdummy
!  END DO
!  READ(unit_in,*) cdummy, varname, cdummy
!  VI = 0
!  IF (TRIM(varname).EQ.'T') THEN
!    VI=0
!  ELSE IF (TRIM(varname).EQ.'v') THEN
!    VI=3
!  ELSE
!    CALL abort(__STAMP__,&
!    'Error in order of variables in BGGdata: '//varname//' is not recognized!')
!  END IF
!  ALLOCATE(VTK_BGGdata_Cells(1:7,ncells))
!  DO icell = 0,INT(ncells/3)-1
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1),VTK_BGGdata_Cells(1+VI:3+VI,3*icell+2) &
!      ,VTK_BGGdata_Cells(1+VI:3+VI,3*icell+3)
!  END DO
!  IF (MOD(ncells,3).EQ.1) THEN
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1)
!  ELSE IF (MOD(ncells,3).EQ.2) THEN
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1),VTK_BGGdata_Cells(1+VI:3+VI,3*icell+2)
!  END IF
!  !var2
!  DO icell=1,2
!    READ(unit_in, '(A)') cdummy
!  END DO
!  READ(unit_in,*) varname, cdummy, cdummy, cdummy
!  IF (TRIM(varname).EQ.'T') THEN
!    VI=0
!  ELSE IF (TRIM(varname).EQ.'v') THEN
!    VI=3
!  ELSE
!    CALL abort(__STAMP__,&
!    'Error in order of variables in BGGdata: '//varname//' is not recognized!')
!  END IF
!  DO icell = 0,INT(ncells/3)-1
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1),VTK_BGGdata_Cells(1+VI:3+VI,3*icell+2) &
!      ,VTK_BGGdata_Cells(1+VI:3+VI,3*icell+3)
!  END DO
!  IF (MOD(ncells,3).EQ.1) THEN
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1)
!  ELSE IF (MOD(ncells,3).EQ.2) THEN
!    READ(unit_in,*) VTK_BGGdata_Cells(1+VI:3+VI,3*icell+1),VTK_BGGdata_Cells(1+VI:3+VI,3*icell+2)
!  END IF
!  !n
!  DO icell=1,2
!    READ(unit_in, '(A)') cdummy
!  END DO
!  DO icell = 0,INT(ncells/9)-1
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4),VTK_BGGdata_Cells(7,9*icell+5),VTK_BGGdata_Cells(7,9*icell+6) &
!      ,VTK_BGGdata_Cells(7,9*icell+7),VTK_BGGdata_Cells(7,9*icell+8),VTK_BGGdata_Cells(7,9*icell+9)
!  END DO
!  SELECT CASE (MOD(ncells,9))
!  CASE(1)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1)
!  CASE(2)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2)
!  CASE(3)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3)
!  CASE(4)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4)
!  CASE(5)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4),VTK_BGGdata_Cells(7,9*icell+5)
!  CASE(6)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4),VTK_BGGdata_Cells(7,9*icell+5),VTK_BGGdata_Cells(7,9*icell+6)
!  CASE(7)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4),VTK_BGGdata_Cells(7,9*icell+5),VTK_BGGdata_Cells(7,9*icell+6) &
!      ,VTK_BGGdata_Cells(7,9*icell+7)
!  CASE(8)
!    READ(unit_in,*) VTK_BGGdata_Cells(7,9*icell+1),VTK_BGGdata_Cells(7,9*icell+2),VTK_BGGdata_Cells(7,9*icell+3) &
!      ,VTK_BGGdata_Cells(7,9*icell+4),VTK_BGGdata_Cells(7,9*icell+5),VTK_BGGdata_Cells(7,9*icell+6) &
!      ,VTK_BGGdata_Cells(7,9*icell+7),VTK_BGGdata_Cells(7,9*icell+8)
!  END SELECT
!  CLOSE(1123)
!
!#ifndef MPI
!  IF (npoints.NE.nNodes) THEN
!    CALL abort(__STAMP__,&
!    'ERROR: wrong number of points in VTK-File')
!  END IF
!  IF (ncells.NE.nElems) THEN
!    CALL abort(__STAMP__,&
!    'ERROR: wrong number of cells in VTK-File')
!  END IF
!#endif /*MPI*/
!
!  ALLOCATE(VTKCellsSP(1:3,ncells))
!  DO icell = 1, ncells
!    Elem_SP(:)=0.
!    DO iNode = 1,8 !SP Nodes_new
!      DO icoord = 1,3
!        Elem_SP(icoord)=Elem_SP(icoord)+VTKNodes(icoord,VTKCells(iNode,icell))
!      END DO
!    END DO
!    VTKCellsSP(:,icell)=Elem_SP(:)/8.
!  END DO
!
!  ALLOCATE(BGGdataAtElem(1:7,nElems))
!  ALLOCATE(IsAssociated1(ncells))
!  ALLOCATE(IsAssociated2(nElems))
!  BGGdataAtElem=0.
!  IsAssociated1=.FALSE.
!  IsAssociated2=.FALSE.
!
!  DO iVTKcell = 1,ncells
!    x = VTKCellsSP(1:3,iVTKcell)
!    CellX = INT((x(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
!    CellY = INT((x(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
!    CellZ = INT((x(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!#ifdef MPI
!    IF ((GEO%FIBGMimax.GE.CellX).AND.(GEO%FIBGMimin.LE.CellX)) THEN  
!    IF ((GEO%FIBGMjmax.GE.CellY).AND.(GEO%FIBGMjmin.LE.CellY)) THEN  
!    IF ((GEO%FIBGMkmax.GE.CellZ).AND.(GEO%FIBGMkmin.LE.CellZ)) THEN  
!#endif /* MPI */
!      DO iElem=1, nElems
!        ! OMG
!        CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(1,iElem)), &
!              GEO%NodeCoords(:,GEO%ElemToNodeID(4,iElem)), &
!              GEO%NodeCoords(:,GEO%ElemToNodeID(3,iElem)), &
!              GEO%NodeCoords(:,GEO%ElemToNodeID(2,iElem)), &
!              VTKCellsSP(:,iVTKcell), &
!              InElementCheck)
!        IF(InElementCheck) THEN
!          CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(3,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(7,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(6,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(2,iElem)), &
!                VTKCellsSP(:,iVTKcell), &
!                InElementCheck)
!        END IF
!        IF(InElementCheck) THEN
!          CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(6,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(5,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(1,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(2,iElem)), &
!                VTKCellsSP(:,iVTKcell), &
!                InElementCheck)
!        END IF
!        IF(InElementCheck) THEN
!          CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(5,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(8,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(4,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(1,iElem)), &
!                VTKCellsSP(:,iVTKcell), &
!                InElementCheck)
!        END IF
!        IF(InElementCheck) THEN
!          CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(8,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(7,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(3,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(4,iElem)), &
!                VTKCellsSP(:,iVTKcell), &
!                InElementCheck)
!        END IF
!        IF(InElementCheck) THEN
!          CALL ParticleInsideQuad3D_DSMC(GEO%NodeCoords(:,GEO%ElemToNodeID(5,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(6,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(7,iElem)), &
!                GEO%NodeCoords(:,GEO%ElemToNodeID(8,iElem)), &
!                VTKCellsSP(:,iVTKcell), &
!                InElementCheck)
!        END IF
!        IF(InElementCheck) THEN
!          IF (IsAssociated1(iVTKcell).OR.IsAssociated2(iElem)) THEN
!            CALL abort(__STAMP__,&
!              'ERROR: Cell is already mapped!')
!          END IF
!          BGGdataAtElem(:,iElem)=VTK_BGGdata_Cells(:,iVTKcell)
!          IsAssociated1(iVTKcell)=.TRUE.
!          IsAssociated2(iElem)=.TRUE.
!          EXIT
!        END IF
!      END DO
!#ifdef MPI
!    END IF
!    END IF 
!    END IF  
!#endif /* MPI */
!  END DO ! iVTKcell
!  !IF (.NOT. ALL(IsAssociated1)) THEN         !only for 1 proc!!!!!!!!
!  !  CALL abort(__STAMP__,&
!  !	'ERROR: Not all VTKcells mapped for BGG!',999,999.)
!  !END IF
!  IF (.NOT. ALL(IsAssociated2)) THEN
!    CALL abort(__STAMP__,&
!              'ERROR: Not all Elems mapped for BGG!')
!  END IF
!  SWRITE(UNIT_stdOut,'(A)')'DONE!'
!  SWRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE ReadinVTKFileBGG


!REAL FUNCTION CalcArea(iLocSide, Element)
!===================================================================================================================================
! Calculates area of mesh element
!===================================================================================================================================
! MODULES
!  USE MOD_Particle_Vars,          ONLY : GEO
!  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
!
! IMPLICIT VARIABLE HANDLING
!  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!  INTEGER, INTENT(IN)         :: iLocSide, Element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!  REAL                        :: xNod1, xNod2, xNod3, xNod4 
!  REAL                        :: yNod1, yNod2, yNod3, yNod4 
!  REAL                        :: zNod1, zNod2, zNod3, zNod4 
!  REAL                        :: Vector1(1:3), Vector2(1:3), Vector3(1:3)
!===================================================================================================================================
!
!CalcArea=1.
!STOP
!!  xNod1 = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
!!  yNod1 = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
!!  zNod1 = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
!!
!!  xNod2 = GEO%NodeCoords(1,GEO%ElemSideNodeID(2,iLocSide,Element))
!!  yNod2 = GEO%NodeCoords(2,GEO%ElemSideNodeID(2,iLocSide,Element))
!!  zNod2 = GEO%NodeCoords(3,GEO%ElemSideNodeID(2,iLocSide,Element))
!!
!!  xNod3 = GEO%NodeCoords(1,GEO%ElemSideNodeID(3,iLocSide,Element))
!!  yNod3 = GEO%NodeCoords(2,GEO%ElemSideNodeID(3,iLocSide,Element))
!!  zNod3 = GEO%NodeCoords(3,GEO%ElemSideNodeID(3,iLocSide,Element))
!!
!!  xNod4 = GEO%NodeCoords(1,GEO%ElemSideNodeID(4,iLocSide,Element))
!!  yNod4 = GEO%NodeCoords(2,GEO%ElemSideNodeID(4,iLocSide,Element))
!!  zNod4 = GEO%NodeCoords(3,GEO%ElemSideNodeID(4,iLocSide,Element))
!!
!!  Vector1(1) = xNod2 - xNod1
!!  Vector1(2) = yNod2 - yNod1
!!  Vector1(3) = zNod2 - zNod1
!!
!!  Vector2(1) = xNod3 - xNod1
!!  Vector2(2) = yNod3 - yNod1
!!  Vector2(3) = zNod3 - zNod1
!!
!!  Vector3(1) = xNod4 - xNod1
!!  Vector3(2) = yNod4 - yNod1
!!  Vector3(3) = zNod4 - zNod1
!!
!!  CalcArea = 0.5*(SQRT((Vector1(2)*Vector2(3)-Vector1(3)*Vector2(2))**2 &
!!         + (-Vector1(1)*Vector2(3)+Vector1(3)*Vector2(1))**2 &
!!         + (Vector1(1)*Vector2(2)-Vector1(2)*Vector2(1))**2) &
!!         + SQRT((Vector3(2)*Vector2(3)-Vector3(3)*Vector2(2))**2 &
!!         + (-Vector3(1)*Vector2(3)+Vector3(3)*Vector2(1))**2 &
!!         + (Vector3(1)*Vector2(2)-Vector3(2)*Vector2(1))**2))
!!
!!  RETURN
!
!END FUNCTION CalcArea


SUBROUTINE FinalizeDSMC() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize dsmc variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(SampDSMC)
SDEALLOCATE(DSMC_RHS)
SDEALLOCATE(PartStateIntEn)
SDEALLOCATE(SpecDSMC)
SDEALLOCATE(DSMC%NumColl)
SDEALLOCATE(DSMC%CalcSurfCollis_SpeciesFlags)
SDEALLOCATE(DSMC%CollProbOut)
SDEALLOCATE(DSMC%CollProbSamp)
SDEALLOCATE(Coll_pData)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(MacroDSMC)
SDEALLOCATE(ChemReac%QKProcedure)
SDEALLOCATE(ChemReac%QKMethod)
SDEALLOCATE(ChemReac%QKCoeff)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReacCount)
SDEALLOCATE(ChemReac%NumReac)
SDEALLOCATE(ChemReac%ReactType)
SDEALLOCATE(ChemReac%DefinedReact)
SDEALLOCATE(ChemReac%ReactCase)
SDEALLOCATE(ChemReac%ReactNum)
SDEALLOCATE(ChemReac%Arrhenius_Prefactor)
SDEALLOCATE(ChemReac%Arrhenius_Powerfactor)
SDEALLOCATE(ChemReac%EActiv)
SDEALLOCATE(ChemReac%EForm)
SDEALLOCATE(ChemReac%MeanEVib_PerIter)
SDEALLOCATE(ChemReac%MeanEVibQua_PerIter)
SDEALLOCATE(ChemReac%CEXa)
SDEALLOCATE(ChemReac%CEXb)
SDEALLOCATE(ChemReac%MEXa)
SDEALLOCATE(ChemReac%MEXb)
SDEALLOCATE(ChemReac%ReactInfo)
SDEALLOCATE(CollInf%Coll_Case)
SDEALLOCATE(CollInf%Coll_CaseNum)
SDEALLOCATE(CollInf%Coll_SpecPartNum)
SDEALLOCATE(CollInf%Cab)
SDEALLOCATE(CollInf%KronDelta)
SDEALLOCATE(CollInf%FracMassCent)
SDEALLOCATE(CollInf%MassRed)
SDEALLOCATE(SampWall)
SDEALLOCATE(MacroSurfaceVal)
SDEALLOCATE(VibQuantsPar)
SDEALLOCATE(XiEq_Surf)
SDEALLOCATE(DSMC_HOSolution)

END SUBROUTINE FinalizeDSMC


END MODULE MOD_DSMC_Init
