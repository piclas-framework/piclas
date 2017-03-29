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
PUBLIC :: InitDSMC, DSMC_SetInternalEnr_LauxVFD, FinalizeDSMC
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDSMC()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc,                    ONLY : PP_N
USE MOD_Mesh_Vars,                  ONLY : nElems, NGEo
USE MOD_Globals_Vars,               ONLY : Pi
USE MOD_ReadInTools
USE MOD_DSMC_ElectronicModel,       ONLY: ReadSpeciesLevel
USE MOD_DSMC_Vars
USE MOD_PARTICLE_Vars,              ONLY: nSpecies, BoltzmannConst, Species, PDM, PartSpecies, useVTKFileBGG
USE MOD_DSMC_Analyze,               ONLY: InitHODSMC
USE MOD_DSMC_ParticlePairing,       ONLY: DSMC_init_octree
USE MOD_DSMC_SteadyState,           ONLY: DSMC_SteadyStateInit
USE MOD_TimeDisc_Vars,              ONLY: TEnd
USE MOD_DSMC_ChemInit,              ONLY: DSMC_chemical_init
USE MOD_DSMC_SurfModelInit,         ONLY: InitDSMCSurfModel
USE MOD_DSMC_ChemReact,             ONLY: CalcBackwardRate, CalcPartitionFunction
USE MOD_DSMC_PolyAtomicModel,       ONLY: InitPolyAtomicMolecs, DSMC_FindFirstVibPick, DSMC_SetInternalEnr_Poly
USE MOD_Particle_Boundary_Sampling, ONLY: InitParticleBoundarySampling
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(32)         :: hilf , hilf2
  INTEGER               :: iCase, iSpec, jSpec, nCase, iPart, iInit, iPolyatMole, iDOF, PartitionArraySize
  INTEGER               :: iInter
  REAL                  :: A1, A2     ! species constant for cross section (p. 24 Laux)
  REAL                  :: JToEv, Temp
  INTEGER,ALLOCATABLE   :: CalcSurfCollis_SpeciesRead(:) !help array for reading surface stuff
  REAL                  :: BGGasEVib, Qtra, Qrot, Qvib, Qelec
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
!  IF (DSMC%OutputMeshInit) THEN
!    SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH...'
!    CALL WriteOutputMesh()
!    SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH DONE!'
!  END IF
! reading and reset general DSMC values
  CollisMode = GETINT('Particles-DSMC-CollisMode','1') !0: no collis, 1:elastic col, 2:elast+rela, 3:chem
  SelectionProc = GETINT('Particles-DSMC-SelectionProcedure','1') !1: Laux, 2:Gimelsheim
  DSMC%RotRelaxProb = GETREAL('Particles-DSMC-RotRelaxProb','0.2')  
  DSMC%VibRelaxProb = GETREAL('Particles-DSMC-VibRelaxProb','0.02')
  DSMC%ElecRelaxProb = GETREAL('Particles-DSMC-ElecRelaxProb','0.01')
  DSMC%GammaQuant   = GETREAL('Particles-DSMC-GammaQuant', '0.5')
!-----------------------------------------------------------------------------------
! Flag for the automatic calculation of the backward reaction rate with the partition functions and equilibrium constant.
! Partition functions are calculated for each species during initialization and stored for values starting with the
! DSMC%PartitionInterval up to DSMC%PartitionMaxTemp, interpolation between the stored values
  DSMC%BackwardReacRate = GETLOGICAL('Particles-DSMC-BackwardReacRate','.FALSE.')
  DSMC%PartitionMaxTemp = GETREAL('Particles-DSMC-PartitionMaxTemp','20000')
  DSMC%PartitionInterval = GETREAL('Particles-DSMC-PartitionInterval','10')
!-----------------------------------------------------------------------------------
  DSMC%TimeFracSamp = GETREAL('Part-TimeFracForSampling','0.0')
  DSMC%NumOutput = GETINT('Particles-NumberForDSMCOutputs','0')
  IF((DSMC%TimeFracSamp.GT.0.0).AND.(DSMC%NumOutput.EQ.0)) DSMC%NumOutput = 1
  DSMC%CalcQualityFactors = GETLOGICAL('Particles-DSMC-CalcQualityFactors','.FALSE.')
  DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim','.FALSE.')
#if (PP_TimeDiscMethod==42)
  DSMC%CalcQualityFactors = .TRUE.
#endif
  DSMC%ReservoirSimuRate = GETLOGICAL('Particles-DSMCReservoirSimRate','.FALSE.')
  DSMC%ReservoirRateStatistic = GETLOGICAL('Particles-DSMCReservoirStatistic','.FALSE.')
  DSMC%VibEnergyModel = GETINT('Particles-ModelForVibrationEnergy','0')
  DSMC%DoTEVRRelaxation = GETLOGICAL('Particles-DSMC-TEVR-Relaxation','.FALSE.')
  DSMC%WallModel = GETINT('Particles-DSMC-WallModel','0') !0: elastic/diffusive reflection, 1:ad-/desorption empiric, 2:chem. ad-/desorption UBI-QEP
  LD_MultiTemperaturMod=GETINT('LD-ModelForMultiTemp','0')
  DSMC%ElectronicModel = GETLOGICAL('Particles-DSMC-ElectronicModel','.FALSE.')
  DSMC%ElectronicModelDatabase = TRIM(GETSTR('Particles-DSMCElectronicDatabase','none'))
  IF ((DSMC%ElectronicModelDatabase .NE. 'none').AND.(CollisMode .GT. 1)) THEN 
    DSMC%EpsElecBin = GETREAL('EpsMergeElectronicState','1E-4')
  ELSEIF(DSMC%ElectronicModel) THEN
    CALL Abort(&
        __STAMP__,&
        'ERROR: Electronic model requires a electronic levels database and CollisMode > 1!')
  END IF
  IF ((DSMC%VibEnergyModel.EQ.1).AND.(CollisMode.EQ.3)) THEN
    CALL Abort(&
     __STAMP__&
    ,'TSHO Model is not working with Chemical Reactions !!!')
  ELSE IF ((DSMC%VibEnergyModel.GT.1).OR.(DSMC%VibEnergyModel.LT.0)) THEN
    CALL Abort(&
     __STAMP__&
    ,'ERROR in ModelForVibrationEnergy Flag!')
  END IF
  IF (DSMC%NumOutput.NE.0) THEN
    DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / REAL(DSMC%NumOutput)
  END IF
  DSMC%NumPolyatomMolecs = 0
  ! Steady - State Detection: Use Q-Criterion or SSD-Alogrithm?
  SamplingActive = .FALSE.
  UseQCrit = GETLOGICAL('Particles-DSMC-UseQCrit','.FALSE.')
  UseSSD = GETLOGICAL('Particles-DSMC-UseSSD','.FALSE.')
  IF(UseQCrit.OR.UseSSD) CALL DSMC_SteadyStateInit()

  ALLOCATE(HValue(nElems))
  HValue(1:nElems) = 0.0

  IF(DSMC%CalcQualityFactors) THEN
    ALLOCATE(DSMC%QualityFacSamp(nElems,3))
    DSMC%QualityFacSamp(1:nElems,1:3) = 0.0
    ALLOCATE(DSMC%QualityFactors(nElems,3))
    DSMC%QualityFactors(1:nElems,1:3) = 0.0
  END IF

! definition of DSMC particle values
  ALLOCATE(DSMC_RHS(PDM%maxParticleNumber,3))
  DSMC_RHS = 0

  IF (nSpecies.LE.0) THEN 
    CALL Abort(&
    __STAMP__&
    ,"ERROR: nSpecies .LE. 0:", nSpecies)
  END IF

  IF (CollisMode.EQ.0) THEN
#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001) || (PP_TimeDiscMethod==42)
    CALL Abort(&
     __STAMP__&
    , "Free Molecular Flow (CollisMode=0) is not supported for LD or DEBUG!")
#endif
  ELSE !CollisMode.GT.0

! reading species data of ini_2
  ALLOCATE(SpecDSMC(nSpecies))
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    SpecDSMC(iSpec)%Name    = TRIM(GETSTR('Part-Species'//TRIM(hilf)//'-SpeciesName','none'))
    SpecDSMC(iSpec)%InterID = GETINT('Part-Species'//TRIM(hilf)//'-InteractionID','0')
    SpecDSMC(iSpec)%TrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceTemp','0')
    SpecDSMC(iSpec)%DrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceDiam','0')
    IF((SpecDSMC(iSpec)%InterID*SpecDSMC(iSpec)%TrefVHS*SpecDSMC(iSpec)%DrefVHS).eq.0) THEN
    CALL Abort(&
    __STAMP__&
    ,"ERROR in species data ini_2")
    END IF
    SpecDSMC(iSpec)%omegaVHS = GETREAL('Part-Species'//TRIM(hilf)//'-omegaVHS','0') ! default case HS 
  ! reading electronic state informations from HDF5 file
    IF((DSMC%ElectronicModelDatabase .NE. 'none').AND.(SpecDSMC(iSpec)%InterID .NE. 4)) THEN
      IF(SpecDSMC(iSpec)%Name.EQ.'none') THEN
        CALL Abort(&
           __STAMP__,&
           "Read-in from electronic database requires the definition of species name! Species:",iSpec)
      END IF
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
! Here, something strange is happen!
!      A1 = 0.5 * SQRT(Pi) * SpecDSMC(iSpec)%DrefVHS*(2*(2-SpecDSMC(iSpec)%omegaVHS) &
!            * BoltzmannConst * SpecDSMC(iSpec)%TrefVHS)**(SpecDSMC(iSpec)%omegaVHS*0.5)
!      A2 = 0.5 * SQRT(Pi) * SpecDSMC(jSpec)%DrefVHS*(2*(2-SpecDSMC(jSpec)%omegaVHS) &
!            * BoltzmannConst * SpecDSMC(jSpec)%TrefVHS)**(SpecDSMC(jSpec)%omegaVHS*0.5)
      A1 = 0.5 * SQRT(Pi) * SpecDSMC(iSpec)%DrefVHS*(2*BoltzmannConst*SpecDSMC(iSpec)%TrefVHS)**(SpecDSMC(iSpec)%omegaVHS*0.5) &
            /SQRT(GAMMA(2.0 - SpecDSMC(iSpec)%omegaVHS))
      A2 = 0.5 * SQRT(Pi) * SpecDSMC(jSpec)%DrefVHS*(2*BoltzmannConst*SpecDSMC(jSpec)%TrefVHS)**(SpecDSMC(jSpec)%omegaVHS*0.5) &
            /SQRT(GAMMA(2.0 - SpecDSMC(jSpec)%omegaVHS))
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
  ! allocate internal energy arrays
    IF ( DSMC%ElectronicModel ) THEN
      ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,3))
    ELSE
      ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,2))
    ENDIF
    ! reading molecular stuff
    SpecDSMC(1:nSpecies)%Xi_Rot = 0
    SpecDSMC(1:nSpecies)%MaxVibQuant = 0
    SpecDSMC(1:nSpecies)%CharaTVib = 0
    SpecDSMC(1:nSpecies)%EZeroPoint = 0.0
    SpecDSMC(1:nSpecies)%PolyatomicMol=.false.
    SpecDSMC(1:nSpecies)%SpecToPolyArray = 0
    ! Check whether calculation of instantaneous translational temperature is require
    IF(DSMC%BackwardReacRate.OR.(SelectionProc.EQ.2)) THEN
      ALLOCATE(DSMC%InstantTransTemp(nSpecies+1))
      DSMC%InstantTransTemp = 0.0
    END IF
    DO iSpec = 1, nSpecies
      IF(.NOT.(SpecDSMC(iSpec)%InterID.EQ.4)) THEN
        WRITE(UNIT=hilf,FMT='(I2)') iSpec
        SpecDSMC(iSpec)%PolyatomicMol=GETLOGICAL('Part-Species'//TRIM(hilf)//'-PolyatomicMol','.FALSE.')
        IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%ElectronicModel)  THEN
          CALL Abort(&
          __STAMP__&
          ,'! Simulation of Polyatomic Molecules and Electronic States are not possible yet!!!')
        END IF
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
          SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
        ELSEIF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          SpecDSMC(iSpec)%Xi_Rot     = 2
          SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib','0.')
          SpecDSMC(iSpec)%CharaTRot  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempRot','0')
          SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')
          IF(SpecDSMC(iSpec)%Ediss_eV*SpecDSMC(iSpec)%CharaTVib.EQ.0) THEN
            CALL Abort(&
            __STAMP__&
            ,'Error! Ediss_eV or CharaTVib is not set or equal to zero!')
          ELSE
            IF (DSMC%VibEnergyModel.EQ.0) THEN
              SpecDSMC(iSpec)%MaxVibQuant = 200
            ELSE
              SpecDSMC(iSpec)%MaxVibQuant = INT(SpecDSMC(iSpec)%Ediss_eV*JToEv/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)) + 1
            END IF
            ! Calculation of the zero-point energy
            SpecDSMC(iSpec)%EZeroPoint = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
          END IF
        END IF
        SpecDSMC(iSpec)%VFD_Phi3_Factor = GETREAL('Part-Species'//TRIM(hilf)//'-VFDPhi3','0.')
        ! Read in species values for rotational relaxation models of Boyd/Zhang if necessary
        IF(DSMC%RotRelaxProb.GT.1.0) THEN
          SpecDSMC(iSpec)%CollNumRotInf = GETREAL('Part-Species'//TRIM(hilf)//'-CollNumRotInf','0.')
          SpecDSMC(iSpec)%TempRefRot    = GETREAL('Part-Species'//TRIM(hilf)//'-TempRefRot','0.')
          IF(SpecDSMC(iSpec)%CollNumRotInf*SpecDSMC(iSpec)%TempRefRot.EQ.0) THEN
            CALL Abort(&
            __STAMP__&
            ,'Error! CollNumRotRef or TempRefRot is not set or equal to zero!')
          END IF
        END IF 
        ! Read in species values for vibrational relaxation models of Milikan-White if necessary
        IF(DSMC%VibRelaxProb.GT.1.0) THEN
          ALLOCATE(SpecDSMC(iSpec)%MW_Const(1:nSpecies))
          DO jSpec = 1, nSpecies
            WRITE(UNIT=hilf2,FMT='(I2)') jSpec
            hilf2=TRIM(hilf)//'-'//TRIM(hilf2)
            SpecDSMC(iSpec)%MW_Const(jSpec)     = GETREAL('Part-Species'//TRIM(hilf)//'-MWConst-'//TRIM(hilf2),'0.')
            IF(SpecDSMC(iSpec)%MW_Const(jSpec).EQ.0) THEN
              CALL Abort(&
              __STAMP__&
              ,'Error! MWConst is not set or equal to zero! Spec-Pair', iSpec)
            END IF
          END DO
          SpecDSMC(iSpec)%CollNumVib     = GETREAL('Part-Species'//TRIM(hilf)//'-CollNumVib','0.')
          SpecDSMC(iSpec)%VibCrossSec    = GETREAL('Part-Species'//TRIM(hilf)//'-VibCrossSection','10E-20')
          IF(SpecDSMC(iSpec)%CollNumVib.EQ.0) THEN
            CALL Abort(&
            __STAMP__&
            ,'Error! CollNumVib not set or equal to zero for Species!', iSpec)
          END IF
        END IF 
        ! Setting the values of Rot-/Vib-RelaxProb to a fix value
        SpecDSMC(iSpec)%RotRelaxProb  = DSMC%RotRelaxProb
        SpecDSMC(iSpec)%VibRelaxProb  = DSMC%VibRelaxProb    !0.02
        SpecDSMC(iSpec)%ElecRelaxProb = DSMC%ElecRelaxProb    !or 0.02 | Bird: somewhere in range 0.01 .. 0.02
        ! multi init stuff
        ALLOCATE(SpecDSMC(iSpec)%Init(0:Species(iSpec)%NumberOfInits))
        DO iInit = 0, Species(iSpec)%NumberOfInits
          IF (iInit .EQ. 0) THEN !0. entry := old style parameter def. (default values if not def., some values might be needed)
            hilf2=TRIM(hilf)
          ELSE ! iInit >0
            WRITE(UNIT=hilf2,FMT='(I2)') iInit
            hilf2=TRIM(hilf)//'-Init'//TRIM(hilf2)
          END IF ! iInit
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            SpecDSMC(iSpec)%Init(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','0.')
            SpecDSMC(iSpec)%Init(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','0.')
            IF (SpecDSMC(iSpec)%Init(iInit)%TRot*SpecDSMC(iSpec)%Init(iInit)%TVib.EQ.0.) THEN
              IF (iInit.EQ.0)THEN
                IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                  CALL Abort(&
                  __STAMP__&
                  ,'Error! TVib and TRot need to be defined in Part-SpeciesXX-TempVib/TempRot for iSpec',iSpec)
                ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                  CALL Abort(&
                  __STAMP__&
                  ,'Error! TVib and TRot need to be defined in Part-SpeciesXX-TempVib/TempRot for BGGas')
                END IF
              ELSE ! iInit >0
                CALL Abort(&
                __STAMP__&
                ,'Error! TVib and TRot need to be defined in Part-SpeciesXX-InitXX-TempVib/TempRot for iSpec, iInit'&
                ,iSpec,REAL(iInit))
              END IF
            END IF
          END IF
          ! read electronic temperature
          IF ( DSMC%ElectronicModel ) THEN
            WRITE(UNIT=hilf,FMT='(I2)') iSpec
            SpecDSMC(iSpec)%Init(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-Tempelec','0.')
            IF (SpecDSMC(iSpec)%Init(iInit)%Telec.EQ.0.) THEN
              IF (iInit.EQ.0)THEN
                IF (Species(iSpec)%StartnumberOfInits.EQ.0)THEN
                  CALL Abort(&
                  __STAMP__&
                  ,' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for Species',iSpec)
                ELSE IF (BGGas%BGGasSpecies.EQ.iSpec) THEN !cases which need values of fixed iInit=0 (indep. from Startnr.OfInits)
                  CALL Abort(&
                  __STAMP__&
                  ,' Error! Telec needs to defined in Part-SpeciesXX-Tempelec for BGGas')
                END IF
              ELSE ! iInit >0
                CALL Abort(&
                __STAMP__&
                ,' Error! Telec needs to defined in Part-SpeciesXX-InitXX-Tempelc for iSpec, iInit',iSpec,REAL(iInit))
              END IF
            END IF
          END IF
        END DO !Inits
        ALLOCATE(SpecDSMC(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
        DO iInit = 1, Species(iSpec)%nSurfacefluxBCs
          IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            WRITE(UNIT=hilf2,FMT='(I2)') iInit
            hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempVib','0.')
            SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot      = GETREAL('Part-Species'//TRIM(hilf2)//'-TempRot','0.')
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%TRot*SpecDSMC(iSpec)%Surfaceflux(iInit)%TVib.EQ.0.) THEN
              CALL Abort(&
              __STAMP__&
              ,'Error! TVib and TRot not def. in Part-SpeciesXX-SurfacefluxXX-TempVib/TempRot for iSpec, iSF',iSpec,REAL(iInit))
            END IF
          END IF
          ! read electronic temperature
          IF ( DSMC%ElectronicModel ) THEN
            WRITE(UNIT=hilf,FMT='(I2)') iSpec
            SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec   = GETREAL('Part-Species'//TRIM(hilf2)//'-Tempelec','0.')
            IF (SpecDSMC(iSpec)%Surfaceflux(iInit)%Telec.EQ.0.) THEN
              CALL Abort(&
              __STAMP__&
              ,' Error! Telec not defined in Part-SpeciesXX-SurfacefluxXX-Tempelc for iSpec, iSF',iSpec,REAL(iInit))
            END IF
          END IF
        END DO !SurfaceFluxBCs
      END IF ! not electron
    END DO !Species

    ! Initialization of polyatomic species and burn-in phase (Metropolis-Hastings) per initialization region
    IF(DSMC%NumPolyatomMolecs.GT.0) THEN
      DSMC%PolySingleMode = GETLOGICAL('Particles-DSMC-PolyRelaxSingleMode','.FALSE.')
      IF((SelectionProc.NE.2).AND.DSMC%PolySingleMode) THEN
        ! Single-mode relaxation of vibrational modes of polyatomic molecules only possible when the prohibiting double relaxation
        ! method is used (Gimelshein, SelectionProc = 2)
        CALL Abort(&
         __STAMP__&
        ,'ERROR: No single-mode polyatomic relaxation possible with chosen selection procedure! SelectionProc:', SelectionProc)
      END IF
      ALLOCATE(VibQuantsPar(PDM%maxParticleNumber))
      ALLOCATE(PolyatomMolDSMC(DSMC%NumPolyatomMolecs))
      DO iSpec = 1, nSpecies
        IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
          CALL InitPolyAtomicMolecs(iSpec)
!          ! Required if the Metropolis-Hastings random-walk is utilized for the initialization of polyatomic molecules (sampling
!          ! of all modes at once)
!          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!          ALLOCATE( &
!              PolyatomMolDSMC(iPolyatMole)%LastVibQuantNums(1:PolyatomMolDSMC(iPolyatMole)%VibDOF, &
!                                                             0:Species(iSpec)%NumberOfInits+Species(iSpec)%nSurfacefluxBCs))
!          DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits 
!            CALL DSMC_FindFirstVibPick(iInit, iSpec, 1)
!          END DO
!          DO iInit = 1,Species(iSpec)%nSurfacefluxBCs
!            CALL DSMC_FindFirstVibPick(iInit, iSpec, 2)
!          END DO
        END IF
      END DO
    END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Comparison with Landau-Teller equation, including a different selection procedure (restricts relaxation to a single mode)
! and a correctional factor, both in dsmc_collis_mode.f90, also the translational temperature is fixed, in timedisc.f90
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==42)
    DSMC%CompareLandauTeller = GETLOGICAL('Particles-DSMC-CompareLandauTeller','.FALSE.')
    IF(DSMC%CompareLandauTeller) THEN
      IF(CollisMode.NE.2) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Comparison with Landau-Teller only available in CollisMode = 2, CollisMode:', CollisMode)
      END IF
      IF(nSpecies.GT.1) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Comparison with Landau-Teller only available for a single species, nSpecies:', nSpecies)
      END IF
    END IF
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==42)
    IF ( DSMC%ElectronicModel ) THEN
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
    ! Setting the internal energy value of every particle
    DO iPart = 1, PDM%ParticleVecLength
      IF (PDM%ParticleInside(iPart)) THEN
        IF (Species(PartSpecies(iPart))%NumberOfInits.EQ.0) THEN
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(PartSpecies(iPart),0,iPart,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),0,iPart,1)
          END IF
        ELSE
          iInit = PDM%PartInit(iPart)
          IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
            CALL DSMC_SetInternalEnr_Poly(PartSpecies(iPart),iInit,iPart,1)
          ELSE
            CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iInit,iPart,1)
          END IF
        END IF
      END IF
    END DO
    
#if ( PP_TimeDiscMethod ==42 )
    ! Debug Output for initialized electronic state
    IF ( DSMC%ElectronicModel ) THEN
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

#if (PP_TimeDiscMethod!=1000) && (PP_TimeDiscMethod!=1001) && (PP_TimeDiscMethod!=300)
    DEALLOCATE(PDM%PartInit)
#endif
  END IF ! CollisMode .EQ. 2 or 3
!-----------------------------------------------------------------------------------------------------------------------------------
! Define chemical reactions (including ionization and backward reaction rate)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (CollisMode.EQ.3) THEN ! perform chemical reactions
    IF(DSMC%BackwardReacRate) THEN
      IF(MOD(DSMC%PartitionMaxTemp,DSMC%PartitionInterval).EQ.0.0) THEN
        PartitionArraySize = INT(DSMC%PartitionMaxTemp / DSMC%PartitionInterval)
      ELSE
        CALL abort(&
        __STAMP__&
        ,'ERROR: Partition temperature limit must be multiple of partition interval!')
      END IF
    END IF
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I2)') iSpec
      SpecDSMC(iSpec)%HeatOfFormation = GETREAL('Part-Species'//TRIM(hilf)//'-HeatOfFormation_K','0.0')*BoltzmannConst
      ! Read-in of species parameters for the partition function calculation -------------------------------------------------------
      IF(DSMC%BackwardReacRate) THEN
        IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          SpecDSMC(iSpec)%SymmetryFactor              = GETINT('Part-Species'//TRIM(hilf)//'-SymmetryFactor','0')
          IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
            IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
              IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
                CALL abort(&
                __STAMP__&
                ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
              END IF
            ELSE
              IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)  &
                 * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
                CALL abort(&
                __STAMP__&
                ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
              END IF
            END IF
          ELSE            
            IF(SpecDSMC(iSpec)%CharaTRot*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
              __STAMP__&
              ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for backward rate!', iSpec)
            END IF
          END IF
        END IF
        IF (DSMC%ElectronicModelDatabase .EQ. 'none') THEN
          SpecDSMC(iSpec)%MaxElecQuant               = GETINT('Part-Species'//TRIM(hilf)//'-NumElectronicLevels','0')
          IF(SpecDSMC(iSpec)%MaxElecQuant.GT.0) THEN
            ALLOCATE(SpecDSMC(iSpec)%ElectronicState(2,0:SpecDSMC(iSpec)%MaxElecQuant-1))
            DO iDOF=1, SpecDSMC(iSpec)%MaxElecQuant
              WRITE(UNIT=hilf2,FMT='(I2)') iDOF
              SpecDSMC(iSpec)%ElectronicState(1,iDOF-1) &
                = GETINT('Part-Species'//TRIM(hilf)//'-ElectronicDegeneracy-Level'//TRIM(hilf2),'0')
              SpecDSMC(iSpec)%ElectronicState(2,iDOF-1) &
                = GETREAL('Part-Species'//TRIM(hilf)//'-ElectronicEnergyLevel-Level'//TRIM(hilf2),'0')
            END DO
          END IF
        END IF
        ALLOCATE(SpecDSMC(iSpec)%PartitionFunction(1:PartitionArraySize))
        DO iInter = 1, PartitionArraySize
          Temp = iInter * DSMC%PartitionInterval
          CALL CalcPartitionFunction(iSpec, Temp, Qtra, Qrot, Qvib, Qelec)
          SpecDSMC(iSpec)%PartitionFunction(iInter) = Qtra * Qrot * Qvib * Qelec
        END DO
      END IF
      !-----------------------------------------------------------------------------------------------------------------------------
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
  ELSE IF (DSMC%WallModel.GT.0 .AND. CollisMode.GT.1) THEN
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilf,FMT='(I2)') iSpec
      IF((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
        SpecDSMC(iSpec)%SymmetryFactor              = GETINT('Part-Species'//TRIM(hilf)//'-SymmetryFactor','0')
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
            IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
              __STAMP__&
              ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
            END IF
          ELSE
            IF(PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)*PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)  &
                * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
              CALL abort(&
              __STAMP__&
              ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
            END IF
          END IF
        ELSE            
          IF(SpecDSMC(iSpec)%CharaTRot*SpecDSMC(iSpec)%SymmetryFactor.EQ.0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR: Char. rotational temperature or symmetry factor not defined properly for Adsorptionmodel!', iSpec)
          END IF
        END IF
      END IF
    END DO
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! Octree + Nearest Neighbour: octree cell splitting is performed until the mean free path is resolved or a certain number of
! particles is reached. Afterwards a nearest neighbour search for the collision partner is performed.
! Source: Pfeiffer, M., Mirza, A. and Fasoulas, S. (2013). A grid-independent particle pairing strategy for DSMC.
! Journal of Computational Physics 246, 28–36. doi:10.1016/j.jcp.2013.03.018
!-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%UseOctree = GETLOGICAL('Particles-DSMC-UseOctree','.FALSE.')
  ! If number of particles is greater than OctreePartNumNode, cell is going to be divided for performance of nearest neighbour
  DSMC%PartNumOctreeNode = GETINT('Particles-OctreePartNumNode','80')
  ! If number of particles is less than OctreePartNumNodeMin, cell is NOT going to be split even if mean free path is not resolved
  ! 50 / 8 -> ca. 6-7 particles per cell
  DSMC%PartNumOctreeNodeMin = GETINT('Particles-OctreePartNumNodeMin','50')
  IF (DSMC%PartNumOctreeNodeMin.LT.20) THEN
    CALL abort(&
    __STAMP__&
    ,'Particles-OctreePartNumNodeMin is less than 20')
  END IF
  IF(DSMC%UseOctree) THEN
    IF(NGeo.GT.PP_N) CALL abort(&
    __STAMP__&
    ,' Set PP_N to NGeo, else, the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set mean VibQua of BGGas for dissoc reaction
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (BGGas%BGGasSpecies.NE.0) THEN
    IF (DSMC%UseOctree) THEN
      CALL abort(__STAMP__,&
              'ERROR: Utilization of the octree and nearest neighbour scheme not possible with the background gas')
    END IF
    IF (BGGas%BGGasDensity.EQ.0) THEN
      CALL abort(__STAMP__,&
        'ERROR: Background gas density is not defined. Set Particles-DSMCBackgroundGasDensity (real number density) !')
    END IF
    IF(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.4) THEN
      CALL abort(__STAMP__,&
        'ERROR: Electrons as background gas are not yet available!!')
    END IF
    IF((SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.2).OR.(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.20)) THEN
      IF(SpecDSMC(BGGas%BGGasSpecies)%PolyatomicMol) THEN
        CALL abort(&
        __STAMP__&
        ,'ERROR: Polyatomic species as background gas are not yet available!')
      ELSE
        BGGasEVib = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib & 
                  + BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib  &
                  /  (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1) &
                  - BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                  / (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                  / SpecDSMC(BGGas%BGGasSpecies)%Init(0)%TVib) - 1)
        BGGasEVib = BGGasEVib/(BoltzmannConst*SpecDSMC(BGGas%BGGasSpecies)%CharaTVib) - DSMC%GammaQuant
        BGGas%BGMeanEVibQua = MIN(INT(BGGasEVib) + 1, SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant)
      END IF
    ELSE
      BGGas%BGMeanEVibQua = 0
    END IF
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate vib collision numbers, according to Boyd & Abe
!-----------------------------------------------------------------------------------------------------------------------------------
  IF(DSMC%VibRelaxProb.EQ.2) THEN
    DO iSpec = 1, nSpecies
      ALLOCATE(SpecDSMC(iSpec)%CharaVelo(1:nSpecies))
      DO jSpec = 1, nSpecies
        iCase = CollInf%Coll_Case(iSpec,jSpec)
        SpecDSMC(iSpec)%CharaVelo(jSpec) = SQRT(  BoltzmannConst / CollInf%MassRed(iCase) &
                                                * (2./3.*SpecDSMC(iSpec)%MW_Const(jSpec))**3.)
        SpecDSMC(iSpec)%CollNumVib = SpecDSMC(iSpec)%CollNumVib * (2.*(2.-SpecDSMC(iSpec)%omegaVHS)*BoltzmannConst &
                                   * SpecDSMC(iSpec)%TrefVHS / CollInf%MassRed(iCase))**SpecDSMC(iSpec)%omegaVHS

      END DO
    END DO
  END IF

  END IF !CollisMode.GT.0
!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing Surface stuff
!-----------------------------------------------------------------------------------------------------------------------------------
  DSMC%CalcSurfaceVal = GETLOGICAL('Particles-DSMC-CalcSurfaceVal','.FALSE.')
  IF (DSMC%CalcSurfaceVal) THEN
    DSMC%CalcSurfaceTime = GETLOGICAL('Particles-DSMC-CalcSurfaceTime','.FALSE.')
    CALL InitParticleBoundarySampling()
    
    DSMC%CalcSurfCollis_OnlySwaps = GETLOGICAL('Particles-DSMC-CalcSurfCollis_OnlySwaps','.FALSE.')
    DSMC%CalcSurfCollis_Only0Swaps = GETLOGICAL('Particles-DSMC-CalcSurfCollis_Only0Swaps','.FALSE.')
    DSMC%CalcSurfCollis_Output = GETLOGICAL('Particles-DSMC-CalcSurfCollis_Output','.FALSE.')
    IF (DSMC%CalcSurfCollis_Only0Swaps) DSMC%CalcSurfCollis_OnlySwaps=.TRUE.
    DSMC%AnalyzeSurfCollis = GETLOGICAL('Particles-DSMC-AnalyzeSurfCollis','.FALSE.')
    IF (DSMC%AnalyzeSurfCollis) THEN
      AnalyzeSurfCollis%maxPartNumber = GETINT('Particles-DSMC-maxSurfCollisNumber','0')
      ALLOCATE(AnalyzeSurfCollis%Data(1:AnalyzeSurfCollis%maxPartNumber,1:9))
      ALLOCATE(AnalyzeSurfCollis%Spec(1:AnalyzeSurfCollis%maxPartNumber))
      ALLOCATE(AnalyzeSurfCollis%Number(1:nSpecies+1))
      !ALLOCATE(AnalyzeSurfCollis%Rate(1:nSpecies+1))
      AnalyzeSurfCollis%Data=0.
      AnalyzeSurfCollis%Spec=0
      AnalyzeSurfCollis%Number=0
      !AnalyzeSurfCollis%Rate=0.
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
      CALL abort(&
      __STAMP__&
      ,'Error in Particles-DSMC-CalcSurfCollis_NbrOfSpecies!')
    END IF
  END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! Initialize surface model (Adsorption/Desorption/Reactions) variables
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (DSMC%WallModel.GT.0 .AND. CollisMode.GT.1) THEN
    IF (.NOT.DSMC%CalcSurfaceVal) THEN
      CALL InitParticleBoundarySampling()
      SWRITE(UNIT_stdOut,'(A)')'WARNING: Particles-DSMC-CalcSurfaceVal == FALSE!'
    END IF
    IF ((DSMC%WallModel.EQ.2) .OR.(DSMC%WallModel.EQ.1)) CALL abort(&
      __STAMP__&
      ,'Error: WallModel 1&2 not working!')
    CALL InitDSMCSurfModel()
  ELSE IF (DSMC%WallModel.GT.0 .AND. CollisMode.LE.1) THEN
    CALL abort(&
        __STAMP__&
        ,'Error in DSMC-Surface model init - wrong collismode!')
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
  SWRITE(UNIT_stdOut,'(A)')' INIT DSMC DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitDSMC


SUBROUTINE DSMC_SetInternalEnr_LauxVFD(iSpecies, iInit, iPart, init_or_sf)
!===================================================================================================================================
! Energy distribution according to dissertation of Laux (diatomic)
!===================================================================================================================================
! MODULES
  USE MOD_Globals,              ONLY : abort
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
  INTEGER, INTENT(IN)           :: iSpecies, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                          :: iRan
  INTEGER                       :: iQuant
  REAL                        :: TVib                       ! vibrational temperature
  REAL                        :: TRot                       ! rotational temperature
  INTEGER                     :: iInitTemp, init_or_sfTemp
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! Set internal energies (vibrational and rotational)
!-----------------------------------------------------------------------------------------------------------------------------------
  IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
    SELECT CASE (init_or_sf)
    CASE(1) !iInit
      TVib=SpecDSMC(iSpecies)%Init(iInit)%TVib
      TRot=SpecDSMC(iSpecies)%Init(iInit)%TRot
    CASE(2) !SurfaceFlux
      TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
      TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,'neither iInit nor Surfaceflux defined as reference!')
    END SELECT
    ! Set vibrational energy
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
    DO WHILE (iQuant.GE.SpecDSMC(iSpecies)%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
    END DO
    !evtl muss partstateinten nochmal geändert werden, mpi, resize etc..
    PartStateIntEn(iPart, 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpecies)%CharaTVib*BoltzmannConst
    ! Set rotational energy
    CALL RANDOM_NUMBER(iRan)
    PartStateIntEn(iPart, 2) = -BoltzmannConst*TRot*LOG(iRan)
  ELSE
    ! Nullify energy for atomic species
    PartStateIntEn(iPart, 1) = 0
    PartStateIntEn(iPart, 2) = 0
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set electronic energy
!-----------------------------------------------------------------------------------------------------------------------------------
  IF ( DSMC%ElectronicModel .and. SpecDSMC(iSpecies)%InterID .ne. 4) THEN
    iInitTemp = iInit
    init_or_sfTemp = init_or_sf
    CALL InitElectronShell(iSpecies,iPart,iInitTemp,init_or_sfTemp)
  ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set internal energy for LD
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1000) || (PP_TimeDiscMethod==1001)
  IF (LD_MultiTemperaturMod .EQ. 3 ) THEN ! no discret vib levels for this LD method
    IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
      PartStateIntEn(iPart, 1) = (BoltzmannConst*SpecDSMC(iSpecies)%CharaTVib) &
                               / (EXP(SpecDSMC(iSpecies)%CharaTVib/TVib) - 1.0) &
                               + DSMC%GammaQuant * BoltzmannConst*SpecDSMC(iSpecies)%CharaTVib
    !set rotational energy
      PartStateIntEn(iPart, 2) = BoltzmannConst*TRot
    ELSE
      PartStateIntEn(iPart, 1) = 0
      PartStateIntEn(iPart, 2) = 0
    END IF
  END IF
#endif

END SUBROUTINE DSMC_SetInternalEnr_LauxVFD


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
!    CALL abort(__STAMP__&
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
!    CALL abort(__STAMP__&
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
!    CALL abort(__STAMP__&
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
!    CALL abort(__STAMP__&
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
!    CALL abort(__STAMP__&
!    'ERROR: wrong number of points in VTK-File')
!  END IF
!  IF (ncells.NE.nElems) THEN
!    CALL abort(__STAMP__&
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
!            CALL abort(__STAMP__&
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
!  !  CALL abort(__STAMP__&
!  !	'ERROR: Not all VTKcells mapped for BGG!',999,999.)
!  !END IF
!  IF (.NOT. ALL(IsAssociated2)) THEN
!    CALL abort(__STAMP__&
!              'ERROR: Not all Elems mapped for BGG!')
!  END IF
!  SWRITE(UNIT_stdOut,'(A)')'DONE!'
!  SWRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE ReadinVTKFileBGG


SUBROUTINE FinalizeDSMC() 
!----------------------------------------------------------------------------------------------------------------------------------!
! finalize dsmc variables
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_DSMC_Vars
USE MOD_Particle_Vars, ONLY:PDM
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
SDEALLOCATE(DSMC%InstantTransTemp)
IF(DSMC%CalcQualityFactors) THEN
  SDEALLOCATE(DSMC%QualityFacSamp)
  SDEALLOCATE(DSMC%QualityFactors)
END IF
SDEALLOCATE(PDM%PartInit)
SDEALLOCATE(Coll_pData)
SDEALLOCATE(SampDSMC)
SDEALLOCATE(MacroDSMC)
SDEALLOCATE(QKBackWard)
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
SDEALLOCATE(ChemReac%ReactNumRecomb)
SDEALLOCATE(ChemReac%Hab)
SDEALLOCATE(CollInf%Coll_Case)
SDEALLOCATE(CollInf%Coll_CaseNum)
SDEALLOCATE(CollInf%Coll_SpecPartNum)
SDEALLOCATE(CollInf%Cab)
SDEALLOCATE(CollInf%KronDelta)
SDEALLOCATE(CollInf%FracMassCent)
SDEALLOCATE(CollInf%MassRed)
SDEALLOCATE(HValue)
!SDEALLOCATE(SampWall)
SDEALLOCATE(MacroSurfaceVal)
SDEALLOCATE(VibQuantsPar)
! SDEALLOCATE(XiEq_Surf)
SDEALLOCATE(DSMC_HOSolution)
SDEALLOCATE(ElemNodeVol)
END SUBROUTINE FinalizeDSMC


END MODULE MOD_DSMC_Init
