#include "boltzplatz.h"

MODULE MOD_DSMC_Init
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

PUBLIC :: InitDSMC, DSMC_SetInternalEnr_LauxVFD
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDSMC()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals!, ONLY: MPIRoot
USE MOD_Mesh_Vars,             ONLY : nElems
USE MOD_ReadInTools
USE MOD_DSMC_ElectronicModel,  ONLY: ReadSpeciesLevel
USE MOD_DSMC_Vars!, ONLY: 
USE MOD_PARTICLE_Vars,         ONLY: nSpecies, BoltzmannConst, Species, PDM, PartSpecies
USE MOD_Equation_Vars,         ONLY: Pi
USE MOD_DSMC_Analyze,          ONLY: WriteOutputMesh
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
  INTEGER               :: iCase, iSpec, jSpec, nCase, iVec, iPart , iINit
  REAL                  :: A1, A2     ! species constant for cross section (p. 24 Laux)
  CHARACTER(LEN=12)     :: rmCommand 
  REAL                  :: JToEv
  REAL                  :: SpecVib(1:nSpecies), SpecRot(1:nSpecies)
#if ( PP_TimeDiscMethod ==42 )
  CHARACTER(LEN=64)     :: DebugElectronicStateFilename
  INTEGER               :: ii
#endif
!#ifdef MPI
!#endif
!===================================================================================================================================
  JToEv = 1.602176565E-19  

  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' DSMC INIT ...'
  
  
! writing od output mesh
  SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH...'
  CALL WriteOutputMesh()
!  write(rmCommand,'(A12)') 'rm DSMCOut_*'
!  CALL SYSTEM(rmCommand)
  SWRITE(UNIT_stdOut,'(A)')' WRITING OUTPUT-MESH DONE!'
! reading and reset general DSMC values
  CollisMode = GETINT('Particles-DSMC-CollisMode','1') !1:elastic col, 2:elast+rela, 3:chem
  DSMC%GammaQuant   = GETREAL('Particles-DSMC-GammaQuant', '0.5')
  DSMC%TimeFracSamp = GETREAL('Part-TimeFracForSampling','0')
  DSMC%NumOutput = GETINT('Particles-NumberForDSMCOutputs','1')
  DSMC%ReservoirSimu = GETLOGICAL('Particles-DSMCReservoirSim','.FALSE.')
  DSMC%ReservoirSimuRate = GETLOGICAL('Particles-DSMCReservoirSimRate','.FALSE.')
  DSMC%ReservoirRateStatistic = GETLOGICAL('Particles-DSMCReservoirStatistic','.FALSE.')
  DSMC%VibEnergyModel = GETINT('Particles-ModelForVibrationEnergy','0')
  DSMC%CollProbMaxOut = GETLOGICAL('Particles-DSMCOutputOfCollisionProb','.FALSE.')
  DSMC%ElectronicStateDatabase = GETSTR('Particles-DSMCElectronicDatabase','none')
  DSMC%ElectronicState = .FALSE.
  IF ( (DSMC%ElectronicStateDatabase .ne. 'none') .AND. (CollisMode .GT. 1)  ) THEN 
   DSMC%ElectronicState = .TRUE.
  END IF
  DSMC%CollMean = 0.0
  DSMC%CollMeanCount = 0.0
  IF ((DSMC%VibEnergyModel.EQ.1).AND.(CollisMode.EQ.3)) THEN
    PRINT*, 'TSHO Model is not working with Chemical Reactions !!!'
    STOP
  ELSE IF ((DSMC%VibEnergyModel.GT.1).OR.(DSMC%VibEnergyModel.LT.0)) THEN
    PRINT*, 'ERROR in ModelForVibrationEnergy Flag!'
    STOP
  END IF
  IF (DSMC%NumOutput.NE.0) THEN
    DSMC%DeltaTimeOutput = (DSMC%TimeFracSamp * TEnd) / DSMC%NumOutput
  END IF
  DSMC%NumPolyatomMolecs = 0
! definition of DSMC sampling values
  DSMC%SampNum = 0
  ALLOCATE(SampDSMC(nElems,nSpecies))
  SampDSMC(1:nElems,1:nSpecies)%PartV(1)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(2)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV(3)  = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(1) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(2) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartV2(3) = 0
  SampDSMC(1:nElems,1:nSpecies)%PartNum   = 0
  SampDSMC(1:nElems,1:nSpecies)%SimPartNum= 0
  SampDSMC(1:nElems,1:nSpecies)%ERot      = 0
  SampDSMC(1:nElems,1:nSpecies)%EVib      = 0
  SampDSMC(1:nElems,1:nSpecies)%EElec     = 0

  ALLOCATE(CollMean(nElems,2))
  CollMean(1:nElems,1:2)=0

! definition of DSMC particle values
  ALLOCATE(DSMC_RHS(PDM%maxParticleNumber,3))
  DSMC_RHS = 0


  IF (nSpecies.LE.0) THEN 
    SWRITE(*,*) "ERROR: nSpecies .LE. 0:", nSpecies
    STOP
  END IF

! reading species data of ini_2
  ALLOCATE(SpecDSMC(nSpecies))
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I2)') iSpec
    SpecDSMC(iSpec)%Name    = GETSTR('Part-Species'//TRIM(hilf)//'-SpeciesName','none')
    SpecDSMC(iSpec)%InterID = GETINT('Part-Species'//TRIM(hilf)//'-InteractionID','0')
    SpecDSMC(iSpec)%TrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceTemp','0')
    SpecDSMC(iSpec)%DrefVHS = GETREAL('Part-Species'//TRIM(hilf)//'-VHSReferenceDiam','0')
    IF((SpecDSMC(iSpec)%InterID*SpecDSMC(iSpec)%TrefVHS*SpecDSMC(iSpec)%DrefVHS).eq.0) THEN
      SWRITE(*,*) "ERROR in species data ini_2"
      STOP
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
  DSMC%CollProbMax = 0   ! Reset of maximal Collision Probability


!-----------------------------------------------------------------------------------------------------------------------------------
! reading/writing molecular stuff
!-----------------------------------------------------------------------------------------------------------------------------------
IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN ! perform relaxation (molecular) reactions
  ! allocate internal enery arrays
  IF ( DSMC%ElectronicState ) THEN
    ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,3))
  ELSE
    ALLOCATE(PartStateIntEn(PDM%maxParticleNumber,2))
  ENDIF
  ! reading molecular stuff
  SpecDSMC(1:nSpecies)%Xi_Rot = 0
  SpecDSMC(1:nSpecies)%MaxVibQuant = 0
  SpecDSMC(1:nSpecies)%CharaTVib = 0
  SpecDSMC(1:nSpecies)%Telec = 0
  SpecDSMC(1:nSpecies)%PolyatomicMol=.false.
  SpecDSMC(1:nSpecies)%SpecToPolyArray = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
     WRITE(UNIT=hilf,FMT='(I2)') iSpec
     
     SpecDSMC(iSpec)%PolyatomicMol=GETLOGICAL('Part-Species'//TRIM(hilf)//'-PolyatomicMol','.FALSE.')
     IF(SpecDSMC(iSpec)%PolyatomicMol.AND.DSMC%ElectronicState)  THEN
        SWRITE(*,*) '! Simulation of Polyatomic Molecules and Electronic States are not possible yet!!!'
        STOP
     END IF
     IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        DSMC%NumPolyatomMolecs = DSMC%NumPolyatomMolecs + 1
        SpecDSMC(iSpec)%SpecToPolyArray = DSMC%NumPolyatomMolecs
     ELSE
    
       SpecDSMC(iSpec)%TVib       = GETREAL('Part-Species'//TRIM(hilf)//'-TempVib','0.')
       SpecDSMC(iSpec)%TRot       = GETREAL('Part-Species'//TRIM(hilf)//'-TempRot','0.')  
       SpecDSMC(iSpec)%Xi_Rot     = 2
       SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib','0.')  
       SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')
       IF(SpecDSMC(iSpec)%Ediss_eV*SpecDSMC(iSpec)%CharaTVib*SpecDSMC(iSpec)%TRot*SpecDSMC(iSpec)%TVib.EQ.0) THEN
         SWRITE(*,*) '! =========================================================================== !'
         SWRITE(*,*) "! ERROR in MolecularData of MolecSpec                                         !", iSpec
         SWRITE(*,*) '! =========================================================================== !'
         IF(SpecDSMC(iSpec)%Ediss_eV*SpecDSMC(iSpec)%CharaTVib.EQ.0) THEN
           STOP
         ELSE
           IF (DSMC%VibEnergyModel.EQ.0) THEN
             SpecDSMC(iSpec)%MaxVibQuant = 200
           ELSE
             SpecDSMC(iSpec)%MaxVibQuant = INT(SpecDSMC(iSpec)%Ediss_eV*JToEv/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)) + 1
           END IF
         END IF
       ELSE
         IF (DSMC%VibEnergyModel.EQ.0) THEN
           SpecDSMC(iSpec)%MaxVibQuant = 200
         ELSE
           SpecDSMC(iSpec)%MaxVibQuant = INT(SpecDSMC(iSpec)%Ediss_eV*JToEv/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)) + 1
         END IF
       END IF
      END IF
      ! read electronic temperature
      IF ( DSMC%ElectronicState ) THEN
        WRITE(UNIT=hilf,FMT='(I2)') iSpec
        SpecDSMC(iSpec)%Telec           = GETREAL('Part-Species'//TRIM(hilf)//'-Tempelec','0.')
      END IF
      SpecDSMC(iSpec)%VFD_Phi3_Factor = GETREAL('Part-Species'//TRIM(hilf)//'-VFDPhi3','0.')
      ! Setting the values of Rot-/Vib-RelaxProb to a fix value!
      ! This should be changed to a calculated value for every coll pair/situation!!!1
      SpecDSMC(iSpec)%RotRelaxProb  = 0.2!0.2
      SpecDSMC(iSpec)%VibRelaxProb  = 0.02!0.02
      SpecDSMC(iSpec)%ElecRelaxProb = 0.01!or 0.02 | Bird: somewhere in range 0.01 .. 0.02
      ! multi init stuff #
      IF(Species(iSpec)%NumberOfInits.NE.0) THEN
        ALLOCATE(SpecDSMC(iSpec)%Init(1:Species(iSpec)%NumberOfInits))
        DO iInit = 1, Species(iSpec)%NumberOfInits
          IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
            WRITE(UNIT=hilf2,FMT='(I2)') iInit
            SpecDSMC(iSpec)%Init(iInit)%TVib = GETREAL('Part-Species'//TRIM(hilf)//'-Init'//TRIM(hilf2)//'-TempVib','0.')
            SpecDSMC(iSpec)%Init(iInit)%TRot = GETREAL('Part-Species'//TRIM(hilf)//'-Init'//TRIM(hilf2)//'-TempRot','0.')
          END IF
        END DO
        ! temp copy
        SpecVib(iSpec) = SpecDSMC(iSpec)%TVib
        SpecRot(iSpec) = SpecDSMC(iSpec)%TRot
      END IF
    END IF
  END DO

  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    ALLOCATE(VibQuantsPar(PDM%maxParticleNumber))
    ALLOCATE(PolyatomMolDSMC(DSMC%NumPolyatomMolecs))
    DO iSpec = 1, nSpecies
      IF (SpecDSMC(iSpec)%PolyatomicMol) CALL InitPolyAtomicMolecs(iSpec)
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
        SpecDSMC(iSpec)%MaxElecQuant         = size( SpecDSMC(iSpec)%ElectronicState,2) ! max Quant + 1
      END IF
    END DO
  END IF
#endif
!CALL DSMC_SetInternalEnr_PolyFast(1,PDM%ParticleVecLength)
  DO iPart = 1, PDM%ParticleVecLength
    IF (PDM%ParticleInside(ipart)) THEN
      IF (Species(PartSpecies(iPart))%NumberOfInits.EQ.0) THEN
        IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
          CALL DSMC_SetInternalEnr_PolyFastPart2(PartSpecies(iPart),iPart)
        ELSE
          CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iPart)
        END IF
      ELSE
        iInit = PDM%PartInit(iPart)
        IF(SpecDSMC(PartSpecies(iPart))%InterID.EQ.2) THEN
          SpecDSMC(PartSpecies(iPart))%TVib = SpecDSMC(PartSpecies(iPart))%Init(iInit)%TVib
          SpecDSMC(PartSpecies(iPart))%TRot = SpecDSMC(PartSpecies(iPart))%Init(iInit)%TRot
          CALL DSMC_SetInternalEnr_LauxVFD(PartSpecies(iPart),iPart)
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
            WRITE(483,'(I3.1,3x,F12.7)') ii, REAL( SpecDSMC(iSpec)%levelcounter(ii) ) / REAL( Species(iSpec)%initialParticleNumber )
          END DO
          close(unit=483)
        END IF
      END IF
    END DO
  END IF
#endif
  DO iSpec = 1, nSpecies
    IF(Species(iSpec)%NumberOfInits.NE.0) THEN
      IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
        ! copy it back
        SpecDSMC(iSpec)%TVib = SpecVib(iSpec)
        SpecDSMC(iSpec)%TRot = SpecRot(iSpec)
      END IF
    END IF
  END DO
  DEALLOCATE(PDM%PartInit)
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
! reading/writing BG Gas stuff
!-----------------------------------------------------------------------------------------------------------------------------------
BGGas%BGGasSpecies  = GETINT('Particles-DSMCBackgroundGas','0')
BGGas%BGGasDensity  = GETREAL('Particles-DSMCBackgroundGasDensity','0')

!Set mean VibQua of BGGas for dissoc reaction
IF (BGGas%BGGasSpecies.NE.0) THEN
  IF(SpecDSMC(BGGas%BGGasSpecies)%InterID.EQ.2) THEN
    BGGas%BGMeanEVibQua = DSMC%GammaQuant * BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib & 
                        + BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib  &
                        /  (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib / SpecDSMC(BGGas%BGGasSpecies)%TVib) - 1) &
                        - BoltzmannConst * SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                        / (EXP(SpecDSMC(BGGas%BGGasSpecies)%CharaTVib * SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant &
                        / SpecDSMC(BGGas%BGGasSpecies)%TVib) - 1)
    BGGas%BGMeanEVibQua = BGGas%BGMeanEVibQua/(BoltzmannConst*SpecDSMC(BGGas%BGGasSpecies)%CharaTVib) - DSMC%GammaQuant
    BGGas%BGMeanEVibQua = MIN(INT(BGGas%BGMeanEVibQua) + 1, SpecDSMC(BGGas%BGGasSpecies)%MaxVibQuant)    
  ELSE
    BGGas%BGMeanEVibQua = 0
  END IF
END IF

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
SWRITE(UNIT_stdOut,'(A)')' INIT DSMC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitDSMC

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


SUBROUTINE DSMC_SetInternalEnr_LauxVFD(iSpecies, iPart)

  USE MOD_DSMC_Vars,            ONLY : PartStateIntEn, SpecDSMC, DSMC
  USE MOD_Particle_Vars,        ONLY : BoltzmannConst
  USE MOD_DSMC_ElectronicModel, ONLY : InitElectronShell
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
  INTEGER, INTENT(IN)           :: iSpecies, iPart
  REAL                          :: iRan, iRan2
  REAL                          :: ElectronicPartitionTemp, ElectronicPartition, summ
  INTEGER                       :: iQuant, ii ! maximal Quant of Species
!#ifdef MPI
!#endif
!===================================================================================================================================

!set vibrational energy
IF (SpecDSMC(iSpecies)%InterID.EQ.2) THEN
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT(-LOG(iRan)*SpecDSMC(iSpecies)%TVib/SpecDSMC(iSpecies)%CharaTVib)
  DO WHILE (iQuant.GE.SpecDSMC(iSpecies)%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*SpecDSMC(iSpecies)%TVib/SpecDSMC(iSpecies)%CharaTVib)
  END DO
  !evtl muß partstateinten nochmal geändert werden, mpi, resize etc..
  PartStateIntEn(iPart, 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpecies)%CharaTVib*BoltzmannConst
!set rotational energy
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn(iPart, 2) = -BoltzmannConst*SpecDSMC(iSpecies)%TRot*LOG(iRan)
ELSE
  PartStateIntEn(iPart, 1) = 0
  PartStateIntEn(iPart, 2) = 0
END IF

! set electronic energy
IF ( DSMC%ElectronicState .and. SpecDSMC(iSpecies)%InterID .ne. 4) THEN
  CALL InitElectronShell(iSpecies,iPart)
ENDIF

END SUBROUTINE DSMC_SetInternalEnr_LauxVFD

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


SUBROUTINE DSMC_BuildSurfaceOutputMapping()

USE MOD_Mesh_Vars,          ONLY:nElems,nBCSides, SideToElem, BC
USE MOD_Particle_Vars,      ONLY:GEO, PartBound
USE MOD_DSMC_Vars,          ONLY:SurfMesh, SampWall
!--------------------------------------------------------------------------------------------------!
! perform Mapping for Surface Output
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
INTEGER                 :: iElem, iLocSide, iSide, iNode, iNode2
INTEGER                 :: SideID
INTEGER, ALLOCATABLE    :: TempBCSurfNodes(:), TempSideSurfNodeMap(:,:)
REAL,ALLOCATABLE        :: TempSurfaceArea(:)
LOGICAL                 :: IsSortedSurfNode
!#ifdef MPI
!#endif
!===================================================================================================================================
ALLOCATE(TempBCSurfNodes(4*nBCSides))
ALLOCATE(TempSideSurfNodeMap(1:4,1:nBCSides))
ALLOCATE(SurfMesh%GlobSideToSurfSideMap(nBCSides))
ALLOCATE(TempSurfaceArea(nBCSides))
SurfMesh%nSurfaceNode=0
SurfMesh%nSurfaceBCSides=0
SurfMesh%GlobSideToSurfSideMap(1:nBCSides)=0

DO iSide=1, nBCSides
  IF (PartBound%Map(BC(iSide)).EQ.PartBound%ReflectiveBC) THEN  
    SurfMesh%nSurfaceBCSides = SurfMesh%nSurfaceBCSides + 1
    SurfMesh%GlobSideToSurfSideMap(iSide) = SurfMesh%nSurfaceBCSides
    iElem = SideToElem(1,iSide)
    IF (iElem.LT.1) THEN
      iElem = SideToElem(2,iSide)
      iLocSide = SideToElem(4,iSide)
    ELSE
      iLocSide = SideToElem(3,iSide)
    END IF
    TempSurfaceArea(SurfMesh%nSurfaceBCSides) = CalcArea(iLocSide, iElem)
    DO iNode2 = 1, 4
    IsSortedSurfNode = .false.
      DO iNode = 1, SurfMesh%nSurfaceNode 
        IF (GEO%ElemSideNodeID(iNode2, iLocSide, iElem).EQ.TempBCSurfNodes(iNode)) THEN
        TempSideSurfNodeMap(iNode2,SurfMesh%nSurfaceBCSides) = iNode
        IsSortedSurfNode = .true.
        EXIT
        END IF
      END DO
      IF(.NOT.IsSortedSurfNode) THEN
        SurfMesh%nSurfaceNode = SurfMesh%nSurfaceNode + 1
        TempBCSurfNodes(SurfMesh%nSurfaceNode) = GEO%ElemSideNodeID(iNode2, iLocSide, iElem)
        TempSideSurfNodeMap(iNode2,SurfMesh%nSurfaceBCSides) = SurfMesh%nSurfaceNode
      END IF
    END DO  
  END IF
END DO

ALLOCATE(SurfMesh%BCSurfNodes(1:SurfMesh%nSurfaceNode))
SurfMesh%BCSurfNodes(1:SurfMesh%nSurfaceNode) = TempBCSurfNodes(1:SurfMesh%nSurfaceNode)
ALLOCATE(SurfMesh%SideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides))
SurfMesh%SideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides) = TempSideSurfNodeMap(1:4,1:SurfMesh%nSurfaceBCSides)
ALLOCATE(SurfMesh%SurfaceArea(1:SurfMesh%nSurfaceBCSides))
SurfMesh%SurfaceArea(1:SurfMesh%nSurfaceBCSides)=TempSurfaceArea(1:SurfMesh%nSurfaceBCSides)
ALLOCATE(SampWall(1:SurfMesh%nSurfaceBCSides))
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(1) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(2) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(3) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(4) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(5) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(6) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(7) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(8) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Energy(9) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Force(1) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Force(2) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Force(3) = 0.0
SampWall(1:SurfMesh%nSurfaceBCSides)%Counter(1) = 0.0
DEALLOCATE(TempBCSurfNodes)
DEALLOCATE(TempSideSurfNodeMap)
DEALLOCATE(TempSurfaceArea)

END SUBROUTINE DSMC_BuildSurfaceOutputMapping

!---------------------------------------------------------------------------------------!
#ifdef MPI
SUBROUTINE DSMC_BuildHaloSurfaceOutputMapping()

USE MOD_Mesh_Vars,          ONLY : nElems,nBCSides, SideToElem, BC
USE MOD_Particle_Vars,      ONLY : GEO, PartBound
#ifdef MPI
USE MOD_DSMC_Vars,          ONLY : SampWallHaloCell
#endif
USE MOD_DSMC_Vars,          ONLY : SurfMesh
USE MOD_part_MPI_Vars,      ONLY : MPIGEO
!--------------------------------------------------------------------------------------------------!
! perform Mapping for Surface Output
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
INTEGER                 :: iSide

!===================================================================================================================================
ALLOCATE(SurfMesh%HaloSideIDToSurfSideMap(SIZE(MPIGEO%BC,2)))
SurfMesh%HaloSideIDToSurfSideMap(1:SIZE(MPIGEO%BC,2))=0
SurfMesh%nHaloSurfaceBCSides = 0
DO iSide=1, SIZE(MPIGEO%BC,2)
  IF ((MPIGEO%BC(1,iSide).NE.0).AND.(MPIGEO%BC(1,iSide).NE.-1).AND.(MPIGEO%BC(1,iSide).NE.424242)) THEN
    IF (PartBound%Map(MPIGEO%BC(1,iSide)).EQ.PartBound%ReflectiveBC) THEN
      SurfMesh%nHaloSurfaceBCSides = SurfMesh%nHaloSurfaceBCSides + 1
      SurfMesh%HaloSideIDToSurfSideMap(iSide) = SurfMesh%nHaloSurfaceBCSides
    END IF
  END IF
END DO

#ifdef MPI
ALLOCATE(SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides))
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(1) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(2) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(3) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(4) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(5) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(6) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(7) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(8) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(9) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(1) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(2) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(3) = 0.0
SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Counter(1) = 0.0
#endif

END SUBROUTINE DSMC_BuildHaloSurfaceOutputMapping
#endif
!---------------------------------------------------------------------------------------!

REAL FUNCTION CalcArea(iLocSide, Element)
  USE MOD_Particle_Vars,          ONLY : GEO
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  INTEGER, INTENT(IN)         :: iLocSide, Element
  REAL                        :: xNod1, xNod2, xNod3, xNod4 
  REAL                        :: yNod1, yNod2, yNod3, yNod4 
  REAL                        :: zNod1, zNod2, zNod3, zNod4 
  REAL                        :: Vector1(1:3), Vector2(1:3), Vector3(1:3)

!--------------------------------------------------------------------------------------------------!

   xNod1 = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNod1 = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNod1 = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

   xNod2 = GEO%NodeCoords(1,GEO%ElemSideNodeID(2,iLocSide,Element))
   yNod2 = GEO%NodeCoords(2,GEO%ElemSideNodeID(2,iLocSide,Element))
   zNod2 = GEO%NodeCoords(3,GEO%ElemSideNodeID(2,iLocSide,Element))

   xNod3 = GEO%NodeCoords(1,GEO%ElemSideNodeID(3,iLocSide,Element))
   yNod3 = GEO%NodeCoords(2,GEO%ElemSideNodeID(3,iLocSide,Element))
   zNod3 = GEO%NodeCoords(3,GEO%ElemSideNodeID(3,iLocSide,Element))

   xNod4 = GEO%NodeCoords(1,GEO%ElemSideNodeID(4,iLocSide,Element))
   yNod4 = GEO%NodeCoords(2,GEO%ElemSideNodeID(4,iLocSide,Element))
   zNod4 = GEO%NodeCoords(3,GEO%ElemSideNodeID(4,iLocSide,Element))

   Vector1(1) = xNod2 - xNod1
   Vector1(2) = yNod2 - yNod1
   Vector1(3) = zNod2 - zNod1

   Vector2(1) = xNod3 - xNod1
   Vector2(2) = yNod3 - yNod1
   Vector2(3) = zNod3 - zNod1

   Vector3(1) = xNod4 - xNod1
   Vector3(2) = yNod4 - yNod1
   Vector3(3) = zNod4 - zNod1

   CalcArea = 0.5*(SQRT((Vector1(2)*Vector2(3)-Vector1(3)*Vector2(2))**2 &
           + (-Vector1(1)*Vector2(3)+Vector1(3)*Vector2(1))**2 &
           + (Vector1(1)*Vector2(2)-Vector1(2)*Vector2(1))**2) &
           + SQRT((Vector3(2)*Vector2(3)-Vector3(3)*Vector2(2))**2 &
           + (-Vector3(1)*Vector2(3)+Vector3(3)*Vector2(1))**2 &
           + (Vector3(1)*Vector2(2)-Vector3(2)*Vector2(1))**2))

  RETURN

END FUNCTION CalcArea


END MODULE MOD_DSMC_Init
