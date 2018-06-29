#include "boltzplatz.h"

PROGRAM Boltzplatz
!===================================================================================================================================
! Control program of the Boltzplatz code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MOD_Globals_vars           ,ONLY: InitializationWallTime
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ParameterFile,ParameterDSMCFile
USE MOD_Commandline_Arguments
USE MOD_ReadInTools            ,ONLY: prms,PrintDefaultparameterFile,ExtractparameterFile
USE MOD_Boltzplatz_Init        ,ONLY: InitBoltzplatz,FinalizeBoltzplatz
USE MOD_Restart_Vars           ,ONLY: RestartFile
USE MOD_Restart                ,ONLY: Restart
USE MOD_Interpolation          ,ONLY: InitInterpolation
USE MOD_IO_HDF5                ,ONLY: InitIO
USE MOD_TimeDisc               ,ONLY: InitTimeDisc,FinalizeTimeDisc,TimeDisc
USE MOD_MPI                    ,ONLY: InitMPI
USE MOD_RecordPoints_Vars      ,ONLY: RP_Data
USE MOD_Mesh_Vars              ,ONLY: DoSwapMesh
USE MOD_Mesh                   ,ONLY: SwapMesh
#ifdef MPI
USE MOD_LoadBalance            ,ONLY: InitLoadBalance,FinalizeLoadBalance
USE MOD_MPI                    ,ONLY: FinalizeMPI
#endif /*MPI*/
USE MOD_Output                 ,ONLY: InitOutput
USE MOD_Define_Parameters_Init ,ONLY: InitDefineParameters
USE MOD_StringTools            ,ONLY: STRICMP, GetFileExtension
#ifdef PARTICLES
USE MOD_Particle_Vars          ,ONLY: DoInitialIonization
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: Time
LOGICAL                 :: userblockFound
!===================================================================================================================================

CALL InitMPI()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')&
 "           ____            ___    __                    ___              __              "
SWRITE(UNIT_stdOut,'(A)')&
 "          /\\  _`\\         /\\_ \\  /\\ \\__                /\\_ \\            /\\ \\__           "
SWRITE(UNIT_stdOut,'(A)')&
 "          \\ \\ \\L\\ \\    ___\\//\\ \\ \\ \\ ,_\\  ____    _____\\//\\ \\       __  \\ \\ ,_\\  ____    "
SWRITE(UNIT_stdOut,'(A)')&
 "           \\ \\  _ <'  / __`\\\\ \\ \\ \\ \\ \\/ /\\_ ,`\\ /\\ '__`\\\\ \\ \\    /'__`\\ \\ \\ \\/ /\\_ ,`\\  "
SWRITE(UNIT_stdOut,'(A)')&
 "            \\ \\ \\L\\ \\/\\ \\L\\ \\\\_\\ \\_\\ \\ \\_\\/_/  /_\\ \\ \\L\\ \\\\_\\ \\_ /\\ \\L\\.\\_\\ \\ \\_\\/_/  /_ "
SWRITE(UNIT_stdOut,'(A)')&
 "             \\ \\____/\\ \\____//\\____\\\\ \\__\\ /\\____\\\\ \\ ,__//\\____\\\\ \\__/.\\_\\\\ \\__\\ /\\____\\"
SWRITE(UNIT_stdOut,'(A)')&
 "              \\/___/  \\/___/ \\/____/ \\/__/ \\/____/ \\ \\ \\/ \\/____/ \\/__/\\/_/ \\/__/ \\/____/"
SWRITE(UNIT_stdOut,'(A)')&
 "                                                    \\ \\_\\                                "
SWRITE(UNIT_stdOut,'(A)')&
 "                                                     \\/_/                                "
SWRITE(UNIT_stdOut,'(A)')&
 ' '
SWRITE(UNIT_stdOut,'(132("="))')

CALL ParseCommandlineArguments()

! Check if the number of arguments is correct
IF ((nArgs.GT.3) .OR. ((nArgs.EQ.0).AND.(doPrintHelp.EQ.0)) ) THEN
  ! Print out error message containing valid syntax
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
    'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
END IF

CALL InitDefineParameters()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

ParameterFile = Args(1)
IF (nArgs.EQ.2) THEN
  ParameterDSMCFile = Args(2)
  IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
      'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "h5")) THEN
    RestartFile = ParameterDSMCFile
    ParameterDSMCFile = '' !'no file found'
  END IF
ELSE IF (nArgs.GT.2) THEN
  ParameterDSMCFile = Args(2)
  RestartFile = Args(3)
  IF (STRICMP(GetFileExtension(ParameterDSMCFile), "h5").OR.STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
    ! Print out error message containing valid syntax
    CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
      'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  END IF
ELSE IF (STRICMP(GetFileExtension(ParameterFile), "h5")) THEN
  ! Print out error message containing valid syntax
  !CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: boltzplatz parameter.ini [DSMC.ini] [restart.h5]'// &
  !  'or boltzplatz --help [option/section name] to print help for a single parameter, parameter sections or all parameters.')
  ParameterFile = ".boltzplatz.ini" 
  CALL ExtractParameterFile(Args(1), ParameterFile, userblockFound)
  IF (.NOT.userblockFound) THEN
    CALL CollectiveStop(__STAMP__, "No userblock found in state file '"//TRIM(Args(1))//"'")
  END IF
  RestartFile = Args(1)
END IF

StartTime=BOLTZPLATZTIME()
CALL prms%read_options(ParameterFile)
! Measure init duration
Time=BOLTZPLATZTIME()
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A,I0,A)') ' READING INI DONE! [',Time-StartTime,' sec ] NOW '&
,prms%count_setentries(),' PARAMETERS ARE SET'
SWRITE(UNIT_stdOut,'(132("="))')
! Check if we want to read in DSMC.ini
IF(nArgs.GE.2)THEN
  IF(STRICMP(GetFileExtension(ParameterDSMCFile), "ini")) THEN
    CALL prms%read_options(ParameterDSMCFile,furtherini=.TRUE.)
    ! Measure init duration
    Time=BOLTZPLATZTIME()
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A,F8.2,A,I0,A)') ' READING FURTHER INI DONE! [',Time-StartTime,' sec ] NOW '&
    ,prms%count_setentries(),' PARAMETERS ARE SET'
    SWRITE(UNIT_stdOut,'(132("="))')
  END IF
END IF

CALL InitOutput()
CALL InitIO()

CALL InitGlobals()
#ifdef MPI
CALL InitLoadBalance()
#endif /*MPI*/
! call init routines
! Measure init duration
!StartTime=BOLTZPLATZTIME()

! Initialization
CALL InitInterpolation()
CALL InitTimeDisc()

CALL InitBoltzplatz(IsLoadBalance=.FALSE.)

! Do SwapMesh
IF(DoSwapMesh)THEN
  ! Measure init duration
  Time=BOLTZPLATZTIME()
  IF(MPIroot)THEN
    Call SwapMesh()
    SWRITE(UNIT_stdOut,'(132("="))')
    SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' SWAPMESH DONE! BOLTZPLATZ DONE! [',Time-StartTime,' sec ]'
    SWRITE(UNIT_stdOut,'(132("="))')
    STOP
  ELSE
  CALL abort(&
  __STAMP__&
  ,'DO NOT CALL SWAPMESH WITH MORE THAN 1 Procs!',iError,999.)
  END IF
END IF

! RESTART
CALL Restart()

#ifdef PARTICLES
! Ionize the current particles
IF(DoInitialIonization) CALL InitialIonization()
#endif /*PARTICLES*/

! Measure init duration
Time=BOLTZPLATZTIME()
InitializationWallTime=Time-StartTime
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' INITIALIZATION DONE! [',InitializationWallTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

! Run Simulation
CALL TimeDisc()


!Finalize
CALL FinalizeBoltzplatz(IsLoadBalance=.FALSE.)

CALL FinalizeTimeDisc()
! mssing arrays to deallocate
SDEALLOCATE(RP_Data)

!Measure simulation duration
Time=BOLTZPLATZTIME()

#ifdef MPI
!! and additional required for restart with load balance
!ReadInDone=.FALSE.
!ParticleMPIInitIsDone=.FALSE.
!ParticlesInitIsDone=.FALSE.
CALL FinalizeLoadBalance()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(&
  __STAMP__&
  ,'MPI finalize error',iError,999.)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)')  ' BOLTZPLATZ FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Boltzplatz


SUBROUTINE InitialIonization() 
!----------------------------------------------------------------------------------------------------------------------------------!
! 1.) assign charges to each atom/molecule using the charge supplied by the user
! 2.) reconstruct the electron phase space using the summed charged per cell for which an electron is 
!     created to achieve an ionization degree supplied by the user
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: ElectronCharge
USE MOD_Particle_Vars ,ONLY: PDM,PEM,PartState,nSpecies,Species,PartSpecies
USE MOD_Eval_xyz      ,ONLY: eval_xyz_elemcheck
USE MOD_Mesh_Vars     ,ONLY: NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_DSMC_Vars     ,ONLY: CollisMode,DSMC,PartStateIntEn
USE MOD_part_emission ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_Vars     ,ONLY: useDSMC
USE MOD_Eval_xyz      ,ONLY: Eval_XYZ_Poly
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: ChargeLower,ChargeUpper
INTEGER :: ElemCharge(1:PP_nElems)
INTEGER :: ElecSpecIndx,iSpec,location,iElem,iPart,ParticleIndexNbr
REAL    :: ChargeProbability
REAL    :: iRan
REAL    :: PartPosRef(1:3)
REAL    :: CellElectronTemperature=300
!REAL    :: MaxElectronTemp_eV
INTEGER :: SpeciesID(3)
INTEGER :: SpeciesCharge(3)
REAL    :: ChargeAverage
!===================================================================================================================================
SpeciesID     = (/1,2,3/)    ! Species IDs for the considered species
SpeciesCharge = (/0,-1,1/)   ! Charge state for each species
ChargeAverage = 0.01         ! Average charge for each atom/molecule in the cell (corresponds to the ionization degree)

! 1.) reconstruct ions and determine charge
SWRITE(UNIT_stdOut,*)'1.) Reconstructing ions and determining charge'
ElemCharge(1:PP_nElems)=0
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart)) THEN
    IF(ANY(PartSpecies(iPart).EQ.SpeciesID(:)))THEN ! 
      ! get the cell charge average and select and upper and lower charge number
      ChargeLower       = ChargeAverage !INT(TTM_Cell_11(PEM%Element(iPart))) ! use first DOF (0,0,0) because the data is const. in each cell
      ChargeUpper       = ChargeLower+1
      ChargeProbability = REAL(ChargeUpper)-ChargeAverage !TTM_Cell_11(PEM%Element(iPart)) ! 2-1,4=0.6 -> 60% probability to get lower charge
      ! distribute the charge using random numbers
      CALL RANDOM_NUMBER(iRan)
      IF(iRan.LT.ChargeProbability)THEN ! select the lower charge number
        location = MINLOC(ABS(SpeciesCharge-ChargeLower),1) !Determines the location of the element in the array with min value
        ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeLower
      ELSE ! select the upper charge number
        location = MINLOC(ABS(SpeciesCharge-ChargeUpper),1) !Determines the location of the element in the array with min value
        ElemCharge(PEM%Element(iPart))=ElemCharge(PEM%Element(iPart))+ChargeUpper
      END IF
      PartSpecies(iPart)=SpeciesID(location) ! set the species ID to atom/singly charged ion/doubly charged ... and so on
    END IF
  END IF
END DO

! 2.) reconstruct electrons
SWRITE(UNIT_stdOut,*)'2.) Reconstructing electrons'
!MaxElectronTemp_eV=MAXVAL(TTM_Cell_2(:))
ElecSpecIndx = -1
DO iSpec = 1, nSpecies
  IF (Species(iSpec)%ChargeIC.GT.0.0) CYCLE
  IF(NINT(Species(iSpec)%ChargeIC/(-1.60217653E-19)).EQ.1) THEN
    ElecSpecIndx = iSpec
    EXIT
  END IF
END DO
IF (ElecSpecIndx.EQ.-1) CALL abort(&
  __STAMP__&
  ,'Electron species not found. Cannot create electrons without the defined species!')

DO iElem=1,PP_nElems
  DO iPart=1,ElemCharge(iElem) ! 1 electron for each charge of each element
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
    ParticleIndexNbr            = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    PDM%ParticleVecLength       = PDM%ParticleVecLength + 1
    !Set new Species of new particle
    PDM%ParticleInside(ParticleIndexNbr) = .true.
    PartSpecies(ParticleIndexNbr) = ElecSpecIndx
     
    CALL RANDOM_NUMBER(PartPosRef(1:3)) ! get random reference space
    PartPosRef(1:3)=PartPosRef(1:3)*2. - 1. ! map (0,1) -> (-1,1)
    CALL Eval_xyz_Poly(PartPosRef(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem) &
                      ,PartState(ParticleIndexNbr,1:3)) !Map into phys. space

    IF ((useDSMC).AND.(CollisMode.GT.1)) THEN
      PartStateIntEn(ParticleIndexNbr, 1) = 0.
      PartStateIntEn(ParticleIndexNbr, 2) = 0.
      IF ( DSMC%ElectronicModel )  PartStateIntEn(ParticleIndexNbr, 3) = 0.
    END IF
    PEM%Element(ParticleIndexNbr) = iElem
    ! IF(TTM_Cell_2(iElem).LE.0.0)THEN ! not enough atoms in FD cell for averaging a temperature: use max value for electrons
    !   CellElectronTemperature=(MaxElectronTemp_eV*ElectronCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    ! ELSE
    !   CellElectronTemperature=(TTM_Cell_2(iElem)*ElectronCharge/BoltzmannConst) ! convert eV to K: 1 [eV] = e/kB [K]
    ! END IF
    CALL CalcVelocity_maxwell_lpn(ElecSpecIndx, PartState(ParticleIndexNbr,4:6),&
                                  Temperature=CellElectronTemperature)
  END DO
END DO


END SUBROUTINE InitialIonization
