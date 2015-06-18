#include "boltzplatz.h"

MODULE MOD_Analyze
!===================================================================================================================================
! Contains DG analyze 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

! GLOBAL VARIABLES 
INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE CalcError
  MODULE PROCEDURE CalcError
END INTERFACE

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE

INTERFACE PerformAnalyze
  MODULE PROCEDURE PerformAnalyze
END INTERFACE

!===================================================================================================================================
PUBLIC:: CalcError, InitAnalyze, FinalizeAnalyze, PerformAnalyze 
!===================================================================================================================================

CONTAINS

SUBROUTINE InitAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Interpolation_Vars,   ONLY: xGP,wBary,InterpolationInitIsDone
USE MOD_Analyze_Vars,         ONLY:Nanalyze,AnalyzeInitIsDone,Analyze_dt
USE MOD_ReadInTools,          ONLY:GETINT,GETREAL
USE MOD_Analyze_Vars,         ONLY:CalcPoyntingInt
USE MOD_AnalyzeField,         ONLY:GetPoyntingIntPlane
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)                :: DefStr
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) THEN
  CALL abort(__STAMP__,'InitAnalyse not ready to be called or already called.',999,999.)
  RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'
WRITE(DefStr,'(i4)') 2*(PP_N+1)
NAnalyze=GETINT('NAnalyze',DefStr) 
CALL InitAnalyzeBasis(PP_N,NAnalyze,xGP,wBary)
Analyze_dt=GETREAL('Analyze_dt','0.')
AnalyzeInitIsDone=.TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

! init Poynting-Integral
IF(CalcPoyntingInt) CALL GetPoyntingIntPlane()

END SUBROUTINE InitAnalyze


SUBROUTINE InitAnalyzeBasis(N_in,Nanalyze_in,xGP,wBary)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars, ONLY:wAnalyze,Vdm_GaussN_NAnalyze
USE MOD_Basis,        ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: N_in,Nanalyze_in
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP,wBary
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL ,DIMENSION(0:Nanalyze_in) :: XiAnalyze
!===================================================================================================================================
  ALLOCATE(wAnalyze(0:NAnalyze_in),Vdm_GaussN_NAnalyze(0:NAnalyze_in,0:N_in))
  CALL LegGaussLobNodesAndWeights(NAnalyze_in,XiAnalyze,wAnalyze)
  CALL InitializeVandermonde(N_in,NAnalyze_in,wBary,xGP,XiAnalyze,Vdm_GaussN_NAnalyze)
END SUBROUTINE InitAnalyzeBasis


SUBROUTINE CalcError(Time)
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars,ONLY:tEnd, iter
USE MOD_Mesh_Vars,ONLY:Elem_xGP,sJ
USE MOD_Equation_Vars,ONLY:IniExactFunc
USE MOD_Analyze_Vars,ONLY:NAnalyze,Vdm_GaussN_NAnalyze,wAnalyze
USE MOD_DG_Vars,ONLY:U
USE MOD_Equation,ONLY:ExactFunc
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
USE MOD_Output_Vars,ONLY:ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                       :: iElem,k,l,m
REAL                          :: CalcTime
REAL                          :: L_Inf_Error(PP_nVar),L_2_Error(PP_nVar),U_exact(PP_nVar),L_2_Error2(PP_nVar),L_Inf_Error2(PP_nVar)
REAL                          :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                          :: Volume,IntegrationWeight
#ifdef MPI
REAL                          :: Volume2
#endif
CHARACTER(LEN=40)             :: formatStr
!===================================================================================================================================
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.
L_Inf_Error2(:)=-1.E10
L_2_Error2(:)=0.
Volume=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,PP_nElems
   ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
   CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
   CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
   ! Interpolate the solution to the analyze grid
   CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
   DO m=0,NAnalyze
     DO l=0,NAnalyze
       DO k=0,NAnalyze
         CALL ExactFunc(IniExactFunc,time,0,Coords_NAnalyze(1:3,k,l,m),U_exact)
         L_Inf_Error = MAX(L_Inf_Error,abs(U_NAnalyze(:,k,l,m) - U_exact))
         IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
         ! To sum over the elements, We compute here the square of the L_2 error
         L_2_Error = L_2_Error+(U_NAnalyze(:,k,l,m) - U_exact)*(U_NAnalyze(:,k,l,m) - U_exact)*IntegrationWeight
         Volume = Volume + IntegrationWeight 
       END DO ! k
     END DO ! l
   END DO ! m
END DO ! iElem=1,PP_nElems
#ifdef MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,L_2_Error,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
    CALL MPI_REDUCE(MPI_IN_PLACE,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
    CALL MPI_REDUCE(MPI_IN_PLACE,L_Inf_Error,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error,L_2_Error2,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
    CALL MPI_REDUCE(volume,volume2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
    CALL MPI_REDUCE(L_Inf_Error,L_Inf_Error2,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
    ! in this case the receive value is not relevant. 
  END IF
#endif /*MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/Volume)

! Graphical output
CalcTime=BOLTZPLATZTIME()
IF(MPIroot) THEN
  WRITE(UNIT_StdOut,'(A13,ES16.7)')' Sim time  : ',Time
  WRITE(formatStr,'(A5,I1,A7)')'(A13,',PP_nVar,'ES16.7)'
  WRITE(UNIT_StdOut,formatStr)' L_2       : ',L_2_Error
  WRITE(UNIT_StdOut,formatStr)' L_inf     : ',L_Inf_Error
  IF (Time.GE.tEnd) CALL AnalyzeToFile(Time,CalcTime,iter,L_2_Error)
  IF (Time.GT.0.) THEN
    WRITE(UNIT_StdOut,'(132("."))')
    WRITE(UNIT_stdOut,'(A,A,A,F8.2,A)') ' BOLTZPLATZ RUNNING ',TRIM(ProjectName),'... [',CalcTime-StartTime,' sec ]'
    WRITE(UNIT_StdOut,'(132("-"))')
    WRITE(UNIT_StdOut,*)
  ELSE
    WRITE(UNIT_StdOut,'(132("="))')
  END IF
END IF
END SUBROUTINE CalcError

SUBROUTINE AnalyzeToFile(Time,CalcTime,iter,L_2_Error)
!===================================================================================================================================
! Writes the L2-error norms to file.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars
USE MOD_Output_Vars  ,ONLY:ProjectName
USE MOD_Mesh_Vars    ,ONLY:nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: TIME                         ! physical time
REAL,INTENT(IN)                :: CalcTime                     ! computational time
INTEGER(KIND=8),INTENT(IN)     :: iter                         ! number of timesteps
REAL,INTENT(IN)                :: L_2_Error(PP_nVar)           ! L2 error norms
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL                           :: Dummyreal(PP_nVar+1),Dummytime  ! Dummy values for file handling
INTEGER                        :: openStat! File IO status
CHARACTER(LEN=50)              :: formatStr                    ! format string for the output and Tecplot header
CHARACTER(LEN=30)              :: L2name(PP_nVar)              ! variable name for the Tecplot header
CHARACTER(LEN=300)             :: Filename                     ! Output filename,
LOGICAL                        :: fileExists                   ! Error handler for file
INTEGER                        :: ioUnit
!===================================================================================================================================
Filename = 'out.'//TRIM(ProjectName)//'.dat'
! Check for file
INQUIRE(FILE = Filename, EXIST = fileExists)
!! File processing starts here open old and extratct information or create new file.
ioUnit=1746
  OPEN(UNIT   = ioUnit       ,&
       FILE   = Filename     ,&
       STATUS = 'Unknown'    ,&
       ACCESS = 'SEQUENTIAL' ,&
       IOSTAT = openStat                 )
  IF (openStat.NE.0) THEN
     WRITE(*,*)'ERROR: cannot open Outfile'
  END IF
  ! Create a new file with the Tecplot header etc. 
  WRITE(ioUnit,*)'TITLE="Analysis,'//TRIM(ProjectName)//'"'
  WRITE(ioUnit,'(A12)')'VARIABLES ='
  ! Fill the formatStr and L2name strings
  CALL getVARformatStr(formatStr,L2name)
  WRITE(ioUnit,formatStr)'"timesteps"',L2name,' "t_sim" "t_CPU" "DOF" "Ncells" "nProcs"'
  WRITE(ioUnit,*) 'ZONE T="Analysis,'//TRIM(ProjectName)//'"'

! Create format string for the variable output
WRITE(formatStr,'(A10,I1,A37)')'(E23.14E5,',PP_nVar,'(1X,E23.14E5),4(1X,E23.14E5),2X,I6.6)'
WRITE(ioUnit,formatstr) REAL(iter),L_2_Error(:),TIME,CalcTime-StartTime, &
                 REAL(nGlobalElems*(PP_N+1)**3),REAL(nGlobalElems),nProcessors

CLOSE(ioUnit) ! outputfile
END SUBROUTINE AnalyzeToFile

SUBROUTINE getVARformatStr(VARformatStr,L2name)
!===================================================================================================================================
! This creates the format string for writeAnalyse2file
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Equation_Vars,ONLY:StrVarNames
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  ! CALC%varName: Name of conservative variables
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  CHARACTER(LEN=30) :: L2name(PP_nVar) ! The name of the Tecplot variables
  CHARACTER(LEN=50) :: VARformatStr ! L2name format string
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
  INTEGER           :: i ! counter
!===================================================================================================================================
DO i=1,PP_nVar
  WRITE(L2name(i),'(A5,A,A2)')' "L2_',TRIM(StrVarNames(i)),'" '
END DO
WRITE(VARformatStr,'(A3)')'(A,'
DO i=1,PP_nVar
  WRITE(VARformatStr,'(A,A1,I2,A1)')TRIM(VARformatStr),'A',LEN_TRIM(L2name(i)),','
END DO
WRITE(VARformatStr,'(A,A2)')TRIM(VARformatStr),'A)'
END SUBROUTINE getVARformatStr

SUBROUTINE FinalizeAnalyze()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars
USE MOD_AnalyzeField,     ONLY:FinalizePoyntingInt
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NAnalyze)
SDEALLOCATE(wAnalyze)
IF(CalcPoyntingInt) CALL FinalizePoyntingInt()
AnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeAnalyze


SUBROUTINE PerformAnalyze(t,iter,tenddiff,forceAnalyze,OutPut,LastIter)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,             ONLY: nElems
USE MOD_Analyze_Vars,          ONLY: CalcPoyntingInt
USE MOD_AnalyzeField,          ONLY: CalcPoyntingIntegral
USE MOD_RecordPoints,          ONLY: RecordPoints
USE MOD_RecordPoints_Vars,     ONLY: RP_onProc
USE MOD_Particle_Analyze_Vars, ONLY: DoAnalyze, PartAnalyzeStep
#if (PP_TimeDiscMethod==42)
USE MOD_TimeDisc_Vars,         ONLY: dt
#else
USE MOD_TimeDisc_Vars,         ONLY: TEnd,dt
#endif
#ifdef PARTICLES
USE MOD_PARTICLE_Vars,         ONLY: WriteMacroValues,MacroValSamplIterNum,nSpecies, Time
USE MOD_Particle_Analyze,      ONLY: AnalyzeParticles
USE MOD_Particle_Analyze_Vars, ONLY: DoAnalyze, PartAnalyzeStep
USE MOD_DSMC_Vars,             ONLY: SampDSMC,nOutput,DSMC,useDSMC, iter_macvalout,SurfMesh,SampWall
USE MOD_DSMC_Analyze,          ONLY: DSMC_output_calc, DSMC_data_sampling, CalcSurfaceValues, WriteOutputMeshSamp
USE MOD_Particle_surfaces_vars, ONLY: ntracks,tTracking,tLocalization,MeassureTrackTime
#ifdef MPI
USE MOD_Particle_MPI_Vars,     ONLY: PartMPI
#endif /*MPI*/
#if (PP_TimeDiscMethod!=42)
!USE MOD_LD_Vars,               ONLY: useLD
#endif
!USE MOD_LD_Analyze,            ONLY: LD_data_sampling, LD_output_calc
#if (PP_TimeDiscMethod==1001)
!USE MOD_LD_DSMC_TOOLS
#endif
#else
USE MOD_AnalyzeField,          ONLY: AnalyzeField
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)            :: t
REAL,INTENT(IN)               :: tenddiff
INTEGER(KIND=8),INTENT(INOUT) :: iter
LOGICAL,INTENT(IN)            :: forceAnalyze,output
LOGICAL,INTENT(IN),OPTIONAL   :: LastIter
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iSampWallAlloc
#ifdef PARTICLES
INTEGER                       :: RECI
REAL                          :: RECR
#endif /*PARTICLES*/

!===================================================================================================================================

! not for first iteration
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
IF((iter.EQ.0).AND.(.NOT.forceAnalyze)) RETURN
!IF(iter.EQ.0) RETURN
#endif

!----------------------------------------------------------------------------------------------------------------------------------
! DG-Solver
!----------------------------------------------------------------------------------------------------------------------------------

! Calculate error norms
IF(forceAnalyze.OR.Output) CALL CalcError(t)

! poynting vector
IF (CalcPoyntingInt) THEN
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
  IF(forceAnalyze)THEN
    CALL CalcPoyntingIntegral(t,doProlong=.TRUE.)
   ELSE
    IF(MOD(iter,PartAnalyzeStep).EQ.0 .AND. .NOT. OutPut) CALL CalcPoyntingIntegral(t,doProlong=.FALSE.)
  END IF ! ForceAnalyze
  IF(PRESENT(LastIter) .AND. LastIter) CALL CalcPoyntingIntegral(t,doProlong=.TRUE.)
#else
  IF(MOD(iter,PartAnalyzeStep).EQ.0) CALL CalcPoyntingIntegral(t)
#endif
END IF

! fill recordpoints buffer
#if (PP_TimeDiscMethod==1) || (PP_TimeDiscMethod==2) || (PP_TimeDiscMethod==6)
IF(RP_onProc) THEN
  CALL RecordPoints(iter,t,forceAnalyze,Output)
END IF
#endif

!----------------------------------------------------------------------------------------------------------------------------------
! PIC & DG-Sovler
!----------------------------------------------------------------------------------------------------------------------------------

#ifdef PARTICLES
! particle analyze
IF (DoAnalyze)  THEN
  IF(forceAnalyze)THEN
    CALL AnalyzeParticles(t) 
  ELSE
    IF(MOD(iter,PartAnalyzeStep).EQ.0 .AND. .NOT. OutPut) CALL AnalyzeParticles(t) 
  END IF
  IF(PRESENT(LastIter) .AND. LastIter) CALL AnalyzeParticles(t) 
END IF

#else /*pure DGSEM */
IF (DoAnalyze)  THEN
  IF(ForceAnalyze)THEN
    CALL AnalyzeField(t) 
  ELSE
    IF(MOD(iter,PartAnalyzeStep).EQ.0 .AND. .NOT. OutPut) CALL AnalyzeField(t) 
  END IF
  IF(PRESENT(LastIter) .AND. LastIter) CALL AnalyzeField(t) 
END IF

  !IF(PartAnalyzeStep.EQ.123456789) CALL AnalyzeParticles(t) 
END IF
#endif /*PARTICLES*/

!----------------------------------------------------------------------------------------------------------------------------------
! DSMC & LD 
!----------------------------------------------------------------------------------------------------------------------------------
! update of time here
#ifdef PARTICLES
Time = t
! write DSMC macroscopic values 
IF ((WriteMacroValues).AND.(.NOT.Output))THEN
#if (PP_TimeDiscMethod==1000)
  CALL LD_data_sampling()  ! Data sampling for output
#elif(PP_TimeDiscMethod==1001)
  CALL LD_DSMC_data_sampling()
#else
  CALL DSMC_data_sampling()
#endif
  iter_macvalout = iter_macvalout + 1
  IF (MacroValSamplIterNum.LE.iter_macvalout) THEN
#if (PP_TimeDiscMethod==1000)
    CALL LD_output_calc()  ! Data sampling for output
#elif(PP_TimeDiscMethod==1001)
    CALL LD_DSMC_output_calc()
#else
    CALL DSMC_output_calc
    IF (DSMC%CalcSurfaceVal) CALL CalcSurfaceValues
    IF (DSMC%OutputMeshSamp) CALL WriteOutputMeshSamp() !EmType6
#endif
    nOutput = nOutput + 1
    iter_macvalout = 0
    DSMC%SampNum = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV(1)  = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV(2)  = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV(3)  = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV2(1) = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV2(2) = 0
    SampDSMC(1:nElems,1:nSpecies)%PartV2(3) = 0
    SampDSMC(1:nElems,1:nSpecies)%PartNum   = 0
    SampDSMC(1:nElems,1:nSpecies)%ERot      = 0
    SampDSMC(1:nElems,1:nSpecies)%EVib      = 0

    IF(DSMC%CalcSurfaceVal) THEN
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
      DO iSampWallAlloc=1,SurfMesh%nSurfaceBCSides
        SampWall(iSampWallAlloc)%Counter(1:nSpecies) = 0.0
      END DO
    END IF

  END IF
END IF

IF(OutPut)THEN
#if (PP_TimeDiscMethod==42)
  IF((dt.EQ.tEndDiff).AND.(useDSMC).AND.(.NOT.DSMC%ReservoirSimu)) THEN
    CALL DSMC_output_calc
  END IF
#else
  IF((dt.EQ.tEndDiff).AND.(useDSMC).AND.(.NOT.WriteMacroValues)) THEN
    nOutput = INT((DSMC%TimeFracSamp * TEnd) / DSMC%DeltaTimeOutput)
    !IF (.NOT. useLD) THEN
      CALL DSMC_output_calc
      IF (DSMC%OutputMeshSamp) CALL WriteOutputMeshSamp() !EmType6
    !END IF
    IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues
  END IF
#endif
END IF

! meassure tracking time for particles // no MPI barrier MPI Wall-time but local CPU time
! allows non-synchronous meassurement of particle tracking
IF(OutPut .AND. MeassureTrackTime)THEN
  IF(MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,nTracks      , 1 ,MPI_INTEGER         ,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,tTracking    , 1 ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE,tLocalization, 1 ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  ELSE ! no Root
    CALL MPI_REDUCE(nTracks      ,RECI,1,MPI_INTEGER         ,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_REDUCE(tTracking    ,RECR,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
    CALL MPI_REDUCE(tLocalization,RECR,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  END IF
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A,I15)')   ' Number of trackings:   ',nTracks
  SWRITE(UNIT_stdOut,'(A,F15.6)') ' Tracking time:         ',tTracking
  SWRITE(UNIT_stdOut,'(A,F15.8)') ' Average Tracking time: ',tTracking/REAL(nTracks)
  SWRITE(UNIT_stdOut,'(A,F15.6)') ' Localization time:     ',tLocalization
  SWRITE(UNIT_StdOut,'(132("-"))')
  nTracks=0
  tTracking=0.
  tLocalization=0.
END IF ! only during output like Doftime
#endif /*PARTICLES*/

END SUBROUTINE PerformAnalyze

END MODULE MOD_Analyze
