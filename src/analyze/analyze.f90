!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

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

INTERFACE FinalizeAnalyze
  MODULE PROCEDURE FinalizeAnalyze
END INTERFACE

INTERFACE PerformAnalyze
  MODULE PROCEDURE PerformAnalyze
END INTERFACE

INTERFACE CalcErrorStateFiles
  MODULE PROCEDURE CalcErrorStateFiles
END INTERFACE

INTERFACE CalcErrorStateFileSigma
  MODULE PROCEDURE CalcErrorStateFileSigma
END INTERFACE

!===================================================================================================================================
PUBLIC:: DefineParametersAnalyze
PUBLIC:: CalcError, InitAnalyze, FinalizeAnalyze, PerformAnalyze, CalcErrorStateFiles, CalcErrorStateFileSigma
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyze()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
!USE MOD_AnalyzeEquation ,ONLY: DefineParametersAnalyzeEquation
IMPLICIT NONE
!==================================================================================================================================
! -------------------------
CALL prms%SetSection("Analyze")
! -------------------------
CALL prms%CreateLogicalOption('DoCalcErrorNorms'     , 'Set true to compute L2 and LInf error norms at analyze step.','.FALSE.')
CALL prms%CreateLogicalOption('OutputErrorNormsToH5' , 'Set true to write the analytical solution, the L2 and LInf error norms at analyze step to .h5 state file.','.FALSE.')
CALL prms%CreateRealOption(   'Analyze_dt'           , 'Specifies time interval at which analysis routines are called.','0.')
CALL prms%CreateIntOption(    'nSkipAnalyze'         , '(Skip Analyze_dt)','1')

CALL prms%CreateRealOption(   'SkipAnalyzeWindow'    , 'Reoccurring time frame when switching between nSkipAnalyze and SkipAnalyzeWindow.','-1.')
CALL prms%CreateRealOption(   'SkipAnalyzeSwitchTime', 'Time within the reoccurring time frame, when using nSkipAnalyzeSwitch instead of nSkipAnalyze.')
CALL prms%CreateIntOption(    'nSkipAnalyzeSwitch'   , 'Skip Analyze_dt with a different values as nSkipAnalyze')
CALL prms%CreateRealOption(   'OutputTimeFixed'      , 'fixed time for writing state to .h5','-1.0')
CALL prms%CreateLogicalOption('DoMeasureAnalyzeTime' , 'Measure time that is spent in analyze routines and count the number of '//&
                                                       'analysis calls (to std out stream)','.FALSE.')
!CALL prms%CreateLogicalOption('AnalyzeToFile',   "Set true to output result of error norms to a file (DoCalcErrorNorms=T)",&
                                                 !'.FALSE.')
!CALL prms%CreateIntOption(    'nWriteData' ,     "Interval as multiple of Analyze_dt at which HDF5 files "//&
                                                 !"(e.g. State,TimeAvg,Fluc) are written.",&
                                                 !'1')
!CALL prms%CreateIntOption(    'AnalyzeExactFunc',"Define exact function used for analyze (e.g. for computing L2 errors). "//&
                                                 !"Default: Same as IniExactFunc")
!CALL prms%CreateIntOption(    'AnalyzeRefState' ,"Define state used for analyze (e.g. for computing L2 errors). "//&
                                                 !"Default: Same as IniRefState")
!CALL DefineParametersAnalyzeEquation()

! -------------------------
CALL prms%SetSection("Analyzefield")
! -------------------------
CALL prms%CreateIntOption(    'Field-AnalyzeStep'   , 'Analyze is performed each Nth time step. Set to 0 to completely skip.','1')

!-- CalcStuff
CALL prms%CreateLogicalOption(  'CalcPotentialEnergy', 'Calculate Potential Energy. Output file is Database.csv','.FALSE.')
CALL prms%CreateLogicalOption(  'CalcPointsPerWavelength', 'Flag to compute the points per wavelength in each cell','.FALSE.')

CALL prms%CreateLogicalOption( 'CalcPoyntingVecIntegral',"Calculate Poynting vector integral, which is the integrated energy density "//&
                                                         "over a plane, perpendicular to PoyntingMainDir axis (default is z-direction)"//&
                                                         ". The plane position must lie on an interface between two adjacent elements.",&
                                                      '.FALSE.')

!-- BoundaryFieldOutput
CALL prms%CreateLogicalOption(  'CalcBoundaryFieldOutput', 'Output the field boundary over time' , '.FALSE.')
CALL prms%CreateIntOption(      'BFO-NFieldBoundaries'   , 'Number of boundaries used for CalcBoundaryFieldOutput')
CALL prms%CreateIntArrayOption( 'BFO-FieldBoundaries'    , 'Vector (length BFO-NFieldBoundaries) with the numbers of each Field-Boundary', no=0)

!-- Poynting Vector
CALL prms%CreateIntOption( 'PoyntingVecInt-Planes', 'Total number of Poynting vector integral planes for measuring the '//&
                                                    'directed power flow (energy flux density: Density and direction of an '//&
                                                    'electromagnetic field.', '0')
CALL prms%CreateRealOption('Plane-Tolerance'  , 'Absolute tolerance for checking the Poynting vector integral plane '//&
                                                'coordinates and normal vectors of the corresponding sides for selecting '//&
                                                'relevant sides', '1E-5')
CALL prms%CreateRealOption('Plane-[$]-x-coord', 'x-coordinate of the n-th Poynting vector plane (when PoyntingMainDir=1)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Plane-[$]-y-coord', 'y-coordinate of the n-th Poynting vector plane (when PoyntingMainDir=2)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Plane-[$]-z-coord', 'z-coordinate of the n-th Poynting vector plane (when PoyntingMainDir=3)', '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption( 'PoyntingMainDir'  , 'Direction in which the Poynting vector integral is to be measured. '//&
                                                   '\n1: x \n2: y \n3: z (default)','3')

#if USE_HDG
!-- AverageElectricPotential
CALL prms%CreateLogicalOption( 'CalcAverageElectricPotential',"Calculate the averaged electric potential at a specific x-coordinate."//&
                                                              " The plane position must lie on an interface between two adjacent elements",&
                                                              '.FALSE.')
CALL prms%CreateRealOption('AvgPotential-Plane-x-coord', 'x-coordinate of the averaged electric potential')
CALL prms%CreateRealOption('AvgPotential-Plane-Tolerance', 'Absolute tolerance for checking the averaged electric potential plane '&
                                                         , '1E-5')
CALL prms%CreateLogicalOption( 'CalcElectricTimeDerivative' ,"Calculate the time derivative of the electric displacement field D=eps*E and output to .h5 and .csv files.",".FALSE.")
#endif /*USE_HDG*/
!-- TimeAverage
CALL prms%CreateLogicalOption( 'CalcTimeAverage'            , 'Flag if time averaging should be performed','.FALSE.')
CALL prms%CreateIntOption(     'nSkipAvg'                   , 'Iter every which CalcTimeAverage is performed')
CALL prms%CreateStringOption(  'VarNameAvg'                 , 'Count of time average variables',multiple=.TRUE.)
CALL prms%CreateStringOption(  'VarNameFluc'                , 'Count of fluctuation variables',multiple=.TRUE.)

!-- Code Analyze
#ifdef CODE_ANALYZE
CALL prms%CreateLogicalOption( 'DoCodeAnalyzeOutput'        , 'print code analyze info to CodeAnalyze.csv','.TRUE.')
#endif /* CODE_ANALYZE */
END SUBROUTINE DefineParametersAnalyze

SUBROUTINE InitAnalyze()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
#if PP_nVar>=6
USE MOD_AnalyzeField          ,ONLY: GetPoyntingIntPlane
USE MOD_Analyze_Vars          ,ONLY: CalcPoyntingInt
#endif /*PP_nVar>=6*/
USE MOD_Analyze_Vars          ,ONLY: AnalyzeInitIsDone,Analyze_dt,DoCalcErrorNorms,OutputErrorNormsToH5
USE MOD_Analyze_Vars          ,ONLY: CalcPointsPerWavelength,PPWCell,OutputTimeFixed,FieldAnalyzeStep
USE MOD_Analyze_Vars          ,ONLY: AnalyzeCount,AnalyzeTime,DoMeasureAnalyzeTime
USE MOD_Analyze_Vars          ,ONLY: doFieldAnalyze,CalcEpot
USE MOD_Analyze_Vars          ,ONLY: CalcBoundaryFieldOutput,BFO
USE MOD_Analyze_Vars          ,ONLY: nSkipAnalyze,SkipAnalyzeWindow,SkipAnalyzeSwitchTime,nSkipAnalyzeSwitch
USE MOD_Interpolation_Vars    ,ONLY: InterpolationInitIsDone,Uex,NAnalyze
USE MOD_IO_HDF5               ,ONLY: AddToElemData,ElementOut
USE MOD_Mesh_Vars             ,ONLY: nElems
USE MOD_ReadInTools           ,ONLY: GETINT,GETREAL,GETLOGICAL,PrintOption,GETINTARRAY
USE MOD_TimeAverage_Vars      ,ONLY: doCalcTimeAverage
USE MOD_TimeAverage           ,ONLY: InitTimeAverage
USE MOD_TimeDisc_Vars         ,ONLY: TEnd
USE MOD_Equation_vars         ,ONLY: Wavelength
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemCharLength_Shared
#if USE_MPI && defined(PARTICLES)
USE MOD_Mesh_Vars             ,ONLY: offSetElem
#endif /*USE_MPI && defined(PARTICLES)*/
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
#if USE_HDG
USE MOD_Analyze_Vars          ,ONLY: CalcAverageElectricPotential,PosAverageElectricPotential,CalcElectricTimeDerivative
USE MOD_AnalyzeField          ,ONLY: GetAverageElectricPotentialPlane
#ifdef PARTICLES
USE MOD_PICInterpolation_Vars ,ONLY: useAlgebraicExternalField,AlgebraicExternalField
#endif /*PARTICLES*/
#endif /*USE_HDG*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars             ,ONLY: BoundaryType,BoundaryName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)   :: DefStr
INTEGER             :: iElem,CNElemID,iBoundary,BCType,BCState,iBC
REAL                :: PPWCellMax,PPWCellMin
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.AnalyzeInitIsDone) THEN
  CALL abort(__STAMP__,'InitAnalyse not ready to be called or already called.')
  RETURN
END IF
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT ANALYZE...'

! Get logical for calculating the error norms L2 and LInf
DoCalcErrorNorms = GETLOGICAL('DoCalcErrorNorms')

IF(DoCalcErrorNorms)THEN
  ! Get logical for writing the analytical solution, the error norms L2 and LInf to .h5
  OutputErrorNormsToH5 = GETLOGICAL('OutputErrorNormsToH5')
  ! Allocate container for exact solution (Gauss-Lobatto nodes)
  ALLOCATE(Uex(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze,1:nElems))
  Uex = 0.
END IF ! DoCalcErrorNorms

! Get the time step for performing analyzes and integer for skipping certain steps
WRITE(DefStr,WRITEFORMAT) TEnd
Analyze_dt        = GETREAL('Analyze_dt',DefStr)
nSkipAnalyze      = GETINT('nSkipAnalyze')

! Get 2nd option for skipping certain steps (within reoccurring time frame SkipAnalyzeWindow when time > SkipAnalyzeSwitchTime)
SkipAnalyzeWindow = GETREAL('SkipAnalyzeWindow')
IF(SkipAnalyzeWindow.GT.0.)THEN
  SkipAnalyzeSwitchTime = GETREAL('SkipAnalyzeSwitchTime')
  nSkipAnalyzeSwitch    = GETINT('nSkipAnalyzeSwitch')
  IF(SkipAnalyzeSwitchTime.GT.SkipAnalyzeWindow) CALL abort(__STAMP__,'SkipAnalyzeSwitchTime must be smaller than SkipAnalyzeWindow')
ELSE
  SkipAnalyzeSwitchTime = HUGE(1.)
  nSkipAnalyzeSwitch    = nSkipAnalyze
END IF ! SkipAnalyzeWindow.GT.0.

OutputTimeFixed   = GETREAL('OutputTimeFixed')
! Time averaged quantises fields (Maxwell/Poisson solver) and deposited particles (PIC)
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6)||(PP_TimeDiscMethod>=501 && PP_TimeDiscMethod<=509)
  doCalcTimeAverage = GETLOGICAL('CalcTimeAverage')
#else
  doCalcTimeAverage = .FALSE.
#endif

IF(doCalcTimeAverage)  CALL InitTimeAverage()

FieldAnalyzeStep  = GETINT('Field-AnalyzeStep')
DoFieldAnalyze    = .FALSE.
#if (USE_HDG)
IF (FieldAnalyzeStep.GT.0) DoFieldAnalyze = .TRUE.
#endif /*USE_HDG*/
IF (FieldAnalyzeStep.EQ.0) FieldAnalyzeStep = HUGE(FieldAnalyzeStep)
CalcEpot          = GETLOGICAL('CalcPotentialEnergy')
IF(CalcEpot)        DoFieldAnalyze = .TRUE.
#if PP_nVar>=6
CalcPoyntingInt = GETLOGICAL('CalcPoyntingVecIntegral') ! PoyntingVecIntegral
IF(CalcPoyntingInt)THEN
  DoFieldAnalyze = .TRUE.
  CALL GetPoyntingIntPlane()
END IF
#endif /*PP_nVar>=6*/

#if USE_HDG
!--- Get averaged electric potential
CalcAverageElectricPotential=.FALSE. ! Initialize
#ifdef PARTICLES
IF(useAlgebraicExternalField.AND.AlgebraicExternalField.EQ.1)THEN
  CalcAverageElectricPotential = .TRUE.
  PosAverageElectricPotential  = 2.4e-2 ! automatic setting for this case
  CALL PrintOption('AlgebraicExternalField.EQ.1, Setting CalcAverageElectricPotential','INFO',LogOpt=CalcAverageElectricPotential)
  CALL PrintOption('AlgebraicExternalField.EQ.1, Setting PosAverageElectricPotential ','INFO',RealOpt=PosAverageElectricPotential)
ELSE
#endif /*PARTICLES*/
  CalcAverageElectricPotential = GETLOGICAL('CalcAverageElectricPotential') ! user-defined activation
  IF(CalcAverageElectricPotential) PosAverageElectricPotential = GETREAL('AvgPotential-Plane-x-coord')    ! user-required input
#ifdef PARTICLES
END IF
#endif /*PARTICLES*/
IF(CalcAverageElectricPotential)THEN
  DoFieldAnalyze = .TRUE.
  CALL GetAverageElectricPotentialPlane()
END IF

!-- Electric displacement current
! Calculate the time derivative of D=eps0*E and output to h5
CalcElectricTimeDerivative = GETLOGICAL('CalcElectricTimeDerivative')
CALL InitCalcElectricTimeDerivativeSurface()

#endif /*USE_HDG*/

!-- BoundaryParticleOutput (after mapping of PartBound on FieldBound and determination of PartBound types = open, reflective etc.)
CalcBoundaryFieldOutput = GETLOGICAL('CalcBoundaryFieldOutput')
IF(CalcBoundaryFieldOutput)THEN
  DoFieldAnalyze = .TRUE.
  BFO%NFieldBoundaries = GETINT('BFO-NFieldBoundaries')
  BFO%FieldBoundaries  = GETINTARRAY('BFO-FieldBoundaries',BFO%NFieldBoundaries)
  DO iBoundary=1,BFO%NFieldBoundaries
    iBC = BFO%FieldBoundaries(iBoundary)
    IF(iBC.GT.SIZE(BoundaryName)) CALL abort(__STAMP__,'BFO-FieldBoundaries BC index is too large: ',IntInfoOpt=iBC)
    BCType  = BoundaryType(iBC,BC_TYPE)
    BCState = BoundaryType(iBC,BC_STATE)
    LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A)')&
       ' Activated BFO of electric potential for ['//TRIM(BoundaryName(iBC))//'] with BCType [',BCType,'] and BCState [',BCState,']'
  END DO
END IF ! CalcBoundaryFieldOutput

! Get logical for measurement of time spent in analyze routines
DoMeasureAnalyzeTime = GETLOGICAL('DoMeasureAnalyzeTime')
! Initialize time and counter for analyze measurement
AnalyzeCount = 0
AnalyzeTime  = 0.0

AnalyzeInitIsDone = .TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT ANALYZE DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

! Points Per Wavelength
CalcPointsPerWavelength = GETLOGICAL('CalcPointsPerWavelength')
IF(CalcPointsPerWavelength)THEN
  ! calculate cell local number excluding neighbor DOFs
  ALLOCATE( PPWCell(1:PP_nElems) )
  PPWCell=0.0
  CALL AddToElemData(ElementOut,'PPWCell',RealArray=PPWCell(1:PP_nElems))
  ! Calculate PPW for each cell
  IF(WaveLength.LT.0.) WaveLength = GETREAL('WaveLength','1.')
  CALL PrintOption('Wavelength for PPWCell','OUTPUT',RealOpt=Wavelength)
  PPWCellMin=HUGE(1.)
  PPWCellMax=-HUGE(1.)
  DO iElem = 1, nElems
    ! In case of MPI=ON and PARTICLES=OFF, no shared array is created and all arrays are processor-local
#if USE_MPI && defined(PARTICLES)
    CNElemID = GetCNElemID(iElem+offSetElem)
#else
    CNElemID = iElem
#endif /*USE_MPI && defined(PARTICLES)*/
    PPWCell(iElem) = (REAL(PP_N)+1.)*Wavelength/ElemCharLength_Shared(CNElemID)
    PPWCellMin     = MIN(PPWCellMin,PPWCell(iElem))
    PPWCellMax     = MAX(PPWCellMax,PPWCell(iElem))
  END DO ! iElem = 1, nElems
#if USE_MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , PPWCellMin , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , PPWCellMax , 1 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(PPWCellMin   , 0          , 1 , MPI_DOUBLE_PRECISION , MPI_MIN , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(PPWCellMax   , 0          , 1 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#endif /*USE_MPI*/
  CALL PrintOption('MIN(PPWCell)','CALCUL.',RealOpt=PPWCellMin)
  CALL PrintOption('MAX(PPWCell)','CALCUL.',RealOpt=PPWCellMax)
END IF
END SUBROUTINE InitAnalyze


#if USE_HDG
SUBROUTINE CalcError(L_2_Error,L_Inf_Error)
#else
SUBROUTINE CalcError(time,L_2_Error,L_Inf_Error)
#endif
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_DG_Vars            ,ONLY: U
USE MOD_Equation           ,ONLY: ExactFunc
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Interpolation_Vars ,ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wAnalyze,Uex
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP,sJ
USE MOD_Particle_Mesh_Vars ,ONLY: MeshVolume
USE MOD_Analyze_Vars       ,ONLY: OutputErrorNormsToH5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#if !(USE_HDG)
REAL,INTENT(IN)               :: time
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: L_2_Error(PP_nVar)   !< L2 error of the solution
REAL,INTENT(OUT)              :: L_Inf_Error(PP_nVar) !< LInf error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,k,l,m
REAL                          :: U_exact(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                          :: IntegrationWeight
!===================================================================================================================================
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,PP_nElems
  ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
  CALL ChangeBasis3D(3,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,Elem_xGP(1:3,:,:,:,iElem),Coords_NAnalyze(1:3,:,:,:))
  ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the jacobian ;-)
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
  ! Interpolate the solution to the analyze grid
  CALL ChangeBasis3D(PP_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,U(1:PP_nVar,:,:,:,iElem),U_NAnalyze(1:PP_nVar,:,:,:))
  DO m=0,NAnalyze
    DO l=0,NAnalyze
      DO k=0,NAnalyze
#if USE_HDG
        CALL ExactFunc(IniExactFunc,Coords_NAnalyze(1:3,k,l,m),U_exact(1:PP_nVar,k,l,m),ElemID=iElem)
#else
        CALL ExactFunc(IniExactFunc,time,0,Coords_NAnalyze(1:3,k,l,m),U_exact(1:PP_nVar,k,l,m))
#endif
        L_Inf_Error = MAX(L_Inf_Error,abs(U_NAnalyze(:,k,l,m) - U_exact(1:PP_nVar,k,l,m)))
        IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
        ! To sum over the elements, We compute here the square of the L_2 error
        L_2_Error = L_2_Error+(U_NAnalyze(:,k,l,m) - U_exact(1:PP_nVar,k,l,m))*&
                              (U_NAnalyze(:,k,l,m) - U_exact(1:PP_nVar,k,l,m))*IntegrationWeight
      END DO ! k
    END DO ! l
  END DO ! m
  ! Output the exact solution, the L2 error and LInf error to .h5 (in NodeTypeGL = 'GAUSS-LOBATTO')
  IF(OutputErrorNormsToH5)THEN
    Uex(1:PP_nVar,:,:,:,iElem) = U_exact(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
  END IF ! OutputErrorNormsToH5
END DO ! iElem=1,PP_nElems
#if USE_MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error   , PP_nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error , PP_nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error   , 0            , PP_nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(L_Inf_Error , 0            , PP_nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/MeshVolume)

END SUBROUTINE CalcError


#ifdef PARTICLES
SUBROUTINE CalcErrorPartSource(PartSource_nVar,L_2_PartSource,L_Inf_PartSource)
!===================================================================================================================================
! Calculates the L2 (particle) source error by integrating the difference between the numerical and analytical solution
! over all elements.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Equation           ,ONLY: ExactFunc
USE MOD_Interpolation_Vars ,ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wAnalyze
USE MOD_Mesh_Vars          ,ONLY: sJ
USE MOD_PICDepo_Vars       ,ONLY: PartSourceOld
USE MOD_Particle_Mesh_Vars ,ONLY: MeshVolume
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: PartSource_nVar
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: L_2_PartSource(PartSource_nVar)   !< L2 error of the source
REAL,INTENT(OUT)              :: L_Inf_PartSource(PartSource_nVar) !< LInf error of the source
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,k,l,m
REAL                          :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                          :: IntegrationWeight

REAL                          :: PartSource_NAnalyze1(1:PartSource_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: PartSource_NAnalyze2(1:PartSource_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
!===================================================================================================================================

! nullify
L_Inf_PartSource(:)=-1.E10
L_2_PartSource(:)=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,PP_nElems
  ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  CALL ChangeBasis3D(1,PP_N,NAnalyze,Vdm_GaussN_NAnalyze,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_NAnalyze(1:1,:,:,:))
  CALL ChangeBasis3D(PartSource_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze &
      ,PartSourceOld(1:PartSource_nVar,1,:,:,:,iElem),PartSource_NAnalyze1(1:PartSource_nVar,:,:,:))
  CALL ChangeBasis3D(PartSource_nVar,PP_N,NAnalyze,Vdm_GaussN_NAnalyze &
      ,PartSourceOld(1:PartSource_nVar,2,:,:,:,iElem),PartSource_NAnalyze2(1:PartSource_nVar,:,:,:))
  PartSourceOld(1:PartSource_nVar,2,:,:,:,iElem)=PartSourceOld(1:PartSource_nVar,1,:,:,:,iElem)
  DO m=0,NAnalyze
    DO l=0,NAnalyze
      DO k=0,NAnalyze
        IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
        L_Inf_PartSource = MAX(L_Inf_PartSource,abs(PartSource_NAnalyze1(:,k,l,m) - PartSource_NAnalyze2(:,k,l,m)))
        L_2_PartSource = L_2_PartSource+(PartSource_NAnalyze1(:,k,l,m) - PartSource_NAnalyze2(:,k,l,m)) &
            *(PartSource_NAnalyze1(:,k,l,m) - PartSource_NAnalyze2(:,k,l,m))*IntegrationWeight
      END DO ! k
    END DO ! l
  END DO ! m
END DO ! iElem=1,PP_nElems
#if USE_MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE  , L_2_PartSource   , PartSource_nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(MPI_IN_PLACE  , L_Inf_PartSource , PartSource_nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
ELSE
  CALL MPI_REDUCE(L_2_PartSource   , 0             , PartSource_nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(L_Inf_PartSource , 0             , PartSource_nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ! in this case the receive value is not relevant.
END IF
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_PartSource = SQRT(L_2_PartSource/MeshVolume)

END SUBROUTINE CalcErrorPartSource
#endif /*PARTICLES*/


SUBROUTINE CalcErrorStateFiles(nVar,N1,N2,U1,U2)
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: sJ
USE MOD_Interpolation_Vars ,ONLY: NAnalyze,wAnalyze
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Basis              ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
INTEGER,INTENT(IN)           :: nVar
INTEGER,INTENT(IN)           :: N1
INTEGER,INTENT(IN)           :: N2
REAL,INTENT(IN)              :: U1(1:nVar,0:N1,0:N1,0:N1,nElems)
REAL,INTENT(IN)              :: U2(1:nVar,0:N2,0:N2,0:N2,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                          :: L_2_Error(nVar)   !< L2 error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE  :: Vdm_GaussN_NAnalyze1(:,:)    ! for interpolation to Analyze points
REAL,ALLOCATABLE  :: Vdm_GaussN_NAnalyze2(:,:)    ! for interpolation to Analyze points
INTEGER                       :: iElem,k,l,m
REAL                          :: L_Inf_Error(nVar),L_2_Error2(nVar),L_Inf_Error2(nVar)
REAL                          :: U1_NAnalyze(1:nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: U2_NAnalyze(1:nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                          :: Volume,IntegrationWeight
#if USE_MPI
REAL                          :: Volume2
#endif
CHARACTER(LEN=40)             :: formatStr
REAL                          :: xGP1(0:N1),xGP2(0:N2),wGP1(0:N1),wGP2(0:N2),wBary1(0:N1),wBary2(0:N2)
REAL ,DIMENSION(0:NAnalyze)   :: XiAnalyze
!===================================================================================================================================
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.
L_Inf_Error2(:)=-1.E10
L_2_Error2(:)=0.
Volume=0.

CALL LegendreGaussNodesAndWeights(N1,xGP1,wGP1)
CALL LegendreGaussNodesAndWeights(N2,xGP2,wGP2)

CALL BarycentricWeights(N1,xGP1,wBary1)
CALL BarycentricWeights(N2,xGP2,wBary2)

ALLOCATE(wAnalyze(0:NAnalyze))
ALLOCATE(Vdm_GaussN_NAnalyze1(0:NAnalyze,0:N1))
ALLOCATE(Vdm_GaussN_NAnalyze2(0:NAnalyze,0:N2))

CALL LegGaussLobNodesAndWeights(NAnalyze,XiAnalyze,wAnalyze)

CALL InitializeVandermonde(N1,NAnalyze,wBary1,xGP1,XiAnalyze,Vdm_GaussN_NAnalyze1)
CALL InitializeVandermonde(N2,NAnalyze,wBary2,xGP2,XiAnalyze,Vdm_GaussN_NAnalyze2)



! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
   ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the Jacobian ;-)
   J_N(1,0:N1,0:N1,0:N1)=1./sJ(:,:,:,iElem)
   CALL ChangeBasis3D(1,N1,NAnalyze,Vdm_GaussN_NAnalyze1,J_N(1:1,0:N1,0:N1,0:N1),J_NAnalyze(1:1,:,:,:))



   ! Interpolate the solution to the analyze grid
   CALL ChangeBasis3D(nVar,N1,NAnalyze,Vdm_GaussN_NAnalyze1,U1(1:nVar,:,:,:,iElem),U1_NAnalyze(1:nVar,:,:,:))
   CALL ChangeBasis3D(nVar,N2,NAnalyze,Vdm_GaussN_NAnalyze2,U2(1:nVar,:,:,:,iElem),U2_NAnalyze(1:nVar,:,:,:))

   DO m=0,NAnalyze
     DO l=0,NAnalyze
       DO k=0,NAnalyze
         L_Inf_Error = MAX(L_Inf_Error,abs(U1_NAnalyze(:,k,l,m) - U2_NAnalyze(:,k,l,m)))
         IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
         ! To sum over the elements, We compute here the square of the L_2 error
         L_2_Error = L_2_Error+ (U1_NAnalyze(:,k,l,m) - U2_NAnalyze(:,k,l,m))*&
                                (U1_NAnalyze(:,k,l,m) - U2_NAnalyze(:,k,l,m))*IntegrationWeight
         Volume = Volume + IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO ! iElem=1,nElems
#if USE_MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error   , nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , volume      , 1    , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error , nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error   , L_2_Error2   , nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(volume      , volume2      , 1    , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(L_Inf_Error , L_Inf_Error2 , nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/Volume)

! Graphical output
IF(MPIroot) THEN
  WRITE(formatStr,'(A5,I1,A7)')'(A13,',nVar,'ES16.7)'
  WRITE(UNIT_StdOut,formatStr)' L_2       : ',L_2_Error
  WRITE(UNIT_StdOut,formatStr)' L_inf     : ',L_Inf_Error
END IF
END SUBROUTINE CalcErrorStateFiles


SUBROUTINE CalcErrorStateFileSigma(nVar,N1,U1)
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of state variables using the Analyze Framework (GL points+weights)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: sJ
USE MOD_Interpolation_Vars ,ONLY: NAnalyze,wAnalyze
USE MOD_Equation           ,ONLY: ExactFunc
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Basis              ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: nVar
INTEGER,INTENT(IN)           :: N1
REAL,INTENT(IN)              :: U1(1:nVar,0:N1,0:N1,0:N1,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                          :: L_2_Error(nVar)   !< L2 error of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE              :: Vdm_GaussN_NAnalyze1(:,:)    ! for interpolation to Analyze points
INTEGER                       :: iElem,k,l,m
REAL                          :: L_Inf_Error(nVar),L_2_Error2(nVar),L_Inf_Error2(nVar)
REAL                          :: U1_NAnalyze(1:nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL                          :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                          :: Volume,IntegrationWeight
#if USE_MPI
REAL                          :: Volume2
#endif
CHARACTER(LEN=40)             :: formatStr
REAL                          :: xGP1(0:N1),wGP1(0:N1),wBary1(0:N1)
REAL ,DIMENSION(0:NAnalyze)   :: XiAnalyze
!===================================================================================================================================
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.
L_Inf_Error2(:)=-1.E10
L_2_Error2(:)=0.
Volume=0.

CALL LegendreGaussNodesAndWeights(N1,xGP1,wGP1)

CALL BarycentricWeights(N1,xGP1,wBary1)

ALLOCATE(Vdm_GaussN_NAnalyze1(0:NAnalyze,0:N1))

CALL LegGaussLobNodesAndWeights(NAnalyze,XiAnalyze,wAnalyze)

CALL InitializeVandermonde(N1,NAnalyze,wBary1,xGP1,XiAnalyze,Vdm_GaussN_NAnalyze1)



! Interpolate values of Error-Grid from GP's
DO iElem=1,nElems
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,0:N1,0:N1,0:N1)=1./sJ(:,:,:,iElem)
   CALL ChangeBasis3D(1,N1,NAnalyze,Vdm_GaussN_NAnalyze1,J_N(1:1,0:N1,0:N1,0:N1),J_NAnalyze(1:1,:,:,:))



   ! Interpolate the solution to the analyze grid
   CALL ChangeBasis3D(nVar,N1,NAnalyze,Vdm_GaussN_NAnalyze1,U1(1:nVar,:,:,:,iElem),U1_NAnalyze(1:nVar,:,:,:))

   DO m=0,NAnalyze
     DO l=0,NAnalyze
       DO k=0,NAnalyze
         L_Inf_Error = MAX(L_Inf_Error,sqrt(abs(U1_NAnalyze(:,k,l,m))))
         IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
         ! To sum over the elements, We compute here the square of the L_2 error (sigma in state is already squared!!!)
         L_2_Error = L_2_Error+ (U1_NAnalyze(:,k,l,m))*IntegrationWeight
         Volume = Volume + IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO ! iElem=1,nElems
#if USE_MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error    , nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , volume       , 1    , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error  , nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error    , L_2_Error2   , nVar , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(volume       , volume2      , 1    , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(L_Inf_Error  , L_Inf_Error2 , nVar , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/Volume)

! Graphical output
IF(MPIroot) THEN
  WRITE(formatStr,'(A5,I1,A7)')'(A13,',nVar,'ES16.7)'
  WRITE(UNIT_StdOut,formatStr)' L2_sigma  : ',L_2_Error
  WRITE(UNIT_StdOut,formatStr)' Linf_sigma: ',L_Inf_Error
END IF
END SUBROUTINE CalcErrorStateFileSigma


SUBROUTINE AnalyzeToFile(time,CalcTime,L_2_Error)
!===================================================================================================================================
! Writes the L2-error norms to file.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_TimeDisc_Vars ,ONLY:iter
USE MOD_Globals_Vars  ,ONLY:ProjectName
USE MOD_Mesh_Vars    ,ONLY:nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: time                         ! physical time
REAL,INTENT(IN)                :: CalcTime                     ! computational time
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
!LOGICAL                        :: fileExists                   ! Error handler for file
INTEGER                        :: ioUnit
!===================================================================================================================================
Filename = 'out.'//TRIM(ProjectName)//'.dat'
! Check for file
! INQUIRE(FILE = Filename, EXIST = fileExists) ! now -> FILEEXISTS(Filename)
! FILEEXISTS(Filename)
!! File processing starts here open old and extract information or create new file.
ioUnit=1746 ! This number must be fixed?
  OPEN(UNIT   = ioUnit       ,&
       FILE   = Filename     ,&
       STATUS = 'Unknown'    ,&
       ACCESS = 'SEQUENTIAL' ,&
       IOSTAT = openStat                 )
  IF (openStat.NE.0) THEN
     WRITE(*,*)'ERROR: cannot open Outfile'
  END IF
  ! Create a new file with the Tecplot (ASCII file, not binary) header etc.
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
USE MOD_Globals
#if PP_nVar>=6
USE MOD_Analyze_Vars       ,ONLY: CalcPoyntingInt
USE MOD_AnalyzeField       ,ONLY: FinalizePoyntingInt
#endif /*PP_nVar>=6*/
USE MOD_Analyze_Vars       ,ONLY: PPWCell,AnalyzeInitIsDone
USE MOD_TimeAverage_Vars   ,ONLY: doCalcTimeAverage
USE MOD_TimeAverage        ,ONLY: FinalizeTimeAverage
#if USE_HDG
USE MOD_Analyze_Vars       ,ONLY: CalcAverageElectricPotential,EDC
USE MOD_AnalyzeField       ,ONLY: FinalizeAverageElectricPotential
USE MOD_Analyze_Vars       ,ONLY: CalcElectricTimeDerivative
#endif /*USE_HDG*/
USE MOD_Interpolation_Vars ,ONLY: UEx
! IMPLICIT VARIABLE HANDLINGDG
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_HDG && USE_MPI
INTEGER :: iEDCBC
#endif /*USE_HDG && USE_MPI*/
!===================================================================================================================================
#if PP_nVar>=6
IF(CalcPoyntingInt) CALL FinalizePoyntingInt()
#endif /*PP_nVar>=6*/
#if USE_HDG
IF(CalcAverageElectricPotential) CALL FinalizeAverageElectricPotential()
! Electric displacement current
IF(CalcElectricTimeDerivative)THEN
  SDEALLOCATE(EDC%Current)
  SDEALLOCATE(EDC%FieldBoundaries)
  SDEALLOCATE(EDC%BCIDToEDCBCID)
#if USE_MPI
  DO iEDCBC = 1, EDC%NBoundaries
    IF(EDC%COMM(iEDCBC)%UNICATOR.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(EDC%COMM(iEDCBC)%UNICATOR,iERROR)
  END DO
  SDEALLOCATE(EDC%COMM)
#endif /*USE_MPI*/
END IF ! CalcElectricTimeDerivative
#endif /*USE_HDG*/
IF(doCalcTimeAverage) CALL FinalizeTimeAverage
SDEALLOCATE(PPWCell)
SDEALLOCATE(UEx)
AnalyzeInitIsDone = .FALSE.
END SUBROUTINE FinalizeAnalyze


SUBROUTINE PerformAnalyze(OutputTime,FirstOrLastIter,OutPutHDF5)
!===================================================================================================================================
! Check if the analyze subroutines are called depending on the input parameters
! Input parameters
! 1) OutputTime
!    * current time of analyze
! 2) FirstOrLastIter
!    * logical flag for first or last iteration
!      This step is required for the correct opening and closing of a *.csv Database. Furthermore it is needed to prevent
!      duplicates from the *.csv Database file
! 3) OutPutHDF5
!    * OutputHDF5 is true if a state file is written
! The perform-analyze routine is called four times within the timedisc
! 1) initialize before the first iteration. call is performed for an initial computation and a restart
! 2) after the time update
! 3) during an analyze step and writing of a state file
! 4) during a load-balance step
! This routine calls all other analyze-subroutines, which write data to a csv file
! Currently these are:
! I)    AnalyzeField             ->  FieldAnalyze.csv
! II)   AnalyzeParticles         ->  PartAnalyze.csv
! III)  AnalyzeSurface           ->  SurfaceDatabase.csv
! IV)   TrackParticlePosition    ->  ParticlePosition.csv
! V)    AnalyticParticleMovement ->  ParticlePositionAnalytic.csv
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars              ,ONLY: DoCalcErrorNorms,OutputErrorNorms,FieldAnalyzeStep
USE MOD_Analyze_Vars              ,ONLY: AnalyzeCount,AnalyzeTime,DoMeasureAnalyzeTime
USE MOD_Restart_Vars              ,ONLY: DoRestart
USE MOD_TimeDisc_Vars             ,ONLY: iter,tEnd
USE MOD_RecordPoints              ,ONLY: RecordPoints
USE MOD_LoadDistribution          ,ONLY: WriteElemTimeStatistics
USE MOD_Globals_Vars              ,ONLY: ProjectName
USE MOD_AnalyzeField              ,ONLY: AnalyzeField
#ifdef PARTICLES
USE MOD_Mesh_Vars                 ,ONLY: MeshFile
USE MOD_Particle_Vars             ,ONLY: WriteMacroVolumeValues,WriteMacroSurfaceValues,MacroValSamplIterNum,ExcitationSampleData
USE MOD_Particle_Vars             ,ONLY: SampleElecExcitation
USE MOD_Particle_Analyze          ,ONLY: AnalyzeParticles
USE MOD_Particle_Analyze_Tools    ,ONLY: CalculatePartElemData
USE MOD_Particle_Analyze_Output   ,ONLY: WriteParticleTrackingData
USE MOD_Particle_Analyze_Vars     ,ONLY: PartAnalyzeStep,DoPartAnalyze,TrackParticlePosition
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfaceAnalyzeStep
USE MOD_SurfaceModel_Analyze      ,ONLY: AnalyzeSurface
USE MOD_DSMC_Vars                 ,ONLY: DSMC, iter_macvalout,iter_macsurfvalout
USE MOD_DSMC_Vars                 ,ONLY: DSMC_Solution
USE MOD_Particle_Tracking_vars    ,ONLY: ntracks,tTracking,tLocalization,MeasureTrackTime
USE MOD_BGK_Vars                  ,ONLY: BGKInitDone, BGK_QualityFacSamp
USE MOD_FPFlow_Vars               ,ONLY: FPInitDone, FP_QualityFacSamp
USE MOD_DSMC_Vars                 ,ONLY: useDSMC
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars    ,ONLY: nComputeNodeSurfTotalSides, CalcSurfaceImpact
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallState,SampWallImpactEnergy,SampWallImpactVector
USE MOD_Particle_Boundary_Vars    ,ONLY: SampWallPumpCapacity,SampWallImpactAngle,SampWallImpactNumber
USE MOD_DSMC_Analyze              ,ONLY: DSMC_data_sampling, WriteDSMCToHDF5
USE MOD_Particle_Boundary_Sampling,ONLY: CalcSurfaceValues
USE MOD_Particle_Vars             ,ONLY: DelayTime
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars    ,ONLY: rTotalBBChecks,rTotalBezierClips,rTotalBezierNewton
!USE MOD_Particle_Surfaces_Vars    ,ONLY: SideBoundingBoxVolume
USE MOD_Particle_Analyze_Code     ,ONLY: AnalyticParticleMovement
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackingMethod
USE MOD_PICInterpolation_Vars     ,ONLY: DoInterpolationAnalytic
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/
#if (PP_nVar>=6)
USE MOD_AnalyzeField              ,ONLY: CalcPoyntingIntegral
#endif /*PP_nVar>=6*/
#if defined(LSERK) || defined(IMPA) || defined(ROS) || USE_HDG
USE MOD_Analyze_Vars              ,ONLY: DoFieldAnalyze
USE MOD_RecordPoints_Vars         ,ONLY: RP_onProc
#endif /*defined(LSERK) ||  defined(IMPA) || defined(ROS) || USE_HDG*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers        ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#ifdef PARTICLES
USE MOD_PICDepo_Vars              ,ONLY: DoDeposition, RelaxDeposition
#endif /*PARTICLES*/
USE MOD_TimeDisc_Vars             ,ONLY: time
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: OutputTime
LOGICAL,INTENT(IN)            :: FirstOrLastIter
LOGICAL,INTENT(IN)            :: OutputHDF5
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef maxwell
LOGICAL                       :: ProlongToFaceNeeded
#endif /*Maxwell*/
#ifdef PARTICLES
INTEGER                       :: iSide
#if USE_MPI
INTEGER                       :: RECI
REAL                          :: RECR
#endif /*USE_MPI*/
!#ifdef CODE_ANALYZE
!REAL                          :: TotalSideBoundingBoxVolume
!#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/
LOGICAL                       :: LastIter
REAL                          :: L_2_Error(PP_nVar)
REAL                          :: L_Inf_Error(PP_nVar)
#if defined(LSERK) || defined(IMPA) || defined(ROS) || USE_HDG
#if USE_LOADBALANCE
REAL                          :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
#endif /* LSERK && IMPA && ROS && USE_HDG*/
REAL                          :: StartAnalyzeTime,EndAnalyzeTime
CHARACTER(LEN=40)             :: formatStr
LOGICAL                       :: DoPerformFieldAnalyze
LOGICAL                       :: DoPerformPartAnalyze
LOGICAL                       :: DoPerformSurfaceAnalyze
LOGICAL                       :: DoPerformErrorCalc
#ifdef PARTICLES
#if ((USE_HDG) && (PP_nVar==1))
INTEGER                       :: PartSource_nVar=1
REAL                          :: L_2_PartSource(1:1)
REAL                          :: L_Inf_PartSource(1:1)
#else
INTEGER                       :: PartSource_nVar=4
REAL                          :: L_2_PartSource(1:4)
REAL                          :: L_Inf_PartSource(1:4)
#endif
#endif /* PARTICLES */
REAL                          :: CurrentTime
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(6))
#endif /*EXTRAE*/
EndAnalyzeTime = -1.0 ! initialize

! Create .csv file for performance analysis and load balance: write header line
CALL WriteElemTimeStatistics(WriteHeader=.TRUE.,iter_opt=iter,time_opt=time)

! check if final/last iteration iteration
LastIter=.FALSE.
IF(OutputHDF5 .AND. FirstOrLastIter) LastIter=.TRUE.

!----------------------------------------------------------------------------------------------------------------------------------
! Determine if an analyze step has to be performed
! selection is identical with/without particles
!----------------------------------------------------------------------------------------------------------------------------------
! The initial selection is performed for DoPerformAnalyze and copied to the other analyzes.
! The iteration dependent steps are performed later

! Prolongtoface in CalcPoyntingIntegral is always needed, because analyze is not any more performed within the first
! RK stage
#ifdef maxwell
ProlongToFaceNeeded=.TRUE.
#endif /*maxwell*/

! 1)
! Initial start of computation before the first iteration
! this check is identical for all time integration methods
! analyze routines are not called for a restart
! PO: not sure if it this check is any longer needed
IF(FirstOrLastIter)THEN
  DoPerformFieldAnalyze   =.TRUE.
  DoPerformPartAnalyze    =.TRUE.
  DoPerformSurfaceAnalyze =.TRUE.
ELSE
  ! nullify
  DoPerformFieldAnalyze   =.FALSE.
  DoPerformPartAnalyze    =.FALSE.
  DoPerformSurfaceAnalyze =.FALSE.
END IF

! 2) check analyze with respect to iteration counter
! selection criterion depending on iteration counter
! * full stage Runge-Kutta schemes
! * LSERK schemes (because analyze is not hidden in RK stage)
! * DSMC

! FieldAnalyzeStep
! 2) normal analyze at analyze step
IF(MOD(iter,FieldAnalyzeStep).EQ.0 .AND. .NOT. OutPutHDF5) DoPerformFieldAnalyze=.TRUE.
! 3) + 4) force analyze during a write-state information and prevent duplicates
IF(MOD(iter,FieldAnalyzeStep).NE.0 .AND. OutPutHDF5)       DoPerformFieldAnalyze=.TRUE.
! Remove analyze during restart or load-balance step
IF(DoRestart .AND. iter.EQ.0) DoPerformFieldAnalyze=.FALSE.
! BUT print analyze info for CODE_ANALYZE and USE_HDG to get the HDG solver statistics (number of iterations, runtime and norm)
#if USE_HDG
#ifdef CODE_ANALYZE
IF(DoRestart .AND. iter.EQ.0) DoPerformFieldAnalyze=.TRUE.
#endif /*CODE_ANALYZE*/
#endif /*USE_HDG*/
! Finally, remove duplicates for last iteration
! This step is needed, because PerformAnalyze is called twice within the iterations
IF(FirstOrLastIter .AND. .NOT.OutPutHDF5 .AND.iter.NE.0) DoPerformFieldAnalyze=.FALSE.

#ifdef PARTICLES
! PartAnalyzeStep
! 2) normal analyze at analyze step
IF(MOD(iter,PartAnalyzeStep).EQ.0 .AND. .NOT. OutPutHDF5) DoPerformPartAnalyze=.TRUE.
! 3) + 4) force analyze during a write-state information and prevent duplicates
IF(MOD(iter,PartAnalyzeStep).NE.0 .AND. OutPutHDF5)       DoPerformPartAnalyze=.TRUE.
! Remove analyze during restart or load-balance step
IF(DoRestart .AND. iter.EQ.0) DoPerformPartAnalyze=.FALSE.
! Finally, remove duplicates for last iteration
! This step is needed, because PerformAnalyze is called twice within the iterations
IF(FirstOrLastIter.AND.(.NOT.OutPutHDF5).AND.(iter.NE.0)) DoPerformPartAnalyze=.FALSE.

! SurfaceAnalyzeStep
! 2) normal analyze at analyze step
IF(MOD(iter,SurfaceAnalyzeStep).EQ.0 .AND. .NOT. OutPutHDF5) DoPerformSurfaceAnalyze=.TRUE.
! 3) + 4) force analyze during a write-state information and prevent duplicates
IF(MOD(iter,SurfaceAnalyzeStep).NE.0 .AND. OutPutHDF5)       DoPerformSurfaceAnalyze=.TRUE.
! Remove analyze during restart or load-balance step
IF(DoRestart .AND. iter.EQ.0) DoPerformSurfaceAnalyze=.FALSE.
! Finally, remove duplicates for last iteration
! This step is needed, because PerformAnalyze is called twice within the iterations
IF(FirstOrLastIter .AND. .NOT.OutPutHDF5 .AND.iter.NE.0) DoPerformSurfaceAnalyze=.FALSE.
#endif /*PARTICLES*/

! selection of error calculation for first iteration, output time but not last iteration
DoPerformErrorCalc=.FALSE.
IF((FirstOrLastIter.OR.OutputHDF5).AND.(.NOT.LastIter)) DoPerformErrorCalc=.TRUE.

IF((DoPerformFieldAnalyze.OR.DoPerformPartAnalyze.OR.DoPerformSurfaceAnalyze).AND.DoMeasureAnalyzeTime)THEN
  StartAnalyzeTime = PICLASTIME()
  AnalyzeCount     = AnalyzeCount + 1
END IF


!----------------------------------------------------------------------------------------------------------------------------------
! DG-Solver
!----------------------------------------------------------------------------------------------------------------------------------
! Calculate error norms
! This computes the error analysis during each dt_Analysis step
IF(DoCalcErrorNorms) THEN
  IF(DoPerformErrorCalc)THEN
    OutputErrorNorms=.TRUE.
#if USE_HDG
    CALL CalcError(L_2_Error,L_Inf_Error)
#else
    CALL CalcError(OutputTime,L_2_Error,L_Inf_Error)
#endif
    IF (OutputTime.GE.tEnd)THEN
      CurrentTime=PICLASTIME()
      CALL AnalyzeToFile(OutputTime,CurrentTime,L_2_Error)
    END IF
#ifdef PARTICLES
    IF (DoDeposition.AND.RelaxDeposition) CALL CalcErrorPartSource(PartSource_nVar,L_2_PartSource,L_Inf_PartSource)
#endif /*PARTICLES*/
  END IF
END IF

! the following analysis are restricted to Runge-Kutta based time-discs and temporal varying electrodynamic fields
#if defined(LSERK) || defined(IMPA) || defined(ROS) || USE_HDG

!----------------------------------------------------------------------------------------------------------------------------------
! Maxwell's equation: Compute Poynting Vector and field energies
!----------------------------------------------------------------------------------------------------------------------------------
IF (DoFieldAnalyze) THEN
  IF(DoPerformFieldAnalyze) CALL AnalyzeField(OutputTime)
END IF

!----------------------------------------------------------------------------------------------------------------------------------
! Recordpoints buffer
! Only usable with
! maxwell or poisson in combination with HDG
!----------------------------------------------------------------------------------------------------------------------------------
! remove duplicate analysis  due to double call of PerformAnalysis
IF(RP_onProc) THEN
  ! for a restart, an initial analysis at t=tRestart is necessary to
  ! compute fill the HDF5 file correctly
  ! Information to time in HDF5-Format:
  ! file1: 0 :t1
  ! file2: t1:tend
  ! Hence, t1 and the fields are stored in both files
  IF(DoPerformFieldAnalyze.OR.iter.EQ.0)THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    CALL RecordPoints(OutputTime,OutputHDF5)
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_DGANALYZE,tLBStart)
#endif /*USE_LOADBALANCE*/
  END IF
END IF

! end the analyzes for all Runge-Kutta based time-discs
#endif /* LSERK && IMPA && ROS && USE_HDG*/

!----------------------------------------------------------------------------------------------------------------------------------
! PIC, DSMC and other Particle-based Solvers
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef PARTICLES
! The following if clause might be simplified
IF((OutPutHDF5.OR.FirstOrLastIter) .AND. .NOT.(FirstOrLastIter.AND.(.NOT.OutPutHDF5).AND.(iter.NE.0))) CALL CalculatePartElemData()
#endif /*PARTICLES*/

#ifdef PARTICLES
IF(DoPartAnalyze.AND.DoPerformPartAnalyze)           CALL AnalyzeParticles(OutputTime)
IF(DoPerformSurfaceAnalyze)                          CALL AnalyzeSurface(OutputTime)
IF(TrackParticlePosition.AND.DoPerformPartAnalyze)   CALL WriteParticleTrackingData(OutputTime,iter) ! new function
#ifdef CODE_ANALYZE
IF(DoInterpolationAnalytic.AND.DoPerformPartAnalyze) CALL AnalyticParticleMovement(OutputTime,iter)
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/


!----------------------------------------------------------------------------------------------------------------------------------
! DSMC
!----------------------------------------------------------------------------------------------------------------------------------
! update of time here
#ifdef PARTICLES

! write volume data for DSMC macroscopic values
IF ((WriteMacroVolumeValues).AND.(.NOT.OutputHDF5))THEN
  CALL DSMC_data_sampling()
  IF (iter.GT.0) iter_macvalout = iter_macvalout + 1
  IF (MacroValSamplIterNum.LE.iter_macvalout) THEN
    CALL WriteDSMCToHDF5(TRIM(MeshFile),OutputTime)
    iter_macvalout = 0
    DSMC%SampNum = 0
    DSMC_Solution = 0.0
    IF(SampleElecExcitation) ExcitationSampleData = 0.0
    IF(DSMC%CalcQualityFactors) THEN
      DSMC%QualityFacSamp(:,:) = 0.
      IF(BGKInitDone) BGK_QualityFacSamp(:,:) = 0.
      IF(FPInitDone) FP_QualityFacSamp(:,:) = 0.
    END IF
  END IF
END IF

! write surface data for DSMC macroscopic values
IF ((WriteMacroSurfaceValues).AND.(.NOT.OutputHDF5))THEN
  IF (iter.GT.0) iter_macsurfvalout = iter_macsurfvalout + 1
  IF (MacroValSamplIterNum.LE.iter_macsurfvalout) THEN
    CALL CalcSurfaceValues()
    DO iSide=1,nComputeNodeSurfTotalSides
      SampWallState(:,:,:,iSide)=0.
      IF(nPorousBC.GT.0) THEN
        SampWallPumpCapacity(iSide)=0.
      END IF
      ! Sampling of impact energy for each species (trans, rot, vib), impact vector (x,y,z) and angle
      IF(CalcSurfaceImpact)THEN
        SampWallImpactEnergy(:,:,:,:,iSide)=0.
        SampWallImpactVector(:,:,:,:,iSide)=0.
        SampWallImpactAngle (:,:,:,  iSide)=0.
        SampWallImpactNumber(:,:,:,  iSide)=0.
      END IF ! CalcSurfaceImpact
    END DO
    iter_macsurfvalout = 0
  END IF
END IF

! Output of macroscopic variables sampled with TimeFracForSampling, which is performed after each DSMC time step
IF(OutPutHDF5)THEN
  IF((LastIter).AND.(useDSMC).AND.(.NOT.WriteMacroVolumeValues).AND.(.NOT.WriteMacroSurfaceValues)) THEN
    IF (DSMC%NumOutput.GT.0) THEN
      CALL WriteDSMCToHDF5(TRIM(MeshFile),OutputTime)
    END IF
    IF ((OutputTime.GE.DelayTime).AND.(DSMC%NumOutput.GT.0)) THEN
      IF(DSMC%CalcSurfaceVal) CALL CalcSurfaceValues
    END IF
  END IF
END IF

! Measure tracking time for particles // no MPI barrier MPI Wall-time but local CPU time
! Allows non-synchronous measurement of particle tracking
IF(OutPutHDF5 .AND. MeasureTrackTime)THEN
#if USE_MPI
  IF(MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE  , nTracks       , 1 , MPI_INTEGER          , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE  , tTracking     , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE  , tLocalization , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  ELSE ! no Root
    CALL MPI_REDUCE(nTracks       , RECI          , 1 , MPI_INTEGER          , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(tTracking     , RECR          , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(tLocalization , RECR          , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  END IF
#endif /*USE_MPI*/
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

!----------------------------------------------------------------------------------------------------------------------------------
! Code Analyze
!----------------------------------------------------------------------------------------------------------------------------------

#ifdef CODE_ANALYZE
! particle analyze
IF(TrackingMethod.NE.TRIATRACKING)THEN
  IF (DoPartAnalyze)  THEN
    IF(DoPerformPartAnalyze) CALL CodeAnalyzeOutput(OutputTime)
    IF(LastIter)THEN
      CALL CodeAnalyzeOutput(OutputTime)
      SWRITE(UNIT_stdOut,'(A51)') 'CODE_ANALYZE: Following output has been accumulated'
      SWRITE(UNIT_stdOut,'(A35,E15.7)') ' rTotalBBChecks    : ' , rTotalBBChecks
      SWRITE(UNIT_stdOut,'(A35,E15.7)') ' rTotalBezierClips : ' , rTotalBezierClips
      SWRITE(UNIT_stdOut,'(A35,E15.7)') ' rTotalBezierNewton: ' , rTotalBezierNewton
!      TotalSideBoundingBoxVolume=SUM(SideBoundingBoxVolume)
!#if USE_MPI
!      IF(MPIRoot) THEN
!        CALL MPI_REDUCE(MPI_IN_PLACE, TotalSideBoundingBoxVolume , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
!      ELSE ! no Root
!        CALL MPI_REDUCE(TotalSideBoundingBoxVolume ,           0 , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
!      END IF
!#endif /*USE_MPI*/
!      SWRITE(UNIT_stdOut,'(A35,E15.7)') ' Total Volume of SideBoundingBox: ' , TotalSideBoundingBoxVolume
    END IF
  END IF
END IF ! TrackingMethod.NE.TRIATRACKING
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/

! Time for analysis
IF((DoPerformFieldAnalyze.OR.DoPerformPartAnalyze.OR.DoPerformSurfaceAnalyze).AND.DoMeasureAnalyzeTime)THEN
  EndAnalyzeTime=PICLASTIME()
  AnalyzeTime = AnalyzeTime + EndAnalyzeTime-StartAnalyzeTime
END IF

!----------------------------------------------------------------------------------------------------------------------------------
! Output info
!----------------------------------------------------------------------------------------------------------------------------------
IF(DoPerformErrorCalc)THEN
  IF(DoCalcErrorNorms) THEN
    ! Graphical output
    IF(MPIroot) THEN
      WRITE(formatStr,'(A5,I1,A7)')'(A13,',PP_nVar,'ES16.7)'
      WRITE(UNIT_StdOut,formatStr)' L_2       : ',L_2_Error
      WRITE(UNIT_StdOut,formatStr)' L_inf     : ',L_Inf_Error
    END IF
  END IF
#ifdef PARTICLES
  IF (DoDeposition .AND. RelaxDeposition) THEN
    ! Graphical output
    IF(MPIroot) THEN
      WRITE(formatStr,'(A5,I1,A7)')'(A13,',PartSource_nVar,'ES16.7)'
      WRITE(UNIT_StdOut,formatStr)' L_2_PartSource  : ',L_2_PartSource
      WRITE(UNIT_StdOut,formatStr)' L_inf_PartSource: ',L_Inf_PartSource
    END IF
  END IF
#endif /* PARTICLES */
  IF(EndAnalyzeTime.LT.0.0) EndAnalyzeTime=PICLASTIME()
  IF(MPIroot) THEN
    ! write out has to be "Sim time" due to analyzes in reggie. Reggie searches for exactly this tag
    WRITE(UNIT_StdOut,'(A13,ES16.7)')' Sim time  : ',OutputTime
    IF(DoMeasureAnalyzeTime)THEN
      WRITE(UNIT_StdOut,'(A17,ES16.7,A9,I11,A)')' Analyze time  : ',AnalyzeTime, ' (called ',AnalyzeCount,' times)'
      AnalyzeCount = 0
      AnalyzeTime  = 0.0
    END IF ! DoMeasureAnalyzeTime
    IF (OutputTime.GT.0.) THEN
      WRITE(UNIT_StdOut,'(132("."))')
      CALL DisplayMessageAndTime(EndAnalyzeTime-StartTime, 'PICLAS RUNNING '//TRIM(ProjectName)//'... ', DisplayDespiteLB=.TRUE.)
    ELSE
      WRITE(UNIT_StdOut,'(132("="))')
    END IF
  END IF
END IF

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
END SUBROUTINE PerformAnalyze

#ifdef PARTICLES
#ifdef CODE_ANALYZE
SUBROUTINE CodeAnalyzeOutput(TIME)
!===================================================================================================================================
! output of code_analyze stuff: costly analyze routines for sanity checks and debugging
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars            ,ONLY:DoCodeAnalyzeOutput
USE MOD_Particle_Analyze_Vars   ,ONLY:IsRestart,DoPartAnalyze
USE MOD_Restart_Vars            ,ONLY:DoRestart
USE MOD_Particle_Surfaces_Vars  ,ONLY:rBoundingBoxChecks,rPerformBezierClip,rTotalBBChecks,rTotalBezierClips,rPerformBezierNewton
USE MOD_Particle_Surfaces_Vars  ,ONLY:rTotalBezierNewton
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: rDummy
LOGICAL             :: isOpen
CHARACTER(LEN=350)  :: outfile
INTEGER             :: unit_index, OutputCounter
!===================================================================================================================================
IF(.NOT.DoCodeAnalyzeOutput) RETURN ! check if the output is to be skipped and return if true

IF ( DoRestart ) THEN
  isRestart = .true.
END IF
IF (DoPartAnalyze) THEN
  !SWRITE(UNIT_StdOut,'(132("-"))')
  !SWRITE(UNIT_stdOut,'(A)') ' PERFORMING PARTICLE ANALYZE...'
  OutputCounter = 2
  unit_index = 555
#if USE_MPI
  IF(MPIROOT)THEN
#endif /*USE_MPI*/
   INQUIRE(UNIT   = unit_index , OPENED = isOpen)
   IF (.NOT.isOpen) THEN
     outfile = 'CodeAnalyze.csv'
     IF (isRestart .and. FILEEXISTS(outfile)) THEN
        OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
        !CALL FLUSH (unit_index)
     ELSE
        OPEN(unit_index,file=TRIM(outfile))
        !CALL FLUSH (unit_index)
        !--- insert header

        WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I3.3,A11)',ADVANCE='NO') OutputCounter,'-nBBChecks     '
          OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A1)',ADVANCE='NO') ','
        WRITE(unit_index,'(I3.3,A12)',ADVANCE='NO') OutputCounter,'-nBezierClips   '
          OutputCounter = OutputCounter + 1
        WRITE(unit_index,'(A14)') ' '
     END IF
   END IF
#if USE_MPI
  END IF
#endif /*USE_MPI*/

 ! MPI Communication
#if USE_MPI
  IF(MPIRoot) THEN
    CALL MPI_REDUCE(MPI_IN_PLACE         , rBoundingBoxChecks   , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE         , rPerformBezierClip   , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(MPI_IN_PLACE         , rPerformBezierNewton , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  ELSE ! no Root
    CALL MPI_REDUCE(rBoundingBoxChecks   , rDummy               , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(rPerformBezierClip   , rDummy               , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
    CALL MPI_REDUCE(rPerformBezierNewton , rDummy               , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
  END IF
#endif /*USE_MPI*/

#if USE_MPI
   IF(MPIROOT)THEN
#endif /*USE_MPI*/
     WRITE(unit_index,104,ADVANCE='NO') Time
     WRITE(unit_index,'(A1)',ADVANCE='NO') ','
     WRITE(unit_index,104,ADVANCE='NO') rBoundingBoxChecks
     WRITE(unit_index,'(A1)',ADVANCE='NO') ','
     WRITE(unit_index,104,ADVANCE='NO') rPerformBezierClip
     WRITE(unit_index,'(A1)',ADVANCE='NO') ','
     WRITE(unit_index,104,ADVANCE='NO') rPerformBezierNewton
     WRITE(unit_index,'(A1)') ' '
#if USE_MPI
   END IF
#endif /*USE_MPI*/

104    FORMAT (e25.14)

!SWRITE(UNIT_stdOut,'(A)')' PARTCILE ANALYZE DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')
ELSE
!SWRITE(UNIT_stdOut,'(A)')' NO PARTCILE ANALYZE TO DO!'
!SWRITE(UNIT_StdOut,'(132("-"))')
END IF ! DoPartAnalyze

! nullify and save total number
rTotalBBChecks=rTotalBBChecks+REAL(rBoundingBoxChecks,16)
rBoundingBoxChecks=0.
rTotalBezierClips=rTotalBezierClips+REAL(rPerformBezierClip,16)
rPerformBezierClip=0.
rTotalBezierNewton=rTotalBezierNewton+REAL(rPerformBezierNewton,16)
rPerformBezierNewton=0.

END SUBROUTINE CodeAnalyzeOutput
#endif /*CODE_ANALYZE*/
#endif /*PARTICLES*/

#if USE_HDG
!===================================================================================================================================
!> Create containers and communicators for each boundary on which the electric displacement current is calculated and agglomerated
!> This is done for all normal BCs except periodic BCs.
!>
!
!> 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
!>     contribute to the total electric displacement current. If yes, then this processor is part of the communicator
!> 2.) Create Mapping from electric displacement current BC index to field BC index
!> 3.) Create Mapping from field BC index to electric displacement current BC index
!> 4.) Check if field BC is on current proc (or MPI root)
!> 5.) Create MPI sub-communicators
!===================================================================================================================================
SUBROUTINE InitCalcElectricTimeDerivativeSurface()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars        ,ONLY: nBCs,BoundaryType
USE MOD_Analyze_Vars     ,ONLY: DoFieldAnalyze,CalcElectricTimeDerivative,EDC
USE MOD_Equation_Vars    ,ONLY: Et
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_MPI
USE MOD_Mesh_Vars        ,ONLY: BoundaryName,nBCSides,BC
#endif /*USE_MPI*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iBC,iEDCBC
#if USE_MPI
LOGICAL,ALLOCATABLE :: BConProc(:)
INTEGER             :: SideID,color
#endif /*USE_MPI*/
!===================================================================================================================================
IF(.NOT.CalcElectricTimeDerivative) RETURN ! Read-in parameter that is set in  InitAnalyze() in analyze.f90

! Allocate temporal derivative for E: No need to nullify as is it overwritten with E the first time it is used
ALLOCATE(Et(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Et = 0.

! 1.) Loop over all field BCs and check if the current processor is either the MPI root or has at least one of the BCs that
! contribute to the total electric displacement current. If yes, then this processor is part of the communicator
EDC%NBoundaries = 0
DO iBC=1,nBCs
  IF(BoundaryType(iBC,BC_ALPHA).NE.0) CYCLE
  EDC%NBoundaries = EDC%NBoundaries + 1
END DO

! If not electric displacement current boundaries exist, no measurement of the current can be performed
IF(EDC%NBoundaries.EQ.0) RETURN

! Automatically activate surface model analyze flag
DoFieldAnalyze = .TRUE.

! 2.) Create Mapping from electric displacement current BC index to field BC index
ALLOCATE(EDC%FieldBoundaries(EDC%NBoundaries))
EDC%NBoundaries = 0
DO iBC=1,nBCs
  IF(BoundaryType(iBC,BC_ALPHA).NE.0) CYCLE
  EDC%NBoundaries = EDC%NBoundaries + 1
  EDC%FieldBoundaries(EDC%NBoundaries) = iBC
END DO

! Allocate the container
ALLOCATE(EDC%Current(1:EDC%NBoundaries))
EDC%Current = 0.

! 3.) Create Mapping from field BC index to electric displacement current BC index
ALLOCATE(EDC%BCIDToEDCBCID(nBCs))
EDC%BCIDToEDCBCID = -1
DO iEDCBC = 1, EDC%NBoundaries
  iBC = EDC%FieldBoundaries(iEDCBC)
  EDC%BCIDToEDCBCID(iBC) = iEDCBC
END DO ! iEDCBC = 1, EDC%NBoundaries

#if USE_MPI
! 4.) Check if field BC is on current proc (or MPI root)
ALLOCATE(BConProc(EDC%NBoundaries))
BConProc = .FALSE.
IF(MPIRoot)THEN
  BConProc = .TRUE.
ELSE
  DO SideID=1,nBCSides
    IF(BoundaryType(BC(SideID),BC_ALPHA).NE.0) CYCLE
    iBC     = BC(SideID)
    iEDCBC  = EDC%BCIDToEDCBCID(iBC)
    BConProc(iEDCBC) = .TRUE.
  END DO ! SideID=1,nBCSides
END IF ! MPIRoot

! 5.) Create MPI sub-communicators
ALLOCATE(EDC%COMM(EDC%NBoundaries))
DO iEDCBC = 1, EDC%NBoundaries
  ! create new communicator
  color = MERGE(iEDCBC, MPI_UNDEFINED, BConProc(iEDCBC))

  ! set communicator id
  EDC%COMM(iEDCBC)%ID=iEDCBC

  ! create new emission communicator for electric displacement current communication. Pass MPI_INFO_NULL as rank to follow the original ordering
  CALL MPI_COMM_SPLIT(MPI_COMM_PICLAS, color, MPI_INFO_NULL, EDC%COMM(iEDCBC)%UNICATOR, iError)

  ! Find my rank on the shared communicator, comm size and proc name
  IF(BConProc(iEDCBC))THEN
    CALL MPI_COMM_RANK(EDC%COMM(iEDCBC)%UNICATOR, EDC%COMM(iEDCBC)%MyRank, iError)
    CALL MPI_COMM_SIZE(EDC%COMM(iEDCBC)%UNICATOR, EDC%COMM(iEDCBC)%nProcs, iError)

    ! inform about size of emission communicator
    IF (EDC%COMM(iEDCBC)%MyRank.EQ.0) THEN
#if USE_LOADBALANCE
      IF(.NOT.PerformLoadBalance)&
#endif /*USE_LOADBALANCE*/
          WRITE(UNIT_StdOut,'(A,I0,A,I0,A)') ' Electric displacement current: Emission-Communicator ',iEDCBC,' on ',&
              EDC%COMM(iEDCBC)%nProcs,' procs for '//TRIM(BoundaryName(EDC%FieldBoundaries(iEDCBC)))
    END IF
  END IF ! BConProc(iEDCBC)
END DO ! iEDCBC = 1, EDC%NBoundaries
DEALLOCATE(BConProc)
#endif /*USE_MPI*/


END SUBROUTINE InitCalcElectricTimeDerivativeSurface
#endif /*USE_HDG*/

END MODULE MOD_Analyze
