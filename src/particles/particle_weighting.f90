!==================================================================================================================================
! Copyright (c) 2024 boltzplatz - numerical plasma dynamics GmbH, Simone Lauterbach, Marcel Pfeiffer
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

MODULE MOD_ParticleWeighting
!===================================================================================================================================
!>
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
PUBLIC :: DefineParametersParticleWeighting, InitParticleWeighting
!===================================================================================================================================

INTEGER,PARAMETER      :: PRM_PARTWEIGHT_CONSTANT   = 0
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_RADIAL     = 1
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_LINEAR     = 2
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_CELLLOCAL  = 3

CONTAINS

!==================================================================================================================================
!> Define parameters for particles weighting
!==================================================================================================================================
SUBROUTINE DefineParametersParticleWeighting()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE

CALL prms%SetSection("Particle Weighting")
CALL prms%CreateIntFromStringOption('Part-Weight-Type', "Particle weighting type: \n"             //&
                                      'constant ('//TRIM(int2strf(PRM_PARTWEIGHT_CONSTANT))//')\n'  //&
                                      'radial ('//TRIM(int2strf(PRM_PARTWEIGHT_RADIAL))//')\n'      //&
                                      'linear ('//TRIM(int2strf(PRM_PARTWEIGHT_LINEAR))//')\n'      //&
                                      'cell_local ('//TRIM(int2strf(PRM_PARTWEIGHT_CELLLOCAL))//')' &
                                      ,'constant')
CALL addStrListEntry('Part-Weight-Type' , 'constant'   , PRM_PARTWEIGHT_CONSTANT)
CALL addStrListEntry('Part-Weight-Type' , 'radial'     , PRM_PARTWEIGHT_RADIAL)
CALL addStrListEntry('Part-Weight-Type' , 'linear'     , PRM_PARTWEIGHT_LINEAR)
CALL addStrListEntry('Part-Weight-Type' , 'cell_local' , PRM_PARTWEIGHT_CELLLOCAL)

! General particle weighting parameters
CALL prms%CreateIntOption(    'Part-Weight-CloneMode',  &
                              'Radial weighting: Select between methods for the delayed insertion of cloned particles:/n'//&
                              '1: Chronological, 2: Random', '2')
CALL prms%CreateIntOption(    'Part-Weight-CloneDelay', &
                              'Radial weighting:  Delay (number of iterations) before the stored cloned particles are inserted '//&
                              'at the position they were cloned', '2')
CALL prms%CreateLogicalOption('Part-Weight-CellAverage', 'Radial/linear only: Enables a cell-average weighting (), '//&
                              'where every particle has the same weighting factor within a cell', '.FALSE.')
CALL prms%CreateIntOption(    'Part-Weight-SurfFluxSubSides', &
                              'Axisymmetric only: Split the surface flux side into the given number of subsides, reduces the '//&
                              'error in the particle distribution across the cell (visible in the number density)', '20')

! Radial weighting parameter
CALL prms%CreateRealOption(   'Part-Weight-Radial-ScaleFactor', 'Radial weighting scale factor, defining '//&
                              'the linear increase of the weighting factor (e.g. factor 2 means that the weighting factor will '//&
                              'be twice as large at the outer radial boundary or at the end point.')

! Linear weighting parameters
CALL prms%CreateIntOption(    'Part-Weight-Linear-nScalePoints', 'Number of coordinates with distinct weighting factors', '2')
CALL prms%CreateIntOption(    'Part-Weight-Linear-CoordinateAxis', '1: x-Axis, 2: y-Axis, 3: z-Axis', '0')
CALL prms%CreateRealArrayOption('Part-Weight-Linear-StartPointForScaling', &
                              'Start coordinate for the scaling along a given vector', no=3)
CALL prms%CreateRealArrayOption('Part-Weight-Linear-EndPointForScaling', &
                              'End coordinate for the scaling along a given vector', no=3)
CALL prms%CreateRealOption(   'Part-Weight-Linear-ScalePoint[$]-Coordinate', &
                              '(Relative ) Coordinate of the respective scale point on the axis', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(   'Part-Weight-Linear-ScalePoint[$]-Factor', &
                              'Weighting factor of the respective scale point', numberedmulti=.TRUE.)

! Cell-local weighting parameters
CALL prms%CreateRealOption(   'Part-Weight-CellLocal-MinParticleNumber', 'Target minimum simulation particle number per cell', '5')
CALL prms%CreateRealOption(   'Part-Weight-CellLocal-MaxParticleNumber', 'Target maximum simulation particle number per cell', '1000')
CALL prms%CreateLogicalOption('Part-Weight-CellLocal-ApplyMedianFilter', 'Applies a median filter to the cell-local distribution  '//&
                              'of the adapted weighting factor', '.FALSE.')
CALL prms%CreateIntOption(    'Part-Weight-CellLocal-RefinementNumber', 'Number of times the median filter is applied', '1')
CALL prms%CreateRealOption(   'Part-Weight-CellLocal-QualityFactor', 'Threshold for the adaption of the weighting based on a quality factor:\n'//&
                              'DSMC: Mean collision separation distance over mean free path (DSMC_MCS_over_MFP)\n'//&
                              'BGK/FP: Maximal relaxation factor (BGK_MaxRelaxationFactor/FP_MaxRelaxationFactor)\n'//&
                              'If these are below the threshold, the weighting factor will be adapted to resolve it', '0.8')
CALL prms%CreateIntOption(    'Part-Weight-CellLocal-SymAxis-MinPartNum', 'Target minimum particle number close to the symmetry axis', '10')
CALL prms%CreateIntOption(    'Part-Weight-CellLocal-Cat-MinPartNum', 'Target minimum particle number close to catalytic boundaries', '10')
CALL prms%CreateLogicalOption('Part-Weight-CellLocal-IncludeMaxPartNum', 'Flag to determine if the maximal particle number should be '//&
                              'included in the refinement process', '.TRUE.')
CALL prms%CreateLogicalOption('Part-Weight-CellLocal-SkipAdaption', 'Flag to skip the adaption process of the weighting factor '//&
                              'and only use previously determined weighting factor values from a state file', '.FALSE.')
CALL prms%CreateRealOption(   'Part-Weight-CellLocal-RefineFactorBGK', 'Ratio between the target BGK and DSMC weighting factor', '1.0')
CALL prms%CreateRealOption(   'Part-Weight-CellLocal-RefineFactorFP', 'Ratio between the target FP and DSMC weighting factor', '1.0')

END SUBROUTINE DefineParametersParticleWeighting


SUBROUTINE InitParticleWeighting()
!===================================================================================================================================
!> Initialize the particle weighting, especially the particle cloning
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools           ,ONLY: GETINT,GETREAL,GETLOGICAL,GETINTFROMSTR
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
USE MOD_DSMC_Vars             ,ONLY: DoRadialWeighting, DoLinearWeighting, DoCellLocalWeighting, ParticleWeighting
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: ParticleWeightType
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT PARTICLE WEIGHTING...'
DoRadialWeighting = .FALSE.
DoLinearWeighting = .FALSE.
DoCellLocalWeighting = .FALSE.
ParticleWeighting%PerformCloning = .FALSE.
ParticleWeighting%EnableOutput = .FALSE.
ParticleWeighting%UseCellAverage = .FALSE.
ParticleWeighting%UseSubdivision = .FALSE.

ParticleWeightType = GETINTFROMSTR('Part-Weight-Type')

SELECT CASE(ParticleWeightType)
Case(PRM_PARTWEIGHT_CONSTANT)
  ParticleWeighting%Type = 'constant'
Case(PRM_PARTWEIGHT_RADIAL)
  ParticleWeighting%Type = 'radial'
  DoRadialWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
  ParticleWeighting%ScaleFactor = GETREAL('Part-Weight-Radial-ScaleFactor')
  IF(ParticleWeighting%ScaleFactor.LT.1.) CALL CollectiveStop(__STAMP__,'Part-Weight-Radial-ScaleFactor has to be greater than 1!')
  IF(.NOT.Symmetry%Axisymmetric) CALL CollectiveStop(__STAMP__,' Part-Weight-Type = radial requires an axisymmetric simulation!')
Case(PRM_PARTWEIGHT_LINEAR)
  ParticleWeighting%Type = 'linear'
  DoLinearWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
Case(PRM_PARTWEIGHT_CELLLOCAL)
  ParticleWeighting%Type = 'cell_local'
  DoCellLocalWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown Part-Weight-Type!' ,IntInfo=ParticleWeightType)
END SELECT

! Cell local radial weighting (all particles have the same weighting factor within a cell)
IF(DoRadialWeighting.OR.DoLinearWeighting) THEN
  ParticleWeighting%UseCellAverage = GETLOGICAL('Part-Weight-CellAverage')
END IF

! Number of subsides to split the surface flux sides into, otherwise a wrong distribution of particles across large cells will be
! inserted, visible in the number density as an increase in the number density closer the axis (e.g. resulting in a heat flux peak)
! (especially when using mortar meshes)
IF(Symmetry%Axisymmetric.AND.ParticleWeighting%PerformCloning.AND..NOT.ParticleWeighting%UseCellAverage) THEN
  ParticleWeighting%UseSubdivision = .TRUE.
  ParticleWeighting%nSubSides=GETINT('Part-Weight-SurfFluxSubSides')
  ALLOCATE(ParticleWeighting%PartInsSide(ParticleWeighting%nSubSides))
  ParticleWeighting%PartInsSide = 0
END IF

! Initialize of particle cloning
IF(ParticleWeighting%PerformCloning) CALL InitParticleCloning()

! Initialize variables for linear weighting
IF(DoLinearWeighting) CALL InitLinearWeighting()

! Initialize variables for cell-local weighting
IF(DoCellLocalWeighting) CALL InitCellLocalWeighting()

END SUBROUTINE InitParticleWeighting


SUBROUTINE InitParticleCloning()
!===================================================================================================================================
!> Read-in and initialize the variables required for the cloning procedures. Two modes with a delayed clone insertion are available:
!> 1: Insert the clones after the delay in the same chronological order as they were created
!> 2: Choose a random list of particles to insert after the delay buffer is full with clones
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_DSMC_Vars               ,ONLY: ParticleWeighting, ClonedParticles
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: DoLoadBalance, UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Clone read-in during load balance is currently only supported via the HDF5 output
#if USE_LOADBALANCE
IF(DoLoadBalance.AND.(.NOT.UseH5IOLoadBalance)) THEN
  CALL abort(__STAMP__,'ERROR: Particle weighting only supports a load balance using an HDF5 output (UseH5IOLoadBalance = T)!')
END IF
#endif /*USE_LOADBALANCE*/

! Cloning parameters
ParticleWeighting%CloneMode = GETINT('Part-Weight-CloneMode')
ParticleWeighting%CloneInputDelay = GETINT('Part-Weight-CloneDelay')

ParticleWeighting%NextClone = 0
ParticleWeighting%CloneVecLengthDelta = 100
ParticleWeighting%CloneVecLength = ParticleWeighting%CloneVecLengthDelta

SELECT CASE(ParticleWeighting%CloneMode)
  CASE(1)
    IF(ParticleWeighting%CloneInputDelay.LT.1) THEN
      CALL Abort(__STAMP__,'ERROR in Particle Weighting: Clone delay should be greater than 0')
    END IF
    ALLOCATE(ParticleWeighting%ClonePartNum(0:(ParticleWeighting%CloneInputDelay-1)))
    ALLOCATE(ClonedParticles(1:ParticleWeighting%CloneVecLength,0:(ParticleWeighting%CloneInputDelay-1)))
    ParticleWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) ParticleWeighting%CloneDelayDiff = 1
  CASE(2)
    IF(ParticleWeighting%CloneInputDelay.LT.2) THEN
      CALL Abort(__STAMP__,'ERROR in Particle Weighting: Clone delay should be greater than 1')
    END IF
    ALLOCATE(ParticleWeighting%ClonePartNum(0:ParticleWeighting%CloneInputDelay))
    ALLOCATE(ClonedParticles(1:ParticleWeighting%CloneVecLength,0:ParticleWeighting%CloneInputDelay))
    ParticleWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) ParticleWeighting%CloneDelayDiff = 0
  CASE DEFAULT
    CALL Abort(__STAMP__,'ERROR in Particle Weighting: The selected cloning mode is not available! Choose between 1 and 2.'//&
        ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

END SUBROUTINE InitParticleCloning


SUBROUTINE InitLinearWeighting()
!===================================================================================================================================
!> Read-in and initialize the variables required for the linear weighting, calculate the average weighting factor for the emission
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_DSMC_Vars               ,ONLY: LinearWeighting, ParticleWeighting
USE MOD_part_tools              ,ONLY: CalcAverageMPF
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)                   :: hilf
INTEGER                         :: iScale, nScalePoints
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT LINEAR PARTICLE WEIGHTING...'
! Linear scaling points for the variable weighting
LinearWeighting%nScalePoints = GETINT('Part-Weight-Linear-nScalePoints')

nScalePoints = LinearWeighting%nScalePoints
ALLOCATE(LinearWeighting%ScalePoint(nScalePoints))
ALLOCATE(LinearWeighting%VarMPF(nScalePoints))

! Variable weights with scaling along one of the coordinate axes
LinearWeighting%ScaleAxis    = GETINT('Part-Weight-Linear-CoordinateAxis')

! Variable weights with scaling not along a coordinate axis
IF (LinearWeighting%ScaleAxis.EQ.0) THEN
  LinearWeighting%StartPointScaling = GETREALARRAY('Part-Weight-Linear-StartPointForScaling',3)
  LinearWeighting%EndPointScaling   = GETREALARRAY('Part-Weight-Linear-EndPointForScaling',3)
  ! Determine the vector along which the scaling is performed
  LinearWeighting%ScalingVector = LinearWeighting%EndPointScaling - LinearWeighting%StartPointScaling
END IF

! Read-In of the scaling points along the chosen direction and the corresponding weighting factor
DO iScale = 1, nScalePoints
  WRITE(UNIT=hilf,FMT='(I0)') iScale
  LinearWeighting%ScalePoint(iScale) = GETREAL('Part-Weight-Linear-ScalePoint'//TRIM(hilf)//'-Coordinate')
  LinearWeighting%VarMPF(iScale)     = GETREAL('Part-Weight-Linear-ScalePoint'//TRIM(hilf)//'-Factor')
END DO

! Sanity check: Accept only relative values for scaling a user-defined vector
IF(LinearWeighting%ScaleAxis.EQ.0) THEN
  IF(ANY(LinearWeighting%ScalePoint.GT.1.0).OR.ANY(LinearWeighting%ScalePoint.LT.0.0)) CALL abort(__STAMP__,&
    'ERROR in InitLinearWeighting: The coordinate of the scale point has to be a relative value along the vector between 0 and 1!')
END IF

! Sanity check: Check whether the coordinate axis has been properly defined or disabled
IF (LinearWeighting%ScaleAxis.LT.0.OR.LinearWeighting%ScaleAxis.GT.3) THEN
  CALL abort(__STAMP__, 'ERROR in InitLinearWeighting: The coordinate axis must be a value between 0 (disabled), 1 (=x), 2 (=y), and 3 (=z)!')
END IF

! Sanity check: Check whether the coordinate axis has been properly defined or disabled
IF (LinearWeighting%nScalePoints.LT.2) THEN
  CALL abort(__STAMP__, 'ERROR in InitLinearWeighting: The number of scaling points must at least be two!')
END IF

! Calculation of the average particle MPF in the simulation domain (utilized for the particle initialization)
ParticleWeighting%ScaleFactor = CalcAverageMPF()

END SUBROUTINE InitLinearWeighting


SUBROUTINE InitCellLocalWeighting()
!===================================================================================================================================
!> Read-in and initialize the variables required for the cell-local weighting
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Restart_Vars            ,ONLY: DoMacroscopicRestart
USE MOD_DSMC_Vars               ,ONLY: CellLocalWeight
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT CELL-LOCAL PARTICLE WEIGHTING...'
! No further adaption process, use of the MPF distribution from the previous adaption process
IF(.NOT.PerformLoadBalance) CellLocalWeight%SkipAdaption       = GETLOGICAL('Part-Weight-CellLocal-SkipAdaption')

IF (.NOT.CellLocalWeight%SkipAdaption) THEN
  IF(.NOT.DoMacroscopicRestart) CALL abort(__STAMP__, 'ERROR: Cell-local weighting adaption process only possible with -DoMacroscopicRestart=T!')
  ! Read-in of the parameter boundaries
  CellLocalWeight%MinPartNum         = GETREAL('Part-Weight-CellLocal-MinParticleNumber')
  CellLocalWeight%MaxPartNum         = GETREAL('Part-Weight-CellLocal-MaxParticleNumber')
  CellLocalWeight%IncludeMaxPartNum  = GETLOGICAL('Part-Weight-CellLocal-IncludeMaxPartNum')
  CellLocalWeight%QualityFactor      = GETREAL('Part-Weight-CellLocal-QualityFactor')
  CellLocalWeight%BGKFactor          = GETREAL('Part-Weight-CellLocal-RefineFactorBGK')
  CellLocalWeight%FPFactor           = GETREAL('Part-Weight-CellLocal-RefineFactorFP')
  CellLocalWeight%SymAxis_MinPartNum = GETINT('Part-Weight-CellLocal-SymAxis-MinPartNum')
  CellLocalWeight%Cat_MinPartNum     = GETINT('Part-Weight-CellLocal-Cat-MinPartNum')
  CellLocalWeight%UseMedianFilter    = GETLOGICAL('Part-Weight-CellLocal-ApplyMedianFilter')
  ! Parameters for the filtering subroutine
  IF (CellLocalWeight%UseMedianFilter) CellLocalWeight%nRefine          = GETINT('Part-Weight-CellLocal-RefinementNumber')
END IF

END SUBROUTINE InitCellLocalWeighting

END MODULE MOD_ParticleWeighting