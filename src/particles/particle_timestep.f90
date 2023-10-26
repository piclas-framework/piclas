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

MODULE MOD_Particle_TimeStep
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
PUBLIC :: DefineParametersVariableTimeStep
PUBLIC :: InitPartTimeStep, GetParticleTimeStep, GetSpeciesTimeStep, VarTimeStep_CalcElemFacs, VarTimeStep_InitDistribution

!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for particles
!==================================================================================================================================
SUBROUTINE DefineParametersVariableTimeStep()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================

CALL prms%SetSection("Variable Timestep")
! === Distribution
CALL prms%CreateLogicalOption('Part-VariableTimeStep-Distribution', &
                              'Utilize a time step distribution, must be available in the particle state file!', '.FALSE.')
CALL prms%CreateLogicalOption('Part-VariableTimeStep-Distribution-Adapt', &
                              'Adapt the time step distribution according to certain parameters (read-in from a DSMC state) '//&
                              'and store it in the particle state file. Requires a macroscopic restart and a DSMC state file:\n'//&
                              'Particles-MacroscopicRestart = T\n'//&
                              'Particles-MacroscopicRestart-Filename = DSMCState.h5', '.FALSE.')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-Distribution-TargetMCSoverMFP', &
                              'DSMC: Target ratio of the mean collision separation distance over the mean free path', '0.25')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-Distribution-TargetMaxCollProb', &
                              'DSMC: Target maximum collision probability', '0.8')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-Distribution-TargetMaxRelaxFactor', &
                              'BGK: Target maximum relaxation factor', '0.8')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-Distribution-MaxFactor', &
                              'Maximum time factor to avoid too large time steps and problems with halo region/particle cloning')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-Distribution-MinFactor', &
                              'Minimum time factor to avoid cells with a large number of particles')
CALL prms%CreateIntOption(    'Part-VariableTimeStep-Distribution-MinPartNum', &
                              'Optional: Define a minimum number of particles per cells to increase the number of particles by '//&
                              'decreasing the time step', '0')
! === Linear Scaling
CALL prms%CreateLogicalOption('Part-VariableTimeStep-LinearScaling', &
                              'Utilize a linearly increasing time step along a given direction (3D) or a linearly increasing '//&
                              'time step in the radial direction (2D)', '.FALSE.')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-ScaleFactor', &
                              'Time step factor f*dt, determines the maximal time step of the linear function')
! 3D: Scaling along a given vector
CALL prms%CreateRealArrayOption('Part-VariableTimeStep-Direction', &
                                'Direction of the vector along which a linear increase is applied to the time step. '//&
                                'Currently only scaling along the x-axis (positive or negative direction) is allowed, '//&
                                'e.g. (/-1.0,0.0,0.0/)', no=3)
CALL prms%CreateRealArrayOption('Part-VariableTimeStep-StartPoint', &
                                'Starting point of the vector, to use the domain border: -99999.',no=3)
CALL prms%CreateRealArrayOption('Part-VariableTimeStep-EndPoint'  , &
                                'End point of the vector, to use the domain border: -99999.', no=3)
! 2D/Axi: Radial and axial scaling towards
CALL prms%CreateLogicalOption('Part-VariableTimeStep-Use2DFunction', &
                              'Only 2D/Axi simulations: Enables the scaling of the time step in the x-direction towards and '//&
                              'away from a user-given stagnation point', '.FALSE.')
CALL prms%CreateLogicalOption('Part-VariableTimeStep-OnlyDecreaseDt', &
                              'Only 2D/Axi simulations: Enables the scaling of the time step in the x-direction towards and '//&
                              'away from a user-given stagnation point', '.FALSE.')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-StagnationPoint', &
                              'Defines the point on the x-axis, towards which the time step is decreased with the factor '//&
                              'ScaleFactor2DFront and away from which the time step is increased with the factor ScaleFactor2DBack')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-ScaleFactor2DFront', &
                              'FRONT: Time step decreases towards the stagnation point')
CALL prms%CreateRealOption(   'Part-VariableTimeStep-ScaleFactor2DBack', &
                              'BACK: Time step increases away from the stagnation points')
! === Species-specific time step (activated through e.g. Species1-TimeStepFactor = 0.1)
CALL prms%CreateLogicalOption('Part-VariableTimeStep-DisableForMCC', &
                              'Disable the variable time step for the MCC routines to perform collisions at the manual '//&
                              'time step (e.g. to accelerate convergence to thermal/chemical equilibrium)', '.FALSE.')

END SUBROUTINE DefineParametersVariableTimeStep


SUBROUTINE InitPartTimeStep()
!===================================================================================================================================
!> Initialization of the variable time step (only in the case of UseLinearScaling or UseDistribution)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_ReadInTools             ,ONLY: GETLOGICAL, GETINT, GETREAL, GETREALARRAY
USE MOD_Particle_Vars           ,ONLY: Symmetry, VarTimeStep
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT VARIABLE TIME STEP...'

IF(VarTimeStep%UseLinearScaling) THEN
  IF(VarTimeStep%UseDistribution) CALL abort(__STAMP__, &
    'ERROR: Cannot use linear time scaling with a given distribution!')
  ! Timestep is varied according to a linear function
  VarTimeStep%ScaleFac = GETREAL('Part-VariableTimeStep-ScaleFactor')
  IF(Symmetry%Order.EQ.2) THEN
    ! Scaling the time step in x- and r-direction with the defintion of a stagnation point
    VarTimeStep%Use2DTimeFunc = GETLOGICAL('Part-VariableTimeStep-Use2DFunction')
    IF(VarTimeStep%Use2DTimeFunc) THEN
      ! FRONT: Time step decreases towards the stagnation point
      ! BACK: Time step increases away from the stagnation points
      VarTimeStep%StagnationPoint = GETREAL('Part-VariableTimeStep-StagnationPoint','0.0')
      VarTimeStep%TimeScaleFac2DFront = GETREAL('Part-VariableTimeStep-ScaleFactor2DFront','1.0')
      VarTimeStep%TimeScaleFac2DBack = GETREAL('Part-VariableTimeStep-ScaleFactor2DBack','1.0')
    END IF
  ELSE IF(Symmetry%Order.EQ.1) THEN
    CALL abort(__STAMP__, &
    'ERROR: 1D and variable timestep is not implemented yet!')
  ELSE
    VarTimeStep%StartPoint = GETREALARRAY('Part-VariableTimeStep-StartPoint',3)
    VarTimeStep%EndPoint = GETREALARRAY('Part-VariableTimeStep-EndPoint',3)
    VarTimeStep%Direction = GETREALARRAY('Part-VariableTimeStep-Direction',3)
    IF(ABS(VarTimeStep%Direction(1)).NE.1.0) CALL abort(&
      __STAMP__&
      ,'ERROR: Currently direction of linear time step scaling must be defined with (/1.0,0.0,0.0/)!')
  END IF
END IF
IF(VarTimeStep%UseDistribution) THEN
  ! Read-in of the maximal collision probability from the DSMC state file and calculate the appropriate time step
  ! Particle time step is utilized for this purpose, although element-wise time step is stored in VarTimeStep%ElemFacs
  ! Flag if the time step distribution should be adapted (else read-in, if array does not exist: abort)
  VarTimeStep%AdaptDistribution = GETLOGICAL('Part-VariableTimeStep-Distribution-Adapt')
  VarTimeStep%OnlyDecreaseDt = GETLOGICAL('Part-VariableTimeStep-OnlyDecreaseDt')
  VarTimeStep%TargetMCSoverMFP = GETREAL('Part-VariableTimeStep-Distribution-TargetMCSoverMFP')
  VarTimeStep%TargetMaxCollProb = GETREAL('Part-VariableTimeStep-Distribution-TargetMaxCollProb')
  ! Read of maximal time factor to avoid too large time steps and problems with halo region/particle cloning
  VarTimeStep%DistributionMaxTimeFactor = GETREAL('Part-VariableTimeStep-Distribution-MaxFactor')
  VarTimeStep%DistributionMinTimeFactor = GETREAL('Part-VariableTimeStep-Distribution-MinFactor')
  ! Optional: Increase number of particles by decreasing the time step
  VarTimeStep%DistributionMinPartNum = GETINT('Part-VariableTimeStep-Distribution-MinPartNum')
  ! BGK/FP: Read-in of the target maximal relaxation factor
  VarTimeStep%TargetMaxRelaxFactor = GETREAL('Part-VariableTimeStep-Distribution-TargetMaxRelaxFactor')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitPartTimeStep


SUBROUTINE VarTimeStep_InitDistribution()
!===================================================================================================================================
!> Calculates/determines the variable time step element-wise
!>-------------------------------------------------------------
!> Every proc does this for the whole domain, time factor is used in readMesh to determine the weight and distribute load
!> accordingly.
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems
USE MOD_HDF5_Input,             ONLY: OpenDataFile,CloseDataFile,DatasetExists,ReadArray,ReadAttribute,GetDataProps
USE MOD_PARTICLE_Vars,          ONLY: VarTimeStep, Symmetry
USE MOD_Restart_Vars,           ONLY: DoRestart, RestartFile, DoMacroscopicRestart, MacroRestartFileName
USE MOD_StringTools             ,ONLY:STRICMP
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                           :: iElem, iVar
LOGICAL                           :: TimeStepExists, QualityExists, TimeStepModified
REAL, ALLOCATABLE                 :: DSMCQualityFactors(:,:), PartNum(:)
REAL                              :: TimeFracTemp
CHARACTER(LEN=255),ALLOCATABLE    :: VarNames_tmp(:)
INTEGER                           :: nVar_HDF5, N_HDF5, nVar_MaxCollProb, nVar_MCSoverMFP, nVar_TotalPartNum, nVar_TimeStep
REAL, ALLOCATABLE                 :: ElemData_HDF5(:,:)
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
INTEGER                           :: nVar_MaxRelaxFac
REAL, ALLOCATABLE                 :: MaxRelaxFactor(:)
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' INIT VARIABLE TIME STEP DISTRIBUTION...'

TimeStepExists = .FALSE.
QualityExists = .FALSE.
nVar_TimeStep = 0

IF(DoRestart) THEN
! Try to get the time step factor distribution directly from state file
  CALL OpenDataFile(TRIM(RestartFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL DatasetExists(File_ID,'ElemTimeStep',TimeStepExists)
  IF(TimeStepExists) THEN
    ! Allocate the array for the element-wise time step factor
    ALLOCATE(VarTimeStep%ElemFac(nGlobalElems))
    VarTimeStep%ElemFac = 1.0
    ! Read-in of the time step
    ASSOCIATE(nGlobalElems    => INT(nGlobalElems,IK))
      CALL ReadArray('ElemTimeStep',2,(/nGlobalElems, 1_IK/),0_IK,1,RealArray=VarTimeStep%ElemFac(1:nGlobalElems))
    END ASSOCIATE
    SWRITE(UNIT_stdOut,*)'Variable Time Step: Read-in of timestep distribution from state file.'
#if USE_MPI
    ! Allocate the array for the element-wise weighting factor
    ALLOCATE(VarTimeStep%ElemWeight(nGlobalElems))
    VarTimeStep%ElemWeight = 1.0
#endif
  ELSEIF(.NOT.VarTimeStep%AdaptDistribution) THEN
    CALL abort(__STAMP__, &
      'ERROR: Variable time step requires a given timestep distribution or -Distribution-Adapt=T!')
  END IF
  CALL CloseDataFile()
ELSE  ! No Restart
  IF(.NOT.VarTimeStep%AdaptDistribution) THEN
    SWRITE(UNIT_stdOut,'(A)') '| Variable Time Step: No restart with given distribution and no adaption of the time step selected.'
    SWRITE(UNIT_stdOut,'(A)') '| Variable Time Step: Time step distribution is initialized uniformly with 1.'
    ! This is performed after each processor knows how many elements it will receive.
  END IF
END IF ! DoRestart = T/F

IF(VarTimeStep%AdaptDistribution) THEN
  SWRITE(UNIT_stdOut,'(A)') &
    ' | Variable Time Step: Adapting the time step according to quality factors in the given DSMC state file.'
  IF((.NOT.DoMacroscopicRestart).OR.(.NOT.DoRestart)) THEN
    CALL abort(__STAMP__,&
    'ERROR: It is required to use a restart and macroscopic restart when adapting the time step distribution!')
  END IF
  ! Open DSMC state file
  CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

  CALL GetDataProps('ElemData',nVar_HDF5,N_HDF5,nGlobalElems)

  ! Arrays might have been allocated if a time step was found in the state file
  IF(.NOT.ALLOCATED(VarTimeStep%ElemFac)) THEN
    ALLOCATE(VarTimeStep%ElemFac(nGlobalElems))
    VarTimeStep%ElemFac = 1.0
  END IF
#if USE_MPI
  IF(.NOT.ALLOCATED(VarTimeStep%ElemWeight)) THEN
    ALLOCATE(VarTimeStep%ElemWeight(nGlobalElems))
    VarTimeStep%ElemWeight = 1.0
  END IF
#endif

  IF(nVar_HDF5.LE.0) THEN
    SWRITE(*,*) 'ERROR: Something is wrong with our MacroscopicRestart file:', TRIM(MacroRestartFileName)
    CALL abort(__STAMP__,&
    'ERROR: Number of variables in the ElemData array appears to be zero!')
  END IF

  ! Get the variable names from the DSMC state and find the position of required quality factors
  ALLOCATE(VarNames_tmp(1:nVar_HDF5))
  CALL ReadAttribute(File_ID,'VarNamesAdd',nVar_HDF5,StrArray=VarNames_tmp(1:nVar_HDF5))

  DO iVar=1,nVar_HDF5
    IF (STRICMP(VarNames_tmp(iVar),"DSMC_MaxCollProb")) THEN
      nVar_MaxCollProb = iVar
    END IF
    IF (STRICMP(VarNames_tmp(iVar),"DSMC_MCS_over_MFP")) THEN
      nVar_MCSoverMFP = iVar
    END IF
    IF (STRICMP(VarNames_tmp(iVar),"Total_SimPartNum")) THEN
      nVar_TotalPartNum = iVar
    END IF
    ! Check if a time step distribution was written out in the DSMC state file
    IF (STRICMP(VarNames_tmp(iVar),"VariableTimeStep")) THEN
      nVar_TimeStep = iVar
    END IF
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
    IF (STRICMP(VarNames_tmp(iVar),"BGK_MaxRelaxationFactor").OR.STRICMP(VarNames_tmp(iVar),"FP_MaxRelaxationFactor")) THEN
      nVar_MaxRelaxFac = iVar
    END IF
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
  END DO

  ALLOCATE(ElemData_HDF5(1:nVar_HDF5,1:nGlobalElems))
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (nVar_HDF5    => INT(nVar_HDF5,IK) ,&
             nGlobalElems => INT(nGlobalElems,IK))
    CALL ReadArray('ElemData',2,(/nVar_HDF5,nGlobalElems/),0_IK,2,RealArray=ElemData_HDF5(:,:))
  END ASSOCIATE

  ALLOCATE(DSMCQualityFactors(nGlobalElems,1:2), PartNum(nGlobalElems))
  DSMCQualityFactors(:,1) = ElemData_HDF5(nVar_MaxCollProb,:)
  DSMCQualityFactors(:,2) = ElemData_HDF5(nVar_MCSoverMFP,:)
  PartNum(:)              = ElemData_HDF5(nVar_TotalPartNum,:)
  ! Check if a time step distribution is available in the DSMC state file and use that instead of the read-in from the state file
  IF(nVar_TimeStep.GT.0) VarTimeStep%ElemFac(:) = ElemData_HDF5(nVar_TimeStep,:)
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
  ALLOCATE(MaxRelaxFactor(nGlobalElems))
  MaxRelaxFactor(:) = ElemData_HDF5(nVar_MaxRelaxFac,:)
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
  DEALLOCATE(ElemData_HDF5)

  ! Calculating the time step per element based on the read-in max collision prob and mean collision separation
  DO iElem = 1, nGlobalElems
    TimeStepModified = .FALSE.
    ! Skipping cells, where less than 2 particles were sampled
    IF(PartNum(iElem).LT.2.0) CYCLE
#if USE_MPI
    ! Storing the old time step factor temporarily
    VarTimeStep%ElemWeight(iElem) = VarTimeStep%ElemFac(iElem)
#endif
    ! Storing either a 1 or the read-in time step factor in a temporary variable
    TimeFracTemp = VarTimeStep%ElemFac(iElem)
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
    ! Adapting the time step in order to achieve a maximal relaxation factor < 0.8
    IF(MaxRelaxFactor(iElem).GT.VarTimeStep%TargetMaxRelaxFactor) THEN
      TimeFracTemp = VarTimeStep%TargetMaxRelaxFactor*VarTimeStep%ElemFac(iElem) / MaxRelaxFactor(iElem)
      TimeStepModified = .TRUE.
    END IF
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
    ! Adapting the time step in order to achieve MaxCollProb < 0.8
    IF(DSMCQualityFactors(iElem,1).GT.VarTimeStep%TargetMaxCollProb) THEN
      TimeFracTemp = VarTimeStep%TargetMaxCollProb*VarTimeStep%ElemFac(iElem) / DSMCQualityFactors(iElem,1)
      TimeStepModified = .TRUE.
    END IF
    ! Adapting the time step in order to reduce the mean collision separation
    IF(DSMCQualityFactors(iElem,2).GT.VarTimeStep%TargetMCSoverMFP) THEN
      IF(Symmetry%Order.EQ.2) THEN
        TimeFracTemp = MIN(TimeFracTemp,VarTimeStep%ElemFac(iElem)*(VarTimeStep%TargetMCSoverMFP/DSMCQualityFactors(iElem,2))**2)
      ELSE
        TimeFracTemp = MIN(TimeFracTemp,VarTimeStep%ElemFac(iElem)*(VarTimeStep%TargetMCSoverMFP/DSMCQualityFactors(iElem,2))**3)
      END IF
      TimeStepModified = .TRUE.
    END IF
    ! Decrease time step according to a given minimal particle number (optional: default is 0)
    IF(VarTimeStep%DistributionMinPartNum.GT.0) THEN
      IF(PartNum(iElem).LT.VarTimeStep%DistributionMinPartNum) THEN
        TimeFracTemp = MIN(TimeFracTemp,PartNum(iElem)/VarTimeStep%DistributionMinPartNum*VarTimeStep%ElemFac(iElem))
        TimeStepModified = .TRUE.
      END IF
    END IF
    ! Limiting the minimal time step factor to the given value
    TimeFracTemp = MAX(TimeFracTemp,VarTimeStep%DistributionMinTimeFactor)
    ! If time step was not adapted due to particle number, collision probability or mean collision separation
    ! Choose appropriate time step to satisfy target MCSoverMFP, MaxCollProb and MinPartNum
    IF(.NOT.TimeStepModified) THEN
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
      IF(MaxRelaxFactor(iElem).GT.0.0) THEN
        TimeFracTemp = VarTimeStep%TargetMaxRelaxFactor*VarTimeStep%ElemFac(iElem) / MaxRelaxFactor(iElem)
      END IF
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
      IF(DSMCQualityFactors(iElem,1).GT.0.0) THEN
        TimeFracTemp = VarTimeStep%TargetMaxCollProb*VarTimeStep%ElemFac(iElem) / DSMCQualityFactors(iElem,1)
      END IF
      IF(DSMCQualityFactors(iElem,2).GT.0.0) THEN
        IF(Symmetry%Order.EQ.2) THEN
          TimeFracTemp = MIN(TimeFracTemp,VarTimeStep%ElemFac(iElem)*(VarTimeStep%TargetMCSoverMFP/DSMCQualityFactors(iElem,2))**2)
        ELSE
          TimeFracTemp = MIN(TimeFracTemp,VarTimeStep%ElemFac(iElem)*(VarTimeStep%TargetMCSoverMFP/DSMCQualityFactors(iElem,2))**3)
        END IF
      END IF
      IF(VarTimeStep%DistributionMinPartNum.GT.0) THEN
        TimeFracTemp = MIN(TimeFracTemp,PartNum(iElem)/VarTimeStep%DistributionMinPartNum*VarTimeStep%ElemFac(iElem))
      END IF
    END IF
    IF (VarTimeStep%OnlyDecreaseDt) THEN
      IF(TimeFracTemp.GT.VarTimeStep%ElemFac(iElem)) TimeFracTemp = VarTimeStep%ElemFac(iElem)
    END IF
    ! Finally, limiting the maximal time step factor to the given value and saving it to the right variable
    VarTimeStep%ElemFac(iElem) = MIN(TimeFracTemp,VarTimeStep%DistributionMaxTimeFactor)
#if USE_MPI
    ! Calculating the weight, multiplied with the particle number from state file during readMesh
    ! (covering the case when a time step distribution is read-in and adapted -> elements have been already once load-balanced with
    ! the old time step, consequently weight should only include difference between old and new time step)
    VarTimeStep%ElemWeight(iElem) = VarTimeStep%ElemWeight(iElem) / VarTimeStep%ElemFac(iElem)
#endif
  END DO
  ! Close the DSMC state file and deallocate not required variables
  CALL CloseDataFile()
  SDEALLOCATE(DSMCQualityFactors)
  SDEALLOCATE(PartNum)
#if (PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400)
  SDEALLOCATE(MaxRelaxFactor)
#endif /*PP_TimeDiscMethod==300 || PP_TimeDiscMethod==400*/
END IF      ! Adapt Distribution

SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE VarTimeStep_InitDistribution


REAL FUNCTION GetParticleTimeStep(xPos, yPos, iElem)
!===================================================================================================================================
!> Calculates/determines the time step
!> a) at a position x/y (only in 2D/Axi) [VarTimeStep%UseLinearScaling]
!> b) of the given element number (3D and VTS distribution) [VarTimeStep%UseDistribution]
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: VarTimeStep, Symmetry
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN), OPTIONAL      :: xPos, yPos
INTEGER, INTENT(IN), OPTIONAL   :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL          :: xFactor
!===================================================================================================================================

GetParticleTimeStep = 1.

IF(VarTimeStep%UseLinearScaling) THEN
  IF(Symmetry%Order.EQ.2) THEN
    IF (VarTimeStep%Use2DTimeFunc) THEN
      IF(.NOT.PRESENT(xPos).OR..NOT.PRESENT(yPos)) CALL abort(__STAMP__,&
        'ERROR: Position in x-direction is required in the call of GetParticleTimeStep for linear scaling in 2D!')
      IF (xPos.LT.VarTimeStep%StagnationPoint) THEN
        xFactor = ABS((VarTimeStep%StagnationPoint-xPos)/(VarTimeStep%StagnationPoint-GEO%xminglob) &
                      * (VarTimeStep%TimeScaleFac2DFront - 1.0))
      ELSE
        xFactor = ABS((xPos-VarTimeStep%StagnationPoint)/(GEO%xmaxglob-VarTimeStep%StagnationPoint) &
                      * (VarTimeStep%TimeScaleFac2DBack - 1.0))
      END IF
      GetParticleTimeStep = (1. + yPos/GEO%ymaxglob*(VarTimeStep%ScaleFac-1.0))*(1.+xFactor)
    ELSE
      IF(.NOT.PRESENT(yPos)) CALL abort(__STAMP__,&
        'ERROR: Position in x-direction is required in the call of GetParticleTimeStep for linear scaling in 2D!')
      GetParticleTimeStep = (1. + yPos/GEO%ymaxglob*(VarTimeStep%ScaleFac-1.0))
    END IF
  ELSE
    IF(.NOT.PRESENT(iElem)) CALL abort(__STAMP__,&
      'ERROR: Element number is required in the call of GetParticleTimeStep for distribution/scaling in 3D!')
    GetParticleTimeStep = VarTimeStep%ElemFac(iElem)
  END IF
ELSE IF(VarTimeStep%UseDistribution) THEN
  IF(.NOT.PRESENT(iElem)) CALL abort(__STAMP__,&
    'ERROR: Element number is required in the call of GetParticleTimeStep for distribution!')
  GetParticleTimeStep = VarTimeStep%ElemFac(iElem)
ELSE
  CALL abort(__STAMP__,'ERROR: GetParticleTimeStep should not be utilized without LinearScaling/Distribution flag!')
END IF

RETURN

END FUNCTION GetParticleTimeStep


PURE REAL FUNCTION GetSpeciesTimeStep(iCase)
!===================================================================================================================================
!> Determines the species-specific time step from the collision case
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: VarTimeStep, Species
USE MOD_DSMC_Vars               ,ONLY: CollInf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iCase
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER               :: iSpec, jSpec, SpecID
!===================================================================================================================================

GetSpeciesTimeStep = 1.

! ManualTimeStep has been utilized for the collisions
IF(VarTimeStep%DisableForMCC) RETURN

! Determine the particle species (one will be the background species with a TimeStepFactor of 1)
iSpec = CollInf%collidingSpecies(iCase,1)
jSpec = CollInf%collidingSpecies(iCase,2)
IF(Species(iSpec)%TimeStepFactor.LT.1.) THEN
  SpecID = iSpec
ELSE IF(Species(jSpec)%TimeStepFactor.LT.1.) THEN
  SpecID = jSpec
END IF
GetSpeciesTimeStep = Species(SpecID)%TimeStepFactor

RETURN

END FUNCTION GetSpeciesTimeStep


SUBROUTINE VarTimeStep_CalcElemFacs()
!===================================================================================================================================
!> Calculates/determines the variable time step element-wise (every proc for his own domain) during the initialization
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Vars          ,ONLY: VarTimeStep
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemMidPoint_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                   :: iElem,CNElemID
!===================================================================================================================================

ALLOCATE(VarTimeStep%ElemFac(nElems))
VarTimeStep%ElemFac = 1.0
IF (VarTimeStep%Direction(1).GT.0.0) THEN
  DO iElem = 1, nElems
    CNElemID = GetCNElemID(iElem + offsetElem)
    IF (ElemMidPoint_Shared(1,CNElemID).LT.VarTimeStep%StartPoint(1)) THEN
      VarTimeStep%ElemFac(iElem)=1.0
    ELSE IF (VarTimeStep%EndPoint(1).EQ.-99999.) THEN
      VarTimeStep%ElemFac(iElem)= 1.0 + (VarTimeStep%ScaleFac-1.0)/(GEO%xmaxglob-VarTimeStep%StartPoint(1)) &
          * (ElemMidPoint_Shared(1,CNElemID)-VarTimeStep%StartPoint(1))
    ELSE
      IF (ElemMidPoint_Shared(1,CNElemID).GT.VarTimeStep%EndPoint(1)) THEN
        VarTimeStep%ElemFac(iElem)=VarTimeStep%ScaleFac
      ELSE
        VarTimeStep%ElemFac(iElem)= 1.0 + (VarTimeStep%ScaleFac-1.0)/(VarTimeStep%EndPoint(1)-VarTimeStep%StartPoint(1)) &
            * (ElemMidPoint_Shared(1,CNElemID)-VarTimeStep%StartPoint(1))
      END IF
    END IF
  END DO
ELSE
  DO iElem = 1, nElems
    CNElemID = GetCNElemID(iElem + offsetElem)
    IF (ElemMidPoint_Shared(1,CNElemID).GT.VarTimeStep%StartPoint(1)) THEN
      VarTimeStep%ElemFac(iElem)=1.0
    ELSE IF (VarTimeStep%EndPoint(1).EQ.-99999.) THEN
      VarTimeStep%ElemFac(iElem)= 1.0 + (VarTimeStep%ScaleFac-1.0)/(VarTimeStep%StartPoint(1)-GEO%xminglob) &
          * (VarTimeStep%StartPoint(1)-ElemMidPoint_Shared(1,CNElemID))
    ELSE
      IF (ElemMidPoint_Shared(1,CNElemID).LT.VarTimeStep%EndPoint(1)) THEN
        VarTimeStep%ElemFac(iElem)=VarTimeStep%ScaleFac
      ELSE
        VarTimeStep%ElemFac(iElem)= 1.0 + (VarTimeStep%ScaleFac-1.0)/(VarTimeStep%StartPoint(1)-VarTimeStep%EndPoint(1)) &
            * (VarTimeStep%StartPoint(1)-ElemMidPoint_Shared(1,CNElemID))
      END IF
    END IF
  END DO
END IF

END SUBROUTINE VarTimeStep_CalcElemFacs

END MODULE MOD_Particle_TimeStep
