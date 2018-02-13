#include "boltzplatz.h"

MODULE MOD_PICInit
!===================================================================================================================================
! Includes PIC Init
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
INTERFACE InitPIC
  MODULE PROCEDURE InitPIC
END INTERFACE
PUBLIC::InitPIC
!===================================================================================================================================
PUBLIC::DefineParametersPIC
CONTAINS

!==================================================================================================================================
!> Define parameters for PIC
!==================================================================================================================================
SUBROUTINE DefineParametersPIC()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("PIC")

CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , 'TODO-DEFINE-PARAMETER', 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'TODO-DEFINE-PARAMETER', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'TODO-DEFINE-PARAMETER', '1.0')
CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateLogicalOption(  'PIC-BG-Field'                , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-BGFileName'              , 'TODO-DEFINE-PARAMETER', 'none')
CALL prms%CreateIntOption(      'PIC-NBG'                     , 'TODO-DEFINE-PARAMETER', '1')
CALL prms%CreateRealOption(     'PIC-BGFieldScaling'          , 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'TODO-DEFINE-PARAMETER', 'none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'TODO-DEFINE-PARAMETER', 'none')

CALL prms%CreateIntOption(      'PIC-nCollectChargesBCs'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntOption(      'PIC-CollectCharges[$]-BC'    , 'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-NumOfRealCharges'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-ChargeDist'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('PIC-NormVecOfWall'  , 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.')
CALL prms%CreateIntOption(      'PIC-DeltaType'      , 'TODO-DEFINE-PARAMETER', '1')
CALL prms%CreateIntOption(      'PIC-DeltaType-N'    , 'TODO-DEFINE-PARAMETER =PP_N', '1')
CALL prms%CreateRealArrayOption('PIC-BGMdeltas'      , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-FactorBGM'      , 'TODO-DEFINE-PARAMETER', '1. , 1. , 1.')
CALL prms%CreateLogicalOption(  'PIC-OutputSource'   , 'TODO-DEFINE-PARAMETER', '.FALSE.')

CALL prms%SetSection("PIC Deposition")

CALL prms%CreateLogicalOption(  'PIC-DoDeposition'         , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Deposition-Type'      , 'TODO-DEFINE-PARAMETER', 'nearest-blurrycenter')

CALL prms%CreateRealOption(     'PIC-epanechnikov-radius'  , 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius' , 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateIntOption(      'PIC-shapefunction-alpha'  , 'TODO-DEFINE-PARAMETER', '2')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-equi'   , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntOption(      'PIC-shapefunction1d-direction'  , 'TODO-DEFINE-PARAMETER', '1')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'TODO-DEFINE-PARAMETER', '0.')
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixes'     , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateLogicalOption(  'PrintSFDepoWarnings'      , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateRealOption(     'PIC-SFdepoFixesEps'       , 'TODO-DEFINE-PARAMETER', '0.')
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Normal'  , 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ChargeMult'  , 'TODO-DEFINE-PARAMETER', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixLinks'  , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntArrayOption( 'PIC-SFdepoFixLink[$]'     , 'TODO-DEFINE-PARAMETER', '1 , 2', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoLayers'  , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Normal'  , 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'PIC-SFdepoLayers[$]-UseFixBounds'   , 'TODO-DEFINE-PARAMETER', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'PIC-SFdepoLayers[$]-Space'  , 'TODO-DEFINE-PARAMETER', 'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector1'  , 'TODO-DEFINE-PARAMETER', '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector2'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-SFdepoLayersRadius', 'TODO-DEFINE-PARAMETER','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-Chargedens'        , 'TODO-DEFINE-PARAMETER','1.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-SFdepoLayers[$]-Spec'        , 'TODO-DEFINE-PARAMETER','1', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'PIC-SFResampleAnalyzeSurfCollis'  , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntArrayOption( 'PIC-SFResampleSurfCollisBC', 'TODO-DEFINE-PARAMETER')
CALL prms%CreateLogicalOption(  'PIC-SFResampleReducePartNumber'   , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntOption(      'PIC-PartNumThreshold'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntOption(      'PIC-SFResampleNumberOfBCs' , 'TODO-DEFINE-PARAMETER', '1')
CALL prms%CreateIntOption(      'PIC-SFResamplePartNumberReduced'  , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateLogicalOption(  'PIC-SFResampleRestart'     , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateStringOption(   'PIC-SFResampleRestartFile' , 'TODO-DEFINE-PARAMETER', 'dummy')
CALL prms%CreateRealOption(     'PIC-SFResample-xmin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-ymin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-xmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-ymax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateLogicalOption(  'PIC-SFResample-UseFixBounds'      , 'TODO-DEFINE-PARAMETER', '.TRUE.')

END SUBROUTINE DefineParametersPIC

SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars,  ONLY: externalField
USE MOD_PIC_Vars ,              ONLY: PICInitIsDone, PIC
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(PICInitIsDone)THEN
   SWRITE(*,*) "InitPIC already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PIC ...'

! So far, nothing to do here...
IF (externalField(6).NE.0) PIC%GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
