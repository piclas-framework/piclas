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

CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , "TODO-DEFINE-PARAMETER Type of Interpolation-Method to calculat"//&
								"e the EM field's value for the particle", 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'TODO-DEFINE-PARAMETER Interpolate with outer iElem-loop (not'//&
								'for many Elems per proc!)', '.TRUE.')
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'TODO-DEFINE-PARAMETER External field is added to the'//&
								'maxwell-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'TODO-DEFINE-PARAMETER', '1.0')
CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , "TODO-DEFINE-PARAMETER Compute the self field's influence "//&
								"on the Particle", '.TRUE.')
CALL prms%CreateLogicalOption(  'PIC-BG-Field'                , 'TODO-DEFINE-PARAMETER BGField data (1:x,0:NBG,0:NBG,0:NBG,'//&
								'1:PP_nElems) DEFAULT=F If PIC-BG-Field=T, define: '//&
								'PIC-BGFilename,PIC-BGFieldScaling,PIC-NBG', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-BGFileName'              , 'TODO-DEFINE-PARAMETER File name for background field '//&
								'([character].h5)', 'none')
CALL prms%CreateIntOption(      'PIC-NBG'                     , 'TODO-DEFINE-PARAMETER Polynomial degree that shall be used '//&
								'for background field during simulation', '1')
CALL prms%CreateRealOption(     'PIC-BGFieldScaling'          , 'TODO-DEFINE-PARAMETER Space scaling of background field', '1.')
CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'TODO-DEFINE-PARAMETER File to curved external field data.','none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'TODO-DEFINE-PARAMETER File containing the external '//&
								'field CSV table', 'none')

CALL prms%CreateIntOption(      'PIC-nCollectChargesBCs'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntOption(      'PIC-CollectCharges[$]-BC'    , 'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-NumOfRealCharges'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-ChargeDist'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('PIC-NormVecOfWall'  , 'TODO-DEFINE-PARAMETER Normal vector for pushTimeStep', '1. , 0. , 0.')
CALL prms%CreateIntOption(      'PIC-DeltaType'      , 'TODO-DEFINE-PARAMETER Flag ', '1')
CALL prms%CreateIntOption(      'PIC-DeltaType-N'    , 'TODO-DEFINE-PARAMETER =PP_N Polynomial degree of delta distribution', '1')
CALL prms%CreateRealArrayOption('PIC-BGMdeltas'      , 'TODO-DEFINE-PARAMETER Dimensions of PIC background mesh', '0. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-FactorBGM'      , 'TODO-DEFINE-PARAMETER Denominator of PIC-BGMdeltas', '1. , 1. , 1.')
CALL prms%CreateLogicalOption(  'PIC-OutputSource'   , 'TODO-DEFINE-PARAMETER Writes the source to hdf5', '.FALSE.')

CALL prms%SetSection("PIC Deposition")

CALL prms%CreateLogicalOption(  'PIC-DoDeposition'         , 'TODO-DEFINE-PARAMETER Switch deposition on/off', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Deposition-Type'      , 'TODO-DEFINE-PARAMETER If Deposition-Type=shape_function, define: '//&
							     'PIC-shapefunction-radius,PIC-shapefunction-alpha.\n(HALOWIKI:) If'//&
							     'Deposition-Type =(cartmesh_volumeweighting/ cartmesh_splines),\n'//&
							     'Define: PIC-BGMdeltas,PIC-FactorBGM', 'nearest-blurrycenter')
CALL prms%CreateStringOption(   'PIC-TimeAverageFile'      , 'TODO-DEFINE-PARAMETER', 'none')

CALL prms%CreateRealOption(     'PIC-epanechnikov-radius'  , 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius' , 'TODO-DEFINE-PARAMETER Radius of shape function', '1.')
CALL prms%CreateIntOption(      'PIC-shapefunction-alpha'  , 'TODO-DEFINE-PARAMETER Exponent of shape function', '2')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-equi'   , 'TODO-DEFINE-PARAMETER Use equidistant points for shapefunction'&
							   , '.FALSE.')
CALL prms%CreateIntOption(      'PIC-shapefunction1d-direction'  , 'TODO-DEFINE-PARAMETER Direction of 1D shape function', '1')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'TODO-DEFINE-PARAMETER Minimal shape function radius', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'TODO-DEFINE-PARAMETER Scaling factor of shape function radius', '0.')
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixes'     , 'TODO-DEFINE-PARAMETER Number of fixes for shape func depo at'//&
							     ' planar BCs', '0')
CALL prms%CreateLogicalOption(  'PrintSFDepoWarnings'      , 'TODO-DEFINE-PARAMETER Print the shapefunction warnings', '.FALSE.')
CALL prms%CreateRealOption(     'PIC-SFdepoFixesEps'       , 'TODO-DEFINE-PARAMETER Epsilon for defined planes', '0.')
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Normal','TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ChargeMult'  , 'TODO-DEFINE-PARAMETER Multiplier for mirrored charges '//&
						 	     '(wall: -1.0, sym: 1.0)', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmin'  , 'TODO-DEFINE-PARAMETER -> SFdepoFixesBounds(:,:,:)     '//&
							     '1:nFixes;1:2(min,max);1:3(x,y,z)?', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixLinks'  , 'TODO-DEFINE-PARAMETER Number of linked SFdepoFixes ', '0')
CALL prms%CreateIntArrayOption( 'PIC-SFdepoFixLink[$]'     , 'TODO-DEFINE-PARAMETER 1:nLinks;1:3 (2 fixes are linked with each'//&
							     ' other!) (:,3 is fraction of 180 deg)', '1 , 2'&
							   , numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoLayers'  ,    'TODO-DEFINE-PARAMETER Number of const. source layer for sf-depo'//&
							      ' at planar BCs', '0')
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Normal', 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmin'  , 'TODO-DEFINE-PARAMETER -> SFdepoLayersBounds(:,:,:)     '//&
							      '1:nFixes;1:2(min,max);1:3(x,y,z)?', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'PIC-SFdepoLayers[$]-UseFixBounds'   , 'TODO-DEFINE-PARAMETER Use alls planes of SFdepoFixes as'//&
							     	       ' additional bounds', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'PIC-SFdepoLayers[$]-Space' ,          'TODO-DEFINE-PARAMETER Name of space (cuboid or cylinder)'&
							    ,          'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector1'  ,   'TODO-DEFINE-PARAMETER Base Vector 1'&
							    ,          '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector2'  ,   'TODO-DEFINE-PARAMETER Base Vector 2', '0. , 0. , 1.'&
							    	   ,    numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-SFdepoLayersRadius', 'TODO-DEFINE-PARAMETER  Radius for cylinder-space'&
									, '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-Chargedens'        , 'TODO-DEFINE-PARAMETER','1.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-SFdepoLayers[$]-Spec'        ,    'TODO-DEFINE-PARAMETER Particle species for respective'//&
								       ' layer','1', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'PIC-SFResampleAnalyzeSurfCollis'  , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntArrayOption( 'PIC-SFResampleSurfCollisBC',        'TODO-DEFINE-PARAMETER BCs to be analyzed (def.: 0 = all)')
CALL prms%CreateLogicalOption(  'PIC-SFResampleReducePartNumber'   , 'TODO-DEFINE-PARAMETER Reduce PartNumberSamp to '//&
								     'PartNumberReduced', '.FALSE.')
CALL prms%CreateIntOption(      'PIC-PartNumThreshold'      , 	     'TODO-DEFINE-PARAMETER Threshold for checking inserted '//&
							      	     'parts per deposition (otherwise abort)', '0')
CALL prms%CreateIntOption(      'PIC-SFResampleNumberOfBCs' ,        'TODO-DEFINE-PARAMETER Number of BC to be analyzed', '1')
CALL prms%CreateIntOption(      'PIC-SFResamplePartNumberReduced'  , 'TODO-DEFINE-PARAMETER Max. allowed number of parts to'//&
								     ' be saved', '0')
CALL prms%CreateLogicalOption(  'PIC-SFResampleRestart'     , 'TODO-DEFINE-PARAMETER Read-in old DSMCSurfCollis-file for restart'&
							    , '.FALSE.')
CALL prms%CreateStringOption(   'PIC-SFResampleRestartFile' , 'TODO-DEFINE-PARAMETER Name of the new DSMCSurfCollis-file to'//&
							      ' read-in by restart', 'dummy')
CALL prms%CreateRealOption(     'PIC-SFResample-xmin'       , 'TODO-DEFINE-PARAMETER ')
CALL prms%CreateRealOption(     'PIC-SFResample-ymin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-xmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-ymax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateLogicalOption(  'PIC-SFResample-UseFixBounds','TODO-DEFINE-PARAMETER Use all planes of SFdepoFixes as '//&
							      'additional bounds?', '.TRUE.')

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
