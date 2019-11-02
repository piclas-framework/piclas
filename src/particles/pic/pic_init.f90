!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

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

#ifdef CODE_ANALYZE
CALL prms%CreateLogicalOption('PIC-DoInterpolationAnalytic'      , "Use an analytic/algebraic function for PIC interpolation "//&
                                                                   "(ifdef CODE_ANALYZE)",&
                                                                   '.FALSE.')

CALL prms%CreateIntOption(    'PIC-AnalyticInterpolation-Type'   , "Type of AnalyticInterpolation-Method for calculating the "//&
                                                                   "EM field's value for the particle (ifdef CODE_ANALYZE)",'0')

CALL prms%CreateIntOption(    'PIC-AnalyticInterpolation-SubType', "SubType of AnalyticInterpolation-Method for calculating the "//&
                                                                   "EM field's value for the particle (ifdef CODE_ANALYZE)",'0')

CALL prms%CreateRealOption(   'PIC-AnalyticInterpolationP'       , "parameter 'p' for AnalyticInterpolationType = 1", '1.')
#endif /*CODE_ANALYZE*/

CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , "TODO-DEFINE-PARAMETER\n"//&
                                                                "Type of Interpolation-Method to calculate"//&
                                                                " the EM field's value for the particle", 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Interpolate with outer iElem-loop (not'//&
                                                                'for many Elems per proc!)', '.TRUE.')
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'External field is added to the'//&
                                                                'maxwell-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'TODO-DEFINE-PARAMETER', '1.0')
CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , "TODO-DEFINE-PARAMETER\n"//&
                                                                "Compute the self field's influence "//&
                                                                "on the Particle", '.TRUE.')
CALL prms%CreateLogicalOption(  'PIC-BG-Field'                , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'BGField data (1:x,0:NBG,0:NBG,0:NBG,'//&
                                                                '1:PP_nElems) \n'//&
                                                                'If PIC-BG-Field=T\n'//&
                                                                'Define:\n'//&
                                                                'PIC-BGFilename\n'//&
                                                                'PIC-BGFieldScaling\n'//&
                                                                'PIC-NBG', '.FALSE.')
CALL prms%CreateLogicalOption(  'PIC-CalcBField'              , 'Calculate the BGField from parameters supplied by user', '.FALSE.')
CALL prms%CreateStringOption(   'PIC-BGFileName'              , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'File name for background field ([character].h5)', 'none')
CALL prms%CreateIntOption(      'PIC-NBG'                     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Polynomial degree that shall be used '//&
                                                                'for background field during simulation', '1')
CALL prms%CreateRealOption(     'PIC-BGFieldScaling'          , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Space scaling of background field', '1.')
CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'File to curved external field data.','none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'File containing the external '//&
                                                                'field CSV table', 'none')

CALL prms%CreateIntOption(      'PIC-nCollectChargesBCs'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntOption(      'PIC-CollectCharges[$]-BC'    , 'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-NumOfRealCharges'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-ChargeDist'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)

CALL prms%CreateRealArrayOption('PIC-NormVecOfWall'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Normal vector for pushTimeStep', '1. , 0. , 0.')
CALL prms%CreateIntOption(      'PIC-DeltaType'      , 'Basis function type.\n'//&
                                                       '1: Lagrange-Polynomial\n'//&
                                                       '2: Bernstein-Polynomial\n'//&
                                                       '3: Uniform B-Spline', '1')
CALL prms%CreateIntOption(      'PIC-DeltaType-N'    , 'Polynomial degree of the delta distribution basis function', '1')
CALL prms%CreateRealArrayOption('PIC-BGMdeltas'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Dimensions of PIC background mesh', '0. , 0. , 0.')
CALL prms%CreateRealArrayOption('PIC-FactorBGM'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Denominator of PIC-BGMdeltas', '1. , 1. , 1.')
CALL prms%CreateLogicalOption(  'PIC-OutputSource'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                       'Writes the source to hdf5', '.FALSE.')

CALL prms%SetSection("PIC Deposition")
CALL prms%CreateLogicalOption(  'PIC-DoDeposition'         , 'Switch deposition of charge (and current density) on/off', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Deposition-Type'      , '1.1)  shape_function\n'                   //&
                                                             '1.2)  shape_function_1d\n'                //&
                                                             '1.3)  shape_function_2d\n'                //&
                                                             '1.4)  shape_function_cylindrical\n'       //&
                                                             '1.5)  shape_function_spherical\n'         //&
                                                             '1.6)  shape_function_simple\n'            //&
                                                             '      1.1) to 1.6) require\n'            //&
                                                             '        PIC-shapefunction-radius\n'//&
                                                             '        PIC-shapefunction-alpha\n' //&
                                                             '      1.2) and 1.3) require\n'            //&
                                                             '        PIC-shapefunction1d-direction\n'  //&
                                                             '      1.4) and 1.5) require\n'            //&
                                                             '        PIC-shapefunction-radius0\n'      //&
                                                             '        PIC-shapefunction-scale\n'        //&
                                                             '2.)   cell_volweight\n'                   //&
                                                             '3.)   epanechnikov\n'                     //&
                                                             '4.)   nearest_gausspoint\n'               //&
                                                             '5.)   delta_distri\n'                     //&
                                                             '      requires PIC-DeltaType\n'           //&
                                                             '               PIC-DeltaType-N\n'         //&
                                                             '6.1)  cartmesh_volumeweighting\n'         //&
                                                             '6.2)  cartmesh_splines\n'                 //&
                                                             '      requires PIC-BGMdeltas\n'           //&
                                                             '               PIC-FactorBGM\n'           //&
                                                             '7.)   nearest-blurrycenter\n'             //&
                                                             '8.)   cell_volweight_mean'                &
                                                           , 'nearest-blurrycenter') ! Default
CALL prms%CreateStringOption(   'PIC-TimeAverageFile'      , 'TODO-DEFINE-PARAMETER', 'none')
CALL prms%CreateLogicalOption(  'PIC-RelaxDeposition'      , 'Relaxation of current PartSource with RelaxFac\n'//&
                                                             'into PartSourceOld', '.FALSE.')
CALL prms%CreateRealOption(     'PIC-RelaxFac'             , 'Relaxation factor of current PartSource with RelaxFac\n'//&
                                                             'into PartSourceOld', '0.001')

CALL prms%CreateRealOption(     'PIC-epanechnikov-radius'  , 'TODO-DEFINE-PARAMETER', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius' , 'Radius of shape function', '1.')
CALL prms%CreateIntOption(      'PIC-shapefunction-alpha'  , 'Exponent of shape function', '2')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-equi'   , 'Use equidistant points for shapefunction deposition' , '.FALSE.')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-local-depo-BC', 'Do not use shape function deposition in elements where a '//&
                                                                   'boundary would truncate the shape function. Use a local '//&
                                                                   'deposition in these elements instead of the shape function.'&
                                                                 , '.FALSE.')
CALL prms%CreateIntOption(      'PIC-shapefunction1d-direction' ,'1D shape function: Deposition direction\n'//&
                                                                 '2D shape function: Perpendicular deposition')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-3D-deposition' ,'Deposite the charge over volume (3D)\n'//&
                                                                   ' or over a line (1D)/area(2D)\n'//&
                                                                   '1D shape function: volume or line\n'//&
                                                                   '2D shape function: volume or area', '.TRUE.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'Minimum shape function radius (for cylindrical and spherical)', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'Scaling factor of shape function radius '//&
                                                             '(for cylindrical and spherical)', '0.')
! Shape Function Deposition Fixes
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixes'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Number of fixes for shape func depo at'//&
                                                             ' planar BCs', '0')
CALL prms%CreateLogicalOption(  'PrintSFDepoWarnings'      , 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Print the shapefunction warnings', '.FALSE.')
CALL prms%CreateRealOption(     'PIC-SFdepoFixesEps'       , 'TODO-DEFINE-PARAMETER\n'//&
                                                             'Epsilon for defined planes', '0.')
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER\n', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoFixes[$]-Normal','TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ChargeMult'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                    'Multiplier for mirrored charges '//&
                                                              '(wall: -1.0, sym: 1.0)', '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmin'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                    '-> SFdepoFixesBounds(:,:,:) 1:nFixes;1:2(min,max);1:3(x,y,z)?'&
                                                           , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoFixes[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoFixLinks'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                    'Number of linked SFdepoFixes ', '0')
CALL prms%CreateIntArrayOption( 'PIC-SFdepoFixLink[$]'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                    '(1:nLinks)\n'//&
                                                             '1:3 (2 fixes are linked with each other!)\n'//&
                                                             ':,3 is fraction of 180 deg', '1 , 2'&
                                                           , numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'PIC-NbrOfSFdepoLayers'  ,    'TODO-DEFINE-PARAMETER\n'//&
                                                                    'Number of const. source layer for sf-depo'//&
                                                              ' at planar BCs', '0')
CALL prms%CreateLogicalOption(  'PIC-ConstantSFdepoLayers'      , 'Do deposition of SFdepoLayers just once', '.FALSE.')
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Basepoint'  , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-Normal', 'TODO-DEFINE-PARAMETER', '1. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmin'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                     '-> SFdepoLayersBounds(:,:,:)\n'//&
                                                              '1:nFixes;1:2(min,max);1:3(x,y,z)?', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmin'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-xmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-ymax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-zmax'  , 'TODO-DEFINE-PARAMETER', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'PIC-SFdepoLayers[$]-UseFixBounds'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                     'Use alls planes of SFdepoFixes as'//&
                                                              ' additional bounds', '.TRUE.', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'PIC-SFdepoLayers[$]-Space' ,          'TODO-DEFINE-PARAMETER\n'//&
                                                                                   'Name of space (cuboid or cylinder)'&
                                                            ,          'cuboid', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector1'  ,   'TODO-DEFINE-PARAMETER\n'//&
                                                                                   'Base Vector 1'&
                                                            ,          '0. , 1. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('PIC-SFdepoLayers[$]-BaseVector2'  ,   'TODO-DEFINE-PARAMETER\n'//&
                                                                                   'Base Vector 2', '0. , 0. , 1.'&
                                                                       ,    numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-SFdepoLayersRadius', 'TODO-DEFINE-PARAMETER\n'//&
                                                                                      'Radius for cylinder-space'&
                                                                        , '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-Chargedens'        , 'TODO-DEFINE-PARAMETER','1.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'PIC-SFdepoLayers[$]-Spec'        ,    'TODO-DEFINE-PARAMETER\n'//&
                                                                                   'Particle species for respective'//&
                                                                       ' layer','1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-SFdepoLayers[$]-MPF'        ,    'MPF for respective'//&
                                                                       ' layer (def.: MPF of resp. species)', numberedmulti=.TRUE.)

CALL prms%CreateLogicalOption(  'PIC-SFResampleAnalyzeSurfCollis'  , 'TODO-DEFINE-PARAMETER', '.FALSE.')
CALL prms%CreateIntArrayOption( 'PIC-SFResampleSurfCollisBC',        'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'BCs to be analyzed (def.: 0 = all)')
CALL prms%CreateLogicalOption(  'PIC-SFResampleReducePartNumber'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                     'Reduce PartNumberSamp to PartNumberReduced', '.FALSE.')
CALL prms%CreateIntOption(      'PIC-PartNumThreshold'      ,              'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'Threshold for checking inserted '//&
                                                                           'parts per deposition (otherwise abort)', '0')
CALL prms%CreateIntOption(      'PIC-SFResampleNumberOfBCs' ,        'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'Number of BC to be analyzed', '1')
CALL prms%CreateIntOption(      'PIC-SFResamplePartNumberReduced'  , 'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'Max. allowed number of parts to be saved', '0')
CALL prms%CreateIntOption(      'PIC-SFResampleNbrOfSpeciesForDtCalc' ,        'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'Number of species used for SFResample-dt', '1')
CALL prms%CreateIntArrayOption( 'PIC-SFResampleSpeciesForDtCalc',        'TODO-DEFINE-PARAMETER\n'//&
                                                                                 'Species used for SFResample-dt (def.: 0 = all)')
CALL prms%CreateLogicalOption(  'PIC-SFResampleRestart'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                     'Read-in old DSMCSurfCollis-file for restart'&
                                                            , '.FALSE.')
CALL prms%CreateStringOption(   'PIC-SFResampleRestartFile' , 'TODO-DEFINE-PARAMETER\n'//&
                                                                     'Name of the new DSMCSurfCollis-file to'//&
                                                              ' read-in by restart', 'dummy')
CALL prms%CreateRealOption(     'PIC-SFResample-xmin'       , 'TODO-DEFINE-PARAMETER ')
CALL prms%CreateRealOption(     'PIC-SFResample-ymin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmin'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-xmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-ymax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateRealOption(     'PIC-SFResample-zmax'       , 'TODO-DEFINE-PARAMETER')
CALL prms%CreateLogicalOption(  'PIC-SFResample-UseFixBounds','TODO-DEFINE-PARAMETER\n'//&
                                                                     'Use all planes of SFdepoFixes as '//&
                                                              'additional bounds?', '.TRUE.')

END SUBROUTINE DefineParametersPIC

SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars,  ONLY: externalField
USE MOD_PICInterpolation       ,ONLY: InitializeParticleInterpolation
USE MOD_PICDepo                ,ONLY: InitializeDeposition
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

CALL InitializeParticleInterpolation()
CALL InitializeDeposition()

! So far, nothing to do here...
IF (externalField(6).NE.0) PIC%GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
