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
USE MOD_ReadInTools    ,ONLY: prms
USE MOD_PICDepo_Method ,ONLY: DefineParametersDepositionMethod
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
CALL prms%CreateRealOption(   'PIC-AnalyticInterpolationPhase'   , "Phase shift angle phi that is used for cos(w*t + phi)", '0.')
#endif /*CODE_ANALYZE*/

CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , "TODO-DEFINE-PARAMETER\n"//&
                                                                "Compute the self field's influence "//&
                                                                "on the Particle", '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Interpolation-Type'      , "Type of Interpolation-Method to calculate the electro(-magnetic)"//&
                                                                " field's value for the particle", 'particle_position')
CALL prms%CreateLogicalOption(  'PIC-InterpolationElemLoop'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'Interpolate with outer iElem-loop (not'//&
                                                                'for many Elems per proc!)', '.TRUE.')
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'External field is added to the'//&
                                                                'maxwell-solver-field', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'TODO-DEFINE-PARAMETER', '1.0')

CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'File to curved external field data.','none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'TODO-DEFINE-PARAMETER\n'//&
                                                                'File containing the external '//&
                                                                'field CSV table', 'none')

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
CALL DefineParametersDepositionMethod() ! Get PIC-DoDeposition and PIC-Deposition-Type
CALL prms%CreateStringOption(   'PIC-TimeAverageFile'      , 'TODO-DEFINE-PARAMETER', 'none')
CALL prms%CreateLogicalOption(  'PIC-RelaxDeposition'      , 'Relaxation of current PartSource with RelaxFac\n'//&
                                                             'into PartSourceOld', '.FALSE.')
CALL prms%CreateRealOption(     'PIC-RelaxFac'             , 'Relaxation factor of current PartSource with RelaxFac\n'//&
                                                             'into PartSourceOld', '0.001')

CALL prms%CreateLogicalOption(  'PIC-shapefunction-charge-conservation'   , 'Enable charge conservation.' , '.FALSE.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius' , 'Radius of shape function', '1.')
CALL prms%CreateIntOption(      'PIC-shapefunction-alpha'  , 'Exponent of shape function', '2')
CALL prms%CreateIntOption(      'PIC-shapefunction-dimension', '1D, 2D or 3D shape function', '3')
CALL prms%CreateIntOption(      'PIC-shapefunction-direction',&
    'Only required for PIC-shapefunction-dimension 1 or 2: Shape function direction for 1D (the direction in which the charge '//&
    'will be distributed) and 2D (the direction in which the charge will be constant)', '1')
!CALL prms%CreateLogicalOption(  'PIC-shapefunction-equi'   , 'Use equidistant points for shapefunction deposition' , '.FALSE.')
!CALL prms%CreateLogicalOption(  'PIC-shapefunction-local-depo-BC', 'Do not use shape function deposition in elements where a '//&
                                                                   !'boundary would truncate the shape function. Use a local '//&
                                                                   !'deposition in these elements instead of the shape function.'&
                                                                 !, '.FALSE.')
!CALL prms%CreateIntOption(      'PIC-shapefunction1d-direction' ,'1D shape function: Deposition direction\n'//&
                                                                 !'2D shape function: Perpendicular deposition')
CALL prms%CreateLogicalOption(  'PIC-shapefunction-3D-deposition' ,'Deposite the charge over volume (3D)\n'//&
                                                                   ' or over a line (1D)/area(2D)\n'//&
                                                                   '1D shape function: volume or line\n'//&
                                                                   '2D shape function: volume or area', '.TRUE.')
CALL prms%CreateRealOption(     'PIC-shapefunction-radius0', 'Minimum shape function radius (for cylindrical and spherical)', '1.')
CALL prms%CreateRealOption(     'PIC-shapefunction-scale'  , 'Scaling factor of shape function radius '//&
                                                             '(for cylindrical and spherical)', '0.')

END SUBROUTINE DefineParametersPIC


SUBROUTINE InitPIC()
!===================================================================================================================================
! PIC Init
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation       ,ONLY: InitializeParticleInterpolation
USE MOD_PICDepo                ,ONLY: InitializeDeposition
USE MOD_PIC_Vars ,              ONLY: PICInitIsDone
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

PICInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PIC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPIC

END MODULE MOD_PICInit
