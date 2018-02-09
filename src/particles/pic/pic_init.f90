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
CALL prms%CreateRealArrayOption('PIC-externalField'           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0. , 0. , 0. , 0.')
CALL prms%CreateRealOption(     'PIC-scaleexternalField'      , 'TODO-DEFINE-PARAMETER', '1.0')
CALL prms%CreateLogicalOption(  'PIC-DoInterpolation'         , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateLogicalOption(  'PIC-BG-Field'                , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-curvedexternalField'     , 'TODO-DEFINE-PARAMETER', 'none')
CALL prms%CreateStringOption(   'PIC-variableexternalField'   , 'TODO-DEFINE-PARAMETER', 'none')

CALL prms%CreateIntOption(      'PIC-nCollectChargesBCs'      , 'TODO-DEFINE-PARAMETER', '0')
CALL prms%CreateIntOption(      'PIC-CollectCharges[$]-BC'    , 'TODO-DEFINE-PARAMETER', '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-NumOfRealCharges'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'PIC-CollectCharges[$]-ChargeDist'  , 'TODO-DEFINE-PARAMETER', '0.', numberedmulti=.TRUE.)

CALL prms%SetSection("PIC Deposition")

CALL prms%CreateLogicalOption(  'PIC-DoDeposition'            , 'TODO-DEFINE-PARAMETER', '.TRUE.')
CALL prms%CreateStringOption(   'PIC-Deposition-Type'     , 'TODO-DEFINE-PARAMETER', 'nearest-blurrycenter')

!  r_sf     = GETREAL('PIC-epanechnikov-radius','1.')
!  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
!  alpha_sf = GETINT('PIC-shapefunction-alpha','2')
!  DoSFEqui = GETLOGICAL('PIC-shapefunction-equi','F')
!  NbrOfSFdepoFixes = GETINT('PIC-NbrOfSFdepoFixes','0')
!  PrintSFDepoWarnings=GETLOGICAL('PrintSFDepoWarnings','.FALSE.')
!    SFdepoFixesEps = GETREAL('PIC-SFdepoFixesEps','0.')
!
!        GETREALARRAY('PIC-SFdepoFixes'//TRIM(hilf)//'-Basepoint',3,'0. , 0. , 0.')
!        GETREALARRAY('PIC-SFdepoFixes'//TRIM(hilf)//'-Normal',3,'1. , 0. , 0.') !directed outwards
!        GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ChargeMult','1.')
!      SFdepoFixesBounds(iSFfix,1,1)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-xmin',TRIM(hilf2))
!      SFdepoFixesBounds(iSFfix,1,2)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ymin',TRIM(hilf2))
!      SFdepoFixesBounds(iSFfix,1,3)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-zmin',TRIM(hilf2))
!      SFdepoFixesBounds(iSFfix,2,1)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-xmax',TRIM(hilf2))
!      SFdepoFixesBounds(iSFfix,2,2)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-ymax',TRIM(hilf2))
!      SFdepoFixesBounds(iSFfix,2,3)   = GETREAL('PIC-SFdepoFixes'//TRIM(hilf)//'-zmax',TRIM(hilf2))
!    NbrOfSFdepoFixLinks = GETINT('PIC-NbrOfSFdepoFixLinks','0')
!          GETINTARRAY('PIC-SFdepoFixLink'//TRIM(hilf),2,'1 , 2')
!
!  NbrOfSFdepoLayers = GETINT('PIC-NbrOfSFdepoLayers','0')
!        GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-Basepoint',3,'0. , 0. , 0.')
!        GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-Normal',3,'1. , 0. , 0.') !directed outwards
!      SFdepoLayersBounds(iSFfix,1,1)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-xmin',TRIM(hilf2)),GEO%xmin-r_sf)
!      SFdepoLayersBounds(iSFfix,1,2)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-ymin',TRIM(hilf2)),GEO%ymin-r_sf)
!      SFdepoLayersBounds(iSFfix,1,3)   = MAX(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-zmin',TRIM(hilf2)),GEO%zmin-r_sf)
!      SFdepoLayersBounds(iSFfix,2,1)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-xmax',TRIM(hilf2)),GEO%xmax+r_sf)
!      SFdepoLayersBounds(iSFfix,2,2)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-ymax',TRIM(hilf2)),GEO%ymax+r_sf)
!      SFdepoLayersBounds(iSFfix,2,3)   = MIN(GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-zmax',TRIM(hilf2)),GEO%zmax+r_sf)
!        SFdepoLayersUseFixBounds(iSFfix) = GETLOGICAL('PIC-SFdepoLayers'//TRIM(hilf)//'-UseFixBounds','.TRUE.')
!        GETSTR('PIC-SFdepoLayers'//TRIM(hilf)//'-Space','cuboid')
!        SFdepoLayersBaseVector(iSFfix,1,1:3)   = GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-BaseVector1',3,'0. , 1. , 0.')
!        SFdepoLayersBaseVector(iSFfix,2,1:3)   = GETREALARRAY('PIC-SFdepoLayers'//TRIM(hilf)//'-BaseVector2',3,'0. , 0. , 1.')
!          GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-SFdepoLayersRadius','1.')
!      SFdepoLayersChargedens = GETREAL('PIC-SFdepoLayers'//TRIM(hilf)//'-Chargedens','1.')
!      SFdepoLayersSpec(iSFFix) = GETINT('PIC-SFdepoLayers'//TRIM(hilf)//'-Spec','1')
!
!  SFResampleAnalyzeSurfCollis = GETLOGICAL('PIC-SFResampleAnalyzeSurfCollis','.FALSE.')
!    LastAnalyzeSurfCollis%ReducePartNumber = GETLOGICAL('PIC-SFResampleReducePartNumber','.FALSE.')
!    LastAnalyzeSurfCollis%PartNumThreshold = GETINT('PIC-PartNumThreshold','0')
!      LastAnalyzeSurfCollis%PartNumberReduced = GETINT('PIC-SFResamplePartNumberReduced',TRIM(hilf)) !def. PartNumThreshold
!    LastAnalyzeSurfCollis%Bounds(1,1)   = MAX(GETREAL('PIC-SFResample-xmin',TRIM(hilf)),GEO%xmin-r_sf)
!    LastAnalyzeSurfCollis%Bounds(1,2)   = MAX(GETREAL('PIC-SFResample-ymin',TRIM(hilf)),GEO%ymin-r_sf)
!    LastAnalyzeSurfCollis%Bounds(1,3)   = MAX(GETREAL('PIC-SFResample-zmin',TRIM(hilf)),GEO%zmin-r_sf)
!    LastAnalyzeSurfCollis%Bounds(2,1)   = MIN(GETREAL('PIC-SFResample-xmax',TRIM(hilf)),GEO%xmax+r_sf)
!    LastAnalyzeSurfCollis%Bounds(2,2)   = MIN(GETREAL('PIC-SFResample-ymax',TRIM(hilf)),GEO%ymax+r_sf)
!    LastAnalyzeSurfCollis%Bounds(2,3)   = MIN(GETREAL('PIC-SFResample-zmax',TRIM(hilf)),GEO%zmax+r_sf)
!      LastAnalyzeSurfCollis%UseFixBounds = GETLOGICAL('PIC-SFResample-UseFixBounds','.TRUE.')
!    LastAnalyzeSurfCollis%NormVecOfWall = GETREALARRAY('PIC-NormVecOfWall',3,'1. , 0. , 0.')  !directed outwards
!    LastAnalyzeSurfCollis%Restart = GETLOGICAL('PIC-SFResampleRestart','.FALSE.')
!      LastAnalyzeSurfCollis%DSMCSurfCollisRestartFile = GETSTR('PIC-SFResampleRestartFile','dummy')
!    LastAnalyzeSurfCollis%NumberOfBCs = GETINT('PIC-SFResampleNumberOfBCs','1')
!      LastAnalyzeSurfCollis%BCs = GETINT('PIC-SFResampleSurfCollisBC','0') ! 0 means all...
!      LastAnalyzeSurfCollis%BCs = GETINTARRAY('PIC-SFResampleSurfCollisBC',LastAnalyzeSurfCollis%NumberOfBCs,hilf2)
!
!  r_sf     = GETREAL('PIC-shapefunction-radius','1.')
!  alpha_sf = GETINT ('PIC-shapefunction-alpha','2')
!  sf1d_dir = GETINT ('PIC-shapefunction1d-direction','1')
!
!  r_sf       = GETREAL   ('PIC-shapefunction-radius','1.')
!  r_sf0      = GETREAL   ('PIC-shapefunction-radius0','1.')
!  r_sf_scale = GETREAL   ('PIC-shapefunction-scale','0.')
!  DoSFEqui   = GETLOGICAL('PIC-shapefunction-equi','F')
!
!  DeltaType = GETINT('PIC-DeltaType','1')
!  WRITE(hilf,'(I2)') PP_N
!  NDepo     = GETINT('PIC-DeltaType-N',hilf)
!
!  BGMdeltas(1:3) = GETREALARRAY('PIC-BGMdeltas',3,'0. , 0. , 0.')
!  FactorBGM(1:3) = GETREALARRAY('PIC-FactorBGM',3,'1. , 1. , 1.')

CALL prms%CreateLogicalOption(  'PIC-OutputSource'            , 'TODO-DEFINE-PARAMETER', '.FALSE.')


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
