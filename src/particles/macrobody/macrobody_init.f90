!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

!===================================================================================================================================
!> Module for initializing macroscopic bodies inside particle domain
!===================================================================================================================================
MODULE MOD_MacroBody_Init
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC :: DefineParametersMacroBody
PUBLIC :: InitMacroBody
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters for macrocopic bodies
!==================================================================================================================================
SUBROUTINE DefineParametersMacroBody()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("MacroParticle")

CALL prms%CreateIntOption(      'MacroPart-nMacroParticle'  &
                           , 'Number of macro particle, which are checked during tracing',  '0')
CALL prms%CreateLogicalOption(  'MacroPart-AccelerationEnabled'  &
                           , 'Enables momentum changes of macro particle',  '.FALSE.')
CALL prms%CreateLogicalOption(  'MacroPart-FluxesEnabled'  &
                           , 'Enables mass and energy changes of macro particle',  '.FALSE.')
CALL prms%CreateLogicalOption(  'MacroPart-WriteElemData'  &
                           , 'Enables write out of elem data for Macro-spheres in state file. e.g. volumeportion','.FALSE.')
CALL prms%CreateRealArrayOption('MacroPart[$]-center'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('MacroPart[$]-velocity'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('MacroPart[$]-rotation'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-radius'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-temp'  &
                           , 'TODO-DEFINE-PARAMETER',  '273.15', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-density'  &
                           , 'TODO-DEFINE-PARAMETER',  '997', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-momentumACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-transACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-vibACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroPart[$]-rotACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersMacroBody


!===================================================================================================================================
!> initialize variables used for macrobodies in domain
!===================================================================================================================================
SUBROUTINE InitMacroBody()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI
USE MOD_ReadInTools            ,ONLY: GETINT, GETLOGICAL, GETREAL, GETREALARRAY
USE MOD_MacroBody_Vars
USE MOD_Particle_Vars          ,ONLY: nPointsMCVolumeEstimate
USE MOD_DSMC_Vars              ,ONLY: ConsiderVolumePortions
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping, TriaTracking
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO, nTotalElems
USE MOD_IO_HDF5                ,ONLY: AddToElemData, ElementOut
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
CHARACTER(32) :: hilf
INTEGER       :: iMP
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==43)
nMacroParticle = GETINT('MacroPart-nMacroParticle')
#else
nMacroParticle = 0
#endif
IF (nMacroparticle.GT.0) THEN
  IF (DoRefMapping.OR.TriaTracking) CALL abort(&
      __STAMP__&
      ,'Macroparticle not possible with dorefmapping or TriaTracking')
  ! if implementation for triatracking intended, fix number of envelopes in halo region build (particle_mpi_halo.f90)
  UseMacroPart=.TRUE.
  MacroPartFluxesEnabled = GETLOGICAL('MacroPart-FluxesEnabled')
  MacroPartFluxesEnabled = GETLOGICAL('MacroPart-AccelerationEnabled')
  ALLOCATE (MacroPart(1:nMacroParticle))
  DO iMP = 1,nMacroParticle
    WRITE(UNIT=hilf,FMT='(I0)') iMP
    MacroPart(iMP)%center=GETREALARRAY('MacroPart'//TRIM(hilf)//'-center',3)
    MacroPart(iMP)%velocity(1:3)=GETREALARRAY('MacroPart'//TRIM(hilf)//'-velocity',3)
    MacroPart(iMP)%velocity(4:6)=GETREALARRAY('MacroPart'//TRIM(hilf)//'-rotation',3)
    MacroPart(iMP)%radius=GETREAL('MacroPart'//TRIM(hilf)//'-radius')
    MacroPart(iMP)%temp=GETREAL('MacroPart'//TRIM(hilf)//'-temp')
    MacroPart(iMP)%density=GETREAL('MacroPart'//TRIM(hilf)//'-density')
    IF (MacroPart(iMP)%density.LE.0.) CALL abort(&
        __STAMP__&
        ,'density must be above 0 for MacroPart',iMP)
    MacroPart(iMP)%mass=4./3.*MacroPart(iMP)%radius**3*PI*MacroPart(iMP)%density
    MacroPart(iMP)%momentumACC=GETREAL('MacroPart'//TRIM(hilf)//'-momentumACC')
    MacroPart(iMP)%transAcc=GETREAL('MacroPart'//TRIM(hilf)//'-transACC')
    MacroPart(iMP)%vibAcc=GETREAL('MacroPart'//TRIM(hilf)//'-vibACC')
    MacroPart(iMP)%rotACC=GETREAL('MacroPart'//TRIM(hilf)//'-rotACC')
    MacroPart(iMP)%RHS(:)=0.
  END DO
  CalcMPVolumePortion=.TRUE.
ELSE
  UseMacroPart=.FALSE.
  MacroPartFluxesEnabled=.FALSE.
  CalcMPVolumePortion=.FALSE.
END IF
ConsiderVolumePortions=.FALSE.
IF (UseMacropart) THEN
  ConsiderVolumePortions=.TRUE.
  ALLOCATE(ElemHasMacroPart(1:nTotalElems, 1:nMacroParticle))
  ElemHasMacroPart(:,:)=.FALSE.
  MacroPartWriteElemData=GETLOGICAL('MacroPart-WriteElemData')
  IF (MacroPartWriteElemData) THEN
    CALL AddToElemData(ElementOut,'ElemHasMacroPart',LogArray=ElemHasMacroPart(:,1))
    CALL AddToElemData(ElementOut,'MPVolumePortion',RealArray=GEO%MPVolumePortion(:))
  END IF
END IF
IF (ConsiderVolumePortions) THEN
  nPointsMCVolumeEstimate = GETINT('Particles-nPointsMCVolumeEstimate')
  IF (nPointsMCVolumeEstimate.LT.1) CALL abort(&
      __STAMP__&
      ,'nPointsMCVolumeEstimate is must be above 0')
END IF

END SUBROUTINE InitMacroBody


END MODULE MOD_MacroBody_Init
