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

CALL prms%CreateIntOption(      'MacroBody-nMacroBody'  &
                           , 'Number of macro particle, which are checked during tracing',  '0')
CALL prms%CreateLogicalOption(  'MacroBody-AccelerationEnabled'  &
                           , 'Enables momentum changes of macro particle',  '.FALSE.')
CALL prms%CreateLogicalOption(  'MacroBody-FluxesEnabled'  &
                           , 'Enables mass and energy changes of macro particle',  '.FALSE.')
CALL prms%CreateLogicalOption(  'MacroBody-WriteElemData'  &
                           , 'Enables write out of elem data for Macro-spheres in state file. e.g. volumeportion','.FALSE.')
CALL prms%CreateRealArrayOption('MacroBody[$]-center'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('MacroBody[$]-velocity'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('MacroBody[$]-rotation'  &
                           , 'TODO-DEFINE-PARAMETER', '0. , 0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-radius'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-temp'  &
                           , 'TODO-DEFINE-PARAMETER',  '273.15', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-density'  &
                           , 'TODO-DEFINE-PARAMETER',  '997', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-momentumACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-transACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-vibACC'  &
                           , 'TODO-DEFINE-PARAMETER',  '1.0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'MacroBody[$]-rotACC'  &
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
INTEGER       :: iMB
!----------------------------------------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==43)
nMacroBody = GETINT('MacroBody-nMacroBody')
#else
nMacroBody = 0
#endif
IF (nMacroBody.GT.0) THEN
  IF (DoRefMapping.OR.TriaTracking) CALL abort(&
      __STAMP__&
      ,'Macroparticle not possible with dorefmapping or TriaTracking')
  ! if implementation for triatracking intended, fix number of envelopes in halo region build (particle_mpi_halo.f90)
  UseMacroBody=.TRUE.
  MacroBodyFluxesEnabled = GETLOGICAL('MacroBody-FluxesEnabled')
  MacroBodyAccelerationEnabled = GETLOGICAL('MacroBody-AccelerationEnabled')
  ALLOCATE (MacroSphere(1:nMacroBody))
  DO iMB = 1,nMacroBody
    WRITE(UNIT=hilf,FMT='(I0)') iMB
    MacroSphere(iMB)%center=GETREALARRAY('MacroBody'//TRIM(hilf)//'-center',3)
    MacroSphere(iMB)%velocity(1:3)=GETREALARRAY('MacroBody'//TRIM(hilf)//'-velocity',3)
    MacroSphere(iMB)%velocity(4:6)=GETREALARRAY('MacroBody'//TRIM(hilf)//'-rotation',3)
    MacroSphere(iMB)%radius=GETREAL('MacroBody'//TRIM(hilf)//'-radius')
    MacroSphere(iMB)%temp=GETREAL('MacroBody'//TRIM(hilf)//'-temp')
    MacroSphere(iMB)%density=GETREAL('MacroBody'//TRIM(hilf)//'-density')
    IF (MacroSphere(iMB)%density.LE.0.) CALL abort(&
        __STAMP__&
        ,'density must be above 0 for MacroBody',iMB)
    MacroSphere(iMB)%mass=4./3.*MacroSphere(iMB)%radius**3*PI*MacroSphere(iMB)%density
    MacroSphere(iMB)%momentumACC=GETREAL('MacroBody'//TRIM(hilf)//'-momentumACC')
    MacroSphere(iMB)%transAcc=GETREAL('MacroBody'//TRIM(hilf)//'-transACC')
    MacroSphere(iMB)%vibAcc=GETREAL('MacroBody'//TRIM(hilf)//'-vibACC')
    MacroSphere(iMB)%rotACC=GETREAL('MacroBody'//TRIM(hilf)//'-rotACC')
    MacroSphere(iMB)%RHS(:)=0.
  END DO
  CalcMPVolumePortion=.TRUE.
ELSE
  UseMacroBody=.FALSE.
  MacroBodyFluxesEnabled=.FALSE.
  MacroBodyAccelerationEnabled=.FALSE.
  CalcMPVolumePortion=.FALSE.
END IF
ConsiderVolumePortions=.FALSE.
IF (UseMacroBody) THEN
  ConsiderVolumePortions=.TRUE.
  ALLOCATE(ElemHasMacroBody(1:nTotalElems, 1:nMacroBody))
  ElemHasMacroBody(:,:)=.FALSE.
  MacroPartWriteElemData=GETLOGICAL('MacroBody-WriteElemData')
  IF (MacroPartWriteElemData) THEN
    CALL AddToElemData(ElementOut,'ElemHasMacroBody',LogArray=ElemHasMacroBody(:,1))
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
