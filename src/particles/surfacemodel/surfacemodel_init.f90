!==================================================================================================================================
! Copyright (c) 2015-2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel_Init
!===================================================================================================================================
!> Module for initialization of surface models
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
PUBLIC :: DefineParametersSurfModel
PUBLIC :: InitSurfaceModel
PUBLIC :: FinalizeSurfaceModel
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for surface model
!==================================================================================================================================
SUBROUTINE DefineParametersSurfModel()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("SurfaceModel")

CALL prms%CreateIntOption( 'Part-Species[$]-PartBound[$]-ResultSpec','Resulting recombination species (one of nSpecies)',&
                           '-1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption('Part-SurfaceModel-SEE-Te','Bulk electron temperature for SEE model by Morozov2004 in Kelvin (default'//&
                           ' corresponds to 50 eV)','5.80226250308285e5')
CALL prms%CreateLogicalOption( 'Part-SurfaceModel-SEE-Te-automatic','Automatically set the bulk electron temperature by using '//&
                               'the global electron temperature for SEE model by Morozov2004', '.FALSE.')
CALL prms%CreateIntOption( 'Part-SurfaceModel-SEE-Te-Spec','Electron species (one of nSpecies) when using the global electron '//&
                           'temperature for SEE model by Morozov2004')

END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: Kelvin2eV
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_ReadInTools            ,ONLY: GETINT,GETREAL,GETLOGICAL
USE MOD_Particle_Boundary_Vars ,ONLY: nPartBound,PartBound
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModResultSpec,SurfModEnergyDistribution,SurfModSEEelectronTemp
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModSEEelectronTempAutoamtic,SurfModSEEelectronTempSpecies
USE MOD_HDF5_Input             ,ONLY: DatasetExists,ReadArray
USE MOD_IO_HDF5                ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_Restart_Vars           ,ONLY: RestartFile,DoRestart
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32) :: hilf, hilf2
INTEGER       :: iSpec, iPartBound
LOGICAL       :: ReadSurfModSEEelectronTemp,SurfModSEEelectronTempExists
REAL          :: TmpArray(1,1)
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%Reactive))) RETURN

ReadSurfModSEEelectronTemp = .FALSE. ! Initialize

ALLOCATE(SurfModResultSpec(1:nPartBound,1:nSpecies))
SurfModResultSpec = 0
ALLOCATE(SurfModEnergyDistribution(1:nPartBound))
SurfModEnergyDistribution = ''
! initialize model specific variables
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  DO iPartBound=1,nPartBound
    IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
    WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
    hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
    SELECT CASE(PartBound%SurfaceModel(iPartBound))
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(5,6,7,8)
      ! 5: SEE by Levko2015
      ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
      ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
      ! 8: SEE-E (bombarding electrons are reflected, e- on dielectric materials is considered for SEE and three different outcomes)
!-----------------------------------------------------------------------------------------------------------------------------------
      SurfModResultSpec(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-ResultSpec')
      IF(PartBound%SurfaceModel(iPartBound).EQ.8)THEN
        SurfModEnergyDistribution  = 'Morozov2004'
        ReadSurfModSEEelectronTemp = .TRUE.
      ELSE
        SurfModEnergyDistribution  = 'deltadistribution'
      END IF ! PartBound%SurfaceModel(iPartBound).EQ.8
!-----------------------------------------------------------------------------------------------------------------------------------
    END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
  END DO
END DO


! If SEE model by Morozov is used, read the additional parameter for the electron bulk temperature
IF(ReadSurfModSEEelectronTemp)THEN
  SurfModSEEelectronTemp = GETREAL('Part-SurfaceModel-SEE-Te') ! default is 50 eV = 5.80226250308285e5 K
  SurfModSEEelectronTemp = SurfModSEEelectronTemp*Kelvin2eV    ! convert to eV to be used in the code
  SurfModSEEelectronTempAutoamtic = GETLOGICAL('Part-SurfaceModel-SEE-Te-automatic')
  IF(SurfModSEEelectronTempAutoamtic)THEN
    SurfModSEEelectronTempSpecies = GETINT('Part-SurfaceModel-SEE-Te-Spec')
  END IF ! SurfModSEEelectronTempAutoamtic

  ! Restart: Only root reads state file to prevent access with a large number of processors
  IF(MPIRoot)THEN
    IF(DoRestart)THEN
      CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
      CALL DatasetExists(File_ID,'SurfModSEEelectronTemp',SurfModSEEelectronTempExists)
      IF(SurfModSEEelectronTempExists)THEN
        CALL ReadArray('SurfModSEEelectronTemp',2,(/1_IK,1_IK/),0_IK,2,RealArray=TmpArray(1,1))
        SurfModSEEelectronTemp = TmpArray(1,1)
        WRITE(UNIT_stdOut,'(1(A,ES10.2E3))') " Read SurfModSEEelectronTemp from restart file ["//TRIM(RestartFile)//"] Te[eV]: ",SurfModSEEelectronTemp
      END IF ! RegionElectronRefExists
      CALL CloseDataFile()
    END IF ! DoRestart
  END IF ! MPIRoot
#if USE_MPI
  ! Broadcast from root to other processors
  CALL MPI_BCAST(SurfModSEEelectronTemp,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iERROR)
#endif /*USE_MPI*/
ELSE
  SurfModSEEelectronTempAutoamtic = .FALSE. ! Default
END IF ! ReadSurfModSEEelectronTemp

END SUBROUTINE InitSurfaceModel


SUBROUTINE FinalizeSurfaceModel()
!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars
USE MOD_SurfaceModel_Analyze_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SurfModelAnalyzeInitIsDone=.FALSE.

SDEALLOCATE(SurfModResultSpec)
SDEALLOCATE(SurfModEnergyDistribution)

! === Surface Analyze Vars
SDEALLOCATE(SurfAnalyzeCount)
SDEALLOCATE(SurfAnalyzeNumOfAds)
SDEALLOCATE(SurfAnalyzeNumOfDes)
IF(CalcBoundaryParticleOutput)THEN
  SDEALLOCATE(BPO%RealPartOut)
  SDEALLOCATE(BPO%PartBoundaries)
  SDEALLOCATE(BPO%BCIDToBPOBCID)
  SDEALLOCATE(BPO%Species)
  SDEALLOCATE(BPO%SpecIDToBPOSpecID)
END IF ! CalcBoundaryParticleOutput

END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
