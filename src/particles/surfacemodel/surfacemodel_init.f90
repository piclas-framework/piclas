!==================================================================================================================================
! Copyright (c) 2015-2019 Wladimir Reschke
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

CALL prms%CreateIntOption(     'Part-Species[$]-PartBound[$]-ResultSpec'    , 'Resulting recombination species (one of nSpecies)' , '-1' , numberedmulti=.TRUE.)
CALL prms%CreateStringOption(  'Part-Boundary[$]-SurfModEnergyDistribution' , 'Energy distribution function for surface emission model (only changable for SurfaceModel=7)' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(    'Part-Boundary[$]-SurfModEmissionEnergy'     , 'Energy of emitted particle for surface emission model (only available for SurfaceModel=7)' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(    'Part-Boundary[$]-SurfModEmissionYield'      , 'Emission yield factor for surface emission model (only changable for SurfaceModel=7)' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(    'Part-SurfaceModel-SEE-Te'                   , 'Bulk electron temperature for SEE model by Morozov2004 in Kelvin (default corresponds to 50 eV)' , '5.80226250308285e5')
CALL prms%CreateLogicalOption( 'Part-SurfaceModel-SEE-Te-automatic'         , 'Automatically set the bulk electron temperature by using the global electron temperature for SEE model by Morozov2004' , '.FALSE.')

CALL prms%CreateRealArrayOption('Part-Boundary[$]-SurfModSEEPowerFit'       , 'SEE Power-fit model (SurfaceModel = 4): coefficients of the form a*E(eV)^b, input as (a,b)', numberedmulti=.TRUE.,no=2)
END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: Kelvin2eV
USE MOD_Particle_Vars          ,ONLY: nSpecies,Species,usevMPF
USE MOD_ReadInTools            ,ONLY: GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY
USE MOD_Particle_Boundary_Vars ,ONLY: nPartBound,PartBound
USE MOD_SurfaceModel_Vars      ,ONLY: BulkElectronTempSEE,SurfModSEEelectronTempAutoamtic
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModResultSpec,SurfModEnergyDistribution,SurfModEmissionEnergy,SurfModEmissionYield
USE MOD_SurfaceModel_Vars      ,ONLY: SurfModSEEPowerFit
USE MOD_Particle_Vars          ,ONLY: CalcBulkElectronTemp,BulkElectronTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)        :: hilf, hilf2, hilf3
INTEGER              :: iSpec, iPartBound
LOGICAL              :: SurfModelElectronTemp
INTEGER, ALLOCATABLE :: SumOfResultSpec(:)
REAL                 :: MPFiSpec,MPFresultSpec
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%Reactive))) RETURN

! Initialize
SurfModelElectronTemp = .FALSE.
CalcBulkElectronTemp  = .FALSE.

ALLOCATE(SurfModResultSpec(1:nPartBound,1:nSpecies))
SurfModResultSpec = 0
ALLOCATE(SurfModEnergyDistribution(1:nPartBound))
SurfModEnergyDistribution = ''
ALLOCATE(SurfModEmissionEnergy(1:nPartBound))
SurfModEmissionEnergy = -2.0
ALLOCATE(SurfModEmissionYield(1:nPartBound))
SurfModEmissionYield = 0.
ALLOCATE(SumOfResultSpec(nPartBound))
SumOfResultSpec = 0

ALLOCATE(SurfModSEEPowerFit(1:2, 1:nPartBound))
SurfModSEEPowerFit = 0

! Loop particle boundaries
DO iPartBound=1,nPartBound
  IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
  WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Read-in parameters for all SEE models
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  CASE(SEE_MODELS_ID)
    ! Loop all species
    DO iSpec = 1,nSpecies
      WRITE(UNIT=hilf,FMT='(I0)') iSpec
      hilf3=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
      SurfModResultSpec(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf3)//'-ResultSpec')
      SumOfResultSpec(iPartBound)         = SumOfResultSpec(iPartBound) + SurfModResultSpec(iPartBound,iSpec)
      IF(SumOfResultSpec(iPartBound).EQ.-nSpecies) CALL abort(__STAMP__,&
          'SEE surface model: All resulting species are -1. Define at least one species that can be created by an SEE event.')
      ! Check that the impacting and SEE particles have the same MPF if vMPF is turned off
      IF(.NOT.usevMPF)THEN
        ! Skip non-initialized values
        IF(SurfModResultSpec(iPartBound,iSpec).NE.-1)THEN
          MPFiSpec      = Species(iSpec)%MacroParticleFactor
          MPFresultSpec = Species(SurfModResultSpec(iPartBound,iSpec))%MacroParticleFactor
          IF(.NOT.(ALMOSTEQUALRELATIVE(MPFiSpec,MPFresultSpec,1e-3)))THEN
            IPWRITE(UNIT_StdOut,*) "Bombarding particle: SpecID =", iSpec
            IPWRITE(UNIT_StdOut,*) "Bombarding particle:    MPF =", MPFiSpec
            IPWRITE(UNIT_StdOut,*) "Secondary electron : SpecID =", SurfModResultSpec(iPartBound,iSpec)
            IPWRITE(UNIT_StdOut,*) "Secondary electron :    MPF =", MPFresultSpec
            CALL abort(__STAMP__,'SEE model: MPF of bomarding particle and secondary electron must be the same.')
          END IF ! .NOT.(ALMOSTEQUALRELATIVE(MPFiSpec,MPFresultSpec,1e-3))
        END IF ! MPFresultSpec.NE.-1
      END IF ! .NOT.usevMPF
    END DO ! iSpec = 1,nSpecies
  END SELECT
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Set the default energy distribution
  SurfModEnergyDistribution(iPartBound)  = 'deltadistribution'
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Read-in and set SEE model specific parameters
  SELECT CASE(PartBound%SurfaceModel(iPartBound))
  ! 4: SEE Power-fit model by Goebel & Katz „Fundamentals of Electric Propulsion - Ion and Hall Thrusters“
  CASE(4)
    CALL abort(__STAMP__,'SEE model power fit: Implementation not yet finished, regression test and documentation is missing!')
    SurfModSEEPowerFit(1:2,iPartBound) = GETREALARRAY('Part-Boundary'//TRIM(hilf2)//'-SurfModSEEPowerFit',2)
  ! 5: SEE by Levko2015
  ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
  ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
  CASE(7)
    ! Note that the define vars help needs to be changed as soon as these parameters are available for other surface models
    SurfModEnergyDistribution(iPartBound) = TRIM(GETSTR('Part-Boundary'//TRIM(hilf2)//'-SurfModEnergyDistribution','deltadistribution'))
    SurfModEmissionEnergy(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilf2)//'-SurfModEmissionEnergy','-1.0')
    IF((SurfModEmissionEnergy(iPartBound).LE.0.).AND.(SurfModEnergyDistribution(iPartBound).EQ.'uniform-energy')) CALL abort(&
        __STAMP__,'SEE model with uniform-energy distribution requires Part-BoundaryX-SurfModEmissionEnergy > 0.')
    SurfModEmissionYield(iPartBound)      = GETREAL('Part-Boundary'//TRIM(hilf2)//'-SurfModEmissionYield' ,'0.13')
  ! 8: SEE-E (e- on dielectric materials is considered for SEE and three different outcomes)
  CASE(8)
    SurfModEnergyDistribution(iPartBound)  = 'Morozov2004'
    SurfModelElectronTemp = .TRUE.
  ! 9: SEE-I when Ar^+ ion bombards surface with 0.01 probability and fixed SEE electron energy of 6.8 eV
  !10: SEE-I (bombarding electrons are removed, Ar+ on copper is considered for SEE)
  END SELECT
  !---------------------------------------------------------------------------------------------------------------------------------
END DO ! iPartBound=1,nPartBound

DEALLOCATE(SumOfResultSpec)

! If SEE model by Morozov is used, read the additional parameter for the electron bulk temperature
IF(SurfModelElectronTemp)THEN
  BulkElectronTempSEE             = GETREAL('Part-SurfaceModel-SEE-Te') ! default is 50 eV = 5.80226250308285e5 K
  BulkElectronTempSEE             = BulkElectronTempSEE*Kelvin2eV       ! convert to eV to be used in the code
  SurfModSEEelectronTempAutoamtic = GETLOGICAL('Part-SurfaceModel-SEE-Te-automatic')
  IF(SurfModSEEelectronTempAutoamtic)THEN
    CalcBulkElectronTemp = .TRUE.
    BulkElectronTemp     = BulkElectronTempSEE
  END IF
END IF ! SurfModelElectronTemp

END SUBROUTINE InitSurfaceModel


SUBROUTINE FinalizeSurfaceModel()
!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(SurfModResultSpec)
SDEALLOCATE(SurfModEnergyDistribution)
SDEALLOCATE(SurfModEmissionEnergy)
SDEALLOCATE(SurfModEmissionYield)
SDEALLOCATE(StickingCoefficientData)
SDEALLOCATE(SurfModSEEPowerFit)
END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
