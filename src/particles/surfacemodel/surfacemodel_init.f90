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

CALL prms%CreateIntOption(     'Part-Species[$]-PartBound[$]-ResultSpec'&
                               ,'Resulting recombination species (one of nSpecies)','-1', numberedmulti=.TRUE.)

END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars             ,ONLY: nSpecies,Species,usevMPF
USE MOD_ReadInTools               ,ONLY: GETINT
USE MOD_Particle_Boundary_Vars    ,ONLY: nPartBound, PartBound
USE MOD_SurfaceModel_Vars         ,ONLY: SurfModResultSpec, SurfModEnergyDistribution
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32) :: hilf,hilf2
REAL          :: MPFiSpec,MPFresultSpec
INTEGER       :: iSpec,iPartBound
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%Reactive))) RETURN

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
    CASE(SEE_MODELS_ID)
      ! 5: SEE by Levko2015
      ! 6: SEE by Pagonakis2016 (originally from Harrower1956)
      ! 7: SEE-I (bombarding electrons are removed, Ar+ on different materials is considered for SEE)
      ! 8: SEE-E (e- on dielectric materials is considered for SEE and three different outcomes)
      ! 9: SEE-I when Ar^+ ion bombards surface with 0.01 probability and fixed SEE electron energy of 6.8 eV
      !10: SEE-I (bombarding electrons are removed, Ar+ on copper is considered for SEE)
!-----------------------------------------------------------------------------------------------------------------------------------
      SurfModResultSpec(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-ResultSpec')
      ! Check that the impacting and SEE particles have the same MPF is vMPF is turned off
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
      ! Set specific distributions functions
      IF(PartBound%SurfaceModel(iPartBound).EQ.8)THEN
        SurfModEnergyDistribution = 'Morozov2004'
      ELSE
        SurfModEnergyDistribution = 'deltadistribution'
      END IF ! PartBound%SurfaceModel(iPartBound).EQ.8
!-----------------------------------------------------------------------------------------------------------------------------------
    END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
  END DO ! iPartBound=1,nPartBound
END DO ! iSpec = 1,nSpecies

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
END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
