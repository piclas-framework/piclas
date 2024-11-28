!==================================================================================================================================
! Copyright (c) 2024 boltzplatz - numerical plasma dynamics GmbH, Simone Lauterbach, Marcel Pfeiffer
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

MODULE MOD_ParticleWeighting
!===================================================================================================================================
!>
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
PUBLIC :: DefineParametersParticleWeighting, InitParticleWeighting
!===================================================================================================================================

INTEGER,PARAMETER      :: PRM_PARTWEIGHT_CONSTANT   = 0
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_RADIAL     = 1
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_LINEAR     = 2
INTEGER,PARAMETER      :: PRM_PARTWEIGHT_CELLLOCAL  = 3

CONTAINS

!==================================================================================================================================
!> Define parameters for particles weighting
!==================================================================================================================================
SUBROUTINE DefineParametersParticleWeighting()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE

CALL prms%SetSection("Particle Weighting")
CALL prms%CreateIntFromStringOption('Part-Weight-Type', "Particle weighting type: \n"             //&
                                      'constant ('//TRIM(int2strf(PRM_PARTWEIGHT_CONSTANT))//')\n'  //&
                                      'radial ('//TRIM(int2strf(PRM_PARTWEIGHT_RADIAL))//')\n'      //&
                                      'linear ('//TRIM(int2strf(PRM_PARTWEIGHT_LINEAR))//')\n'      //&
                                      'cell_local ('//TRIM(int2strf(PRM_PARTWEIGHT_CELLLOCAL))//')' &
                                      ,'constant')
CALL addStrListEntry('Part-Weight-Type' , 'constant'   , PRM_PARTWEIGHT_CONSTANT)
CALL addStrListEntry('Part-Weight-Type' , 'radial'     , PRM_PARTWEIGHT_RADIAL)
CALL addStrListEntry('Part-Weight-Type' , 'linear'     , PRM_PARTWEIGHT_LINEAR)
CALL addStrListEntry('Part-Weight-Type' , 'cell_local' , PRM_PARTWEIGHT_CELLLOCAL)

END SUBROUTINE DefineParametersParticleWeighting


SUBROUTINE InitParticleWeighting()
!===================================================================================================================================
!> Initialize the particle weighting, especially the particle cloning
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools           ,ONLY: GETLOGICAL,GETINTFROMSTR
USE MOD_Symmetry_Vars         ,ONLY: Symmetry
USE MOD_DSMC_Vars             ,ONLY: DoRadialWeighting, DoLinearWeighting, DoCellLocalWeighting, ParticleWeighting
USE MOD_DSMC_Symmetry         ,ONLY: InitParticleCloning
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: ParticleWeightType
!===================================================================================================================================

DoRadialWeighting = .FALSE.
DoLinearWeighting = .FALSE.
DoCellLocalWeighting = .FALSE.
ParticleWeighting%PerformCloning = .FALSE.

ParticleWeightType = GETINTFROMSTR('Part-Weight-Type')

SELECT CASE(ParticleWeightType)
Case(PRM_PARTWEIGHT_CONSTANT)
  ParticleWeighting%Type = 'constant'
Case(PRM_PARTWEIGHT_RADIAL)
  ParticleWeighting%Type = 'radial'
  DoRadialWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
  IF(.NOT.Symmetry%Axisymmetric) CALL CollectiveStop(__STAMP__,' Part-Weight-Type = radial requires an axisymmetric simulation!')
Case(PRM_PARTWEIGHT_LINEAR)
  ParticleWeighting%Type = 'linear'
  DoLinearWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
Case(PRM_PARTWEIGHT_CELLLOCAL)
  ParticleWeighting%Type = 'cell_local'
  DoCellLocalWeighting = .TRUE.
  ParticleWeighting%PerformCloning = .TRUE.
  ParticleWeighting%EnableOutput = .TRUE.
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown Part-Weight-Type!' ,IntInfo=ParticleWeightType)
END SELECT

IF(ParticleWeighting%PerformCloning) CALL InitParticleCloning()

END SUBROUTINE InitParticleWeighting

END MODULE MOD_ParticleWeighting