!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_Analyze_Pure
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalcEkinPart,CalcEkinPart2
!===================================================================================================================================

CONTAINS


PPURE REAL FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2, c2_inv, RelativisticLimit
USE MOD_Particle_Vars ,ONLY: PartState,PartSpecies,Species,usevMPF,PartLorentzType
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting
USE MOD_part_tools    ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: partV2, gamma1, WeightingFactor
!===================================================================================================================================

IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
  WeightingFactor = GetParticleWeight(iPart)
ELSE
  WeightingFactor = GetParticleWeight(iPart) * Species(PartSpecies(iPart))%MacroParticleFactor
END IF

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+DOTPRODUCT(PartState(4:6,iPart))*c2_inv)
  CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
ELSE
  partV2 = DOTPRODUCT(PartState(4:6,iPart))
  IF (partV2.LT.RelativisticLimit)THEN ! |v| < 1000000 when speed of light is 299792458
    CalcEkinPart= 0.5 * Species(PartSpecies(iPart))%MassIC * partV2 * WeightingFactor
  ELSE
    gamma1=partV2*c2_inv
    ! Sanity check: Lorentz factor must be below 1.0
    IF(gamma1.GE.1.0)THEN
      CalcEkinPart=-1.0
    ELSE
      gamma1=1.0/SQRT(1.-gamma1)
      CalcEkinPart=(gamma1-1.0)*Species(PartSpecies(iPart))%MassIC*c2 * WeightingFactor
    END IF ! gamma1.GE.1.0
  END IF ! ipartV2
END IF
END FUNCTION CalcEkinPart


PPURE REAL FUNCTION CalcEkinPart2(velocity,Species_IN,WeightingFactor)
!===================================================================================================================================
! computes the kinetic energy of one particle given its velocity, species and weighting factor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2,c2_inv,RelativisticLimit
USE MOD_Particle_Vars ,ONLY: Species,PartLorentzType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: velocity(1:3)
INTEGER,INTENT(IN)              :: Species_IN
REAL,INTENT(IN)                 :: WeightingFactor
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: partV2, gamma1
!===================================================================================================================================
partV2 = DOT_PRODUCT(velocity,velocity)

IF (PartLorentzType.EQ.5)THEN
  ! gamma v is pushed instead of gamma, therefore, only the relativistic kinetic energy is computed
  ! compute gamma
  gamma1=SQRT(1.0+partV2*c2_inv)
  CalcEkinPart2=(gamma1-1.0)*Species(Species_IN)%MassIC*c2 * WeightingFactor
ELSE
  IF (partV2.LT.RelativisticLimit)THEN ! |v| < 1000000 when speed of light is 299792458
    CalcEkinPart2= 0.5 * Species(Species_IN)%MassIC * partV2 * WeightingFactor
  ELSE
    gamma1=partV2*c2_inv
    ! Sanity check: Lorentz factor must be below 1.0
    IF(gamma1.GE.1.0)THEN
      CalcEkinPart2=-1.0
    ELSE
      gamma1=1.0/SQRT(1.-gamma1)
      CalcEkinPart2=(gamma1-1.0)*Species(Species_IN)%MassIC*c2 * WeightingFactor
    END IF ! gamma1.GE.1.0
  END IF ! ipartV2
END IF
END FUNCTION CalcEkinPart2

#endif /*PARTICLES*/
END MODULE MOD_Particle_Analyze_Pure
