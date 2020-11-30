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

MODULE MOD_Particle_Analyze_Tools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CalcEkinPart
  MODULE PROCEDURE CalcEkinPart
END INTERFACE

INTERFACE CalcEkinPart2
  MODULE PROCEDURE CalcEkinPart2
END INTERFACE

INTERFACE CalcNumPartsOfSpec
  MODULE PROCEDURE CalcNumPartsOfSpec
END INTERFACE

PUBLIC :: CalcEkinPart,CalcEkinPart2
PUBLIC :: CalcNumPartsOfSpec
!===================================================================================================================================

CONTAINS


PURE REAL FUNCTION CalcEkinPart(iPart)
!===================================================================================================================================
! computes the kinetic energy of one particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2, c2_inv
USE MOD_Particle_Vars ,ONLY: PartState, PartSpecies, Species
USE MOD_PARTICLE_Vars ,ONLY: usevMPF
USE MOD_Particle_Vars ,ONLY: PartLorentzType
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
  IF (partV2.LT.1E12)THEN
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


PURE REAL FUNCTION CalcEkinPart2(velocity,Species_IN,WeightingFactor)
!===================================================================================================================================
! computes the kinetic energy of one particle given its velocity, species and weighting factor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Globals_Vars  ,ONLY: c2,c2_inv
USE MOD_Particle_Vars ,ONLY: Species
USE MOD_Particle_Vars ,ONLY: PartLorentzType
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
  IF (partV2.LT.1E12)THEN
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


SUBROUTINE CalcNumPartsOfSpec(NumSpec,SimNumSpec,CalcNumSpec_IN,CalcSimNumSpec_IN)
!===================================================================================================================================
! Computes the number of simulated particles AND number of real particles within the domain
! Last section of the routine contains the MPI-communication
! CAUTION: SimNumSpec equals NumSpec only for constant weighting factor
! NOTE: Background gas particles are not considered
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Particle_Vars         ,ONLY: PDM,PartSpecies
USE MOD_Particle_Analyze_Vars ,ONLY: nSpecAnalyze
USE MOD_DSMC_Vars             ,ONLY: BGGas
USE MOD_Particle_Vars         ,ONLY: nSpecies
USE MOD_part_tools            ,ONLY: GetParticleWeight
#if USE_MPI
USE MOD_Particle_MPI_Vars     ,ONLY: PartMPI
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: NumSpec(nSpecAnalyze)
INTEGER(KIND=IK),INTENT(OUT)       :: SimNumSpec(nSpecAnalyze)
LOGICAL,INTENT(IN)                 :: CalcNumSpec_IN,CalcSimNumSpec_IN ! Flags for performing MPI reduce
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iPart, iSpec
!===================================================================================================================================

NumSpec    = 0.
SimNumSpec = 0
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    ! NumSpec = real particle number
    NumSpec(PartSpecies(iPart))    = NumSpec(PartSpecies(iPart)) + GetParticleWeight(iPart)
    ! SimNumSpec = simulated particle number
    SimNumSpec(PartSpecies(iPart)) = SimNumSpec(PartSpecies(iPart)) + 1
  END IF
END DO
IF(BGGas%NumberOfSpecies.GT.0) THEN
  DO iSpec = 1, nSpecies
    IF(BGGas%BackgroundSpecies(iSpec)) THEN
      NumSpec(iSpec) = 0.
      SimNumSpec(iSpec) = 0
    END IF
  END DO
END IF
IF(nSpecAnalyze.GT.1)THEN
  NumSpec(nSpecAnalyze)    = SUM(NumSpec(1:nSpecies))
  SimNumSpec(nSpecAnalyze) = SUM(SimNumSpec(1:nSpecies))
END IF

#if USE_MPI
IF (PartMPI%MPIRoot) THEN
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(MPI_IN_PLACE,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(MPI_IN_PLACE,SimNumSpec ,nSpecAnalyze,MPI_INTEGER_INT_KIND,MPI_SUM,0,PartMPI%COMM,IERROR)
ELSE
  IF(CalcNumSpec_IN)    CALL MPI_REDUCE(NumSpec     ,NumSpec    ,nSpecAnalyze,MPI_DOUBLE_PRECISION,MPI_SUM,0,PartMPI%COMM,IERROR)
  IF(CalcSimNumSpec_IN) CALL MPI_REDUCE(SimNumSpec  ,SimNumSpec ,nSpecAnalyze,MPI_INTEGER_INT_KIND,MPI_SUM,0,PartMPI%COMM,IERROR)
END IF
#endif /*USE_MPI*/

! Set global number of particles (info for std out)
IF(CalcSimNumSpec_IN)THEN
  GlobalNbrOfParticlesUpdated = .TRUE.
#if USE_MPI
  IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
    nGlobalNbrOfParticles = INT(SimNumSpec(nSpecAnalyze),KIND=IK)
#if USE_MPI
  END IF ! PartMPI%MPIRoot
#endif /*USE_MPI*/
END IF ! CalcSimNumSpec_IN

END SUBROUTINE CalcNumPartsOfSpec


#endif /*PARTICLES*/


END MODULE MOD_Particle_Analyze_Tools
