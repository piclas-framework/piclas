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

MODULE MOD_PICDepo_HDG
!===================================================================================================================================
! MOD PIC Depo HDG
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!===================================================================================================================================
#if USE_HDG
PUBLIC :: DepositVirtualDielectricLayerParticles
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Loop over all particles and find that ones that have been flagged during particle-boundary interaction and have hit a VDL
!> boundary. They are flagged there as they might move to another process during that interaction.
!> Here, after MPI communication, they can be deleted and deposited at the target position to form a surface charge on a (virtual)
!> dielectric layer.
!===================================================================================================================================
SUBROUTINE DepositVirtualDielectricLayerParticles()
! MODULES
USE MOD_Particle_Vars   ,ONLY: PEM, PDM, Species, PartSpecies, usevmpf, PartMPF, PartState
USE MOD_Particle_Vars   ,ONLY: IsVDLSpecID, SpeciesOffsetVDL
USE MOD_PICDepo_Tools   ,ONLY: DepositParticleOnNodes
USE MOD_part_operations ,ONLY: RemoveParticle
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iPart
REAL    :: charge,SignSwitch
!===================================================================================================================================
! Loop over all particles
DO iPart=1,PDM%ParticleVecLength
  ! Only consider un-deleted particles
  IF (PDM%ParticleInside(iPart)) THEN
    ! Check particle index for VDL particles
    IF(IsVDLSpecID(iPart))THEN
      ! Check for negative sign
      IF(PartSpecies(iPart).LT.0)THEN
        ! If negative sign is found in the species index, invert the deposited charge
        SignSwitch = -1
        ! Invert the species index so it is meaningful again
        PartSpecies(iPart) = -PartSpecies(iPart)
      ELSE
        ! Use same sign for charge deposition
        SignSwitch =  1
      END IF ! PartSpecies(iPart).LT.0
      ! Reset to original species index
      PartSpecies(iPart) = PartSpecies(iPart) - SpeciesOffsetVDL
      ! Check if vMPF is active
      IF(usevMPF)THEN
        ! Calculate the charge considering the MPF of the specific particle
        charge = Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
      ELSE
        ! Calculate the charge considering the MPF of the species
        charge = Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
      END IF
      ! Deposit the charge on the corner nodes of the element
      CALL DepositParticleOnNodes(SignSwitch*charge, PartState(1:3,iPart), PEM%GlobalElemID(iPart))
      ! After deposition, delete the particle from existence
      CALL RemoveParticle(iPart)
    END IF ! IsVDLSpecID(iPart)
  END IF !PDM%ParticleInside(iPart)
END DO ! iPart=1,PDM%ParticleVecLength

END SUBROUTINE DepositVirtualDielectricLayerParticles
#endif /*USE_HDG*/

#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
END MODULE MOD_PICDepo_HDG
