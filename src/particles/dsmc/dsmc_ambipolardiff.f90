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

MODULE MOD_DSMC_AmbipolarDiffusion
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE AD_InsertParticles
  MODULE PROCEDURE AD_InsertParticles
END INTERFACE

!INTERFACE BGGas_DeleteParticles
!  MODULE PROCEDURE BGGas_DeleteParticles
!END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: AD_InsertParticles,
!===================================================================================================================================

CONTAINS

SUBROUTINE AD_InsertParticles(iPartIndx_Node,nPart, vBulk)
!===================================================================================================================================
!> Creating electrons for each actual ion simulation particle
!===================================================================================================================================
! MODULES
USE MOD_Globals                ,ONLY: Abort
USE MOD_DSMC_Init              ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_Vars              ,ONLY: BGGas, SpecDSMC, CollisMode
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_PARTICLE_Vars          ,ONLY: PDM, PartSpecies, PartState, PEM, PartPosRef
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF,CalcVelocity_maxwell_lpn
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefmapping
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr, iSpec, LocalElemID, iPartIndx_Nodetmp(nPart)
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

iNewPart=0
PositionNbr = 0
nNewElectrons = 0
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) nNewElectrons = nNewElectrons + 1
END DO
IF (nNewElectrons.EQ.0) RETURN
iPartIndx_Nodetmp = iPartIndx_Node
DEALLOCATE(iPartIndx_Node)
ALLOCATE(iPartIndx_Node(nPart + nNewElectrons))
iPartIndx_Node(1:nPart) = iPartIndx_Nodetmp(1:nPart)

!DO iPart = 1, PDM%ParticleVecLength
!  IF (PDM%ParticleInside(iPart)) THEN
!    ! Skip background particles that have been created within this loop
!    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) CYCLE
!    iNewPart = iNewPart + 1
!    PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
!    IF (PositionNbr.EQ.0) THEN
!      CALL Abort(&
!__STAMP__&
!,'ERROR in BGGas: MaxParticleNumber should be twice the expected number of particles, to account for the BGG particles!')
!    END IF
!    PartState(1:3,PositionNbr) = PartState(1:3,iPart)
!    IF(DoRefMapping)THEN ! here Nearst-GP is missing
!      PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
!    END IF
!    iSpec = BGGas_GetSpecies()
!    PartSpecies(PositionNbr) = iSpec
!    IF(CollisMode.GT.1) THEN
!      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
!        CALL DSMC_SetInternalEnr_Poly(iSpec,1,PositionNbr,1)
!      ELSE
!        CALL DSMC_SetInternalEnr_LauxVFD(iSpec,1,PositionNbr,1)
!      END IF
!    END IF
!    PEM%GlobalElemID(PositionNbr) = PEM%GlobalElemID(iPart)
!    LocalElemID = PEM%LocalElemID(PositionNbr)
!    PDM%ParticleInside(PositionNbr) = .true.
!    PEM%pNext(PEM%pEnd(LocalElemID)) = PositionNbr     ! Next Particle of same Elem (Linked List)
!    PEM%pEnd(LocalElemID) = PositionNbr
!    PEM%pNumber(LocalElemID) = PEM%pNumber(LocalElemID) + 1
!    BGGas%PairingPartner(iPart) = PositionNbr
!    CALL CalcVelocity_maxwell_lpn(FractNbr=iSpec, Vec3D=PartState(4:6,PositionNbr), iInit=1)
!  END IF
!END DO
!PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
!PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DSMC,tLBStart)
#endif /*USE_LOADBALANCE*/
END SUBROUTINE AD_InsertParticles



!SUBROUTINE BGGas_DeleteParticles()
!!===================================================================================================================================
!! Deletes all background gas particles and updates the particle index list
!!===================================================================================================================================
!! MODULES
!USE MOD_DSMC_Vars,          ONLY : BGGas
!USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies
!USE MOD_part_tools,         ONLY : UpdateNextFreePosition
!! IMPLICIT VARIABLE HANDLING
!  IMPLICIT NONE
!!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                     :: iPart
!!===================================================================================================================================

!DO iPart = 1, PDM%ParticleVecLength
!  IF (PDM%ParticleInside(iPart)) THEN
!    IF(BGGas%BackgroundSpecies(PartSpecies(iPart))) PDM%ParticleInside(iPart) = .FALSE.
!  END IF
!END DO
!BGGas%PairingPartner = 0
!CALL UpdateNextFreePosition()

!END SUBROUTINE BGGas_DeleteParticles

END MODULE MOD_DSMC_AmbipolarDiffusion
