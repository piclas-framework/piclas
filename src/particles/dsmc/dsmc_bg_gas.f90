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

MODULE MOD_DSMC_BGGas
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_InitBGGas
  MODULE PROCEDURE DSMC_InitBGGas
END INTERFACE

INTERFACE DSMC_pairing_bggas
  MODULE PROCEDURE DSMC_pairing_bggas
END INTERFACE

INTERFACE DSMC_FinalizeBGGas
  MODULE PROCEDURE DSMC_FinalizeBGGas
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_InitBGGas, DSMC_pairing_bggas, DSMC_FinalizeBGGas
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_InitBGGas()
!===================================================================================================================================
! Initialization of background gas
!===================================================================================================================================
! MODULES
  USE MOD_Globals,            ONLY : Abort
  USE MOD_DSMC_Init,          ONLY : DSMC_SetInternalEnr_LauxVFD
  USE MOD_DSMC_Vars,          ONLY : BGGas, SpecDSMC
  USE MOD_DSMC_PolyAtomicModel, ONLY : DSMC_SetInternalEnr_Poly
  USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies, PartState, PEM, PartPosRef
  USE MOD_part_emission,      ONLY : SetParticleChargeAndMass, SetParticleVelocity, SetParticleMPF
  USE MOD_part_tools,         ONLY : UpdateNextFreePosition
  USE MOD_Particle_Tracking_Vars, ONLY:DoRefmapping
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER           :: iNewPart, iPart, PositionNbr
!===================================================================================================================================
  iNewPart=0
  PositionNbr = 0
  DO iPart = 1, PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart).AND.(PartSpecies(iPart).NE.BGGas%BGGasSpecies)) THEN
      iNewPart = iNewPart + 1
      PositionNbr = PDM%nextFreePosition(iNewPart+PDM%CurrentNextFreePosition)
      IF (PositionNbr.EQ.0) THEN
        CALL Abort(&
__STAMP__&
,'ERROR in BGGas: Too many Particles!')
      END IF
      PartState(PositionNbr,1:3) = PartState(iPart,1:3)
      IF(DoRefMapping)THEN ! here Nearst-GP is missing
        PartPosRef(1:3,PositionNbr)=PartPosRef(1:3,iPart)
      END IF
      PartSpecies(PositionNbr) = BGGas%BGGasSpecies
      IF(SpecDSMC(BGGas%BGGasSpecies)%PolyatomicMol) THEN
        CALL DSMC_SetInternalEnr_Poly(BGGas%BGGasSpecies,0,PositionNbr,1)
      ELSE
        CALL DSMC_SetInternalEnr_LauxVFD(BGGas%BGGasSpecies,0,PositionNbr,1)
      END IF
      PEM%Element(PositionNbr) = PEM%Element(iPart)
      PDM%ParticleInside(PositionNbr) = .true.
      PEM%pNext(PEM%pEnd(PEM%Element(PositionNbr))) = PositionNbr     ! Next Particle of same Elem (Linked List)
      PEM%pEnd(PEM%Element(PositionNbr)) = PositionNbr
      PEM%pNumber(PEM%Element(PositionNbr)) = &                       ! Number of Particles in Element
      PEM%pNumber(PEM%Element(PositionNbr)) + 1
    END IF
  END DO
  CALL SetParticleVelocity(BGGas%BGGasSpecies,0,iNewPart,1) ! Properties of BG gas are stored in iInit=0
  PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,PositionNbr)
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + iNewPart 

END SUBROUTINE DSMC_InitBGGas

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_pairing_bggas(iElem)
!===================================================================================================================================
! Building of pairs for the background gas
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, BGGas, CollisMode, ChemReac, PartStateIntEn
  USE MOD_Particle_Vars,          ONLY : PEM,PartSpecies,nSpecies,PartState,Species,usevMPF,PartMPF,Species
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem          
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                       !
  INTEGER                       :: nPair, iPair, iPart, iLoop, cPart1, cPart2, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase
  INTEGER, ALLOCATABLE          :: iPartIndx(:) ! List of particles in the cell nec for stat pairing
!===================================================================================================================================
  nPart = PEM%pNumber(iElem)
  nPair = INT(nPart/2)
  
  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  ALLOCATE(iPartIndx(nPart))
  Coll_pData%Ec=0
  iPartIndx = 0
  IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iPartIndx(iLoop) = iPart
    CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + 1 
          ! counter for part num of spec per cell   
    IF (CollisMode.EQ.3) ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) = &
                        ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) &
                        + PartStateIntEn(iPart,1) !Calulation of mean evib per cell and iter, necessary for disso prob 
    iPart = PEM%pNext(iPart)    
  END DO
  
  ! Setting Number of BGGas Particles per Cell
  IF (Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensityFileID.GT.0) THEN
    BGGas%BGColl_SpecPartNum = Species(BGGas%BGGasSpecies)%Init(0)%ElemPartDensity(iElem) * GEO%Volume(iElem)      &
                                               / Species(BGGas%BGGasSpecies)%MacroParticleFactor
  ELSE
    BGGas%BGColl_SpecPartNum = BGGas%BGGasDensity * GEO%Volume(iElem)      &
                                               / Species(BGGas%BGGasSpecies)%MacroParticleFactor
  END IF

  CollInf%Coll_SpecPartNum(BGGas%BGGasSpecies) = BGGas%BGColl_SpecPartNum

  DO iPair = 1, nPair 
    DO cPart1 = 1, nPart ! Searching the first non BG particle
      IF (PartSpecies(iPartIndx(cPart1)).NE.BGGas%BGGasSpecies) THEN
        Coll_pData(iPair)%iPart_p1 = iPartIndx(cPart1)
        iPartIndx(cPart1) = iPartIndx(nPart)
        EXIT
      END IF
    END DO
    nPart = nPart - 1

    DO cPart2 = 1, nPart ! Searching the first BG particle
      IF (PartSpecies(iPartIndx(cPart2)).EQ.BGGas%BGGasSpecies) THEN
        Coll_pData(iPair)%iPart_p2 = iPartIndx(cPart2)
        iPartIndx(cPart2) = iPartIndx(nPart)
        EXIT
      END IF
    END DO
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
    IF (usevMPF) PartMPF(Coll_pData(iPair)%iPart_p2) = PartMPF(Coll_pData(iPair)%iPart_p1) 
    iCase = CollInf%Coll_Case(cSpec1, cSpec2) 
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2 
    Coll_pData(iPair)%PairType = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
  END DO
  
  DEALLOCATE(iPartIndx)
END SUBROUTINE DSMC_pairing_bggas


SUBROUTINE DSMC_FinalizeBGGas()
!===================================================================================================================================
! Kills all BG Particles
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,          ONLY : BGGas
USE MOD_PARTICLE_Vars,      ONLY : PDM, PartSpecies
USE MOD_part_tools,         ONLY : UpdateNextFreePosition
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  WHERE (PartSpecies(1:PDM%ParticleVecLength).EQ.BGGas%BGGasSpecies) &
         PDM%ParticleInside(1:PDM%ParticleVecLength) = .FALSE.
  CALL UpdateNextFreePosition()

END SUBROUTINE DSMC_FinalizeBGGas

END MODULE MOD_DSMC_BGGas
