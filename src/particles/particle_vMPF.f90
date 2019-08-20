!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer and Asim Mirza
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

MODULE MOD_vMPF
!===================================================================================================================================
! Module controlling particle number by merge and split routines
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
PUBLIC :: SplitMerge_main
!===================================================================================================================================

CONTAINS

SUBROUTINE SplitMerge_main()
!===================================================================================================================================
!> Main routine for split and merge particles
!> Loop over all elements:
!> 1.) build partindx list
!> 2.) Call split or merge routine
!===================================================================================================================================
! MODULES
USE MOD_PARTICLE_Vars         ,ONLY: vMPFNewPartNum, PEM
USE MOD_Mesh_Vars             ,ONLY: nElems

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem, nPart, iLoop, iPart
INTEGER, ALLOCATABLE  :: iPartIndx_Node(:)
!===================================================================================================================================

DO iElem = 1, nElems
  nPart = PEM%pNumber(iElem)
  IF (nPart.LT.vMPFNewPartNum) CYCLE
! 1.) build partindx list
  ALLOCATE(iPartIndx_Node(nPart))
  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    iPartIndx_Node(iLoop) = iPart
    iPart = PEM%pNext(iPart)
  END DO
! 2.) Call split or merge routine
  CALL MergeParticles(iPartIndx_Node, nPart, vMPFNewPartNum)
  DEALLOCATE(iPartIndx_Node)
END DO

END SUBROUTINE SplitMerge_main

SUBROUTINE MergeParticles(iPartIndx_Node, nPart, nPartNew)
!===================================================================================================================================
!> Routine for merge particles
!> 1.) Calc bulkvelocity (for momentum conservation)
!> 2.) Calc energy (for energy conservation)
!> 3.) Delete particles randomly (until nPartNew is reached)
!> 4.) Calc bulkvelocity after deleting
!> 5.) Calc energy after deleting
!> 6.) Ensuring momentum and energy conservation
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PartMPF
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                  :: nPart, nPartNew
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: V_rel(3), vmag2, iRan, vBulk(3), Energy
INTEGER               :: iLoop,fillMa1, fillMa2, nDelete, nTemp, iPart, iPartIndx_NodeTMP(nPart), joa(2)
REAL                  :: partWeight, totalWeight, vBulkTmp(3), ENew, alpha
!===================================================================================================================================
vBulk = 0.0; totalWeight = 0.0; Energy = 0.

! 1.) calc bulkvelocity (for momentum conservation)
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
END DO
vBulk(1:3) = vBulk(1:3)/ totalWeight

! 2.) calc energy (for energy conservation)
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  Energy = Energy + 0.5 * vmag2 * partWeight
  ! sample inner energies here!
END DO

! 3.) delete particles randomly (until nPartNew is reached)
iPartIndx_NodeTMP = iPartIndx_Node
nTemp = nPart
nDelete = nPart - nPartNew
DO iLoop = 1, nDelete
  CALL RANDOM_NUMBER(iRan)
  iPart = INT(iRan*nTemp) + 1
  PDM%ParticleInside(iPartIndx_Node(iPart)) = .FALSE.
  iPartIndx_Node(iPart) = iPartIndx_Node(nTemp)
  nTemp = nTemp - 1
END DO

! 4.) calc bulkvelocity after deleting
vBulkTmp = 0.
DO iLoop = 1, nPartNew
  PartMPF(iPartIndx_Node(iLoop)) = totalWeight / nPartNew
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  vBulkTmp(1:3) = vBulkTmp(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
END DO
vBulkTmp(1:3) = vBulkTmp(1:3) / totalWeight

! 5.) calc energy after deleting
ENew = 0.
totalWeight=0.0
DO iLoop = 1, nPartNew
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulkTmp(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  ENew = ENew + 0.5 * vmag2 * partWeight
  ! sample inner energies here!
END DO

! 6.) ensuring momentum and energy conservation
alpha = SQRT(Energy/ENew)
DO iLoop = 1, nPartNew
  PartState(iPartIndx_Node(iLoop),4:6) = vBulk(1:3) + alpha*(PartState(iPartIndx_Node(iLoop),4:6)-vBulkTmp(1:3))
END DO

END SUBROUTINE MergeParticles

SUBROUTINE CalculateDistMoments(iPartIndx_Node, nPart, vBulk, Vtherm2, PressTens, HeatVec, Energy)
!===================================================================================================================================
!> Calculation of distribution moments
!> 1.) Calc bulkvelocity
!> 2.) Summing up the relative velocities and their square to calculate the moments (PressTens, HeatVec)
!> 3.) Fill missing entries in PressTens
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars         ,ONLY: PartState
USE MOD_part_tools            ,ONLY: GetParticleWeight
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT)                  :: nPart
INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
REAL, INTENT(INOUT)                       :: vBulk(3), Energy, Vtherm2, PressTens(3,3), HeatVec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: V_rel(3), vmag2
INTEGER               :: iLoop,fillMa1, fillMa2
REAL                  :: partWeight, totalWeight
!===================================================================================================================================
Vtherm2 = 0.0; PressTens = 0.0; HeatVec = 0.0
vBulk = 0.0; totalWeight = 0.0; Energy = 0.

! 1.) calc bulkvelocity
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
END DO
vBulk(1:3) = vBulk(1:3)/ totalWeight

! 2.) Summing up the relative velocities and their square to calculate the moments (PressTens, HeatVec)
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  Vtherm2 = Vtherm2 + vmag2 * partWeight
  DO fillMa1 =1, 3
    DO fillMa2 =fillMa1, 3
      PressTens(fillMa1, fillMa2)= PressTens(fillMa1, fillMa2) + V_rel(fillMa1)*V_rel(fillMa2) * partWeight
    END DO
  END DO
  HeatVec(1:3) = HeatVec(1:3) + V_rel(1:3)*vmag2 * partWeight
  Energy = Energy + 0.5 * vmag2 * partWeight
  ! sample inner energies here!
END DO
IF(nPart.GT.2) THEN
  HeatVec = HeatVec*nPart*nPart/((nPart-1.)*(nPart-2.)*totalWeight)
ELSE
  HeatVec = 0.0
END IF
Vtherm2 = Vtherm2*nPart/((nPart-1.)*totalWeight)
! 3.) Fill missing entries in PressTens
PressTens(2,1)=PressTens(1,2)
PressTens(3,1)=PressTens(1,3)
PressTens(3,2)=PressTens(2,3)
PressTens = PressTens/totalWeight

END SUBROUTINE CalculateDistMoments

SUBROUTINE ARChapEnsk(nPart, iRanPart, Vtherm, HeatVec, PressTens)
!===================================================================================================================================
!> Acceptance rejection method for reconstruct distribution moments (according Chapman Enskog)
!===================================================================================================================================
! MODULES
USE Ziggurat
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: nPart
REAL, INTENT(IN)              :: HeatVec(3), Vtherm, PressTens(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: iRanPart(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Vheat, V2, iRan, OldProb, Envelope, Envelope2, cMat, cPress
INTEGER                        :: iPart, fillMa1, fillMa2
!===================================================================================================================================
Envelope = MAX(ABS(HeatVec(1)),ABS(HeatVec(2)),ABS(HeatVec(3)))/Vtherm**(3./2.)
Envelope2 = MAX(ABS(PressTens(1,2)),ABS(PressTens(1,3)),ABS(PressTens(2,3)))/Vtherm
Envelope =  1.+4.*MAX(Envelope, Envelope2)

DO iPart = 1, nPart
  iRanPart(1,iPart) = rnor()
  iRanPart(2,iPart) = rnor()
  iRanPart(3,iPart) = rnor()
  cMat = 0.0
  cPress = 0.0
  cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
  cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
  cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
  cPress=cPress + (PressTens(1,1)-Vtherm)*(iRanPart(1,iPart)*iRanPart(1,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
  cPress=cPress + (PressTens(2,2)-Vtherm)*(iRanPart(2,iPart)*iRanPart(2,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
  V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
  Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
  OldProb =  (1. + cMat/Vtherm + cPress/(2.*Vtherm) + VHeat/(2.*Vtherm**(3./2.))*(V2/5.-1.))
  CALL RANDOM_NUMBER(iRan)
  DO WHILE (Envelope*iRan.GT.OldProb)
    iRanPart(1,iPart) = rnor()
    iRanPart(2,iPart) = rnor()
    iRanPart(3,iPart) = rnor()
    cMat = 0.0
    cPress = 0.0
    cMat=cMat + iRanPart(1,iPart)*iRanPart(2,iPart)*PressTens(1,2)
    cMat=cMat + iRanPart(1,iPart)*iRanPart(3,iPart)*PressTens(1,3)
    cMat=cMat + iRanPart(2,iPart)*iRanPart(3,iPart)*PressTens(2,3)
    cPress=cPress + (PressTens(1,1)-Vtherm)*(iRanPart(1,iPart)*iRanPart(1,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
    cPress=cPress + (PressTens(2,2)-Vtherm)*(iRanPart(2,iPart)*iRanPart(2,iPart)-iRanPart(3,iPart)*iRanPart(3,iPart))
    V2 = iRanPart(1,iPart)*iRanPart(1,iPart) + iRanPart(2,iPart)*iRanPart(2,iPart) + iRanPart(3,iPart)*iRanPart(3,iPart)
    Vheat = iRanPart(1,iPart)*HeatVec(1) + iRanPart(2,iPart)*HeatVec(2) + iRanPart(3,iPart)*HeatVec(3)
    OldProb =  (1. + cMat/Vtherm + cPress/(2.*Vtherm) + VHeat/(2.*Vtherm**(3./2.))*(V2/5.-1.))
    CALL RANDOM_NUMBER(iRan)
  END DO
END DO

END SUBROUTINE ARChapEnsk


END MODULE MOD_vMPF
