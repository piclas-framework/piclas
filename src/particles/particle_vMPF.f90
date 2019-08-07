!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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
! Module approximating the collision operator using the Bhatnagar-Gross-Krook model
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
PUBLIC :: 
!===================================================================================================================================

CONTAINS


SUBROUTINE MergeParticles(iPartIndx_Node, nPart, nPartNew)
!===================================================================================================================================
!
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
REAL, INTENT(OUT)                       :: vBulk(3), Energy, Vtherm2, PressTens(3,3), HeatVec(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: V_rel(3), vmag2, iRan
INTEGER               :: iLoop,fillMa1, fillMa2, nDelete, nTemp, iPart
REAL                  :: partWeight, totalWeight
!===================================================================================================================================
Vtherm2 = 0.0; PressTens = 0.0; HeatVec = 0.0
vBulk = 0.0; totalWeight = 0.0; Energy = 0.
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
END DO
vBulk(1:3) = vBulk(1:3)/ totalWeight

! 1.) Summing up the relative velocities and their square to calculate the moments
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  V_rel(1:3)=PartState(iPartIndx_Node(iLoop),4:6)-vBulk(1:3)
  vmag2 = V_rel(1)**2 + V_rel(2)**2 + V_rel(3)**2
  Energy = Energy + 0.5 * vmag2 * partWeight
  ! sample inner energies here!
END DO

nTemp = nPart
nDelete = nPart - nPartNew
DO iLoop = 1, nDelete
  CALL RANDOM_NUMBER(iRan)
  iPart = INT(iRan*nTemp) + 1
  PDM%ParticleInside(iPartIndx_Node(iPart)) = .FALSE.
  iPartIndx_Node(iPart) = iPartIndx_Node(nTemp)
  nTemp = nTemp - 1
END DO

END SUBROUTINE MergeParticles

SUBROUTINE CalculateDistMomements(iPartIndx_Node, nPart, vBulk, Vtherm2, PressTens, HeatVec, Energy)
!===================================================================================================================================
!
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
REAL, INTENT(OUT)                       :: vBulk(3), Energy, Vtherm2, PressTens(3,3), HeatVec(3)
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
DO iLoop = 1, nPart
  partWeight = GetParticleWeight(iPartIndx_Node(iLoop))
  totalWeight = totalWeight + partWeight
  vBulk(1:3) = vBulk(1:3) + PartState(iPartIndx_Node(iLoop),4:6) * partWeight
END DO
vBulk(1:3) = vBulk(1:3)/ totalWeight

! 1.) Summing up the relative velocities and their square to calculate the moments
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
PressTens(2,1)=PressTens(1,2)
PressTens(3,1)=PressTens(1,3)
PressTens(3,2)=PressTens(2,3)
PressTens = PressTens/totalWeight

END SUBROUTINE CalculateDistMomements

SUBROUTINE ARChapEnsk(nPart, iRanPart, Vtherm, HeatVec, PressTens)
!===================================================================================================================================
!> description
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
