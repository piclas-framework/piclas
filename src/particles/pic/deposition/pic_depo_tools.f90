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

MODULE  MOD_PICDepo_Tools
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
INTERFACE DepositParticleOnNodes
MODULE PROCEDURE DepositParticleOnNodes
END INTERFACE

PUBLIC:: DepositParticleOnNodes
!===================================================================================================================================

CONTAINS


SUBROUTINE DepositParticleOnNodes(Charge,PartPos,ElemID) 
!----------------------------------------------------------------------------------------------------------------------------------!
! Deposit the charge of a single particle on the nodes corresponding to the deposition method 'cell_volweight_mean', where the
! charge density is stored in NodeSourceExt, which is added to NodeSource in the standard deposition procedure.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PICDepo_Vars             ,ONLY: NodeSourceExt
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO
#if ((USE_HDG) && (PP_nVar==1))
USE MOD_TimeDisc_Vars            ,ONLY: dt,tAnalyzeDiff,tEndDiff
#endif
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBElemPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Tracking_Vars   ,ONLY: TrackInfo
USE MOD_Particle_Vars            ,ONLY: LastPartPos,PartSpecies,Species
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)     :: Charge
REAL,INTENT(IN)     :: PartPos(1:3)
INTEGER,INTENT(IN)  :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3)
REAL                             :: TSource(1:4)
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER               :: SourceDim=1
LOGICAL, PARAMETER               :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL                          :: doCalculateCurrentDensity
INTEGER                          :: SourceDim
#endif
!===================================================================================================================================
! Check whether charge and current density have to be computed or just the charge density
#if ((USE_HDG) && (PP_nVar==1))
IF(ALMOSTEQUAL(dt,tAnalyzeDiff).OR.ALMOSTEQUAL(dt,tEndDiff))THEN
doCalculateCurrentDensity=.TRUE.
SourceDim=1
ELSE ! do not calculate current density
doCalculateCurrentDensity=.FALSE.
SourceDim=4
END IF
#endif
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

CALL GetPositionInRefElem(PartPos,TempPartPos(1:3),ElemID,ForceMode=.TRUE.)
TSource(:) = 0.0
IF(doCalculateCurrentDensity)THEN
TSource(1:3) = PartPos*Charge
END IF
TSource(4) = Charge

alpha1=0.5*(TempPartPos(1)+1.0)
alpha2=0.5*(TempPartPos(2)+1.0)
alpha3=0.5*(TempPartPos(3)+1.0)
ASSOCIATE( NodeID     => GEO%ElemToNodeID(:,ElemID) )
NodeSourceExt(:,NodeID(1)) = NodeSourceExt(:,NodeID(1))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(2)) = NodeSourceExt(:,NodeID(2))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(3)) = NodeSourceExt(:,NodeID(3))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(4)) = NodeSourceExt(:,NodeID(4))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(5)) = NodeSourceExt(:,NodeID(5))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(6)) = NodeSourceExt(:,NodeID(6))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(7)) = NodeSourceExt(:,NodeID(7))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(8)) = NodeSourceExt(:,NodeID(8))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*  (alpha3))
END ASSOCIATE
#if USE_LOADBALANCE
CALL LBElemPauseTime(ElemID,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to ElemID
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositParticleOnNodes


END MODULE MOD_PICDepo_Tools
