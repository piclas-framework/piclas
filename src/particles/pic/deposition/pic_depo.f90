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

MODULE  MOD_PICDepo
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
 IMPLICIT NONE
 PRIVATE
!===================================================================================================================================
INTERFACE Deposition
  MODULE PROCEDURE Deposition
END INTERFACE

INTERFACE Deboor
  MODULE PROCEDURE Deboor
END INTERFACE

PUBLIC::Deposition
PUBLIC::Deboor
!===================================================================================================================================

CONTAINS

SUBROUTINE Deposition(doInnerParts,doParticle_In)
!===================================================================================================================================
! This subroutine performes the deposition of the particle charge and current density to the grid
! following list of distribution methods are implemted
! - nearest blurrycenter (barycenter of hexahedra)
! - nearest Gauss Point  (only volome of IP - higher resolution than nearest blurrycenter )
! - shape function       (only one type implemented)
! - delta distributio
! useVMPF added, therefore, this routine contains automatically the use of variable mpfs
!===================================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Particle_Vars
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY:PI
USE MOD_Mesh_Vars,              ONLY:nElems, Elem_xGP, sJ, nNodes
USE MOD_ChangeBasis,            ONLY:ChangeBasis3D
USE MOD_Interpolation_Vars,     ONLY:wGP
USE MOD_PICInterpolation_Vars,  ONLY:InterpolationType
USE MOD_Eval_xyz,               ONLY:GetPositionInRefElem
USE MOD_Basis,                  ONLY:LagrangeInterpolationPolys,BernSteinPolynomial
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,casematrix, NbrOfCases
!USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNGeo,ElemRadiusNGeo,ElemRadius2NGeo
#ifdef MPI
! only required for shape function??
USE MOD_Particle_MPI_Vars,      ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,ExtPartToFIBGM,NbrOfExtParticles
USE MOD_Particle_MPI_Vars,      ONLY:PartMPIExchange
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
USE MOD_Particle_MPI,           ONLY:AddHaloNodeData
#endif  /*MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools,      ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
#endif /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT variable declaration
LOGICAL,INTENT(IN)               :: doInnerParts
LOGICAL,INTENT(IN),OPTIONAL      :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                          :: doPartInExists
INTEGER                          :: firstPart,lastPart
INTEGER                          :: i,j, k, l, m, iElem, iPart, iPart2, iSFfix
LOGICAL                          :: chargedone(1:nElems)!, SAVE_GAUSS
LOGICAL                          :: SAVE_GAUSS
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: ElemID, iCase, ind
REAL                             :: radius2, S, S1, Fac(4)!, Fac2(4)
REAL                             :: dx,dy,dz
!REAL                             :: GaussDistance(0:PP_N,0:PP_N,0:PP_N)
REAL, ALLOCATABLE                :: BGMSourceCellVol(:,:,:,:,:), tempsource(:,:,:), tempgridsource(:)
REAL, ALLOCATABLE                :: NodeSource(:,:), tempNodeSource(:,:)
REAL                             :: Vec1(1:3), Vec2(1:3), Vec3(1:3), ShiftedPart(1:3)!, caseShiftedPart(1:3)
INTEGER                          :: a,b, ii, expo
REAL                             :: ElemSource(nElems,1:4)
REAL                             :: Charge, TSource(1:4), auxiliary(0:3),weight(1:3,0:3), locweight
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3), alpha
INTEGER                          :: PosInd(3),r,ss,t,u,v,w, dir, weightrun
INTEGER                          :: iLayer, layerParts
!LOGICAL                          :: DoCycle
REAL,DIMENSION(3,0:NDepo)        :: L_xi
REAL                             :: DeltaIntCoeff,prefac!, SFfixDistance
REAL                             :: local_r_sf, local_r2_sf, local_r2_sf_inv
REAL                             :: RandVal, RandVal2(2), layerPartPos(3), PartRadius, FractPush(3), SFfixDistance
LOGICAL                          :: DoCycle
#if USE_LOADBALANCE
REAL                             :: tLBStart ! load balance
#endif /*USE_LOADBALANCE*/
!============================================================================================================================
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.

IF(doInnerParts)THEN
  IF(.NOT.DoDeposition) RETURN
  PartSource=0.0
  firstPart=1
  lastPart =PDM%ParticleVecLength
  !IF(firstPart.GT.lastPart) RETURN
ELSE
  IF(.NOT.DoDeposition) RETURN
#ifdef MPI
  firstPart=PDM%ParticleVecLength-PartMPIExchange%nMPIParticles+1
  lastPart =PDM%ParticleVecLength
#else
  firstPart=1
  lastPart =0
#endif /*MPI*/
END IF

!IF((firstPart.GT.lastPart).AND.(DepositionType.NE.'delta_distri').AND.(DepositionType.NE.'shape_function')&
!                          .AND.(DepositionType.NE.'nearest_blurrycenter')) RETURN

SELECT CASE(TRIM(DepositionType))
CASE('nearest_blurrycenter')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  ElemSource=0.0
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    DO iPart=firstPart,lastPart
      IF(doPartInExists)THEN
        IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
      ELSE
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
      END IF
      IF(PEM%Element(iPart).EQ.iElem)THEN
        IF(usevMPF)THEN
!#if (PP_nVar==8)
         ElemSource(iElem,1:3) = ElemSource(iElem,1:3)+ &
                PartState(iPart,4:6)* Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
!#endif
         ElemSource(iElem,4) = ElemSource(iElem,4) + &
              Species(PartSpecies(iPart))%ChargeIC* PartMPF(iPart)
        ELSE
!#if (PP_nVar==8)
         ElemSource(iElem,1:3) = ElemSource(iElem,1:3)+ &
                PartState(iPart,4:6)* Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
!#endif
         ElemSource(iElem,4) = ElemSource(iElem,4) + &
              Species(PartSpecies(iPart))%ChargeIC* Species(PartSpecies(iPart))%MacroParticleFactor
        END IF ! usevMPF
      END IF ! Element(iPart).EQ.iElem
    END DO ! iPart
!#if (PP_nVar==8)
    PartSource(1,:,:,:,iElem) = PartSource(1,:,:,:,iElem)+ElemSource(iElem,1)
    PartSource(2,:,:,:,iElem) = PartSource(2,:,:,:,iElem)+ElemSource(iElem,2)
    PartSource(3,:,:,:,iElem) = PartSource(3,:,:,:,iElem)+ElemSource(iElem,3)
!#endif
    PartSource(4,:,:,:,iElem) = PartSource(4,:,:,:,iElem)+ElemSource(iElem,4)
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem=1,PP_nElems
  IF(.NOT.doInnerParts)THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    DO iElem=1,PP_nElems
!#if (PP_nVar==8)
      PartSource(1:4,:,:,:,iElem) = PartSource(1:4,:,:,:,iElem) / GEO%Volume(iElem)
!#else
!      PartSource(4,:,:,:,iElem) = PartSource(4,:,:,:,iElem) / GEO%Volume(iElem)
!#endif
#if USE_LOADBALANCE
      CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! iElem=1,PP_nElems
  END IF ! .NOT. doInnerParts
CASE('cell_volweight')
  ALLOCATE(BGMSourceCellVol(0:1,0:1,0:1,1:nElems,1:4))
  BGMSourceCellVol(:,:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
  IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    IF (usevMPF) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
    ELSE
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF ! usevMPF
    iElem = PEM%Element(iPart)
    IF(DoRefMapping)THEN
      TempPartPos(1:3)=PartPosRef(1:3,iPart)
    ELSE
      CALL GetPositionInRefElem(PartState(iPart,1:3),TempPartPos,iElem,ForceMode=.TRUE.)
    END IF
    TSource(:) = 0.0
!#if (PP_nVar==8)
    TSource(1) = PartState(iPart,4)*Charge
    TSource(2) = PartState(iPart,5)*Charge
    TSource(3) = PartState(iPart,6)*Charge
!#endif
    TSource(4) = Charge
    alpha1=(TempPartPos(1)+1.0)/2.0
    alpha2=(TempPartPos(2)+1.0)/2.0
    alpha3=(TempPartPos(3)+1.0)/2.0
    BGMSourceCellVol(0,0,0,iElem,1:4) = BGMSourceCellVol(0,0,0,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
    BGMSourceCellVol(0,0,1,iElem,1:4) = BGMSourceCellVol(0,0,1,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(alpha3))
    BGMSourceCellVol(0,1,0,iElem,1:4) = BGMSourceCellVol(0,1,0,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(1-alpha3))
    BGMSourceCellVol(0,1,1,iElem,1:4) = BGMSourceCellVol(0,1,1,iElem,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(alpha3))
    BGMSourceCellVol(1,0,0,iElem,1:4) = BGMSourceCellVol(1,0,0,iElem,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(1-alpha3))
    BGMSourceCellVol(1,0,1,iElem,1:4) = BGMSourceCellVol(1,0,1,iElem,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(alpha3))
    BGMSourceCellVol(1,1,0,iElem,1:4) = BGMSourceCellVol(1,1,0,iElem,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(1-alpha3))
    BGMSourceCellVol(1,1,1,iElem,1:4) = BGMSourceCellVol(1,1,1,iElem,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(alpha3))
#if USE_LOADBALANCE
    CALL LBElemPauseTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1, nElems
    BGMSourceCellVol(0,0,0,iElem,:) = BGMSourceCellVol(0,0,0,iElem,1:4)/CellVolWeight_Volumes(0,0,0,iElem)
    BGMSourceCellVol(0,0,1,iElem,:) = BGMSourceCellVol(0,0,1,iElem,1:4)/CellVolWeight_Volumes(0,0,1,iElem)
    BGMSourceCellVol(0,1,0,iElem,:) = BGMSourceCellVol(0,1,0,iElem,1:4)/CellVolWeight_Volumes(0,1,0,iElem)
    BGMSourceCellVol(0,1,1,iElem,:) = BGMSourceCellVol(0,1,1,iElem,1:4)/CellVolWeight_Volumes(0,1,1,iElem)
    BGMSourceCellVol(1,0,0,iElem,:) = BGMSourceCellVol(1,0,0,iElem,1:4)/CellVolWeight_Volumes(1,0,0,iElem)
    BGMSourceCellVol(1,0,1,iElem,:) = BGMSourceCellVol(1,0,1,iElem,1:4)/CellVolWeight_Volumes(1,0,1,iElem)
    BGMSourceCellVol(1,1,0,iElem,:) = BGMSourceCellVol(1,1,0,iElem,1:4)/CellVolWeight_Volumes(1,1,0,iElem)
    BGMSourceCellVol(1,1,1,iElem,:) = BGMSourceCellVol(1,1,1,iElem,1:4)/CellVolWeight_Volumes(1,1,1,iElem)
  END DO
#if USE_LOADBALANCE
  CALL LBElemPauseTime_avg(tLBStart) ! average over the number of elems
#endif /*USE_LOADBALANCE*/

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         alpha1 = CellVolWeightFac(kk)
         alpha2 = CellVolWeightFac(ll)
         alpha3 = CellVolWeightFac(mm)
         PartSource(1:4,kk,ll,mm,iElem) =PartSource(1:4,kk,ll,mm,iElem) +&
              BGMSourceCellVol(0,0,0,iElem,1:4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSourceCellVol(0,0,1,iElem,1:4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSourceCellVol(0,1,0,iElem,1:4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSourceCellVol(0,1,1,iElem,1:4) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSourceCellVol(1,0,0,iElem,1:4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSourceCellVol(1,0,1,iElem,1:4) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSourceCellVol(1,1,0,iElem,1:4) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSourceCellVol(1,1,1,iElem,1:4) * (alpha1) * (alpha2) * (alpha3)
       END DO !mm
     END DO !ll
   END DO !kk
#if USE_LOADBALANCE
   CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
 END DO !iEle
 DEALLOCATE(BGMSourceCellVol)
CASE('cell_volweight_mean','cell_volweight_mean2')
  ALLOCATE(NodeSource(1:4,1:nNodes))
  NodeSource = 0.0

  DO iPart=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(iPart)) THEN
!#if (PP_TimeDiscMethod==440) || (PP_TimeDiscMethod==441) || (PP_TimeDiscMethod==442) || (PP_TimeDiscMethod==443) || (PP_TimeDiscMethod==445)
!      Charge = Species(PartSpecies(iPart))%MacroParticleFactor
!#else
      IF (usevMPF) THEN
        Charge = Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
      ELSE
        Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
      END IF
!#endif
      iElem = PEM%Element(iPart)
      CALL GetPositionInRefElem(PartState(iPart,1:3),TempPartPos(1:3),iElem,ForceMode=.TRUE.)
      !CALL GeoCoordToMap(PartState(iPart,1:3), TempPartPos(1:3), iElem)
      TSource(:) = 0.0
!#if (PP_TimeDiscMethod==440) || (PP_TimeDiscMethod==441) || (PP_TimeDiscMethod==442) || (PP_TimeDiscMethod==443) || (PP_TimeDiscMethod==445)
!      IF (PartSpecies(iPart).EQ.1) THEN
!        TSource(4) = Charge
!      ELSE
!        TSource(1) = Charge
!      END IF
!#else
#if !(defined (PP_HDG) && (PP_nVar==1))
      TSource(1) = PartState(iPart,4)*Charge
      TSource(2) = PartState(iPart,5)*Charge
      TSource(3) = PartState(iPart,6)*Charge
#endif
      TSource(4) = Charge
!#endif

      alpha1=(TempPartPos(1)+1.0)/2.0
      alpha2=(TempPartPos(2)+1.0)/2.0
      alpha3=(TempPartPos(3)+1.0)/2.0
      NodeSource(1:4,GEO%ElemToNodeID(1,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(1,iElem)) &
        + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(2,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(2,iElem)) &
        + (TSource(1:4)*(alpha1)*(1-alpha2)*(1-alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(3,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(3,iElem)) &
        + (TSource(1:4)*(alpha1)*(alpha2)*(1-alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(4,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(4,iElem)) &
        + (TSource(1:4)*(1-alpha1)*(alpha2)*(1-alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(5,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(5,iElem)) &
        + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(6,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(6,iElem)) &
        + (TSource(1:4)*(alpha1)*(1-alpha2)*(alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(7,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(7,iElem)) &
        + (TSource(1:4)*(alpha1)*(alpha2)*(alpha3))
      NodeSource(1:4,GEO%ElemToNodeID(8,iElem)) = NodeSource(1:4,GEO%ElemToNodeID(8,iElem)) &
        + (TSource(1:4)*(1-alpha1)*(alpha2)*(alpha3))
    END IF
  END DO
#ifdef MPI
  CALL AddHaloNodeData(NodeSource(1,:))
  CALL AddHaloNodeData(NodeSource(2,:))
  CALL AddHaloNodeData(NodeSource(3,:))
  CALL AddHaloNodeData(NodeSource(4,:))
#endif /*MPI*/

  DO iElem=1, nNodes
    NodeSource(1:4,iElem) = NodeSource(1:4,iElem)/CellLocNodes_Volumes(iElem)
  END DO

  IF (TRIM(DepositionType).EQ.'cell_volweight_mean2') THEN
    ALLOCATE(tempNodeSource(1:4,1:nNodes))
    tempNodeSource = 0.0
    DO iElem=1, nNodes
      tempNodeSource(1:4,iElem) = NodeSource(1:4,iElem)
      DO kk =1, GEO%NeighNodesOnNode(iElem)
        tempNodeSource(1:4,iElem) = tempNodeSource(1:4,iElem) + NodeSource(1:4,GEO%NodeToNeighNode(iElem)%ElemID(kk))
      END DO
      tempNodeSource(1:4,iElem) = tempNodeSource(1:4,iElem) / (GEO%NeighNodesOnNode(iElem) + 1.0)
    END DO
    NodeSource = tempNodeSource
  END IF


  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         alpha1 = CellVolWeightFac(kk)
         alpha2 = CellVolWeightFac(ll)
         alpha3 = CellVolWeightFac(mm)
         Partsource(1:4,kk,ll,mm,iElem) = NodeSource(1:4,GEO%ElemToNodeID(1,iElem)) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(2,iElem)) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(3,iElem)) * (alpha1) * (alpha2) * (1-alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(4,iElem)) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(5,iElem)) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(6,iElem)) * (alpha1) * (1-alpha2) * (alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(7,iElem)) * (alpha1) * (alpha2) * (alpha3) + &
              NodeSource(1:4,GEO%ElemToNodeID(8,iElem)) * (1-alpha1) * (alpha2) * (alpha3)
         END DO !mm
       END DO !ll
     END DO !kk
   END DO !iEle
   DEALLOCATE(NodeSource)
CASE('epanechnikov')
  ALLOCATE(tempsource(0:PP_N,0:PP_N,0:PP_N))
  IF(DoInnerParts)  tempcharge= 0.0
  DO iPart = firstPart, lastPart
  IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
!      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
    IF (usevMPF) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
    ELSE
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF ! usevMPF
    iElem = PEM%Element(iPart)
    alpha = 0.0
    tempcharge(iElem) = tempcharge(iElem) + Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
          radius2 = (PartState(iPart,1) - Elem_xGP(1,kk,ll,mm,iElem)) * (PartState(iPart,1) - Elem_xGP(1,kk,ll,mm,iElem)) &
                  + (PartState(iPart,2) - Elem_xGP(2,kk,ll,mm,iElem)) * (PartState(iPart,2) - Elem_xGP(2,kk,ll,mm,iElem)) &
                  + (PartState(iPart,3) - Elem_xGP(3,kk,ll,mm,iElem)) * (PartState(iPart,3) - Elem_xGP(3,kk,ll,mm,iElem))
         IF (radius2 .LT. r2_sf) THEN
           tempsource(kk,ll,mm) = r2_sf - radius2
           alpha = alpha + tempsource(kk,ll,mm)
         ELSE
           tempsource(kk,ll,mm) = 0.0
         END IF
       END DO !mm
     END DO !ll
    END DO !kk
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         PartSource(1:3,kk,ll,mm,iElem) = PartSource(1:3,kk,ll,mm,iElem)  + 1./alpha*tempsource(kk,ll,mm)*PartState(iPart,4:6) &
                                        * Species(PartSpecies(iPart))%ChargeIC &
                                        * Species(PartSpecies(iPart))%MacroParticleFactor
         PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem)  + 1./alpha*tempsource(kk,ll,mm) &
                                        * Species(PartSpecies(iPart))%ChargeIC &
                                        * Species(PartSpecies(iPart))%MacroParticleFactor
       END DO !mm
     END DO !ll
    END DO !kk
#if USE_LOADBALANCE
    CALL LBElemPauseTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO

  IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    ALLOCATE(tempgridsource(1:nElems))
    tempgridsource= 0.0
    ! seemps to be finalize
    DO iElem=1,nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
            alpha = wGP(kk)*wGP(ll)*wGP(mm)/sJ(kk,ll,mm,iElem)
            PartSource(1:4,kk,ll,mm,iElem) = 1./SQRT(alpha) * PartSource(1:4,kk,ll,mm,iElem)*(PP_N+1)**3/ (GEO%Volume(iElem))
            tempgridsource(iElem) = tempgridsource(iElem) + PartSource(4,kk,ll,mm,iElem)*alpha
         END DO !mm
       END DO !ll
      END DO !kk
      ! possible ABS???
      !IF (tempgridsource(iElem).GT.0.0) THEN
      IF (ABS(tempgridsource(iElem)).GT.0.0) THEN
        alpha = tempcharge(iElem)/tempgridsource(iElem)
        PartSource(1:4,:,:,:,iElem) = PartSource(1:4,:,:,:,iElem)*alpha
      END IF
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO
    DEALLOCATE(tempgridsource)
  END IF
  DEALLOCATE(tempsource)

CASE('shape_function','shape_function_simple')
  !-- "normal" particles
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  IF (usevMPF) THEN
    DO iPart=firstPart,LastPart
      IF(doPartInExists)THEN
        IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
      ELSE
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
      END IF
      CALL calcSfSource(4,Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)*w_sf &
        ,Vec1,Vec2,Vec3,PartState(iPart,1:3),iPart,PartVelo=PartState(iPart,4:6))
    END DO ! iPart
  ELSE
    DO iPart=firstPart,LastPart
      IF(doPartInExists)THEN
        IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
      ELSE
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
      END IF
      CALL calcSfSource(4,Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor*w_sf &
        ,Vec1,Vec2,Vec3,PartState(iPart,1:3),iPart,PartVelo=PartState(iPart,4:6))
    END DO ! iPart
  END IF ! usevMPF
  IF(.NOT.DoInnerParts)THEN

    !-- layer particles (only once, i.e., during call with .NOT.DoInnerParts)
    DO iLayer=1,NbrOfSFdepoLayers
      IF (SFdepoLayersAlreadyDone) EXIT
      CALL RANDOM_NUMBER(RandVal)
      IF (SFdepoLayersPartNum(iLayer).GT.0.) THEN
        layerParts=INT(SFdepoLayersPartNum(iLayer)+RandVal)
      ELSE
        layerParts=0
      END IF
      DO iPart=1,layerParts
        SELECT CASE (TRIM(SFdepoLayersSpace(iLayer)))
        CASE('cuboid')
          CALL RANDOM_NUMBER(RandVal2)
          layerPartPos = SFdepoLayersGeo(iLayer,1,:) + &
            (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
        CASE('cylinder')
          PartRadius = SFdepoLayersRadius(iLayer) + 1
          DO WHILE(PartRadius.GT.SFdepoLayersRadius(iLayer))
            CALL RANDOM_NUMBER(RandVal2)
            RandVal2 = RandVal2 * 2. - 1.
            layerPartPos = SFdepoLayersGeo(iLayer,1,:) + SFdepoLayersRadius(iLayer) * &
              (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
            PartRadius = SQRT( (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) * &
              (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) + &
              (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) * &
              (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) + &
              (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) * &
              (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) )
          END DO
        CASE DEFAULT
          CALL abort(__STAMP__, &
            ' Wrong Space for SFdepoLayer: only cuboid and cylinder implemented!')
        END SELECT
        CALL RANDOM_NUMBER(RandVal)
        layerPartPos = layerPartPos + RandVal*SFdepoLayersGeo(iLayer,2,:)
        IF ( SFdepoLayersUseFixBounds(iLayer) ) THEN
          DoCycle=.FALSE.
          DO iSFfix=1,NbrOfSFdepoFixes
            SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
              + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
              + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
            IF (SFfixDistance .GT. 0.) THEN !outside of plane
              DoCycle=.TRUE.
              EXIT
            END IF
          END DO
          IF (DoCycle) CYCLE
        ELSE IF ( SFdepoLayersBounds(iLayer,1,1).GT.layerPartPos(1) .OR. layerPartPos(1).GT.SFdepoLayersBounds(iLayer,2,1) .OR. &
          SFdepoLayersBounds(iLayer,1,2).GT.layerPartPos(2) .OR. layerPartPos(2).GT.SFdepoLayersBounds(iLayer,2,2) .OR. &
          SFdepoLayersBounds(iLayer,1,3).GT.layerPartPos(3) .OR. layerPartPos(3).GT.SFdepoLayersBounds(iLayer,2,3) ) THEN
          CYCLE !outside of bounds
        END IF
        CALL calcSfSource(1 &
          ,Species(SFdepoLayersSpec(iLayer))%ChargeIC*SFdepoLayersMPF(iLayer)*w_sf &
          ,Vec1,Vec2,Vec3,layerPartPos,iPart,const_opt=ConstantSFdepoLayers)
      END DO ! iPart
      IF (iLayer.EQ.NbrOfSFdepoLayers .AND. ConstantSFdepoLayers) THEN
        SFdepoLayersAlreadyDone=.TRUE.
      END IF
    END DO ! iLayer=1,NbrOfSFdepoLayers

    !--SFResampleAnalyzeSurfCollis
    IF (SFResampleAnalyzeSurfCollis) THEN
      iPart=0
      DO iPart2=1,LastAnalyzeSurfCollis%PartNumberDepo
        !get random (equal!) position between [1,PartNumberSamp]
        CALL RANDOM_NUMBER(RandVal)
        iPart=MIN(1+INT(RandVal*REAL(LastAnalyzeSurfCollis%PartNumberSamp)),LastAnalyzeSurfCollis%PartNumberSamp)
        !perform surfaceflux-like push into sf-layer outside of mesh
        CALL RANDOM_NUMBER(RandVal)
        FractPush = RandVal*LastAnalyzeSurfCollis%pushTimeStep*LastAnalyzeSurfCollis%WallState(4:6,iPart)
        IF ( DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,FractPush).LE.r_SF  ) THEN
          layerPartPos = LastAnalyzeSurfCollis%WallState(1:3,iPart) + FractPush
          IF ( LastAnalyzeSurfCollis%UseFixBounds ) THEN
            DoCycle=.FALSE.
            DO iSFfix=1,NbrOfSFdepoFixes
              SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
                + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
                + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
              IF (SFfixDistance .GT. 0.) THEN !outside of plane
                DoCycle=.TRUE.
                EXIT
              END IF
            END DO
            IF (DoCycle) CYCLE
          ELSE IF ( LastAnalyzeSurfCollis%Bounds(1,1).GT.layerPartPos(1) .OR. &
            layerPartPos(1).GT.LastAnalyzeSurfCollis%Bounds(2,1) .OR. &
            LastAnalyzeSurfCollis%Bounds(1,2).GT.layerPartPos(2) .OR. &
            layerPartPos(2).GT.LastAnalyzeSurfCollis%Bounds(2,2) .OR. &
            LastAnalyzeSurfCollis%Bounds(1,3).GT.layerPartPos(3) .OR. &
            layerPartPos(3).GT.LastAnalyzeSurfCollis%Bounds(2,3) ) THEN
            CYCLE !outside of bounds
          END IF
        ELSE
          CYCLE !outside of r_SF
        END IF
        CALL calcSfSource(4 &
          ,Species(LastAnalyzeSurfCollis%Species(iPart))%ChargeIC &
          *Species(LastAnalyzeSurfCollis%Species(iPart))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,layerPartPos,iPart2,PartVelo=LastAnalyzeSurfCollis%WallState(4:6,iPart))
      END DO ! iPart2
    END IF !SFResampleAnalyzeSurfCollis

    !-- external particles
#ifdef MPI
    IF (usevMPF) THEN
      DO iPart=1,NbrOfextParticles  !external Particles
        CALL calcSfSource(4,Species(ExtPartSpecies(iPart))%ChargeIC*ExtPartMPF(iPart)*w_sf &
          ,Vec1,Vec2,Vec3,ExtPartState(iPart,1:3),iPart,PartVelo=ExtPartState(iPart,4:6))
      END DO
    ELSE
      DO iPart=1,NbrOfextParticles  !external Particles
        CALL calcSfSource(4 &
          ,Species(ExtPartSpecies(iPart))%ChargeIC*Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf &
          ,Vec1,Vec2,Vec3,ExtPartState(iPart,1:3),iPart,PartVelo=ExtPartState(iPart,4:6))
      END DO
    END IF ! usevMPF
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
#endif /*MPI*/

    !-- add const. PartSource (only once, i.e., during call with .NOT.DoInnerParts)
    IF (PartSourceConstExists) THEN
      DO iElem = 1, nElems
        DO kk = 0, PP_N
          DO ll = 0, PP_N
            DO mm = 0, PP_N
              PartSource(1:4,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem) + PartSourceConst(1:4,mm,ll,kk,iElem)
            END DO !mm
          END DO !ll
        END DO !kk
      END DO !iElem
    END IF !PartSourceConstExists
  END IF !.NOT.DoInnerParts

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map PartSource from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('shape_function_1d')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    IF (usevMPF) THEN
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
    ELSE
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
    END IF ! usevMPF
    !IF(fac(4).GT.0.) print*,'charge pos'
    Fac(1:3) = PartState(iPart,4:6)*Fac(4)
    !-- determine which background mesh cells (and interpolation points within) need to be considered
    chargedone(:) = .FALSE.
    DO iCase = 1, NbrOfCases
      DO ind = 1,3
        ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
             casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
      END DO
      IF(sf1d_dir.EQ.1)THEN
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = GEO%FIBGMjmax
        lmin = GEO%FIBGMjmin
        mmax = GEO%FIBGMkmax
        mmin = GEO%FIBGMkmin
      ELSEIF(sf1d_dir.EQ.2)THEN
        kmax = GEO%FIBGMimax
        kmin = GEO%FIBGMimin
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = GEO%FIBGMkmax
        mmin = GEO%FIBGMkmin
      ELSE
        kmax = GEO%FIBGMimax
        kmin = GEO%FIBGMimin
        lmax = GEO%FIBGMjmax
        lmin = GEO%FIBGMjmin
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
      END IF
      !-- go through all these cells
      DO kk = kmin,kmax
        DO ll = lmin, lmax
          DO mm = mmin, mmax
            !--- go through all mapped elements not done yet
            DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
              IF(ElemID.GT.nElems) CYCLE
              IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                !--- go through all gauss points
                !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  !-- calculate distance between gauss and particle
                  radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                          * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                  !-- calculate charge and current density at ip point using a shape function
                  !-- currently only one shapefunction available, more to follow (including structure change)
                  IF (radius2 .LT. r2_sf) THEN
                    S = 1. - r2_sf_inv * radius2
                  !radius2=GaussDistance(k,l,m)
                  !IF (radius2 .LT. 1.0) THEN
                  !  S = 1 -  radius2
                    S1 = S*S
                    DO expo = 3, alpha_sf
                      S1 = S*S1
                    END DO
                    PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                    PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                  END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END DO ! iCase (periodicity)
  END DO ! i
#ifdef MPI
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF

    DO iPart=1,NbrOfextParticles  !external Particles
      IF (usevMPF) THEN
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
      ELSE
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
      END IF ! usevMPF
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      chargedone(:) = .FALSE.
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        IF(sf1d_dir.EQ.1)THEN
          kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
          kmax = MIN(kmax,GEO%FIBGMimax)
          kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
          kmin = MAX(kmin,GEO%FIBGMimin)
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSEIF(sf1d_dir.EQ.2)THEN
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
          lmax = MIN(lmax,GEO%FIBGMjmax)
          lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
          lmin = MAX(lmin,GEO%FIBGMjmin)
          mmax = GEO%FIBGMkmax
          mmin = GEO%FIBGMkmin
        ELSE
          kmax = GEO%FIBGMimax
          kmin = GEO%FIBGMimin
          lmax = GEO%FIBGMjmax
          lmin = GEO%FIBGMjmin
          mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
          mmax = MIN(mmax,GEO%FIBGMkmax)
          mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
          mmin = MAX(mmin,GEO%FIBGMkmin)
        END IF
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain)
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                              * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. r2_sf) THEN
                        S = 1. - r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/

  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map source from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('shape_function_cylindrical','shape_function_spherical')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (GEO%nPeriodicVectors.EQ.1) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  END IF
  IF (GEO%nPeriodicVectors.EQ.2) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  END IF
  IF (GEO%nPeriodicVectors.EQ.3) THEN
    Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
  END IF
  DO iPart=firstPart,LastPart
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! compute local radius
    local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:SfRadiusInt),PartState(iPart,1:SfRadiusInt)))
    local_r2_sf=local_r_sf*local_r_sf
    local_r2_sf_inv=1./local_r2_sf
    IF (usevMPF) THEN
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
    ELSE
      Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
    END IF ! usevMPF
    !IF(fac(4).GT.0.) print*,'charge pos'
    Fac(1:3) = PartState(iPart,4:6)*Fac(4)
    !-- determine which background mesh cells (and interpolation points within) need to be considered
    DO iCase = 1, NbrOfCases
      chargedone(:) = .FALSE.
      DO ind = 1,3
        ShiftedPart(ind) = PartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
             casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
      END DO
      kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      !-- go through all these cells
      DO kk = kmin,kmax
        DO ll = lmin, lmax
          DO mm = mmin, mmax
            !--- go through all mapped elements not done yet
            DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
              ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
              IF(ElemID.GT.nElems) CYCLE
              IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                !--- go through all gauss points
                !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  !-- calculate distance between gauss and particle
                  dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                  IF(dX.GT.local_r_sf) CYCLE
                  dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                  IF(dY.GT.local_r_sf) CYCLE
                  dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                  IF(dZ.GT.local_r_sf) CYCLE
                  radius2 = dX*dX+dY*dY+dZ*dZ
                  !-- calculate charge and current density at ip point using a shape function
                  !-- currently only one shapefunction available, more to follow (including structure change)
                  IF (radius2 .LT. local_r2_sf) THEN
                    S = 1. - local_r2_sf_inv * radius2
                  !radius2=GaussDistance(k,l,m)
                  !IF (radius2 .LT. 1.0) THEN
                  !  S = 1 -  radius2
                    S1 = S*S
                    DO expo = 3, alpha_sf
                      S1 = S*S1
                    END DO
                    PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                    PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                  END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END DO ! iCase (periodicity)
  END DO ! i
#ifdef MPI
  IF(.NOT.DoInnerParts)THEN
    Vec1(1:3) = 0.
    Vec2(1:3) = 0.
    Vec3(1:3) = 0.
    IF (NbrOfextParticles .GT. 0) THEN
      IF (GEO%nPeriodicVectors.EQ.1) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      END IF
      IF (GEO%nPeriodicVectors.EQ.2) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      END IF
      IF (GEO%nPeriodicVectors.EQ.3) THEN
        Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
        Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
        Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
      END IF
    END IF

    DO iPart=1,NbrOfextParticles  !external Particles
      local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(iPart,1:SfRadiusInt),PartState(iPart,1:SfRadiusInt)))
      local_r2_sf=local_r_sf*local_r_sf
      local_r2_sf_inv=1./local_r2_sf
      IF (usevMPF) THEN
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
      ELSE
        Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC &
              * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
      END IF ! usevMPF
      Fac(1:3) = ExtPartState(iPart,4:6)*Fac(4)
      !-- determine which background mesh cells (and interpolation points within) need to be considered
      DO iCase = 1, NbrOfCases
        chargedone(:) = .FALSE.
        DO ind = 1,3
          ShiftedPart(ind) = ExtPartState(iPart,ind) + casematrix(iCase,1)*Vec1(ind) + &
               casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
        END DO
        kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
        DO kk = kmin,kmax
          DO ll = lmin, lmax
            DO mm = mmin, mmax
              !--- go through all mapped elements not done yet
              DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
                ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
                IF(ElemID.GT.nElems) CYCLE
                IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
                  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
                  !--- go through all gauss points
                  !CALL ComputeGaussDistance(PP_N,r2_sf_inv,ShiftedPart,ElemDepo_xGP(:,:,:,:,ElemID),GaussDistance)
                  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                    !-- calculate distance between gauss and particle
                      dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                      IF(dX.GT.local_r_sf) CYCLE
                      dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                      IF(dY.GT.local_r_sf) CYCLE
                      dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                      IF(dZ.GT.local_r_sf) CYCLE
                      radius2 = dX*dX+dY*dY+dZ*dZ
                      !-- calculate charge and current density at ip point using a shape function
                      !-- currently only one shapefunction available, more to follow (including structure change)
                      IF (radius2 .LT. local_r2_sf) THEN
                        S = 1. - local_r2_sf_inv * radius2
                        !IF(S.LT.0.) print*,'dist neg '
                      !radius2=GaussDistance(k,l,m)
                      !IF (radius2 .LT. 1.0) THEN
                      !  S = 1 -  radius2
                        S1 = S*S
                        DO expo = 3, alpha_sf
                          S1 = S*S1
                        END DO
                        PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                        PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                      END IF
                  END DO; END DO; END DO
                  chargedone(ElemID) = .TRUE.
                END IF
              END DO ! ppp
            END DO ! mm
          END DO ! ll
        END DO ! kk
      END DO
    END DO
    ! deallocate external state
    SDEALLOCATE(ExtPartState)
    SDEALLOCATE(ExtPartSpecies)
    SDEALLOCATE(ExtPartToFIBGM)
    SDEALLOCATE(ExtPartMPF)
    NbrOfExtParticles=0
  END IF
#endif /*MPI*/
  IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
    ! map PartSource from Equististant points on Gauss-Points
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
    END DO ! iElem=1,PP_nElems
  END IF

CASE('delta_distri')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    DO iPart=firstPart,LastPart
      IF(doPartInExists)THEN
        IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
      ELSE
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
      END IF
      IF(PEM%Element(iPart).EQ.iElem)THEN
        IF (usevMPF) THEN
          prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
        ELSE
          prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
        END IF ! usevMPF
        ! Map Particle to -1|1 space (re-used in interpolation)
        IF(.NOT.DoRefMapping)THEN
          CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem)
        END IF
        ! get value of test function at particle position
        SELECT CASE(DeltaType)
        CASE(1)
          ! xi   -direction
          CALL LagrangeInterpolationPolys(PartPosRef(1,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(1,:))
          ! eta  -direction
          CALL LagrangeInterpolationPolys(PartPosRef(2,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(2,:))
          ! zeta -direction
          CALL LagrangeInterpolationPolys(PartPosRef(3,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(3,:))
        CASE(2)
          DO i=0,NDepo
            ! xi   -direction
             CALL BernsteinPolynomial(NDepo,i,PartPosRef(1,iPart),L_xi(1,i),NDepoChooseK)
            ! eta  -direction
             CALL BernsteinPolynomial(NDepo,i,PartPosRef(2,iPart),L_xi(2,i),NDepoChooseK)
            ! zeta  -direction
             CALL BernsteinPolynomial(NDepo,i,PartPosRef(3,iPart),L_xi(3,i),NDepoChooseK)
          END DO ! i
        CASE(3)
          ! xi - direction
          CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(1,iPart),L_xi(1,:))
          ! eta - direction
          CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(2,iPart),L_xi(2,:))
          ! zeta - direction
          CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(3,iPart),L_xi(3,:))
        END SELECT
        DO k=0,NDepo
          DO j=0,NDepo
            DO i=0,NDepo
         !     print*,'i,j,k,L',i,j,k,L_xi(1,i)* L_xi(2,j)* L_xi(3,k)
              DeltaIntCoeff = L_xi(1,i)* L_xi(2,j)* L_xi(3,k)*prefac
              PartSource(1:3,i,j,k,iElem) = PartSource(1:3,i,j,k,iElem) + DeltaIntCoeff*PartState(iPart,4:6)
              PartSource( 4 ,i,j,k,iElem) = PartSource( 4 ,i,j,k,iElem) + DeltaIntCoeff
            END DO ! i
          END DO ! j
        END DO ! k
      END IF ! Particle in Element
    END DO ! ParticleVecLength
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem
  IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    DO iElem=1,PP_nElems
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
            PartSource( : ,i,j,k,iElem) = PartSource( : ,i,j,k,iElem) * DDMassInv(i,j,k,iElem)
          END DO ! i
        END DO ! j
      END DO ! k
      IF(DeltaDistriChangeBasis)THEN
        CALL ChangeBasis3D(4,NDepo,PP_N,Vdm_NDepo_GaussN,PartSource(1:4,0:NDepo,0:NDepo,0:NDepo,iElem)&
                                                        ,PartSource(1:4,0:PP_N ,0:PP_N ,0:PP_N, iElem))
      END IF
#if USE_LOADBALANCE
      CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! loop over all elems
  END IF ! DoInnerParts
CASE('nearest_gausspoint')
  IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
  SAVE_GAUSS = .FALSE.
  IF(TRIM(InterpolationType).EQ.'nearest_gausspoint') SAVE_GAUSS = .TRUE.
  IF(MOD(PP_N,2).EQ.0) THEN
    a = PP_N/2
    b = a
  ELSE
    a = (PP_N+1)/2
    b = a-1
  END IF
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    DO iPart=firstPart,LastPart
      IF(doPartInExists)THEN
        IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
      ELSE
        IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
      END IF
      IF(PEM%Element(iPart).EQ.iElem)THEN
        IF (usevMPF) THEN
          prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
        ELSE
          prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
        END IF ! usevMPF
        ! Map Particle to -1|1 space (re-used in interpolation)
        !IF(.NOT.DoRefMapping)THEN
        !  CALL GetPositionInRefElem(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
        !END IF
        ! Find out which gausspoint is closest and add up charges and currents
        !! x-direction
        IF(.NOT.SAVE_GAUSS) THEN
          k = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(1,iPart)).GE.GaussBorder(PP_N-ii))THEN
              k = PP_N-ii
              EXIT
            END IF
          END DO
          k = NINT((PP_N+SIGN(2.0*k-PP_N,PartPosRef(1,iPart)))/2)
          !! y-direction
          l = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(2,iPart)).GE.GaussBorder(PP_N-ii))THEN
              l = PP_N-ii
              EXIT
            END IF
          END DO
          l = NINT((PP_N+SIGN(2.0*l-PP_N,PartPosRef(2,iPart)))/2)
          !! z-direction
          m = a
          DO ii = 0,b-1
            IF(ABS(PartPosRef(3,iPart)).GE.GaussBorder(PP_N-ii))THEN
              m = PP_N-ii
              EXIT
            END IF
          END DO
          m = NINT((PP_N+SIGN(2.0*m-PP_N,PartPosRef(3,iPart)))/2)
        END IF
!#if (PP_nVar==8)
        PartSource(1:3,k,l,m,iElem) = PartSource(1:3,k,l,m,iElem) + PartState(iPart,4:6) * prefac
!#endif
        PartSource( 4 ,k,l,m,iElem) = PartSource( 4 ,k,l,m,iElem) + prefac
        !IF (SAVE_GAUSS) THEN
        !  PartPosGauss(iPart,1) = k
        !  PartPosGauss(iPart,2) = l
        !  PartPosGauss(iPart,3) = m
        !END IF
      END IF ! Element .EQ. iElem
    END DO ! iPart
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem=1,PP_nElems
  IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
    DO iElem=1,PP_nElems
      DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
      ! get densities by dividing by gauss volume
!#if (PP_nVar==8)
        PartSource(1:4,k,l,m,iElem) = PartSource(1:4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#else
!        PartSource(4,k,l,m,iElem) = PartSource(4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#endif
      END DO; END DO; END DO
#if USE_LOADBALANCE
      CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! iElem=1,PP_nElems
  END IF
CASE('cartmesh_volumeweighting')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  ! IF(DoInnerParts) BGMSource=0.0 ! not possible due to periodic stuff --> two communications
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    IF (usevMPF) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
    ELSE
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF ! usevMPF
    !Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
    k = FLOOR(PartState(iPart,1)/BGMdeltas(1))
    l = FLOOR(PartState(iPart,2)/BGMdeltas(2))
    m = FLOOR(PartState(iPart,3)/BGMdeltas(3))
    alpha1 = (PartState(iPart,1) / BGMdeltas(1)) - k
    alpha2 = (PartState(iPart,2) / BGMdeltas(2)) - l
    alpha3 = (PartState(iPart,3) / BGMdeltas(3)) - m
    TSource(:) = 0.0
!#if (PP_nVar==8)
    TSource(1) = PartState(iPart,4)*Charge
    TSource(2) = PartState(iPart,5)*Charge
    TSource(3) = PartState(iPart,6)*Charge
!#endif
    TSource(4) = Charge

    BGMSource(k,l,m,1:4)       = BGMSource(k,l,m,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(1-alpha3))
    BGMSource(k,l,m+1,1:4)     = BGMSource(k,l,m+1,1:4) + (TSource * (1-alpha1)*(1-alpha2)*(alpha3))
    BGMSource(k,l+1,m,1:4)     = BGMSource(k,l+1,m,1:4) + (TSource * (1-alpha1)*(alpha2)*(1-alpha3))
    BGMSource(k,l+1,m+1,1:4)   = BGMSource(k,l+1,m+1,1:4) + (TSource * (1-alpha1)*(alpha2)*(alpha3))
    BGMSource(k+1,l,m,1:4)     = BGMSource(k+1,l,m,1:4) + (TSource * (alpha1)*(1-alpha2)*(1-alpha3))
    BGMSource(k+1,l,m+1,1:4)   = BGMSource(k+1,l,m+1,1:4) + (TSource * (alpha1)*(1-alpha2)*(alpha3))
    BGMSource(k+1,l+1,m,1:4)   = BGMSource(k+1,l+1,m,1:4) + (TSource * (alpha1)*(alpha2)*(1-alpha3))
    BGMSource(k+1,l+1,m+1,1:4) = BGMSource(k+1,l+1,m+1,1:4) + (TSource * (alpha1)*(alpha2)*(alpha3))
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_LOADBALANCE
  CALL LBPauseTime(LB_CARTMESHDEPO,tLBStart)
#endif /*USE_LOADBALANCE*/
#ifdef MPI
  ! should be treated in this way, unforunately, we would negelct the periodic stuff
  !IF(.NOT.DoInnerParts)
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif

  ! Step 2: Interpolation of densities onto grid
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
         k = GaussBGMIndex(1,kk,ll,mm,iElem)
         l = GaussBGMIndex(2,kk,ll,mm,iElem)
         m = GaussBGMIndex(3,kk,ll,mm,iElem)
         alpha1 = GaussBGMFactor(1,kk,ll,mm,iElem)
         alpha2 = GaussBGMFactor(2,kk,ll,mm,iElem)
         alpha3 = GaussBGMFactor(3,kk,ll,mm,iElem)
!#if (PP_nVar==8)
         DO i = 1,3
           PartSource(i,kk,ll,mm,iElem) = PartSource(i,kk,ll,mm,iElem)            + &
                BGMSource(k,l,m,i) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k,l,m+1,i) * (1-alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k,l+1,m,i) * (1-alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k,l+1,m+1,i) * (1-alpha1) * (alpha2) * (alpha3) + &
                BGMSource(k+1,l,m,i) * (alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k+1,l,m+1,i) * (alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k+1,l+1,m,i) * (alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k+1,l+1,m+1,i) * (alpha1) * (alpha2) * (alpha3)
         END DO
!#endif
           PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem)          + &
              BGMSource(k,l,m,4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k,l,m+1,4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k,l+1,m,4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k,l+1,m+1,4) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSource(k+1,l,m,4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k+1,l,m+1,4) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k+1,l+1,m,4) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k+1,l+1,m+1,4) * (alpha1) * (alpha2) * (alpha3)
       END DO !mm
     END DO !ll
   END DO !kk
#if USE_LOADBALANCE
   CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
 END DO !iElem
 !DEALLOCATE(BGMSource)
CASE('cartmesh_splines')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  !ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
  ! IF(DoInnerParts) BGMSource=0. not possible due to periodic stuff
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  BGMSource(:,:,:,:) = 0.0
  DO iPart = firstPart, lastPart
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
!      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
    IF (usevMPF) THEN
      Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
    ELSE
      Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
    END IF ! usevMPF
    PosInd(1) = FLOOR(PartState(iPart,1)/BGMdeltas(1))
    PosInd(2) = FLOOR(PartState(iPart,2)/BGMdeltas(2))
    PosInd(3) = FLOOR(PartState(iPart,3)/BGMdeltas(3))
    !print*,'posind(1:3),charge',posInd,charge
    DO dir = 1,3               ! x,y,z direction
      DO weightrun = 0,3
        DO mm = 0, 3
          IF (mm.EQ.weightrun) then
            auxiliary(mm) = 1.0
          ELSE
            auxiliary(mm) = 0.0
          END IF
        END DO
        CALL DeBoor(PosInd(dir),auxiliary,PartState(iPart,dir),weight(dir,weightrun),dir)
      END DO
    END DO
    DO k = PosInd(1)-1, PosInd(1)+2
      kk = abs(k - PosInd(1) - 2)
      DO l = PosInd(2)-1, PosInd(2)+2
        ll = abs(l - PosInd(2) - 2)
        DO m = PosInd(3)-1, PosInd(3)+2
          mm = abs(m - PosInd(3) - 2)
          locweight = weight(1,kk)*weight(2,ll)*weight(3,mm)*charge
!#if (PP_nVar==8)
          BGMSource(k,l,m,1) = BGMSource(k,l,m,1) + PartState(iPart,4)* locweight
          BGMSource(k,l,m,2) = BGMSource(k,l,m,2) + PartState(iPart,5)* locweight
          BGMSource(k,l,m,3) = BGMSource(k,l,m,3) + PartState(iPart,6)* locweight
!#endif
          BGMSource(k,l,m,4) = BGMSource(k,l,m,4) + locweight
       !   print*,'BMGSOURCE4',BGMSOURCE(k,l,m,4)
        END DO
      END DO
    END DO
  END DO
  BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_LOADBALANCE
  CALL LBPauseTime(LB_CARTMESHDEPO,tLBStart)
#endif /*USE_LOADBALANCE*/
#ifdef MPI
  !IF(.NOT.DoInnerParts)THEN has to be communicated each time :(
  CALL MPISourceExchangeBGM()
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif

  ! Step 2: Interpolation of densities onto grid
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem = 1, nElems
    DO kk = 0, PP_N
      DO ll = 0, PP_N
        DO mm = 0, PP_N
          k = GaussBGMIndex(1,kk,ll,mm,iElem)
          l = GaussBGMIndex(2,kk,ll,mm,iElem)
          m = GaussBGMIndex(3,kk,ll,mm,iElem)
          DO r = k-1,k+2
            u = r-k+2
            DO ss = l-1,l+2
              v = ss-l+2
              DO t = m-1,m+2
                w = t-m+2
!#if (PP_nVar==8)
                PartSource(1:4,kk,ll,mm,iElem) = PartSource(1:4,kk,ll,mm,iElem) + BGMSource(r,ss,t,1:4) * GPWeight(iElem,kk,ll,mm,u,v,w)
                !DO i = 1,3
                !  PartSource(i,kk,ll,mm,iElem) = PartSource(i,kk,ll,mm,iElem) + BGMSource(r,ss,t,i) * GPWeight(iElem,kk,ll,mm,u,v,w)
                !END DO
!#endif
                !PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem) + BGMSource(r,ss,t,4) * GPWeight(iElem,kk,ll,mm,u,v,w)
              END DO !t
            END DO !s
          END DO !r
        END DO !mm
      END DO !ll
    END DO !kk
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO !iElem
 !DEALLOCATE(BGMSource)
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Unknown DepositionType in pic_depo.f90')
END SELECT

RETURN
END SUBROUTINE Deposition


#ifndef MPI
SUBROUTINE PeriodicSourceExchange()
!============================================================================================================================
! Exchange sources in periodic case
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
USE MOD_Particle_Mesh_Vars,  ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)         :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i,k,l,m,k2,l2,m2
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,GEO%nPeriodicVectors
  DO k = BGMminX, BGMmaxX
    k2 = k + GEO%PeriodicBGMVectors(1,i)
    DO l = BGMminY, BGMmaxY
      l2 = l + GEO%PeriodicBGMVectors(2,i)
      DO m = BGMminZ, BGMmaxZ
        m2 = m + GEO%PeriodicBGMVectors(3,i)
        IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
          IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
            IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
              BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
              BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO

END SUBROUTINE PeriodicSourceExchange
#else /*MPI*/
SUBROUTINE MPISourceExchangeBGM()
!=================================================================================================================================
! Exchange sources in periodic case for MPI
!==================================================================================================================================
! use MODULES
USE MOD_Particle_MPI_Vars,  ONLY: PartMPI,tMPIMEssage
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_PICDepo_Vars
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT)        :: BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)
TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)
INTEGER                     :: send_request(0:PartMPI%nProcs-1)
INTEGER                     :: recv_request(0:PartMPI%nProcs-1)
INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)
INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)
INTEGER                     :: iProc,k,l,m,n, ppp, Counter, MsgLength(0:PartMPI%nProcs-1), iPer
INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)
INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)
INTEGER                     :: allocStat, Counter2
INTEGER                     :: messageCounterS, messageCounterR
INTEGER                     :: myRealKind, k2,l2,m2
REAL                        :: myRealTestValue
!-----------------------------------------------------------------------------------------------------------------------------------

myRealKind = KIND(myRealTestValue)
IF (myRealKind.EQ.4) THEN
 myRealKind = MPI_REAL
ELSE IF (myRealKind.EQ.8) THEN
 myRealKind = MPI_DOUBLE_PRECISION
ELSE
 myRealKind = MPI_REAL
END IF

!--- Determine which Sources actually need to be sent (<> 0) and build corresponding true/false list
!    One list per process
DO iProc = 0,PartMPI%nProcs-1
   MsgLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      MsgLength(iProc) = (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2) + 1) * &
                         (PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3) - PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3) + 1)
   END IF
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO k = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         MsgLength(iProc) = MsgLength(iProc) + &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,1) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,2) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
              (PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(2,3) -&
               PartMPI%DepoBGMConnect(iProc)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
      END DO
   END IF
   IF((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor))THEN
      ALLOCATE(send_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
      ALLOCATE(recv_message(iProc)%content_log(1:MsgLength(iProc)), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
      END IF
   END IF
   !--- check which sources are <> 0
   Counter = 0
   SourceLength(iProc) = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
               send_message(iProc)%content_log(Counter) = .TRUE.
               SourceLength(iProc) = SourceLength(iProc) + 1
            ELSE
               send_message(iProc)%content_log(Counter) = .FALSE.
            END IF
          END DO
        END DO
      END DO
   END IF
   !--- same for periodic
   IF(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
          DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                 PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
            DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                   PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter = Counter + 1
               IF (ANY(BGMSource(k,l,m,:).NE.0.0)) THEN
                  send_message(iProc)%content_log(Counter) = .TRUE.
                  SourceLength(iProc) = SourceLength(iProc) + 1
               ELSE
                  send_message(iProc)%content_log(Counter) = .FALSE.
               END IF
            END DO
          END DO
         END DO
      END DO
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%Content_log,MsgLength(iProc),MPI_LOGICAL,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF

!--- Assemble actual sources to send
 DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ALLOCATE(send_message(iProc)%content(1:SourceLength(iProc)*4), STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
      END IF
   END IF
   Counter = 0
   Counter2 = 0
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
      DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
        DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
          DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
             Counter = Counter + 1
             IF (send_message(iProc)%content_log(Counter)) THEN
                Counter2 = Counter2 + 1
                DO n = 1,4
                   send_message(iProc)%content((Counter2-1)*4 +n) = BGMSource(k,l,m,n)
                END DO
             END IF
          END DO
        END DO
      END DO
   END IF
   IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
      DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
         DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
           DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                Counter = Counter + 1
                IF (send_message(iProc)%content_log(Counter)) THEN
                   Counter2 = Counter2 + 1
                   DO ppp = 1,4
                      send_message(iProc)%content((Counter2-1)*4 +ppp) = BGMSource(k,l,m,ppp)
                   END DO
                END IF
             END DO
           END DO
         END DO
      END DO
   END IF
END DO

!--- allocate actual PartSource receive buffer
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      Counter = 0
      DO k = 1, MsgLength(iProc)
         IF (recv_message(iProc)%Content_log(k)) THEN
            Counter = Counter + 1
         END IF
      END DO
      RecvLength(iProc) = Counter
      IF (RecvLength(iProc).GT.0) THEN
         ALLOCATE(recv_message(iProc)%content(1:Counter*4), STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
         END IF
      END IF
   END IF
END DO
!--- communicate
messageCounterS = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(SourceLength(iProc).GT.0)) THEN
      ! MPI_ISEND true/false list for all border BGM points
      messageCounterS = messageCounterS + 1
      CALL MPI_ISEND(send_message(iProc)%content,SourceLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     send_request(messageCounterS), IERROR)
   END IF
END DO
messageCounterR = 0
DO iProc = 0,PartMPI%nProcs-1
   IF (((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.&
        (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)).AND.(RecvLength(iProc).GT.0)) THEN
      ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
      messageCounterR = messageCounterR + 1
      CALL MPI_IRECV(recv_message(iProc)%content,RecvLength(iProc)*4,myRealKind,iProc,1,PartMPI%COMM, &
                     recv_request(messageCounterR), IERROR)
   END IF
END DO
! MPI_WAITALL for the non-blocking MPI-communication to be finished
IF (messageCounterS .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
END IF
IF (messageCounterR .GE. 1) THEN
   CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
END IF
!--- Deallocate Send Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor)) THEN
      DEALLOCATE(send_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
      END IF
      IF (SourceLength(iProc).GT.0) THEN
         DEALLOCATE(send_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
         END IF
      END IF
   END IF
END DO

!--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
!---           received sources have been added!
IF ((GEO%nPeriodicVectors.GT.0).AND.(GEO%SelfPeriodic)) THEN
   DO iPer = 1, GEO%nPeriodicVectors
      DO k = BGMminX, BGMmaxX
         k2 = k + GEO%PeriodicBGMVectors(1,iPer)
        DO l = BGMminY, BGMmaxY
          l2 = l + GEO%PeriodicBGMVectors(2,iPer)
          DO m = BGMminZ, BGMmaxZ
             m2 = m + GEO%PeriodicBGMVectors(3,iPer)
             IF ((k2.GE.BGMminX).AND.(k2.LE.BGMmaxX)) THEN
             IF ((l2.GE.BGMminY).AND.(l2.LE.BGMmaxY)) THEN
             IF ((m2.GE.BGMminZ).AND.(m2.LE.BGMmaxZ)) THEN
                BGMSource(k,l,m,:) = BGMSource(k,l,m,:) + BGMSource(k2,l2,m2,:)
                BGMSource(k2,l2,m2,:) = BGMSource(k,l,m,:)
             END IF
             END IF
             END IF
          END DO
        END DO
      END DO
   END DO
END IF

!--- Add Sources and Deallocate Receive Message Buffers
DO iProc = 0,PartMPI%nProcs-1
   IF (RecvLength(iProc).GT.0) THEN
      Counter = 0
      Counter2 = 0
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor) THEN
         DO k = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,1), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,1)
         DO l = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,2), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,2)
         DO m = PartMPI%DepoBGMConnect(iProc)%BGMBorder(1,3), PartMPI%DepoBGMConnect(iProc)%BGMBorder(2,3)
            Counter = Counter + 1
            IF(recv_message(iProc)%content_log(Counter))THEN
               Counter2 = Counter2 + 1
               DO n = 1,4
                 BGMSource(k,l,m,n) = BGMSource(k,l,m,n) + recv_message(iProc)%content((Counter2-1)*4+n)
               END DO
            END IF
         END DO
         END DO
         END DO
      END IF
      IF (PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor) THEN
         DO n = 1, PartMPI%DepoBGMConnect(iProc)%BGMPeriodicBorderCount
           DO k = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,1),&
                  PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,2)
               DO m = PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(1,3),&
                      PartMPI%DepoBGMConnect(iProc)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter = Counter + 1
                  IF(recv_message(iProc)%content_log(Counter))THEN
                     Counter2 = Counter2 + 1
                     DO ppp = 1,4
                        BGMSource(k,l,m,ppp) = BGMSource(k,l,m,ppp) + recv_message(iProc)%content((Counter2-1)*4+ppp)
                     END DO
                  END IF
               END DO
             END DO
           END DO
         END DO
      END IF
      IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
         DEALLOCATE(recv_message(iProc)%content, STAT=allocStat)
         IF (allocStat .NE. 0) THEN
            CALL abort(&
            __STAMP__&
            ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
         END IF
      END IF
   END IF
   IF ((PartMPI%DepoBGMConnect(iProc)%isBGMPeriodicNeighbor).OR.(PartMPI%DepoBGMConnect(iProc)%isBGMNeighbor)) THEN
      DEALLOCATE(recv_message(iProc)%content_log, STAT=allocStat)
      IF (allocStat .NE. 0) THEN
         CALL abort(&
         __STAMP__&
         ,'ERROR in MPISourceExchangeBGM: cannot deallocate recv_message')
      END IF
   END IF
END DO

END SUBROUTINE MPISourceExchangeBGM
#endif /*MPI*/

SUBROUTINE DeBoor(PosInd, aux, coord, results, dir)
!============================================================================================================================
! recursive function for evaluating a b-spline basis function
!============================================================================================================================
! use MODULES
   USE MOD_PICDepo_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   INTEGER, INTENT(IN)                          :: PosInd, dir
   REAL, INTENT(IN)                             :: coord
   REAL, INTENT(INOUT)                          :: aux(0:3), results
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,k,jL,jR
   REAL                              :: hlp1,hlp2
!-----------------------------------------------------------------------------------------------------------------------------------
    DO i = 0, 2
       DO k = 0, 2-i
          jL = PosInd - k
          jR = jL + 3 - i

          hlp1 = jR * BGMdeltas(dir) - coord
          hlp2 = coord - jL * BGMdeltas(dir)

          aux(k) = (hlp1 * aux(k+1) + hlp2 * aux(k)) / (hlp1+hlp2)
       ENDDO
    ENDDO
    results = aux(0)

    RETURN
END SUBROUTINE DeBoor


SUBROUTINE DeBoorRef(N_in,NKnots,Knots,Xi_in,UBspline)
!===================================================================================================================================
! DeBoor algorithms for uniform B-Splines
! rule of DeBoor algorithm
! N_i^0 (x) = 1 if x in [u_i, u_i+1); else 0
! N_i^n (x) = (x-u_i) /(u_i+n-u_i) N_i^n-1(x) + (u_i+n+1 - x)/(u_i+n+1 - u_i+1) N_i+1^n-1(x)
! this algorithm evaluates the complete 1D basis fuction, because certain knots can be reused
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: N_in,Nknots
REAL,INTENT(IN)         :: Xi_in
REAL,INTENT(OUT)        :: UBspline(0:N_in)
REAL,INTENT(IN)         :: knots(0:Nknots) ! range of parameter
!REAL,INTENT(IN)         :: DXi is 2 for [-1,1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,n, last!,first
REAL                    :: tmpArray(0:NKnots-1,0:NKnots-1)
REAL                    :: sDxiN!,DxiN1
!===================================================================================================================================


! init, first layer (constant)
! zero order
n=0
tmpArray=0.
last=nknots-1
DO i=0,last
  IF((knots(i).LE.Xi_in).AND.(Xi_in.LT.knots(i+1)))THEN
    tmpArray(i,n)=1.0
  END IF ! select
END DO ! i

DO n=1,N_in
  last=last-1
  DO i=0,last
    ! standard
!    tmpArray(i,n) = (Xi_in-knots(i))/(knots(i+n)-knots(i))*tmpArray(i,n-1) &
!                  + (knots(i+n+1)-Xi_in)/(knots(i+n+1)-knots(i+1))*tmpArray(i+1,n-1)
    ! optimized
    sDxiN=0.5/REAL(n)
    tmpArray(i,n) = sDxiN*( (Xi_in-knots(i) )*tmparray(i,n-1)+(knots(i+n+1)-Xi_in)*tmpArray(i+1,n-1))
  END DO ! i
END DO ! n

! move back to correct range
UBSpline(0:N_in)=tmpArray(0:N_in,N_in)

END SUBROUTINE DeBoorRef


#ifdef donotcompilethis
SUBROUTINE ComputeGaussDistance(N_In,scaleR,X_in,Elem_xGP,GaussDistance)
!----------------------------------------------------------------------------------------------------------------------------------!
! compute all distance between X_in and given array
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: N_In
REAL,INTENT(IN)         :: X_in(1:3)
REAL,INTENT(IN)         :: scaleR
!REAL,INTENT(IN)         :: Elem_xGP(1:3,1:N_In)
REAL,INTENT(IN)         :: Elem_xGP(1:3,0:N_in,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!REAL,INTENT(OUT)        :: GaussDistance(1:N_In)
REAL,INTENT(OUT)        :: GaussDistance(0:N_in,0:N_in,0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k
!===================================================================================================================================

DO k=0,N_in
  DO j=0,N_in
    DO i=0,N_in
      !tmp=X_in-Elem_xGP(:,i,j,k)
      GaussDistance(i,j,k) =((X_in(1)-Elem_xGP(1,i,j,k))*(X_in(1)-Elem_xGP(1,i,j,k)) &
                            +(X_in(2)-Elem_xGP(2,i,j,k))*(X_in(2)-Elem_xGP(2,i,j,k)) &
                            +(X_in(3)-Elem_xGP(3,i,j,k))*(X_in(3)-Elem_xGP(3,i,j,k)))*scaleR
    END DO ! i=0,N_in
  END DO !  j=0,N_in
END DO ! k=0,N_in


!DO i=1,N_in
!  tmp=X_in-Elem_xGP(:,i)
!  GaussDistance(i) = DOT_PRODUCT(tmp,tmp)*scaleR
!END DO ! i = 1,N_in

END SUBROUTINE ComputeGaussDistance
#endif


SUBROUTINE calcSfSource(SourceSize_in,ChargeMPF,Vec1,Vec2,Vec3,PartPos,PartIdx,PartVelo,const_opt)
!============================================================================================================================
! deposit charges on DOFs via shapefunction including periodic displacements and mirroring with SFdepoFixes
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars,           ONLY:r_sf,DepositionType
USE MOD_PICDepo_Vars,           ONLY:NbrOfSFdepoFixes,SFdepoFixesGeo,SFdepoFixesBounds,SFdepoFixesChargeMult
USE MOD_PICDepo_Vars,           ONLY:SFdepoFixesPartOfLink,SFdepoFixesEps,NbrOfSFdepoFixLinks,SFdepoFixLinks
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:casematrix,NbrOfCases
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: SourceSize_in,PartIdx
REAL, INTENT(IN)                 :: ChargeMPF,PartPos(3),Vec1(3),Vec2(3),Vec3(3)
REAL, INTENT(IN), OPTIONAL       :: PartVelo(3)
LOGICAL, INTENT(IN), OPTIONAL    :: const_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if (defined (PP_HDG) && (PP_nVar==1))
!yes, PartVelo and SourceSize_in are not used, but the subroutine-call and -head would be ugly with the preproc-flags...
INTEGER, PARAMETER               :: SourceSize=1
REAL                             :: Fac(4:4), Fac2(4:4)
#else
INTEGER                          :: SourceSize
REAL                             :: Fac(4-SourceSize_in+1:4), Fac2(4-SourceSize_in+1:4)
#endif
INTEGER                          :: iCase, ind
REAL                             :: ShiftedPart(1:3), caseShiftedPart(1:3), n_loc(1:3)
INTEGER                          :: iSFfix, LinkLoopEnd(2), iSFfixLink, iTwin, iLinkRecursive, SFfixIdx, SFfixIdx2
LOGICAL                          :: DoCycle, DoNotDeposit
REAL                             :: SFfixDistance, SFfixDistance2
LOGICAL , ALLOCATABLE            :: SFdepoFixDone(:)
LOGICAL                          :: const
!----------------------------------------------------------------------------------------------------------------------------------
IF (PRESENT(const_opt)) THEN
  const=const_opt
ELSE
  const=.FALSE.
END IF
#if !(defined (PP_HDG) && (PP_nVar==1))
SourceSize=SourceSize_in
#endif
IF (SourceSize.EQ.1) THEN
  Fac2= ChargeMPF
#if !(defined (PP_HDG) && (PP_nVar==1))
ELSE IF (SourceSize.EQ.4) THEN
  Fac2(1:3) = PartVelo*ChargeMPF
  Fac2(4)= ChargeMPF
#endif
ELSE
  CALL abort(&
__STAMP__ &
,'SourceSize has to be either 1 or 4!',SourceSize)
END IF
IF (NbrOfSFdepoFixes.EQ.0) THEN
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      ShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    Fac = Fac2
    SELECT CASE(TRIM(DepositionType))
    CASE('shape_function')
      CALL depoChargeOnDOFs_sf(ShiftedPart,SourceSize,Fac,const)
    CASE('shape_function_simple')
      CALL depoChargeOnDOFs_sf_simple(ShiftedPart,SourceSize,Fac,const)
    END SELECT
  END DO ! iCase (periodicity)
ELSE ! NbrOfSFdepoFixes.NE.0
  ALLOCATE(SFdepoFixDone(0:NbrOfSFdepoFixes))
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      caseShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    SFdepoFixDone=.FALSE.
    DO iSFfix=0,NbrOfSFdepoFixes
      IF (SFdepoFixesPartOfLink(iSFfix)) CYCLE !this SFfix will be already covered as part of a FixLink
      IF (iSFfix.EQ.0) THEN
        LinkLoopEnd(1)=NbrOfSFdepoFixLinks !non-mirrored position: consider FixLinks
      ELSE
        LinkLoopEnd(1)=0
      END IF
      DO iSFfixLink=0,LinkLoopEnd(1)
        DoCycle=.FALSE.
        !-- strategy for SFdepoFixLink:
        !-- 1.: create 2 identities by iTwin for iLinkRecursive=0 (non-mirror and mirror at SFdepoFixLinks(*,1))
        !-- 2.: mirror recursevely further for iLinkRecursive>0 (alternating at SFdepoFixLinks(*,2) and SFdepoFixLinks(*,1))
        DO iTwin=1,2
          IF (iSFfixLink.EQ.0 .AND. iTwin.EQ.2) EXIT !no SFdepoFixLink
          Fac = Fac2
          ShiftedPart=caseShiftedPart
          IF (iSFfixLink.EQ.0) THEN
            LinkLoopEnd(2)=0 !no SFdepoFixLink
          ELSE
            LinkLoopEnd(2)=ABS(SFdepoFixLinks(iSFfixLink,3))-1 !might be negative as flag for special case
          END IF
          DO iLinkRecursive=0,LinkLoopEnd(2)
            IF (iLinkRecursive.EQ.0) THEN !(first!)
              IF (iTwin.EQ.1) THEN
                SFfixIdx = iSFfix
              ELSE
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            ELSE
              IF (MOD(iLinkRecursive,2).NE.0) THEN !uneven (second and later: 1, 3, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,2)
              ELSE !even (third and later: 2, 4, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            END IF
            DoNotDeposit=.FALSE.
            IF (iLinkRecursive.EQ.0 .OR. (iLinkRecursive.EQ.1 .AND. iTwin.EQ.1)) THEN !single or one-time-mirrored identity
              IF (SFdepoFixDone(SFfixIdx)) DoNotDeposit=.TRUE. !do not deposit this charge (but save position for further mirroring)
              SFdepoFixDone(SFfixIdx) = .TRUE.
            ELSE IF (iSFfixLink.GT.0) THEN
              IF (SFdepoFixLinks(iSFfixLink,3).LT.0) EXIT !skip double mirroring for 120 deg case
            END IF
            IF (SFfixIdx.GT.0) THEN
              IF ( SFdepoFixesBounds(SFfixIdx,1,1).GT.ShiftedPart(1) .OR. ShiftedPart(1).GT.SFdepoFixesBounds(SFfixIdx,2,1) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,2).GT.ShiftedPart(2) .OR. ShiftedPart(2).GT.SFdepoFixesBounds(SFfixIdx,2,2) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,3).GT.ShiftedPart(3) .OR. ShiftedPart(3).GT.SFdepoFixesBounds(SFfixIdx,2,3) ) THEN
                DoCycle=.TRUE.
                EXIT !do not shift this particle (-> go to next single SFfixIdx -> iSFfixLink-loop)
              END IF
              SFfixDistance = SFdepoFixesGeo(SFfixIdx,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx,1,1)) &
                + SFdepoFixesGeo(SFfixIdx,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx,1,2)) &
                + SFdepoFixesGeo(SFfixIdx,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx,1,3))
              SFfixIdx2=0 !init
              IF (SFfixDistance .GT. SFdepoFixesEps) THEN
                IPWRITE(*,'(I4,A,3(x,E12.5))')' original case-pos:',caseShiftedPart
                IPWRITE(*,'(I4,A,3(x,E12.5))')' current pos:',ShiftedPart
                IPWRITE(*,'(I4,7(A,I0))') &
                  ' iCase: ',iCase,', iSFfix:',iSFfix,', iSFfixLink:',iSFfixLink,', iTwin:',iTwin,', iLinkRec:',iLinkRecursive,&
                  ', SFfixIdx:',SFfixIdx,', PartIdx:',PartIdx
                CALL abort(&
                  __STAMP__ &
                  ,'Particle is outside of SF-Fix-Plane! (For Layer-/Resample-Parts: try -UseFixBounds)',SFfixIdx,SFfixDistance)
              ELSE IF ( (SFfixDistance.LT.-r_sf) ) THEN !.OR. (SFfixDistance.GT.0.) ) THEN !SFfixDistance>0 are particle within eps
                DoNotDeposit=.TRUE. !too far inside so that mirrored SF would not reach any DOF
              ELSE IF (iSFfixLink.GT.0) THEN !check which is the other SFfixIdx of current link
                IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,1)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,2)
                ELSE IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,2)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,1)
                ELSE
                  CALL abort(&
                    __STAMP__ &
                    ,'Something is wrong with iSFfixLink',iSFfixLink)
                END IF
              END IF
              ShiftedPart(1:3) = ShiftedPart(1:3) - 2.*SFfixDistance*SFdepoFixesGeo(SFfixIdx,2,1:3)
              Fac = Fac * SFdepoFixesChargeMult(SFfixIdx)
#if !(defined (PP_HDG) && (PP_nVar==1))
              IF (SourceSize.EQ.4) THEN
                ! change velocity
                n_loc = SFdepoFixesGeo(SFfixIdx,2,1:3)
                Fac(1:3) = Fac2(1:3) -2.*DOT_PRODUCT(Fac2(1:3),n_loc)*n_loc
              END IF
#endif
              IF (SFfixIdx2.NE.0) THEN !check if new position would not reach a dof because of the other plane
                SFfixDistance2 = SFdepoFixesGeo(SFfixIdx2,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx2,1,1)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx2,1,2)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx2,1,3))
                IF (SFfixDistance2 .GT. r_sf) DoNotDeposit=.TRUE. !too far outside of plane
              END IF
            END IF
            IF (DoNotDeposit) CYCLE !(-> do not deposit but save position for possible further recursive mirroring)
            !------------- actual deposition:
            SELECT CASE(TRIM(DepositionType))
            CASE('shape_function')
              CALL depoChargeOnDOFs_sf(ShiftedPart,SourceSize,Fac,const)
            CASE('shape_function_simple')
              CALL depoChargeOnDOFs_sf_simple(ShiftedPart,SourceSize,Fac,const)
            END SELECT
          END DO ! iLinkRecursive
          IF (DoCycle) EXIT
        END DO ! iTwin
      END DO ! iSFfixLink
    END DO ! iSFfix
  END DO ! iCase (periodicity)
END IF !NbrOfSFdepoFixes

END SUBROUTINE calcSfSource


SUBROUTINE depoChargeOnDOFs_sf(Position,SourceSize,Fac,const)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars,           ONLY:PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf, ElemDepo_xGP, PartSourceConst
USE MOD_Mesh_Vars,              ONLY:nElems
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_PreProc,                ONLY:PP_N
#ifdef MPI
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
#if (defined (PP_HDG) && (PP_nVar==1))
REAL, INTENT(IN)                 :: Fac(4:4)
#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
#endif
LOGICAL, INTENT(IN)              :: const
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
LOGICAL                          :: chargedone(1:nElems)
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: ElemID
REAL                             :: radius2, S, S1
REAL                             :: dx,dy,dz
INTEGER                          :: expo
!----------------------------------------------------------------------------------------------------------------------------------

chargedone(:) = .FALSE.
!-- determine which background mesh cells (and interpolation points within) need to be considered
kmax = CEILING((Position(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
kmax = MIN(kmax,GEO%FIBGMimax)
kmin = FLOOR((Position(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
kmin = MAX(kmin,GEO%FIBGMimin)
lmax = CEILING((Position(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
lmax = MIN(lmax,GEO%FIBGMjmax)
lmin = FLOOR((Position(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
lmin = MAX(lmin,GEO%FIBGMjmin)
mmax = CEILING((Position(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
mmax = MIN(mmax,GEO%FIBGMkmax)
mmin = FLOOR((Position(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
mmin = MAX(mmin,GEO%FIBGMkmin)
DO kk = kmin,kmax
  DO ll = lmin, lmax
    DO mm = mmin, mmax
      !--- go through all mapped elements not done yet
      DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
        ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
        IF(ElemID.GT.nElems) CYCLE
        IF (.NOT.chargedone(ElemID)) THEN
#ifdef MPI
          nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
          !--- go through all gauss points
          DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
            !-- calculate distance between gauss and particle
            dX = ABS(Position(1) - ElemDepo_xGP(1,k,l,m,ElemID))
            IF(dX.GT.r_sf) CYCLE
            dY = ABS(Position(2) - ElemDepo_xGP(2,k,l,m,ElemID))
            IF(dY.GT.r_sf) CYCLE
            dZ = ABS(Position(3) - ElemDepo_xGP(3,k,l,m,ElemID))
            IF(dZ.GT.r_sf) CYCLE
            radius2 = dX*dX+dY*dY+dZ*dZ
            !-- calculate charge and current density at ip point using a shape function
            !-- currently only one shapefunction available, more to follow (including structure change)
            IF (radius2 .LT. r2_sf) THEN
              S = 1. - r2_sf_inv * radius2
              S1 = S*S
              DO expo = 3, alpha_sf
                S1 = S*S1
              END DO
              IF (const) THEN
                IF (SourceSize.EQ.1) THEN
                  PartSourceConst(4,k,l,m,ElemID) = PartSourceConst(4,k,l,m,ElemID) + Fac(4) * S1
#if !(defined (PP_HDG) && (PP_nVar==1))
                ELSE IF (SourceSize.EQ.4) THEN
                  PartSourceConst(1:4,k,l,m,ElemID) = PartSourceConst(1:4,k,l,m,ElemID) + Fac(1:4) * S1
#endif
                END IF
              ELSE !.NOT.const
                IF (SourceSize.EQ.1) THEN
                  PartSource(4,k,l,m,ElemID) = PartSource(4,k,l,m,ElemID) + Fac(4) * S1
#if !(defined (PP_HDG) && (PP_nVar==1))
                ELSE IF (SourceSize.EQ.4) THEN
                  PartSource(1:4,k,l,m,ElemID) = PartSource(1:4,k,l,m,ElemID) + Fac(1:4) * S1
#endif
                END IF
              END IF !const.
            END IF
          END DO; END DO; END DO
          chargedone(ElemID) = .TRUE.
        END IF
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

END SUBROUTINE depoChargeOnDOFs_sf


SUBROUTINE depoChargeOnDOFs_sf_simple(Position,SourceSize,Fac,const)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction_simple (i.e. loop through all elems instead of part-dependency: efficient for small elem-nbr!)
!============================================================================================================================
! use MODULES
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
USE MOD_PICDepo_Vars,           ONLY:PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf, ElemDepo_xGP, ElemRadius2_sf, PartSourceConst
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo
USE MOD_PreProc,                ONLY:PP_N, PP_nElems
#ifdef MPI
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
#if (defined (PP_HDG) && (PP_nVar==1))
REAL, INTENT(IN)                 :: Fac(4:4)
#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
#endif
LOGICAL, INTENT(IN)              :: const
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: ElemID
REAL                             :: radius2, S, S1
REAL                             :: dx,dy,dz
INTEGER                          :: expo
!----------------------------------------------------------------------------------------------------------------------------------

DO ElemID=1,PP_nElems
  dX = ABS(Position(1) - ElemBaryNgeo(1,ElemID))
  IF(dX.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  dY = ABS(Position(2) - ElemBaryNgeo(2,ElemID))
  IF(dY.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  dZ = ABS(Position(3) - ElemBaryNgeo(3,ElemID))
  IF(dZ.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  radius2 = dX*dX+dY*dY+dZ*dZ
  IF(radius2.GT.ElemRadius2_sf(ElemID)) CYCLE
#ifdef MPI
  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*MPI*/
  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    !-- calculate distance between gauss and particle
    dX = ABS(Position(1) - ElemDepo_xGP(1,k,l,m,ElemID))
    IF(dX.GT.r_sf) CYCLE
    dY = ABS(Position(2) - ElemDepo_xGP(2,k,l,m,ElemID))
    IF(dY.GT.r_sf) CYCLE
    dZ = ABS(Position(3) - ElemDepo_xGP(3,k,l,m,ElemID))
    IF(dZ.GT.r_sf) CYCLE
    radius2 = dX*dX+dY*dY+dZ*dZ
    !-- calculate charge and current density at ip point using a shape function
    !-- currently only one shapefunction available, more to follow (including structure change)
    IF (radius2 .GT. r2_sf) CYCLE
    S = 1. - r2_sf_inv * radius2
    S1 = S*S
    DO expo = 3, alpha_sf
      S1 = S*S1
    END DO
    IF (const) THEN
      IF (SourceSize.EQ.1) THEN
        PartSourceConst(4,k,l,m,ElemID) = PartSourceConst(4,k,l,m,ElemID) + Fac(4) * S1
#if !(defined (PP_HDG) && (PP_nVar==1))
      ELSE IF (SourceSize.EQ.4) THEN
        PartSourceConst(1:4,k,l,m,ElemID) = PartSourceConst(1:4,k,l,m,ElemID) + Fac(1:4) * S1
#endif
      END IF
    ELSE !.NOT.const
      IF (SourceSize.EQ.1) THEN
        PartSource(4,k,l,m,ElemID) = PartSource(4,k,l,m,ElemID) + Fac(4) * S1
#if !(defined (PP_HDG) && (PP_nVar==1))
      ELSE IF (SourceSize.EQ.4) THEN
        PartSource(1:4,k,l,m,ElemID) = PartSource(1:4,k,l,m,ElemID) + Fac(1:4) * S1
#endif
      END IF
    END IF !const
  END DO; END DO; END DO
END DO !ElemID=1,PP_nElems

END SUBROUTINE depoChargeOnDOFs_sf_simple


END MODULE MOD_PICDepo
