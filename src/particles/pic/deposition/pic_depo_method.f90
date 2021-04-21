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

MODULE MOD_PICDepo_Method
!===================================================================================================================================
! Module containing the different deposition methods (NGP, linear (inter-cell) weighting, shape function
!===================================================================================================================================
IMPLICIT NONE
PRIVATE

INTERFACE PartRHS
  PROCEDURE DepositionMethod
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DepositionMethod
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE DepositionMethodInterface(doParticle_In)
    USE MOD_Particle_Vars ,ONLY: PDM
    LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
  END SUBROUTINE
END INTERFACE

PROCEDURE(DepositionMethodInterface),POINTER :: DepositionMethod    !< pointer defining the standard inner Riemann solver

INTEGER,PARAMETER      :: PRM_DEPO_SF   = 0  ! shape_function
INTEGER,PARAMETER      :: PRM_DEPO_SF_CC= 1  ! shape_function_cc
INTEGER,PARAMETER      :: PRM_DEPO_SF_ADAPTIVE= 2  ! shape_function_adaptive
INTEGER,PARAMETER      :: PRM_DEPO_CVW  = 6  ! cell_volweight
INTEGER,PARAMETER      :: PRM_DEPO_CVWM = 12 ! cell_volweight_mean

INTERFACE InitDepositionMethod
  MODULE PROCEDURE InitDepositionMethod
END INTERFACE

PUBLIC :: InitDepositionMethod
!==================================================================================================================================

PUBLIC :: DefineParametersDepositionMethod
CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("PIC Deposition")
CALL prms%CreateLogicalOption('PIC-DoDeposition', 'Switch deposition of charge (and current density) on/off', '.TRUE.')

CALL prms%CreateIntFromStringOption('PIC-Deposition-Type', "Type/Method used in the deposition step: \n"           //&
                                    '1.1)  shape_function ('//TRIM(int2strf(PRM_DEPO_SF))//')\n'                   //&
                                    '1.2)  shape_function_cc ('//TRIM(int2strf(PRM_DEPO_SF_CC))//')\n'             //&
                                    '1.3)  shape_function_adaptive ('//TRIM(int2strf(PRM_DEPO_SF_ADAPTIVE))//')\n' //&
                                    '2.)   cell_volweight ('//TRIM(int2strf(PRM_DEPO_CVW))//')\n'                  //&
                                    '3.)   cell_volweight_mean ('//TRIM(int2strf(PRM_DEPO_CVWM))//')'                &
                                    ,'cell_volweight')

CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function'             , PRM_DEPO_SF)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_cc'          , PRM_DEPO_SF_CC)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_adaptive'    , PRM_DEPO_SF_ADAPTIVE)
CALL addStrListEntry('PIC-Deposition-Type' , 'cell_volweight'             , PRM_DEPO_CVW)
CALL addStrListEntry('PIC-Deposition-Type' , 'cell_volweight_mean'        , PRM_DEPO_CVWM)
END SUBROUTINE DefineParametersDepositionMethod


!==================================================================================================================================!
!> Initialize deposition method function pointer
!==================================================================================================================================!
SUBROUTINE InitDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools            ,ONLY: GETINTFROMSTR
USE MOD_PICDepo_Vars           ,ONLY: DepositionType,r_sf,dim_sf,dim_sf_dir,SFAdaptiveSmoothing,alpha_sf,sfDepo3D
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_ReadInTools            ,ONLY: GETREAL,PrintOption,GETINT,GETLOGICAL
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: DepositionType_loc
!==================================================================================================================================
r_sf=-1.0 ! default
DepositionType_loc = GETINTFROMSTR('PIC-Deposition-Type')
! check for interpolation type incompatibilities (cannot be done at interpolation_init
! because DepositionType_loc is not known yet)
IF((DepositionType_loc.EQ.PRM_DEPO_CVWM).AND.(TrackingMethod.NE.TRIATRACKING)) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in pic_depo.f90: PIC-Deposition-Type = cell_volweight_mean only allowed with TriaTracking!')
END IF

! Select the deposition method function pointer
SELECT CASE(DepositionType_loc)
Case(PRM_DEPO_SF) ! shape_function
  DepositionMethod => DepositionMethod_SF
  DepositionType   = 'shape_function'
Case(PRM_DEPO_SF_CC) ! shape_function
  DepositionMethod => DepositionMethod_SF
  DepositionType   = 'shape_function_cc'
Case(PRM_DEPO_SF_ADAPTIVE) ! shape_function
  DepositionMethod => DepositionMethod_SF
  DepositionType   = 'shape_function_adaptive'
Case(PRM_DEPO_CVW) ! cell_volweight
  DepositionType   = 'cell_volweight'
  DepositionMethod => DepositionMethod_CVW
Case(PRM_DEPO_CVWM) ! cell_volweight_mean
  DepositionType   = 'cell_volweight_mean'
  DepositionMethod => DepositionMethod_CVWM
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
      'Unknown DepositionMethod!' ,IntInfo=DepositionType_loc)
END SELECT

! If shape function is used, the radius must be read here as it is used for the BGM setup
IF(StringBeginsWith(DepositionType,'shape_function'))THEN
  IF(TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
    ! When using shape function adaptive, the radius is scaled as such that only the direct element neighbours are considered for
    ! deposition (all corner node connected elements) and each element has a separate shape function radius. Therefore, the global
    ! radius is set to zero
    r_sf = 0.
    CALL PrintOption('Global shape function radius is set to zero: PIC-shapefunction-radius' , 'INFO.' , RealOpt=r_sf)
    SFAdaptiveSmoothing = GETLOGICAL('PIC-shapefunction-adaptive-smoothing')
  ELSE
    r_sf = GETREAL('PIC-shapefunction-radius')
  END IF ! TRIM(DepositionType).EQ.'shape_function_adaptive'

  ! Get dimension of shape function kernel (distributes in 1, 2 or 3 dimensions)
  dim_sf   = GETINT('PIC-shapefunction-dimension')

  ! Get shape function direction for 1D (the direction in which the charge will be distributed) and 2D (the direction in which the
  ! charge will be constant)
  IF(dim_sf.NE.3) dim_sf_dir = GETINT('PIC-shapefunction-direction')

  ! Get shape function exponent and dimension (1D, 2D or 3D). This parameter is required in InitShapeFunctionDimensionalty()
  alpha_sf = GETINT('PIC-shapefunction-alpha')

  ! Get deposition parameter, the default is TRUE (3D), that distributes the charge over
  !  FALSE: line (1D) / area (2D)
  !   TRUE: volume (3D)
  sfDepo3D = GETLOGICAL('PIC-shapefunction-3D-deposition')
  IF((dim_sf.EQ.3).AND.(.NOT.sfDepo3D)) &
      CALL abort(__STAMP__,'PIC-shapefunction-dimension=F and PIC-shapefunction-3D-deposition=T is not allowed')
END IF ! StringBeginsWith(DepositionType,'shape_function')

! Suppress compiler warnings
RETURN
CALL DepositionMethod_SF()
CALL DepositionMethod_CVW()
CALL DepositionMethod_CVWM()

END SUBROUTINE InitDepositionMethod


SUBROUTINE DepositionMethod_CVW(doParticle_In)
!===================================================================================================================================
! 'cell_volweight'
! Linear charge density distribution within a cell (discontinuous across cell interfaces)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies,PDM,PEM,PartPosRef,usevMPF,PartMPF
USE MOD_Particle_Vars          ,ONLY: PartState
USE MOD_PICDepo_Vars           ,ONLY: PartSource,CellVolWeight_Volumes,CellVolWeightFac
USE MOD_Part_Tools             ,ONLY: isDepositParticle
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers     ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars              ,ONLY: nElems, offSetElem
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
#if ((USE_HDG) && (PP_nVar==1))
USE MOD_TimeDisc_Vars          ,ONLY: dt,tAnalyzeDiff,tEndDiff
#endif
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: PartSource_Shared_Win
USE MOD_Globals            ,ONLY: IERROR
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE  :: BGMSourceCellVol(:,:,:,:,:)
REAL               :: Charge, TSource(1:4)
REAL               :: alpha1, alpha2, alpha3, TempPartPos(1:3)
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER :: SourceDim=1
LOGICAL, PARAMETER :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL            :: doCalculateCurrentDensity
INTEGER            :: SourceDim
#endif
INTEGER            :: kk, ll, mm
INTEGER            :: iPart,iElem
!===================================================================================================================================

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

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

ALLOCATE(BGMSourceCellVol(SourceDim:4,0:1,0:1,0:1,1:nElems))
BGMSourceCellVol(:,:,:,:,:) = 0.0
DO iPart = 1,PDM%ParticleVecLength
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(PRESENT(doParticle_In))THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  IF (usevMPF) THEN
    Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
  ELSE
    Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF ! usevMPF
  IF(TrackingMethod.EQ.REFMAPPING)THEN
    TempPartPos(1:3)=PartPosRef(1:3,iPart)
  ELSE
    CALL GetPositionInRefElem(PartState(1:3,iPart),TempPartPos,PEM%GlobalElemID(iPart),ForceMode=.TRUE.)
  END IF
  IF(doCalculateCurrentDensity)THEN
    TSource(1:3) = PartState(4:6,iPart)*Charge
  ELSE
    TSource(1:3) = 0.0
  END IF
  iElem = PEM%LocalElemID(iPart)
  TSource(4) = Charge
  alpha1=(TempPartPos(1)+1.0)/2.0
  alpha2=(TempPartPos(2)+1.0)/2.0
  alpha3=(TempPartPos(3)+1.0)/2.0
  BGMSourceCellVol(:,0,0,0,iElem) = BGMSourceCellVol(:,0,0,0,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
  BGMSourceCellVol(:,0,0,1,iElem) = BGMSourceCellVol(:,0,0,1,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(alpha3))
  BGMSourceCellVol(:,0,1,0,iElem) = BGMSourceCellVol(:,0,1,0,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(alpha2)*(1-alpha3))
  BGMSourceCellVol(:,0,1,1,iElem) = BGMSourceCellVol(:,0,1,1,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(alpha2)*(alpha3))
  BGMSourceCellVol(:,1,0,0,iElem) = BGMSourceCellVol(:,1,0,0,iElem) + (TSource(SourceDim:4)*(alpha1)*(1-alpha2)*(1-alpha3))
  BGMSourceCellVol(:,1,0,1,iElem) = BGMSourceCellVol(:,1,0,1,iElem) + (TSource(SourceDim:4)*(alpha1)*(1-alpha2)*(alpha3))
  BGMSourceCellVol(:,1,1,0,iElem) = BGMSourceCellVol(:,1,1,0,iElem) + (TSource(SourceDim:4)*(alpha1)*(alpha2)*(1-alpha3))
  BGMSourceCellVol(:,1,1,1,iElem) = BGMSourceCellVol(:,1,1,1,iElem) + (TSource(SourceDim:4)*(alpha1)*(alpha2)*(alpha3))
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to iElem
#endif /*USE_LOADBALANCE*/
END DO

DO iElem=1, nElems
  BGMSourceCellVol(:,0,0,0,iElem) = BGMSourceCellVol(:,0,0,0,iElem)/CellVolWeight_Volumes(0,0,0,iElem)
  BGMSourceCellVol(:,0,0,1,iElem) = BGMSourceCellVol(:,0,0,1,iElem)/CellVolWeight_Volumes(0,0,1,iElem)
  BGMSourceCellVol(:,0,1,0,iElem) = BGMSourceCellVol(:,0,1,0,iElem)/CellVolWeight_Volumes(0,1,0,iElem)
  BGMSourceCellVol(:,0,1,1,iElem) = BGMSourceCellVol(:,0,1,1,iElem)/CellVolWeight_Volumes(0,1,1,iElem)
  BGMSourceCellVol(:,1,0,0,iElem) = BGMSourceCellVol(:,1,0,0,iElem)/CellVolWeight_Volumes(1,0,0,iElem)
  BGMSourceCellVol(:,1,0,1,iElem) = BGMSourceCellVol(:,1,0,1,iElem)/CellVolWeight_Volumes(1,0,1,iElem)
  BGMSourceCellVol(:,1,1,0,iElem) = BGMSourceCellVol(:,1,1,0,iElem)/CellVolWeight_Volumes(1,1,0,iElem)
  BGMSourceCellVol(:,1,1,1,iElem) = BGMSourceCellVol(:,1,1,1,iElem)/CellVolWeight_Volumes(1,1,1,iElem)
END DO

DO iElem = 1, nElems
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
        alpha1 = CellVolWeightFac(kk)
        alpha2 = CellVolWeightFac(ll)
        alpha3 = CellVolWeightFac(mm)
        PartSource(SourceDim:4,kk,ll,mm,GetCNElemID(iElem+offSetElem)) = &
            PartSource(SourceDim:4,kk,ll,mm,GetCNElemID(iElem+offSetElem)) + &
            BGMSourceCellVol(:,0,0,0,iElem) * (1-alpha1) * (1-alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,0,0,1,iElem) * (1-alpha1) * (1-alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,0,1,0,iElem) * (1-alpha1) *   (alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,0,1,1,iElem) * (1-alpha1) *   (alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,1,0,0,iElem) *   (alpha1) * (1-alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,1,0,1,iElem) *   (alpha1) * (1-alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,1,1,0,iElem) *   (alpha1) *   (alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,1,1,1,iElem) *   (alpha1) *   (alpha2) *   (alpha3)
      END DO ! mm
    END DO ! ll
  END DO ! kk
END DO ! iElem
#if USE_LOADBALANCE
CALL LBElemSplitTime_avg(tLBStart) ! Average over the number of elems (and Start again)
#endif /*USE_LOADBALANCE*/
DEALLOCATE(BGMSourceCellVol)
#if USE_MPI
CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
#endif /*USE_MPI*/
END SUBROUTINE DepositionMethod_CVW


SUBROUTINE DepositionMethod_CVWM(doParticle_In)
!===================================================================================================================================
! 'cell_volweight_mean'
! Linear charge density distribution within a cell (continuous across cell interfaces)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars    ,ONLY: DoDielectricSurfaceCharge
USE MOD_Eval_xyz           ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars          ,ONLY: nElems,OffsetElem
USE MOD_Particle_Vars      ,ONLY: Species,PartSpecies,PDM,PEM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared, nUniqueGlobalNodes, NodeInfo_Shared
USE MOD_PICDepo_Vars       ,ONLY: PartSource,CellVolWeightFac,NodeSourceExt,NodeVolume,NodeSource
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Part_Tools         ,ONLY: isDepositParticle
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceLoc, NodeMapping, NodeSource_Shared_Win
USE MOD_MPI_Shared_Vars    ,ONLY: MPI_COMM_LEADERS_SHARED, MPI_COMM_SHARED, myComputeNodeRank, myLeaderGroupRank
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeProcessors, nLeaderGroupProcs
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmpLoc
USE MOD_PICDepo_Vars       ,ONLY: PartSource_Shared_Win
#else
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmp
#endif  /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime,LBElemSplitTime,LBElemPauseTime_avg
#endif /*USE_LOADBALANCE*/
#if ((USE_HDG) && (PP_nVar==1))
USE MOD_TimeDisc_Vars      ,ONLY: dt,tAnalyzeDiff,tEndDiff
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Charge, TSource(1:4)
REAL               :: alpha1, alpha2, alpha3, TempPartPos(1:3)
INTEGER            :: kk, ll, mm, iPart, iElem, iProc
INTEGER            :: NodeID(1:8), firstElem, lastElem, firstNode, lastNode, iNode
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER :: SourceDim=1
LOGICAL, PARAMETER :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL            :: doCalculateCurrentDensity
INTEGER            :: SourceDim
#endif
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
#if USE_MPI
INTEGER            :: RecvRequest(0:nLeaderGroupProcs-1),SendRequest(0:nLeaderGroupProcs-1)
INTEGER            :: MessageSize
#endif
!===================================================================================================================================
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

! Check whether charge and current density have to be computed or just the charge density
! For HDG the current density is only required for output to HDF5, i.e., analysis reasons
#if ((USE_HDG) && (PP_nVar==1))
IF(ALMOSTEQUAL(dt,tAnalyzeDiff).OR.ALMOSTEQUAL(dt,tEndDiff))THEN
  doCalculateCurrentDensity=.TRUE.
  SourceDim=1
ELSE ! do not calculate current density
  doCalculateCurrentDensity=.FALSE.
  SourceDim=4
END IF

! Quick-fix for multi-node
doCalculateCurrentDensity=.TRUE.
SourceDim=1
#endif

#if USE_MPI
ASSOCIATE(NodeSource       => NodeSourceLoc       ,&
          NodeSourceExtTmp => NodeSourceExtTmpLoc )
  firstNode = INT(REAL( myComputeNodeRank   *nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))+1
  lastNode  = INT(REAL((myComputeNodeRank+1)*nUniqueGlobalNodes)/REAL(nComputeNodeProcessors))
#else
  firstNode = 1
  lastNode = nUniqueGlobalNodes
#endif

  NodeSource=0.0
  DO iPart=1,PDM%ParticleVecLength
    IF(PRESENT(doParticle_In))THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    IF (isDepositParticle(iPart)) THEN
      IF (usevMPF) THEN
        Charge = Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
      ELSE
        Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
      END IF
      CALL GetPositionInRefElem(PartState(1:3,iPart),TempPartPos(1:3),PEM%GlobalElemID(iPart),ForceMode=.TRUE.)
      TSource(:) = 0.0
      IF(doCalculateCurrentDensity)THEN
        TSource(1:3) = PartState(4:6,iPart)*Charge
      END IF
      TSource(4) = Charge

      alpha1=0.5*(TempPartPos(1)+1.0)
      alpha2=0.5*(TempPartPos(2)+1.0)
      alpha3=0.5*(TempPartPos(3)+1.0)

      NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,PEM%CNElemID(iPart)))
      NodeSource(SourceDim:4,NodeID(1)) = NodeSource(SourceDim:4,NodeID(1)) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
      NodeSource(SourceDim:4,NodeID(2)) = NodeSource(SourceDim:4,NodeID(2)) + (TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*(1-alpha3))
      NodeSource(SourceDim:4,NodeID(3)) = NodeSource(SourceDim:4,NodeID(3)) + (TSource(SourceDim:4)*  (alpha1)*  (alpha2)*(1-alpha3))
      NodeSource(SourceDim:4,NodeID(4)) = NodeSource(SourceDim:4,NodeID(4)) + (TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*(1-alpha3))
      NodeSource(SourceDim:4,NodeID(5)) = NodeSource(SourceDim:4,NodeID(5)) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*  (alpha3))
      NodeSource(SourceDim:4,NodeID(6)) = NodeSource(SourceDim:4,NodeID(6)) + (TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*  (alpha3))
      NodeSource(SourceDim:4,NodeID(7)) = NodeSource(SourceDim:4,NodeID(7)) + (TSource(SourceDim:4)*  (alpha1)*  (alpha2)*  (alpha3))
      NodeSource(SourceDim:4,NodeID(8)) = NodeSource(SourceDim:4,NodeID(8)) + (TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*  (alpha3))
#if USE_LOADBALANCE
    CALL LBElemSplitTime(PEM%LocalElemID(iPart),tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to iElem
#endif /*USE_LOADBALANCE*/
    END IF
  END DO

#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  ! 1/2 Add the local non-synchronized surface charge contribution (does not consider the charge contribution from restart files) from
  ! NodeSourceExtTmp. This contribution accumulates over time, but remains locally to each processor as it is communicated via the
  ! normal NodeSource container. The synchronized part is added after communication.
  IF(DoDielectricSurfaceCharge)THEN
      NodeSource(4,:) = NodeSource(4,:) + NodeSourceExtTmp(:)
  END IF ! DoDielectricSurfaceCharge
#if USE_LOADBALANCE
  CALL LBElemPauseTime_avg(tLBStart) ! Average over the number of elems
#endif /*USE_LOADBALANCE*/

#if USE_MPI
END ASSOCIATE
MessageSize = (5-SourceDim)*nUniqueGlobalNodes
IF(myComputeNodeRank.EQ.0)THEN
  CALL MPI_REDUCE(NodeSourceLoc(SourceDim:4,:),NodeSource(SourceDim:4,:),MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(NodeSourceLoc(SourceDim:4,:),0                        ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
END IF ! myrank.eq.0
CALL MPI_WIN_SYNC(NodeSource_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

! Multi-node communication
IF(nLeaderGroupProcs.GT.1)THEN
  IF(myComputeNodeRank.EQ.0)THEN

    ! 1) Send/Receive charge density
    DO iProc = 0, nLeaderGroupProcs - 1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE

      ! Open receive buffer
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        CALL MPI_IRECV( NodeMapping(iProc)%RecvNodeSourceCharge(:) &
                  , NodeMapping(iProc)%nRecvUniqueNodes            &
                  , MPI_DOUBLE_PRECISION                           &
                  , iProc                                          &
                  , 666                                            &
                  , MPI_COMM_LEADERS_SHARED                        &
                  , RecvRequest(iProc)                             &
                  , IERROR)
      END IF
      ! Send message (non-blocking)
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        DO iNode = 1, NodeMapping(iProc)%nSendUniqueNodes
          NodeMapping(iProc)%SendNodeSourceCharge(iNode) = NodeSource(4,NodeMapping(iProc)%SendNodeUniqueGlobalID(iNode))
        END DO
        CALL MPI_ISEND( NodeMapping(iProc)%SendNodeSourceCharge(:) &
                      , NodeMapping(iProc)%nSendUniqueNodes        &
                      , MPI_DOUBLE_PRECISION                       &
                      , iProc                                      &
                      , 666                                        &
                      , MPI_COMM_LEADERS_SHARED                    &
                      , SendRequest(iProc)                         &
                      , IERROR)
      END IF
    END DO

    ! Finish communication
    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
      IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
        CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
        IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
    END DO

    ! 2) Send/Receive current density
    IF(doCalculateCurrentDensity)THEN

      DO iProc = 0, nLeaderGroupProcs - 1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE

        ! Open receive buffer
        IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
          CALL MPI_IRECV( NodeMapping(iProc)%RecvNodeSourceCurrent(1:3,:) &
              , 3*NodeMapping(iProc)%nRecvUniqueNodes                     &
              , MPI_DOUBLE_PRECISION                                      &
              , iProc                                                     &
              , 666                                                       &
              , MPI_COMM_LEADERS_SHARED                                   &
              , RecvRequest(iProc)                                        &
              , IERROR)
        END IF
        ! Send message (non-blocking)
        IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
          DO iNode = 1, NodeMapping(iProc)%nSendUniqueNodes
            NodeMapping(iProc)%SendNodeSourceCurrent(1:3,iNode) = NodeSource(1:3,NodeMapping(iProc)%SendNodeUniqueGlobalID(iNode))
          END DO
          CALL MPI_ISEND( NodeMapping(iProc)%SendNodeSourceCurrent(1:3,:) &
              , 3*NodeMapping(iProc)%nSendUniqueNodes                     &
              , MPI_DOUBLE_PRECISION                                      &
              , iProc                                                     &
              , 666                                                       &
              , MPI_COMM_LEADERS_SHARED                                   &
              , SendRequest(iProc)                                        &
              , IERROR)
        END IF
      END DO

      ! Finish communication
      DO iProc = 0,nLeaderGroupProcs-1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (NodeMapping(iProc)%nSendUniqueNodes.GT.0) THEN
          CALL MPI_WAIT(SendRequest(iProc),MPISTATUS,IERROR)
          IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF
        IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
          CALL MPI_WAIT(RecvRequest(iProc),MPISTATUS,IERROR)
          IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
        END IF
      END DO
    END IF ! doCalculateCurrentDensity

    ! 3) Extract messages
    IF(doCalculateCurrentDensity)THEN! SourceDim=1
      DO iProc = 0, nLeaderGroupProcs - 1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
          DO iNode = 1, NodeMapping(iProc)%nRecvUniqueNodes
            ASSOCIATE( NS => NodeSource(SourceDim:4,NodeMapping(iProc)%RecvNodeUniqueGlobalID(iNode)))
              NS = NS + (/NodeMapping(iProc)%RecvNodeSourceCurrent(1:3,iNode), NodeMapping(iProc)%RecvNodeSourceCharge(iNode)/)
            END ASSOCIATE
          END DO
        END IF
      END DO
    ELSE
      DO iProc = 0, nLeaderGroupProcs - 1
        IF (iProc.EQ.myLeaderGroupRank) CYCLE
        IF (NodeMapping(iProc)%nRecvUniqueNodes.GT.0) THEN
          DO iNode = 1, NodeMapping(iProc)%nRecvUniqueNodes
            ASSOCIATE( NS => NodeSource(4,NodeMapping(iProc)%RecvNodeUniqueGlobalID(iNode)))
              NS = NS + NodeMapping(iProc)%RecvNodeSourceCharge(iNode)
            END ASSOCIATE
          END DO
        END IF
      END DO
    END IF ! doCalculateCurrentDensity
  END IF ! myComputeNodeRank.EQ.0
  CALL MPI_WIN_SYNC(NodeSource_Shared_Win,IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
END IF ! nLeaderGroupProcs.GT.1
#endif

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

! 2/2 Add the global, synchronized surface charge contribution (considers the charge contribution from restart files) from
! NodeSourceExt. The container NodeSourceExt is updated when it is written to .h5, where, additionally, the container
! NodeSourceExtTmp is nullified
IF(DoDielectricSurfaceCharge)THEN
  DO iNode=firstNode, lastNode
    NodeSource(4,iNode) = NodeSource(4,iNode) + NodeSourceExt(iNode)
  END DO
END IF ! DoDielectricSurfaceCharge

! Currently also "Nodes" are included in time measurement that is averaged across all elements. Can this be improved?
DO iNode=firstNode, lastNode
  IF(NodeVolume(iNode).GT.0.)THEN
    NodeSource(SourceDim:4,iNode) = NodeSource(SourceDim:4,iNode)/NodeVolume(iNode)
  END IF ! NodeVolume(iNode).GT.0.
END DO
#if USE_LOADBALANCE
CALL LBElemPauseTime_avg(tLBStart) ! Average over the number of elems
#endif /*USE_LOADBALANCE*/

#if USE_MPI
CALL MPI_WIN_SYNC(NodeSource_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)
firstElem = offSetElem+1
lastElem  = offSetElem+nElems
#else
firstElem = 1
lastElem = nElems
#endif

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
! Interpolate node source values to volume polynomial
DO iElem = firstElem, lastElem
  ! Get UniqueNodeID from NonUniqueNodeID = ElemNodeID_Shared(:,GetCNElemID(iElem))
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem)))
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
        alpha1 = CellVolWeightFac(kk)
        alpha2 = CellVolWeightFac(ll)
        alpha3 = CellVolWeightFac(mm)
        Partsource(SourceDim:4,kk,ll,mm,GetCNElemID(iElem)) = &
             NodeSource(SourceDim:4,NodeID(1)) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(2)) * (alpha1)   * (1-alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(3)) * (alpha1)   * (alpha2)   * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(4)) * (1-alpha1) * (alpha2)   * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(5)) * (1-alpha1) * (1-alpha2) * (alpha3)   + &
             NodeSource(SourceDim:4,NodeID(6)) * (alpha1)   * (1-alpha2) * (alpha3)   + &
             NodeSource(SourceDim:4,NodeID(7)) * (alpha1)   * (alpha2)   * (alpha3)   + &
             NodeSource(SourceDim:4,NodeID(8)) * (1-alpha1) * (alpha2)   * (alpha3)
      END DO !mm
    END DO !ll
  END DO !kk
END DO !iEle
#if USE_LOADBALANCE
CALL LBElemPauseTime_avg(tLBStart) ! Average over the number of elems
#endif /*USE_LOADBALANCE*/
#if USE_MPI
CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
#endif /*USE_MPI*/
END SUBROUTINE DepositionMethod_CVWM


SUBROUTINE DepositionMethod_SF(doParticle_In)
!===================================================================================================================================
! 'shape_function'
! Smooth polynomial deposition via "shape functions" of various order in 3D
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_globals
USE MOD_Particle_Vars               ,ONLY: Species, PartSpecies,PDM,PartMPF,usevMPF
USE MOD_Particle_Vars               ,ONLY: PartState
USE MOD_PICDepo_Vars                ,ONLY: PartSource
USE MOD_PICDepo_Shapefunction_Tools ,ONLY: calcSfSource
USE MOD_Mesh_Tools                  ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_PICDepo_Vars                ,ONLY: PartSourceProc
USE MOD_PICDepo_Vars                ,ONLY: PartSource_Shared_Win, ShapeMapping, nSendShapeElems, SendShapeElemID
USE MOD_PICDepo_Vars                ,ONLY: CNShapeMapping
USE MOD_MPI_Shared_Vars             ,ONLY: MPI_COMM_SHARED, myComputeNodeRank, nComputeNodeProcessors
USE MOD_MPI_Shared_Vars             ,ONLY: MPI_COMM_LEADERS_SHARED, myLeaderGroupRank, nLeaderGroupProcs
#endif /*USE_MPI*/
USE MOD_Part_Tools                  ,ONLY: isDepositParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Charge
INTEGER            :: iElem, iPart
#if USE_MPI
INTEGER            :: SendRequest, RecvRequest(nComputeNodeProcessors-1), iProc, CNElemID
INTEGER            :: RecvRequestCN(0:nLeaderGroupProcs-1), SendRequestCN(0:nLeaderGroupProcs-1)
#endif
!===================================================================================================================================
#if USE_MPI
PartSourceProc = 0.
#endif
!Vec1(1:3) = 0.
!Vec2(1:3) = 0.
!Vec3(1:3) = 0.
!IF (GEO%nPeriodicVectors.EQ.1) THEN
!  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!END IF
!IF (GEO%nPeriodicVectors.EQ.2) THEN
!  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!END IF
!IF (GEO%nPeriodicVectors.EQ.3) THEN
!  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
!  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
!  Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
!END IF
!CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
DO iPart=1,PDM%ParticleVecLength
  IF(PRESENT(doParticle_In))THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  IF (.NOT.isDepositParticle(iPart)) CYCLE
  IF (usevMPF) THEN
    Charge = Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
  ELSE
    Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
  END IF
  ! Fill PartSourceProc and deposit charge in local part of PartSource(CNElem(1:nElems + offset))
  CALL calcSfSource(4,Charge,PartState(1:3,iPart),iPart,PartVelo=PartState(4:6,iPart))
END DO
#if USE_MPI
! Communication
CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)

! 1 of 2: Inner-Node Communication
IF (myComputeNodeRank.EQ.0) THEN
  DO iProc = 1,nComputeNodeProcessors-1
      IF (ShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE
      CALL MPI_IRECV( ShapeMapping(iProc)%RecvBuffer(1:4,0:PP_N,0:PP_N,0:PP_N,1:ShapeMapping(iProc)%nRecvShapeElems)&
                    , ShapeMapping(iProc)%nRecvShapeElems*4*(PP_N+1)**3                                   &
                    , MPI_DOUBLE_PRECISION                &
                    , iProc                               &
                    , 2001                                &
                    , MPI_COMM_SHARED                     &
                    , RecvRequest(iProc)                  &
                    , IERROR)
  END DO

  ! Add contributions of node slaves
  DO iProc = 1,nComputeNodeProcessors-1
    IF (ShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE
    CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
    DO iElem = 1, ShapeMapping(iProc)%nRecvShapeElems
      ASSOCIATE( ShapeID => ShapeMapping(iProc)%RecvShapeElemID(iElem))
        PartSource(:,:,:,:,ShapeID) = PartSource(:,:,:,:,ShapeID) + ShapeMapping(iProc)%RecvBuffer(:,:,:,:,iElem)
      END ASSOCIATE
    END DO
  END DO

  ! Add contribution of node root
  DO iElem = 1, nSendShapeElems
    PartSource(:,:,:,:,SendShapeElemID(iElem)) = PartSource(:,:,:,:,SendShapeElemID(iElem)) + PartSourceProc(:,:,:,:,iElem)
  END DO
ELSE
  IF (nSendShapeElems.GT.1) THEN
    CALL MPI_ISEND( PartSourceProc                         &
                  , nSendShapeElems*4*(PP_N+1)**3          &
                  , MPI_DOUBLE_PRECISION                   &
                  , 0                                      &
                  , 2001                                   &
                  , MPI_COMM_SHARED                        &
                  , SendRequest                            &
                  , IERROR)

    CALL MPI_WAIT(SendRequest,MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END IF
END IF
!CALL MPI_WIN_FLUSH(0,PartSource_Shared_Win, IERROR)
CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)

! 2 of 2: Multi-node communication
IF(nLeaderGroupProcs.GT.1)THEN
  IF(myComputeNodeRank.EQ.0)THEN
    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (CNShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE

      CALL MPI_IRECV( CNShapeMapping(iProc)%RecvBuffer   &
                    , CNShapeMapping(iProc)%nRecvShapeElems*4*(PP_N+1)**3   &
                    , MPI_DOUBLE_PRECISION                                  &
                    , iProc                                                 &
                    , 2002                                                  &
                    , MPI_COMM_LEADERS_SHARED                               &
                    , RecvRequestCN(iProc)                                  &
                    , IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (CNShapeMapping(iProc)%nSendShapeElems.EQ.0) CYCLE

      DO iElem=1, CNShapeMapping(iProc)%nSendShapeElems
        CNElemID = GetCNElemID(CNShapeMapping(iProc)%SendShapeElemID(iElem))
        CNShapeMapping(iProc)%SendBuffer(:,:,:,:,iElem) = PartSource(:,:,:,:,CNElemID)
      END DO

      CALL MPI_ISEND( CNShapeMapping(iProc)%SendBuffer   &
                    , CNShapeMapping(iProc)%nSendShapeElems*4*(PP_N+1)**3   &
                    , MPI_DOUBLE_PRECISION                                  &
                    , iProc                                                 &
                    , 2002                                                  &
                    , MPI_COMM_LEADERS_SHARED                               &
                    , SendRequestCN(iProc)                                  &
                    , IERROR)
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE

      IF (CNShapeMapping(iProc)%nRecvShapeElems.NE.0) THEN
        CALL MPI_WAIT(RecvRequestCN(iProc),MPIStatus,IERROR)
        IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF

      IF (CNShapeMapping(iProc)%nSendShapeElems.NE.0) THEN
        CALL MPI_WAIT(SendRequestCN(iProc),MPIStatus,IERROR)
        IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
      END IF
    END DO

    DO iProc = 0,nLeaderGroupProcs-1
      IF (iProc.EQ.myLeaderGroupRank) CYCLE
      IF (CNShapeMapping(iProc)%nRecvShapeElems.EQ.0) CYCLE
      DO iElem=1, CNShapeMapping(iProc)%nRecvShapeElems
        CNElemID = GetCNElemID(CNShapeMapping(iProc)%RecvShapeElemID(iElem))
        PartSource(:,:,:,:,CNElemID) = PartSource(:,:,:,:,CNElemID) + CNShapeMapping(iProc)%RecvBuffer(:,:,:,:,iElem)
      END DO
    END DO
  END IF ! myComputeNodeRank.EQ.0
  CALL MPI_WIN_SYNC(PartSource_Shared_Win, IERROR)
  CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
END IF ! nLeaderGroupProcs.GT.1
#endif
END SUBROUTINE DepositionMethod_SF

END MODULE MOD_PICDepo_Method
