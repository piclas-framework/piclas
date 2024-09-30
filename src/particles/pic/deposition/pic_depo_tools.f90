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

MODULE MOD_PICDepo_Tools
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
INTERFACE DepositParticleOnNodes
  MODULE PROCEDURE DepositParticleOnNodes
END INTERFACE

INTERFACE CalcCellLocNodeVolumes
  MODULE PROCEDURE CalcCellLocNodeVolumes
END INTERFACE

INTERFACE ReadTimeAverage
  MODULE PROCEDURE ReadTimeAverage
END INTERFACE

INTERFACE beta
  MODULE PROCEDURE beta
END INTERFACE

INTERFACE DepositPhotonSEEHoles
  MODULE PROCEDURE DepositPhotonSEEHoles
END INTERFACE

PUBLIC:: DepositParticleOnNodes,CalcCellLocNodeVolumes,ReadTimeAverage,beta,DepositPhotonSEEHoles
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Deposit surface charge of positive charges (electron holes) due to SEE from a surface
!> Check if the current iBC is connected to a dielectric region (surface charge currently only for dielectrics)
!===================================================================================================================================
SUBROUTINE DepositPhotonSEEHoles(iBC,NbrOfParticle)
! MODULES
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_PICDepo_Vars           ,ONLY: DoDeposition
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_Particle_Vars          ,ONLY: PEM, PartSpecies, PartState, Species, usevMPF, PartMPF
USE MOD_Part_Tools             ,ONLY: GetNextFreePosition
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iBC           !< BC of emitted particle (only if defined and >0)
INTEGER,INTENT(IN) :: NbrOfParticle !< Number of newly inserted electrons
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: iPart          !< i-th inserted electron
INTEGER  :: ParticleIndex  !< index of i-th inserted electron in particle array
REAL     :: MPF            !< macro-particle factor
REAL     :: ChargeHole     !< Charge of SEE electrons holes
!===================================================================================================================================
! Only continue when Species(i)%Init(iInit)%PartBCIndex.GT.0
IF(iBC.LE.0) RETURN
! Check if deposition and surface charge is active and if the current BC is a charged BC
IF(DoDeposition.AND.DoDielectricSurfaceCharge.AND.PartBound%Dielectric(iBC))THEN
  DO iPart = 1, NbrOfParticle
    ! Get index from next free position array
    ParticleIndex = GetNextFreePosition(iPart)

    ! Get charge
    IF(usevMPF)THEN
      MPF = PartMPF(ParticleIndex)
    ELSE
      MPF = Species(PartSpecies(ParticleIndex))%MacroParticleFactor
    END IF ! usevMPF
    ChargeHole = -Species(PartSpecies(ParticleIndex))%ChargeIC*MPF

    ! Create electron hole (i.e. positive surface charge)
    CALL DepositParticleOnNodes(ChargeHole, PartState(1:3,ParticleIndex), PEM%GlobalElemID(ParticleIndex))
  END DO
END IF

END SUBROUTINE DepositPhotonSEEHoles

!===================================================================================================================================
!> Deposit the charge of a single particle on the nodes corresponding to the deposition method 'cell_volweight_mean', where the
!> charge is stored in NodeSourceExtTmp, which is added to NodeSource in the standard deposition procedure.
!> Note that the corresponding volumes are not accounted for yet. The volumes are applied in the deposition routine.
!===================================================================================================================================
SUBROUTINE DepositParticleOnNodes(Charge,PartPos,GlobalElemID)
! MODULES
USE MOD_Globals
USE MOD_Globals            ,ONLY: VECNORM,ElementOnProc
USE MOD_Globals_Vars       ,ONLY: ElementaryCharge
USE MOD_Eval_xyz           ,ONLY: GetPositionInRefElem
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeCoords_Shared
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBElemPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars ,ONLY: NodeInfo_Shared
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmp
#else
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: Charge        !< Charge that is deposited on nodes
REAL,INTENT(IN)                  :: PartPos(1:3)
INTEGER,INTENT(IN)               :: GlobalElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3)
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER                          :: NodeID(1:8),iNode
LOGICAL                          :: SucRefPos
REAL                             :: norm,PartDistDepo(8),DistSum
!===================================================================================================================================

! Skip for neutral particles and reflected particles or species swapped particles where impacting and reflecting particle carry the
! same charge. Deposit only particles that are deleted on the surface or change their charge on contact (e.g. neutralization)
IF(ABS(Charge).LE.0.0) RETURN

#if USE_LOADBALANCE
! Only measure time if particle is deposited on local proc
IF(ElementOnProc(GlobalElemID)) CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

CALL GetPositionInRefElem(PartPos, TempPartPos(1:3), GlobalElemID, ForceMode = .TRUE., isSuccessful = SucRefPos)

#if USE_MPI
ASSOCIATE( NodeSourceExt => NodeSourceExtTmp )
#endif
  ! Check if GetPositionInRefElem was able to find the reference position (via ref. mapping), else use distance-based deposition
  IF(SucRefPos)THEN
    alpha1=0.5*(TempPartPos(1)+1.0)
    alpha2=0.5*(TempPartPos(2)+1.0)
    alpha3=0.5*(TempPartPos(3)+1.0)

    ! Apply charge to nodes (note that the volumes are not accounted for yet here!)
    NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(GlobalElemID)))
    NodeSourceExt(NodeID(1)) = NodeSourceExt(NodeID(1)) + (Charge*(1-alpha1)*(1-alpha2)*(1-alpha3))
    NodeSourceExt(NodeID(2)) = NodeSourceExt(NodeID(2)) + (Charge*  (alpha1)*(1-alpha2)*(1-alpha3))
    NodeSourceExt(NodeID(3)) = NodeSourceExt(NodeID(3)) + (Charge*  (alpha1)*  (alpha2)*(1-alpha3))
    NodeSourceExt(NodeID(4)) = NodeSourceExt(NodeID(4)) + (Charge*(1-alpha1)*  (alpha2)*(1-alpha3))
    NodeSourceExt(NodeID(5)) = NodeSourceExt(NodeID(5)) + (Charge*(1-alpha1)*(1-alpha2)*  (alpha3))
    NodeSourceExt(NodeID(6)) = NodeSourceExt(NodeID(6)) + (Charge*  (alpha1)*(1-alpha2)*  (alpha3))
    NodeSourceExt(NodeID(7)) = NodeSourceExt(NodeID(7)) + (Charge*  (alpha1)*  (alpha2)*  (alpha3))
    NodeSourceExt(NodeID(8)) = NodeSourceExt(NodeID(8)) + (Charge*(1-alpha1)*  (alpha2)*  (alpha3))
  ELSE
     NodeID = ElemNodeID_Shared(:,GetCNElemID(GlobalElemID))
     DO iNode = 1, 8
       norm = VECNORM(NodeCoords_Shared(1:3, NodeID(iNode)) - PartPos(1:3))
       IF(norm.GT.0.)THEN
         PartDistDepo(iNode) = 1./norm
       ELSE
         PartDistDepo(:) = 0.
         PartDistDepo(iNode) = 1.0
         EXIT
       END IF ! norm.GT.0.
     END DO
     DistSum = SUM(PartDistDepo(1:8))
     DO iNode = 1, 8
       NodeSourceExt(NodeInfo_Shared(NodeID(iNode))) = NodeSourceExt(NodeInfo_Shared(NodeID(iNode)))  &
         +  PartDistDepo(iNode)/DistSum*Charge
     END DO
  END IF ! SucRefPos
#if USE_MPI
END ASSOCIATE
#endif

#if USE_LOADBALANCE
! Only measure time if particle is deposited on local proc
IF(ElementOnProc(GlobalElemID)) CALL LBElemPauseTime(GlobalElemID-offsetElem,tLBStart)
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositParticleOnNodes


SUBROUTINE CalcCellLocNodeVolumes()
!===================================================================================================================================
!> Initialize sub-cell volumes around nodes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis              ,ONLY: InitializeVandermonde
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars ,ONLY: wGP, xGP, wBary
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems, nComputeNodeProcessors, myComputeNodeRank, MPI_COMM_SHARED
USE MOD_PICDepo_Vars       ,ONLY: NodeVolume_Shared, NodeVolume_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif
USE MOD_PICDepo_Vars       ,ONLY: NodeVolume,Periodic_nNodes,Periodic_offsetNode,Periodic_Nodes
USE MOD_Particle_Mesh_Vars ,ONLY: ElemsJ, ElemNodeID_Shared, nUniqueGlobalNodes, NodeInfo_Shared,GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: Vdm_loc(0:1,0:PP_N),wGP_loc,xGP_loc(0:1),DetJac(1,0:1,0:1,0:1)
REAL                             :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N)
INTEGER                          :: j,k,l,iElem, firstElem, lastElem, iNode, jNode
REAL                             :: NodeVolumeLoc(1:nUniqueGlobalNodes)
#if USE_MPI
INTEGER                          :: MessageSize
#if USE_DEBUG
INTEGER                          :: I
#endif /*USE_DEBUG*/
#endif
INTEGER                          :: NodeID(1:8)
!===================================================================================================================================
NodeVolumeLoc = 0.
#if USE_MPI
CALL Allocate_Shared((/nUniqueGlobalNodes/),NodeVolume_Shared_Win,NodeVolume_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeVolume_Shared_Win,IERROR)
NodeVolume => NodeVolume_Shared
firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
! only CN root nullifies
IF (myComputeNodeRank.EQ.0) NodeVolume = 0.0
! This sync/barrier is required as it cannot be guaranteed that the zeros have been written to memory by the time the MPI_REDUCE
! is executed (see MPI specification). Until the Sync is complete, the status is undefined, i.e., old or new value or utter nonsense.
CALL BARRIER_AND_SYNC(NodeVolume_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(NodeVolume(1:nUniqueGlobalNodes))
NodeVolume = 0.0
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

IF (PP_N.NE.1) THEN
  xGP_loc(0) = -0.5
  xGP_loc(1) = 0.5
  wGP_loc = 1.
  CALL InitializeVandermonde(PP_N,1,wBary,xGP,xGP_loc, Vdm_loc)
END IF
! ElemNodeID and ElemsJ use compute node elems
DO iElem = firstElem, lastElem
  IF (PP_N.EQ.1) THEN
    wGP_loc = wGP(0)
    DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
      DetJac(1,j,k,l)=1./ElemsJ(j,k,l,iElem)
    END DO; END DO; END DO
  ELSE
    DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
      DetLocal(1,j,k,l)=1./ElemsJ(j,k,l,iElem)
    END DO; END DO; END DO
    CALL ChangeBasis3D(1,PP_N, 1, Vdm_loc, DetLocal(:,:,:,:),DetJac(:,:,:,:))
  END IF
#if USE_MPI
  ASSOCIATE( NodeVolume => NodeVolumeLoc )
#endif /*USE_MPI*/
    ! Get UniqueNodeIDs
    NodeID = NodeInfo_Shared(ElemNodeID_Shared(1:8,iElem))
    NodeVolume(NodeID(1)) = NodeVolume(NodeID(1)) + DetJac(1,0,0,0)
    NodeVolume(NodeID(2)) = NodeVolume(NodeID(2)) + DetJac(1,1,0,0)
    NodeVolume(NodeID(3)) = NodeVolume(NodeID(3)) + DetJac(1,1,1,0)
    NodeVolume(NodeID(4)) = NodeVolume(NodeID(4)) + DetJac(1,0,1,0)
    NodeVolume(NodeID(5)) = NodeVolume(NodeID(5)) + DetJac(1,0,0,1)
    NodeVolume(NodeID(6)) = NodeVolume(NodeID(6)) + DetJac(1,1,0,1)
    NodeVolume(NodeID(7)) = NodeVolume(NodeID(7)) + DetJac(1,1,1,1)
    NodeVolume(NodeID(8)) = NodeVolume(NodeID(8)) + DetJac(1,0,1,1)
#if USE_MPI
  END ASSOCIATE
#endif /*USE_MPI*/
END DO

#if USE_MPI
! collect the information from the proc-local shadow arrays in the compute-node shared array
MessageSize = nUniqueGlobalNodes
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(NodeVolumeLoc , NodeVolume , MessageSize , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
ELSE
  CALL MPI_REDUCE(NodeVolumeLoc , 0          , MessageSize , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_SHARED , IERROR)
END IF
CALL BARRIER_AND_SYNC(NodeVolume_Shared_Win,MPI_COMM_SHARED)
#endif
IF (GEO%nPeriodicVectors.GT.0) THEN
#if USE_MPI
  IF (myComputeNodeRank.EQ.0) THEN
#endif /*USE_MPI*/
    ! Root node acquires the periodic contribution
    NodeVolumeLoc = 0.
    DO iNode = 1,nUniqueGlobalNodes
      IF (Periodic_nNodes(iNode).GT.0) THEN
        DO jNode = Periodic_offsetNode(iNode)+1,Periodic_offsetNode(iNode)+Periodic_nNodes(iNode)
          NodeVolumeLoc(iNode) = NodeVolumeLoc(iNode) + NodeVolume(Periodic_Nodes(jNode))
        END DO ! jNode
      END IF ! Periodic_nNodes(iNode).GT.0
    END DO ! iNode = nUniqueGlobalNodes
    ! Only node root adds periodic contribution to shared array
    NodeVolume = NodeVolume + NodeVolumeLoc
#if USE_MPI
  END IF ! myComputeNodeRank.EQ.0
  CALL BARRIER_AND_SYNC(NodeVolume_Shared_Win,MPI_COMM_SHARED)
#endif
END IF

#if USE_MPI
#if USE_DEBUG
! Sanity Check: Only check UniqueGlobalNodes that are on the compute node (total)
DO iElem = firstElem, lastElem
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(1:8,iElem))
  DO I = 1, 8
    IF(NodeVolume(NodeID(I)).LE.0.0)THEN
      IPWRITE(UNIT_StdOut,'(I0,A,I0,A,ES25.17E3)') " NodeVolume(NodeID(",I,")) =", NodeVolume(NodeID(I))
      CALL abort(__STAMP__,'NodeVolume(NodeID(I)) <= 0.0 for NodeID(I) = ',IntInfoOpt=NodeID(I))
    END IF ! NodeVolume(NodeID(1)).LE.0.0
  END DO
END DO ! I = 1, nUniqueGlobalNodes
#endif /*USE_DEBUG*/
#endif /*USE_MPI*/

END SUBROUTINE CalcCellLocNodeVolumes


SUBROUTINE ReadTimeAverage(FileName)
!===================================================================================================================================
! Read in ChargeDensity and save to PartSource
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Input,              ONLY:ReadArray,ReadAttribute,File_ID,OpenDataFile,CloseDataFile,DatasetExists
USE MOD_Particle_Vars,           ONLY:nSpecies
USE MOD_PICDepo_Vars,            ONLY:PartSource
USE MOD_Mesh_Vars,               ONLY:OffsetElem,nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE         :: U(:,:,:,:,:)
INTEGER                  :: iSpec, iElem, kk, ll, mm
INTEGER                  :: Rank
INTEGER                  :: nVars, iVar, N_HDF5
INTEGER,ALLOCATABLE      :: PartSourceToVar(:)
INTEGER(HID_T)           :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
LOGICAL                  :: SolutionExists
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
CHARACTER(LEN=10)        :: strhelp
#if USE_MPI
REAL                     :: StartT,EndT
#endif /*USE_MPI*/
!===================================================================================================================================

SWRITE(UNIT_stdOut,*)'Using TimeAverage as constant PartSource(4) from file:',TRIM(FileName)
#if USE_MPI
StartT=MPI_WTIME()
#endif

IF(MPIRoot) THEN
  IF(.NOT.FILEEXISTS(FileName))  CALL abort(__STAMP__, &
        'TimeAverage-File "'//TRIM(FileName)//'" does not exist',999,999.)
END IF
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

! get attributes
CALL DatasetExists(File_ID,'DG_Solution',SolutionExists)
IF(MPIRoot)THEN
  IF(.NOT.SolutionExists)  CALL abort(&
    __STAMP__&
    ,'DG_Solution in TimeAverage-File "'//TRIM(FileName)//'" does not exist!')
END IF
CALL H5DOPEN_F(File_ID, 'DG_Solution', Dset_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)
IF(MPIRoot)THEN
  IF(INT(Dims(Rank),4).NE.nGlobalElems)  CALL abort(&
    __STAMP__&
    ,' MeshSize and Size of TimeAverage-File "'//TRIM(FileName)//'" does not match!')
END IF
nVars=INT(Dims(1),4)
ALLOCATE(VarNames(nVars))
CALL ReadAttribute(File_ID,'VarNames',nVars,StrArray=VarNames)

ALLOCATE(PartSourceToVar(nSpecies))
PartSourceToVar=0
DO iSpec=1,nSpecies
  WRITE(strhelp,'(I2.2)') iSpec
  DO iVar=1,nVars
    IF (VarNames(iVar).EQ.TRIM('ChargeDensity-Spec')//TRIM(strhelp)) THEN
      PartSourceToVar(iSpec)=iVar
      EXIT
    END IF
  END DO
END DO
IF (.NOT.ANY(PartSourceToVar.NE.0)) CALL abort(__STAMP__, &
  'No PartSource found in TimeAverage-File "'//TRIM(FileName)//'"!!!',999,999.)
DEALLOCATE(VarNames)

!-- read state
ALLOCATE(U(nVars,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
CALL ReadAttribute(File_ID,'N',1,IntScalar=N_HDF5)
IF(N_HDF5.EQ.PP_N)THEN! No interpolation needed, read solution directly from file
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVars       => INT(nVars,IK)     ,&
        PP_nElems   => INT(PP_nElems,IK) ,&
        OffsetElem  => INT(OffsetElem,IK) )
        CALL ReadArray('DG_Solution',5,(/nVars,INT(PP_N,IK)+1_IK,INT(PP_N,IK)+1_IK,INT(PP_N,IK)+1_IK,PP_nElems/),OffsetElem,5,RealArray=U)
  END ASSOCIATE
ELSE
  CALL abort(__STAMP__, &
        'N_HDF5.NE.PP_N !',999,999.)
END IF
CALL CloseDataFile()

!-- save to PartSource
PartSource(4,:,:,:,:)=0.
DO iSpec=1,nSpecies
  IF (PartSourceToVar(iSpec).NE.0) THEN
    DO iElem=1,PP_nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
            PartSource(4,mm,ll,kk,iElem)=PartSource(4,mm,ll,kk,iElem)+U(PartSourceToVar(iSpec),mm,ll,kk,iElem)
          END DO
        END DO
      END DO
    END DO
  END IF
END DO
DEALLOCATE(U)
DEALLOCATE(PartSourceToVar)

#if USE_MPI
EndT=MPI_WTIME()
SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' Readin took  [',EndT-StartT,'s].'
#endif
SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' DONE!'

END SUBROUTINE ReadTimeAverage


FUNCTION beta(z,w)
!============================================================================================================================
! calculates the beta function
!============================================================================================================================
! use MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                             :: w, z
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                          :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)
END FUNCTION beta


END MODULE MOD_PICDepo_Tools
