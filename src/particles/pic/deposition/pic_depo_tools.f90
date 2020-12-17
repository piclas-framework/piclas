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

PUBLIC:: DepositParticleOnNodes,CalcCellLocNodeVolumes,ReadTimeAverage,beta
!===================================================================================================================================

CONTAINS


SUBROUTINE DepositParticleOnNodes(iPart,PartPos,GlobalElemID)
!----------------------------------------------------------------------------------------------------------------------------------!
! Deposit the charge of a single particle on the nodes corresponding to the deposition method 'cell_volweight_mean', where the
! charge is stored in NodeSourceExtTmp, which is added to NodeSource in the standard deposition procedure.
! Note that the corresponding volumes are not accounted for yet. The volumes are applied in the deposition routine.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Eval_xyz           ,ONLY: GetPositionInRefElem
USE MOD_Particle_Vars      ,ONLY: PartSpecies,Species
USE MOD_Particle_Vars      ,ONLY: usevMPF,PartMPF,PDM
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_part_tools         ,ONLY: isChargedParticle
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBElemPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars ,ONLY: NodeInfo_Shared
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmpLoc
#else
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExtTmp
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: iPart
REAL,INTENT(IN)                  :: PartPos(1:3)
INTEGER,INTENT(IN)               :: GlobalElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3), Charge
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER                          :: NodeID(1:8)
!===================================================================================================================================

! Sanity checks
IF(.NOT.PDM%ParticleInside(iPart))THEN
  IPWRITE (*,*) "iPart         :", iPart
  IPWRITE (*,*) "global ElemID :", GlobalElemID
  CALL abort(&
      __STAMP__&
      ,'DepositParticleOnNodes(): Particle not inside element.')
ELSEIF(PartSpecies(iPart).LT.0)THEN
  IPWRITE (*,*) "iPart         :", iPart
  IPWRITE (*,*) "global ElemID :", GlobalElemID
  CALL abort(&
      __STAMP__&
      ,'DepositParticleOnNodes(): Negative speciesID')
END IF ! PartSpecies(iPart)

! Skip for neutral particles
IF(isChargedParticle(iPart)) RETURN

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
IF (usevMPF) THEN
  Charge = Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
ELSE
  Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
END IF

CALL GetPositionInRefElem(PartPos,TempPartPos(1:3),GlobalElemID,ForceMode=.TRUE.)

alpha1=0.5*(TempPartPos(1)+1.0)
alpha2=0.5*(TempPartPos(2)+1.0)
alpha3=0.5*(TempPartPos(3)+1.0)

#if USE_MPI
ASSOCIATE( NodeSourceExtTmp => NodeSourceExtTmpLoc )
#endif
  ! Apply charge to nodes (note that the volumes are not accounted for yet here!)
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(GlobalElemID)))
  NodeSourceExtTmp(NodeID(1)) = NodeSourceExtTmp(NodeID(1)) + (Charge*(1-alpha1)*(1-alpha2)*(1-alpha3))
  NodeSourceExtTmp(NodeID(2)) = NodeSourceExtTmp(NodeID(2)) + (Charge*  (alpha1)*(1-alpha2)*(1-alpha3))
  NodeSourceExtTmp(NodeID(3)) = NodeSourceExtTmp(NodeID(3)) + (Charge*  (alpha1)*  (alpha2)*(1-alpha3))
  NodeSourceExtTmp(NodeID(4)) = NodeSourceExtTmp(NodeID(4)) + (Charge*(1-alpha1)*  (alpha2)*(1-alpha3))
  NodeSourceExtTmp(NodeID(5)) = NodeSourceExtTmp(NodeID(5)) + (Charge*(1-alpha1)*(1-alpha2)*  (alpha3))
  NodeSourceExtTmp(NodeID(6)) = NodeSourceExtTmp(NodeID(6)) + (Charge*  (alpha1)*(1-alpha2)*  (alpha3))
  NodeSourceExtTmp(NodeID(7)) = NodeSourceExtTmp(NodeID(7)) + (Charge*  (alpha1)*  (alpha2)*  (alpha3))
  NodeSourceExtTmp(NodeID(8)) = NodeSourceExtTmp(NodeID(8)) + (Charge*(1-alpha1)*  (alpha2)*  (alpha3))
#if USE_MPI
END ASSOCIATE
#endif

#if USE_LOADBALANCE
CALL LBElemPauseTime(GlobalElemID-offsetElem,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to ElemID
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
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems, nComputeNodeProcessors, myComputeNodeRank, MPI_COMM_SHARED
USE MOD_MPI_Shared!        ,ONLY: Allocate_Shared
USE MOD_PICDepo_Vars       ,ONLY: NodeVolume_Shared, NodeVolume_Shared_Win
#else
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif
USE MOD_PICDepo_Vars       ,ONLY: NodeVolume
USE MOD_Particle_Mesh_Vars ,ONLY: ElemsJ, ElemNodeID_Shared, nUniqueGlobalNodes, NodeInfo_Shared
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
INTEGER                          :: j,k,l,iElem, firstElem, lastElem
#if USE_MPI
INTEGER(KIND=MPI_ADDRESS_KIND)   :: MPISharedSize
INTEGER                          :: MessageSize
REAL                             :: NodeVolumeLoc(1:nUniqueGlobalNodes)
#endif
#if USE_DEBUG
INTEGER                          :: I
#endif /*USE_DEBUG*/
INTEGER                          :: NodeID(1:8)
!===================================================================================================================================
#if USE_MPI
MPISharedSize = INT((nUniqueGlobalNodes),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/nUniqueGlobalNodes/),NodeVolume_Shared_Win,NodeVolume_Shared)
CALL MPI_WIN_LOCK_ALL(0,NodeVolume_Shared_Win,IERROR)
NodeVolume => NodeVolume_Shared
firstElem = INT(REAL( myComputeNodeRank   *nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
lastElem  = INT(REAL((myComputeNodeRank+1)*nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
NodeVolumeLoc = 0.
#else
ALLOCATE(NodeVolume(1:nUniqueGlobalNodes))
firstElem = 1
lastElem  = nElems
#endif /*USE_MPI*/

! only CN root nullifies
#if USE_MPI
IF (myComputeNodeRank.EQ.0) THEN
#endif /* USE_MPI*/
  NodeVolume = 0.0
#if USE_MPI
END IF
#endif /* USE_MPI*/

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
MessageSize =  nUniqueGlobalNodes
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(NodeVolumeLoc,NodeVolume,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
ELSE
  CALL MPI_REDUCE(NodeVolumeLoc,0         ,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_SHARED,IERROR)
END IF
CALL MPI_WIN_SYNC(NodeVolume_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

#if USE_DEBUG
! Sanity Check
DO I = 1, nUniqueGlobalNodes
  IF(NodeVolume(I).LE.0.0)THEN
    IPWRITE(UNIT_StdOut,*) "NodeVolume(",I,") =", NodeVolume(I)
    CALL abort(&
        __STAMP__&
        ,'NodeVolume(I) <= 0.0 for I = ',IntInfoOpt=I)
  END IF ! NodeVolume(NodeID(1)).LE.0.0
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
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

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
  CALL ReadAttribute(File_ID,'N',1,IntegerScalar=N_HDF5)
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
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!'
#else
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!'
#endif

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
