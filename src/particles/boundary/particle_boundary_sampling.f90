#include "boltzplatz.h"

MODULE MOD_Particle_Boundary_Sampling
!===================================================================================================================================
!! Determines how particles interact with a given boundary condition 
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
INTERFACE InitParticleBoundarySampling
  MODULE PROCEDURE InitParticleBoundarySampling
END INTERFACE

INTERFACE FinalizeParticleBoundarySampling
  MODULE PROCEDURE FinalizeParticleBoundarySampling
END INTERFACE
 
INTERFACE WriteSurfSampleToHDF5
  MODULE PROCEDURE WriteSurfSampleToHDF5
END INTERFACE

#ifdef MPI
INTERFACE ExchangeSurfData
  MODULE PROCEDURE ExchangeSurfData
END INTERFACE
#endif /*MPI*/

PUBLIC::InitParticleBoundarySampling
PUBLIC::WriteSurfSampleToHDF5
PUBLIC::FinalizeParticleBoundarySampling
#ifdef MPI
PUBLIC::ExchangeSurfData
#endif /*MPI*/
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleBoundarySampling() 
!===================================================================================================================================
! init of particle boundary sampling
! default: use for sampling same polynomial degree as NGeo
! 1) mark sides for sampling
! 2) build special MPI communicator
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY:NGeo,BC,nSides,nBCSides,nBCs,BoundaryName
USE MOD_ReadInTools             ,ONLY:GETINT
USE MOD_Particle_Boundary_Vars  ,ONLY:nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample,SurfMesh,SampWall,nSurfBC,SurfBCName
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides
USE MOD_Particle_Vars           ,ONLY:nSpecies
USE MOD_DSMC_Vars               ,ONLY:useDSMC,DSMC
USE MOD_Basis                   ,ONLY:LegendreGaussNodesAndWeights
USE MOD_Particle_Surfaces       ,ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY:BezierControlPoints3D,BezierSampleN
USE MOD_Particle_Mesh_Vars      ,ONLY:PartBCSideList
USE MOD_Particle_Tracking_Vars  ,ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI_Vars       ,ONLY:PartMPI
#else
USE MOD_Particle_Boundary_Vars  ,ONLY:offSetSurfSide,SurfCOMM
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q,iSide,SurfSideID,SideID
INTEGER                                :: iSample,jSample, iBC
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,area,tmpI2,tmpJ2
CHARACTER(2)                           :: hilf
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
!===================================================================================================================================
 
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'
WRITE(UNIT=hilf,FMT='(I2)') NGeo
nSurfSample = GETINT('DSMC-nSurfSample',hilf)
! IF (NGeo.GT.nSurfSample) THEN
!   nSurfSample = NGeo
! END IF
IF (useDSMC) THEN
  IF (DSMC%WallModel.GT.0) THEN
    IF (nSurfSample.NE.BezierSampleN) THEN
!       nSurfSample = BezierSampleN
    CALL abort(&
__STAMP__&
,'Error: nSurfSample not equal to BezierSampleN. Problem for Desorption + Surfflux')
    END IF
  END IF
END IF
 
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))

dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1. 
END DO

! create boundary name mapping for surfaces SurfaceBC number mapping
nSurfBC = 0
ALLOCATE(BCName(1:nBCs))
DO iBC=1,nBCs
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(iBC)).EQ.PartBound%ReflectiveBC) THEN
  nSurfBC = nSurfBC + 1
  BCName(nSurfBC) = BoundaryName(iBC)
  END IF
END DO
IF (nSurfBC.GE.1) THEN
ALLOCATE(SurfBCName(1:nSurfBC))
  DO iBC=1,nSurfBC
    SurfBCName(iBC) = BCName(iBC)
  END DO
END IF
DEALLOCATE(BCName)

! get number of BC-Sides
ALLOCATE(SurfMesh%SideIDToSurfID(1:nTotalSides))
SurfMesh%SideIDToSurfID(1:nTotalSides)=-1
! first own sides
SurfMesh%nSides=0
DO iSide=1,nBCSides
  IF(BC(iSide).EQ.0) CYCLE
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
    SurfMesh%nSides = SurfMesh%nSides + 1
    SurfMesh%SideIDToSurfID(iSide)=SurfMesh%nSides
    !SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nSides
  END IF
END DO

! halo sides
SurfMesh%nTotalSides=SurfMesh%nSides
DO iSide=nSides+1,nTotalSides
  IF(BC(iSide).EQ.0) CYCLE
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
    SurfMesh%nTotalSides = SurfMesh%nTotalSides + 1
    SurfMesh%SideIDToSurfID(iSide)=SurfMesh%nTotalSides
    !SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
  END IF
END DO

SurfMesh%SurfOnProc=.FALSE.
!IF(SurfMesh%nSides.GT.0) 
IF(SurfMesh%nTotalSides.GT.0) SurfMesh%SurfOnProc=.TRUE.

#ifdef MPI
CALL MPI_ALLREDUCE(SurfMesh%nSides,SurfMesh%nGlobalSides,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError)
#else
SurfMesh%nGlobalSides=SurfMesh%nSides
#endif


SWRITE(UNIT_stdOut,'(A,I8)') ' nGlobalSurfSides ', SurfMesh%nGlobalSides

SurfMesh%SampSize=9+3+nSpecies ! Energy + Force + nSpecies
#ifdef MPI
! split communitator
CALL InitSurfCommunicator()
IF(SurfMesh%SurfOnProc) THEN
  CALL GetHaloSurfMapping()
END IF
#else
SurfCOMM%MyRank=0
SurfCOMM%MPIRoot=.TRUE.
SurfCOMM%nProcs=1
SurfCOMM%MyOutputRank=0
SurfCOMM%nOutputProcs = 1
SurfCOMM%nOutputProcs=1
! get correct offsets
OffSetSurfSide=0
#endif /*MPI*/

IF(.NOT.SurfMesh%SurfOnProc) RETURN


! allocate everything
ALLOCATE(SampWall(1:SurfMesh%nTotalSides))


DO iSide=1,SurfMesh%nTotalSides ! caution: iSurfSideID
  ALLOCATE(SampWall(iSide)%State(1:SurfMesh%SampSize,1:nSurfSample,1:nSurfSample))
  SampWall(iSide)%State=0.
  !ALLOCATE(SampWall(iSide)%Energy(1:9,0:nSurfSample,0:nSurfSample)         &
  !        ,SampWall(iSide)%Force(1:9,0:nSurfSample,0:nSurfSample)          &
  !        ,SampWall(iSide)%Counter(1:nSpecies,0:nSurfSample,0:nSurfSample) )
  ! nullify
  !SampWall(iSide)%Energy(1:9,0:nSurfSample,0:nSurfSample)         = 0.
  !SampWall(iSide)%Force(1:9,0:nSurfSample,0:nSurfSample)          = 0.
  !SampWall(iSide)%Counter(1:nSpecies,0:nSurfSample,0:nSurfSample) = 0.
END DO

ALLOCATE(SurfMesh%SurfaceArea(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides)) 
SurfMesh%SurfaceArea=0.


ALLOCATE(Xi_NGeo( 0:NGeo)  &
        ,wGP_NGeo(0:NGeo) )
CALL LegendreGaussNodesAndWeights(NGeo,Xi_NGeo,wGP_NGeo)

! compute area of sub-faces
tmp1=dXiEQ_SurfSample/2.0 !(b-a)/2
DO iSide=1,nTotalSides
  SurfSideID=SurfMesh%SideIDToSurfID(iSide)
  IF(SurfSideID.EQ.-1) CYCLE
  IF(DoRefMapping)THEN
    SideID=PartBCSideList(iSide)
  ELSE
    SideID=iSide
  END IF
  ! call here stephens algorithm to compute area 
  DO jSample=1,nSurfSample
    DO iSample=1,nSurfSample
      area=0.
      tmpI2=(XiEQ_SurfSample(iSample-1)+XiEQ_SurfSample(iSample))/2. ! (a+b)/2
      tmpJ2=(XiEQ_SurfSample(jSample-1)+XiEQ_SurfSample(jSample))/2. ! (a+b)/2
      DO q=0,NGeo
        DO p=0,NGeo
          XiOut(1)=tmp1*Xi_NGeo(p)+tmpI2
          XiOut(2)=tmp1*Xi_NGeo(q)+tmpJ2
          CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
                                                  ,Gradient=gradXiEta3D)
          ! calculate first fundamental form
          E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
          F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
          G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
          D=SQRT(E*G-F*F)
          area = area+tmp1*tmp1*D*wGP_NGeo(p)*wGP_NGeo(q)      
        END DO
      END DO
      SurfMesh%SurfaceArea(iSample,jSample,SurfSideID) = area 
    END DO ! iSample=1,nSurfSample
  END DO ! jSample=1,nSurfSample
END DO ! iSide=1,nTotalSides

DEALLOCATE(Xi_NGeo,wGP_NGeo)

SWRITE(UNIT_stdOut,'(A)') ' ... DONE.'

END SUBROUTINE InitParticleBoundarySampling

#ifdef MPI
SUBROUTINE InitSurfCommunicator()
!===================================================================================================================================
! Creates two new subcommunicators. 
! SurfCOMM%COMM contains all MPI-Ranks which have reflective boundary faces in their halo-region and process which have reflective
! boundary faces in their origin region. This communicator is used to communicate the wall-sampled values of halo-faces to the
! origin face
! SurfCOMM%OutputCOMM is another subset. This communicator contains only the processes with origin surfaces. It is used to perform
! collective writes of the surf-sampled values.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars   ,ONLY:SurfMesh
USE MOD_Particle_Boundary_Vars   ,ONLY:SurfCOMM
USE MOD_Particle_MPI_Vars        ,ONLY:PartMPI
USE MOD_Particle_Boundary_Vars   ,ONLY:OffSetSurfSideMPI,OffSetSurfSide
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color,iProc
INTEGER                   :: noSurfrank,Surfrank
LOGICAL                   :: hasSurf
INTEGER,ALLOCATABLE       :: countSurfSideMPI(:)
LOGICAL                   :: OutputOnProc
!===================================================================================================================================
color=MPI_UNDEFINED
IF(SurfMesh%SurfonProc) color=2
! THEN
! color=2
! ELSE
! color=1
! END IF
! create ranks for RP communicator
IF(PartMPI%MPIRoot) THEN
  Surfrank=-1
  noSurfrank=-1
  SurfCOMM%Myrank=0
  IF(SurfMesh%SurfonProc) THEN
    Surfrank=0
  ELSE 
    noSurfrank=0
  END IF
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    IF(hasSurf) THEN
      SurfRank=SurfRank+1
      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    ELSE
      noSurfRank=noSurfRank+1
      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    END IF
  END DO
ELSE
  CALL MPI_SEND(SurfMesh%SurfOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
  CALL MPI_RECV(SurfCOMM%MyRank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
END IF

! create new SurfMesh communicator for SurfMesh communiction
CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyRank, SurfCOMM%COMM,iError)
IF(SurfMesh%SurfOnPRoc) THEN
  CALL MPI_COMM_SIZE(SurfCOMM%COMM, SurfCOMM%nProcs,iError)
ELSE
  SurfCOMM%nProcs = 0
END IF
SurfCOMM%MPIRoot=.FALSE.
IF(SurfCOMM%MyRank.EQ.0 .AND. SurfMesh%SurfOnProc) THEN
  SurfCOMM%MPIRoot=.TRUE.
  WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF COMM:        ',SurfCOMM%nProcs,' procs'
END IF

! now, create output communicator
OutputOnProc=.FALSE.
color=MPI_UNDEFINED
IF(SurfMesh%nSides.GT.0) THEN
  OutputOnProc=.TRUE.
  color=4
END IF

IF(PartMPI%MPIRoot) THEN
  Surfrank=-1
  noSurfrank=-1
  SurfCOMM%MyOutputRank=0
  IF(SurfMesh%nSides.GT.0) THEN
    Surfrank=0
  ELSE 
    noSurfrank=0
  END IF
  DO iProc=1,nProcessors-1
    CALL MPI_RECV(hasSurf,1,MPI_LOGICAL,iProc,0,MPI_COMM_WORLD,MPIstatus,iError)
    IF(hasSurf) THEN
      SurfRank=SurfRank+1
      CALL MPI_SEND(SurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    ELSE
      noSurfRank=noSurfRank+1
      CALL MPI_SEND(noSurfRank,1,MPI_INTEGER,iProc,0,MPI_COMM_WORLD,iError)
    END IF
  END DO
ELSE
  CALL MPI_SEND(OutputOnProc,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,iError)
  CALL MPI_RECV(SurfCOMM%MyOutputRank,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,MPIstatus,iError)
END IF

! create new SurfMesh Output-communicator 
CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyOutputRank, SurfCOMM%OutputCOMM,iError)
IF(OutputOnPRoc)THEN
  CALL MPI_COMM_SIZE(SurfCOMM%OutputCOMM, SurfCOMM%nOutputProcs,iError)
ELSE
  SurfCOMM%nOutputProcs = 0
END IF
SurfCOMM%MPIOutputRoot=.FALSE.
IF(SurfCOMM%MyOutputRank.EQ.0 .AND. OutputOnProc) THEN
  SurfCOMM%MPIOutputRoot=.TRUE.
  WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF OUTPUT-COMM: ',SurfCOMM%nOutputProcs,' procs'
END IF

IF(SurfMesh%nSides.EQ.0) RETURN

! get correct offsets
ALLOCATE(offsetSurfSideMPI(0:SurfCOMM%nOutputProcs))
offsetSurfSideMPI=0
ALLOCATE(countSurfSideMPI(0:SurfCOMM%nOutputProcs-1))
countSurfSideMPI=0

CALL MPI_GATHER(SurfMesh%nSides,1,MPI_INTEGER,countSurfSideMPI,1,MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)

IF (SurfCOMM%MPIOutputRoot) THEN
  DO iProc=1,SurfCOMM%nOutputProcs-1
    offsetSurfSideMPI(iProc)=SUM(countSurfSideMPI(0:iProc-1))
  END DO
  offsetSurfSideMPI(SurfCOMM%nOutputProcs)=SUM(countSurfSideMPI(:))
END IF
CALL MPI_BCAST (offsetSurfSideMPI,size(offsetSurfSideMPI),MPI_INTEGER,0,SurfCOMM%OutputCOMM,iError)
offsetSurfSide=offsetSurfSideMPI(SurfCOMM%MyOutputRank)

END SUBROUTINE InitSurfCommunicator


SUBROUTINE GetHaloSurfMapping() 
!===================================================================================================================================
! build all missing stuff for surface-sampling communicator, like
! offSetMPI
! MPI-neighbor list
! PartHaloSideToProc
! only receiving process knows to which local side the sending information is going, the sending process does not know the final 
! sideid 
! if only a processes has to send his data to another, the whole structure is build, but not filled. only if the nsidesrecv,send>0
! communication is performed
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars                   ,ONLY:nSides,nBCSides
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars           ,ONLY:PartHaloSideToProc,PartHaloElemToProc,SurfSendBuf,SurfRecvBuf,SurfExchange
USE MOD_Particle_Mesh_Vars          ,ONLY:nTotalSides,PartSideToElem,PartElemToSide
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: isMPINeighbor(0:SurfCOMM%nProcs-1)
LOGICAL                           :: RecvMPINeighbor(0:SurfCOMM%nProcs-1)
INTEGER                           :: nDOF,ALLOCSTAT,SideID
INTEGER                           :: iProc, GlobalProcID,iSide,ElemID,SurfSideID,LocalProcID,iSendSide,iRecvSide,iPos
INTEGER,ALLOCATABLE               :: recv_status_list(:,:) 
INTEGER                           :: RecvRequest(0:SurfCOMM%nProcs-1),SendRequest(0:SurfCOMM%nProcs-1)
INTEGER                           :: SurfToGlobal(0:SurfCOMM%nProcs-1)
!===================================================================================================================================

nDOF=(nSurfSample)*(nSurfSample)


! get mapping from local rank to global rank...
CALL MPI_ALLGATHER(MyRank,1, MPI_INTEGER, SurfToGlobal(0:SurfCOMM%nProcs-1), 1, MPI_INTEGER, SurfCOMM%COMM, IERROR)

! get list of mpi surf-comm neighbors in halo-region
isMPINeighbor=.FALSE.
IF(SurfMesh%nTotalSides.GT.SurfMesh%nSides)THEN
  ! PartHaloSideToProc has the mapping from the local sideid to the corresponding 
  ! element-id, localsideid, native-proc-id and local-proc-id
  ! caution: 
  ! native-proc-id is id in global list || PartMPI%COMM
  ! local-proc-id is the nth neighbor in SurfCOMM%COMM
  ! caution2:
  !  this mapping is done only for reflective bcs, thus, each open side or non-reflective side 
  !  points to -1
  ! get all MPI-neighbors to communicate with
  DO iProc=0,SurfCOMM%nProcs-1
    IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
    GlobalProcID=SurfToGlobal(iProc)
    DO iSide=nSides+1,nTotalSides
      SurfSideID=SurfMesh%SideIDToSurfID(iSide)
      IF(SurfSideID.EQ.-1) CYCLE
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(ElemID.LE.PP_nElems)THEN
        IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
      END IF
      IF(ElemID.LE.1)THEN
        IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
      END IF
      IF(GlobalProcID.EQ.PartHaloElemToProc(NATIVE_PROC_ID,ElemID))THEN
        IF(.NOT.isMPINeighbor(iProc))THEN
          isMPINeighbor(iProc)=.TRUE.
        END IF
      END IF
    END DO ! iSide=nSides+1,nTotalSides
  END DO ! iProc = 0, SurfCOMM%nProcs-1
END IF

! Make sure SurfCOMM%MPINeighbor is consistent
! 1) communication of found is neighbor
ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,0:SurfCOMM%nProcs-1))
! receive and send connection information
DO iProc=0,SurfCOMM%nProcs-1
  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
  CALL MPI_IRECV( RecvMPINeighbor(iProc)                    &
                , 1                                         &
                , MPI_LOGICAL                               &
                , iProc                                     &
                , 1001                                      &
                , SurfCOMM%COMM                             &
                , RecvRequest(iProc)                        &  
                , IERROR )
END DO ! iProc
DO iProc=0,SurfCOMM%nProcs-1
  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
  CALL MPI_ISEND( isMPINeighbor(iProc)                      &
                , 1                                         &
                , MPI_LOGICAL                               &
                , iProc                                     &
                , 1001                                      &
                , SurfCOMM%COMM                             &
                , SendRequest(iProc)                        & 
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! finish communication
DO iProc=1,SurfCOMM%nProcs-1
  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
  CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc
DEALLOCATE(RECV_STATUS_LIST)

! 2) finalize MPI consistency check.
DO iProc=0,SurfCOMM%nProcs-1
  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
  IF (RecvMPINeighbor(iProc).NEQV.isMPINeighbor(iProc)) THEN
    IF(.NOT.isMPINeighbor(iProc))THEN
      isMPINeighbor(iProc)=.TRUE.
      ! it is a debug output
      !IPWRITE(UNIT_stdOut,*) ' Found missing mpi-neighbor'
    END IF
  END IF
END DO ! iProc

! build the mapping
SurfCOMM%nMPINeighbors=0
IF(ANY(IsMPINeighbor))THEN
  ! PartHaloSideToProc has the mapping from the local sideid to the corresponding 
  ! element-id, localsideid, native-proc-id and local-proc-id
  ! caution: 
  ! native-proc-id is id in global list || PartMPI%COMM
  ! local-proc-id is the nth neighbor in SurfCOMM%COMM
  ! caution2:
  !  this mapping is done only for reflective bcs, thus, each open side or non-reflective side 
  !  points to -1
  IF(SurfMesh%nTotalSides.GT.SurfMesh%nSides)THEN
    ALLOCATE(PartHaloSideToProc(1:4,nSides+1:nTotalSides))
    PartHaloSideToProc=-1
    ! get all MPI-neighbors to communicate with
    DO iProc=0,SurfCOMM%nProcs-1
      IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
      IF(isMPINeighbor(iProc))THEN
        SurfCOMM%nMPINeighbors=SurfCOMM%nMPINeighbors+1
      ELSE
        CYCLE
      END IF
      GlobalProcID=SurfToGlobal(iProc)
      DO iSide=nSides+1,nTotalSides
        SurfSideID=SurfMesh%SideIDToSurfID(iSide)
        IF(SurfSideID.EQ.-1) CYCLE
        ! get elemid
        ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
        IF(ElemID.LE.PP_nElems)THEN
          IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
        END IF
        IF(ElemID.LE.1)THEN
          IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
        END IF
        IF(GlobalProcID.EQ.PartHaloElemToProc(NATIVE_PROC_ID,ElemID))THEN
          ! caution: 
          ! native-proc-id is id in global list || PartMPI%COMM
          ! local-proc-id is the nth neighbor in SurfCOMM%COMM
          PartHaloSideToProc(NATIVE_PROC_ID,iSide)=GlobalProcID
          PartHaloSideToProc(LOCAL_PROC_ID,iSide) =SurfCOMM%nMPINeighbors
        END IF
      END DO ! iSide=nSides+1,nTotalSides
    END DO ! iProc = 0, SurfCOMM%nProcs-1
  END IF
END IF


! build SurfMesh exchange information
! next, allocate SurfCOMM%MPINeighbor
ALLOCATE(SurfCOMM%MPINeighbor(1:SurfCOMM%nMPINeighbors))
! set native proc-id of each SurfCOMM-MPI-Neighbor
LocalProcID=0
DO iProc = 0,SurfCOMM%nProcs-1
  IF(iProc.EQ.SurfCOMM%MyRank) CYCLE
  IF(isMPINeighbor(iProc))THEN
    LocalProcID=LocalProcID+1
    ! map from local proc id to global
    SurfCOMM%MPINeighbor(LocalProcID)%NativeProcID=iProc !PartMPI%MPINeighbor(iProc)
  END IF
END DO ! iProc=1,PartMPI%nMPINeighbors

! array how many data has to be communicated
! number of Sides
ALLOCATE(SurfExchange%nSidesSend(1:SurfCOMM%nMPINeighbors) &
        ,SurfExchange%nSidesRecv(1:SurfCOMM%nMPINeighbors) &
        ,SurfExchange%SendRequest(SurfCOMM%nMPINeighbors)  &
        ,SurfExchange%RecvRequest(SurfCOMM%nMPINeighbors)  )


SurfExchange%nSidesSend=0
SurfExchange%nSidesRecv=0
! loop over all neighbors  
DO iProc=1,SurfCOMM%nMPINeighbors
  ! proc-id in SurfCOMM%nProcs
  DO iSide=nSides+1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    IF(SurfSideID.EQ.-1) CYCLE
    ! get elemid
    IF(iProc.EQ.PartHaloSideToProc(LOCAL_PROC_ID,iSide))THEN
      SurfExchange%nSidesSend(iProc)=SurfExchange%nSidesSend(iProc)+1
      PartHaloSideToProc(LOCAL_SEND_ID,iSide) =SurfExchange%nSidesSend(iProc)
    END IF
  END DO ! iSide=nSides+1,nTotalSides
END DO ! iProc=1,SurfCOMM%nMPINeighbors

! open receive number of send particles
ALLOCATE(RECV_STATUS_LIST(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  CALL MPI_IRECV( SurfExchange%nSidesRecv(iProc)            &
                , 1                                         &
                , MPI_INTEGER                               &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID  &
                , 1001                                      &
                , SurfCOMM%COMM                             &
                , SurfExchange%RecvRequest(iProc)           & 
                , IERROR )
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%SendList(SurfExchange%nSidesSend(iProc)))
  SurfCOMM%MPINeighbor(iProc)%SendList=0
  CALL MPI_ISEND( SurfExchange%nSidesSend(iProc)           &
                , 1                                        &
                , MPI_INTEGER                              &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID &
                , 1001                                     &
                , SurfCOMM%COMM                            &
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )
END DO ! iProc


! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
END DO ! iProc

! allocate send and receive buffer
ALLOCATE(SurfSendBuf(SurfCOMM%nMPINeighbors))
ALLOCATE(SurfRecvBuf(SurfCOMM%nMPINeighbors))
DO iProc=1,SurfCOMM%nMPINeighbors
  ALLOCATE(SurfSendBuf(iProc)%content(2*SurfExchange%nSidesSend(iProc)),STAT=ALLOCSTAT)
  SurfSendBuf(iProc)%Content=-1
  ALLOCATE(SurfRecvBuf(iProc)%content(2*SurfExchange%nSidesRecv(iProc)),STAT=ALLOCSTAT)
END DO ! iProc=1,PartMPI%nMPINeighbors
 
! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                , 2*SurfExchange%nSidesRecv(iProc)             &
                , MPI_DOUBLE_PRECISION                         &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1004                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              & 
                , IERROR )
END DO ! iProc

! build message 
! after this message, the receiving process knows to which of his sides the sending process will send the 
  ! surface data
DO iProc=1,SurfCOMM%nMPINeighbors
  iSendSide=0
  iPos=1
  DO iSide=nSides+1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    IF(SurfSideID.EQ.-1) CYCLE
    ! get elemid
    IF(iProc.EQ.PartHaloSideToProc(LOCAL_PROC_ID,iSide))THEN
      iSendSide=iSendSide+1
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(ElemID.LE.PP_nElems)THEN
        IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
      END IF
      IF(ElemID.LE.1)THEN
        IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
      END IF
      SurfCOMM%MPINeighbor(iProc)%SendList(iSendSide)=SurfSideID
      SurfSendBuf(iProc)%content(iPos  )= REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID))
      SurfSendBuf(iProc)%content(iPos+1)= REAL(PartSideToElem(S2E_LOC_SIDE_ID,iSide))
      iPos=iPos+2
    END IF
  END DO ! iSide=nSides+1,nTotalSides
  IF(iSendSide.NE.SurfExchange%nSidesSend(iProc)) CALL abort(&
__STAMP__&
          ,' Message too short!',iProc)
  IF(ANY(SurfSendBuf(iProc)%content.LE.0))THEN  
    CALL abort(&
__STAMP__&
          ,' Sent NATIVE_ELEM_ID or LOCSIDEID is zero!')
  END IF
END DO

DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
                , 2*SurfExchange%nSidesSend(iProc)         & 
                , MPI_DOUBLE_PRECISION                     &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                , 1004                                     &
                , SurfCOMM%COMM                            &   
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )                                     
END DO ! iProc                                                

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! fill list with received side ids
! store the receiving data
DO iProc=1,SurfCOMM%nMPINeighbors
  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList(SurfExchange%nSidesRecv(iProc)))
  iPos=1
  DO iRecvSide=1,SurfExchange%nSidesRecv(iProc)
    SideID=PartElemToSide(E2S_SIDE_ID,INT(SurfRecvBuf(iProc)%content(iPos+1)),INT(SurfRecvBuf(iProc)%content(iPos)))
    IF(SideID.GT.nSides)THEN
      IPWRITE(UNIT_stdOut,*) ' Received wrong sideid ', SideID
      IPWRITE(UNIT_stdOut,*) ' Sending process has error in halo-region! ', SurfCOMM%MPINeighbor(iProc)%NativeProcID
    END IF
    IF(SideID.GT.nBCSides)THEN
      IPWRITE(UNIT_stdOut,*) ' Received wrong sideid. Is not a BC side! ', SideID
      IPWRITE(UNIT_stdOut,*) ' Sending process has error in halo-region! ', SurfCOMM%MPINeighbor(iProc)%NativeProcID
    END IF
    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
    IF(SurfSideID.EQ.-1)THEN
      IPWRITE(UNIT_stdOut,*) ' Received wrong sideid. SurfSideID is corrupted! '
    END IF
    SurfCOMM%MPINeighbor(iProc)%RecvList(iRecvSide)=SurfSideID
    iPos=iPos+2
  END DO ! RecvSide=1,SurfExchange%nSidesRecv(iProc)-1,2
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
  DEALLOCATE(SurfSendBuf(iProc)%content)
  DEALLOCATE(SurfRecvBuf(iProc)%content)
  ALLOCATE(SurfSendBuf(iProc)%content(SurfMesh%SampSize*nDOF*SurfExchange%nSidesSend(iProc)))
  ALLOCATE(SurfRecvBuf(iProc)%content(SurfMesh%SampSize*nDOF*SurfExchange%nSidesRecv(iProc)))
  SurfSendBuf(iProc)%content=0.
  SurfRecvBuf(iProc)%content=0.
END DO ! iProc
DEALLOCATE(recv_status_list)

CALL MPI_BARRIER(SurfCOMM%Comm,iError)

END SUBROUTINE GetHaloSurfMapping


SUBROUTINE ExchangeSurfData() 
!===================================================================================================================================
! exchange the surface data
! only processes with samling sides in their halo region and the original process participate on the communication
! structure is similar to particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample,SampWall
USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID
INTEGER                         :: iPos,p,q,iProc
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nValues=SurfMesh%SampSize*(nSurfSample)**2
!
! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                , MessageSize                                  &
                , MPI_DOUBLE_PRECISION                         &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1009                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              & 
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SurfSendBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)= SampWall(SurfSideID)%State(:,p,q)
        iPos=iPos+SurfMesh%SampSize
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
    SampWall(SurfSideID)%State(:,:,:)=0.
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
                , MessageSize                              & 
                , MPI_DOUBLE_PRECISION                     &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                , 1009                                     &
                , SurfCOMM%COMM                            &   
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )                                     
END DO ! iProc                                                

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SampWall(SurfSideID)%State(:,p,q)=SampWall(SurfSideID)%State(:,p,q) &
                                         +SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)
        iPos=iPos+SurfMesh%SampSize
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
    SurfRecvBuf(iProc)%content = 0.
    SurfSendBuf(iProc)%content = 0.
 END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO ! iProc

END SUBROUTINE ExchangeSurfData
#endif /*MPI*/

SUBROUTINE WriteSurfSampleToHDF5(MeshFileName,OutputTime) 
!===================================================================================================================================
! write the final values of the surface sampling to a HDF5 state file
! additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Globals_Vars,               ONLY:ProjectName
USE MOD_Particle_Boundary_Vars,     ONLY:nSurfSample,SurfMesh,offSetSurfSide
USE MOD_DSMC_Vars,                  ONLY:MacroSurfaceVal , CollisMode
USE MOD_Particle_Vars,              ONLY:nSpecies
USE MOD_HDF5_Output,                ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_Particle_Boundary_Vars,     ONLY:SurfCOMM,nSurfBC,SurfBCName
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255),ALLOCATABLE      :: StrVarNames(:)
INTEGER                             :: nVar
REAL                                :: tstart,tend
!===================================================================================================================================

#ifdef MPI
CALL MPI_BARRIER(SurfCOMM%COMM,iERROR)
IF(SurfMesh%nSides.EQ.0) RETURN
#endif /*MPI*/
IF(SurfCOMM%MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMCSurfSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString=TRIM(FileName)//'.h5'


! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
IF(SurfCOMM%MPIOutputRoot)THEN
#ifdef MPI
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.)
#else
  CALL OpenDataFile(FileString,create=.TRUE.)
#endif
  Statedummy = 'DSMCSurfState'  
  
  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  ! Write dataset properties "Time","MeshFile","DSMC_nSurfSampl","DSMC_nSpecies","DSMC_CollisMode"
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
  CALL WriteAttributeToHDF5(File_ID,'BC_Surf',nSurfBC,StrArray=SurfBCName)

  ! Create dataset attribute "VarNames" and write to file
  nVar=5
  ALLOCATE(StrVarNames(1:nVar))
  StrVarnames(1)='ForceX'
  StrVarnames(2)='ForceY'
  StrVarnames(3)='ForceZ'
  StrVarnames(4)='HeatFlux'
  StrVarnames(5)='Counter'

  CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarnames)

  CALL CloseDataFile()
  DEALLOCATE(StrVarNames)
END IF

#ifdef MPI
CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
#endif /*MPI*/

#ifdef MPI
!   IF(SurfCOMM%nProcs.GT.1)THEN
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
!   ELSE
!     CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.)
!   END IF
#else
  CALL OpenDataFile(FileString,create=.FALSE.)
#endif

CALL WriteArrayToHDF5(DataSetName='DSMC_SurfaceSampling', rank=4,&
                    nValGlobal=(/5,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/5,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/0,          0,          0,offsetSurfSide/),&
                    collective=.TRUE., RealArray=MacroSurfaceVal)
CALL CloseDataFile()

IF(SurfCOMM%MPIOutputROOT)THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfSampleToHDF5


SUBROUTINE FinalizeParticleBoundarySampling() 
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars
#ifdef MPI
USE MOD_Particle_MPI_Vars           ,ONLY:SurfSendBuf,SurfRecvBuf,SurfExchange,PartHaloSideToProc
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iProc,iSurfSide
!===================================================================================================================================

SDEALLOCATE(XiEQ_SurfSample)
SDEALLOCATE(SurfMesh%SurfaceArea)
SDEALLOCATE(SurfMesh%SideIDToSurfID)
!SDALLOCATE(SampWall%Energy)
!SDEALLOCATE(SampWall%Force)
!SDEALLOCATE(SampWall%Counter)
DO iSurfSide=1,SurfMesh%nSides
  SDEALLOCATE(SampWall(iSurfSide)%State)
END DO
SDEALLOCATE(SurfBCName)
SDEALLOCATE(SampWall)
#ifdef MPI
SDEALLOCATE(PartHaloSideToProc)
SDEALLOCATE(SurfExchange%nSidesSend )
SDEALLOCATE(SurfExchange%nSidesRecv )
DO iProc=1,SurfCOMM%nMPINeighbors
  SDEALLOCATE(SurfSendBuf(iProc)%content)
  SDEALLOCATE(SurfRecvBuf(iProc)%content)
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%SendList)
  SDEALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList)
END DO ! iProc=1,PartMPI%nMPINeighbors
SDEALLOCATE(SurfSendBuf)
SDEALLOCATE(SurfRecvBuf)
SDEALLOCATE(OffSetSurfSideMPI)
#endif /*MPI*/
END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
