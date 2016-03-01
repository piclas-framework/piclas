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
USE MOD_Mesh_Vars               ,ONLY:NGeo,BC,nSides,nBCSides
USE MOD_ReadInTools             ,ONLY:GETINT
USE MOD_Particle_Boundary_Vars  ,ONLY:nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample,SurfMesh,SampWall
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides
USE MOD_Particle_Vars           ,ONLY:nSpecies
USE MOD_Basis                   ,ONLY:LegendreGaussNodesAndWeights
USE MOD_Particle_Surfaces       ,ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY:BezierControlPoints3D
USE MOD_Particle_Mesh_Vars      ,ONLY:PartBCSideList
#ifdef MPI
USE MOD_Particle_MPI_Vars       ,ONLY:PartMPI
#else
USE MOD_Particle_Boundary_Vars  ,ONLY:offSetSurfSide
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
INTEGER                                :: iSample,jSample
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
REAL                                   :: tmpI1,tmpJ1,tmpI2,tmpJ2,XiOut(1:2),E,F,G,D,tmp1
!CHARACTER(2)                :: hilf
!===================================================================================================================================
 
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'
nSurfSample=nGeo
!WRITE(UNIT=hilf,FMT='(I2)') NGeo
!nSurfSample = GETINT('DSMC-nSurfSample',hilf)
 
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))

dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1. 
END DO

! get number of BC-Sides
ALLOCATE(SurfMesh%SideIDToSurfID(1:nTotalSides)) 
! first own sides
SurfMesh%nSides=0
DO iSide=1,nBCSides
  IF(BC(iSide).LE.1) CYCLE
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
    SurfMesh%nSides = SurfMesh%nSides + 1
    SurfMesh%SideIDToSurfID(iSide)=SurfMesh%nSides
    SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nSides
  END IF
END DO

! halo sides
SurfMesh%nTotalSides=SurfMesh%nSides
DO iSide=nSides+1,nTotalSides
  IF(BC(iSide).LE.1) CYCLE
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
    SurfMesh%nTotalSides = SurfMesh%nTotalSides + 1
    SurfMesh%SideIDToSurfID(iSide)=SurfMesh%nTotalSides
    SurfMesh%SideIDToSurfID(iSide) = SurfMesh%nTotalSides
  END IF
END DO

SurfMesh%SurfOnProc=.FALSE.
IF(SurfMesh%nSides.GT.0) SurfMesh%SurfOnProc=.TRUE.

#ifdef MPI
CALL MPI_ALLREDUCE(SurfMesh%nSides,SurfMesh%nGlobalSides,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,iError)
#else
SurfMesh%nGlobalSides=SurfMesh%nSides
#endif


SWRITE(UNIT_stdOut,'(A,I8)') ' nGlobalSurfSides ', SurfMesh%nGlobalSides


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
! get correct offsets
OffSetSurfSide=0
#endif /*MPI*/

IF(.NOT.SurfMesh%SurfOnProc) RETURN


! allocate everything
ALLOCATE(SampWall(1:SurfMesh%nTotalSides))

SurfMesh%SampSize=9+3+nSpecies ! Energy + Force + nSpecies
DO iSide=1,SurfMesh%nTotalSides
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
  SideID=PartBCSideList(iSide)
  ! call here stephens algorithm to compute area 
  DO jSample=1,nSurfSample
    DO iSample=1,nSurfSample
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
          SurfMesh%SurfaceArea(iSample,jSample,SurfSideID)=tmp1*tmp1*D*wGP_NGeo(p)*wGP_NGeo(q)      
        END DO
      END DO
    END DO ! iSample=1,nSurfSample
  END DO ! jSample=1,nSurfSample
END DO ! iSide=1,nTotalSides

DEALLOCATE(Xi_NGeo,wGP_NGeo)

SWRITE(UNIT_stdOut,'(A)') ' ... DONE.'

END SUBROUTINE InitParticleBoundarySampling

#ifdef MPI
SUBROUTINE InitSurfCommunicator()
!===================================================================================================================================
! Read RP parameters from ini file and RP definitions from HDF5 
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
!===================================================================================================================================
color=MPI_UNDEFINED
IF(SurfMesh%SurfonProc) color=2

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

! create new RP communicator for RP output
CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyRank, SurfCOMM%COMM,iError)
IF(SurfMesh%SurfOnPRoc) CALL MPI_COMM_SIZE(SurfCOMM%COMM, SurfCOMM%nProcs,iError)
SurfCOMM%MPIRoot=.FALSE.
IF(SurfCOMM%MyRank.EQ.0 .AND. SurfMesh%SurfOnProc) THEN
  SurfCOMM%MPIRoot=.TRUE.
  WRITE(UNIT_stdout,'(A11,I5,A6)') 'SURF COMM: ',SurfCOMM%nProcs,' procs'
END IF

! get correct offsets
ALLOCATE(offsetSurfSideMPI(0:SurfCOMM%nProcs))
offsetSurfSideMPI=0
ALLOCATE(countSurfSideMPI(0:SurfCOMM%nProcs-1))
countSurfSideMPI=0

CALL MPI_GATHER(SurfMesh%nSides,1,MPI_INTEGER,countSurfSideMPI,1,MPI_INTEGER,0,SurfCOMM%COMM,iError)

IF (SurfCOMM%MPIroot) THEN
  DO iProc=1,SurfCOMM%nProcs-1
    offsetSurfSideMPI(iProc)=SUM(countSurfSideMPI(0:iProc-1))
  END DO
  offsetSurfSideMPI(SurfCOMM%nProcs)=SUM(countSurfSideMPI(:))
END IF
CALL MPI_BCAST (offsetSurfSideMPI,size(offsetSurfSideMPI),MPI_INTEGER,0,SurfCOMM%COMM,iError)
offsetSurfSide=offsetSurfSideMPI(SurfCOMM%MyRank)

END SUBROUTINE InitSurfCommunicator


SUBROUTINE GetHaloSurfMapping() 
!===================================================================================================================================
! build all missing stuff for surface-sampling communicator, like
! offSetMPI
! MPI-neighbor list
! PartHaloSideToProc
! only receiving process nows to which local side the sending information is going, the sending process does not know the final 
! sideid 
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars                   ,ONLY:nSides
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfMesh,SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars           ,ONLY:PartHaloSideToProc,PartMPI,PartHaloElemToProc,SurfSendBuf,SurfRecvBuf,SurfExchange
USE MOD_Particle_Mesh_Vars          ,ONLY:nTotalSides,PartSideToElem,PartElemToSide
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: isMPINeighbor(0:PartMPI%nMPINeighbors)
INTEGER                           :: nDOF,ALLOCSTAT,SideID
INTEGER                           :: iProc, GlobalProcID,iSide,ElemID,SurfSideID,LocalProcID,iSendSide,iRecvSide,iPos
INTEGER,ALLOCATABLE               :: recv_status_list(:,:)
!===================================================================================================================================

nDOF=(nSurfSample)*(nSurfSample)

IF(SurfMesh%nTotalSides.GT.SurfMesh%nSides)THEN
  ALLOCATE(PartHaloSideToProc(1:4,SurfMesh%nSides+1:SurfMesh%nTotalSides))
  PartHaloSideToProc=-1
  ! get all MPI-neighbors to communicate with
  ! loop over all mpi neighbors
  isMPINeighbor=.FALSE.
  SurfCOMM%nMPINeighbors=0
  DO iProc=1,PartMPI%nMPINeighbors
    GlobalProcID=PartMPI%MPINeighbor(iProc)
    DO iSide=nSides+1,nTotalSides
      SurfSideID=SurfMesh%SideIDToSurfID(iSide)
      IF(SurfSideID.EQ.-1) CYCLE
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(GlobalProcID.EQ.PartHaloElemToProc(NATIVE_PROC_ID,ElemID))THEN
        IF(.NOT.isMPINeighbor(iProc))THEN
          isMPINeighbor(iProc)=.TRUE.
          SurfCOMM%nMPINeighbors=SurfCOMM%nMPINeighbors+1
        END IF
        PartHaloSideToProc(NATIVE_PROC_ID,iSide)=GlobalProcID
        PartHaloSideToProc(LOCAL_PROC_ID,iSide) =PartHaloElemToProc(LOCAL_PROC_ID,ElemID)
      END IF
    END DO ! iSide=nSides+1,nTotalSides
  END DO
  ! next, allocate SurfCOMM%MPINeighbor
  ALLOCATE(SurfCOMM%MPINeighbor(1:SurfCOMM%nMPINeighbors))
  LocalProcID=0
  DO iProc = 1,PartMPI%nMPINeighbors
    IF(isMPINeighbor(iProc))THEN
      LocalProcID=LocalProcID+1
      ! map from local proc id to global
      SurfCOMM%MPINeighbor(LocalProcID)%NativeProcID=PartMPI%MPINeighbor(iProc)
    END IF
  END DO ! iProc=1,PartMPI%nMPINeighbors
  ALLOCATE(SurfExchange%nSidesSend(1:SurfCOMM%nMPINeighbors) &
          ,SurfExchange%nSidesRecv(1:SurfCOMM%nMPINeighbors) &
          ,SurfExchange%SendRequest(SurfCOMM%nMPINeighbors)  &
          ,SurfExchange%RecvRequest(SurfCOMM%nMPINeighbors)  )

  ALLOCATE(recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors))

  SurfExchange%nSidesSend=0
  SurfExchange%nSidesRecv=0
  ! loop over all neighbors  
  DO iProc=1,PartMPI%nMPINeighbors
    GlobalProcID=PartMPI%MPINeighbor(iProc)
    DO iSide=nSides+1,nTotalSides
      SurfSideID=SurfMesh%SideIDToSurfID(iSide)
      IF(SurfSideID.EQ.-1) CYCLE
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(iProc.EQ.PartHaloElemToProc(LOCAL_PROC_ID,ElemID))THEN
        SurfExchange%nSidesSend(iProc)=SurfExchange%nSidesSend(iProc)+1
        PartHaloSideToProc(LOCAL_SEND_ID,SideID) =SurfExchange%nSidesSend(iProc)
      END IF
    END DO ! iSide=nSides+1,nTotalSides
  END DO
  ! open receive number of send particles
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
    ALLOCATE(SurfRecvBuf(iProc)%content(2*SurfExchange%nSidesRecv(iProc)),STAT=ALLOCSTAT)
  END DO ! iProc=1,PartMPI%nMPINeighbors
 
  ! open receive buffer
  DO iProc=1,SurfCOMM%nMPINeighbors
    CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                  , 2*SurfExchange%nSidesRecv(iProc)             &
                  , MPI_INTEGER                                  &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                  , 1004                                         &
                  , SurfCOMM%COMM                                &
                  , SurfExchange%RecvRequest(iProc)              & 
                  , IERROR )
  END DO ! iProc

  ! build message
  DO iProc=1,SurfCOMM%nMPINeighbors
    iSendSide=0
    iPos=1
    DO iSide=nSides+1,nTotalSides
      SurfSideID=SurfMesh%SideIDToSurfID(iSide)
      IF(SurfSideID.EQ.-1) CYCLE
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
      IF(iProc.EQ.PartHaloElemToProc(LOCAL_PROC_ID,ElemID))THEN
        iSendSide=iSendSide+1
        SurfCOMM%MPINeighbor(iProc)%SendList(iSendSide)=iSide
        SurfSendBuf(iProc)%content(iPos  )= PartHaloElemToProc(NATIVE_ELEM_ID,ElemID)
        SurfSendBuf(iProc)%content(iPos+1)= PartSideToElem(S2E_LOC_SIDE_ID,iSide)
      END IF
    END DO ! iSide=nSides+1,nTotalSides
  END DO


  DO iProc=1,SurfCOMM%nMPINeighbors
    CALL MPI_ISEND( SurfSendBuf(iProc)%content               &
                  , 2*SurfExchange%nSidesSend(iProc)         & 
                  , MPI_INTEGER                              &
                  , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                  , 1004                                     &
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

  ! fill list with received side ids
  DO iProc=1,SurfCOMM%nMPINeighbors
    ALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList(SurfExchange%nSidesRecv(iProc)))
    iPos=1
    DO iRecvSide=1,SurfExchange%nSidesRecv(iProc)
      SideID=PartElemToSide(E2S_SIDE_ID,INT(SurfRecvBuf(iProc)%content(iPos+1)),INT(SurfRecvBuf(iProc)%content(iPos)))
      SurfSideID=SurfMesh%SideIDToSurfID(SideID)
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

END IF

END SUBROUTINE GetHaloSurfMapping


SUBROUTINE ExchangeSurfData() 
!===================================================================================================================================
! exchange the surface data like particles
! has to be added to the already computed values 
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
  CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
      __STAMP__&
          ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
      __STAMP__&
          ,' MPI Communication error', IERROR)
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SampWall(SurfSideID)%State(:,p,q)=SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)
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
USE MOD_Output_Vars,                ONLY:ProjectName
USE MOD_Particle_Boundary_Vars,     ONLY:nSurfSample,SurfMesh,offSetSurfSide
USE MOD_DSMC_Vars,                  ONLY:MacroSurfaceVal , CollisMode
USE MOD_Particle_Vars,              ONLY:nSpecies
USE MOD_HDF5_Output,                ONLY:WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
#ifdef MPI
USE MOD_Particle_Boundary_Vars,     ONLY:SurfCOMM
#endif
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
!===================================================================================================================================

SWRITE(*,*) ' WRITE DSMCSurfSTATE TO HDF5 FILE...'
FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString=TRIM(FileName)//'.h5'
#ifdef MPI
  CALL OpenDataFile(FileString,create=.TRUE.,single=.FALSE.,communicatorOpt=SurfCOMM%COMM)
#else
  CALL OpenDataFile(FileString,create=.TRUE.)
#endif

Statedummy = 'DSMCSurfState'
CALL WriteHDF5Header(Statedummy,File_ID)


CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=nSurfSample-1)
CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)


nVar=5
ALLOCATE(StrVarNames(1:nVar))
StrVarnames(1)='ForceX'
StrVarnames(2)='ForceY'
StrVarnames(3)='ForceZ'
StrVarnames(4)='HeatFlux'
StrVarnames(5)='Counter'

CALL WriteArrayToHDF5(DataSetName='DSMC_SurfaceSampling', rank=4,&
                    nValGlobal=(/5,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/5,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/0,          0,          0,offsetSurfSide/),&
                    collective=.TRUE., RealArray=MacroSurfaceVal)
CALL CloseDataFile()

DEALLOCATE(StrVarNames)

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
