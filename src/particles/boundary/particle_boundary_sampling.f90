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
USE MOD_ReadInTools
USE MOD_Mesh_Vars               ,ONLY:NGeo,BC,nSides,nBCSides,nBCs,BoundaryName
USE MOD_ReadInTools             ,ONLY:GETINT
USE MOD_Particle_Boundary_Vars  ,ONLY:nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample,SurfMesh,SampWall,nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars  ,ONLY:SurfCOMM,CalcSurfCollis,AnalyzeSurfCollis
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides
USE MOD_Particle_Vars           ,ONLY:nSpecies
USE MOD_DSMC_Vars               ,ONLY:useDSMC,DSMC,Adsorption
USE MOD_Basis                   ,ONLY:LegendreGaussNodesAndWeights
USE MOD_Particle_Surfaces       ,ONLY:EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Surfaces_Vars  ,ONLY:BezierControlPoints3D,BezierSampleN
USE MOD_Particle_Mesh_Vars      ,ONLY:PartBCSideList
USE MOD_Particle_Tracking_Vars  ,ONLY:DoRefMapping
#ifdef MPI
USE MOD_Particle_MPI_Vars       ,ONLY:PartMPI
#else
USE MOD_Particle_Boundary_Vars  ,ONLY:offSetSurfSide
#endif /*MPI*/
USE MOD_PICDepo_Vars            ,ONLY:SFResampleAnalyzeSurfCollis
USE MOD_PICDepo_Vars            ,ONLY:LastAnalyzeSurfCollis
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: p,q,iSide,SurfSideID,SideID
INTEGER                                :: iSample,jSample, iBC, iSpec
REAL,DIMENSION(2,3)                    :: gradXiEta3D
REAL,ALLOCATABLE,DIMENSION(:)          :: Xi_NGeo,wGP_NGeo
REAL                                   :: XiOut(1:2),E,F,G,D,tmp1,area,tmpI2,tmpJ2
CHARACTER(2)                           :: hilf, hilf2
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
INTEGER,ALLOCATABLE                    :: CalcSurfCollis_SpeciesRead(:) !help array for reading surface stuff
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
  BCName=''
END DO
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

ALLOCATE(SurfMesh%SurfSideToGlobSideMap(1:SurfMesh%nTotalSides))
SurfMesh%SurfSideToGlobSideMap(:) = -1
DO iSide = 1,nTotalSides
  IF (SurfMesh%SideIDToSurfID(iSide).LE.0) CYCLE
  SurfMesh%SurfSideToGlobSideMap(SurfMesh%SideIDToSurfID(iSide)) = iSide
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

! get the full area of the BC's of all faces
Area=0.
DO iSide=1,nBCSides
  SurfSideID=SurfMesh%SideIDToSurfID(iSide)
  IF(SurfSideID.EQ.-1) CYCLE
  Area = Area + SUM(SurfMesh%SurfaceArea(:,:,SurfSideID))
END DO ! iSide=1,nTotalSides

#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Area,1,MPI_DOUBLE_PRECISION,MPI_SUM,SurfCOMM%COMM,iError)
#endif /*MPI*/

SWRITE(UNIT_stdOut,'(A,E25.14E3)') ' Surface-Area: ', Area

DEALLOCATE(Xi_NGeo,wGP_NGeo)

! Initialize surface collision sampling and analyze
CalcSurfCollis%OnlySwaps = GETLOGICAL('Particles-CalcSurfCollis_OnlySwaps','.FALSE.')
CalcSurfCollis%Only0Swaps = GETLOGICAL('Particles-CalcSurfCollis_Only0Swaps','.FALSE.')
CalcSurfCollis%Output = GETLOGICAL('Particles-CalcSurfCollis_Output','.FALSE.')
IF (CalcSurfCollis%Only0Swaps) CalcSurfCollis%OnlySwaps=.TRUE.
CalcSurfCollis%AnalyzeSurfCollis = GETLOGICAL('Particles-AnalyzeSurfCollis','.FALSE.')
AnalyzeSurfCollis%NumberOfBCs = 1 !initialize for ifs (BCs=0 means all)
ALLOCATE(AnalyzeSurfCollis%BCs(1))
AnalyzeSurfCollis%BCs = 0
IF (.NOT.CalcSurfCollis%AnalyzeSurfCollis .AND. SFResampleAnalyzeSurfCollis) THEN
  CALL abort(__STAMP__,&
    'ERROR: SFResampleAnalyzeSurfCollis was set without CalcSurfCollis%AnalyzeSurfCollis!')
END IF
IF (CalcSurfCollis%AnalyzeSurfCollis) THEN
  AnalyzeSurfCollis%maxPartNumber = GETINT('Particles-DSMC-maxSurfCollisNumber','0')
  AnalyzeSurfCollis%NumberOfBCs = GETINT('Particles-DSMC-NumberOfBCs','1')
  IF (AnalyzeSurfCollis%NumberOfBCs.EQ.1) THEN !already allocated
    AnalyzeSurfCollis%BCs = GETINT('Particles-DSMC-SurfCollisBC','0') ! 0 means all...
  ELSE
    DEALLOCATE(AnalyzeSurfCollis%BCs)
    ALLOCATE(AnalyzeSurfCollis%BCs(1:AnalyzeSurfCollis%NumberOfBCs)) !dummy
    hilf2=''
    DO iBC=1,AnalyzeSurfCollis%NumberOfBCs !build default string: 0,0,0,...
      WRITE(UNIT=hilf,FMT='(I0)') 0
      hilf2=TRIM(hilf2)//TRIM(hilf)
      IF (iBC.NE.AnalyzeSurfCollis%NumberOfBCs) hilf2=TRIM(hilf2)//','
    END DO
    AnalyzeSurfCollis%BCs = GETINTARRAY('Particles-SurfCollisBC',AnalyzeSurfCollis%NumberOfBCs,hilf2)
  END IF
  ALLOCATE(AnalyzeSurfCollis%Data(1:AnalyzeSurfCollis%maxPartNumber,1:9))
  ALLOCATE(AnalyzeSurfCollis%Spec(1:AnalyzeSurfCollis%maxPartNumber))
  ALLOCATE(AnalyzeSurfCollis%BCid(1:AnalyzeSurfCollis%maxPartNumber))
  ALLOCATE(AnalyzeSurfCollis%Number(1:nSpecies+1))
  IF (LastAnalyzeSurfCollis%Restart) THEN
    CALL ReadAnalyzeSurfCollisToHDF5()
  END IF
  !ALLOCATE(AnalyzeSurfCollis%Rate(1:nSpecies+1))
  AnalyzeSurfCollis%Data=0.
  AnalyzeSurfCollis%Spec=0
  AnalyzeSurfCollis%BCid=0
  AnalyzeSurfCollis%Number=0
  !AnalyzeSurfCollis%Rate=0.
END IF
! Species-dependent calculations
ALLOCATE(CalcSurfCollis%SpeciesFlags(1:nSpecies))
CalcSurfCollis%NbrOfSpecies = GETINT('Particles-DSMC-CalcSurfCollis_NbrOfSpecies','0')
IF ( (CalcSurfCollis%NbrOfSpecies.GT.0) .AND. (CalcSurfCollis%NbrOfSpecies.LE.nSpecies) ) THEN
  ALLOCATE(CalcSurfCollis_SpeciesRead(1:CalcSurfCollis%NbrOfSpecies))
  hilf2=''
  DO iSpec=1,CalcSurfCollis%NbrOfSpecies !build default string: 1 - CSC_NoS
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    hilf2=TRIM(hilf2)//TRIM(hilf)
    IF (ispec.NE.CalcSurfCollis%NbrOfSpecies) hilf2=TRIM(hilf2)//','
  END DO
  CalcSurfCollis_SpeciesRead = GETINTARRAY('Particles-CalcSurfCollis_Species',CalcSurfCollis%NbrOfSpecies,hilf2)
  CalcSurfCollis%SpeciesFlags(:)=.FALSE.
  DO iSpec=1,CalcSurfCollis%NbrOfSpecies
    CalcSurfCollis%SpeciesFlags(CalcSurfCollis_SpeciesRead(ispec))=.TRUE.
  END DO
  DEALLOCATE(CalcSurfCollis_SpeciesRead)
ELSE IF (CalcSurfCollis%NbrOfSpecies.EQ.0) THEN !default
  CalcSurfCollis%SpeciesFlags(:)=.TRUE.
ELSE
  CALL abort(&
  __STAMP__&
  ,'Error in Particles-CalcSurfCollis_NbrOfSpecies!')
END IF

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
! Sets also used for communication of adsorption variables
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

! create new SurfMesh communicator for SurfMesh communication
CALL MPI_COMM_SPLIT(PartMPI%COMM, color, SurfCOMM%MyRank, SurfCOMM%COMM,iError)
IF(SurfMesh%SurfOnPRoc) THEN
  CALL MPI_COMM_SIZE(SurfCOMM%COMM, SurfCOMM%nProcs,iError)
ELSE
  SurfCOMM%nProcs = 0
END IF
SurfCOMM%MPIRoot=.FALSE.
IF(SurfCOMM%MyRank.EQ.0 .AND. SurfMesh%SurfOnProc) THEN
  SurfCOMM%MPIRoot=.TRUE.
!   WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF COMM:        ',SurfCOMM%nProcs,' procs'
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
!   WRITE(UNIT_stdout,'(A18,I5,A6)') 'SURF OUTPUT-COMM: ',SurfCOMM%nOutputProcs,' procs'
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
USE MOD_Particle_Vars               ,ONLY:nSpecies
USE MOD_Particle_MPI_Vars           ,ONLY:PartHaloSideToProc,PartHaloElemToProc,SurfSendBuf,SurfRecvBuf,SurfExchange
USE MOD_Particle_Mesh_Vars          ,ONLY:nTotalSides,PartSideToElem,PartElemToSide
USE MOD_DSMC_Vars                   ,ONLY:Adsorption
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
INTEGER                           :: NativeElemID, NativeLocSideID
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
        CALL abort(&
__STAMP__&
,' Error in PartSideToElem. Halo-Side cannot be connected to local element', ElemID  )
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
DO iProc=0,SurfCOMM%nProcs-1
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
        CALL abort(&
__STAMP__&
,' Error in PartSideToElem. Halo-Side cannot be connected to local element', ElemID  )
        END IF
        IF(GlobalProcID.EQ.PartHaloElemToProc(NATIVE_PROC_ID,ElemID))THEN
          ! caution: 
          ! native-proc-id is id in global list || PartMPI%COMM
          ! local-proc-id is the nth neighbor in SurfCOMM%COMM
          PartHaloSideToProc(NATIVE_PROC_ID,iSide)=GlobalProcID
          PartHaloSideToProc(LOCAL_PROC_ID ,iSide)=SurfCOMM%nMPINeighbors
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
  IF(SurfExchange%nSidesSend(iProc).GT.0)THEN
    ALLOCATE(SurfSendBuf(iProc)%content(2*SurfExchange%nSidesSend(iProc)),STAT=ALLOCSTAT)
    SurfSendBuf(iProc)%Content=-1
  END IF
  IF(SurfExchange%nSidesRecv(iProc).GT.0)THEN
    ALLOCATE(SurfRecvBuf(iProc)%content(2*SurfExchange%nSidesRecv(iProc)),STAT=ALLOCSTAT)
  END IF
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
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iSendSide=0
  iPos=1
  DO iSide=nSides+1,nTotalSides
    SurfSideID=SurfMesh%SideIDToSurfID(iSide)
    IF(SurfSideID.EQ.-1) CYCLE
    ! get elemid
    IF(iProc.EQ.PartHaloSideToProc(LOCAL_PROC_ID,iSide))THEN
      iSendSide=iSendSide+1
      ! get elemid
      ElemID=PartSideToElem(S2E_ELEM_ID,iSide)
!       IF(ElemID.LE.PP_nElems)THEN
!         IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
!       END IF
!       IF(ElemID.LE.1)THEN
!         IPWRITE(UNIT_stdOut,*) ' Error in PartSideToElem'
!       END IF
      SurfCOMM%MPINeighbor(iProc)%SendList(iSendSide)=SurfSideID
!       IPWRITE(*,*) 'negative elem id',PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),PartSideToElem(S2E_LOC_SIDE_ID,iSide)
      SurfSendBuf(iProc)%content(iPos  )= REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID))
      SurfSendBuf(iProc)%content(iPos+1)= REAL(PartSideToElem(S2E_LOC_SIDE_ID,iSide))
      iPos=iPos+2
    END IF
  END DO ! iSide=nSides+1,nTotalSides
  IF(iSendSide.NE.SurfExchange%nSidesSend(iProc)) CALL abort(&
__STAMP__&
          ,' Message too short!',iProc)
  IF(ANY(SurfSendBuf(iProc)%content.LE.0))THEN  
    IPWRITE(UNIT_stdOut,*) ' nSendSides', SurfExchange%nSidesSend(iProc), ' to Proc ', iProc
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
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  ALLOCATE(SurfCOMM%MPINeighbor(iProc)%RecvList(SurfExchange%nSidesRecv(iProc)))
  iPos=1
  DO iRecvSide=1,SurfExchange%nSidesRecv(iProc)
    NativeElemID   =INT(SurfRecvBuf(iProc)%content(iPos))
    NativeLocSideID=INT(SurfRecvBuf(iProc)%content(iPos+1))
!     IPWRITE(*,*) 'received- elemid, locsideid',NativeElemID,NativeLocSideID
    IF(NativeElemID.GT.PP_nElems)THEN
     CALL abort(&
__STAMP__&
          ,' Cannot send halo-data to other progs. big error! ', ElemID, REAL(PP_nElems))
    END IF
    SideID=PartElemToSide(E2S_SIDE_ID,NativeLocSideID,NativeElemID)
    IF(SideID.GT.nBCSides)THEN
      IPWRITE(UNIT_stdOut,*) ' Received wrong sideid. Is not a BC side! '
      IPWRITE(UNIT_stdOut,*) ' SideID, nBCSides, nSides ', SideID, nBCSides, nSides
      IPWRITE(UNIT_stdOut,*) ' ElemID, locsideid        ', NativeElemID, NativeLocSideID
      IPWRITE(UNIT_stdOut,*) ' Sending process has error in halo-region! ', SurfCOMM%MPINeighbor(iProc)%NativeProcID
     CALL abort(&
__STAMP__&
          ,' Big error in halo region! NativeLocSideID ', NativeLocSideID )
    END IF
    SurfSideID=SurfMesh%SideIDToSurfID(SideID)
    IF(SurfSideID.EQ.-1)THEN
     CALL abort(&
__STAMP__&
          ,' Side is not even a reflective BC side! ', SurfSideID )
    END IF
    SurfCOMM%MPINeighbor(iProc)%RecvList(iRecvSide)=SurfSideID
    iPos=iPos+2
  END DO ! RecvSide=1,SurfExchange%nSidesRecv(iProc)-1,2
END DO ! iProc

DO iProc=1,SurfCOMM%nMPINeighbors
  SDEALLOCATE(SurfSendBuf(iProc)%content)
  SDEALLOCATE(SurfRecvBuf(iProc)%content)
  IF(SurfExchange%nSidesSend(iProc).GT.0) THEN
    ALLOCATE(SurfSendBuf(iProc)%content((SurfMesh%SampSize)*nDOF*SurfExchange%nSidesSend(iProc)))
    SurfSendBuf(iProc)%content=0.
  END IF
  IF(SurfExchange%nSidesRecv(iProc).GT.0) THEN
    ALLOCATE(SurfRecvBuf(iProc)%content((SurfMesh%SampSize)*nDOF*SurfExchange%nSidesRecv(iProc)))
    SurfRecvBuf(iProc)%content=0.
  END IF
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
USE MOD_Particle_Vars               ,ONLY:nSpecies
USE MOD_DSMC_Vars                   ,ONLY:Adsorption,DSMC
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
INTEGER                         :: iPos,p,q,iProc,iReact
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

IF (DSMC%WallModel.GT.0)THEN
  ! additional array entries for Coverage, Accomodation and recombination coefficient
  nValues = (SurfMesh%SampSize+(nSpecies+1)+nSpecies+(Adsorption%RecombNum*nSpecies))*(nSurfSample)**2
ELSE
  nValues = SurfMesh%SampSize*nSurfSample**2
END IF
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
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  SurfSendBuf(iProc)%content = 0.
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SurfSendBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)= SampWall(SurfSideID)%State(:,p,q)
        iPos=iPos+SurfMesh%SampSize
        IF (DSMC%WallModel.GT.0)THEN
          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies+1)= SampWall(SurfSideID)%Adsorption(:,p,q)
          iPos=iPos+nSpecies+1
          SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%Accomodation(:,p,q)
          iPos=iPos+nSpecies
          DO iReact=1,Adsorption%RecombNum
            SurfSendBuf(iProc)%content(iPos+1:iPos+nSpecies)= SampWall(SurfSideID)%Reaction(iReact,:,p,q)
            iPos=iPos+nSpecies
          END DO
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
    SampWall(SurfSideID)%State(:,:,:)=0.
    IF (DSMC%WallModel.GT.0)THEN
      SampWall(SurfSideID)%Adsorption(:,:,:)=0.
      SampWall(SurfSideID)%Accomodation(:,:,:)=0.
      SampWall(SurfSideID)%Reaction(:,:,:,:)=0.
    END IF
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
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        SampWall(SurfSideID)%State(:,p,q)=SampWall(SurfSideID)%State(:,p,q) &
                                         +SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfMesh%SampSize)
        iPos=iPos+SurfMesh%SampSize
        IF (DSMC%WallModel.GT.0)THEN
          SampWall(SurfSideID)%Adsorption(:,p,q)=SampWall(SurfSideID)%Adsorption(:,p,q) &
                                                +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies+1)
          iPos=iPos+nSpecies+1
          SampWall(SurfSideID)%Accomodation(:,p,q)=SampWall(SurfSideID)%Accomodation(:,p,q) &
                                                  +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
          iPos=iPos+nSpecies
          DO iReact=1,Adsorption%RecombNum
            SampWall(SurfSideID)%Reaction(iReact,:,p,q)=SampWall(SurfSideID)%Reaction(iReact,:,p,q) &
                                                       +SurfRecvBuf(iProc)%content(iPos+1:iPos+nSpecies)
            iPos=iPos+nSpecies
          END DO
        END IF
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  SurfRecvBuf(iProc)%content = 0.
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
USE MOD_DSMC_Vars,                  ONLY:MacroSurfaceVal,MacroSurfaceSpecVal, CollisMode, DSMC
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
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255),ALLOCATABLE      :: StrVarNames(:), StrOutNames(:)
INTEGER                             :: nVar, Species_nOut, iSpec
REAL                                :: tstart,tend
!===================================================================================================================================

#ifdef MPI
CALL MPI_BARRIER(SurfCOMM%COMM,iERROR)
IF(SurfMesh%nSides.EQ.0) RETURN
#endif /*MPI*/
IF(SurfCOMM%MPIOutputRoot)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMCSurfSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString=TRIM(FileName)//'.h5'


! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#ifdef MPI
IF(SurfCOMM%MPIOutputRoot)THEN
  CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
#else
  CALL OpenDataFile(FileString,create=.TRUE.,readOnly=.FALSE.)
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
  NodeTypeTemp='VISU'
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeTemp/))

  ! Create dataset attribute "VarNames" and write to file
  nVar=5
  ALLOCATE(StrVarNames(1:nVar))
  StrVarNames(1)='ForceX'
  StrVarNames(2)='ForceY'
  StrVarNames(3)='ForceZ'
  StrVarNames(4)='HeatFlux'
  StrVarNames(5)='Counter_Total'
  
  CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)
  
! Create dataset attribute "Species_Varnames" and write to file
  IF (DSMC%WallModel.GT.0) THEN
    Species_nOut=4
    ALLOCATE(StrOutNames(1:Species_nOut))
    StrOutNames(1)='Spec_Counter'
    StrOutNames(2)='Accomodation'
    StrOutNames(3)='Coverage'
    StrOutNames(4)='Recomb_Coeff'
  ELSE
    Species_nOut=1
    ALLOCATE(StrOutNames(1:Species_nOut))
    StrOutNames(1)='Counter'
  END IF
  
  CALL WriteAttributeToHDF5(File_ID,'Species_Varnames',Species_nOut,StrArray=StrOutNames)

  CALL CloseDataFile()
  DEALLOCATE(StrVarNames)
  DEALLOCATE(StrOutNames)
#ifdef MPI
END IF

CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
#endif /*MPI*/

#ifdef MPI
!   IF(SurfCOMM%nProcs.GT.1)THEN
  CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
!   ELSE
!     CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.)
!   END IF
#else
  CALL OpenDataFile(FileString,create=.FALSE.,readOnly=.FALSE.)
#endif

CALL WriteArrayToHDF5(DataSetName='DSMC_SurfaceSampling', rank=4,&
                    nValGlobal=(/5,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                    nVal=      (/5,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                    offset=    (/0,          0,          0,offsetSurfSide/),&
                    collective=.TRUE., RealArray=MacroSurfaceVal)
IF (DSMC%WallModel.GT.0) THEN
  DO iSpec = 1,nSpecies
      WRITE(H5_Name,'(A,I3.3,A)') 'DSMC_Spec',iSpec,'_SurfaceSampling'
      CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=4,&
                      nValGlobal=(/4,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                      nVal=      (/4,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                      offset=    (/0,          0,          0,offsetSurfSide/),&
                      collective=.TRUE.,  RealArray=MacroSurfaceSpecVal(:,:,:,:,iSpec))
  END DO
ELSE
  DO iSpec = 1,nSpecies
      WRITE(H5_Name,'(A,I3.3,A)') 'DSMC_Spec',iSpec,'_SurfaceSampling'
      CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=4,&
                      nValGlobal=(/1,nSurfSample,nSurfSample,SurfMesh%nGlobalSides/),&
                      nVal=      (/1,nSurfSample,nSurfSample,SurfMesh%nSides/),&
                      offset=    (/0,          0,          0,offsetSurfSide/),&
                      collective=.TRUE.,  RealArray=MacroSurfaceSpecVal(:,:,:,:,iSpec))
  END DO
END IF
CALL CloseDataFile()

IF(SurfCOMM%MPIOutputROOT)THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WriteSurfSampleToHDF5


SUBROUTINE ReadAnalyzeSurfCollisToHDF5()
!===================================================================================================================================
! Reading AnalyzeSurfCollis-Data from hdf5 file for restart
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_input,         ONLY: OpenDataFile,CloseDataFile,ReadArray,File_ID,GetDataSize,nDims,HSize,ReadAttribute
USE MOD_Particle_Vars,      ONLY: nSpecies
USE MOD_Particle_Boundary_Vars,ONLY: AnalyzeSurfCollis
USE MOD_PICDepo_Vars,       ONLY: LastAnalyzeSurfCollis, r_SF !, SFResampleAnalyzeSurfCollis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: Filename, H5_Name
INTEGER                        :: PartDataSize, iSpec
LOGICAL                        :: fileExists
REAL, ALLOCATABLE              :: PartSpecData(:,:,:)
INTEGER                        :: TotalNumberMPF, counter2, counter
REAL                           :: TotalFlowrateMPF, RandVal
LOGICAL,ALLOCATABLE            :: PartDone(:)
!===================================================================================================================================
  FileName = TRIM(LastAnalyzeSurfCollis%DSMCSurfCollisRestartFile)
  SWRITE(UNIT_stdOut,*)'Reading Particles from DSMCSurfCollis-File:',TRIM(FileName)

  !-- initialize data (check if file exists and determine size of arrays)
  PartDataSize=9
  IF(MPIRoot) THEN
    INQUIRE (FILE=TRIM(FileName), EXIST=fileExists)
    IF(.NOT.FileExists)  CALL abort(__STAMP__, &
          'DSMCSurfCollis-File "'//TRIM(FileName)//'" does not exist',999,999.)
#ifdef MPI
    CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
#else
    CALL OpenDataFile(TRIM(FileName),create=.FALSE.,readOnly=.TRUE.)
#endif
    DO iSpec=1,nSpecies
      WRITE(H5_Name,'(A,I3.3)') 'SurfCollisData_Spec',iSpec
      CALL GetDataSize(File_ID,TRIM(H5_Name),nDims,HSize)
      AnalyzeSurfCollis%Number(iSpec)=INT(HSize(1),4) !global number of particles
      IF ( INT(HSize(nDims),4) .NE. PartDataSize ) THEN
        CALL Abort(&
        __STAMP__,&
        'Error in ReadAnalyzeSurfCollisToHDF5. Array has size of ',nDims,REAL(INT(HSize(nDims),4)))
      END IF
      DEALLOCATE(HSize)
    END DO !iSpec
    CALL CloseDataFile()
    AnalyzeSurfCollis%Number(nSpecies+1) = SUM( AnalyzeSurfCollis%Number(1:nSpecies) )
  END IF !MPIRoot
#ifdef MPI
  CALL MPI_BCAST(AnalyzeSurfCollis%Number(:),nSpecies+1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#endif
  TotalNumberMPF=AnalyzeSurfCollis%Number(nSpecies+1)
  ALLOCATE( PartSpecData(nSpecies,MAXVAL(AnalyzeSurfCollis%Number(1:nSpecies)),PartDataSize) )

  !-- open file for actual read-in
#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
#else
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,readOnly=.TRUE.)
#endif
  ! Read in state
  DO iSpec=1,nSpecies
    WRITE(H5_Name,'(A,I3.3)') 'SurfCollisData_Spec',iSpec
    IF (AnalyzeSurfCollis%Number(iSpec).GT.0) CALL ReadArray(TRIM(H5_Name),2,(/AnalyzeSurfCollis%Number(iSpec),PartDataSize/), &
      0,1,RealArray=PartSpecData(iSpec,1:AnalyzeSurfCollis%Number(iSpec),1:PartDataSize))
  END DO !iSpec
  CALL ReadAttribute(File_ID,'TotalFlowrateMPF',1,RealScalar=TotalFlowrateMPF)
  SWRITE(UNIT_stdOut,*)'DONE!' 
  CALL CloseDataFile() 

!--save data
  IF (TotalNumberMPF.GT.0) THEN
    IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts to MaxPartNumber
      LastAnalyzeSurfCollis%PartNumberSamp=MIN(TotalNumberMPF,LastAnalyzeSurfCollis%PartNumberReduced)
      ALLOCATE(PartDone(1:TotalNumberMPF))
      PartDone(:)=.FALSE.
    ELSE
      LastAnalyzeSurfCollis%PartNumberSamp=TotalNumberMPF
    END IF
    SWRITE(*,*) 'Number of saved particles for SFResampleAnalyzeSurfCollis: ',LastAnalyzeSurfCollis%PartNumberSamp
    SDEALLOCATE(LastAnalyzeSurfCollis%WallState)
    SDEALLOCATE(LastAnalyzeSurfCollis%Species)
    ALLOCATE(LastAnalyzeSurfCollis%WallState(6,LastAnalyzeSurfCollis%PartNumberSamp))
    ALLOCATE(LastAnalyzeSurfCollis%Species(LastAnalyzeSurfCollis%PartNumberSamp))
    LastAnalyzeSurfCollis%pushTimeStep = HUGE(LastAnalyzeSurfCollis%pushTimeStep)

    ! Add particle to list
    counter2 = 0
    DO counter = 1, LastAnalyzeSurfCollis%PartNumberSamp
      IF (LastAnalyzeSurfCollis%ReducePartNumber) THEN !reduce saved number of parts (differently for each proc. Could be changed)
        DO !get random (equal!) position between [1,TotalNumberMPF] and accept if .NOT.PartDone
          CALL RANDOM_NUMBER(RandVal)
          counter2 = MIN(1+INT(RandVal*REAL(TotalNumberMPF)),TotalNumberMPF)
          IF (.NOT.PartDone(counter2)) THEN
            PartDone(counter2)=.TRUE.
            EXIT
          END IF
        END DO
      ELSE
        counter2 = counter
      END IF

      iSpec=nSpecies
      DO !determine in which species-"batch" the counter is located (use logical since it is used for ReducePartNumber anyway)
        IF (iSpec.EQ.1) THEN
          IF ( counter2 .GE. 1 ) THEN
            EXIT
          ELSE
            CALL Abort(&
              __STAMP__, &
              'Error in SFResampleAnalyzeSurfCollis. Could not determine iSpec for counter2 ',counter2)
          END IF
        ELSE IF ( counter2 - SUM(AnalyzeSurfCollis%Number(1:iSpec-1)) .GE. 1 ) THEN
          EXIT
        ELSE
          iSpec = iSpec - 1
        END IF
      END DO
      IF (iSpec.GT.1) THEN
        counter2 = counter2 - SUM(AnalyzeSurfCollis%Number(1:iSpec-1))
      END IF
      IF (counter2 .GT. AnalyzeSurfCollis%Number(iSpec)) THEN
        CALL Abort(&
          __STAMP__, &
          'Error in SFResampleAnalyzeSurfCollis. Determined iSpec is wrong for counter2 ',counter2)
      END IF

      LastAnalyzeSurfCollis%WallState(:,counter) = PartSpecData(iSpec,counter2,1:6)
      LastAnalyzeSurfCollis%Species(counter) = iSpec
      LastAnalyzeSurfCollis%pushTimeStep = MIN( LastAnalyzeSurfCollis%pushTimeStep &
        , DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,LastAnalyzeSurfCollis%WallState(4:6,counter)) )
    END DO

    IF (LastAnalyzeSurfCollis%pushTimeStep .LE. 0.) THEN
      CALL Abort(&
        __STAMP__,&
        'Error with SFResampleAnalyzeSurfCollis. Something is wrong with velocities or NormVecOfWall!')
    ELSE
      LastAnalyzeSurfCollis%pushTimeStep = r_SF / LastAnalyzeSurfCollis%pushTimeStep !dt required for smallest projected velo to cross r_SF
      LastAnalyzeSurfCollis%PartNumberDepo = NINT(TotalFlowrateMPF * LastAnalyzeSurfCollis%pushTimeStep)
      SWRITE(*,'(A,E12.5,x,I0)') 'Total Flowrate and to be inserted number of MP for SFResampleAnalyzeSurfCollis: ' &
        ,TotalFlowrateMPF, LastAnalyzeSurfCollis%PartNumberDepo
      IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumberSamp) THEN
        SWRITE(*,*) 'WARNING: PartNumberDepo .GT. PartNumberSamp!'
      END IF
      IF (LastAnalyzeSurfCollis%PartNumberDepo .GT. LastAnalyzeSurfCollis%PartNumThreshold) THEN
        CALL Abort(&
          __STAMP__,&
          'Error with SFResampleAnalyzeSurfCollis: PartNumberDepo .gt. PartNumThreshold',LastAnalyzeSurfCollis%PartNumberDepo)
      END IF
    END IF
  END IF !TotalNumberMPF.GT.0

END SUBROUTINE ReadAnalyzeSurfCollisToHDF5


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
SDEALLOCATE(SurfMesh%SurfSideToGlobSideMap)
!SDALLOCATE(SampWall%Energy)
!SDEALLOCATE(SampWall%Force)
!SDEALLOCATE(SampWall%Counter)
DO iSurfSide=1,SurfMesh%nSides
  SDEALLOCATE(SampWall(iSurfSide)%State)
  SDEALLOCATE(SampWall(iSurfSide)%Adsorption)
  SDEALLOCATE(SampWall(iSurfSide)%Accomodation)
  SDEALLOCATE(SampWall(iSurfSide)%Reaction)
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
SDEALLOCATE(CalcSurfCollis%SpeciesFlags)
SDEALLOCATE(AnalyzeSurfCollis%Data)
SDEALLOCATE(AnalyzeSurfCollis%Spec)
SDEALLOCATE(AnalyzeSurfCollis%BCid)
SDEALLOCATE(AnalyzeSurfCollis%Number)
!SDEALLOCATE(AnalyzeSurfCollis%Rate)
SDEALLOCATE(AnalyzeSurfCollis%BCs)
END SUBROUTINE FinalizeParticleBoundarySampling

END MODULE MOD_Particle_Boundary_Sampling
