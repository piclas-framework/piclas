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

MODULE MOD_Particle_MPI
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE InitParticleMPI
  MODULE PROCEDURE InitParticleMPI
END INTERFACE

#if USE_MPI
INTERFACE IRecvNbOfParticles
  MODULE PROCEDURE IRecvNbOfParticles
END INTERFACE

INTERFACE SendNbOfParticles
  MODULE PROCEDURE SendNbOfParticles
END INTERFACE

INTERFACE FinalizeParticleMPI
  MODULE PROCEDURE FinalizeParticleMPI
END INTERFACE

INTERFACE MPIParticleSend
  MODULE PROCEDURE MPIParticleSend
END INTERFACE

INTERFACE MPIParticleRecv
  MODULE PROCEDURE MPIParticleRecv
END INTERFACE

INTERFACE InitHaloMesh
  MODULE PROCEDURE InitHaloMesh
END INTERFACE

INTERFACE InitParticleCommSize
  MODULE PROCEDURE InitParticleCommSize
END INTERFACE

INTERFACE InitEmissionComm
  MODULE PROCEDURE InitEmissionComm
END INTERFACE

INTERFACE BoxInProc
  MODULE PROCEDURE BoxInProc
END INTERFACE

INTERFACE ExchangeBezierControlPoints3D
  MODULE PROCEDURE ExchangeBezierControlPoints3D
END INTERFACE

INTERFACE AddHaloNodeData
  MODULE PROCEDURE AddHaloNodeData
END INTERFACE

PUBLIC :: InitParticleMPI,FinalizeParticleMPI,InitHaloMesh, InitParticleCommSize, IRecvNbOfParticles, MPIParticleSend
PUBLIC :: MPIParticleRecv
PUBLIC :: InitEmissionComm
PUBLIC :: ExchangeBezierControlPoints3D
PUBLIC :: AddHaloNodeData
#else
PUBLIC :: InitParticleMPI
#endif /*USE_MPI*/

!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                             :: myRealTestValue
INTEGER                         :: color
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI ... '
IF(ParticleMPIInitIsDone) &
  CALL abort(&
    __STAMP__&
  ,' Particle MPI already initialized!')

#if USE_MPI
PartMPI%myRank = myRank
color = 999
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,PartMPI%MyRank,PartMPI%COMM,iERROR)
CALL MPI_COMM_SIZE (PartMPI%COMM,PartMPI%nProcs ,iError)
!PartMPI%COMM   = MPI_COMM_WORLD
IF(PartMPI%nProcs.NE.nProcessors) CALL abort(&
    __STAMP__&
    ,' MPI Communicater-size does not match!', IERROR)
PartCommSize   = 0
IF(PartMPI%MyRank.EQ.0) THEN
  PartMPI%MPIRoot=.TRUE.
ELSE
  PartMPI%MPIRoot=.FALSE.
END IF
iMessage=0
#else
PartMPI%myRank = 0
PartMPI%nProcs = 1
PartMPI%MPIRoot=.TRUE.
#endif  /*USE_MPI*/
!! determine datatype length for variables to be sent
!myRealKind = KIND(myRealTestValue)
!IF (myRealKind.EQ.4) THEN
!  myRealKind = MPI_REAL
!ELSE IF (myRealKind.EQ.8) THEN
!  myRealKind = MPI_DOUBLE_PRECISION
!ELSE
!  myRealKind = MPI_REAL
!END IF

ParticleMPIInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE MPI DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleMPI


#if USE_MPI
SUBROUTINE InitParticleCommSize()
!===================================================================================================================================
! get size of Particle-MPI-Message. Unfortunately, this subroutine have to be called after particle_init because
! all required features have to be read from the ini-File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars
USE MOD_DSMC_Vars,              ONLY:useDSMC, CollisMode, DSMC
USE MOD_Particle_Vars,          ONLY:usevMPF, PDM
USE MOD_Particle_Tracking_vars, ONLY:DoRefMapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: ALLOCSTAT
!===================================================================================================================================

PartCommSize   = 0
! PartState: position and velocity
PartCommSize   = PartCommSize + 6
! Tracking: Include Reference coordinates
IF(DoRefMapping) PartCommSize=PartCommSize+3
! Species-ID
PartCommSize   = PartCommSize + 1
! id of element
PartCommSize   = PartCommSize + 1

IF (useDSMC.AND.(CollisMode.GT.1)) THEN
  IF (usevMPF .AND. DSMC%ElectronicModel) THEN
    ! vib. , rot and electronic energy and macroparticle factor for each particle
    PartCommSize = PartCommSize + 4
  ELSE IF (usevMPF ) THEN
    ! vib. and rot energy and macroparticle factor for each particle
    PartCommSize = PartCommSize + 3
  ELSE IF ( DSMC%ElectronicModel ) THEN
    ! vib., rot. and electronic energy
    PartCommSize = PartCommSize + 3
  ELSE
    ! vib. and rot. energy
    PartCommSize = PartCommSize + 2
  END IF
ELSE
  ! PIC simulation with MPI
  IF (usevMPF) PartCommSize = PartCommSize+1
END IF

! time integration
#if defined(LSERK)
! Pt_tmp for pushing: Runge-Kutta derivative of position and velocity
PartCommSize   = PartCommSize + 6
! IsNewPart for RK-Reconstruction
PartCommSize   = PartCommSize + 1
#endif

! additional stuff for full RK schemes, e.g. implicit and imex RK
#if defined(IMPA) || defined(ROS)
! PartXK
PartCommSize   = PartCommSize + 6
! R_PartXK
PartCommSize   = PartCommSize + 6
! PartQ
PartCommSize   = PartCommSize + 6
! PartDtFrac
PartCommSize   = PartCommSize + 1
! IsNewPart
PartCommSize   = PartCommSize + 1
! FieldAtParticle
PartCommSize   = PartCommSize + 6
! NormVec
PartCommSize   = PartCommSize + 3
! ElementN
PartCommSize   = PartCommSize + 1
! PeriodicMoved
PartCommSize   = PartCommSize + 1
#endif /*IMPA or ROS*/
#if defined(IMPA)
! PartDeltaX
PartCommSize   = PartCommSize + 6
! PartLambdaAccept
PartCommSize   = PartCommSize + 1
! F_PartX0
PartCommSize   = PartCommSize + 6
! F_PartXk
PartCommSize   = PartCommSize + 6
! Norm_F_PartX0, Norm2_F_PartXK, Norm_F_PartXK_old
PartCommSize   = PartCommSize + 3
! DoPartInNewton
PartCommSize   = PartCommSize + 1
! and PartisImplicit
PartCommSize   = PartCommSize + 1
#endif
! if iStage=0, then the PartStateN is not communicated
PartCommSize0  = PartCommSize

ALLOCATE( PartMPIExchange%nPartsSend(2,PartMPI%nMPINeighbors)  &
        , PartMPIExchange%nPartsRecv(2,PartMPI%nMPINeighbors)  &
        , PartRecvBuf(1:PartMPI%nMPINeighbors)                 &
        , PartSendBuf(1:PartMPI%nMPINeighbors)                 &
        , PartMPIExchange%SendRequest(2,PartMPI%nMPINeighbors) &
        , PartMPIExchange%RecvRequest(2,PartMPI%nMPINeighbors) &
        , PartTargetProc(1:PDM%MaxParticleNumber)              &
        , PartMPIDepoSend(1:PDM%MaxParticleNumber)             &
        , STAT=ALLOCSTAT                                       )

IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,' Cannot allocate Particle-MPI-Variables! ALLOCSTAT',ALLOCSTAT)

PartMPIExchange%nPartsSend=0
PartMPIExchange%nPartsRecv=0

IF(DoExternalParts)THEN
  ExtPartCommSize=7
  IF(usevMPF) ExtPartCommSize=8
END IF

END SUBROUTINE InitParticleCommSize


SUBROUTINE IRecvNbOfParticles()
!===================================================================================================================================
! Open Recv-Buffer for number of received particles
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI,PartMPIExchange
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iProc
!===================================================================================================================================

PartMPIExchange%nPartsRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_IRECV( PartMPIExchange%nPartsRecv(:,iProc)                        &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(1,iProc)                       &
                , IERROR )
 ! IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__&
 !         ,' MPI Communication error', IERROR)
END DO ! iProc

END SUBROUTINE IRecvNbOfParticles


SUBROUTINE SendNbOfParticles(doParticle_In)
!===================================================================================================================================
! This routine sends the number of send particles, for which the following steps are performed:
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! The remaining steps are performed in SendParticles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloElemToProc, PartTargetProc
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,PEM,PDM,Species,PartPosRef
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,PartMPIDepoSend
USE MOD_Particle_MPI_Vars,        ONLY:PartShiftVector
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL   :: doParticle_In(1:PDM%ParticleVecLength)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                       :: doPartInExists
INTEGER                       :: iPart,ElemID,iProc
! shape function
INTEGER                       :: CellX,CellY,CellZ!, iPartShape
INTEGER                       :: PartDepoProcs(1:PartMPI%nProcs+1), nDepoProcs, ProcID,LocalProcID
INTEGER                       :: nPartShape
REAL                          :: ShiftedPart(1:3)
LOGICAL                       :: PartInBGM
!===================================================================================================================================
doPartInExists=.FALSE.
IF(PRESENT(DoParticle_IN)) doPartInExists=.TRUE.

! 1) get number of send particles
PartMPIExchange%nPartsSend=0
!ALLOCATE(PartTargetProc(1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
!IF (ALLOCSTAT.NE.0) CALL abort(&
!    __STAMP__ &
!    ' Cannot allocate PartMPIDepoSend!')
PartTargetProc=-1
DO iPart=1,PDM%ParticleVecLength
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ElemID=PEM%Element(iPart)
  IF(ElemID.GT.PP_nElems) THEN
    PartMPIExchange%nPartsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,ElemID))=             &
                                        PartMPIExchange%nPartsSend(1,PartHaloElemToProc(LOCAL_PROC_ID,ElemID))+1
    PartTargetProc(iPart)=PartHaloElemToProc(LOCAL_PROC_ID,ElemID)
  END IF
END DO ! iPart

!DO CellX=GEO%FIBGMimin,GEO%FIBGMimax
!  DO Celly=GEO%FIBGMjmin,GEO%FIBGMjmax
!    DO Cellz=GEO%FIBGMkmin,GEO%FIBGMkmax
!      IF(ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) IPWRITE(*,*) GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs
!    END DO
!  END DO
!END DO


! external particles for communication
IF(DoExternalParts)THEN
  !ALLOCATE(PartMPIDepoSend(1:PDM%ParticleVecLength),STAT=ALLOCSTAT)
  !IF (ALLOCSTAT.NE.0) CALL abort(&
  !    __STAMP__ &
  !    ' Cannot allocate PartMPIDepoSend!')
  PartMPIDepoSend=.FALSE.
  nPartShape=0
  DO iPart=1,PDM%ParticleVecLength
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! Don't deposit neutral external particles!
    IF(.NOT.CHARGEDPARTICLE(iPart)) CYCLE
    ! Don't deposit external shape function particles in cells where local deposition is used (only when DoSFLocalDepoAtBounds=T)
    IF(SkipExternalSFParticles(iPart)) CYCLE
    ! Get indices of background mesh cells
    CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
    CellX = MIN(GEO%FIBGMimax,CellX)
    CellX = MAX(GEO%FIBGMimin,CellX)
    CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
    CellY = MIN(GEO%FIBGMjmax,CellY)
    CellY = MAX(GEO%FIBGMjmin,CellY)
    CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
    CellZ = MIN(GEO%FIBGMkmax,CellZ)
    CellZ = MAX(GEO%FIBGMkmin,CellZ)
    IF(ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
      IF(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1) .GT. 0) THEN
        nPartShape=nPartShape+1
        !PartMPIDepoSend(nPartShape) =  !iPart
        PartMPIDepoSend(iPart) = .TRUE.
      END IF
    END IF
  END DO ! iPart=1,PDM%ParticleVecLength
  ! now, get correct BGM cell for particle
  ! including periodic displacement or BGM element without mpi neighbors
  ! shape-padding could be modified for all other deposition methods? reuse?
  DO iPart=1,PDM%ParticleVecLength
    !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    !IF(PartTargetProc(iPart).EQ.-1) CYCLE
    IF(.NOT.PartMPIDepoSend(iPart)) CYCLE
    IF (Species(PartSpecies(iPart))%ChargeIC.EQ.0) CYCLE ! get BMG cell
    CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
    CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
    CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
    PartInBGM = .TRUE.
    ! check if particle is in range of my FIBGM
    ! first check is outside
    IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
        (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
        (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
      PartInBGM = .FALSE.
    ELSE
      ! particle inside, check if particle is not moved by periodic BC
      ! if this is the case, then the ShapeProcs is not allocated!
      IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
      END IF
    END IF
    IF (.NOT.PartInBGM) THEN
      ! it is possible that the particle has been moved over a periodic side
      IF (GEO%nPeriodicVectors.GT.0) THEN
        ShiftedPart(1:3) = PartState(iPart,1:3) + partShiftVector(1:3,iPart)
        CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
            (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
            (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
          CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellX = MIN(GEO%FIBGMimax,CellX)
          CellX = MAX(GEO%FIBGMimin,CellX)
          CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellY = MIN(GEO%FIBGMjmax,CellY)
          CellY = MAX(GEO%FIBGMjmin,CellY)
          CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          CellZ = MIN(GEO%FIBGMkmax,CellZ)
          CellZ = MAX(GEO%FIBGMkmin,CellZ)
        ELSE
          IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
            IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM! Err2'
            IF(doPartInExists)THEN
              IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))
            ELSE
              IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
            END IF
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
            IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
            IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',ShiftedPart(1),',PartY=',ShiftedPart(2),',PartZ=',&
                    ShiftedPart(3)
            IF(DoRefMapping)THEN
              IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartXi=',PartPosRef(1,iPart)   &
                                                    ,',PartEta=',PartPosRef(2,iPart) &
                                                    ,',PartZeta=',PartPosRef(3,iPart)
            END IF
            CALL Abort(&
            __STAMP__&
            ,'Particle outside BGM! Err2')
          END IF
        END IF
      ELSE
        IPWRITE(UNIT_errOut,*)'Warning in SendNbOfParticles: Particle outside BGM!'
        IF(doPartInExists)THEN
          IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',(PDM%ParticleInside(iPart).AND.DoParticle_In(iPart))
        ELSE
          IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
        END IF
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
        IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',PartState(iPart,1),',PartY=',PartState(iPart,2),',PartZ=',&
                PartState(iPart,3)
        IF(DoRefMapping)THEN
          IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartXi=',PartPosRef(1,iPart),',PartEta=',PartPosRef(2,iPart),',PartZeta=',&
                  PartPosRef(3,iPart)
        END IF
        IPWRITE(UNIT_errOut,*)'Remap particle!'

        CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
        CellX = MIN(GEO%FIBGMimax,CellX)
        CellX = MAX(GEO%FIBGMimin,CellX)
        CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
        CellY = MIN(GEO%FIBGMjmax,CellY)
        CellY = MAX(GEO%FIBGMjmin,CellY)
        CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
        CellZ = MIN(GEO%FIBGMkmax,CellZ)
        CellZ = MAX(GEO%FIBGMkmin,CellZ)
        IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'New-CellX=',CellX,',New-CellY=',CellY,',New-CellZ=',CellZ
        ! nothing to do, because of tolerance, particle could be outside
        !IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
        !    (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
        !    (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN

        !  CALL Abort(&
        !       __STAMP__&
        !      'Particle outside BGM!')
        !END IF
      END IF ! GEO%nPeriodicVectors
    END IF ! PartInBGM
    nDepoProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
    PartDepoProcs(1:nDepoProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nDepoProcs+1)
    DO iProc=1,nDepoProcs
      ProcID=PartDepoProcs(iProc)
      ! particle shall not be send to MyRank, is fixed without MPI communication in MPIParticleSend
      IF(ProcID.EQ.PartMPI%MyRank) CYCLE
      ! if shapeproc is target proc, to net send
      !ElemID=PEM%Element(iPart)
      !IF(PartHaloElemToProc(NATIVE_PROC_ID,ElemID).EQ.ProcID) CYCLE
      ! short version
      LocalProcID=PartMPI%GlobalToLocal(ProcID)
      IF(PartTargetProc(iPart).EQ.LocalProcID) CYCLE
      PartMPIExchange%nPartsSend(2,LocalProcID)= PartMPIExchange%nPartsSend(2,LocalProcID) +1
    END DO ! iProc=1,nDepoProcs
  END DO ! iPart=1,PDM%ParticleVecLength
END IF


! 2) send number of send particles
DO iProc=1,PartMPI%nMPINeighbors
  !IPWRITE(UNIT_stdOut,*) 'Target Number of send particles',PartMPI%MPINeighbor(iProc),PartMPIExchange%nPartsSend(:,iProc)
  CALL MPI_ISEND( PartMPIExchange%nPartsSend(:,iProc)                      &
                , 2                                                          &
                , MPI_INTEGER                                                &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1001                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(1,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)

!  CALL MPI_ISEND( PartMPIExchange%nPartsSend(iProc)                          &
!                , 1                                                          &
!                , MPI_INTEGER                                                &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1001                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%SendRequest(1,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

END SUBROUTINE SendNbOfParticles


SUBROUTINE MPIParticleSend()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! first steps are performed in SendNbOfParticles
! 1) Compute number of Send Particles
! 2) Perform MPI_ISEND with number of particles
! Starting Here:
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
! CAUTION: If particles are sent for deposition, PartTargetProc has the information, if a particle is sent
!          and after the build and wait for number of particles reused to build array with external parts
!          information in PartState,.. can be reused, because they are not overwritten
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartHaloElemToProc, PartCommSize,PartSendBuf, PartRecvBuf &
                                      ,PartTargetProc
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, Species,PartPosRef
#if defined(LSERK)
USE MOD_Particle_Vars,            ONLY:Pt_temp
#endif
USE MOD_DSMC_Vars,                ONLY:useDSMC, CollisMode, DSMC, PartStateIntEn, SpecDSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
USE MOD_LD_Vars,                  ONLY:useLD,PartStateBulkValues
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,PartMPIDepoSend,PartShiftVector, ExtPartCommSize, PartMPIDepoSend
USE MOD_Particle_MPI_Vars,        ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF,  NbrOfExtParticles
#if defined(ROS) || defined(IMPA)
USE MOD_Particle_Vars,            ONLY:PartStateN,PartStage,PartDtFrac,PartQ
USE MOD_Particle_MPI_Vars,        ONLY:PartCommSize0
USE MOD_Timedisc_Vars,            ONLY:iStage
USE MOD_LinearSolver_Vars,        ONLY:PartXK,R_PartXK
USE MOD_Particle_Mesh_Vars,       ONLY:ElemToGlobalElemID
USE MOD_PICInterpolation_Vars,    ONLY:FieldAtParticle
#endif /*ROS or IMPLICIT*/
#if defined(IMPA)
USE MOD_Particle_Vars,            ONLY:F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton &
                                     ,PartDeltaX,PartLambdaAccept
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
#endif /*IMPA*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart,ElemID,iPos,iProc,jPos
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nSendParticles, nSendExtParticles, nRecvExtParticles
INTEGER                       :: ALLOCSTAT
! shape function
INTEGER                       :: CellX,CellY,CellZ!, iPartShape
INTEGER                       :: PartDepoProcs(1:PartMPI%nProcs+1), nDepoProcs, ProcID, jProc, iExtPart, LocalProcID
REAL                          :: ShiftedPart(1:3)
LOGICAL                       :: PartInBGM
#if defined(ROS) || defined(IMPA)
INTEGER                       :: iCounter
#endif /*ROS or IMPA*/
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, MsgRecvLengthPoly
INTEGER                       :: MsgLengthPoly(1:PartMPI%nMPINeighbors), pos_poly(1:PartMPI%nMPINeighbors)
!===================================================================================================================================

#if defined(ROS)
PartCommSize=PartCommSize0+iStage*6
#endif /*ROS*/
#if defined (IMPA)
PartCommSize=PartCommSize0+iStage*6
#endif /*IMPA*/

! ! 1) get number of send particles
! PartMPIExchange%nPartsSend=0
! DO iPart=1,PDM%ParticleVecLength
!   IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
!   ElemID=PEM%Element(iPart)
!   IF(ElemID.GT.PP_nElems)  PartMPIExchange%nPartsSend(PartHaloElemToProc(LOCAL_PROC_ID,ElemID))=             &
!                                         PartMPIExchange%nPartsSend(PartHaloElemToProc(LOCAL_PROC_ID,ElemID))+1
! END DO ! iPart

! external particles add to same message
! IF(DoExternalParts)THEN
!   ALLOCATE( shape_indices(1:PDM%ParticleVecLength), stat=allocstat)
!   shape_indices(:) = 0
!   iPartShape=1
!   DO iPart=1,PDM%ParticleVecLength
!     IF(PDM%ParticleInside(iPart))THEN
!       IF(ALMOSTZERO(Species(PartSpecies(iPart))%ChargeIC) CYCLE        ! Don't deposite neutral particles!
!       CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
!       CellX = MIN(GEO%FIBGMimax,CellX)
!       CellX = MAX(GEO%FIBGMimin,CellX)
!       CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
!       CellY = MIN(GEO%FIBGMkmax,CellY)
!       CellY = MAX(GEO%FIBGMkmin,CellY)
!       CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
!       CellZ = MIN(GEO%FIBGMlmax,CellZ)
!       CellZ = MAX(GEO%FIBGMlmin,CellZ)
!       IF(ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
!         IF(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1) .GT. 0) THEN
!           shape_indices(iPartShape) = i
!           iPartShape = iPartShape + 1
!         END IF
!       END IF
!     END IF ! ParticleInside
!   END DO ! iPart
!   iPartShape = iPartShape - 1
!
! END IF ! DoExternalParts

! ! 2) send number of send particles
! DO iProc=1,PartMPI%nMPINeighbors
!   !IPWRITE(UNIT_stdOut,*) 'Target Number of send particles',PartMPI%MPINeighbor(iProc),PartMPIExchange%nPartsSend(iProc)
!   CALL MPI_ISEND( PartMPIExchange%nPartsSend(iProc)                          &
!                 , 1                                                          &
!                 , MPI_INTEGER                                                &
!                 , PartMPI%MPINeighbor(iProc)                                 &
!                 , 1001                                                       &
!                 , PartMPI%COMM                                               &
!                 , PartMPIExchange%SendRequest(1,iProc)                       &
!                 , IERROR )
!   IF(IERROR.NE.MPI_SUCCESS) CALL abort(__STAMP__&
!           ,' MPI Communication error', IERROR)
!
! !  CALL MPI_ISEND( PartMPIExchange%nPartsSend(iProc)                          &
! !                , 1                                                          &
! !                , MPI_INTEGER                                                &
! !                , PartMPI%MPINeighbor(iProc)                                 &
! !                , 1001                                                       &
! !                , PartMPI%COMM                                               &
! !                , PartMPIExchange%SendRequest(1,PartMPI%MPINeighbor(iProc))  &
! !                , IERROR )
! END DO ! iProc

!DO iProc=1,PartMPI%nMPINeighbors
!  IPWRITE(UNIT_stdOut,*) 'Number of send  particles',  PartMPIExchange%nPartsSend(iProc)&
!                        , 'target proc', PartMPI%MPINeighbor(iProc)
!END DO

!--- Determining the number of additional variables due to VibQuantsPar of polyatomic particles
!--- (size varies depending on the species of particle)
MsgLengthPoly(:) = 0
IF(DSMC%NumPolyatomMolecs.GT.0) THEN
  DO iPart=1,PDM%ParticleVecLength
    IF(PartTargetProc(iPart).EQ.-1) CYCLE
    IF(SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
      MsgLengthPoly(PartTargetProc(iPart)) = MsgLengthPoly(PartTargetProc(iPart)) + PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

! 3) Build Message
DO iProc=1, PartMPI%nMPINeighbors
  ! allocate SendBuf
  nSendParticles=PartMPIExchange%nPartsSend(1,iProc)
  iPos=0
  IF(DoExternalParts)THEN
    nSendExtParticles=PartMPIExchange%nPartsSend(2,iProc)
    IF((nSendExtParticles.EQ.0).AND.(nSendParticles.EQ.0)) CYCLE
    MessageSize=nSendParticles*PartCommSize &
               +nSendExtParticles*ExtPartCommSize
  ELSE
    IF(nSendParticles.EQ.0) CYCLE
    MessageSize=nSendParticles*PartCommSize
  END IF
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    pos_poly(iProc) = MessageSize
    MessageSize = MessageSize + MsgLengthPoly(iProc)
  END IF

  ALLOCATE(PartSendBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'  Cannot allocate PartSendBuf, local ProcId, ALLOCSTAT',iProc,REAL(ALLOCSTAT))
  ! fill message
  DO iPart=1,PDM%ParticleVecLength
    !IF(.NOT.PDM%ParticleInside(iPart)) CYCLE
    !IF(ElemID.GT.PP_nElems) THEN
    !  !IF(PartHaloElemToProc(NATIVE_PROC_ID,ElemID).NE.PartMPI%MPINeighbor(iProc))THEN
    !  IF(PartHaloElemToProc(LOCAL_PROC_ID,ElemID).NE.iProc) CYCLE
    IF(nSendParticles.EQ.0) EXIT
    IF(PartTargetProc(iPart).EQ.iProc) THEN
      !iPos=iPos+1
      ! fill content
      ElemID=PEM%Element(iPart)
      PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(iPart,1:6)
      IF(DoRefMapping) THEN ! + deposition type....
        PartSendBuf(iProc)%content(7+iPos:9+iPos) = PartPosRef(1:3,iPart)
        jPos=iPos+9
      ELSE
        jPos=iPos+6
      END IF
      PartSendBuf(iProc)%content(       1+jPos) = REAL(PartSpecies(iPart),KIND=8)
      jPos=jPos+1
#if defined(LSERK)
      PartSendBuf(iProc)%content(1+jPos:6+jPos) = Pt_temp(iPart,1:6)
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(7+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(7+jPos) = 0.
      END IF
      jPos=jPos+7
#endif
#if defined(ROS) || defined(IMPA)
      ! send iStage - 1 messages
      IF(iStage.GT.0)THEN ! should GT.1
        IF(iStage.EQ.1) CALL Abort(&
               __STAMP__&
         ,' You should never send particles now!')
        PartSendBuf(iProc)%content(1+jpos:6+jpos)        = PartStateN(iPart,1:6)
        DO iCounter=1,iStage-1
          PartSendBuf(iProc)%content(jpos+7+(iCounter-1)*6:jpos+6+iCounter*6) = PartStage(iPart,1:6,iCounter)
        END DO
        jPos=jPos+iStage*6
      ENDIF
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartXK(1:6,iPart)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = R_PartXK(1:6,iPart)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartQ(1:6,iPart)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1) = PartDtFrac(iPart)
      jPos=jPos+1
      IF (PDM%IsNewPart(iPart)) THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.
      END IF
      jPos=jPos+1
      ! fieldatparticle
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = FieldAtParticle(iPart,1:6)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1:jPos+3) = PEM%NormVec(iPart,1:3)
      PEM%NormVec(iPart,1:3)=0.
      jPos=jPos+3
      PartSendBuf(iProc)%content(jPos+1) =REAL(ElemToGlobalElemID(PEM%ElementN(iPart)))
      jPos=jPos+1
      IF (PEM%PeriodicMoved(iPart)) THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.
      END IF
      PEM%PeriodicMoved(iPart)=.FALSE.
      jPos=jPos+1
#endif /*ROS or IMEX */
#if defined(IMPA)
      ! required for particle newton && closed particle description
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = PartDeltaX(1:6,iPart)
      IF (PartLambdaAccept(iPart)) THEN
        PartSendBuf(iProc)%content(7+jPos) = 1.
      ELSE
        PartSendBuf(iProc)%content(7+jPos) = 0.
      END IF
      jPos=jPos+7
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = F_PartX0(1:6,iPart)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1:jPos+6) = F_PartXk(1:6,iPart)
      jPos=jPos+6
      PartSendBuf(iProc)%content(jPos+1)  = Norm_F_PartX0    (iPart)
      PartSendBuf(iProc)%content(jPos+2)  = Norm_F_PartXK    (iPart)
      PartSendBuf(iProc)%content(jPos+3)  = Norm_F_PartXK_old(iPart)
      IF(DoPartInNewton(iPart))THEN
        PartSendBuf(iProc)%content(jPos+4) = 1.0
      ELSE
        PartSendBuf(iProc)%content(jPos+4) = 0.0
      END IF
      jPos=jPos+4
      IF(PartIsImplicit(iPart))THEN
        PartSendBuf(iProc)%content(jPos+1) = 1.0
      ELSE
        PartSendBuf(iProc)%content(jPos+1) = 0.0
      END IF
      jPos=jPos+1
#endif /*IMPA*/
      !PartSendBuf(iProc)%content(       14+jPos) = REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),KIND=8)
      PartSendBuf(iProc)%content(    1+jPos) = REAL(PartHaloElemToProc(NATIVE_ELEM_ID,ElemID),KIND=8)
      jPos=jPos+1
      IF(.NOT.UseLD) THEN
        IF (useDSMC.AND.(CollisMode.GT.1)) THEN
          IF (usevMPF .AND. DSMC%ElectronicModel) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos) = PartMPF(iPart)
            PartSendBuf(iProc)%content(4+jPos) = PartStateIntEn(iPart, 3)
            jPos=jPos+4
          ELSE IF (usevMPF) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos) = PartMPF(iPart)
            jPos=jPos+3
          ELSE IF ( DSMC%ElectronicModel ) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos) = PartStateIntEn(iPart, 3)
            jPos=jPos+3
          ELSE
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            jPos=jPos+2
          END IF
        ELSE
          IF (usevMPF) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartMPF(iPart)
            jPos=jPos+1
          END IF
        END IF
      ELSE ! UseLD == true      =>      useDSMC == true
        IF (CollisMode.GT.1) THEN
          IF (usevMPF .AND. DSMC%ElectronicModel) THEN
            PartSendBuf(iProc)%content( 1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content( 2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content( 3+jPos) = PartMPF(iPart)
            PartSendBuf(iProc)%content( 4+jPos) = PartStateIntEn(iPart, 3)
            PartSendBuf(iProc)%content( 5+jPos:9+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+9
          ELSE IF (usevMPF) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos) = PartMPF(iPart)
            PartSendBuf(iProc)%content(4+jPos:8+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+8
          ELSE IF ( DSMC%ElectronicModel ) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos) = PartStateIntEn(iPart, 3)
            PartSendBuf(iProc)%content(4+jPos:8+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+8
          ELSE
            PartSendBuf(iProc)%content(1+jPos) = PartStateIntEn(iPart, 1)
            PartSendBuf(iProc)%content(2+jPos) = PartStateIntEn(iPart, 2)
            PartSendBuf(iProc)%content(3+jPos:7+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+7
          END IF
        ELSE
          IF (usevMPF) THEN
            PartSendBuf(iProc)%content(1+jPos) = PartMPF(iPart)
            PartSendBuf(iProc)%content(2+jPos:6+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+6
          ELSE
            PartSendBuf(iProc)%content(1+jPos:5+jPos) = PartStateBulkValues(iPart,1:5)
            jPos=jPos+5
          END IF
        END IF
      END IF
      !--- put the polyatomic vibquants per particle at the end of the message
      IF (DSMC%NumPolyatomMolecs.GT.0) THEN
        IF(SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
          PartSendBuf(iProc)%content(pos_poly(iProc)+1:pos_poly(iProc)+PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                                                          = VibQuantsPar(iPart)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
          pos_poly(iProc) = pos_poly(iProc) + PolyatomMolDSMC(iPolyatMole)%VibDOF
        END IF
      END IF
      IF(MOD(jPos,PartCommSize).NE.0) THEN
#if defined(ROS) || defined(IMPA)
        IPWRITE(UNIT_stdOut,*)  'iStage',iStage
        IPWRITE(UNIT_stdOut,*)  'PartCommSize0',PartCommSize0,PartCommSize
#else
        IPWRITE(UNIT_stdOut,*)  'PartCommSize',PartCommSize
#endif
        IPWRITE(UNIT_stdOut,*)  'jPos',jPos
        CALL Abort(&
            __STAMP__&
            ,' Particle-wrong sending message size!')
      END IF
      ! endif timedisc
      ! here iPos because PartCommSize contains DoRefMapping
      iPos=iPos+PartCommSize
      ! particle is ready for send, now it can deleted
      PDM%ParticleInside(iPart) = .FALSE.
#ifdef IMPA
      DoPartInNewton(iPart)   = .FALSE.
      PartLambdaAccept(iPart) = .TRUE.
      PartIsImplicit(iPart)     = .FALSE.
#endif /*IMPA*/
    END IF ! Particle is particle with target proc-id equals local proc id
  END DO  ! iPart
  ! next, external particles has to be handled for deposition
  IF(DoExternalParts)THEN
    IF(nSendExtParticles.EQ.0) CYCLE
    iPos=nSendParticles*PartCommSize
    DO iPart=1,PDM%ParticleVecLength
      IF(PartTargetProc(iPart).EQ.iProc) CYCLE
      IF(.NOT.PartMPIDepoSend(iPart)) CYCLE
      IF (Species(PartSpecies(iPart))%ChargeIC.EQ.0) CYCLE
      ! get BMG cell
      CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      PartInBGM = .TRUE.
      ! check if particle is in range of my FIBGM
      ! first check is outside
      IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
          (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
          (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
        PartInBGM = .FALSE.
      ELSE
        ! particle inside, check if particle is not moved by periodic BC
        ! if this is the case, then the ShapeProcs is not allocated!
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
          PartInBGM = .FALSE.
        END IF
      END IF
      IF (.NOT.PartInBGM) THEN
        ! it is possible that the particle has been moved over a periodic side
        IF (GEO%nPeriodicVectors.GT.0) THEN
          ShiftedPart(1:3) = PartState(iPart,1:3) + partShiftVector(1:3,iPart)
          CellX = INT((ShiftedPart(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellY = INT((ShiftedPart(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellZ = INT((ShiftedPart(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
              (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
              (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
            CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
            CellX = MIN(GEO%FIBGMimax,CellX)
            CellX = MAX(GEO%FIBGMimin,CellX)
            CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
            CellY = MIN(GEO%FIBGMjmax,CellY)
            CellY = MAX(GEO%FIBGMjmin,CellY)
            CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
            CellZ = MIN(GEO%FIBGMkmax,CellZ)
            CellZ = MAX(GEO%FIBGMkmin,CellZ)
          ELSE
            IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) THEN
              IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM! Err2'
              IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
              IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
              IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',ShiftedPart(1),',PartY=',ShiftedPart(2),',PartZ=',&
                      ShiftedPart(3)
              CALL Abort(&
              __STAMP__&
              ,'Particle outside BGM! Err2')
            END IF
          END IF
        ELSE
          IPWRITE(UNIT_errOut,*)'Remap particle!'

          CellX = INT((PartState(iPart,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
          CellX = MIN(GEO%FIBGMimax,CellX)
          CellX = MAX(GEO%FIBGMimin,CellX)
          CellY = INT((PartState(iPart,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
          CellY = MIN(GEO%FIBGMjmax,CellY)
          CellY = MAX(GEO%FIBGMjmin,CellY)
          CellZ = INT((PartState(iPart,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
          CellZ = MIN(GEO%FIBGMkmax,CellZ)
          CellZ = MAX(GEO%FIBGMkmin,CellZ)
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'New-CellX=',CellX,',New-CellY=',CellY,',New-CellZ=',CellZ
          IPWRITE(UNIT_errOut,*)'ERROR in SendNbOfParticles: Particle outside BGM!'
          IPWRITE(UNIT_errOut,*)'iPart =',iPart,',ParticleInside =',PDM%ParticleInside(iPart)
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'minX =',GEO%FIBGMimin,',minY =',GEO%FIBGMjmin,',minZ =',GEO%FIBGMkmin
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'CellX=',CellX,',CellY=',CellY,',CellZ=',CellZ
          IPWRITE(UNIT_errOut,'(I4,3(A,I4))')'maxX =',GEO%FIBGMimax,',maxY =',GEO%FIBGMjmax,',maxZ =',GEO%FIBGMkmax
          IPWRITE(UNIT_errOut,'(I4,3(A,ES13.5))')'PartX=',PartState(iPart,1),',PartY=',PartState(iPart,2),',PartZ=',&
                  PartState(iPart,3)
          !CALL Abort(&
          !     __STAMP__&
          !    'Particle outside BGM!')
        END IF ! GEO%nPeriodicVectors
      END IF ! PartInBGM
      nDepoProcs=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)
      PartDepoProcs(1:nDepoProcs)=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(2:nDepoProcs+1)
      DO jProc=1,nDepoProcs
        !IF(PartMPI%GlobalToLocal(PartDepoProcs(jProc)).EQ.iProc) CYCLE
  !      IPWRITE(*,*) PartDepoProcs(jProc)
        ProcID=PartDepoProcs(jProc)
        LocalProcID=PartMPI%GlobalToLocal(ProcID)
  !      IPWRITE(*,*) 'localprocid,iproc',iProc,localprocid
        IF(PartTargetProc(iPart).EQ.LocalProcID) CYCLE
        IF(LocalProcID.NE.iProc) CYCLE
        PartSendBuf(iProc)%content(1+iPos:6+iPos) = PartState(iPart,1:6)
        PartSendBuf(iProc)%content(       7+iPos) = REAL(PartSpecies(iPart),KIND=8)
        IF (usevMPF) PartSendBuf(iProc)%content( 8+iPos) = PartMPF(iPart)
        ! count only, if particle is sent
        iPos=iPos+ExtPartCommSize
      END DO ! jProc=1,nDepoProcs
    END DO ! iPart=1,PDM%ParticleVecLength
  END IF ! DoExternalParts
  IF(iPos.NE.(MessageSize-MsgLengthPoly(iProc))) IPWRITE(*,*) ' error message size', iPos,(MessageSize-MsgLengthPoly(iProc))
END DO ! iProc

! 4) Finish Received number of particles
DO iProc=1,PartMPI%nMPINeighbors
  CALL MPI_WAIT(PartMPIExchange%SendRequest(1,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(1,iProc),recv_status_list(:,iProc),IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! total number of received particles
PartMPIExchange%nMPIParticles=SUM(PartMPIExchange%nPartsRecv(1,:))
!IPWRITE(UNIT_stdOut,*) 'Number of received particles',SUM(PartMPIExchange%nPartsRecv(1,:))
!IPWRITE(UNIT_stdOut,*) 'Number of received extparticles',SUM(PartMPIExchange%nPartsRecv(2,:))


! caution, fancy trick, particles are sent, but information is not deleted
! temporary storage
IF(DoExternalParts) THEN
  NbrOfExtParticles =SUM(PartMPIExchange%nPartsSend(1,:))+SUM(PartMPIExchange%nPartsRecv(2,:))
  ALLOCATE(ExtPartState  (1:NbrOfExtParticles,1:6) &
          ,ExtPartSpecies(1:NbrOfExtParticles)     &
          !,ExtPartToFIBGM(1:6,1:NbrOfExtParticles) &
          ,STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(&
  __STAMP__&
  ,'  Cannot allocate ExtPartState on Rank',PartMPI%MyRank,REAL(ALLOCSTAT))
  IF (usevMPF) THEN
    ALLOCATE(ExtPartMPF (1:NbrOfExtParticles) &
            ,STAT=ALLOCSTAT)
    IF (ALLOCSTAT.NE.0) CALL abort(&
    __STAMP__&
    ,'  Cannot allocate ExtPartState on Rank',PartMPI%MyRank)
  END IF
  ! map alt state to ext
  iExtPart=0
  DO iPart=1,PDM%ParticleVecLength
    IF(PartTargetProc(iPart).EQ.-1) CYCLE
    iExtPart=iExtPart+1
    ExtPartState(iExtPart,1:6)        = PartState(iPart,1:6)
    ExtPartSpecies(iExtPart)          = PartSpecies(iPart)
    IF (usevMPF) ExtPartMPF(iExtPart) = PartMPF(iPart)
  END DO ! iPart=1,PDM%ParticleVecLength
END IF

DO iPart=1,PDM%ParticleVecLength
  IF(PartTargetProc(iPart).EQ.-1) CYCLE
  PartState(iPart,1:6)=0.
  PartSpecies(iPart)=0
#if defined(LSERK)
  Pt_temp(iPart,1:6)=0.
#endif
END DO ! iPart=1,PDM%ParticleVecLength


! 5) Allocate received buffer and open MPI_IRECV
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE
  nRecvParticles=PartMPIExchange%nPartsRecv(1,iProc)
  MessageSize=nRecvParticles*PartCommSize
  IF(DoExternalParts) THEN
    nRecvExtParticles=PartMPIExchange%nPartsRecv(2,iProc)
    MessageSize=MessageSize   &
               +nRecvExtParticles*ExtPartCommSize
  END IF
  ! determine the maximal possible polyatomic addition to the regular recv message
  IF (DSMC%NumPolyatomMolecs.GT.0) THEN
    MsgRecvLengthPoly = MAXVAL(PolyatomMolDSMC(:)%VibDOF)*nRecvParticles
    MessageSize = MessageSize + MsgRecvLengthPoly
  END IF
  ALLOCATE(PartRecvBuf(iProc)%content(MessageSize),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    IPWRITE(*,*) 'sum of total received particles            ', SUM(PartMPIExchange%nPartsRecv(1,:))
    IPWRITE(*,*) 'sum of total received deposition particles ', SUM(PartMPIExchange%nPartsRecv(2,:))
    CALL abort(&
    __STAMP__&
    ,'  Cannot allocate PartRecvBuf, local source ProcId, Allocstat',iProc,REAL(ALLOCSTAT))
  END IF
  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%RecvRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)

!  CALL MPI_IRECV( PartRecvBuf(iProc)%content                                 &
!                , MessageSize                                                &
!                , MPI_DOUBLE_PRECISION                                       &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1002                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%RecvRequest(2,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

! 6) Send Particles
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE
  nSendParticles=PartMPIExchange%nPartsSend(1,iProc)
  MessageSize=nSendParticles*PartCommSize
  IF(DoExternalParts)THEN
    nSendExtParticles=PartMPIExchange%nPartsSend(2,iProc)
    MessageSize=MessageSize &
               +nSendExtParticles*ExtPartCommSize
  END IF
  IF(DSMC%NumPolyatomMolecs.GT.0) THEN
    MessageSize = MessageSize + MsgLengthPoly(iProc)
  END IF
  !IF(nSendParticles.EQ.0) CYCLE
  !MessageSize=nSendParticles*PartCommSize
  !ALLOCATE(SendBuf(iProc)%content(PartCommSize,nSendParticles))
  CALL MPI_ISEND( PartSendBuf(iProc)%content                                 &
                , MessageSize                                                &
                , MPI_DOUBLE_PRECISION                                       &
                , PartMPI%MPINeighbor(iProc)                                 &
                , 1002                                                       &
                , PartMPI%COMM                                               &
                , PartMPIExchange%SendRequest(2,iProc)                       &
                , IERROR )
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
!  CALL MPI_ISEND( SendBuf(iProc)%content                                     &
!                , MessageSize                                                &
!                , MPI_DOUBLE_PRECISION                                       &
!                , PartMPI%MPINeighbor(iProc)                                 &
!                , 1002                                                       &
!                , PartMPI%COMM                                               &
!                , PartMPIExchange%SendRequest(2,PartMPI%MPINeighbor(iProc))  &
!                , IERROR )
END DO ! iProc

! deallocate send buffer, deallocated to early, first MPI_WAIT??
!DO iProc=1,PartMPI%nMPINeighbors
!  SDEALLOCATE(SendBuf(iProc)%content)
!END DO ! iProc

! and not deallocate, because global, fixed variable
! SDEALLOCATE(PartTargetProc)
! SDEALLOCATE(PartMPIDepoSend)

END SUBROUTINE MPIParticleSend


SUBROUTINE MPIParticleRecv()
!===================================================================================================================================
! this routine sends the particles. Following steps are performed
! 1) Compute number of Send Particles
! 2) Performe MPI_ISEND with number of particles
! 3) Build Message
! 4) MPI_WAIT for number of received particles
! 5) Open Receive-Buffer for particle message -> MPI_IRECV
! 6) Send Particles -> MPI_ISEND
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Tracking_vars,   ONLY:DoRefMapping
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI,PartMPIExchange,PartCommSize, PartRecvBuf,PartSendBuf!,iMessage
USE MOD_Particle_Vars,            ONLY:PartState,PartSpecies,usevMPF,PartMPF,PEM,PDM, PartPosRef
#if defined(LSERK)
USE MOD_Particle_Vars,            ONLY:Pt_temp
#endif
USE MOD_DSMC_Vars,                ONLY:useDSMC, CollisMode, DSMC, PartStateIntEn, SpecDSMC, PolyatomMolDSMC, VibQuantsPar
USE MOD_LD_Vars,                  ONLY:useLD,PartStateBulkValues
! variables for parallel deposition
USE MOD_Particle_MPI_Vars,        ONLY:DoExternalParts,ExtPartCommSize
USE MOD_Particle_MPI_Vars,        ONLY:ExtPartState,ExtPartSpecies,ExtPartMPF
#if defined(ROS) || defined(IMPA)
USE MOD_Particle_Vars,            ONLY:PartStateN,PartStage,PartDtFrac,PartQ
USE MOD_Particle_MPI_Vars,        ONLY:PartCommSize0
USE MOD_Timedisc_Vars,            ONLY:iStage
USE MOD_LinearSolver_Vars,        ONLY:PartXK,R_PartXK
USE MOD_PICInterpolation_Vars,    ONLY:FieldAtParticle
USE MOD_Particle_Mesh_Vars,       ONLY:nTotalElems
USE MOD_Particle_Mesh_Vars,       ONLY:ElemToGlobalElemID
USE MOD_Mesh_Vars,                ONLY:OffSetElem
#endif /*ROS or IMPA*/
#if defined(IMPA)
USE MOD_Particle_Vars,           ONLY:F_PartX0,F_PartXk,Norm_F_PartX0,Norm_F_PartXK,Norm_F_PartXK_old,DoPartInNewton &
                                     ,PartDeltaX,PartLambdaAccept
USE MOD_Particle_Vars,           ONLY:PartIsImplicit
#endif /*IMPA*/
USE MOD_DSMC_Vars,               ONLY: RadialWeighting
USE MOD_DSMC_Symmetry2D,         ONLY: DSMC_2D_RadialWeighting
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iProc, iPos, nRecv, PartID,jPos, iPart, TempNextFreePosition
INTEGER                       :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINeighbors)
INTEGER                       :: MessageSize, nRecvParticles, nRecvExtParticles
!INTEGER,ALLOCATABLE           :: RecvArray(:,:), RecvArray_glob(:,:,:)
!CHARACTER(LEN=64)             :: filename,hilf
! shape function
INTEGER                       :: iExtPart
#if defined(ROS) || defined(IMPA)
INTEGER                       :: iCounter, LocElemID,iElem
#endif /*ROS or IMPA*/
! Polyatomic Molecules
INTEGER                       :: iPolyatMole, pos_poly, MsgLengthPoly
! INTEGER , ALLOCATABLE         :: MsgLengthPoly(:)
!===================================================================================================================================

!IPWRITE(UNIT_stdOut,*) 'exchange',PartMPIExchange%nMPIParticles

!DO iProc=1,PartMPI%nMPINeighbors
!  IPWRITE(UNIT_stdOut,*) 'Number of received  particles',  PartMPIExchange%nPartsRecv(iProc) &
!                        , 'source proc', PartMPI%MPINeighbor(iProc)
!END DO

#if defined(ROS)
PartCommSize=PartCommSize0+iStage*6
#endif /*ROS*/
#if defined (IMPA)
PartCommSize=PartCommSize0+iStage*6
#endif /*MPA*/
! IF (DSMC%NumPolyatomMolecs.GT.0) ALLOCATE(MsgLengthPoly(1:PartMPI%nMPINeighbors))

DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsSend(:,iProc)).EQ.0) CYCLE
  CALL MPI_WAIT(PartMPIExchange%SendRequest(2,iProc),MPIStatus,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
    __STAMP__&
    ,' MPI Communication error', IERROR)
END DO ! iProc

! old number of already filled ExtParticles
IF(DoExternalParts) iExtPart=SUM(PartMPIExchange%nPartsSend(1,:))

nRecv=0
DO iProc=1,PartMPI%nMPINeighbors
  IF(SUM(PartMPIExchange%nPartsRecv(:,iProc)).EQ.0) CYCLE
  nRecvParticles=PartMPIExchange%nPartsRecv(1,iProc)
  IF(DoExternalParts) THEN
    nRecvExtParticles=PartMPIExchange%nPartsRecv(2,iProc)
  END IF
  ! determine the maximal possible polyatomic addition to the regular message
  IF (DSMC%NumPolyatomMolecs.GT.0) THEN
    MsgLengthPoly = MAXVAL(PolyatomMolDSMC(:)%VibDOF)*nRecvParticles
  ELSE
    MsgLengthPoly = 0
  END IF
  !IF(nRecvParticles.EQ.0) CYCLE
  MessageSize=nRecvParticles*PartCommSize
  pos_poly = MessageSize
  IF (DSMC%NumPolyatomMolecs.GT.0) MessageSize = MessageSize + MsgLengthPoly
  ! finish communication with iproc
  CALL MPI_WAIT(PartMPIExchange%RecvRequest(2,iProc),recv_status_list(:,iProc),IERROR)
  ! correct loop shape
  ! DO iPart=1,nRecvParticles
  ! nParts 1 Pos=1..17
  ! nPart2 2 Pos=1..17,18..34
  DO iPos=0,MessageSize-1-MsgLengthPoly,PartCommSize
    IF(nRecvParticles.EQ.0) EXIT
    nRecv=nRecv+1
    PartID = PDM%nextFreePosition(nRecv+PDM%CurrentNextFreePosition)
    IF(PartID.EQ.0)  CALL abort(&
      __STAMP__&
      ,' Error in ParticleExchange_parallel. Corrupted list: PIC%nextFreePosition', nRecv)
    PartState(PartID,1:6)   = PartRecvBuf(iProc)%content( 1+iPos: 6+iPos)
    IF(DoRefMapping)THEN
      PartPosRef(1:3,PartID) = PartRecvBuf(iProc)%content(7+iPos: 9+iPos)
      jPos=iPos+9
    ELSE
      jPos=iPos+6
    END IF
    PartSpecies(PartID)     = INT(PartRecvBuf(iProc)%content( 1+jPos),KIND=4)
    jPos=jPos+1
#if defined(LSERK)
    Pt_temp(PartID,1:6)     = PartRecvBuf(iProc)%content( 1+jPos:6+jPos)
    IF ( INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL Abort(&
        __STAMP__&
        ,'Error with IsNewPart in MPIParticleRecv!')
    END IF
    jPos=jPos+7
#endif
#if defined(ROS) || defined(IMPA)
    IF(iStage.GT.0)THEN
      IF(iStage.EQ.1) CALL Abort(&
             __STAMP__&
       ,' You should never receive particle now!')
      PartStateN(PartID,1:6)     = PartRecvBuf(iProc)%content(jpos+1:jpos+6)
      DO iCounter=1,iStage-1
        PartStage(PartID,1:6,iCounter) = PartRecvBuf(iProc)%content(jpos+7+(iCounter-1)*6:jpos+6+iCounter*6)
      END DO
      jPos=jPos+iStage*6
    END IF
    PartXK(1:6,PartID)         = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    R_PartXK(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    PartQ(1:6,PartID)          = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    PartDtFrac(PartID) = PartRecvBuf(iProc)%content(jPos+1)
    jPos=jPos+1
    IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PDM%IsNewPart(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PDM%IsNewPart(PartID)=.FALSE.
    ELSE
      CALL Abort(&
        __STAMP__&
        ,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content( 1+jPos))
    END IF
    jPos=jPos+1
    ! fieldatparticle
    FieldAtParticle(PartID,1:6)  = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    PEM%NormVec(PartID,1:3)  = PartRecvBuf(iProc)%content(jPos+1:jPos+3)
    jPos=jPos+3
    LocElemID=INT(PartRecvBuf(iProc)%content(jPos+1),KIND=4)
    IF((LocElemID-OffSetElem.GE.1).AND.(LocElemID-OffSetElem.LE.PP_nElems))THEN
      PEM%ElementN(PartID)=LocElemID-OffSetElem
    ELSE
      PEM%ElementN(PartID)=0
      DO iElem=PP_nElems+1,nTotalElems
        IF(ElemToGlobalElemID(iElem).EQ.LocElemID)THEN
          PEM%ElementN(PartID)=iElem
          EXIT
        END IF
      END DO ! iElem=1,nTotalElems
      IF(PEM%ElementN(PartID).EQ.0)THEN
        CALL Abort(&
          __STAMP__&
          ,'Error with IsNewPart in MPIParticleRecv: PEM%ElementN(PartID) = 0!')
      END IF
    END IF
    jPos=jPos+1
    IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 1) THEN
      PEM%PeriodicMoved(PartID)=.TRUE.
    ELSE IF ( INT(PartRecvBuf(iProc)%content( 1+jPos)) .EQ. 0) THEN
      PEM%PeriodicMoved(PartID)=.FALSE.
    ELSE
      CALL Abort(&
        __STAMP__&
        ,'Error with IsNewPart in MPIParticleRecv!',1,PartRecvBuf(iProc)%content( 1+jPos))
    END IF
    jPos=jPos+1
#endif /*ROS or IMPA*/
#if defined(IMPA)
    PartDeltaX(1:6,PartID)     = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    IF ( INT(PartRecvBuf(iProc)%content( 7+jPos)) .EQ. 1) THEN
      PartLambdaAccept(PartID)=.TRUE.
    ELSE ! IF ( INT(PartRecvBuf(iProc)%content( 14+jPos)) .EQ. 0) THEN
      PartLambdaAccept(PartID)=.FALSE.
    END IF
    jPos=jPos+7
    F_PartX0(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    F_PartXk(1:6,PartID)       = PartRecvBuf(iProc)%content(jPos+1:jPos+6)
    jPos=jPos+6
    Norm_F_PartX0    (PartID) = PartRecvBuf(iProc)%content(jPos+1)
    Norm_F_PartXk    (PartID) = PartRecvBuf(iProc)%content(jPos+2)
    Norm_F_PartXk_Old(PartID) = PartRecvBuf(iProc)%content(jPos+3)
    IF(PartRecvBuf(iProc)%content(jPos+4).EQ.1.0)THEN
      DoPartInNewton(PartID) = .TRUE.
    ELSE
      DoPartInNewton(PartID) = .FALSE.
    END IF
    jPos=jPos+4
    IF(PartRecvBuf(iProc)%content(jPos+1).EQ.1.0)THEN
        PartIsImplicit(PartID) = .TRUE.
    ELSE
        PartIsImplicit(PartID) = .FALSE.
    END IF
    jPos=jPos+1
#endif /*IMPA*/
    PEM%Element(PartID)     = INT(PartRecvBuf(iProc)%content(1+jPos),KIND=4)
    jPos=jPos+1
    IF(.NOT.UseLD) THEN
      IF (useDSMC.AND.(CollisMode.GT.1)) THEN
        IF (usevMPF .AND. DSMC%ElectronicModel) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(2+jPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(3+jPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(4+jPos)
          jPos=jPos+4
        ELSE IF ( usevMPF ) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(2+jPos)
          PartMPF(PartID)           = PartRecvBuf(iProc)%content(3+jPos)
          jPos=jPos+3
        ELSE IF ( DSMC%ElectronicModel) THEN
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(2+jPos)
          PartStateIntEn(PartID, 3) = PartRecvBuf(iProc)%content(3+jPos)
          jPos=jPos+3
        ELSE
          PartStateIntEn(PartID, 1) = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2) = PartRecvBuf(iProc)%content(2+jPos)
          jPos=jPos+2
        END IF
      ELSE
        IF (usevMPF)THEN
          PartMPF(PartID) = PartRecvBuf(iProc)%content(1+jPos)
          jPos=jPos+1
        END IF
      END IF
    ELSE ! UseLD == true      =>      useDSMC == true
      IF (CollisMode.GT.1) THEN
        IF (usevMPF .AND. DSMC%ElectronicModel) THEN
          PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(2+jPos)
          PartMPF(PartID)                 = PartRecvBuf(iProc)%content(3+jPos)
          PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(4+jPos)
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(5+jPos:9+jPos)
          jPos=jPos+9
        ELSE IF ( usevMPF) THEN
          PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(2+jPos)
          PartMPF(PartID)                 = PartRecvBuf(iProc)%content(3+jPos)
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(4+jPos:8+jPos)
          jPos=jPos+8
        ELSE IF ( DSMC%ElectronicModel ) THEN
          PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(2+jPos)
          PartStateIntEn(PartID, 3)       = PartRecvBuf(iProc)%content(3+jPos)
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(4+jPos:8+jPos)
          jPos=jPos+8
        ELSE
          PartStateIntEn(PartID, 1)       = PartRecvBuf(iProc)%content(1+jPos)
          PartStateIntEn(PartID, 2)       = PartRecvBuf(iProc)%content(2+jPos)
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(3+jPos:7+jPos)
          jPos=jPos+7
        END IF
      ELSE
        IF (usevMPF) THEN
          PartMPF(PartID)                 = PartRecvBuf(iProc)%content(1+jPos)
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(2+jPos:6+jPos)
          jPos=jPos+6
        ELSE
          PartStateBulkValues(PartID,1:5) = PartRecvBuf(iProc)%content(1+jPos:5+jPos)
          jPos=jPos+5
        END IF
      END IF
    END IF
    IF(MOD(jPos,PartCommSize).NE.0)THEN
#if defined(ROS) || defined(IMPA)
      IPWRITE(UNIT_stdOut,*)  'iStage',iStage
      IPWRITE(UNIT_stdOut,*)  'PartCommSize0',PartCommSize0,PartCommSize
#else
      IPWRITE(UNIT_stdOut,*)  'jPos',jPos
#endif
      CALL Abort(&
          __STAMP__&
          ,' Particle-wrong receiving message size!')
    END IF
    IF (DSMC%NumPolyatomMolecs.GT.0) THEN
      IF(SpecDSMC(PartSpecies(PartID))%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(PartSpecies(PartID))%SpecToPolyArray
        IF(ALLOCATED(VibQuantsPar(PartID)%Quants)) DEALLOCATE(VibQuantsPar(PartID)%Quants)
        ALLOCATE(VibQuantsPar(PartID)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        VibQuantsPar(PartID)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                            = NINT(PartRecvBuf(iProc)%content(pos_poly+1:pos_poly+PolyatomMolDSMC(iPolyatMole)%VibDOF))
        pos_poly = pos_poly + PolyatomMolDSMC(iPolyatMole)%VibDOF
      END IF
    END IF
    ! Set Flag for received parts in order to localize them later
    PDM%ParticleInside(PartID) = .TRUE.
#if defined(IMPA) || defined(ROS)
    ! only for fully implicit
    !IF(PartIsImplicit(PartID))THEN
    !    ! get LastElement in local system
    !    ! element position is verified in armijo if it is needed. This prevents a
    !    ! circular definition.
    !    PEM%LastElement(PartID)=-1
    ! END IF
!   ! PEM%StageElement(PartID)= PEM%Element(PartID)
!    StagePartPos(PartID,1)  = PartState(PartID,1)
!    StagePartPos(PartID,2)  = PartState(PartID,2)
!    StagePartPos(PartID,3)  = PartState(PartID,3)
#else
    PEM%lastElement(PartID) = -888
#endif
  END DO
  IF(DoExternalParts)THEN
    jPos=MessageSize
    IF(nRecvExtParticles.EQ.0) CYCLE
    MessageSize=nRecvExtParticles*ExtPartCommSize+jPos
    DO iPos=jPos,MessageSize-1,ExtPartCommSize
      iExtPart=iExtPart+1
      ExtPartState(iExtPart,1:6) = PartRecvBuf(iProc)%content( 1+iPos: 6+iPos)
      ExtPartSpecies(iExtPart)   = INT(PartRecvBuf(iProc)%content( 7+iPos),KIND=4)
      IF (usevMPF) ExtPartMPF(iExtPart) = PartRecvBuf(iProc)%content( 8+iPos)
    END DO ! iPos
  END IF ! DoExternalParts
  ! be nice: deallocate the receive buffer
  ! deallocate non used array
END DO ! iProc

TempNextFreePosition        = PDM%CurrentNextFreePosition
PDM%ParticleVecLength       = PDM%ParticleVecLength + PartMPIExchange%nMPIParticles
PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + PartMPIExchange%nMPIParticles
IF(PDM%ParticleVecLength.GT.PDM%MaxParticleNumber) CALL abort(&
    __STAMP__&
    ,' ParticleVecLegnth>MaxParticleNumber due to MPI-communication!')

IF(RadialWeighting%DoRadialWeighting) THEN
  ! Checking whether received particles have to be cloned or deleted
  DO iPart = 1,nrecv
    PartID = PDM%nextFreePosition(iPart+TempNextFreePosition)
    CALL DSMC_2D_RadialWeighting(PartID,PEM%Element(PartID))
  END DO
END IF

! validate solution and check
! debug
! ProcID=PartMPI%MyRank
! iMessage=iMessage+1
! ALLOCATE(RecvArray(1:2,0:PartMPI%nProcs-1))
! IF(PartMPI%MPIROOT)THEN
!   ALLOCATE(RecvArray_glob(1:2,0:PartMPI%nProcs-1,0:PartMPI%nProcs-1))
! ELSE
!   ALLOCATE(RecvArray_glob(1:1,1:1,1:1))
! END IF
! RecvArray=0
! DO iProc=1,PartMPI%nMPINeighbors
!   RecvArray(1,PartMPI%MPINeighbor(iProc)) = PartMPIExchange%nPartsSend(iProc)
!   RecvArray(2,PartMPI%MPINeighbor(iProc)) = PartMPIExchange%nPartsRecv(iProc)
! END DO ! iProc
!
! ! mpi gather
! CALL MPI_GATHER(RecvArray,2*PartMPI%nProcs,MPI_INTEGER,RecvArray_glob,PartMPI%nProcs*2,MPI_INTEGER,0,PartMPI%COMM,iError)
! IF(PartMPI%MPIROOT) THEN
!   WRITE( hilf, '(I5.5)') iMessage
!   filename = 'particle_send_'//TRIM(hilf)//'.csv'
!   OPEN(unit=63,FILE=filename,status='UNKNOWN')
!   !  WRITE(63,'(A6)',ADVANCE='NO') ' Procs'
!   !  DO iProc=0,PartMPI%nProcs-1
!   !    WRITE(63, '(A3,I5)',ADVANCE='NO') ',  ',iProc
!   !  END DO ! iProc
!   !  WRITE(63, '(A)') ''
!   DO iProc=0,PartMPI%nProcs-1
!   !    WRITE(63,'(I5)',ADVANCE='NO') iProc
!      DO ProcID=0,PartMPI%nProcs-1
!        WRITE(63,'(3x,I5)',ADVANCE='NO') RecvArray_glob(1,ProcID,iProc)
!      END DO ! ProcID
!      WRITE(63,'(A)') ''
!    END DO ! iProc
!   CLOSE(63)
!
!   filename = 'particle_recv_'//TRIM(hilf)//'.csv'
!   OPEN(unit=63,FILE=filename,status='UNKNOWN')
! !  WRITE(63,'(A6)',ADVANCE='NO') ' Procs'
! !  WRITE(63,'(A6)',ADVANCE='NO') ' Procs'
! !  DO iProc=0,PartMPI%nProcs-1
! !    WRITE(63, '(A3,I5)',ADVANCE='NO') ',  ',iProc
! !  END DO ! iProc
! !  WRITE(63, '(A)') ''
!   DO iProc=0,PartMPI%nProcs-1
! !    WRITE(63,'(I5)',ADVANCE='NO') iProc
!     DO ProcID=0,PartMPI%nProcs-1
!       WRITE(63,'(3x,I5)',ADVANCE='NO') RecvArray_glob(2,ProcID,iProc)
!     END DO ! ProcID
!     WRITE(63,'(A)') ''
!   END DO ! iProc
!
!   CLOSE(63)
! END IF


! deallocate send,receive buffer
DO iProc=1,PartMPI%nMPINeighbors
  SDEALLOCATE(PartRecvBuf(iProc)%content)
  SDEALLOCATE(PartSendBuf(iProc)%content)
END DO ! iProc


! final
PartMPIExchange%nPartsRecv=0
PartMPIExchange%nPartsSend=0


END SUBROUTINE MPIParticleRecv


!==================================================================================================================================
!> Check whether a particle is inside a cell where a local deposition method is used instead of the shape function and
!> returns true or false (when true, the particle is not communicated via MPI to the neighboring ranks)
!==================================================================================================================================
PURE FUNCTION SkipExternalSFParticles(PartID)
! MODULES
USE MOD_PICDepo_Vars       ,ONLY: DoSFLocalDepoAtBounds
USE MOD_Particle_Vars      ,ONLY: PEM
USE MOD_Particle_Mesh_Vars ,ONLY: IsLocalDepositionBCElem
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: PartID                  !< Particle ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL             :: SkipExternalSFParticles !< returns true is the particle must not be sent to other ranks via MPI
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Only check when shape functions are used in combination with local deposition
IF(.NOT.DoSFLocalDepoAtBounds)THEN
  SkipExternalSFParticles=.FALSE.
  RETURN
END IF

! Check if particle is inside a local deposition element
IF(IsLocalDepositionBCElem(PEM%Element(PartID)))THEN
  SkipExternalSFParticles=.TRUE.
ELSE
  SkipExternalSFParticles=.FALSE.
END IF

END FUNCTION SkipExternalSFParticles


SUBROUTINE FinalizeParticleMPI()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_MPI_Vars
USE MOD_Particle_Vars,            ONLY:Species,nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: nInitRegions,iInitRegions,iSpec
!===================================================================================================================================

nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits+(1-Species(iSpec)%StartnumberOfInits)
END DO ! iSpec
IF(nInitRegions.GT.0) THEN
  DO iInitRegions=1,nInitRegions
    IF(PartMPI%InitGroup(iInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(PartMPI%InitGroup(iInitRegions)%Comm,iERROR)
    END IF
  END DO ! iInitRegions
END IF
CALL MPI_COMM_FREE(PartMPI%COMM,iERROR)


SDEALLOCATE( PartHaloElemToProc)
SDEALLOCATE( PartHaloNodeToProc)
SDEALLOCATE( PartMPI%isMPINeighbor)
SDEALLOCATE( PartMPI%MPINeighbor )
SDEALLOCATE( PartMPI%GlobalToLocal )
SDEALLOCATE( PartMPIExchange%nPartsSend)
SDEALLOCATE( PartMPIExchange%nPartsRecv)
SDEALLOCATE( PartMPIExchange%RecvRequest)
SDEALLOCATE( PartMPIExchange%SendRequest)
SDEALLOCATE( PartMPIExchange%Send_message)
SDEALLOCATE( PartMPI%isMPINodeNeighbor)
SDEALLOCATE( PartMPI%MPINodeNeighbor)
SDEALLOCATE( PartMPI%InitGroup)
SDEALLOCATE( PartSendBuf)
SDEALLOCATE( PartRecvBuf)
SDEALLOCATE( NodeSendBuf)
SDEALLOCATE( NodeRecvBuf)
SDEALLOCATE( NodeExchange%nNodesSend)
SDEALLOCATE( NodeExchange%nNodesRecv)
SDEALLOCATE( NodeExchange%RecvRequest)
SDEALLOCATE( NodeExchange%SendRequest)
SDEALLOCATE( ExtPartState)
SDEALLOCATE( ExtPartSpecies)

! and for communication
SDEALLOCATE( PartTargetProc )
SDEALLOCATE( PartMPIDepoSend )

ParticleMPIInitIsDone=.FALSE.
END SUBROUTINE FinalizeParticleMPI


SUBROUTINE ExchangeBezierControlPoints3D()
!===================================================================================================================================
! exchange all beziercontrolpoints at MPI interfaces
! maybe extended to periodic sides, to be tested
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_Mesh_Vars,                  ONLY:NGeo,NGeoElevated,nSides,firstMPISide_YOUR,MortarSlave2MasterInfo&
                                        ,firstMortarMPISide,lastMortarMPISide,lastMPISide_YOUR,SideToElem
USE MOD_Particle_Surfaces,          ONLY:GetSideSlabNormalsAndIntervals
USE MOD_Particle_Surfaces_vars,     ONLY:BezierControlPoints3D,SideSlabIntervals,BezierControlPoints3DElevated &
                                        ,SideSlabIntervals,SideSlabNormals,BoundingBoxIsEmpty

!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 ::BezierSideSize,SendID, iSide,SideID,ElemID
!===================================================================================================================================

! funny: should not be required, as sides are built for master and slave sides??
! communicate the MPI Master Sides to Slaves
! all processes have now filled sides and can compute the particles inside the proc region

SendID=1
BezierSideSize=3*(NGeo+1)*(NGeo+1)
DO iNbProc=1,nNbProcs
  ! Start receive face data
  IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
    nRecVal     =BezierSideSize*nMPISides_rec(iNbProc,SendID)
    SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
    CALL MPI_IRECV(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nRecVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest_Flux(iNbProc),iError)
  END IF
  ! Start send face data
  IF(nMPISides_send(iNbProc,SendID).GT.0)THEN
    nSendVal    =BezierSideSize*nMPISides_send(iNbProc,SendID)
    SideID_start=OffsetMPISides_send(iNbProc-1,SendID)+1
    SideID_end  =OffsetMPISides_send(iNbProc,SendID)
    CALL MPI_ISEND(BezierControlPoints3D(:,:,:,SideID_start:SideID_end),nSendVal,MPI_DOUBLE_PRECISION,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest_Flux(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs

DO iNbProc=1,nNbProcs
  IF(nMPISides_rec(iNbProc,SendID).GT.0) CALL MPI_WAIT(RecRequest_Flux(iNbProc) ,MPIStatus,iError)
  IF(iERROR.NE.0) CALL abort(&
  __STAMP__&
  ,' MPI-Error during BezierControlPoint-exchange. iError', iERROR)
END DO !iProc=1,nNBProcs
! Check send operations
DO iNbProc=1,nNbProcs
  IF(nMPISides_send(iNbProc,SendID).GT.0) CALL MPI_WAIT(SendRequest_Flux(iNbProc),MPIStatus,iError)
  IF(iERROR.NE.0) CALL abort(&
  __STAMP__&
  ,' MPI-Error during BezierControlPoint-exchange. iError', iERROR)
END DO !iProc=1,nNBProcs

! build the bounding box for YOUR-MPI-sides without mortar sides
DO iSide=firstMPISide_YOUR,lastMPISide_YOUR
  ! elevation occurs within this routine
  CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                         &
                                     ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide) &
                                     ,SideSlabNormals(1:3,1:3,iSide)                                         &
                                     ,SideSlabInterVals(1:6,iSide)                                           &
                                     ,BoundingBoxIsEmpty(iSide)                                              )
END DO

! build the bounding box for missing MPI-mortar sides, or YOUR mortar sides
! actually, I do not know, if this is requried
DO iSide=firstMortarMPISide,lastMortarMPISide
  ElemID=SideToElem(S2E_ELEM_ID,iSide)
  SideID=MortarSlave2MasterInfo(iSide)
  IF(ElemID.NE.-1) CYCLE
  IF(SideID.EQ.-1) CYCLE
  ! elevation occurs within this routine
  CALL GetSideSlabNormalsAndIntervals(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,iSide)                         &
                                     ,BezierControlPoints3DElevated(1:3,0:NGeoElevated,0:NGeoElevated,iSide) &
                                     ,SideSlabNormals(1:3,1:3,iSide)                                         &
                                     ,SideSlabInterVals(1:6,iSide)                                           &
                                     ,BoundingBoxIsEmpty(iSide)                                              )
END DO

DO iSide=1,nSides
  !IF(MortarSlave2MasterInfo(iSide).NE.-1)CYCLE
  IF(SUM(ABS(SideSlabIntervals(:,iSide))).EQ.0)THEN
    CALL abort(&
    __STAMP__&
    ,'  Zero bounding box found!, iSide',iSide)
  END IF
END DO

END SUBROUTINE ExchangeBezierControlPoints3D


SUBROUTINE InitHaloMesh()
!===================================================================================================================================
! communicate all direct neighbor sides from master to slave
! has to be called after GetSideType and MPI_INIT of DG solver
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars
USE MOD_PreProc
!USE MOD_Mesh_Vars,                  ONLY:nSides
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_vars,     ONLY:DoRefMapping
USE MOD_Particle_Mesh_Vars,         ONLY:nTotalBCSides,nTotalSides
#endif /*CODE_ANALYZE*/
USE MOD_Particle_MPI_Vars,          ONLY:PartMPI,PartHaloElemToProc,printMPINeighborWarnings
USE MOD_Particle_MPI_Halo,          ONLY:IdentifyHaloMPINeighborhood,ExchangeHaloGeometry
USE MOD_Particle_Mesh_Vars,         ONLY:nPartSides,nTotalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 ::iElem
INTEGER                 ::iProc,ALLOCSTAT,iMPINeighbor
LOGICAL                 ::TmpNeigh
INTEGER,ALLOCATABLE     ::SideIndex(:),ElemIndex(:)
!===================================================================================================================================

ALLOCATE(SideIndex(1:nPartSides),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate SideIndex!')
SideIndex=0
ALLOCATE(ElemIndex(1:PP_nElems),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL abort(&
__STAMP__&
,'  Cannot allocate ElemIndex!')
ElemIndex=0

! check epsilondistance
DO iProc=0,PartMPI%nProcs-1
  IF(iProc.EQ.PartMPI%MyRank) CYCLE
  LOGWRITE(*,*)'  - Identify non-immediate MPI-Neighborhood...'
  !--- AS: identifies which of my node have to be sent to iProc w.r.t. to
  !        eps vicinity region.
  CALL IdentifyHaloMPINeighborhood(iProc,SideIndex,ElemIndex)
  LOGWRITE(*,*)'    ...Done'

  LOGWRITE(*,*)'  - Exchange Geometry of MPI-Neighborhood...'
  CALL ExchangeHaloGeometry(iProc,ElemIndex)
  LOGWRITE(*,*)'    ...Done'
  SideIndex(:)=0
  ElemIndex(:)=0
END DO
DEALLOCATE(SideIndex,STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) THEN
  CALL abort(&
__STAMP__&
,'Could not deallocate SideIndex')
END IF

#ifdef CODE_ANALYZE
IF(DoRefMapping) CALL CheckArrays(nTotalSides,nTotalElems,nTotalBCSides)
#endif /*CODE_ANALYZE*/

! Make sure PMPIVAR%MPINeighbor is consistent
DO iProc=0,PartMPI%nProcs-1
  IF (PartMPI%MyRank.EQ.iProc) CYCLE
  IF (PartMPI%MyRank.LT.iProc) THEN
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
  ELSE IF (PartMPI%MyRank.GT.iProc) THEN
    CALL MPI_RECV(TmpNeigh,1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
    CALL MPI_SEND(PartMPI%isMPINeighbor(iProc),1,MPI_LOGICAL,iProc,1101,PartMPI%COMM,IERROR)
  END IF
  !IPWRITE(UNIT_stdOut,*) 'check',tmpneigh,PartMPI%isMPINeighbor(iProc)
  IF (TmpNeigh.NEQV.PartMPI%isMPINeighbor(iProc)) THEN
    IF(printMPINeighborWarnings)THEN
      WRITE(*,*) 'WARNING: MPINeighbor set to TRUE',PartMPI%MyRank,iProc
    END IF
    IF(.NOT.PartMPI%isMPINeighbor(iProc))THEN
      PartMPI%isMPINeighbor(iProc) = .TRUE.
      PartMPI%nMPINeighbors=PartMPI%nMPINeighbors+1
    END IF
  END IF
END DO


! fill list with neighbor proc id and add local neighbor id to PartHaloElemToProc
ALLOCATE( PartMPI%MPINeighbor(PartMPI%nMPINeighbors) &
        , PartMPI%GlobalToLocal(0:PartMPI%nProcs-1)  )
iMPINeighbor=0
PartMPI%GlobalToLocal=-1
DO iProc=0,PartMPI%nProcs-1
  IF(PartMPI%isMPINeighbor(iProc))THEN
    iMPINeighbor=iMPINeighbor+1
    PartMPI%MPINeighbor(iMPINeighbor)=iProc
    PartMPI%GlobalToLocal(iProc)     =iMPINeighbor
    DO iElem=PP_nElems+1,nTotalElems
      IF(iProc.EQ.PartHaloElemToProc(NATIVE_PROC_ID,iElem)) PartHaloElemToProc(LOCAL_PROC_ID,iElem)=iMPINeighbor
    END DO ! iElem
  END IF
END DO

IF(iMPINeighbor.NE.PartMPI%nMPINeighbors) CALL abort(&
  __STAMP__&
  , ' Found number of mpi neighbors does not match! ', iMPINeighbor,REAL(PartMPI%nMPINeighbors))

IF(PartMPI%nMPINeighbors.GT.0)THEN
  IF(ANY(PartHaloElemToProc(LOCAL_PROC_ID,:).EQ.-1)) IPWRITE(UNIT_stdOut,*) ' Local proc id not found'
  IF(MAXVAL(PartHaloElemToProc(LOCAL_PROC_ID,:)).GT.PartMPI%nMPINeighbors) IPWRITE(UNIT_stdOut,*) ' Local proc id too high.'
  IF(MINVAL(PartHaloElemToProc(NATIVE_ELEM_ID,:)).LT.1) IPWRITE(UNIT_stdOut,*) ' native elem id too low'
  IF(MINVAL(PartHaloElemToProc(NATIVE_PROC_ID,:)).LT.0) IPWRITE(UNIT_stdOut,*) ' native proc id not found'
  IF(MAXVAL(PartHaloElemToProc(NATIVE_PROC_ID,:)).GT.PartMPI%nProcs-1) IPWRITE(UNIT_stdOut,*) ' native proc id too high.'
END IF

END SUBROUTINE InitHaloMesh


SUBROUTINE InitEmissionComm()
!===================================================================================================================================
! build emission communicators for particle emission regions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartMPI
USE MOD_Particle_Vars,          ONLY:Species,nSpecies
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
#ifndef PP_HDG
USE MOD_CalcTimeStep,           ONLY:CalcTimeStep
#endif /*PP_HDG*/
USE MOD_Particle_MPI_Vars,      ONLY:halo_eps
!USE MOD_Particle_Mesh,          ONLY:BoxInProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec,iInit,iNode,iRank
INTEGER                         :: nInitRegions
LOGICAL                         :: RegionOnProc
REAL                            :: xCoords(3,8),lineVector(3),radius,height
REAL                            :: xlen,ylen,zlen
REAL                            :: dt
INTEGER                         :: color,iProc
INTEGER                         :: noInitRank,InitRank
!INTEGER,ALLOCATABLE             :: DummyRank(:)
LOGICAL                         :: hasRegion
!===================================================================================================================================

! get number of total init regions
nInitRegions=0
DO iSpec=1,nSpecies
  nInitRegions=nInitRegions+Species(iSpec)%NumberOfInits+(1-Species(iSpec)%StartnumberOfInits)
END DO ! iSpec
IF(nInitRegions.EQ.0) RETURN

! allocate communicators
ALLOCATE( PartMPI%InitGroup(1:nInitRegions))
!ALLOCATE( DummyRank(0:PartMPI%nProcs-1) )
!DO iRank=0,PartMPI%nProcs-1
!  DummyRank(iRank)=iRank
!END DO ! iRank

nInitRegions=0
DO iSpec=1,nSpecies
  RegionOnProc=.FALSE.
  DO iInit=Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
    nInitRegions=nInitRegions+1
    SELECT CASE(TRIM(Species(iSpec)%Init(iInit)%SpaceIC))
    CASE ('point')
       xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
       RegionOnProc=PointInProc(xCoords(1:3,1))
    CASE ('line_with_equidistant_distribution')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE ('line')
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      RegionOnProc=BoxInProc(xCoords(1:3,1:2),2)
    CASE('disc')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('gyrotron_circle')
      Radius=Species(iSpec)%Init(iInit)%RadiusIC+Species(iSpec)%Init(iInit)%RadiusICGyro
      !xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      !ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
      !     SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      xlen=Radius
      ylen=Radius
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      IF(Species(iSpec)%Init(iInit)%initialParticleNumber.NE.0)THEN
        lineVector(1:3)=(/0.,0.,Species(iSpec)%Init(iInit)%CuboidHeightIC/)
      ELSE
#ifndef PP_HDG
        dt = CALCTIMESTEP()
#endif /*PP_HDG*/
        lineVector(1:3)= dt* Species(iSpec)%Init(iInit)%VeloIC/Species(iSpec)%Init(iInit)%alpha
        zlen=0.
      END IF
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+lineVector+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('circle_equidistant')
      xlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(1)*Species(iSpec)%Init(iInit)%NormalIC(1))
      ylen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(2)*Species(iSpec)%Init(iInit)%NormalIC(2))
      zlen=Species(iSpec)%Init(iInit)%RadiusIC * &
           SQRT(1.0 - Species(iSpec)%Init(iInit)%NormalIC(3)*Species(iSpec)%Init(iInit)%NormalIC(3))
      ! all 8 edges
      xCoords(1:3,1) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,-zlen/)
      xCoords(1:3,2) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,-zlen/)
      xCoords(1:3,3) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,-zlen/)
      xCoords(1:3,4) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,-zlen/)
      xCoords(1:3,5) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,-ylen,+zlen/)
      xCoords(1:3,6) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,-ylen,+zlen/)
      xCoords(1:3,7) = Species(iSpec)%Init(iInit)%BasePointIC+(/-xlen,+ylen,+zlen/)
      xCoords(1:3,8) = Species(iSpec)%Init(iInit)%BasePointIC+(/+xlen,+ylen,+zlen/)
      RegionOnProc=BoxInProc(xCoords(1:3,1:8),8)
    CASE('cuboid')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC&
                                                           +Species(iSpec)%Init(iInit)%BaseVector2IC

      IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        height = halo_eps
      ELSE
        height= Species(iSpec)%Init(iInit)%CuboidHeightIC
      END IF
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('sphere')
      ASSOCIATE ( radius => Species(iSpec)%Init(iInit)%RadiusIC        ,&
                  origin => Species(iSpec)%Init(iInit)%BasePointIC(1:3) )
        ! Set the 8 bounding box coordinates depending on the origin and radius
        xCoords(1:3,1)=origin + (/ radius  , -radius , -radius/)
        xCoords(1:3,2)=origin + (/ radius  , radius  , -radius/)
        xCoords(1:3,3)=origin + (/ -radius , radius  , -radius/)
        xCoords(1:3,4)=origin + (/ -radius , -radius , -radius/)
        xCoords(1:3,5)=origin + (/ radius  , -radius , radius /)
        xCoords(1:3,6)=origin + (/ radius  , radius  , radius /)
        xCoords(1:3,7)=origin + (/ -radius , radius  , radius /)
        xCoords(1:3,8)=origin + (/ -radius , -radius , radius /)
      END ASSOCIATE
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      radius = Species(iSpec)%Init(iInit)%RadiusIC
      ! here no radius, already inclueded
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC-Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                           -Species(iSpec)%Init(iInit)%BaseVector2IC

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*Species(iSpec)%Init(iInit)%BaseVector1IC&
                                   +2.0*Species(iSpec)%Init(iInit)%BaseVector2IC

      IF (Species(iSpec)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        height = halo_eps
      ELSE
        height= Species(iSpec)%Init(iInit)%CylinderHeightIC
      END IF
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cuboid_vpi')

      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC
      xCoords(1:3,2)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=Species(iSpec)%Init(iInit)%BasePointIC+Species(iSpec)%Init(iInit)%BaseVector1IC&
                                                           +Species(iSpec)%Init(iInit)%BaseVector2IC

      height = halo_eps
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)
    CASE('cylinder_vpi')
      lineVector(1) = Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(3) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(iSpec)%Init(iInit)%BaseVector1IC(3) * Species(iSpec)%Init(iInit)%BaseVector2IC(1) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(iSpec)%Init(iInit)%BaseVector1IC(1) * Species(iSpec)%Init(iInit)%BaseVector2IC(2) - &
        Species(iSpec)%Init(iInit)%BaseVector1IC(2) * Species(iSpec)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
         CALL abort(&
         __STAMP__&
         ,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      radius = Species(iSpec)%Init(iInit)%RadiusIC

      xCoords(1:3,1)=Species(iSpec)%Init(iInit)%BasePointIC-radius*Species(iSpec)%Init(iInit)%BaseVector1IC &
                                                           -radius*Species(iSpec)%Init(iInit)%BaseVector2IC

      xCoords(1:3,2)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector1IC
      xCoords(1:3,3)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector2IC
      xCoords(1:3,4)=xCoords(1:3,1)+2.0*radius*Species(iSpec)%Init(iInit)%BaseVector1IC&
                                   +2.0*radius*Species(iSpec)%Init(iInit)%BaseVector2IC

      height = halo_eps
      DO iNode=1,4
        xCoords(1:3,iNode+4)=xCoords(1:3,iNode)+lineVector*height
      END DO ! iNode
      RegionOnProc=BoxInProc(xCoords,8)


    CASE('LD_insert')
      RegionOnProc=.TRUE.
    CASE('cell_local')
      RegionOnProc=.TRUE.
    CASE('cuboid_equal')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL abort(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)

     !~j CALL abort(&
     !~j __STAMP__&
     !~j ,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
    CASE ('cuboid_with_equidistant_distribution')
       xlen = SQRT(Species(iSpec)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(iSpec)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(iSpec)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(iSpec)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(iSpec)%Init(iInit)%BaseVector1IC(1)).OR. &
           (ylen.NE.Species(iSpec)%Init(iInit)%BaseVector2IC(2)).OR. &
           (zlen.NE.Species(iSpec)%Init(iInit)%CuboidHeightIC)) THEN
          CALL abort(&
          __STAMP__&
          ,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
       END IF
       DO iNode=1,8
        xCoords(1:3,iNode) = Species(iSpec)%Init(iInit)%BasePointIC(1:3)
       END DO
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE('sin_deviation')
       IF(Species(iSpec)%Init(iInit)%initialParticleNumber.NE. &
            (Species(iSpec)%Init(iInit)%maxParticleNumberX * Species(iSpec)%Init(iInit)%maxParticleNumberY &
            * Species(iSpec)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',iSpec,' does not match number of particles in each direction!'
         CALL abort(&
         __STAMP__&
         ,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       xCoords(1:3,1) = (/GEO%xminglob,GEO%yminglob,GEO%zminglob/)
       xCoords(1:3,2) = xCoords(1:3,1) + (/xlen,0.,0./)
       xCoords(1:3,3) = xCoords(1:3,1) + (/0.,ylen,0./)
       xCoords(1:3,4) = xCoords(1:3,1) + (/xlen,ylen,0./)
       xCoords(1:3,5) = xCoords(1:3,1) + (/0.,0.,zlen/)
       xCoords(1:3,6) = xCoords(1:3,5) + (/xlen,0.,0./)
       xCoords(1:3,7) = xCoords(1:3,5) + (/0.,ylen,0./)
       xCoords(1:3,8) = xCoords(1:3,5) + (/xlen,ylen,0./)
       RegionOnProc=BoxInProc(xCoords,8)
    CASE ('IMD')
       RegionOnProc=.TRUE.
    CASE DEFAULT
      IPWRITE(*,*) 'ERROR: Species ', iSpec, 'of', iInit, 'is using an unknown SpaceIC!'
      CALL abort(&
      __STAMP__&
      ,'ERROR: Given SpaceIC is not implemented!')
    END SELECT
    ! create new communicator
    color=MPI_UNDEFINED
    !IPWRITE(UNIT_stdOut,*) RegionOnProc
    IF(RegionOnProc) color=nInitRegions!+1
    ! set communicator id
    Species(iSpec)%Init(iInit)%InitComm=nInitRegions

    ! create ranks for RP communicator
    IF(PartMPI%MPIRoot) THEN
      InitRank=-1
      noInitRank=-1
      !myRPRank=0
      iRank=0
      PartMPI%InitGroup(nInitRegions)%MyRank=0
      IF(RegionOnProc) THEN
        InitRank=0
      ELSE
        noInitRank=0
      END IF
      DO iProc=1,PartMPI%nProcs-1
        CALL MPI_RECV(hasRegion,1,MPI_LOGICAL,iProc,0,PartMPI%COMM,MPIstatus,iError)
        IF(hasRegion) THEN
          InitRank=InitRank+1
          CALL MPI_SEND(InitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        ELSE
          noInitRank=noInitRank+1
          CALL MPI_SEND(noInitRank,1,MPI_INTEGER,iProc,0,PartMPI%COMM,iError)
        END IF
      END DO
    ELSE
      CALL MPI_SEND(RegionOnProc,1,MPI_LOGICAL,0,0,PartMPI%COMM,iError)
      CALL MPI_RECV(PartMPI%InitGroup(nInitRegions)%MyRank,1,MPI_INTEGER,0,0,PartMPI%COMM,MPIstatus,iError)
    END IF

    ! create new emission communicator
!    IPWRITE(UNIT_stdOut,*) 'color',color
!    SWRITE(*,*) 'comm null',MPI_COMM_NULL
    CALL MPI_COMM_SPLIT(PartMPI%COMM, color, PartMPI%InitGroup(nInitRegions)%MyRank, PartMPI%InitGroup(nInitRegions)%COMM,iError)
    IF(RegionOnProc) CALL MPI_COMM_SIZE(PartMPI%InitGroup(nInitRegions)%COMM,PartMPI%InitGroup(nInitRegions)%nProcs ,iError)
!    IPWRITE(UNIT_stdOut,*) 'loc rank',PartMPI%InitGroup(nInitRegions)%MyRank
!    IPWRITE(UNIT_stdOut,*) 'loc comm',PartMPI%InitGroup(nInitRegions)%COMM
    IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0 .AND. RegionOnProc) &
    WRITE(UNIT_StdOut,*) ' Emission-Region,Emission-Communicator:',nInitRegions,PartMPI%InitGroup(nInitRegions)%nProcs,' procs'
    IF(PartMPI%InitGroup(nInitRegions)%COMM.NE.MPI_COMM_NULL) THEN
      IF(PartMPI%InitGroup(nInitRegions)%MyRank.EQ.0) THEN
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.TRUE.
      ELSE
        PartMPI%InitGroup(nInitRegions)%MPIRoot=.FALSE.
      END IF
      !ALLOCATE(PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1))
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%GroupToComm(PartMPI%InitGroup(nInitRegions)%MyRank)= PartMPI%MyRank
      CALL MPI_ALLGATHER(PartMPI%MyRank,1,MPI_INTEGER&
                        ,PartMPI%InitGroup(nInitRegions)%GroupToComm(0:PartMPI%InitGroup(nInitRegions)%nProcs-1)&
                       ,1,MPI_INTEGER,PartMPI%InitGroup(nInitRegions)%COMM,iERROR)
      ALLOCATE(PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1))
      PartMPI%InitGroup(nInitRegions)%CommToGroup(0:PartMPI%nProcs-1)=-1
      DO iRank=0,PartMPI%InitGroup(nInitRegions)%nProcs-1
        PartMPI%InitGroup(nInitRegions)%CommToGroup(PartMPI%InitGroup(nInitRegions)%GroupToComm(iRank))=iRank
      END DO ! iRank
      !IPWRITE(*,*) 'CommToGroup',PartMPI%InitGroup(nInitRegions)%GroupToComm
    END IF
  END DO ! iniT
END DO ! iSpec

!DEALLOCATE(DummyRank)

END SUBROUTINE InitEmissionComm


FUNCTION BoxInProc(CartNodes,nNodes)
!===================================================================================================================================
! check if bounding box is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNodes(1:3,1:nNodes)
INTEGER,INTENT(IN):: nNodes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: BoxInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

BoxInProc=.FALSE.
! get background of nodes
xmin=HUGE(1)
xmax=-HUGE(1)
ymin=HUGE(1)
ymax=-HUGE(1)
zmin=HUGE(1)
zmax=-HUGE(1)
testval = CEILING((MINVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((MAXVAL(CartNodes(1,:))-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((MINVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((MAXVAL(CartNodes(2,:))-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((MINVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((MAXVAL(CartNodes(3,:))-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) BoxInProc=.TRUE.

END FUNCTION BoxInProc


FUNCTION PointInProc(CartNode)
!===================================================================================================================================
! check if point is on proc
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars,       ONLY:GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: CartNode(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL           :: PointInProc
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: xmin,xmax,ymin,ymax,zmin,zmax,testval
!===================================================================================================================================

PointInProc=.FALSE.
! get background of nodes
xmin=HUGE(1)
xmax=-HUGE(1)
ymin=HUGE(1)
ymax=-HUGE(1)
zmin=HUGE(1)
zmax=-HUGE(1)
testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmin    = MIN(xmin,testval)
testval = CEILING((CartNode(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
xmax    = MAX(xmax,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymin    = MIN(ymin,testval)
testval = CEILING((CartNode(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
ymax    = MAX(ymax,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmin    = MIN(zmin,testval)
testval = CEILING((CartNode(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
zmax    = MAX(zmax,testval)

IF(    ((xmin.LE.GEO%FIBGMimax).AND.(xmax.GE.GEO%FIBGMimin)) &
  .AND.((ymin.LE.GEO%FIBGMjmax).AND.(ymax.GE.GEO%FIBGMjmin)) &
  .AND.((zmin.LE.GEO%FIBGMkmax).AND.(zmax.GE.GEO%FIBGMkmin)) ) PointInProc=.TRUE.

END FUNCTION PointInProc


SUBROUTINE CheckArrays(nTotalSides,nTotalElems,nTotalBCSides)
!===================================================================================================================================
! check if any entry of the checked arrays exists and if the entry is not NAN
! instead of using IEEE standard, the infamous nan-check a(i).NE.a(i) is used
! Sanity check for refmapping and mpi-communication.
! PO: not sure if it is required any more.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_MPI_Vars,      ONLY:PartHaloElemToProc
USE MOD_Mesh_Vars,              ONLY:BC,nGeo,XCL_NGeo,DXCL_NGEO
USE MOD_Particle_Mesh_Vars,     ONLY:SidePeriodicType,PartBCSideList
USE MOD_Particle_Mesh_Vars,     ONLY:PartElemToSide,PartSideToElem!,PartElemToElemGlob
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
USE MOD_Particle_Surfaces_Vars, ONLY:SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                 :: nTotalSides,nTotalElems
INTEGER,INTENT(IN),OPTIONAL        :: nTotalBCSides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iElem,iVar,iVar2,i,j,k
INTEGER                            :: ilocSide,iSide
!===================================================================================================================================

DO iElem=1,nTotalElems
  ! PartElemToSide
  DO ilocSide=1,6
    DO iVar=1,2
      IF(PartElemToSide(iVar,ilocSide,iElem).NE.PartElemToSide(iVar,ilocSide,iElem)) CALL abort(&
__STAMP__&
, ' Error in PartElemToSide')
    END DO ! iVar=1,2
  END DO ! ilocSide=1,6
  IF(DoRefMapping)THEN
    ! XCL_NGeo & dXCL_NGeo
    DO k=0,NGeo
      DO j=0,NGeo
        DO i=0,NGeo
          DO iVar=1,3
            IF(XCL_NGeo(iVar,i,j,k,iElem).NE.XCL_NGeo(iVar,i,j,k,iElem)) CALL abort(&
__STAMP__&
, ' Error in XCL_NGeo')
            DO iVar2=1,3
              IF(dXCL_NGeo(iVar2,iVar,i,j,k,iElem).NE.dXCL_NGeo(iVar2,iVar,i,j,k,iElem)) CALL abort(&
__STAMP__&
, ' Error in dXCL_NGeo')
            END DO ! iVar2=1,3
          END DO ! iVar=1,3
        END DO ! i=0,NGeo
      END DO ! j=0,NGeo
    END DO ! k=0,NGeo
  END IF ! DoRefMapping
!  ! PartElemToElem
!  DO ilocSide=1,6
!    IF(PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem).NE.PartElemToElem(E2E_NB_ELEM_ID,ilocSide,iElem)) CALL abort(&
!       __STAMP__&
!       , ' Error in PartElemToElem')
!    IF(PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem).NE.PartElemToElem(E2E_NB_LOC_SIDE_ID,ilocSide,iElem)) CALL abort(&
!       __STAMP__&
!       , ' Error in PartElemToElem')
!  END DO ! ilocSide=1,6
END DO ! iElem=1,nTotalElems
IF(DoRefMapping)THEN
  ! PartBCSideList
  DO iSide = 1,nTotalSides
    IF(PartBCSideList(iSide).NE.PartBCSideList(iSide)) CALL abort(&
__STAMP__&
        , ' Error in dXCL_NGeo')
  END DO ! iSide=1,nTotalSides
  ! BezierControlPoints3D
  DO iSide=1,nTotalBCSides
    DO k=0,NGeo
      DO j=0,NGeo
        DO iVar=1,3
          IF(BezierControlPoints3D(iVar,j,k,iSide) &
         .NE.BezierControlPoints3D(iVar,j,k,iSide)) CALL abort(&
__STAMP__&
, ' Error in dXCL_NGeo')
        END DO ! iVar=1,3
      END DO ! j=0,nGeo
    END DO ! k=0,nGeo
    ! Slabnormals & SideSlabIntervals
    DO iVar=1,3
      DO iVar2=1,3
        IF(SideSlabNormals(iVar2,iVar,iSide).NE.SideSlabNormals(iVar2,iVar,iSide)) CALL abort(&
__STAMP__&
, ' Error in PartHaloElemToProc')
      END DO ! iVar2=1,PP_nVar
    END DO ! iVar=1,PP_nVar
    DO ilocSide=1,6
      IF(SideSlabIntervals(ilocSide,iSide).NE.SideSlabIntervals(ilocSide,iSide)) CALL abort(&
__STAMP__&
, ' Error in SlabInvervalls')
    END DO ! ilocSide=1,6
    IF(BoundingBoxIsEmpty(iSide).NEQV.BoundingBoxIsEmpty(iSide)) CALL abort(&
__STAMP__&
, ' Error in BoundingBoxIsEmpty')
  END DO ! iSide=1,nTotalBCSides
END IF ! DoRefMapping
! PartHaloElemToProc
DO iElem=PP_nElems+1,nTotalElems
  DO iVar=1,3
    IF(PartHaloElemToProc(iVar,iElem).NE.PartHaloElemToProc(iVar,iElem)) CALL abort(&
__STAMP__&
, ' Error in PartHaloElemToProc')
  END DO ! iVar=1,3
END DO ! iElem=PP_nElems+1,nTotalElems
DO iSide = 1,nTotalSides
  DO iVar=1,5
    IF(PartSideToElem(iVar,iSide).NE.PartSideToElem(iVar,iSide)) CALL abort(&
        __STAMP__&
        , ' Error in PartSideToElem')
  END DO ! iVar=1,5
  IF(SidePeriodicType(iSide).NE.SidePeriodicType(iSide)) CALL abort(&
      __STAMP__&
      , ' Error in BCSideType')
  IF(BC(iSide).NE.BC(iSide)) CALL abort(&
      __STAMP__&
      , ' Error in BC')
END DO ! iSide=1,nTotalSides


END SUBROUTINE CheckArrays


SUBROUTINE AddHaloNodeData(DataInReal)
!===================================================================================================================================
!> Add the cell node data of halo nodes to local nodes at same position and send local node data to halo procs
!> Input Array of proc local nodes (REAL)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars         ,ONLY: nNodes
USE MOD_Particle_MPI_Vars ,ONLY: PartMPI
USE MOD_Particle_MPI_Vars ,ONLY: NodeSendBuf, NodeRecvBuf, NodeExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(INOUT) :: DataInReal(1:nNodes)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iProc, iPos, iSendNode, iRecvNode
INTEGER :: recv_status_list(1:MPI_STATUS_SIZE,1:PartMPI%nMPINodeNeighbors)
!===================================================================================================================================

! open receive buffer
DO iProc=1,PartMPI%nMPINodeNeighbors
  IF(NodeExchange%nNodesRecv(iProc).EQ.0) CYCLE
  CALL MPI_IRECV( NodeRecvBuf(iProc)%content                &
                , NodeExchange%nNodesRecv(iProc)            &
                , MPI_DOUBLE_PRECISION                      &
                , PartMPI%MPINodeNeighbor(iProc)%COMMProcID &
                , 1415                                      &
                , PartMPI%COMM                              &
                , NodeExchange%RecvRequest(iProc)           &
                , IERROR )
END DO ! iProc

! build message
! after this message, the receiving process knows to which of his nodes it receives and the sending process will know which nodes to
! send
DO iProc=1,PartMPI%nMPINodeNeighbors
  IF(NodeExchange%nNodesSend(iProc).EQ.0) CYCLE
  iPos=1
  ! zero send content buffers
  NodeSendBuf(iProc)%content(:)= 0.
  ! fill send content buffers
  DO iSendNode=1,NodeExchange%nNodesSend(iProc)
    NodeSendBuf(iProc)%content(iPos)=DataInReal(PartMPI%MPINodeNeighbor(iProc)%SendList(iSendNode))
    iPos=iPos+1
  END DO ! iSendNode=1,NodeExchange%nNodesSend(iProc)
END DO

DO iProc=1,PartMPI%nMPINodeNeighbors
  IF(NodeExchange%nNodesSend(iProc).EQ.0) CYCLE
  CALL MPI_ISEND( NodeSendBuf(iProc)%content                &
                , NodeExchange%nNodesSend(iProc)            &
                , MPI_DOUBLE_PRECISION                      &
                , PartMPI%MPINodeNeighbor(iProc)%COMMProcID &
                , 1415                                      &
                , PartMPI%COMM                              &
                , NodeExchange%SendRequest(iProc)           &
                , IERROR )
END DO ! iProc

! 4) Finish Received indexing of received nodes
DO iProc=1,PartMPI%nMPINodeNeighbors
  IF(NodeExchange%nNodesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(NodeExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error (Node data)', IERROR)
  END IF
  IF(NodeExchange%nNodesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(NodeExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error (Node data)', IERROR)
  END IF
END DO ! iProc

! fill list with received Node-IDs
DO iProc=1,PartMPI%nMPINodeNeighbors
  IF(NodeExchange%nNodesRecv(iProc).EQ.0) CYCLE
  iPos=1
  ! fill data-array with data from recv content buffers
  DO iRecvNode=1,NodeExchange%nNodesRecv(iProc)
    DataInReal(PartMPI%MPINodeNeighbor(iProc)%RecvList(iRecvNode))=DataInReal(PartMPI%MPINodeNeighbor(iProc)%RecvList(iRecvNode)) &
        + NodeRecvBuf(iProc)%content(iPos)
    iPos=iPos+1
  END DO ! RecvNode=1,NodeExchange%nNodesRecv(iProc)
  ! zero recv content buffers
  NodeRecvBuf(iProc)%content(:)= 0.
END DO ! iProc

END SUBROUTINE AddHaloNodeData
#endif /*USE_MPI*/

END MODULE MOD_Particle_MPI
