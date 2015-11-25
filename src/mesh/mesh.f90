#include "boltzplatz.h"

MODULE MOD_Mesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

INTERFACE Flip_S2M
  MODULE PROCEDURE Flip_S2M
END INTERFACE

INTERFACE Flip_M2S
  MODULE PROCEDURE Flip_M2S
END INTERFACE

INTERFACE CGNS_SideToVol
  MODULE PROCEDURE CGNS_SideToVol
END INTERFACE

INTERFACE CGNS_VolToSide
  MODULE PROCEDURE CGNS_VolToSide
END INTERFACE

INTERFACE SideToVol
  MODULE PROCEDURE SideToVol
END INTERFACE

INTERFACE VolToSide
  MODULE PROCEDURE VolToSide
END INTERFACE

INTERFACE VolToVol
  MODULE PROCEDURE VolToVol
END INTERFACE

INTERFACE ElemToNBElem
  MODULE PROCEDURE ElemToNBElem
END INTERFACE

PUBLIC::InitMesh
PUBLIC::Flip_S2M
PUBLIC::Flip_M2S
PUBLIC::CGNS_SideToVol
PUBLIC::CGNS_VolToSide
PUBLIC::SideToVol
PUBLIC::VolToSide
PUBLIC::VolToVol
PUBLIC::ElemToNBElem
PUBLIC::FinalizeMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars, ONLY:xGP,InterpolationInitIsDone
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_Metrics,            ONLY:CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Analyze_Vars,       ONLY:CalcPoyntingInt
#ifdef PARTICLES
USE MOD_Particle_Mesh,          ONLY:InitParticleMesh,InitElemVolumes ! new
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,BezierControlPoints3DElevated,SideSlabNormals,SideSlabIntervals
USE MOD_Particle_Surfaces_Vars, ONLY:BoundingBoxIsEmpty,BezierElevation,ElemSlabNormals,ElemSlabIntervals
USE MOD_Particle_Tracking_Vars, ONLY:DoRefMapping
#endif
#ifdef MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
USE MOD_MPI_Vars,           ONLY:offsetSurfElemMPI
#endif
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE  :: NodeCoords(:,:,:,:,:)
REAL              :: x(3),PI,meshScale
INTEGER           :: iElem,i,j,k,f,s,p,q, Flip(3), ijk(3),pq(3)
LOGICAL           :: debugmesh
INTEGER           :: iSide,countSurfElem,iProc
INTEGER,ALLOCATABLE :: countSurfElemMPI(:)
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL abort(__STAMP__,'InitMesh not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

! prepare pointer structure (get nElems, etc.)
MeshFile = GETSTR('MeshFile')

useCurveds=GETLOGICAL('useCurveds','.TRUE.')
#ifdef MPI
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.)
#else
CALL OpenDataFile(MeshFile,create=.FALSE.)
#endif
CALL ReadAttribute(File_ID,'Ngeo',1,IntegerScalar=NGeo)
SWRITE(UNIT_stdOut,'(A67,I2.0)') ' |                           NGeo |                                ', NGeo

CALL CloseDataFile()
IF(useCurveds)THEN
  IF(PP_N.LT.NGeo)THEN
    SWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                             & can cause problems on periodic boundaries! Set N>NGeo'
  END IF
ELSE
  IF(NGeo.GT.1) THEN
    SWRITE(*,*) ' WARNING: Using linear elements although NGeo>1!'
  END IF
END IF

CALL readMesh(MeshFile,NodeCoords) !set nElems

!schmutz fink
PP_nElems=nElems

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
CALL setLocalSideIDs()

ALLOCATE(XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems))
XCL_NGeo = 0.

#ifdef MPI
! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
CALL exchangeFlip()
#endif

! fill ElemToSide, SideToElem,BC
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(SideToElem2(4,2*nInnerSides+nBCSides+nMPISides))
ALLOCATE(BC(1:nSides))
!ALLOCATE(AnalyzeSide(1:nSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
SideToElem2 = -1   !mapping side to elem, sorted by elem ID (for ProlongToFace) 
BC          = 0
!AnalyzeSide = 0

!lower and upper index of U/gradUx/y/z _plus
!lower and upper index of U/gradUx/y/z _plus
sideID_minus_lower = 1
sideID_minus_upper = nBCSides+nInnerSides+nMPISides_MINE
sideID_plus_lower  = nBCSides+1
sideID_plus_upper  = nBCSides+nInnerSides+nMPISides

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
CALL fillMeshInfo()

#ifdef PARTICLES
! save geometry information for particle tracking
CALL InitParticleMesh()
!CALL InitElemVolumes()
#endif

! calculating offset of surface elements for DSMC surface output

#ifdef MPI

IF(ALLOCATED(offsetSurfElemMPI))DEALLOCATE(offsetSurfElemMPI)
ALLOCATE(offsetSurfElemMPI(0:nProcessors))
offsetSurfElemMPI=0

countSurfElem=0

DO iSide=1,nBCSides
  IF (BoundaryType(BC(iSide),1).EQ.4) THEN
    countSurfElem = countSurfElem + 1
  END IF
END DO

!IF (MPIroot) THEN
ALLOCATE(countSurfElemMPI(0:nProcessors-1))
countSurfElemMPI=0
!END IF

CALL MPI_GATHER(countSurfElem,1,MPI_INTEGER,countSurfElemMPI,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)

IF (MPIroot) THEN
DO iProc=1,nProcessors-1
  offsetSurfElemMPI(iProc)=SUM(countSurfElemMPI(0:iProc-1))
END DO
offsetSurfElemMPI(nProcessors)=SUM(countSurfElemMPI(:))
END IF

CALL MPI_BCAST (offsetSurfElemMPI,size(offsetSurfElemMPI),MPI_INTEGER,0,MPI_COMM_WORLD,iError)

offsetSurfElem=offsetSurfElemMPI(myRank)

DEALLOCATE(countSurfElemMPI)
#else /* MPI */
offsetSurfElem=0          ! offset is the index of first entry, hdf5 array starts at 0-.GT. -1 
#endif /* MPI */


! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----
! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale','1.0')
IF(ABS(meshScale-1.).GT.1e-14)&
  NodeCoords = NodeCoords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  Pi = ACOS(-1.) 
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=NodeCoords(:,i,j,k,iElem)
      NodeCoords(:,i,j,k,iElem) = x+ 0.1*SIN(Pi*x(1))*SIN(Pi*x(2))*SIN(Pi*x(3))
    END DO; END DO; END DO;
  END DO
END IF


ALLOCATE(Elem_xGP(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(BCFace_xGP(3,0:PP_N,0:PP_N,1:nBCSides))

! PoyntingVecIntegral
CalcPoyntingInt = GETLOGICAL('CalcPoyntingVecIntegral','.FALSE.')
IF (CalcPoyntingInt) ALLOCATE(Face_xGP(3,0:PP_N,0:PP_N,1:nSides))

ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(       NormVec(3,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)) 
ALLOCATE(      TangVec1(3,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper)) 
ALLOCATE(      TangVec2(3,0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))  
ALLOCATE(      SurfElem(  0:PP_N,0:PP_N,sideID_minus_lower:sideID_minus_upper))  

! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
! assign 1/detJ (sJ)
! assign normal and tangential vectors and surfElems on faces
#ifdef PARTICLES
ALLOCATE(BezierControlPoints3D(1:3,0:NGeo,0:NGeo,1:nSides) ) 
BezierControlPoints3D=0.

ALLOCATE(SideSlabNormals(1:3,1:3,1:nSides),SideSlabIntervals(1:6,nSides),BoundingBoxIsEmpty(1:nSides) )
SideSlabNormals=0.
SideSlabIntervals=0.
BoundingBoxIsEmpty=.TRUE.
ALLOCATE(ElemSlabNormals(1:3,1:3,1:nElems),ElemSlabIntervals(1:6,nElems) )
ElemSlabNormals=0.
ElemSlabIntervals=0.
#endif /*PARTICLES*/

crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
CALL InitMeshBasis(NGeo,PP_N,xGP)
CALL CalcMetrics(NodeCoords)
DEALLOCATE(NodeCoords)

debugmesh=GETLOGICAL('debugmesh','.FALSE.')
IF(debugmesh)THEN
  CALL  WriteDebugMesh()
END IF !/*debugmesh*/

ALLOCATE(VolToSideA(3,0:PP_N,0:PP_N,0:PP_N,0:4,1:6))
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO f=0,4
    DO s=1,6
      VolToSideA(:,i,j,k,f,s) = VolToSide(i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

ALLOCATE(CGNS_VolToSideA(3,0:PP_N,0:PP_N,0:PP_N,1:6))
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO s=1,6
    CGNS_VolToSideA(:,i,j,k,s) = CGNS_VolToSide(i,j,k,s)
  END DO
END DO; END DO; END DO


ALLOCATE(SideToVolA(3,0:PP_N,0:PP_N,0:PP_N,0:4,1:6))
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO f=0,4
    DO s=1,6
      SideToVolA(:,i,j,k,f,s) = SideToVol(i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

ALLOCATE(SideToVol2A(2,0:PP_N,0:PP_N,0:4,1:6))
DO j=0,PP_N; DO i=0,PP_N
  DO f=0,4
    DO s=1,6
      SideToVol2A(:,i,j,f,s) = SideToVol2(i,j,f,s)
    END DO
  END DO
END DO; END DO

ALLOCATE(CGNS_SideToVol2A(2,0:PP_N,0:PP_N,1:6))
DO j=0,PP_N; DO i=0,PP_N
  DO s=1,6
    CGNS_SideToVol2A(:,i,j,s) = CGNS_SideToVol2(i,j,s)
  END DO
END DO; END DO

DO f = 0, 4
  DO s = 1, 6
    DO q = 0,PP_N; DO p = 0,PP_N
      ijk = SideToVolA(:,0,p,q,f,s)
      pq = VolToSideA(:,ijk(1),ijk(2),ijk(3),f,s)
      IF ((pq(1).NE.p).OR.(pq(2).NE.q)) THEN
        CALL abort(__STAMP__,'SideToVol does not fit to VolToSideA')
      END IF
    END DO; END DO 
  END DO ! s = 1, 6
END DO ! f = 0, 4


#ifndef PARTICLES
DEALLOCATE(XCL_NGeo)
#else
! init element volume
CALL InitElemVolumes()
#endif

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


FUNCTION Flip_S2M(p, q, flip)
!===================================================================================================================================
! Transforms Coordinates from RHS of Slave to RHS of Master 
!    input: p,q in Slave-RHS, flip;  
!   output: indices in Master-RHS
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: Flip_S2M
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(flip)
  CASE(0) 
    Flip_S2M = (/     p,     q/)
  CASE(1) 
    Flip_S2M = (/     q,     p/)
  CASE(2) 
    Flip_S2M = (/PP_N-p,     q/)
  CASE(3) 
    Flip_S2M = (/PP_N-q,PP_N-p/)
  CASE(4) 
    Flip_S2M = (/     p,PP_N-q/)
END SELECT
END FUNCTION Flip_S2M


FUNCTION Flip_M2S(p, q, flip)
!===================================================================================================================================
! Transforms Coordinates from RHS of Master to RHS of Slave
!   actualy this is the same function as Flip_S2M
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: Flip_M2S
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Flip_M2S=Flip_S2M(p,q,flip)
END FUNCTION Flip_M2S


FUNCTION CGNS_VolToSide(i,j,k, locSideID)
!===================================================================================================================================
! Transforms Volume-Coordinates into RHS of the Side (uses CGNS-Notation)
! input: i,j,k, locSideID 
!   where: i,j,k = volume-indices 
! output: indices in Master-RHS  +  volume-index which is not used (depending on locSideID)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: i,j,k,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(3) :: CGNS_VolToSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(locSideID)
  CASE(XI_MINUS)   
    CGNS_VolToSide = (/k,j,i/)
  CASE(XI_PLUS)    
    CGNS_VolToSide = (/j,k,PP_N-i/)
  CASE(ETA_MINUS)  
    CGNS_VolToSide = (/i,k,j/)
  CASE(ETA_PLUS)   
    CGNS_VolToSide = (/PP_N-i,k,PP_N-j/)
  CASE(ZETA_MINUS) 
    CGNS_VolToSide = (/j,i,k/)
  CASE(ZETA_PLUS)  
    CGNS_VolToSide = (/i,j,PP_N-k/)
END SELECT
END FUNCTION CGNS_VolToSide


FUNCTION CGNS_SideToVol(l, p, q, locSideID)
!===================================================================================================================================
! Transforms RHS-Coordinates of Side (CGNS-Notation) into Volume-Coordinates 
! input: l, p,q, locSideID
!   where: p,q are in Master-RHS;
!          l is the xi-,eta- or zeta-index in 0:PP_N corresponding to locSideID
! output: volume-indicies
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: l,p,q,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(3) :: CGNS_SideToVol
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(locSideID)
  CASE(XI_MINUS)   
    CGNS_SideToVol = (/l,q,p/)
  CASE(XI_PLUS)    
    CGNS_SideToVol = (/PP_N-l,p,q/)
  CASE(ETA_MINUS)  
    CGNS_SideToVol = (/p,l,q/)
  CASE(ETA_PLUS)   
    CGNS_SideToVol = (/PP_N-p,PP_N-l,q/)
  CASE(ZETA_MINUS) 
    CGNS_SideToVol = (/q,p,l/)
  CASE(ZETA_PLUS)  
    CGNS_SideToVol = (/p,q,PP_N-l/)
END SELECT
END FUNCTION CGNS_SideToVol


FUNCTION CGNS_SideToVol2(p, q, locSideID)
!===================================================================================================================================
! Transforms RHS-Coordinates of Side (CGNS-Notation) into Volume-Coordinates 
! input: l, p,q, locSideID
!   where: p,q are in Master-RHS;
!          l is the xi-,eta- or zeta-index in 0:PP_N corresponding to locSideID
! output: volume-indicies
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: CGNS_SideToVol2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(locSideID)
  CASE(XI_MINUS)   
    CGNS_SideToVol2 = (/q,p/)
  CASE(XI_PLUS)    
    CGNS_SideToVol2 = (/p,q/)
  CASE(ETA_MINUS)  
    CGNS_SideToVol2 = (/p,q/)
  CASE(ETA_PLUS)   
    CGNS_SideToVol2 = (/PP_N-p,q/)
  CASE(ZETA_MINUS) 
    CGNS_SideToVol2 = (/q,p/)
  CASE(ZETA_PLUS)  
    CGNS_SideToVol2 = (/p,q/)
END SELECT
END FUNCTION CGNS_SideToVol2


FUNCTION VolToSide(i,j,k, flip, locSideID)
!===================================================================================================================================
! Transform Volume-Coordinates to RHS-Coordinates of Master. This is: VolToSide = Flip_S2M(CGNS_VolToSide(...))
! input: i,j,k, flip, locSideID 
!   where: i,j,k = volume-indices 
! output: indices in Master-RHS
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: i,j,k,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(3) :: VolToSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
!===================================================================================================================================
pq = CGNS_VolToSide(i,j,k,locSideID)
VolToSide(1:2) = Flip_S2M(pq(1),pq(2),flip)
VolToSide(3) = pq(3)
END FUNCTION VolToSide


FUNCTION SideToVol(l, p, q, flip, locSideID)
!===================================================================================================================================
! Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol = CGNS_SideToVol(Flip_M2S(...))
! input: l, p,q, flip, locSideID
!     where: p,q are in Master-RHS;
!            l is the xi-,eta- or zeta-index in 0:PP_N corresponding to locSideID
! output: volume-indicies
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: l,p,q,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(3) :: SideToVol
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!===================================================================================================================================
pq = Flip_M2S(p,q,flip)
SideToVol = CGNS_SideToVol(l,pq(1),pq(2),locSideID)
END FUNCTION SideToVol


FUNCTION SideToVol2(p, q, flip, locSideID)
!===================================================================================================================================
! Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol = CGNS_SideToVol(Flip_M2S(...))
! input: l, p,q, flip, locSideID
!     where: p,q are in Master-RHS;
!            l is the xi-,eta- or zeta-index in 0:PP_N corresponding to locSideID
! output: volume-indicies
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: SideToVol2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!===================================================================================================================================
pq = Flip_M2S(p,q,flip)
SideToVol2 = CGNS_SideToVol2(pq(1),pq(2),locSideID)
END FUNCTION SideToVol2


FUNCTION ElemToNBElem(locSideID,iElem) 
!===================================================================================================================================
! get index of neighboring Elem, return -1 if none exists
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: locSideID,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: ElemToNBElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: flip, SideID, i
!===================================================================================================================================
SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
IF ((SideID.GE.nBCSides+1).AND.(SideID.LE.nBCSides+nInnerSides)) THEN
  flip  =ElemToSide(E2S_FLIP,locSideID,iElem)
  IF (flip.EQ.0) THEN
    ElemToNBElem = SideToElem(S2E_NB_ELEM_ID,SideID)
  ELSE
    ElemToNBElem = SideToElem(S2E_ELEM_ID,SideID)
  END IF
ELSE 
  ElemToNBElem = -1
END IF
END FUNCTION ElemToNBElem


FUNCTION VolToVol(i,j,k,locSideID,iElem,neighbor_iElem)
!===================================================================================================================================
! Transform Volume-Coordinates to neighboring Volume-Coordinates.  This is: VolToVol = SideToVol(VolToSide(...))
! input: i,j,k, iElem, locSideID
!     where: i,j,k  are Volume-Indizices of element iElem;
!            locSideID  side to the neighboring element 
! output: volume-indicies of neighboring elemnt, that are next to ijk in direction of locSideID
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: i,j,k,locSideID,iElem,neighbor_iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(3) :: VolToVol
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
INTEGER              :: flip, l, SideID, neighbor_flip, neighbor_locSideID
!===================================================================================================================================
SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
flip  =ElemToSide(E2S_FLIP,locSideID,iElem)
pq = VolToSide(i,j,k, flip, locSideID)
l = pq(3)
IF (flip.EQ.0) THEN
  neighbor_locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  neighbor_flip      = SideToElem(S2E_FLIP,SideID)
ELSE 
  neighbor_locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  neighbor_flip      = 0
END IF
VolToVol = SideToVol(l, pq(1), pq(2), neighbor_flip, neighbor_locSideID)
END FUNCTION VolToVol


SUBROUTINE InitMeshBasis(NGeo_in,N_in,xGP)
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,               ONLY: Xi_NGeo,Vdm_CLN_GaussN,Vdm_CLNGeo_CLN,Vdm_CLNGeo_GaussN,Vdm_NGeo_CLNGeo,DCL_NGeo,DCL_N&
                                       ,wBaryCL_NGeo,XiCL_NGeo,DeltaXi_NGeo
USE MOD_Basis,                   ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,                   ONLY: ChebyGaussLobNodesAndWeights,PolynomialDerivativeMatrix,InitializeVandermonde
#ifdef PARTICLES
USE MOD_Mesh_Vars,               ONLY: wBaryCL_NGeo1,Vdm_CLNGeo1_CLNGeo,XiCL_NGeo1
USE MOD_Particle_Surfaces_Vars,  ONLY: Vdm_Bezier,sVdm_Bezier
USE MOD_Basis,                   ONLY: BuildBezierVdm
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: NGeo_in,N_in
REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in)                     :: XiCL_N,wBaryCL_N
REAL,DIMENSION(0:NGeo_in)                  :: wBary_NGeo!: XiCL_NGeo,!,wBaryCL_NGeo,wBary_NGeo
#ifdef PARTICLES
!REAL,DIMENSION(0:NGeo_in)                  :: XiEquiPartCurved
#endif
INTEGER                                    :: i
!===================================================================================================================================
ALLOCATE(DCL_N(0:N_in,0:N_in),Vdm_CLN_GaussN(0:N_in,0:N_in))
ALLOCATE(Xi_NGeo(0:NGeo_in))
ALLOCATE(DCL_NGeo(0:NGeo_in,0:NGeo_in))
ALLOCATE(Vdm_CLNGeo_GaussN(0:N_in,0:NGeo_in))
ALLOCATE(Vdm_CLNGeo_CLN(0:N_in,0:NGeo_in))
ALLOCATE(Vdm_NGeo_CLNGeo(0:NGeo_in,0:NGeo_in))

ALLOCATE(wBaryCL_NGeo(0:NGeo_In))
ALLOCATE(XiCL_NGeo(0:NGeo_In))
! Chebyshev-Lobatto N
CALL ChebyGaussLobNodesAndWeights(N_in,XiCL_N)
CALL BarycentricWeights(N_in,XiCL_N,wBaryCL_N)
CALL PolynomialDerivativeMatrix(N_in,XiCL_N,DCL_N)
CALL InitializeVandermonde(N_in,N_in,wBaryCL_N,XiCL_N,xGP,Vdm_CLN_GaussN)
!equidistant-Lobatto NGeo
DO i=0,NGeo_in
  Xi_NGeo(i) = 2./REAL(NGeo_in) * REAL(i) - 1. 
END DO
DeltaXi_NGeo=2./NGeo_in
CALL BarycentricWeights(NGeo_in,Xi_NGeo,wBary_NGeo)

! Chebyshev-Lobatto NGeo
CALL ChebyGaussLobNodesAndWeights(NGeo_in,XiCL_NGeo)
CALL BarycentricWeights(NGeo_in,XiCL_NGeo,wBaryCL_NGeo)
CALL PolynomialDerivativeMatrix(NGeo_in,XiCL_NGeo,DCL_NGeo)

CALL InitializeVandermonde(NGeo_in,N_in   ,wBaryCL_NGeo,XiCL_NGeo,xGP      ,Vdm_CLNGeo_GaussN)
CALL InitializeVandermonde(NGeo_in,N_in   ,wBaryCL_NGeo,XiCL_NGeo,XiCL_N   ,Vdm_CLNGeo_CLN   )
CALL InitializeVandermonde(NGeo_in,NGeo_in,wBary_NGeo  ,Xi_NGeo  ,XiCL_NGeo,Vdm_NGeo_CLNGeo  )
#ifdef PARTICLES
! small wBaryCL_NGEO
ALLOCATE(wBaryCL_NGeo1(0:1),XiCL_NGeo1(0:1))
CALL ChebyGaussLobNodesAndWeights(1,XiCL_NGeo1)
CALL BarycentricWeights(1,XiCL_NGeo1,wBaryCL_NGeo1)
ALLOCATE(Vdm_CLNGeo1_CLNGeo(0:NGeo_In,0:1) )
CALL InitializeVandermonde(1, NGeo_in,wBaryCL_NGeo1,XiCL_NGeo1,XiCL_NGeo ,Vdm_CLNGeo1_CLNGeo)

! new for curved particle sides
ALLOCATE(Vdm_Bezier(0:NGeo_in,0:NGeo_in),sVdm_Bezier(0:NGeo_in,0:NGeo_in))
! initialize vandermonde for super-sampled surfaces (particle tracking with curved elements)
!DO i=0,NGeo_in
!  XiEquiPartCurved(i) = 2./REAL(NGeo_in) * REAL(i) - 1. 
!END DO
! initialize vandermonde for bezier basis surface representation (particle tracking with curved elements)
CALL BuildBezierVdm(NGeo_in,XiCL_NGeo,Vdm_Bezier,sVdm_Bezier) !CHANGETAG
#endif

END SUBROUTINE InitMeshBasis


SUBROUTINE FinalizeMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
#ifdef PARTICLES
!USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(Xi_NGeo)
SDEALLOCATE(DCL_N)
SDEALLOCATE(DCL_NGeo)
SDEALLOCATE(VdM_CLN_GaussN)
SDEALLOCATE(VdM_CLNGeo_GaussN)
SDEALLOCATE(Vdm_CLNGeo_CLN)
SDEALLOCATE(Vdm_NGeo_CLNgeo)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(SideToElem2)
SDEALLOCATE(BC)
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(BCFace_xGP)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(sJ)
SDEALLOCATE(NormVec) 
SDEALLOCATE(TangVec1) 
SDEALLOCATE(TangVec2)  
SDEALLOCATE(SurfElem)  
SDEALLOCATE(Face_xGP)
SDEALLOCATE(XCL_NGeo)
SDEALLOCATE(dXCL_NGeo)
SDEALLOCATE(Face_xGP)
SDEALLOCATE(wbaryCL_NGeo)
SDEALLOCATE(XiCL_NGeo)
SDEALLOCATE(VolToSideA)
SDEALLOCATE(CGNS_VolToSideA)
SDEALLOCATE(SideToVolA)
SDEALLOCATE(SideToVol2A)
SDEALLOCATE(CGNS_SideToVol2A)
SDEALLOCATE(Vdm_CLNGeo1_CLNGeo)
SDEALLOCATE(wBaryCL_NGeo1)
SDEALLOCATE(XiCL_NGeo1)
SDEALLOCATE(CurvedElem)
MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
