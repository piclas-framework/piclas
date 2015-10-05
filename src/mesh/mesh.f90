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

PUBLIC::InitMesh
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
USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
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
INTEGER           :: iElem,i,j,k
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
ALLOCATE(SlabNormals(1:3,1:3,1:nSides),SlabIntervalls(1:6,nSides),BoundingBoxIsEmpty(1:nSides) )
SlabNormals=0.
SlabIntervalls=0.
BoundingBoxIsEmpty=.TRUE.
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
!USE MOD_Particle_Surfaces_Vars, ONLY:BezierControlPoints3D,SlabNormals,SlabIntervalls,BoundingBoxIsEmpty
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

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
