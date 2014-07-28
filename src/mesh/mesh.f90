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
USE MOD_Interpolation_Vars, ONLY:xGP,InterpolationInitIsDone
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo,fillElemGeo,getVolumeMapping
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETINT,GETINTARRAY,CNTSTR,GETSTR
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_Metrics,            ONLY:CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Analyze_Vars,       ONLY:CalcPoyntingInt
#ifdef PARTICLES
USE MOD_ParticleInit,       ONLY:InitParticleGeometry,InitElemVolumes
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
INTEGER           :: iBC
LOGICAL           :: debugmesh
REAL,ALLOCATABLE  :: XCL_NGeo(:,:,:,:,:)
REAL              :: x(3),PI
INTEGER           :: iElem,i,j,k,iSide,countSurfElem,iProc
INTEGER,ALLOCATABLE :: countSurfElemMPI(:)
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL abort(__STAMP__,'InitMesh not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

NGeo=GETINT('GeometricNGeo','1') 
CALL initMeshBasis(NGeo,PP_N,xGP)
MeshType=GETINT('MeshType','1')

useCurveds=GETLOGICAL('useCurveds','.FALSE.')
IF(usecurveds)THEN
  IF(PP_N.LT.NGeo) THEN
    SWRITE(*,*)'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                                can cause problems on periodic boundaries! Set N>NGeo'
  END IF
END IF

IF(.NOT.(usecurveds))THEN
  IF(NGeo.GT.1) THEN
    CALL abort(__STAMP__,'NGeo > 1 is only allowed for curved hexahedra.',999,999.)
  END IF
END IF
! read in boundary conditions, will overwrite BCs from meshfile!
nUserBCs = CNTSTR('BoundaryName','0')
IF(nUserBCs.GT.0)THEN
  ALLOCATE(BoundaryName(1:nUserBCs))
  ALLOCATE(BoundaryType(1:nUserBCs,2))
  DO iBC=1,nUserBCs
    BoundaryName(iBC)   = GETSTR('BoundaryName')
    BoundaryType(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
  END DO
END IF !nUserBCs>0

! prepare pointer structure (get nElems, etc.)
MeshFile = GETSTR('MeshFile')

CALL readMesh(MeshFile) !set nElems

!schmutz fink
PP_nElems=nElems

SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
CALL setLocalSideIDs()

ALLOCATE(XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems))
XCL_NGeo = 0.

! map pointer structure to XCL_NGeo
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillElemGeo..."
CALL fillElemGeo(XCL_NGeo)

#ifdef MPI
! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
CALL exchangeFlip()
#endif

! fill ElemToSide, SideToElem,BC
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(SideToElem2(4,2*nInnerSides+nBCSides+nMPISides))
ALLOCATE(BC(1:nBCSides))
ElemToSide  = 0
SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
SideToElem2 = -1   !mapping side to elem, sorted by elem ID (for ProlongToFace) 
BC         = 0

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
CALL InitParticleGeometry()
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

! IF(NGeo.GT.1) CALL getVolumeMapping(XCL_NGeo)

IF(GETLOGICAL('deform','.FALSE.'))THEN
  Pi = ACOS(-1.) 
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=XCL_Ngeo(:,i,j,k,iElem)
      XCL_Ngeo(:,i,j,k,iElem) = x+ 0.1*SIN(Pi*x(1))*SIN(Pi*x(2))*SIN(Pi*x(3))
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

!SWRITE(*,*) "InterpolateToGaussPoints"
!CALL InterpolateToGaussPoints(XCL_NGeo)

! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
! assign 1/detJ (sJ)
! assign normal and tangential vectors and surfElems on faces
crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
CALL CalcMetrics(XCL_NGeo) 
#ifdef PARTICLES
! save geometry information for particle tracking
CALL InitElemVolumes()
#endif

debugmesh=GETLOGICAL('debugmesh','.FALSE.')
IF(debugmesh)THEN
  CALL  WriteDebugMesh()
END IF !/*debugmesh*/

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


SUBROUTINE InitMeshBasis(NGeo_in,N_in,xGP)
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY: Xi_NGeo,Vdm_CLN_GaussN,Vdm_CLNGeo_CLN,Vdm_CLNGeo_GaussN,Vdm_NGeo_CLNGeo,DCL_NGeo,DCL_N
USE MOD_Basis,     ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,BarycentricWeights
USE MOD_Basis,     ONLY: ChebyGaussLobNodesAndWeights,PolynomialDerivativeMatrix,InitializeVandermonde
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
REAL,DIMENSION(0:NGeo_in)                  :: XiCL_NGeo,wBaryCL_NGeo,wBary_NGeo
INTEGER                                    :: i
!===================================================================================================================================
ALLOCATE(DCL_N(0:N_in,0:N_in),Vdm_CLN_GaussN(0:N_in,0:N_in))
ALLOCATE(Xi_NGeo(0:NGeo_in))
ALLOCATE(DCL_NGeo(0:NGeo_in,0:NGeo_in))
ALLOCATE(Vdm_CLNGeo_GaussN(0:N_in,0:NGeo_in))
ALLOCATE(Vdm_CLNGeo_CLN(0:N_in,0:NGeo_in))
ALLOCATE(Vdm_NGeo_CLNGeo(0:NGeo_in,0:NGeo_in))
! Chebyshev-Lobatto N
CALL ChebyGaussLobNodesAndWeights(N_in,XiCL_N)
CALL BarycentricWeights(N_in,XiCL_N,wBaryCL_N)
CALL PolynomialDerivativeMatrix(N_in,XiCL_N,DCL_N)
CALL InitializeVandermonde(N_in,N_in,wBaryCL_N,XiCL_N,xGP,Vdm_CLN_GaussN)
!equidistant-Lobatto NGeo
DO i=0,NGeo_in
  Xi_NGeo(i) = 2./REAL(NGeo_in) * REAL(i) - 1. 
END DO
CALL BarycentricWeights(NGeo_in,Xi_NGeo,wBary_NGeo)

! Chebyshev-Lobatto NGeo
CALL ChebyGaussLobNodesAndWeights(NGeo_in,XiCL_NGeo)
CALL BarycentricWeights(NGeo_in,XiCL_NGeo,wBaryCL_NGeo)
CALL PolynomialDerivativeMatrix(NGeo_in,XiCL_NGeo,DCL_NGeo)

CALL InitializeVandermonde(NGeo_in,N_in   ,wBaryCL_NGeo,XiCL_NGeo,xGP      ,Vdm_CLNGeo_GaussN)
CALL InitializeVandermonde(NGeo_in,N_in   ,wBaryCL_NGeo,XiCL_NGeo,XiCL_N   ,Vdm_CLNGeo_CLN   )
CALL InitializeVandermonde(NGeo_in,NGeo_in,wBary_NGeo  ,Xi_NGeo  ,XiCL_NGeo,Vdm_NGeo_CLNGeo  )
END SUBROUTINE InitMeshBasis


SUBROUTINE InterpolateToGaussPoints(XCL_NGeo)
!===================================================================================================================================
! interpolate from 3D Chebyshev points (NGeo+1)x(NGeo+1)x(NGeo+1) onto inner 3D Gauss (or GL) points 
! and to 2D Gauss Points (N+1)x(N+1) on BC faces
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY:NGeo,Vdm_CLNGeo_GaussN
USE MOD_Mesh_Vars,          ONLY:nElems,nBCSides
USE MOD_Mesh_Vars,          ONLY:ElemToSide
USE MOD_Mesh_Vars,          ONLY:Elem_xGP,BCFace_xGP,Face_xGP
USE MOD_Analyze_Vars,       ONLY:CalcPoyntingInt
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_ChangeBasis,        ONLY:ChangeBasis2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo,nElems) ! XYZ postions of the cheb points
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iLocSide,iElem,p,q
INTEGER          :: SideID
REAL             :: tmp(3,0:PP_N,0:PP_N)
!===================================================================================================================================
! assign coordinates to 3D Element interpolation points (points for computation, either Gauss or Gauss-Lob)
DO iElem=1,nElems
  CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,:,iElem),Elem_xGP(:,:,:,:,iElem))
  ! assign coordinates to 2D BC face interpolation points (points for computation, either Gauss or Gauss-Lob)
  DO iLocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    IF(SideID .LE. nBCSides)THEN  ! BC side
      SELECT CASE(iLocSide)
      CASE(XI_MINUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,   0,:,:,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            BCFace_xGP(:,p,q,SideID)=tmp(:,q,p)
          END DO !p
        END DO !q
      CASE(XI_PLUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,NGeo,:,:,iElem),BCFace_xGP(:,:,:,SideID))
      CASE(ETA_MINUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,   0,:,iElem),BCFace_xGP(:,:,:,SideID))
      CASE(ETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,NGeo,:,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            BCFace_xGP(:,p,q,SideID)=tmp(:,PP_N-p,q)
          END DO !p
        END DO !q
      CASE(ZETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,   0,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            BCFace_xGP(:,p,q,SideID)=tmp(:,q,p)
          END DO !p
        END DO !q
      CASE(ZETA_PLUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,NGeo,iElem),BCFace_xGP(:,:,:,SideID))
      END SELECT !LocSideID
    END IF !(SideID .LE. nBCSides)
  END DO ! iLocSide
END DO !iElem=1,nElems

IF(CalcPoyntingInt) THEN
  ! neccessary to get coordinates of all faces for later check
  ! assign coordinates to 3D Element interpolation points (points for computation, either Gauss or Gauss-Lob)
  DO iElem=1,nElems
    CALL ChangeBasis3D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,:,iElem),Elem_xGP(:,:,:,:,iElem))
    ! assign coordinates to 2D BC face interpolation points (points for computation, either Gauss or Gauss-Lob)
    DO iLocSide=1,6
      SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
      SELECT CASE(iLocSide)
      CASE(XI_MINUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,   0,:,:,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            Face_xGP(:,p,q,SideID)=tmp(:,q,p)
          END DO !p
        END DO !q
      CASE(XI_PLUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,NGeo,:,:,iElem),Face_xGP(:,:,:,SideID))
      CASE(ETA_MINUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,   0,:,iElem),Face_xGP(:,:,:,SideID))
      CASE(ETA_PLUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,NGeo,:,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            Face_xGP(:,p,q,SideID)=tmp(:,PP_N-p,q)
          END DO !p
        END DO !q
      CASE(ZETA_MINUS)
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,   0,iElem),tmp)
        ! turn into local right hand system of side
        DO q=0,PP_N
          DO p=0,PP_N
            Face_xGP(:,p,q,SideID)=tmp(:,q,p)
          END DO !p
        END DO !q
      CASE(ZETA_PLUS) 
        CALL ChangeBasis2D(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,XCL_NGeo(:,:,:,NGeo,iElem),Face_xGP(:,:,:,SideID))
      END SELECT !LocSideID
    END DO ! iLocSide
  END DO !iElem=1,nElems
END IF ! CalcPoyntingInt = .TRUE.
END SUBROUTINE InterpolateToGaussPoints


SUBROUTINE FinalizeMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Analyze_Vars,       ONLY:CalcPoyntingInt
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

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
