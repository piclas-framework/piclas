#include "boltzplatz.h"

MODULE MOD_Particle_Surfaces
!===================================================================================================================================
! Contains subroutines to build the requiered data to track particles on (curviilinear) meshes, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE GetSideType
  MODULE PROCEDURE GetSideType
END INTERFACE

INTERFACE GetBCSideType
  MODULE PROCEDURE GetBCSideType
END INTERFACE

INTERFACE InitParticleSurfaces
  MODULE PROCEDURE InitParticleSurfaces
END INTERFACE

INTERFACE FinalizeParticleSurfaces
  MODULE PROCEDURE FinalizeParticleSurfaces
END INTERFACE

INTERFACE CalcBiLinearNormVec
  MODULE PROCEDURE CalcBiLinearNormVec
END INTERFACE

!INTERFACE GetSuperSampledSurface
!  MODULE PROCEDURE GetSuperSampledSurface
!END INTERFACE

INTERFACE GetBezierControlPoints3D
  MODULE PROCEDURE GetBezierControlPoints3D
END INTERFACE

INTERFACE CalcNormVecBezier
  MODULE PROCEDURE CalcNormVecBezier
END INTERFACE

INTERFACE CalcBiLinearNormVecBezier
  MODULE PROCEDURE CalcBiLinearNormVecBezier
END INTERFACE

INTERFACE CalcBiLinearNormAndTang
  MODULE PROCEDURE CalcBiLinearNormAndTang
END INTERFACE

INTERFACE CalcNormAndTangBezier
  MODULE PROCEDURE CalcNormAndTangBezier
END INTERFACE

INTERFACE GetSlabNormalsAndIntervalls
  MODULE PROCEDURE GetSlabNormalsAndIntervalls
END INTERFACE


PUBLIC::GetSideType, InitParticleSurfaces, FinalizeParticleSurfaces, CalcBiLinearNormVec, &!GetSuperSampledSurface, &
        GetBezierControlPoints3D,CalcBiLinearNormVecBezier,CalcNormVecBezier, GetSlabNormalsAndIntervalls

PUBLIC::GetBCSideType,CalcBiLinearNormAndTang, CalcNormAndTangBezier

!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_vars
USE MOD_Preproc
USE MOD_Mesh_Vars,                  ONLY:nSides,ElemToSide,SideToElem,NGeo,nBCSides,BC,nSides
USE MOD_ReadInTools,                ONLY:GETREAL,GETINT,GETLOGICAL
USE MOD_Particle_Vars,              ONLY:PDM
USE MOD_Particle_Mesh_Vars,         ONLY:PartBCSideList,nTotalBCSides
USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
!USE MOD_Particle_SFC_Vars,          ONLY:whichBoundBox
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: tmp,iSide!,iBCSide
CHARACTER(LEN=2)                :: dummy                         
!===================================================================================================================================

IF(ParticleSurfaceInitIsDone) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES ...!'

epsilontol      = GETREAL('epsOne','1e-12')
hitepsbi        = GETREAL('hitepsbi','1e-12')
hitepsbi=1.0+hitepsbi
!hitepsbi=1.0+SQRT(EPSILON(0.0))
!hitepsbi=1.000800
Mepsilontol     = -epsilontol
epsilonOne      = 1.0 + epsilontol
OneMepsilon     = 1.0 - epsilontol
ClipTolerance   = GETREAL('ClipTolerance','1e-4')
SplitLimit      = GETREAL('SplitLimit','0.4')
SplitLimit      =2.*SplitLimit
ClipMaxIter     = GETINT('ClipMaxIter','10')
ClipHit         = GETREAL('ClipHit','1e-7')
ClipHit =1.0+ClipHit
tmp=2*(NGeo+1)
WRITE(dummy,'(I2.2)') tmp
ClipMaxInter    = GETINT('ClipMaxInter',dummy)

IF(DoRefMapping)THEN
  !MultipleBCs    = GETLOGICAL('MultibleBCs',".FALSE.")
  ALLOCATE(PartBCSideList(1:nSides))
  PartBCSideList(:) = -1
  DO iSide=1,nBCSides
    PartBCSideList(iSide)=iSide
  END DO 
  nTotalBCSides=nBCSides
 ! iBCSide=nBCSides
 ! DO iSide=nBCSides+1,nSides
 !   IF(BC(iSide).EQ.1) THEN
 !     iBCSide=iBCSide+1
 !     PartBCSideList(iSide)=iBCSide
 !   END IF
 ! END DO 
 ! nTotalBCSides=iBCSide
END IF

!! ElemBaryNGeo are required for particle mapping| SingleParticleToExactElem
!IF(.NOT.DoRefMapping)THEN
!!   ALLOCATE(XiEtaZetaBasis(1:3,1:6,1:PP_nElems) &
!!           ,slenXiEtaZetaBasis(1:6,1:PP_nElems) &
!!           ,ElemRadiusNGeo(1:PP_nElems)         &
!!           ,ElemBaryNGeo(1:3,1:PP_nElems)       )
!!   CALL BuildElementBasis()
!ELSE
!  !whichBoundBox = GETINT('PartSFC-BoundBox','1')
!END IF

ALLOCATE( locAlpha(1:ClipMaxInter) &
        , locXi   (1:ClipMaxInter) &
        , locEta  (1:ClipMaxInter) )
ALLOCATE( XiArray (1:2,1:ClipMaxIter) &
        , EtaArray(1:2,1:ClipMaxIter) )


! moved into mesh init
! construct connections to neighbor elems
!ALLOCATE( neighborElemID    (1:6,1:PP_nElems) &
!        , neighborlocSideID (1:6,1:PP_nElems) )
!neighborElemID=-1
!neighborlocSideID=-1

!DO iElem=1,PP_nElems
!  DO ilocSide=1,6
!    flip = ElemToSide(E2S_FLIP,ilocSide,iElem)
!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(flip.EQ.0)THEN
!      ! SideID of slave
!      neighborlocSideID(ilocSide,iElem)=SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
!      neighborElemID   (ilocSide,iElem)=SideToElem(S2E_NB_ELEM_ID,SideID)
!    ELSE
!      ! SideID of master
!      neighborlocSideID(ilocSide,iElem)=SideToElem(S2E_LOC_SIDE_ID,SideID)
!      neighborElemID   (ilocSide,iElem)=SideToElem(S2E_ELEM_ID,SideID)
!    END IF
!  END DO ! ilocSide
!END DO ! Elem

ParticleSurfaceInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitParticleSurfaces

SUBROUTINE FinalizeParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_vars
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
!===================================================================================================================================

SDEALLOCATE(SideType)
!SDEALLOCATE(BiLinearCoeff)
SDEALLOCATE(SideNormVec)
SDEALLOCATE(SideDistance)
SDEALLOCATE(BezierControlPoints3D)
!SDEALLOCATE(SuperSampledBiLinearCoeff)
SDEALLOCATE(SlabNormals)
SDEALLOCATE(SlabIntervalls)
SDEALLOCATE(BoundingBoxIsEmpty)
SDEALLOCATE(locAlpha)
SDEALLOCATE(locXi)
SDEALLOCATE(locEta)
SDEALLOCATE(XiArray)
SDEALLOCATE(EtaArray)
SDEALLOCATE(Vdm_Bezier)
SDEALLOCATE(sVdm_Bezier)
SDEALLOCATE(arrayNChooseK)
SDEALLOCATE(FacNchooseK)
SDEALLOCATE(SideType)
!SDEALLOCATE(gElemBCSide)
ParticleSurfaceInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleSurfaces

! obsolete
!SUBROUTINE GetBilinearPlane()
!===================================================================================================================================
! computes the required coefficients for a bi-linear plane and performs the decision between planar and bi-linear planes
!===================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_Mesh_Vars,                ONLY:nSides,ElemToSide
!USE MOD_Particle_Vars,            ONLY:GEO
!USE MOD_Particle_Surfaces_Vars,   ONLY:epsilonbilinear, SideIsPlanar,BiLinearCoeff, SideNormVec, SideDistance
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!LOGICAL,ALLOCATABLE               :: SideIsDone(:)
!REAL                              :: Displacement,nlength
!INTEGER                           :: iElem,ilocSide, SideID,iNode,iSide
!! debug information
!INTEGER                           :: nBilinear,nPlanar
!===================================================================================================================================
!
!ALLOCATE( SideIsDone(1:nSides) )
!SideIsDone=.FALSE.
!nBiLinear=0
!nPlanar=0
!
!DO iElem=1,PP_nElems ! caution, if particles are not seeded in the whole domain
!  DO ilocSide=1,6
!    SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem) 
!    IF(.NOT.SideIsDone(SideID))THEN
!
!      ! for ray-bi-linear patch intersection see. ramsay
!      ! compute the bi-linear coefficients for this side
!      ! caution: the parameter space is [-1;1] x [-1;1] instead of [0,1]x[0,2] 
!      ! the numbering of the nodes should be counterclockwise 
!      ! DEBUGGGG!!!!!!!!!!!!!!!!
!      ! check if the nodes are in correct numbering
!      ! for ray-bi-linear patch intersection see. ramsay
!      BiLinearCoeff(:,1,SideID) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
!                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
!                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))
!
!      BiLinearCoeff(:,2,SideID) =-GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
!                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))
!
!      BiLinearCoeff(:,3,SideID) =-GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
!                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))
!
!      BiLinearCoeff(:,4,SideID) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
!                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))
!      BiLinearCoeff(:,:,SideID) = 0.25*BiLinearCoeff(:,:,SideID)
!      ! compute displacement vector (is displacement form planar plane)
!      Displacement = BiLinearCoeff(1,1,SideID)*BiLinearCoeff(1,1,SideID) &
!                   + BiLinearCoeff(2,1,SideID)*BiLinearCoeff(2,1,SideID) &
!                   + BiLinearCoeff(3,1,SideID)*BiLinearCoeff(3,1,SideID) 
!      IF(Displacement.LT.epsilonbilinear)THEN
!        SideIsPlanar(SideID)=.TRUE.
!print*,"SideIsPlanar(",SideID,")=",SideIsPlanar(SideID)
!        nPlanar=nPlanar+1
!      ELSE
!        nBilinear=nBilinear+1
!      END IF
!      SideIsDone(SideID)=.TRUE.
!    ELSE
!      CYCLE  
!    END IF ! SideID
!  END DO ! ilocSide
!END DO ! iElem
!
!! get number of bc-sides of each element
!! nElemBCSides=0
!! DO iElem=1,PP_nelems
!!   DO ilocSide=1,6
!!     SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
!!     IF(SideID.LT.nBCSides) nElemBCSides=nElemBCSides+1
!!   END DO ! ilocSide
!! END DO ! nElemBCSides
!
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar    surfaces: ', nPlanar
!SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear surfaces: ', nBilinear
!ALLOCATE(SideNormVec(1:3,nSides))
!SideNormVec=0.
!! compute normal vector of planar sides
!DO iSide=1,nSides
!  IF(SideIsPlanar(SideID))THEN
!    SideNormVec(:,iSide)=CROSS(BiLinearCoeff(:,2,iSide),BiLinearCoeff(:,3,iSide))
!    nlength=SideNormVec(1,iSide)*SideNormVec(1,iSide) &
!           +SideNormVec(2,iSide)*SideNormVec(2,iSide) &
!           +SideNormVec(3,iSide)*SideNormVec(3,iSide) 
!    SideNormVec(:,iSide) = SideNormVec(:,iSide)/SQRT(nlength)
!    SideDistance(iSide)  = DOT_PRODUCT(SideNormVec(:,iSide),BiLinearCoeff(:,4,iSide))
!  END IF
!END DO ! iSide
!
!DEALLOCATE( SideIsDone)
!
!END SUBROUTINE GetBilinearPlane

FUNCTION CalcBiLinearNormVec(xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Globals,                              ONLY:CROSS
USE MOD_Particle_Surfaces_Vars,               ONLY:BiLinearCoeff
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(3)                      :: CalcBiLinearNormVec
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: a,b,nVec
REAL                                   :: nlength
!================================================================================================================================

a=xi* BiLinearCoeff(:,1,SideID)+BiLinearCoeff(:,2,SideID)
b=eta*BiLinearCoeff(:,1,SideID)+BiLinearCoeff(:,3,SideID)

nVec=CROSS(a,b)
nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2)+nVec(3)*nVec(3)
nlength=SQRT(nlength)
CalcBiLinearNormVec=nVec/nlength

END FUNCTION CalcBiLinearNormVec

FUNCTION CalcBiLinearNormVecBezier(xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Globals,                              ONLY:CROSS
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,DIMENSION(3)                      :: CalcBilinearNormVecBezier
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: a,b,nVec
REAL                                   :: nlength
!================================================================================================================================


b=xi*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)  & 
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
   +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)   &
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)+BezierControlPoints3D(:,0   ,NGeo,SideID) )

a=eta*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
    +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) )


nVec=CROSS(a,b)
nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2)+nVec(3)*nVec(3)
nlength=SQRT(nlength)
CalcBiLinearNormVecBezier=nVec/nlength
END FUNCTION CalcBiLinearNormVecBezier


SUBROUTINE CalcBiLinearNormAndTang(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Globals,                              ONLY:CROSSNORM
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,INTENT(OUT)                       :: nVec(3), tang1(3), tang2(3)
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: a,b
REAL                                   :: nlength
!================================================================================================================================


b=xi*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)  & 
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
   +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)   &
          +BezierControlPoints3D(:,NGeo,NGeo,SideID)+BezierControlPoints3D(:,0   ,NGeo,SideID) )

a=eta*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
    +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
           +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) )

tang1=a/DOT_PRODUCT(a,a)
tang2=b/DOT_PRODUCT(b,b)
nVec=CROSSNORM(tang1,tang2)

END SUBROUTINE CalcBiLinearNormAndTang


SUBROUTINE CalcNormAndTangBezier(nVec,tang1,tang2,xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Globals,                              ONLY:CROSSNORM,CROSS
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D,facNchooseK,ArrayNchooseK,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,               ONLY:SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,DIMENSION(3),INTENT(OUT)          :: nVec,tang1,tang2
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: v,u!,nVec
INTEGER                                :: p,q,m
!REAL                                   :: nlength
REAL                                   :: MinusXi,PlusXI,MinusEta,PlusEta
REAL                                   :: xiup(0:NGeo),etaup(0:NGeo),xidown(0:NGeo),etadown(0:NGeo)
!================================================================================================================================

IF(BoundingBoxIsEmpty(SideID))THEN
  nVec=SideNormVec(1:3,SideID)
  u=xi*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)  & 
            +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
     +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0  ,SideID)   &
            +BezierControlPoints3D(:,NGeo,NGeo,SideID)+BezierControlPoints3D(:,0   ,NGeo,SideID) )
  
  v=eta*0.25*(BezierControlPoints3D(:,0   ,0   ,SideID)-BezierControlPoints3D(:,NGeo,0   ,SideID)   &
             +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) ) &
      +0.25*(-BezierControlPoints3D(:,0   ,0   ,SideID)+BezierControlPoints3D(:,NGeo,0   ,SideID)   &
             +BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0   ,NGeo,SideID) )
  
  tang1=u/DOT_PRODUCT(u,u)
  tang2=v/DOT_PRODUCT(v,v)
ELSE ! no flat side
  ! compute norm vec
  ! caution we require the formula in [0;1]
  M=nGeo-1
  MinusXi=1.0-xi
  PlusXi=1.0+xi
  !PlusXi=xi
  MinusEta=1.0-eta
  PlusEta=1.0+eta
  !PlusEta=eta
  
  !! compute the required stuff
  xiup(0)=1.0
  etaup(0)=1.0
  ! caution, here the indicies are switched from n-j to j 
  xidown(NGeo)=1.0
  etadown(NGeo)=1.0
  DO p=1,NGeo
    xiup(p)=xiup(p-1)*PlusXi
    etaup(p)=etaup(p-1)*PlusEta
    xidown(NGeo-p)=xidown(NGeo-p+1)*MinusXi
    etadown(NGeo-p)=etadown(NGeo-p+1)*MinusEta
  END DO ! p
  
  u=0
  v=0
  DO q=0,NGeo
    DO p=0,M
      ! derivative in xi
      u=u+(BezierControlPoints3D(:,p+1,q,SideID)-BezierControlPoints3D(:,p,q,SideID)) &
          *facNchooseK(M,p)*xiup(p)*xidown(p)     &
          *facNChooseK(NGeo,q)*etaup(q)*etadown(q)
      ! derivative in eta ! caution - exchange indicies
      v=v+(BezierControlPoints3D(:,q,p+1,SideID)-BezierControlPoints3D(:,q,p,SideID)) &
          *facNchooseK(NGeo,q)*xiup(q)*xidown(q)  &
          *facNChooseK(M,p)*etaup(p)*etadown(p)
  
    END DO ! p
  END DO ! q

  tang1=u/DOT_PRODUCT(u,u)
  tang2=v/DOT_PRODUCT(v,v)
  nVec =CROSSNORM(u,v)
END IF ! BoundingBoxIsEmpty

END SUBROUTINE CalcNormAndTangBezier


FUNCTION CalcNormVecBezier(xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Globals,                              ONLY:CROSSNORM,CROSS
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D,facNchooseK,ArrayNchooseK,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,               ONLY:SideNormVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                        :: xi,eta
INTEGER,INTENT(IN)                     :: SideID
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
REAL,DIMENSION(3)                      :: CalcNormVecBezier
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)                      :: v,u!,nVec
INTEGER                                :: p,q,m
!REAL                                   :: nlength
REAL                                   :: MinusXi,PlusXI,MinusEta,PlusEta
REAL                                   :: xiup(0:NGeo),etaup(0:NGeo),xidown(0:NGeo),etadown(0:NGeo)
!================================================================================================================================

IF(BoundingBoxIsEmpty(SideID))THEN
  CalcNormVecBezier=SideNormVec(1:3,SideID)
ELSE ! no flat side
  ! compute norm vec
  ! caution we require the formula in [0;1]
  M=nGeo-1
  MinusXi=1.0-xi
  PlusXi=1.0+xi
  !PlusXi=xi
  MinusEta=1.0-eta
  PlusEta=1.0+eta
  !PlusEta=eta
  
  !! compute the required stuff
  xiup(0)=1.0
  etaup(0)=1.0
  ! caution, here the indicies are switched from n-j to j 
  xidown(NGeo)=1.0
  etadown(NGeo)=1.0
  DO p=1,NGeo
    xiup(p)=xiup(p-1)*PlusXi
    etaup(p)=etaup(p-1)*PlusEta
    xidown(NGeo-p)=xidown(NGeo-p+1)*MinusXi
    etadown(NGeo-p)=etadown(NGeo-p+1)*MinusEta
  END DO ! p
  
  ! B = (1./(2**N_in))*REAL(CHOOSE(N_in,j))*((x+1.)**j)*((1.-x)**(N_in-j))
  ! complete form
  !u=0.
  !v=0.
  !DO q=0,NGeo
  !  DO p=0,M
  !    ! derivative in xi
  !!   BezierControlPoints2D_temp(:,q,p)=&
  !!   BezierControlPoints2D_temp(:,q,p)+&
  !!   !BezierControlPoints2D(:,l,q)*B(p,l,Smax)
  !!   BezierControlPoints2D     (:,q,l)*(1./(2.**p))       &
  !!                         *arrayNchooseK(p,l) &
  !!                         *(1.+Etamax)**l       &
  !!                         *(1.-Etamax)**(p-l)
  !
  !    u=u+(BezierControlPoints3D(:,p+1,q,SideID)-BezierControlPoints3D(:,p,q,SideID))        &
  !        *facNchooseK(M,p)*(PlusXi**p)*(MinusXi**(M-p)) &
  !        *facNChooseK(NGeo,q)*(PlusEta**q)*(MinusEta**(NGeo-q))
  !
  !    v=v+(BezierControlPoints3D(:,q,p+1,SideID)-BezierControlPoints3D(:,q,p,SideID))        &
  !        *facNchooseK(NGeo,q)*(PlusXi**q)*(MinusXi**(NGeo-q)) &
  !                                            *facNChooseK(M,p)*(PlusEta**p)*(MinusEta**(M-p))
  !
  !!    u=u+(BezierControlPoints3D(:,p+1,q,SideID)-BezierControlPoints3D(:,p,q,SideID))        &
  !!        *BezierControlPoints3D(:,p,q,SideID)*facNchooseK(M,p)*(PlusXi**p)*(MinusXi**(M-p)) &
  !!                                            *facNChooseK(NGeo,q)*(PlusEta**q)*(MinusEta**(NGeo-q))
  !
  !!    v=v+(BezierControlPoints3D(:,q,p+1,SideID)-BezierControlPoints3D(:,q,p,SideID))        &
  !!        *BezierControlPoints3D(:,q,p,SideID)*facNchooseK(NGeo,q)*(PlusXi**q)*(MinusXi**(NGeo-q)) &
  !!                                            *facNChooseK(M,p)*(PlusEta**p)*(MinusEta**(M-p))
  !
  !  END DO ! p
  !END DO ! q
  !u=u!*NGeo
  !v=v!*NGeo
  
  u=0
  v=0
  DO q=0,NGeo
    DO p=0,M
      ! derivative in xi
      u=u+(BezierControlPoints3D(:,p+1,q,SideID)-BezierControlPoints3D(:,p,q,SideID)) &
          *facNchooseK(M,p)*xiup(p)*xidown(p)     &
          *facNChooseK(NGeo,q)*etaup(q)*etadown(q)
      ! derivative in eta ! caution - exchange indicies
      v=v+(BezierControlPoints3D(:,q,p+1,SideID)-BezierControlPoints3D(:,q,p,SideID)) &
          *facNchooseK(NGeo,q)*xiup(q)*xidown(q)  &
          *facNChooseK(M,p)*etaup(p)*etadown(p)
  
    END DO ! p
  END DO ! q
  !u=u*NGeo
  !v=v*NGeo
  
  CalcNormVecBezier=CROSSNORM(u,v)
  !nVec=CROSS(u,v)
  !nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2)+nVec(3)*nVec(3)
  !nlength=SQRT(nlength)
  !CalcNormVecBezier=nVec/nlength
END IF ! BoundingBoxIsEmpty

END FUNCTION CalcNormVecBezier


SUBROUTINE GetBezierControlPoints3D(XCL_NGeo,iElem)
!===================================================================================================================================
! computes the nodes for Bezier Control Points for [P][I][C] [A]daptive [S]uper [S]ampled Surfaces [O]perations
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:nSides,ElemToSide,SideToElem,NGeo
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,sVdm_Bezier
USE MOD_Particle_Tracking_Vars,   ONLY:DoRefMapping
USE MOD_Mesh_Vars,                ONLY:nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_ChangeBasis,              ONLY:ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: lastSideID,SideID
INTEGER                           :: p,q
REAL                              :: tmp(3,0:NGeo,0:NGeo)  

!===================================================================================================================================
!print*,"SIZE(BezierControlPoints)"
!print*,SIZE(BezierControlPoints)
!print*,"SHAPE(BezierControlPoints)"
!print*,SHAPE(BezierControlPoints)
! BCSides, InnerSides and MINE MPISides are filled
lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
!IF(DoRefMapping) lastSideID  = nBCSides
!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) XI_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,XI_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
    CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,0,:,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) XI_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,XI_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:),tmp)
    !print*,'ixi'
    BezierControlPoints3D(:,:,:,SideID)=tmp
    CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,NGeo,:,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) ETA_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:),BezierControlPoints3D(1:3,:,:,sideID))
      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
   END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 4.) ETA_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ETA_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,NGeo-p,q)
      END DO !p
    END DO !q
      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,NGeo,:),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 5.) ZETA_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ZETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0),tmp)
    ! turn into right hand system of side
    DO q=0,NGeo
      DO p=0,NGeo
        BezierControlPoints3D(1:3,p,q,sideID)=tmp(:,q,p)
      END DO !p
    END DO !q
      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,0),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 6.) ZETA_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo),BezierControlPoints3D(1:3,:,:,sideID))
    CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
  END IF !flip=0
END IF
!ELSE ! no master, here has to come the suff with the slave
!  CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo),tmp)
!  flip= SideToElem(S2E_FLIP,SideID)
!  SELECT CASE(flip)
!    CASE(1) ! slave side, SideID=q,jSide=p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,q,p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSlabNormalsAndIntervalls(NGeo,SideID)
!  END SELECT
!END IF
END SUBROUTINE GetBezierControlPoints3D

SUBROUTINE GetSlabNormalsAndIntervalls(NGeo,SideID)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervalls)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,    ONLY:EpsMach
USE MOD_Preproc
USE MOD_Particle_Surfaces_Vars,   ONLY:SlabNormals,SlabIntervalls,BezierControlPoints3D,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilonbilinear
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: SideID,NGeo
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER                           :: lastSideID,flip,SideID
INTEGER            :: p,q
!REAL                              :: tmp(3,0:NGeo,0:NGeo)  
REAL               :: skalprod(3),dx,dy,dz,dMax,dMin,w,h,l
LOGICAL            :: SideIsPlanar
LOGICAL            :: SideIsCritical
!===================================================================================================================================


!BezierControlPoints(:,:,:,SideID)
!SlabNormals( x y z,1 2 3 , SideID)

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------
!SideIsPlanar=.FALSE.

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors
!-----------------------------------------------------------------------------------------------------------------------------------
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
SlabNormals(:,1,SideID)=BezierControlPoints3D(:,NGeo,0,SideID)   -BezierControlPoints3D(:,0,0,SideID)+&
                        BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,0,NGeo,SideID)
SlabNormals(:,1,SideID)=SlabNormals(:,1,SideID)/SQRT(DOT_PRODUCT(SlabNormals(:,1,SideID),SlabNormals(:,1,SideID)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
SlabNormals(:,2,SideID)=BezierControlPoints3D(:,0,NGeo,SideID)   -BezierControlPoints3D(:,0,0,SideID)+&
                        BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,NGeo,0,SideID)
SlabNormals(:,2,SideID)=CROSSNORM(SlabNormals(:,1,SideID),SlabNormals(:,2,SideID))


!b                      =BezierControlPoints3D(:,0,NGeo,SideID)   -BezierControlPoints3D(:,0,0,SideID)+&
                        !BezierControlPoints3D(:,NGeo,NGeo,SideID)-BezierControlPoints3D(:,NGeo,0,SideID)
!SlabNormals(2,1,SideID)=SlabNormals(1,2,SideID)*b(3) - SlabNormals(1,3,SideID)*b(2)
!SlabNormals(2,2,SideID)=SlabNormals(1,3,SideID)*b(1) - SlabNormals(1,1,SideID)*b(3)
!SlabNormals(2,3,SideID)=SlabNormals(1,1,SideID)*b(2) - SlabNormals(1,2,SideID)*b(1)
!SlabNormals(2,:,SideID)=SlabNormals(2,:,SideID)/SQRT(DOT_PRODUCT(SlabNormals(2,:,SideID),SlabNormals(2,:,SideID)))

!SlabNormals(1,2,SideID)=SlabNormals(2,1,SideID)*(                                                     &
!                        BezierControlPoints3D(3,0,NGeo,SideID)   -BezierControlPoints3D(3,0,0,SideID)+    &
!                        BezierControlPoints3D(3,NGeo,NGeo,SideID)-BezierControlPoints3D(3,NGeo,0,SideID)) &
!                       -SlabNormals(3,1,SideID)*(                                                     &
!                        BezierControlPoints3D(2,0,NGeo,SideID)   -BezierControlPoints3D(2,0,0,SideID)+    &
!                        BezierControlPoints3D(2,NGeo,NGeo,SideID)-BezierControlPoints3D(2,NGeo,0,SideID)) 
!SlabNormals(2,2,SideID)=SlabNormals(3,1,SideID)*(                                                     &
!                        BezierControlPoints3D(1,0,NGeo,SideID)   -BezierControlPoints3D(1,0,0,SideID)+    &
!                        BezierControlPoints3D(1,NGeo,NGeo,SideID)-BezierControlPoints3D(1,NGeo,0,SideID)) &
!                       -SlabNormals(1,1,SideID)*(                                                     &
!                        BezierControlPoints3D(3,0,NGeo,SideID)   -BezierControlPoints3D(3,0,0,SideID)+    &
!                        BezierControlPoints3D(3,NGeo,NGeo,SideID)-BezierControlPoints3D(3,NGeo,0,SideID)) 
!SlabNormals(3,2,SideID)=SlabNormals(1,1,SideID)*(                                                     &
!                        BezierControlPoints3D(2,0,NGeo,SideID)   -BezierControlPoints3D(2,0,0,SideID)+    &
!                        BezierControlPoints3D(2,NGeo,NGeo,SideID)-BezierControlPoints3D(2,NGeo,0,SideID)) &
!                       -SlabNormals(2,1,SideID)*(                                                     &
!                        BezierControlPoints3D(1,0,NGeo,SideID)   -BezierControlPoints3D(1,0,0,SideID)+    &
!                        BezierControlPoints3D(1,NGeo,NGeo,SideID)-BezierControlPoints3D(1,NGeo,0,SideID))
!SlabNormals(:,2,SideID)=SlabNormals(:,2,SideID)/SQRT(DOT_PRODUCT(SlabNormals(:,2,SideID),SlabNormals(:,2,SideID)))
! n_3=n_1 x n_2
SlabNormals(:,3,SideID)=CROSSNORM(SlabNormals(:,2,SideID),SlabNormals(:,1,SideID))
!SlabNormals(1,3,SideID)=SlabNormals(2,2,SideID)*SlabNormals(3,2,SideID) - SlabNormals(3,2,SideID)*SlabNormals(2,1,SideID)
!SlabNormals(2,3,SideID)=SlabNormals(3,2,SideID)*SlabNormals(1,2,SideID) - SlabNormals(1,2,SideID)*SlabNormals(3,1,SideID)
!SlabNormals(3,3,SideID)=SlabNormals(1,2,SideID)*SlabNormals(2,2,SideID) - SlabNormals(2,2,SideID)*SlabNormals(1,1,SideID)
!SlabNormals(3,:,SideID)=SlabNormals(3,:,SideID)/SQRT(DOT_PRODUCT(SlabNormals(3,:,SideID),SlabNormals(3,:,SideID)))
!print*,"slab normal vector length: ",SQRT(DOT_PRODUCT(SlabNormals(1,:,SideID),SlabNormals(1,:,SideID))),&
!                                     SQRT(DOT_PRODUCT(SlabNormals(2,:,SideID),SlabNormals(2,:,SideID))),&
!                                     SQRT(DOT_PRODUCT(SlabNormals(3,:,SideID),SlabNormals(3,:,SideID)))
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervalls beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!SlabIntervalls(x- x+ y- y+ z- z+, SideID)
SideIsPlanar=.FALSE.
SideIsCritical=.FALSE.
! Intervall beta_1
!print*,"SideID",SideID
SlabIntervalls(:, SideID)=0.
! DEBUGGG output
!IF(SideID.EQ.2)THEN
!print*,'SideID',SideID
!DO q=0,NGeo
!  DO p=0,NGeo
!    IF((p.EQ.0).AND.(q.EQ.0))CYCLE
!    skalprod(1)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,1,SideID))
!    skalprod(2)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,2,SideID))
!    skalprod(3)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,3,SideID))
!    print*,'skalarprod1 ', skalprod(1)
!    print*,'skalarprod2 ', skalprod(2)
!    print*,'skalarprod3 ', skalprod(3)
!    IF    (skalprod(1).LT.0.)THEN
!      SlabIntervalls(1, SideID)=MIN(SlabIntervalls(1,SideID),skalprod(1))
!      print*,'SlabIntervals(1)',SlabIntervalls(1,SideID)
!    ELSEIF(skalprod(1).GT.0.)THEN
!      SlabIntervalls(2, SideID)=MAX(SlabIntervalls(2,SideID),skalprod(1))
!      print*,'SlabIntervals(2)',SlabIntervalls(2,SideID)
!    END IF
!    IF    (skalprod(2).LT.0.)THEN
!      SlabIntervalls(3, SideID)=MIN(SlabIntervalls(3,SideID),skalprod(2))
!      print*,'SlabIntervals(3)',SlabIntervalls(3,SideID)
!    ELSEIF(skalprod(2).GT.0.)THEN
!      SlabIntervalls(4, SideID)=MAX(SlabIntervalls(4,SideID),skalprod(2))
!      print*,'SlabIntervals(4)',SlabIntervalls(4,SideID)
!    END IF
!    IF    (skalprod(3).LT.0.)THEN
!      SlabIntervalls(5, SideID)=MIN(SlabIntervalls(5,SideID),skalprod(3))
!      print*,'SlabIntervals(5)',SlabIntervalls(5,SideID)
!    ELSEIF(skalprod(3).GT.0.)THEN
!      SlabIntervalls(6, SideID)=MAX(SlabIntervalls(6,SideID),skalprod(3))
!      print*,'SlabIntervals(6)',SlabIntervalls(6,SideID)
!    END IF
!    read*
!  END DO !p
!END DO !q
!END IF

DO q=0,NGeo
  DO p=0,NGeo
    IF((p.EQ.0).AND.(q.EQ.0))CYCLE
    skalprod(1)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,1,SideID))
    skalprod(2)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,2,SideID))
    skalprod(3)=DOT_PRODUCT(BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,0,SideID),SlabNormals(:,3,SideID))
    IF    (skalprod(1).LT.0.)THEN
      SlabIntervalls(1, SideID)=MIN(SlabIntervalls(1,SideID),skalprod(1))
    ELSEIF(skalprod(1).GT.0.)THEN
      SlabIntervalls(2, SideID)=MAX(SlabIntervalls(2,SideID),skalprod(1))
    END IF
    IF    (skalprod(2).LT.0.)THEN
      SlabIntervalls(3, SideID)=MIN(SlabIntervalls(3,SideID),skalprod(2))
    ELSEIF(skalprod(2).GT.0.)THEN
      SlabIntervalls(4, SideID)=MAX(SlabIntervalls(4,SideID),skalprod(2))
    END IF
    IF    (skalprod(3).LT.0.)THEN
      SlabIntervalls(5, SideID)=MIN(SlabIntervalls(5,SideID),skalprod(3))
    ELSEIF(skalprod(3).GT.0.)THEN
      SlabIntervalls(6, SideID)=MAX(SlabIntervalls(6,SideID),skalprod(3))
    END IF
  END DO !p
END DO !q
!print*,"SideID",SideID," is planar:",SideIsPlanar(SideID)
dx=ABS(ABS(SlabIntervalls(2, SideID))-ABS(SlabIntervalls(1, SideID)))
dy=ABS(ABS(SlabIntervalls(4, SideID))-ABS(SlabIntervalls(3, SideID)))
dz=ABS(ABS(SlabIntervalls(6, SideID))-ABS(SlabIntervalls(5, SideID)))
dMax=MAX(dx,dy,dz)
dMin=MIN(dx,dy,dz)
!print*,dx
!IF(dx/dMax.LT.1.E-3)THEN
IF(dx/dMax.LT.epsilonbilinear)THEN
  SlabIntervalls(1:2, SideID)=0.
  dx=0.
END IF
!print*,dy
!IF(dy/dMax.LT.1.E-3)THEN
IF(dy/dMax.LT.epsilonbilinear)THEN
  SlabIntervalls(3:4, SideID)=0.
  dy=0.
END IF
!print*,dz
!IF(dz/dMax.LT.1.E-3)THEN
IF(dz/dMax.LT.epsilonbilinear)THEN
  SlabIntervalls(5:6, SideID)=0.
  dz=0.
END IF

IF(dx*dy*dz.LT.0) THEN
  IPWRITE(UNIT_stdOut,*) ' Warning, no bounding box'
END IF

!IF(dx*dy*dz.EQ.0.)THEN
IF(ALMOSTZERO(dx*dy*dz))THEN
!IF(ABS(dx*dy*dz).LT.2.0*epsMach)THEN
  SideIsPlanar=.TRUE.
  BoundingBoxIsEmpty(SideID)=.TRUE.
ELSE
  BoundingBoxIsEmpty(SideID)=.FALSE.
END IF

! print*,"SideIsPlanar(",SideID,")",SideIsPlanar
! print*,"BezierControlPoints3D(:,0,0,SideID)",BezierControlPoints3D(:,0,0,SideID)
! print*,"SlabNormals(:,1,SideID)",SlabNormals(:,1,SideID),"beta1 ",SlabIntervalls(1, SideID),SlabIntervalls(2,SideID),&
! "delta",ABS(SlabIntervalls(2, SideID))-ABS(SlabIntervalls(1, SideID))
! print*,"SlabNormals(:,2,SideID)",SlabNormals(:,2,SideID),"beta2 ",SlabIntervalls(3, SideID),SlabIntervalls(4, SideID),&
! "delta",ABS(SlabIntervalls(4, SideID))-ABS(SlabIntervalls(3, SideID))
! print*,"SlabNormals(:,3,SideID)",SlabNormals(:,3,SideID),"beta3 ",SlabIntervalls(5, SideID),SlabIntervalls(6, SideID),&
! "delta",ABS(SlabIntervalls(6, SideID))-ABS(SlabIntervalls(5, SideID))
!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) Is Side critical? (particle path parallel to the larger surface, therefore numerous intersections are possilbe)
!-----------------------------------------------------------------------------------------------------------------------------------
RETURN


IF(.NOT.SideIsPlanar)THEN
  IF(dx.EQ.dMin)THEN
    h=dx
    IF(dy.GE.dz)THEN
      l=dy
      w=dz
    ELSE
      l=dz
      w=dy
    END IF
  END IF

  IF(dy.EQ.dMin)THEN
    h=dy
    IF(dx.GE.dz)THEN
      l=dx
      w=dz
    ELSE
      l=dz
      w=dx
    END IF
  END IF

  IF(dz.EQ.dMin)THEN
  h=dz
    IF(dy.GE.dz)THEN
      l=dy
      w=dz
    ELSE
      l=dz
      w=dy
    END IF
  END IF
END IF

END SUBROUTINE GetSlabNormalsAndIntervalls


SUBROUTINE GetBCSideType()
!================================================================================================================================
! select the side type for each side 
! check if points on edges are linear. if linear, the cross product of the vector between two vertices and a vector between a 
! vercites and a edge point has to be zero
! SideType
! 0 - planar
! 1 - bilinear
! 2 - curved
!================================================================================================================================
USE MOD_Globals!,                  ONLY:CROSS
USE MOD_Mesh_Vars,                ONLY:nSides,NGeo,Xi_NGeo,Sideid_minus_upper,nBCSides
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,BoundingBoxIsEmpty,epsilontol,SideType,SideNormVec,SideDistance
USE MOD_Particle_Mesh_Vars,       ONLY:nTotalSides,IsBCElem,nTotalBCSides,nTotalElems,nTotalBCElems
USE MOD_Particle_Mesh_Vars,       ONLY:PartElemToSide,BCElem,PartSideToElem
USE MOD_Particle_Mesh_Vars,       ONLY:PartBCSideList,nTotalBCSides
USE MOD_Particle_MPI_Vars,        ONLY:halo_eps,halo_eps2
#ifdef MPI
USE MOD_Particle_MPI_HALO,        ONLY:WriteParticleMappingPartitionInformation
#endif /*MPI*/
#if ((PP_TimeDiscMethod!=1) && (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6))  /* RK3 and RK4 only */
USE MOD_Mesh_Vars,                ONLY:XCL_NGeo
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iSide,p,q, nPlanar,nBilinear,nCurved,nDummy,SideID,iElem,ilocSide,nBCElems
INTEGER                     :: nSideCount, BCSideID, BCSideID2, s,r
INTEGER,ALLOCATABLE         :: SideIndex(:)
REAL,DIMENSION(1:3)         :: v1,v2,NodeX
REAL                        :: length,eps
LOGICAL                     :: isLinear,leave
#ifdef MPI  
INTEGER                     :: nPlanarTot,nBilinearTot,nCurvedTot,nBCElemsTot
#endif /*MPI*/
#if ((PP_TimeDiscMethod!=1) && (PP_TimeDiscMethod!=2) && (PP_TimeDiscMethod!=6))  /* RK3 and RK4 only */
REAL,DIMENSION(1:3,0:NGeo,0:NGeo) :: xNodes
#endif
!================================================================================================================================

! allocated here,!!! should be moved?
!ALLOCATE( SideType(nSides)        &
!        , SideDistance(nSides)    &
!        , SideNormVec(1:3,nSides) )

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Get Side Type incl. HALO-Sides...'

ALLOCATE( SideType(nTotalBCSides)        &
        , SideDistance(nTotalBCSides)    &
        , isBCElem(nTotalElems)          &
        , SideNormVec(1:3,nTotalBCSides) )

SideDistance=0.
SideNormVec=0.

eps=1e-8
nPlanar=0
nBilinear=0
nCurved=0
nBCElems=0
#ifdef MPI
nPlanarTot=0
nBilinearTot=0
nCurvedTot=0
nBCElemsTot=0
#endif /*MPI*/
DO iSide=1,nTotalSides
  isLinear=.TRUE.
  SideID  =PartBCSideList(iSide)
  IF(SideID.EQ.-1) CYCLE
  ! all four edges
  !IF(iSide.GT.nSides) IPWRITE(UNIT_stdOut,*) BezierControlPOints3D(:,:,:,iSide)
  IF(SUM(ABS(BezierControlPoints3D(:,:,:,SideID))).LT.1e-10) IPWRITE(UNIT_stdOut,*) 'missing side',SideID
  q=0
  v1=BezierControlPoints3D(:,NGeo,q,SideID)-BezierControlPoints3D(:,0,q,SideID)
  DO p=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,q,SideID)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! p
  q=NGeo
  v1=BezierControlPoints3D(:,NGeo,q,SideID)-BezierControlPoints3D(:,0,q,SideID)
  DO p=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,0,q,SideID)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! p
  p=0
  v1=BezierControlPoints3D(:,p,NGeo,SideID)-BezierControlPoints3D(:,p,0,SideID)
  DO q=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,p,0,SideID)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! q
  p=NGeo
  v1=BezierControlPoints3D(:,p,NGeo,SideID)-BezierControlPoints3D(:,p,0,SideID)
  DO q=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,SideID)-BezierControlPoints3D(:,p,0,SideID)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! q
  IF(isLinear)THEN
    IF(BoundingBoxIsEmpty(SideID))THEN
      SideType(SideID)=PLANAR
      IF(SideID.LE.SideID_Minus_Upper) nPlanar=nPlanar+1
#ifdef MPI
      IF(SideID.GT.SideID_Minus_Upper) nPlanartot=nPlanartot+1
#endif /*MPI*/
      ! compute the norm vec of side and distance from origin
      v1=BezierControlPoints3D(:,NGeo,0,SideID)-BezierControlPoints3D(:,0,0,SideID)
      v2=BezierControlPoints3D(:,0,NGeo,SideID)-BezierControlPoints3D(:,0,0,SideID)
      SideNormVec(:,SideID) = CROSSNORM(v1,v2)
!      length=SQRT(v2(1)*v2(1)+v1(2)*v1(2)+v1(3)*v1(3))
!      SideNormVec(:,iSide) =SideNormVec(:,iSide)/length
      !      xNodes(:,1)=SuperSampledNodes(1:3,p  ,q  ,SideID)
      !      xNodes(:,2)=SuperSampledNodes(1:3,p+1,q  ,SideID)
      !      xNodes(:,3)=SuperSampledNodes(1:3,p+1,q+1,SideID)
      !      xNodes(:,4)=SuperSampledNodes(1:3,p  ,q+1,SideID)
      v1=0.25*(BezierControlPoints3D(:,0,0,SideID)     &
              +BezierControlPoints3D(:,NGeo,0,SideID)  &
              +BezierControlPoints3D(:,0,NGeo,SideID)  &
              +BezierControlPoints3D(:,NGeo,NGeo,SideID))
      SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
!      IF(SideID.EQ.8)THEN
!        print*,'v1',v1
!        print*,'SideDistance',SideDistance(SideID)
!        read*
!      END IF
    ELSE
      SideType(SideID)=BILINEAR
      IF(SideID.LE.SideID_Minus_Upper) nBiLinear=nBiLinear+1
#ifdef MPI
      IF(SideID.GT.SideID_Minus_Upper) nBilinearTot=nBilinearTot+1
#endif /*MPI*/
    END IF ! BoundingBoxIsEmpty
  ELSE ! non-linear sides
    SideType(SideID)=CURVED
    IF(SideID.LE.SideID_Minus_Upper) nCurved=nCurved+1
#ifdef MPI
    IF(SideID.GT.SideID_Minus_Upper) nCurvedTot=nCurvedTot+1
#endif /*MPI*/
    IF(BoundingBoxIsEmpty(SideID))THEN
      v1=BezierControlPoints3D(:,NGeo,0,SideID)-BezierControlPoints3D(:,0,0,SideID)
      v2=BezierControlPoints3D(:,0,NGeo,SideID)-BezierControlPoints3D(:,0,0,SideID)
      SideNormVec(:,SideID) = CROSSNORM(v1,v2)
      v1=0.25*(BezierControlPoints3D(:,0,0,SideID)     &
              +BezierControlPoints3D(:,NGeo,0,SideID)  &
              +BezierControlPoints3D(:,0,NGeo,SideID)  &
              +BezierControlPoints3D(:,NGeo,NGeo,SideID))
      SideDistance(SideID)=DOT_PRODUCT(v1,SideNormVec(:,SideID))
    END IF ! BoundingBoxIsEmpty
  END IF ! isLinear
END DO ! iSide

! mark elements as bc element
IsBCElem=.FALSE.
nTotalBCElems=0
DO iElem=1,nTotalElems
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    IF (SideID.EQ.-1) CYCLE
    IF((SideID.LE.nBCSides).OR.(SideID.GT.nSides))THEN
      IF(.NOT.isBCElem(iElem))THEN
        IsBCElem(iElem)=.TRUE.
        nTotalBCElems=nTotalBCElems+1
        IF(SideID.LE.nBCSides)THEN
          nBCElems=nBCElems+1
        END IF
      END IF ! count only single
    END IF
  END DO ! ilocSide
END DO ! iElem

! build list with elements in halo-eps vicinity around bc-elements
ALLOCATE( BCElem(1:nTotalElems) )
ALLOCATE( SideIndex(1:nTotalSides) )
! number of element local BC-Sides
DO iElem=1,nTotalElems
  BCElem(iElem)%nInnerSides=0
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* only LSERK */
  IF(.NOT.isBCElem(iElem)) CYCLE
#endif
  DO ilocSide=1,6
    SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    IF(SideID.EQ.-1) CYCLE
    IF(PartBCSideList(SideID).EQ.-1) CYCLE
    BCElem(iElem)%nInnerSides = BCElem(iElem)%nInnerSides+1
    !END IF
  END DO ! ilocSide
  BCElem(iElem)%lastSide=BCElem(iElem)%nInnerSides
  ! loop over all sides
  SideIndex=0
  DO iSide=1,nTotalSides
    ! ignore sides of the same element
    IF(PartSideToElem(S2E_ELEM_ID,iSide).EQ.iElem) CYCLE
    BCSideID  =PartBCSideList(iSide)
    IF(BCSideID.EQ.-1) CYCLE
    ! next, get all sides in halo-eps vicinity
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(SideID.EQ.-1) CYCLE
#if ((PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==6))  /* only LSERK */
      BCSideID2 =PartBCSideList(SideID)
      IF(BCSideID2.EQ.-1) CYCLE
      leave=.FALSE.
      ! all points of bc side
      DO q=0,NGeo
        DO p=0,NGeo
          NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)
          ! all nodes of current side
          DO s=0,NGeo
            DO r=0,NGeo
              IF(SQRT(DOT_Product(BezierControlPoints3D(:,r,s,BCSideID2)-NodeX &
                                 ,BezierControlPoints3D(:,r,s,BCSideID2)-NodeX )).LE.halo_eps)THEN
                IF(SideIndex(iSide).EQ.0)THEN
                  BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
                  SideIndex(iSide)=BCElem(iElem)%lastSide
                  leave=.TRUE.
                  EXIT
                END IF
              END IF
            END DO ! r
            IF(leave) EXIT
          END DO ! s
          IF(leave) EXIT
        END DO ! p
        IF(leave) EXIT
      END DO ! q
#else /* no LSERK */
      SELECT CASE(ilocSide)
      CASE(XI_MINUS)
        xNodes=XCL_NGeo(1:3,0,0:NGeo,0:NGeo,iElem)
      CASE(XI_PLUS)
        xNodes=XCL_NGeo(1:3,NGeo,0:NGeo,0:NGeo,iElem)
      CASE(ETA_MINUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0,0:NGeo,iElem)
      CASE(ETA_PLUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,NGeo,0:NGeo,iElem)
      CASE(ZETA_MINUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,0,iElem)
      CASE(ZETA_PLUS)
        xNodes=XCL_NGeo(1:3,0:NGeo,0:NGeo,NGeo,iElem)
      END SELECT
      leave=.FALSE.
      ! all points of bc side
      DO q=0,NGeo
        DO p=0,NGeo
          NodeX(:) = BezierControlPoints3D(:,p,q,BCSideID)
          ! all nodes of current side
          DO s=0,NGeo
            DO r=0,NGeo
              IF(SQRT(DOT_Product(xNodes(:,r,s)-NodeX &
                                 ,xNodes(:,r,s)-NodeX )).LE.halo_eps)THEN
                IF(SideIndex(iSide).EQ.0)THEN
                  BCElem(iElem)%lastSide=BCElem(iElem)%lastSide+1
                  SideIndex(iSide)=BCElem(iElem)%lastSide
                  leave=.TRUE.
                  EXIT
                END IF
              END IF
            END DO ! r
            IF(leave) EXIT
          END DO ! s
          IF(leave) EXIT
        END DO ! p
        IF(leave) EXIT
      END DO ! q
#endif
      IF(leave) EXIT
    END DO ! ilocSide
  END DO ! iSide
  ! finally, allocate the bc side list
  IF(BCElem(iElem)%lastSide.EQ.0) CYCLE
  ! set true, only required for elements without an own bc side
  IF(.NOT.isBCElem(iElem)) nBCElems=nBCElems+1
  isBCElem(iElem)=.TRUE.
  ! allocate complete side list
  ALLOCATE( BCElem(iElem)%BCSideID(BCElem(iElem)%lastSide) )
  ! 1) inner sides
  nSideCount=0
  IF(BCElem(iElem)%nInnerSides.GT.0)THEN
    DO ilocSide=1,6
      SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(SideID.EQ.-1) CYCLE
      BCSideID=PartBCSideList(SideID)
      IF(BCSideID.EQ.-1) CYCLE
     ! IF((SideID.LE.nBCSides).OR.(SideID.GT.nSides))THEN
      nSideCount=nSideCount+1
      !BCElem(iElem)%BCSideID(nSideCount)= BCSideID
      BCElem(iElem)%BCSideID(nSideCount)= SideID
      !END IF
    END DO ! ilocSide
  END IF ! nInnerSides.GT.0
  ! 2) outer sides
  DO iSide=1,nTotalSides
    IF(SideIndex(iSide).GT.0)THEN
      nSideCount=nSideCount+1
      BCSideID=PartBCSideList(iSide)
      !BCElem(iElem)%BCSideID(nSideCount)=BCSideID
      BCElem(iElem)%BCSideID(nSideCount)=iSide
    END IF
  END DO  ! iSide
END DO ! iElem

#ifdef MPI
nPlanarTot=nPlanar+nPlanarTot
nBilinearTot=nBilinear+nBilinearTot
nCurvedTot=nCurved+nCurvedTot
IF(MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nPlanar  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nBilinear,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nCurved  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nBCElems ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
ELSE ! no Root
  CALL MPI_REDUCE(nPlanar  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nBilinear,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurved  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nBCElems  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
END IF
#endif /*MPI*/


SWRITE(UNIT_StdOut,'(A,I8)') ' Number of BC-adjoined elems: ', nBCElems
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar      faces: ', nPlanar
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear   faces: ', nBilinear
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved      faces: ', nCurved

#ifdef MPI
CALL  WriteParticleMappingPartitionInformation(nPlanarTot,nBilinearTot,nCurvedTot,nTotalBCElems)
#endif

END SUBROUTINE GetBCSideType


SUBROUTINE GetSideType()
!================================================================================================================================
! select the side type for each side 
! check if points on edges are linear. if linear, the cross product of the vector between two vertices and a vector between a 
! vercites and a edge point has to be zero
! SideType
! 0 - planar
! 1 - bilinear
! 2 - curved
!================================================================================================================================
USE MOD_Globals!,                  ONLY:CROSS
USE MOD_Mesh_Vars,                ONLY:nSides,NGeo,Xi_NGeo,Sideid_minus_upper
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,BoundingBoxIsEmpty,epsilontol,SideType,SideNormVec,SideDistance
!USE MOD_Particle_Surfaces_Vars,   ONLY:epsilonbilinear
USE MOD_Particle_Mesh_Vars,       ONLY:nTotalSides
USE MOD_Particle_MPI_Vars,        ONLY:PartMPI
#ifdef MPI
USE MOD_Particle_MPI_HALO,        ONLY:WriteParticlePartitionInformation
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
!OUTPUT VARIABLES
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iSide,p,q, nPlanar,nBilinear,nCurved,nDummy
REAL,DIMENSION(1:3)         :: v1,v2
REAL                        :: length,eps
LOGICAL                     :: isLinear
#ifdef MPI  
INTEGER                     :: nPlanarTot,nBilinearTot,nCurvedTot
#endif /*MPI*/
!================================================================================================================================

! allocated here,!!! should be moved?
!ALLOCATE( SideType(nSides)        &
!        , SideDistance(nSides)    &
!        , SideNormVec(1:3,nSides) )

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A)') ' Get Side Type incl. HALO-Sides...'

ALLOCATE( SideType(nTotalSides)        &
        , SideDistance(nTotalSides)    &
        , SideNormVec(1:3,nTotalSides) )

SideDistance=0.
SideNormVec=0.

eps=1e-8
nPlanar=0
nBilinear=0
nCurved=0
#ifdef MPI
nPlanarTot=0
nBilinearTot=0
nCurvedTot=0
#endif /*MPI*/
DO iSide=1,nTotalSides
  isLinear=.TRUE.
  ! all four edges
  !IF(iSide.GT.nSides) IPWRITE(UNIT_stdOut,*) BezierControlPOints3D(:,:,:,iSide)
  IF(SUM(ABS(BezierControlPoints3D(:,:,:,iSide))).LT.1e-10) &
   CALL abort(__STAMP__, &
        ' no BezierControlPoints',PartMPI%MyRank)
  q=0
  v1=BezierControlPoints3D(:,NGeo,q,iSide)-BezierControlPoints3D(:,0,q,iSide)
  DO p=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,iSide)-BezierControlPoints3D(:,0,q,iSide)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! p
  q=NGeo
  v1=BezierControlPoints3D(:,NGeo,q,iSide)-BezierControlPoints3D(:,0,q,iSide)
  DO p=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,iSide)-BezierControlPoints3D(:,0,q,iSide)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! p
  p=0
  v1=BezierControlPoints3D(:,p,NGeo,iSide)-BezierControlPoints3D(:,p,0,iSide)
  DO q=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,iSide)-BezierControlPoints3D(:,p,0,iSide)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! q
  p=NGeo
  v1=BezierControlPoints3D(:,p,NGeo,iSide)-BezierControlPoints3D(:,p,0,iSide)
  DO q=1,NGeo-1
    v2=BezierControlPoints3D(:,p,q,iSide)-BezierControlPoints3D(:,p,0,iSide)
    v2=CROSS(v1,v2)
    length=SQRT(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
    IF(length.GT.eps) isLinear=.FALSE.
  END DO ! q
  IF(isLinear)THEN
    IF(BoundingBoxIsEmpty(iSide))THEN
    !v0= BezierControlPoints3D(:,0,0,iSide)-BezierControlPoints3D(:,NGeo,0,iSide)     &
    !   -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide)  
    !IF(VECNORM(v0).LT.epsilontol)THEN
    !IF(VECNORM(v0).LT.epsilonbilinear)THEN
      SideType(iSide)=PLANAR
      IF(iSide.LE.SideID_Minus_Upper) nPlanar=nPlanar+1
#ifdef MPI
      IF(iSide.GT.SideID_Minus_Upper) nPlanartot=nPlanartot+1
#endif /*MPI*/
      ! compute the norm vec of side and distance from origin
      !v1=BezierControlPoints3D(:,NGeo,0,iSide)-BezierControlPoints3D(:,0,0,iSide)
      !v2=BezierControlPoints3D(:,0,NGeo,iSide)-BezierControlPoints3D(:,0,0,iSide)
      v1=(-BezierControlPoints3D(:,0,0   ,iSide)+BezierControlPoints3D(:,NGeo,0   ,iSide)   &
          -BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )
      
      v2=(-BezierControlPoints3D(:,0,0   ,iSide)-BezierControlPoints3D(:,NGeo,0   ,iSide)   &
          +BezierControlPoints3D(:,0,NGeo,iSide)+BezierControlPoints3D(:,NGeo,NGeo,iSide) )

      SideNormVec(:,iSide) = CROSSNORM(v1,v2)
!      length=SQRT(v2(1)*v2(1)+v1(2)*v1(2)+v1(3)*v1(3))
!      SideNormVec(:,iSide) =SideNormVec(:,iSide)/length
      !      xNodes(:,1)=SuperSampledNodes(1:3,p  ,q  ,SideID)
      !      xNodes(:,2)=SuperSampledNodes(1:3,p+1,q  ,SideID)
      !      xNodes(:,3)=SuperSampledNodes(1:3,p+1,q+1,SideID)
      !      xNodes(:,4)=SuperSampledNodes(1:3,p  ,q+1,SideID)
      v1=0.25*(BezierControlPoints3D(:,0,0,iSide)     &
              +BezierControlPoints3D(:,NGeo,0,iSide)  &
              +BezierControlPoints3D(:,0,NGeo,iSide)  &
              +BezierControlPoints3D(:,NGeo,NGeo,iSide))
      SideDistance(iSide)=DOT_PRODUCT(v1,SideNormVec(:,iSide))
!      IF(iSide.EQ.8)THEN
!        print*,'v1',v1
!        print*,'SideDistance',SideDistance(iSide)
!        read*
!      END IF
    ELSE
      !IPWRITE(UNIT_stdOut,*) 'Boundingboxisempty',boundingboxisempty(iside)
      SideType(iSide)=BILINEAR
      IF(iSide.LE.SideID_Minus_Upper) nBiLinear=nBiLinear+1
#ifdef MPI
      IF(iSide.GT.SideID_Minus_Upper) nBilinearTot=nBilinearTot+1
#endif /*MPI*/
    END IF ! BoundingBoxIsEmpty
  ELSE ! non-linear sides
    SideType(iSide)=CURVED
    IF(iSide.LE.SideID_Minus_Upper) nCurved=nCurved+1
#ifdef MPI
    IF(iSide.GT.SideID_Minus_Upper) nCurvedTot=nCurvedTot+1
#endif /*MPI*/
    IF(BoundingBoxIsEmpty(iSide))THEN
      v1=BezierControlPoints3D(:,NGeo,0,iSide)-BezierControlPoints3D(:,0,0,iSide)
      v2=BezierControlPoints3D(:,0,NGeo,iSide)-BezierControlPoints3D(:,0,0,iSide)
      SideNormVec(:,iSide) = CROSSNORM(v1,v2)
      v1=0.25*(BezierControlPoints3D(:,0,0,iSide)     &
              +BezierControlPoints3D(:,NGeo,0,iSide)  &
              +BezierControlPoints3D(:,0,NGeo,iSide)  &
              +BezierControlPoints3D(:,NGeo,NGeo,iSide))
      SideDistance(iSide)=DOT_PRODUCT(v1,SideNormVec(:,iSide))
    END IF ! BoundingBoxIsEmpty
  END IF ! isLinear
END DO ! iSide

#ifdef MPI
nPlanarTot=nPlanar+nPlanarTot
nBilinearTot=nBilinear+nBilinearTot
nCurvedTot=nCurved+nCurvedTot
IF(MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nPlanar  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nBilinear,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,nCurved  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
ELSE ! no Root
  CALL MPI_REDUCE(nPlanar  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nBilinear,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
  CALL MPI_REDUCE(nCurved  ,nDummy,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,IERROR)
END IF
#endif /*MPI*/

SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar    faces: ', nPlanar
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear faces: ', nBilinear
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of curved    faces: ', nCurved
SWRITE(UNIT_StdOut,'(132("-"))')


#ifdef MPI
CALL  WriteParticlePartitionInformation(nPlanarTot,nBilinearTot,nCurvedTot)
#endif

END SUBROUTINE GetSideType



END MODULE MOD_Particle_Surfaces
