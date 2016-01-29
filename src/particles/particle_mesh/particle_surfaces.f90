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

INTERFACE GetSideSlabNormalsAndIntervals
  MODULE PROCEDURE GetSideSlabNormalsAndIntervals
END INTERFACE

INTERFACE GetElemSlabNormalsAndIntervals
  MODULE PROCEDURE GetElemSlabNormalsAndIntervals
END INTERFACE

INTERFACE EvaluateBezierPolynomialAndGradient
  MODULE PROCEDURE EvaluateBezierPolynomialAndGradient
END INTERFACE
  

PUBLIC::InitParticleSurfaces, FinalizeParticleSurfaces, CalcBiLinearNormVec, &!GetSuperSampledSurface, &
        GetBezierControlPoints3D,CalcBiLinearNormVecBezier,CalcNormVecBezier,GetSideSlabNormalsAndIntervals,&
        GetElemSlabNormalsAndIntervals,EvaluateBezierPolynomialAndGradient

PUBLIC::CalcBiLinearNormAndTang, CalcNormAndTangBezier

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
USE MOD_Mesh_Vars,                  ONLY:nSides,ElemToSide,NGeo,nBCSides,nSides
USE MOD_ReadInTools,                ONLY:GETREAL,GETINT,GETLOGICAL
USE MOD_Particle_Mesh_Vars,         ONLY:PartBCSideList
USE MOD_Particle_Tracking_Vars,     ONLY:DoRefMapping
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,     ONLY:rBoundingBoxChecks,rPerformBezierClip,rPerformBezierNewton
#endif /*CODE_ANALYZE*/
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

BezierNewtonAngle     = GETREAL('BezierNewtonAngle','1.570796326')! 1°=0.01754 (in rad)

epsilontol            = GETREAL('epsilontol','1e-10')
BezierHitEpsBi        = GETREAL('BezierHitEpsBi','1e-12')
BezierHitEpsBi        = 1.0+BezierHitEpsBi
!BezierHitEpsBi=1.0+SQRT(EPSILON(0.0))
!BezierHitEpsBi=1.000800
MinusEps              = -epsilontol
OnePlusEps            = 1.0 + epsilontol
OneMinusEps           = 1.0 - epsilontol
BezierClipTolerance   = GETREAL('BezierClipTolerance','1e-8')
BezierSplitLimit      = GETREAL('BezierSplitLimit','0.6')
BezierSplitLimit      = 2.*BezierSplitLimit
BezierClipMaxIter     = GETINT('BezierClipMaxIter','100')

BezierClipHit         = GETREAL('BezierClipHit','1e-6')
BezierClipHit         = 1.0+BezierClipHit
tmp=2*(NGeo+1)
WRITE(dummy,'(I2.2)') tmp
BezierClipMaxIntersec = GETINT('BezierClipMaxIntersec',dummy)

IF(DoRefMapping)THEN
  !MultipleBCs    = GETLOGICAL('MultibleBCs',".FALSE.")
  ALLOCATE(PartBCSideList(1:nSides))
  PartBCSideList(:) = -1
  DO iSide=1,nBCSides
    PartBCSideList(iSide)=iSide
  END DO 
  !nTotalBCSides=nBCSides
 ! iBCSide=nBCSides
 ! DO iSide=nBCSides+1,nSides
 !   IF(BC(iSide).EQ.1) THEN
 !     iBCSide=iBCSide+1
 !     PartBCSideList(iSide)=iBCSide
 !   END IF
 ! END DO 
 ! nTotalBCSides=iBCSide
END IF

#ifdef CODE_ANALYZE
rBoundingBoxChecks=0.
rPerformBezierClip=0.
rPerformBezierNewton=0.
rTotalBBChecks    =0.
rTotalBezierClips =0.
#endif /*CODE_ANALYZE*/
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

ALLOCATE( locAlpha(1:BezierClipMaxIntersec) &
        , locXi   (1:BezierClipMaxIntersec) &
        , locEta  (1:BezierClipMaxIntersec) )
ALLOCATE( XiArray (1:2,1:BezierClipMaxIter) &
        , EtaArray(1:2,1:BezierClipMaxIter) )


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
SDEALLOCATE(SideSlabNormals)
SDEALLOCATE(SideSlabIntervals)
SDEALLOCATE(ElemSlabNormals)
SDEALLOCATE(ElemSlabIntervals)
SDEALLOCATE(BoundingBoxIsEmpty)
SDEALLOCATE(ElevationMatrix)
SDEALLOCATE(BezierControlPoints3DElevated)
SDEALLOCATE(locAlpha)
SDEALLOCATE(locXi)
SDEALLOCATE(locEta)
SDEALLOCATE(XiArray)
SDEALLOCATE(EtaArray)
SDEALLOCATE(Vdm_Bezier)
SDEALLOCATE(sVdm_Bezier)
SDEALLOCATE(D_Bezier)
SDEALLOCATE(arrayNChooseK)
SDEALLOCATE(BezierControlPoints3DElevated)
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
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierEpsilonBilinear, BiLinearCoeff, SideNormVec, SideDistance
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
!      IF(Displacement.LT.BezierEpsilonBilinear)THEN
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
! book:
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! 1) computation of Bezier Polynomial is from Farin, Chap. 5, p. 57, EQ. 1
!    precompute t^i (1-t)^(n-i) by means of Horner Schema and store values direcly sorted 
! 2) evaluate gradients
!================================================================================================================================
USE MOD_Mesh_Vars,                            ONLY:NGeo
USE MOD_Globals,                              ONLY:CROSSNORM,CROSS
USE MOD_Particle_Surfaces_Vars,               ONLY:BezierControlPoints3D,facNchooseK,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,               ONLY:SideNormVec,D_Bezier
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
REAL,DIMENSION(2,3)                    :: gradXiEta
REAL,DIMENSION(2)                      :: XiOut(2)
!================================================================================================================================

IF(BoundingBoxIsEmpty(SideID))THEN
  CalcNormVecBezier=SideNormVec(1:3,SideID)
ELSE ! no flat side
  ! compute norm vec
  XiOut(1)=Xi
  XiOut(2)=Eta

  CALL EvaluateBezierPolynomialAndGradient(XiOut,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Gradient=gradXiEta)
  CalcNormVecBezier=CROSSNORM(gradXiEta(1,1:3),gradXiEta(2,1:3))
  !nVec=CROSS(u,v)
  !nlength=nVec(1)*nVec(1)+nVec(2)*nVec(2)+nVec(3)*nVec(3)
  !nlength=SQRT(nlength)
  !CalcNormVecBezier=nVec/nlength
END IF ! BoundingBoxIsEmpty

END FUNCTION CalcNormVecBezier


SUBROUTINE EvaluateBezierPolynomialAndGradient(Xi,N_in,iSize,BezierControlPoints,dBezierControlPoints,Point,Gradient) 
!===================================================================================================================================
! evaluate the BezierPolynomial and/or gradient by means of direct evaluation and help of a pseudo-horner scheme
! faster for the evaluation of gradients
! results depend on mode
! mode=1 compute/interpolate point in patch at xi
! mode=2 compute gradient in xi at xi
! mode=3 point and gradient || both
! range: [-1,1]
! info for gradient
! grad(1,:) are all gradients in u direction,...
! book:
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
! 1) computation of Bezier Polynomial is from Farin, Chap. 5, p. 57, EQ. 1
!    precompute t^i (1-t)^(n-i) by means of Horner Schema and store values direcly sorted 
! 2) point value
! 3) evaluate gradients
! optimized in terms of floating point operations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars,               ONLY:facNchooseK,D_Bezier
USE MOD_TimeDisc_Vars,                        ONLY: iter
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
REAL,INTENT(IN)              :: Xi(2)                 ! value of the plane (u,v)
INTEGER,INTENT(IN)           :: N_in
INTEGER,INTENT(IN)           :: iSize                 ! depending if two or three dimensional
REAL,INTENT(IN)              :: BezierControlPoints(1:iSize,0:N_in,0:N_in)
REAL,INTENT(IN),OPTIONAL     :: dBezierControlPoints(1:2,1:iSize,0:N_in,0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(OUT),OPTIONAL    :: Point(1:iSize)        ! point value
REAL,INTENT(OUT),OPTIONAL    :: Gradient(1:2,1:iSize) ! gradient in u,v
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: Mode,iMode
INTEGER                      :: p,q,m
REAL                         :: MinusXi,PlusXI,MinusEta,PlusEta
REAL                         :: xiup(0:N_in),etaup(0:N_in),xidown(0:N_in),etadown(0:N_in)
REAL                         :: xitild(0:N_in), etatild(0:N_in)
!===================================================================================================================================

Mode=0
IF(PRESENT(Point)) Mode=1
iMode=0
IF(Present(Gradient)) iMode=2
Mode=Mode+iMode

M=N_in-1
MinusXi =1.0-xi(1)
PlusXi  =1.0+xi(1)
MinusEta=1.0-xi(2)
PlusEta =1.0+xi(2)

! compute the required stuff || pseudo Horner or precomputation
xiup(0)=1.0
etaup(0)=1.0
! caution, here the indicies are switched from n-j to j  for **down
xidown(N_in)=1.0
etadown(N_in)=1.0
DO p=1,N_in
  xiup   (p)     =xiup   (p-1)*PlusXi
  etaup  (p)     =etaup  (p-1)*PlusEta
  xidown (N_in-p)=xidown (N_in-p+1)*MinusXi
  etadown(N_in-p)=etadown(N_in-p+1)*MinusEta
END DO ! p

! only once
DO p=0,N_in
  xitild (p)=xiup (p)*xidown (p)
  etatild(p)=etaup(p)*etadown(p)
END DO ! p=0,N_in

IF((Mode.EQ.1).OR.(Mode.EQ.3))THEN
  Point=0.
  DO q=0,N_in
    DO p=0,N_in
      ! derivative in xi: d/dxi = u
      !Point(:)=Point(:)+BezierControlPoints(:,p,q)              &
      !                 *facNchooseK(N_In,p)*xiup(p) *xidown(p)  &
      !                 *facNChooseK(N_In,q)*etaup(q)*etadown(q)
      Point(:)=Point(:)+BezierControlPoints(:,p,q)     &
                       *facNchooseK(N_In,p)*xitild (p) &
                       *facNChooseK(N_In,q)*etatild(q)
    END DO ! p
  END DO ! q
  !Point=Point
END IF

IF(Mode.GE.2)THEN ! gradient
  gradient=0.
  IF(PRESENT(dBezierControlPoints))THEN
    DO q=0,N_in
      DO p=0,N_in
        gradient(:,:)=gradient(:,:)+dBezierControlPoints(:,:,p,q)  &
                                   *facNchooseK(N_In,p)*xitild (p) &
                                   *facNChooseK(N_In,q)*etatild(q)
      END DO ! p
    END DO ! q
  ELSE
    DO q=0,N_in
      DO p=0,M
        ! derivative in xi: d/dxi = u
        gradient(1,:)=gradient(1,:)+(BezierControlPoints(:,p+1,q)-BezierControlPoints(:,p,q)) &
                                   *facNchooseK(M,p)*xitild(p)                                &
                                   *facNChooseK(N_in,q)*etatild(q)
        ! derivative in eta: d/deta = v ! caution - exchange indicies
        gradient(2,:)=gradient(2,:)+(BezierControlPoints(:,q,p+1)-BezierControlPoints(:,q,p)) &
                                   *facNchooseK(N_in,q)*xitild(q)                             &
                                   *facNChooseK(M,p)*etatild(p)
      END DO ! p
    END DO ! q
    ! caution, maybe one 1/(b-a) missing??
    gradient=REAL(N_in)*gradient
  END IF
  ! apply gradient
  !dBezierControlPoints=0.
  !DO nn=1,iSize
  !  DO dd=1,2 
  !    DO q=0,N_In
  !      CycPQ(2)=q
  !      DO p=0,N_In
  !        CycPQ(1)=p
  !        ! Matrix-vector multiplication
  !        Cyc=CycPQ
  !        DO l=0,N_in
  !          Cyc(dd)=l  ! d/dxi_dd
  !          dBezierControlPoints(dd,nn,p,q) = dBezierControlPoints(dd,nn,p,q) + &
  !                                            D_Bezier(CycPQ(dd),l)*BezierControlPoints(nn,Cyc(1),Cyc(2))
  !        END DO ! l=0,N_in
  !      END DO ! p=0,N_in
  !    END DO ! p=0,N_in
  !  END DO ! dd=1,iSize
  !END DO ! nn=1,iSize
END IF

END SUBROUTINE EvaluateBezierPolynomialAndGradient


SUBROUTINE GetBezierControlPoints3D(XCL_NGeo,iElem)
!===================================================================================================================================
! computes the nodes for Bezier Control Points for [P][I][C] [A]daptive [S]uper [S]ampled Surfaces [O]perations
! the control points (coeffs for bezier basis) are calculated using the change basis subroutine that interpolates the points 
! from the curved lagrange basis geometry (pre-computed inverse of vandermonde is required)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:ElemToSide,NGeo
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,sVdm_Bezier
USE MOD_Mesh_Vars,                ONLY:nBCSides,nInnerSides,nMPISides_MINE
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
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) ETA_MINUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ETA_MINUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,0,:),BezierControlPoints3D(1:3,:,:,sideID))
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! 6.) ZETA_PLUS
!-----------------------------------------------------------------------------------------------------------------------------------
SideID=ElemToSide(E2S_SIDE_ID,ZETA_PLUS,iElem)
IF(SideID.LE.lastSideID)THEN
  IF(ElemToSide(E2S_FLIP,ZETA_PLUS,iElem).EQ.0) THEN !if flip=0, master side!!
    CALL ChangeBasis2D(3,NGeo,NGeo,sVdm_Bezier,XCL_NGeo(1:3,:,:,NGeo),BezierControlPoints3D(1:3,:,:,sideID))
    CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
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
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(2) ! slave side, SideID=N-p,jSide=q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-p,q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(3) ! slave side, SideID=N-q,jSide=N-p
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,NGeo-q,NGeo-p)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!    CASE(4) ! slave side, SideID=p,jSide=N-q
!      DO q=0,NGeo
!        DO p=0,NGeo
!          BezierControlPoints3D(:,p,q,SideID)=tmp(:,p,NGeo-q)
!        END DO ! p
!      END DO ! q
!      CALL GetSideSlabNormalsAndIntervals(NGeo,SideID)
!  END SELECT
!END IF
END SUBROUTINE GetBezierControlPoints3D

SUBROUTINE GetSideSlabNormalsAndIntervals(NGeo,SideID)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervalls)
! see article:  
!    author = {Shyue-wu Wang and Zen-chung Shih and Ruei-chuan Chang},                    
!    title = {An Efficient and Stable Ray Tracing Algorithm for Parametric Surfaces},     
!    year = {2001},
! original article: oriented-slab box (cartesian bounding box)
!   author = {Yen, Jonathan and Spach, Susan and Smith, Mark and Pulleyblank, Ron},       
!   title = {Parallel Boxing in B-Spline Intersection},                                   
!   issue_date = {January 1991},
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_Globals_Vars,    ONLY:EpsMach
USE MOD_Preproc
USE MOD_Particle_Surfaces_Vars,   ONLY: SideSlabNormals,SideSlabIntervals,BoundingBoxIsEmpty
USE MOD_Particle_Surfaces_Vars,   ONLY: BezierControlPoints3D,BezierControlPoints3DElevated,BezierElevation
USE MOD_Particle_Surfaces_Vars,   ONLY: BezierEpsilonBilinear
#ifdef CODE_ANALYZE
USE MOD_Particle_Surfaces_Vars,   ONLY: SideBoundingBoxVolume
#endif /*CODE_ANALYZE*/
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
LOGICAL            :: SideIsCritical
!===================================================================================================================================


IF(BezierElevation.EQ.0)THEN
  BezierControlPoints3DElevated=BezierControlPoints3D
ELSE
  CALL GetBezierControlPoints3DElevated(NGeo,SideID)
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors
!-----------------------------------------------------------------------------------------------------------------------------------
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
SideSlabNormals(:,1,SideID)=BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideID)                      &
                           -BezierControlPoints3DElevated(:,0,0,SideID)                                         &
                           +BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideID)   &
                           -BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideID)
SideSlabNormals(:,1,SideID)=SideSlabNormals(:,1,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
SideSlabNormals(:,2,SideID)=BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideID)                      &
                           -BezierControlPoints3DElevated(:,0,0,SideID)                                         &
                           +BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideID)   &
                           -BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideID)

!print*,"length(SideSlabNormals(:,1,SideID))",DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID))
!fehlt das?
!SideSlabNormals(:,2,SideID)=SideSlabNormals(:,2,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
SideSlabNormals(:,2,SideID)=CROSSNORM(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID))
!print*,"length(SideSlabNormals(:,2,SideID))",DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID))

!SideSlabNormals(:,2,SideID)=SideSlabNormals(:,2,SideID)/SQRT(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
! n_3=n_1 x n_2
SideSlabNormals(:,3,SideID)=CROSSNORM(SideSlabNormals(:,2,SideID),SideSlabNormals(:,1,SideID))

! check vector length=1
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 1 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,1,SideID)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 2 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,2,SideID)))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,3,SideID),SideSlabNormals(:,3,SideID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 3 does not have the length 1 .',1,DOT_PRODUCT(SideSlabNormals(:,3,SideID),SideSlabNormals(:,3,SideID)))

! check perpendicularity
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 1 and 2 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,2,SideID))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,3,SideID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 1 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,1,SideID),SideSlabNormals(:,3,SideID))))
IF((ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,3,SideID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Side slab normal 2 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(SideSlabNormals(:,2,SideID),SideSlabNormals(:,3,SideID))))
!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervalls beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!SideSlabIntervals(x- x+ y- y+ z- z+, SideID)


! Intervall beta_1
!print*,"SideID",SideID
SideSlabIntervals(:, SideID)=0.

DO q=0,NGeo+BezierElevation
  DO p=0,NGeo+BezierElevation
    IF((p.EQ.0).AND.(q.EQ.0))CYCLE
    skalprod(1)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,1,SideID))
    skalprod(2)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,2,SideID))
    skalprod(3)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                            BezierControlPoints3DElevated(:,0,0,SideID),SideSlabNormals(:,3,SideID))
    IF    (skalprod(1).LT.0.)THEN
      SideSlabIntervals(1, SideID)=MIN(SideSlabIntervals(1,SideID),skalprod(1))
    ELSEIF(skalprod(1).GT.0.)THEN
      SideSlabIntervals(2, SideID)=MAX(SideSlabIntervals(2,SideID),skalprod(1))
    END IF
    IF    (skalprod(2).LT.0.)THEN
      SideSlabIntervals(3, SideID)=MIN(SideSlabIntervals(3,SideID),skalprod(2))
    ELSEIF(skalprod(2).GT.0.)THEN
      SideSlabIntervals(4, SideID)=MAX(SideSlabIntervals(4,SideID),skalprod(2))
    END IF
    IF    (skalprod(3).LT.0.)THEN
      SideSlabIntervals(5, SideID)=MIN(SideSlabIntervals(5,SideID),skalprod(3))
    ELSEIF(skalprod(3).GT.0.)THEN
      SideSlabIntervals(6, SideID)=MAX(SideSlabIntervals(6,SideID),skalprod(3))
    END IF
  END DO !p
END DO !q
!print*,"SideID",SideID," is planar:",SideIsPlanar(SideID)
dx=ABS(ABS(SideSlabIntervals(2, SideID))-ABS(SideSlabIntervals(1, SideID)))
dy=ABS(ABS(SideSlabIntervals(4, SideID))-ABS(SideSlabIntervals(3, SideID)))
dz=ABS(ABS(SideSlabIntervals(6, SideID))-ABS(SideSlabIntervals(5, SideID)))

!-----------------------------------------------------------------------------------------------------------------------------------
! 3.) Is Side critical? (particle path parallel to the larger surface, therefore numerous intersections are possilbe)
! from Wang et al. 2001 An Efficient and Stable ray tracing algorithm for parametric surfaces
! critical sides: h << w <= l
! h: height
! w: width (short side)
! l: length (long side)
!-----------------------------------------------------------------------------------------------------------------------------------

h=dy
w=MIN(dx,dz)
l=MAX(dx,dz)

IF((h.LT.w*1.E-2))THEN!.AND.(w.LE.l))THEN ! critical case possible: second condition is obsolete
  SideIsCritical=.TRUE.
ELSE
  SideIsCritical=.FALSE.
END IF


!-----------------------------------------------------------------------------------------------------------------------------------
! 4.) determine is side is planar -> "BoundingBoxIsEmpty(SideID)"
!     this results also in the decision whether a side is also considered flat or bilinear! 
!-----------------------------------------------------------------------------------------------------------------------------------

dMax=MAX(dx,dy,dz)
dMin=MIN(dx,dy,dz)
IF(dx/dMax.LT.BezierEpsilonBilinear)THEN
  CALL Abort(__STAMP__,&
  'Bezier side length is degenerated. dx/dMax.LT.BezierEpsilonBilinear ->',0,dx/dMax)
END IF
IF(dy/dMax.LT.BezierEpsilonBilinear)THEN
  SideSlabIntervals(3:4, SideID)=0.
  dy=0.
END IF
IF(dz/dMax.LT.BezierEpsilonBilinear)THEN
  CALL Abort(__STAMP__,&
  'Bezier side length is degenerated. dz/dMax.LT.BezierEpsilonBilinear ->',0,dz/dMax)
END IF

IF(dx*dy*dz.LT.0) THEN
  IPWRITE(UNIT_stdOut,*) ' Warning, no bounding box'
  IF(dx*dy*dz.LT.0) CALL Abort(&
  __STAMP__,&
  'A bounding box (for sides) is negative!?. dx*dy*dz.LT.0 ->',0,(dx*dy*dz))
END IF

IF(ALMOSTZERO(dx*dy*dz))THEN ! bounding box volume is approx zeros
  BoundingBoxIsEmpty(SideID)=.TRUE.
ELSE
  BoundingBoxIsEmpty(SideID)=.FALSE.
END IF

#ifdef CODE_ANALYZE
SideBoundingBoxVolume(SideID)=dx*dy*dz
#endif /*CODE_ANALYZE*/
END SUBROUTINE GetSideSlabNormalsAndIntervals


SUBROUTINE GetElemSlabNormalsAndIntervals(NGeo,ElemID)
!===================================================================================================================================
! computes the oriented-slab box for each bezier basis surface (i.e. 3 slab normals + 3 intervalls)
! of each element. This routine must be called after GetSideSlabNormalsAndIntervals(...), because the elevation takes place there
! see article:  
!    author = {Shyue-wu Wang and Zen-chung Shih and Ruei-chuan Chang},                                                              
!    title = {An Efficient and Stable Ray Tracing Algorithm for Parametric Surfaces},                                               
!    year = {2001},
! original article: oriented-slab box (cartesian bounding box)
!   author = {Yen, Jonathan and Spach, Susan and Smith, Mark and Pulleyblank, Ron},                                                 
!   title = {Parallel Boxing in B-Spline Intersection},                                                                             
!   issue_date = {January 1991},
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Mesh_Vars,       ONLY:PartElemToSide,GEO,PartBCSideList
USE MOD_Particle_Surfaces_Vars,   ONLY:ElemSlabNormals,ElemSlabIntervals,BezierControlPoints3DElevated,BezierElevation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID,NGeo
!REAL,INTENT(IN)    :: XCL_NGeo(3,0:NGeo,0:NGeo,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,iLocSide,SideID,SideIDOrigin,tmpSideID
REAL               :: skalprod(3),dx,dy,dz
!===================================================================================================================================

!BezierControlPoints(:,:,:,ElemID)
!ElemSlabNormals( x y z,1 2 3 , ElemID)
IF(GEO%nPeriodicVectors.GT.0)THEN
  CALL  Abort(&
  __STAMP__,&
  ' computation of wrong bounding box!')
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
! 0.) check if side is planar
!-----------------------------------------------------------------------------------------------------------------------------------
!ElemIsPlanar=.FALSE.

!-----------------------------------------------------------------------------------------------------------------------------------
! 1.) slab normal vectors (use the first local element side)
!-----------------------------------------------------------------------------------------------------------------------------------
SideIDOrigin=PartElemToSide(E2S_SIDE_ID,1,ElemID)
! n_1=V_1+V_2 (V: corner vectors in xi-direction)
ElemSlabNormals(:,1,ElemID)= &
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideIDOrigin)- &
                BezierControlPoints3DElevated(:,0,0,SideIDOrigin)+&
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideIDOrigin)-&
                BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideIDOrigin)
ElemSlabNormals(:,1,ElemID)=ElemSlabNormals(:,1,ElemID)/SQRT(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID)))
! n_2=n_1 x (U_1+U_2) (U: corner vectors in eta-direction)
ElemSlabNormals(:,2,ElemID)= &
                BezierControlPoints3DElevated(:,0,NGeo+BezierElevation,SideIDOrigin)   &
                -BezierControlPoints3DElevated(:,0,0,SideIDOrigin)+&
                BezierControlPoints3DElevated(:,NGeo+BezierElevation,NGeo+BezierElevation,SideIDOrigin)&
                -BezierControlPoints3DElevated(:,NGeo+BezierElevation,0,SideIDOrigin)
ElemSlabNormals(:,2,ElemID)=CROSSNORM(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID))
! n_3=n_1 x n_2
ElemSlabNormals(:,3,ElemID)=CROSSNORM(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,1,ElemID))

! check vector length=1
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 1 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,1,ElemID)))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,2,ElemID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 2 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,2,ElemID)))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,3,ElemID),ElemSlabNormals(:,3,ElemID))-1.)).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 3 does not have the length 1 .',1,DOT_PRODUCT(ElemSlabNormals(:,3,ElemID),ElemSlabNormals(:,3,ElemID)))

! check perpendicularity
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 1 and 2 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,2,ElemID))))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,3,ElemID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 1 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,1,ElemID),ElemSlabNormals(:,3,ElemID))))
IF((ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,3,ElemID)))).GT.1.E-6) CALL Abort(&
  __STAMP__,&
  'Element slab normal 2 and 3 are not perpendicular.',0,ABS(DOT_PRODUCT(ElemSlabNormals(:,2,ElemID),ElemSlabNormals(:,3,ElemID))))

!-----------------------------------------------------------------------------------------------------------------------------------
! 2.) slab box intervalls beta_1, beta_2, beta_3
!-----------------------------------------------------------------------------------------------------------------------------------
!ElemSlabIntervals(x- x+ y- y+ z- z+, ElemID)

ElemSlabIntervals(:,ElemID)=0.
DO iLocSide=1,6
  SideID=PartElemToSide(E2S_SIDE_ID,ilocSide,ElemID)
  DO q=0,NGeo+BezierElevation
    DO p=0,NGeo+BezierElevation
      IF((p.EQ.0).AND.(q.EQ.0))CYCLE
      skalprod(1)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,1,ElemID))
      skalprod(2)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,2,ElemID))
      skalprod(3)=DOT_PRODUCT(BezierControlPoints3DElevated(:,p,q,SideID)-&
                              BezierControlPoints3DElevated(:,0,0,SideIDOrigin),ElemSlabNormals(:,3,ElemID))
      IF    (skalprod(1).LT.0.)THEN
        ElemSlabIntervals(1, ElemID)=MIN(ElemSlabIntervals(1,ElemID),skalprod(1))
      ELSEIF(skalprod(1).GT.0.)THEN
        ElemSlabIntervals(2, ElemID)=MAX(ElemSlabIntervals(2,ElemID),skalprod(1))
      END IF
      IF    (skalprod(2).LT.0.)THEN
        ElemSlabIntervals(3, ElemID)=MIN(ElemSlabIntervals(3,ElemID),skalprod(2))
      ELSEIF(skalprod(2).GT.0.)THEN
        ElemSlabIntervals(4, ElemID)=MAX(ElemSlabIntervals(4,ElemID),skalprod(2))
      END IF
      IF    (skalprod(3).LT.0.)THEN
        ElemSlabIntervals(5, ElemID)=MIN(ElemSlabIntervals(5,ElemID),skalprod(3))
      ELSEIF(skalprod(3).GT.0.)THEN
        ElemSlabIntervals(6, ElemID)=MAX(ElemSlabIntervals(6,ElemID),skalprod(3))
      END IF
    END DO !p
  END DO !q
END DO !iLocSide=1:6
dx=ABS(ABS(ElemSlabIntervals(2, ElemID))-ABS(ElemSlabIntervals(1, ElemID)))
dy=ABS(ABS(ElemSlabIntervals(4, ElemID))-ABS(ElemSlabIntervals(3, ElemID)))
dz=ABS(ABS(ElemSlabIntervals(6, ElemID))-ABS(ElemSlabIntervals(5, ElemID)))
IF(dx*dy*dz.LT.0) CALL Abort(&
  __STAMP__,&
  'A bounding box (for elements) is negative!?. dx*dy*dz.LT.0 ->',0,(dx*dy*dz))
!IF((dx*dy*dz).LT.GEO%Volume(ElemID))THEN
!  IPWRITE(*,*) 'Volume', dx*dy*dz
!  IPWRITE(*,*) 'DG-Volume', GEO%Volume(ElemID)
!  CALL Abort(&
!  __STAMP__,&
!  'The bounding box is smaller than element! BoundingBox-volume:')
!END IF
END SUBROUTINE GetElemSlabNormalsAndIntervals


SUBROUTINE GetBezierControlPoints3DElevated(NGeo,SideID)
!===================================================================================================================================
! compute the elevated bezier control points (use pre-computed ElevationMatrix with binomial coefficicents)
! uses the tensor-product basis and combines two 1-D elevation steps
!===================================================================================================================================
! MODULES
!USE MOD_Globals
!USE MOD_Globals_Vars,    ONLY:EpsMach
!USE MOD_Preproc
!USE MOD_Particle_Surfaces_Vars,   ONLY:SideSlabNormals,SideSlabIntervals,BezierControlPoints3D,BoundingBoxIsEmpty
!USE MOD_Particle_Surfaces_Vars,   ONLY:BezierEpsilonBilinear
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierControlPoints3D,BezierControlPoints3DElevated,BezierElevation
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
!REAL               :: skalprod(3),dx,dy,dz,dMax,dMin,w,h,l
!LOGICAL            :: SideIsCritical
REAL               :: temp(1:3,0:NGeo+BezierElevation,0:NGeo)
!===================================================================================================================================
temp=0.
! p-direction
DO q=0,NGeo
  temp(:,:,q)=ElevateBezierPolynomial(NGeo,BezierControlPoints3D(:,:,q,SideID))
END DO
! q-direction
DO p=0,NGeo+BezierElevation
  BezierControlPoints3DElevated(:,p,:,SideID)=ElevateBezierPolynomial(NGeo,temp(:,p,:))
END DO
END SUBROUTINE GetBezierControlPoints3DElevated


FUNCTION ElevateBezierPolynomial(NGeo,BezierPolynomial)
!===================================================================================================================================
! this function creates a new equidistantly distributed set of control points (bézier polynomial basis coefficients) based on the 
! control points "BezierPolynomial" with order "N" and elevates them to "ElevateBezierPolynomial" on "Xi_NGeo_elevated" with order 
! "p+BezierElevation"
! book: single step procedure (1D)
!   author = {Piegl, Les and Tiller, Wayne},                                                                       
!   title = {The NURBS Book (2Nd Ed.)},                                                                            
!   year = {1997},
! book: see also for general remarks
! author = {Farin, Gerald},
! title = {Curves and Surfaces for CAGD: A Practical Guide},
! year = {2002},
!===================================================================================================================================
! MODULES
USE MOD_Particle_Surfaces_Vars,   ONLY:BezierElevation,ElevationMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: NGeo
REAL,INTENT(IN)    :: BezierPolynomial(1:3,0:NGeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL               :: ElevateBezierPolynomial(1:3,0:NGeo+BezierElevation)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL            :: length
INTEGER            :: i,j,jStart,jEnd
!===================================================================================================================================
! algorithm originally from "The NURBS Book" by Les Piegl, Wayne Tiller (p.205)
ElevateBezierPolynomial = 0.
! the first and last points remain the same!
! edge points remain: P_0 and P_p
ElevateBezierPolynomial(:,0)                    = BezierPolynomial(:,0)
ElevateBezierPolynomial(:,NGeo+BezierElevation) = BezierPolynomial(:,NGeo)
! inner points change: P_1,...,P_p-1
DO i=1,NGeo+BezierElevation-1
  jStart = MAX(0,i-BezierElevation)
  jEnd   = MIN(NGeo,i)
  DO j=jStart,jEnd 
    ElevateBezierPolynomial(:,i)=ElevateBezierPolynomial(:,i)+ElevationMatrix(i,j)*BezierPolynomial(:,j)
  END DO
END DO
END FUNCTION ElevateBezierPolynomial


END MODULE MOD_Particle_Surfaces
