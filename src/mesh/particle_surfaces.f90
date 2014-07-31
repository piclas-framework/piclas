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

INTERFACE GetBiLinearPlane
  MODULE PROCEDURE GetBiLinearPlane
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

PUBLIC::GetBiLinearPlane, InitParticleSurfaces, FinalizeParticleSurfaces, CalcBiLinearNormVec

!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleSurfaces()
!===================================================================================================================================
! read required parameters
!===================================================================================================================================
! MODULES
USE MOD_Particle_Surfaces_vars
USE MOD_Mesh_Vars,                  ONLY:nSides
USE MOD_ReadInTools,                ONLY:GetReal
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(ParticleSurfaceInitIsDone) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT PARTICLE SURFACES ...!'

epsilonbilinear = GETREAL('eps-bilinear','1e-6')
epsilontol      = GETREAL('epsOne',1e-12)
epsilonOne      = 1.0 + epsilontol

ALLOCATE( SideIsPlanar(nSides)            &
        , SideDistance(nSides)            &
        , nElemBCSides(PP_nElems)         &
        , BiLinearCoeff(1:3,1:4,1:nSides) )
SideIsPlanar=.FALSE.

CALL GetBiLinearPlane()

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
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( SideIsPlanar, BiLinearCoeff,SideNormVec,nElemBCSides,SideDistance)
ParticleSurfaceInitIsDone=.FALSE.

END SUBROUTINE FinalizeParticleSurfaces

SUBROUTINE GetBilinearPlane()
!===================================================================================================================================
! computes the required coefficients for a bi-linear plane and performs the decision between planar and bi-linear planes
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,                ONLY:nSides,ElemToSide,nBCSides
USE MOD_Particle_Vars,            ONLY:GEO
USE MOD_Particle_Surfaces_Vars,   ONLY:epsilonbilinear, SideIsPlanar,BiLinearCoeff, SideNormVec
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE               :: SideIsDone(:)
REAL                              :: Displacement,nlength
INTEGER                           :: iElem,ilocSide, SideID,iNode,iSide
! debug information
INTEGER                           :: nBilinear,nPlanar
!===================================================================================================================================

ALLOCATE( SideIsDone(1:nSides) )
SideIsDone=.FALSE.
nBiLinear=0
nPlanar=0

DO iElem=1,PP_nElems ! caution, if particles are not seeded in the whole domain
  DO ilocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem) 
    IF(.NOT.SideIsDone(SideID))THEN
      ! compute the bi-linear coefficients for this side
      ! caution: the parameter space is [-1;1] x [-1;1] instead of [0,1]x[0,2] 
      ! the numbering of the nodes should be counterclockwise 
      ! DEBUGGGG!!!!!!!!!!!!!!!!
      ! check if the nodes are in correct numbering
      ! for ray-bi-linear patch intersection see. ramsay
      BiLinearCoeff(:,1,SideID) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))

      BiLinearCoeff(:,2,SideID) =-GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))

      BiLinearCoeff(:,3,SideID) =-GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
                                - GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))

      BiLinearCoeff(:,4,SideID) = GEO%NodeCoords(:,GEO%ElemSideNodeID(1,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(2,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(3,ilocSide,iElem)) &
                                + GEO%NodeCoords(:,GEO%ElemSideNodeID(4,ilocSide,iElem))
      BiLinearCoeff(:,:,SideID) = 0.25*BiLinearCoeff(:,:,SideID)
      ! compute displacement vector (is displacement form planar plane)
      Displacement = BiLinearCoeff(1,1,SideID)*BiLinearCoeff(1,1,SideID) &
                   + BiLinearCoeff(2,1,SideID)*BiLinearCoeff(2,1,SideID) &
                   + BiLinearCoeff(3,1,SideID)*BiLinearCoeff(3,1,SideID) 
      IF(Displacement.LT.epsilonbilinear)THEN
        SideIsPlanar(SideID)=.TRUE.
        nPlanar=nPlanar+1
      ELSE
        nBilinear=nBilinear+1
      END IF
      SideIsDone(SideID)=.TRUE.
    ELSE
      CYCLE  
    END IF ! SideID
  END DO ! ilocSide
END DO ! iElem

! get number of bc-sides of each element
nElemBCSides=0
DO iElem=1,PP_nelems
  DO ilocSide=1,6
    SideID=ElemToSide(E2S_SIDE_ID,ilocSide,ElemID) 
    IF(SideID.LT.nBCSides) nElemBCSides=nElemBCSides+1
  END DO ! ilocSide
END DO ! nElemBCSides

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of planar    surfaces: ', nPlanar
SWRITE(UNIT_StdOut,'(A,I8)') ' Number of bi-linear surfaces: ', nBilinear

ALLOCATE(SideNormVec(nPlanar))
! compute normal vector of planar sides
DO iSide=1,nSides
  IF(SideIsPlanar(SideID))THEN
    SideNormVec(:,SideID)=CROSS(BiLinearCoeff(:,2,SideID),BiLinearCoeff(:,3,SideID))
    nlength=SideNormVec(1,SideID)*SideNormVec(1,SideID) &
           +SideNormVec(2,SideID)*SideNormVec(2,SideID) &
           +SideNormVec(3,SideID)*SideNormVec(3,SideID) 
    SideNormVec(:,SideID) = SideNormVec(:,SideID)/SQRT(nlength)
    SideDisctance(SideID) = DOT_PRODUCT(SideNormVec(:,SideID),BiLinearCoeff(:,4,SideID))
  END IF
END DO

DEALLOCATE( SideIsDone)

END SUBROUTINE GetBilinearPlane

FUNCTION CalcBiLinearNormVec(xi,eta,SideID)
!================================================================================================================================
! function to compute the normal vector of a bi-linear surface
!================================================================================================================================
USE MOD_Globals,                              ONLY:CROSS
USE MOD_Particle_Surfaces_Vars,               ONLY:BiLinearCoeff
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

n=CROSS(a,b)
length=n(1)*n(1)+n(2)*n(2)+n(3)*n(3)
length=SQRT(length)
CalcBiLinearNormVec=n/length

END FUNCTION CalcBiLinearNormVec

END MODULE MOD_Particle_Surfaces
