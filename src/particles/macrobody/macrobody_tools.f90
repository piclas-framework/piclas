!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
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

!===================================================================================================================================
!> Tools for macroscopic bodies inside particle domain
!===================================================================================================================================
MODULE MOD_MacroBody_Tools
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE INSIDEMACROBODY
  MODULE PROCEDURE INSIDEMACROBODY
END INTERFACE

INTERFACE ComputeMacroSphereIntersection
  MODULE PROCEDURE ComputeMacroSphereIntersection
END INTERFACE

INTERFACE GetInteractionWithMacroBody
  MODULE PROCEDURE GetInteractionWithMacroBody
END INTERFACE

PUBLIC :: MarkMacroBodyElems
PUBLIC :: INSIDEMACROBODY
PUBLIC :: ComputeMacroSphereIntersection
PUBLIC :: GetInteractionWithMacroBody
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> 1: check if MacroParticle are inside of Elements and add a safetyfactor to guarantee waterproof tracing
!> 2: calcualte volume portion of each cell, which is occupied with macroparticles
!===================================================================================================================================
SUBROUTINE MarkMacroBodyElems()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: epsMach
USE MOD_Particle_Vars          ,ONLY: ManualTimeStep, nPointsMCVolumeEstimate
USE MOD_MacroBody_Vars         ,ONLY: ElemHasMacroBody, CalcMPVolumePortion
USE MOD_MacroBody_Vars         ,ONLY: MacroSphere, nMacroBody, UseMacroBody
USE MOD_Particle_Mesh_Vars     ,ONLY: nTotalElems, GEO
USE MOD_Mesh_Vars              ,ONLY: XCL_NGeo!, wBaryCL_NGeo, XiCL_NGeo
USE MOD_Mesh_Vars              ,ONLY: NGeo, nElems
USE MOD_TimeDisc_Vars          ,ONLY: dt, RKdtFrac, iter
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem!,TensorProductInterpolation
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iElem,iMB,ii,jj,kk
INTEGER :: kBGM,jBGM,iBGM
REAL    :: MPBounds(1:2,1:3)
INTEGER :: BGMCellXmax,BGMCellXmin
INTEGER :: BGMCellYmax,BGMCellYmin
INTEGER :: BGMCellZmax,BGMCellZmin
REAL    :: DistVec(1:3)
REAL    :: MacroBodyTrajectory(1:3), LengthMacroBodyTrajectory, eps, safetyFac
REAL    :: DistVecLength, BoundsDiagonal, BoundsDiagonalVec(1:3)
INTEGER :: nodesInside, matchedParts
!INTEGER :: nInsertPartsX!,nInsertPartsY,nInsertPartsZ
!INTEGER :: xref,yref,zref
INTEGER :: iPart
REAL    :: refPos(1:3),physPos(1:3)
REAL    :: dtLocal
!===================================================================================================================================
IF (.NOT.UseMacroBody) RETURN

!set a factor for volume around spheres to be marked (necessary for tracking)
safetyFac=1.1

! only if macroparticles are moving or volume portion has to be calculated
! also done in initial iteration
IF (MAXVAL(ABS(MacroSphere(:)%velocity(1))).GT.0. .OR.MAXVAL(ABS(MacroSphere(:)%velocity(2))).GT.0. &
    .OR. MAXVAL(ABS(MacroSphere(:)%velocity(3))).GT.0. .OR. CalcMPVolumePortion) THEN
  CalcMPVolumePortion=.TRUE.
  GEO%MPVolumePortion(:)=0.
  !--- 1: first coarse check if element has MP via BGM with Halo cells included
  DO iMB=1,nMacroBody
    IF (iter.EQ.0) THEN
      IF (ManualTimeStep.EQ.0.0) THEN
         CALL abort(&
              __STAMP__&
              , 'ManualTimeStep.EQ.0.0 -> ManualTimeStep needs to be defined for Macroparticles!'// &
' Particles-ManualTimeStep = ',RealInfoOpt=ManualTimeStep)
      ELSE
        dtLocal=ManualTimeStep
      END IF
    ELSE
      dtLocal=dt*RKdtfrac
    END IF
    MacroBodyTrajectory(1:3)=MacroSphere(iMB)%velocity(1:3)*dtLocal
    LengthMacroBodyTrajectory=SQRT(DOT_PRODUCT(MacroBodyTrajectory,MacroBodyTrajectory))
    IF (LengthMacroBodyTrajectory.GT.0) MacroBodyTrajectory=MacroBodyTrajectory/LengthMacroBodyTrajectory

    ! define a box with BGM around sphere center
    MPBounds(1,1:3)=MacroSphere(iMB)%center(1:3)-(MacroSphere(iMB)%radius*safetyFac+LengthMacroBodyTrajectory+epsMach)
    MPBounds(2,1:3)=MacroSphere(iMB)%center(1:3)+(MacroSphere(iMB)%radius*safetyFac+LengthMacroBodyTrajectory+epsMach)
    BGMCellXmin = MAX(GEO%TFIBGMimin,CEILING((MPBounds(1,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)))
    BGMCellXmax = MIN(GEO%TFIBGMimax,CEILING((MPBounds(2,1)-GEO%xminglob)/GEO%FIBGMdeltas(1)))
    BGMCellYmin = MAX(GEO%TFIBGMjmin,CEILING((MPBounds(1,2)-GEO%xminglob)/GEO%FIBGMdeltas(2)))
    BGMCellYmax = MIN(GEO%TFIBGMjmax,CEILING((MPBounds(2,2)-GEO%xminglob)/GEO%FIBGMdeltas(2)))
    BGMCellZmin = MAX(GEO%TFIBGMkmin,CEILING((MPBounds(1,3)-GEO%xminglob)/GEO%FIBGMdeltas(3)))
    BGMCellZmax = MIN(GEO%TFIBGMkmax,CEILING((MPBounds(2,3)-GEO%xminglob)/GEO%FIBGMdeltas(3)))
    !IF (BGMCellXmin.LT.GEO%TFIBGMimin) BGMCellXmin=GEO%TFIBGMimin
    !BGMCellXmax = CEILING((MPBounds(2,1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
    !IF (BGMCellXmax.GT.GEO%TFIBGMimax) BGMCellXmax=GEO%TFIBGMimax
    !BGMCellYmin = CEILING((MPBounds(1,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
    !IF (BGMCellYmin.LT.GEO%TFIBGMjmin) BGMCellYmin=GEO%TFIBGMjmin
    !BGMCellYmax = CEILING((MPBounds(2,2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
    !IF (BGMCellYmax.GT.GEO%TFIBGMjmax) BGMCellYmax=GEO%TFIBGMjmax
    !BGMCellZmin = CEILING((MPBounds(1,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
    !IF (BGMCellZmin.LT.GEO%TFIBGMkmin) BGMCellZmin=GEO%TFIBGMkmin
    !BGMCellZmax = CEILING((MPBounds(2,3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
    !IF (BGMCellZmax.GT.GEO%TFIBGMkmax) BGMCellZmax=GEO%TFIBGMkmax

    ! add current Element to BGM-Elem
    DO kBGM = BGMCellZmin,BGMCellZmax
      DO jBGM = BGMCellYmin,BGMCellYmax
        DO iBGM = BGMCellXmin,BGMCellXmax
          DO iElem = 1,GEO%TFIBGM(iBGM,jBGM,kBGM)%nElem
            ElemHasMacroBody(GEO%TFIBGM(iBGM,jBGM,kBGM)%Element(iElem),iMB) = .TRUE.
          END DO
        END DO ! kBGM
      END DO ! jBGM
    END DO ! iBGM

    ! loop over all elements that have a MP in BGM and find those which are in a certain safety radius
    DO iElem=1,nTotalElems
      IF (ElemHasMacroBody(iElem,iMB)) THEN
        ! check with all 3 element diagonals wether macroparticle is smaller than element
        BoundsDiagonalVec(1:3)=GEO%BoundsOfElem(2,1:3,iElem)-GEO%BoundsOfElem(1,1:3,iElem)
        ! choose the greater length
        BoundsDiagonal=MAX( SQRT(DOT_PRODUCT(BoundsDiagonalVec,BoundsDiagonalVec)) , MacroSphere(iMB)%radius)*safetyFac
        ElemHasMacroBody(iElem,iMB)=.FALSE.
        nodesInside=0
        ! check distance of each node if it is within the defined distance (Distvec.LE.Boundsdiagonal) of te sphere center
        DO kk = 0,NGeo
          IF (.NOT.ElemHasMacroBody(iElem,iMB) .OR. CalcMPVolumePortion) THEN
            DO ii=0,NGeo
              IF (.NOT.ElemHasMacroBody(iElem,iMB) .OR. CalcMPVolumePortion) THEN
                DO jj=0,NGeo
                  IF (.NOT.ElemHasMacroBody(iElem,iMB) .OR. CalcMPVolumePortion) THEN
                    DistVec(1:3)=XCL_NGeo(1:3,ii,jj,kk,iElem)-MacroSphere(iMB)%center(1:3)
                    DistVecLength=SQRT(DOT_PRODUCT(DistVec,DistVec))
                    IF (DistVecLength.GT.0) DistVec=DistVec/DistVecLength
                    eps = LengthMacroBodyTrajectory*MAX(0.,DOT_PRODUCT(DistVec,MacroBodyTrajectory))+epsMach
                    IF (DistVecLength.LE.BoundsDiagonal+eps) THEN
                      ElemHasMacroBody(iElem,iMB)=.TRUE.
                    END IF
                    IF (DistVecLength.LE.MacroSphere(iMB)%radius) THEN
                      nodesInside=nodesInside+1
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
        IF (CalcMPVolumePortion .AND. nodesInside.EQ.(NGeo+1)**3 .AND. iElem.LE.nElems) GEO%MPVolumePortion(iElem)=1.0
      END IF
    END DO
  END DO

!--- 2: calculate volume portions using monte carlo inserting and rejections
  IF (CalcMPVolumePortion) THEN
    !nInsertPartsX=INT(nPointsMCVolumeEstimate**(1./3.))
    !nInsertPartsY=nInsertPartsX
    !nInsertPartsZ=nInsertPartsX
    DO iElem=1,nElems
      IF (ANY(ElemHasMacroBody(iElem,:)) .AND. GEO%MPVolumePortion(iElem).LT.1.0) THEN
        matchedParts=0
        !! if element is a cube equidistant inserting in reference space
        !DO xref=1,nInsertPartsX
        !  DO yref=1,nInsertPartsY
        !    DO zref=1,nInsertPartsZ
        !      ! refpos needs to be between -1:+1 not 0:1
        !      refPos(1) = (REAL(xref)-0.5) * 1./REAL(nInsertPartsX)
        !      refPos(2) = (REAL(yref)-0.5) * 1./REAL(nInsertPartsY)
        !      refPos(3) = (REAL(zref)-0.5) * 1./REAL(nInsertPartsZ)
        !      CALL TensorProductInterpolation(refPos(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem)&
        !                ,physPos(1:3)) !Map into phys. space
        !      IF (INSIDEMACROBODY(physPos)) matchedParts=matchedParts+1
        !    END DO
        !  END DO
        !END DO
        ! for arbitrary elements, positions are chosen randomly and checked if they are in element first
        DO iPart=1,nPointsMCVolumeEstimate
          DO
            CALL RANDOM_NUMBER(physPos)
            PhysPos = GEO%BoundsOfElem(1,:,iElem) + physPos*(GEO%BoundsOfElem(2,:,iElem)-GEO%BoundsOfElem(1,:,iElem))
            CALL GetPositionInRefElem(physPos,refPos,iElem)
            IF(ALL(ABS(refPos).LE.1.0)) EXIT ! particle inside of element
          END DO
          IF (INSIDEMACROBODY(physPos)) matchedParts=matchedParts+1
        END DO
        GEO%MPVolumePortion(iElem)=REAL(matchedParts)/REAL(nPointsMCVolumeEstimate)
      END IF
    END DO
  END IF ! CalcMPVolumePortion
END IF ! Velo_Macropart > 0 | CalcMPVolumePortion

CalcMPVolumePortion=.FALSE.

END SUBROUTINE MarkMacroBodyElems


!===================================================================================================================================
!> Function for checking if particle position would be inside of any macro-particle in the local domain
!===================================================================================================================================
LOGICAL FUNCTION INSIDEMACROBODY(Particle_pos)
! MODULES
USE MOD_MacroBody_Vars ,ONLY: MacroSphere, nMacroBody
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: refPosSphere(1:3), distance
INTEGER :: iMB
!===================================================================================================================================
INSIDEMACROBODY = .FALSE.
DO iMB=1,nMacroBody
  refPosSphere(1:3) = MacroSphere(iMB)%center(1:3)
  distance=SQRT(DOT_PRODUCT((Particle_pos-refPosSphere),(Particle_pos-refPosSphere)))
  IF (distance.LE.MacroSphere(iMB)%radius) THEN
    INSIDEMACROBODY=.TRUE.
    RETURN
  END IF
END DO
END FUNCTION INSIDEMACROBODY


!===================================================================================================================================
!> Calculates intersection of particle path with defined spherical, solid, moving macroparticle
!===================================================================================================================================
SUBROUTINE ComputeMacroSphereIntersection(isHit,PartTrajectory,lengthPartTrajectory,macroBodyID &
                                         ,alpha,alphaSphere,alphaDoneRel,partID,alpha2)
! MODULES
USE MOD_Globals
USE MOD_Utils                  ,ONLY: QuadraticSolver
USE MOD_Particle_Vars          ,ONLY: LastPartPos
USE MOD_MacroBody_Vars         ,ONLY: MacroSphere
USE MOD_TimeDisc_Vars          ,ONLY: dt,RKdtFrac
#ifdef CODE_ANALYZE
USE MOD_Particle_Vars          ,ONLY: PartState
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
REAL,INTENT(IN),DIMENSION(1:3)    :: PartTrajectory
REAL,INTENT(IN)                   :: lengthPartTrajectory
INTEGER,INTENT(IN)                :: partID
INTEGER,INTENT(IN)                :: macroBodyID
REAL,INTENT(IN)                   :: alphaDoneRel
REAL,INTENT(IN),OPTIONAL          :: alpha2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                  :: alpha
REAL,INTENT(OUT)                  :: alphaSphere
LOGICAL,INTENT(OUT)               :: isHit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: A,B,C
REAL                              :: t(2), scaleFac
INTEGER                           :: InterType,nRoot
REAL                              :: PosSphere(1:3),P2_rel(1:3),relPartTrajectory(1:3)
REAL                              :: relLengthPartTrajectory
REAL                              :: MacroBodyTrajectory(1:3)
#ifdef CODE_ANALYZE
REAL                              :: distance, distance2
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
! set alpha to minus one // no intersection
alpha=-1.0
isHit=.FALSE.
#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(64("-"))')
      WRITE(UNIT_stdout,'((A,G0))')  '     | partID: ',PartID
      WRITE(UNIT_stdout,'((A,G0))')  '     | macroBodyID: ',macroBodyID
      WRITE(UNIT_stdout,'((A,G0))')  '     | alphaDoneRel: ',alphaDoneRel
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! transform particle trajectory to sphere system v2_rel=v2-v1
MacroBodyTrajectory(1:3)=MacroSphere(macroBodyID)%velocity(1:3)*dt*RKdtFrac*(1.-alphaDoneRel)
relPartTrajectory(1:3)=PartTrajectory(1:3)*LengthPartTrajectory-MacroBodyTrajectory(1:3)
! Transform particle position to sphere system
PosSphere(1:3) = MacroSphere(macroBodyID)%center(1:3)+MacroSphere(macroBodyID)%velocity(1:3)*dt*RKdtFrac*alphaDoneRel
P2_rel(1:3)=LastPartPos(PartID,1:3)-PosSphere(1:3)

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      CALL OutputTrajectory(PartID,PartState(PartID,1:3),PartTrajectory,lengthPartTrajectory)
      WRITE(UNIT_stdOut,'(A,3(X,E15.8))') '     | Macrosphere center Pos:       ', PosSphere(1:3)
      WRITE(UNIT_stdOut,'(A,3(X,E15.8))') '     | Macrosphere Trajectory:       ', MacroBodyTrajectory(1:3)
      WRITE(UNIT_stdOut,'(A,3(X,E15.8))') '     | relative Position:       ', P2_rel(1:3)
      WRITE(UNIT_stdOut,'(A,3(X,E15.8))') '     | relative Trajectory:       ', relPartTrajectory(1:3)
      WRITE(UNIT_stdOut,'(A,L)') '     | particle moves away from sphere:  ',(DOT_PRODUCT(P2_rel,relPartTrajectory).GT.0)
      WRITE(UNIT_stdout,'(12("-"))')
    END IF
  END IF
#endif /*CODE_ANALYZE*/

! particle moves away from sphere in reference system
IF (DOT_PRODUCT(P2_rel,relPartTrajectory).GT.0) RETURN
! calculate lenght of v2_rel
relLengthPartTrajectory=SQRT(DOT_PRODUCT(relPartTrajectory,relPartTrajectory))
IF (relLengthPartTrajectory.EQ.0) RETURN

A = DOT_PRODUCT(relPartTrajectory,relPartTrajectory)
B = 2*DOT_PRODUCT(P2_rel,relPartTrajectory)
C = DOT_PRODUCT(P2_rel,P2_rel) - MacroSphere(macroBodyID)%radius**2

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | Quadratic equation constants in MacroParticle intersection: '
      WRITE(UNIT_stdout,'(3(A,G0))') '     | A: ',A,' | B: ',B,' | C: ',C
    END IF
  END IF
#endif /*CODE_ANALYZE*/

scaleFac = relLengthPartTrajectory**2 ! * MacroSphere(macroBodyID)%radius !<...>^2 * cell-scale
scaleFac = 1./scaleFac
A = A * scaleFac
B = B * scaleFac
C = C * scaleFac

#ifdef CODE_ANALYZE
  IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
    IF(PartID.EQ.PARTOUT)THEN
      WRITE(UNIT_stdout,'(A)') '     | Quadratic equation constants in MacroParticle intersection (after scaling): '
      WRITE(UNIT_stdout,'(3(A,G0))') '     | A: ',A,' | B: ',B,' | C: ',C
    END IF
  END IF
#endif /*CODE_ANALYZE*/

CALL QuadraticSolver(A,B,C,nRoot,t(1),t(2))

#ifdef CODE_ANALYZE
distance = SQRT(DOT_PRODUCT(P2_rel,P2_rel))
distance2 = SQRT(DOT_PRODUCT(P2_rel+relPartTrajectory,P2_rel+relParttrajectory))
IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
  IF(PartID.EQ.PARTOUT)THEN
    WRITE(UNIT_stdout,'(2(A,G0))') '     | distance last partpos: ',distance,' | distance partstate: ',distance2
    WRITE(UNIT_stdout,'(2(A,G0))') '     | rel_alpha1 [t(1)]: ',t(1),' | rel_alpha2 [t(2)]: ',t(2)
    WRITE(UNIT_stdout,'(2(A,G0))') '     | alpha(1): ',t(1)*LengthPartTrajectory,' | alpha(2): ',t(2)*LengthPartTrajectory
    WRITE(UNIT_stdout,'(64("-"))')
  END IF
END IF
IF (distance.LE.MacroSphere(macroBodyID)%radius*0.99 .AND. distance2.LE.MacroSphere(macroBodyID)%radius*0.99) THEN
  WRITE(UNIT_stdout,'(A)') '     | Particle is inside of Macro particle (sphere): '
  WRITE(UNIT_stdout,'((A,G0))')  '     | partID: ',PartID
  WRITE(UNIT_stdout,'((A,G0))')  '     | macroBodyID: ',macroBodyID
  WRITE(UNIT_stdout,'(2(A,G0))') '     | distance last partpos: ',distance,' | distance partstate: ',distance2
  CALL abort(&
      __STAMP__&
      ,'Particle is inside of Macro particle (sphere)')
END IF
#endif /*CODE_ANALYZE*/

IF(nRoot.EQ.0)THEN
  RETURN
END IF

IF (nRoot.EQ.1) THEN
  IF(t(1).LE.1.0 .AND. t(1).GE.0.) THEN
    alpha=t(1)*LengthPartTrajectory
    alphaSphere=t(1)
    IF (PRESENT(alpha2)) THEN
      IF (alpha2.GT.-1.0) THEN
        IF (ALMOSTEQUAL(alpha,alpha2)) THEN
          alpha=-1.0
          RETURN
        END IF
      END IF
    END IF
    isHit=.TRUE.
    RETURN
  ELSE
    RETURN
  END IF
ELSE
  InterType=0

  IF(t(1).LE.1.0 .AND. t(1).GE.0.) THEN
    InterType=InterType+1
    IF (PRESENT(alpha2)) THEN
      IF (alpha2.GT.-1.0) THEN
        IF (ALMOSTEQUAL(t(1)*LengthPartTrajectory,alpha2)) THEN
          t(1)=-1.0
          InterType=InterType-1
        END IF
      END IF
    END IF
  ELSE
    t(1)=-1.
  END IF

  IF(t(2).LE.1.0 .AND. t(2).GE.0.) THEN
    InterType=InterType+2
    IF (PRESENT(alpha2)) THEN
      IF (alpha2.GT.-1.0) THEN
        IF (ALMOSTEQUAL(t(2)*LengthPartTrajectory,alpha2)) THEN
          t(2)=-1.0
          InterType=InterType-2
        END IF
      END IF
    END IF
  ELSE
    t(2)=-1.
  END IF

  IF(InterType.EQ.0) THEN
    RETURN
  END IF
  isHit=.TRUE.
  SELECT CASE(InterType)
  CASE(1)
    alpha=t(1)*LengthPartTrajectory
    alphaSphere=t(1)
  CASE(2)
    alpha=t(2)*LengthPartTrajectory
    alphaSphere=t(2)
  CASE DEFAULT
   ! two intersections
    IF(t(1).LT.t(2))THEN
      alpha=t(1)*LengthPartTrajectory
      alphaSphere=t(1)
    ELSE
      alpha=t(2)*LengthPartTrajectory
      alphaSphere=t(2)
    END IF
  END SELECT
  RETURN
END IF

END SUBROUTINE ComputeMacroSphereIntersection


!===================================================================================================================================
!> Computes the post boundary state of a particle that interacts with an spherical macro body
!===================================================================================================================================
SUBROUTINE GetInteractionWithMacroBody(PartTrajectory,lengthPartTrajectory &
                                      ,alpha,alphaSphere,alphaDoneRel,macroBodyID,partID,opt_Reflected)
! MODULES
USE MOD_PreProc
USE MOD_Globals                 ,ONLY: CROSSNORM,abort,UNITVECTOR,CROSS
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars           ,ONLY: PartSpecies
USE MOD_Particle_Vars           ,ONLY: PartState,LastPartPos,Species
USE MOD_MacroBody_Vars          ,ONLY: MacroSphere
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Boundary_Tools ,ONLY: SurfaceToPartEnergyInternal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: PartID,macroBodyID
REAL   ,INTENT(IN)                   :: alpha
REAL   ,INTENT(IN)                   :: alphaSphere
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                   :: alphaDoneRel
REAL,INTENT(INOUT)                   :: PartTrajectory(1:3),lengthPartTrajectory
LOGICAL,INTENT(OUT),OPTIONAL         :: opt_Reflected
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                 :: intersectPoint(1:3),nLoc(1:3),relVeloPart(1:3),sphereIntersectionPos(1:3)
REAL                                 :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
REAL                                 :: EtraWall, EtraNew
REAL                                 :: WallVelo(1:3), WallTemp, TransACC
REAL                                 :: tang1(1:3),tang2(1:3), NewVelo(3)
REAL                                 :: POI_fak!,TildTrajectory(3)
REAL                                 :: Phi, Cmr, VeloCx, VeloCy, VeloCz
REAL                                 :: relPartTrajectory(1:3)!,relLengthPartTrajectory
REAL                                 :: PosSphere(1:3)!,P2_rel(1:3)
REAL                                 :: MacroBodyTrajectory(1:3)
REAL                                 :: force(1:3), moment(1:3), inertiaMoment
!===================================================================================================================================
! transform particle trajectory to sphere system v2_rel=v2-v1
MacroBodyTrajectory(1:3)=MacroSphere(macroBodyID)%velocity(1:3)*dt*RKdtFrac*(1.-alphaDoneRel)
relPartTrajectory(1:3)=PartTrajectory(1:3)*lengthPartTrajectory-MacroBodyTrajectory(1:3)
! Transform particle posiiton to sphere system
PosSphere(1:3) = MacroSphere(macroBodyID)%center(1:3)+MacroSphere(macroBodyID)%velocity(1:3)*dt*RKdtFrac*alphaDoneRel
!P2_rel(1:3)=LastPartPos(PartID,1:3)-PosSphere(1:3)

! already checked during intersection
!IF (DOT_PRODUCT(P2_rel,relPartTrajectory).GT.0) THEN
!  ! particle moves away from sphere in reference system
!  IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
!  RETURN
!END IF
! calculate lenght of v2_rel (already checked in intersection)
!relLengthPartTrajectory=SQRT(DOT_PRODUCT(relPartTrajectory,relPartTrajectory))
!IF (relLengthPartTrajectory.EQ.0) THEN
!  ! particle and sphere move in the same direction with same velocities
!  IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
!  RETURN
!END IF
!relPartTrajectory=relPartTrajectory/relLengthPartTrajectory

intersectPoint(1:3) = LastPartPos(PartID,1:3) + alpha*PartTrajectory(1:3)
sphereIntersectionPos(1:3) = PosSphere(1:3) + alphaSphere*MacroBodyTrajectory(1:3)
nLoc = UNITVECTOR(intersectPoint - sphereIntersectionPos)
! nLoc points outwards of sphere
IF(DOT_PRODUCT(nLoc,relPartTrajectory).GT.0.)  THEN
  IF(PRESENT(opt_Reflected)) opt_Reflected=.FALSE.
  RETURN
ELSE IF(DOT_PRODUCT(nLoc,relPartTrajectory).LE.0.) THEN
  IF(PRESENT(opt_Reflected)) opt_Reflected=.TRUE.
END IF
! change nLoc to point inwards of sphere
nLoc=-nLoc

! transform velocity in reference to sphere
relVeloPart(1:3)=PartState(PartID,4:6)-MacroSphere(macroBodyID)%velocity(1:3)
! perfect reflection on sphere
CALL RANDOM_NUMBER(RanNum)
IF (RanNum.GE.MacroSphere(macroBodyID)%momentumACC) THEN
  NewVelo(1:3) =  relVeloPart(1:3) - 2.*DOT_PRODUCT(relVeloPart(1:3),nLoc)*nLoc
ELSE
  WallVelo = CROSS(MacroSphere(macroBodyID)%velocity(4:6),-nLoc*MacroSphere(macroBodyID)%radius)
  relVeloPart(1:3)=relVelopart(1:3) - WallVelo
  IF (nLoc(3).NE.0.) THEN
    tang1(1) = 1.0
    tang1(2) = 1.0
    tang1(3) = -(nLoc(1)+nLoc(2))/nLoc(3)
  ELSE
    IF (nLoc(2).NE.0.) THEN
      tang1(1) = 1.0
      tang1(3) = 1.0
      tang1(2) = -(nLoc(1)+nLoc(3))/nLoc(2)
    ELSE
      IF (nLoc(1).NE.0.) THEN
        tang1(2) = 1.0
        tang1(3) = 1.0
        tang1(1) = -(nLoc(2)+nLoc(3))/nLoc(1)
      ELSE
        CALL abort(&
            __STAMP__&
            ,'Error in GetInteractionWithMacroPart, n_vec is zero for macro particle',macroBodyID)
      END IF
    END IF
  END IF
  tang1=UNITVECTOR(tang1)
  tang2=CROSSNORM(nLoc,tang1)

  ! diffuse reflection
  WallTemp   = MacroSphere(macroBodyID)%temp
  TransAcc   = MacroSphere(macroBodyID)%transACC
  CALL RANDOM_NUMBER(RanNum)
  VeloCrad    = SQRT(-LOG(RanNum))
  CALL RANDOM_NUMBER(RanNum)
  VeloCz      = SQRT(-LOG(RanNum))
  Fak_D       = VeloCrad**2 + VeloCz**2
  EtraWall    = BoltzmannConst * WallTemp * Fak_D
  VeloReal    = SQRT(DOT_PRODUCT(relVeloPart,relVeloPart))
  EtraOld     = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
  EtraNew     = EtraOld + TransACC * (EtraWall - EtraOld)
  Cmr         = SQRT(2.0 * EtraNew / (Species(PartSpecies(PartID))%MassIC * Fak_D))
  CALL RANDOM_NUMBER(RanNum)
  Phi     = 2.0 * PI * RanNum
  VeloCx  = Cmr * VeloCrad * COS(Phi) ! tang1
  VeloCy  = Cmr * VeloCrad * SIN(Phi) ! tang2
  VeloCz  = Cmr * VeloCz
  NewVelo(1:3) = VeloCx*tang1-tang2*VeloCy-VeloCz*nLoc + WallVelo

  ! Adding the energy that is transferred from the surface onto the internal energies of the particle
  CALL SurfaceToPartEnergyInternal(PartID,WallTemp)
END IF
! transform velocity back to mesh reference
PartState(PartID,4:6) = NewVelo(1:3) + MacroSphere(macroBodyID)%velocity(1:3)

! intersection point with surface
LastPartPos(PartID,1:3) = intersectPoint(1:3)
POI_fak = alphaDoneRel+(1.-alphaDoneRel)*alphaSphere
PartState(PartID,1:3)   = LastPartPos(PartID,1:3) + (1.0 - POI_fak) * dt*RKdtFrac * (PartState(PartID,4:6))

PartTrajectory=PartState(PartID,1:3) - LastPartPos(PartID,1:3)
lengthPartTrajectory=SQRT(DOT_PRODUCT(PartTrajectory,PartTrajectory))
IF (lengthPartTrajectory.GT.0.) PartTrajectory=PartTrajectory/lengthPartTrajectory

!----  Sampling Forces at MacroPart and calculating velocity change of macroparticle due to impule change
force(1:3) = Species(PartSpecies(PartID))%MassIC &
           * (relVeloPart(1:3) - NewVelo(1:3)) * Species(PartSpecies(PartID))%MacroParticleFactor / dt
! delta velocity
MacroSphere(macroBodyID)%RHS(1:3) = MacroSphere(macroBodyID)%RHS(1:3) + force(1:3)*dt/MacroSphere(macroBodyID)%mass

! moment and moment of inertia
moment(1:3) = CROSS(-nLoc*MacroSphere(macroBodyID)%radius,Force(1:3))
inertiaMoment = 2./5.*MacroSphere(macroBodyID)%mass*MacroSphere(macroBodyID)%radius**2
! delta rot velo
MacroSphere(macroBodyID)%RHS(4:6) = MacroSphere(macroBodyID)%RHS(4:6) + moment(1:3)*dt/inertiaMoment

END SUBROUTINE GetInteractionWithMacroBody


END MODULE MOD_MacroBody_Tools
