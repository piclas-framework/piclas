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

MODULE MOD_Particle_Localization
!===================================================================================================================================
! Contains global variables provided by the particle surfaces routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
INTERFACE SinglePointToElement
  MODULE PROCEDURE SinglePointToElement
END INTERFACE

INTERFACE LocateParticleInElement
  MODULE PROCEDURE LocateParticleInElement
END INTERFACE

PUBLIC :: SinglePointToElement,LocateParticleInElement

CONTAINS


SUBROUTINE LocateParticleInElement(PartID,doHALO)
!----------------------------------------------------------------------------------------------------------------------------------!
! Finds a single particle in its host element
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: PDM,PEM,PartState,PartPosRef
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: PartID
LOGICAL,INTENT(IN) :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ElemID
!===================================================================================================================================
ElemID = SinglePointToElement(PartState(PartID,1:3),doHALO=doHALO)
PEM%Element(PartID) = ElemID
IF(ElemID.EQ.-1)THEN
  PDM%ParticleInside(PartID)=.FALSE.
ELSE
  PDM%ParticleInside(PartID)=.TRUE.
  IF(DoRefMapping)THEN
    CALL GetPositionInRefElem(PartState(PartID,1:3),PartPosRef(1:3,PartID),ElemID)
  END IF ! DoRefMapping
END IF ! ElemID.EQ.-1
END SUBROUTINE LocateParticleInElement


!===================================================================================================================================
!> This function maps a 3D point to an element
!> returns elementID or -1 in case no element was found
!===================================================================================================================================
INTEGER FUNCTION SinglePointToElement(Pos3D,doHALO)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemRadius2NGeo
USE MOD_Particle_Mesh_Vars     ,ONLY: Geo
USE MOD_Utils                  ,ONLY: InsertionSort
USE MOD_Particle_Tracking_Vars ,ONLY: Distance, ListDistance, TriaTracking
USE MOD_Mesh_Vars              ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D,GetCNElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: FIBGM_nElems, FIBGM_offsetElem, FIBGM_Element
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)    :: Pos3D(1:3)
LOGICAL,INTENT(IN) :: doHalo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iBGMElem,nBGMElems, ElemID, iBGM,jBGM,kBGM
REAL    :: Distance2, RefPos(1:3)
REAL    :: Det(6,2)
LOGICAL :: InElementCheck
!===================================================================================================================================
SinglePointToElement = -1

!IF ( (Pos3D(1).LT.GEO%xminglob).OR.(Pos3D(1).GT.GEO%xmaxglob).OR. &
!     (Pos3D(2).LT.GEO%yminglob).OR.(Pos3D(2).GT.GEO%ymaxglob).OR. &
!     (Pos3D(3).LT.GEO%zminglob).OR.(Pos3D(3).GT.GEO%zmaxglob)) RETURN

! --- get background mesh cell of point
iBGM = CEILING((Pos3D(1)-GEO%xminglob)/GEO%FIBGMdeltas(1))
iBGM = MAX(MIN(GEO%FIBGMimax,iBGM),GEO%FIBGMimin)
jBGM = CEILING((Pos3D(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
jBGM = MAX(MIN(GEO%FIBGMjmax,jBGM),GEO%FIBGMjmin)
kBGM = CEILING((Pos3D(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
kBGM = MAX(MIN(GEO%FIBGMkmax,kBGM),GEO%FIBGMkmin)

!--- check all cells associated with this beckground mesh cell
nBGMElems = FIBGM_nElems(iBGM,jBGM,kBGM)

! get closest element barycenter
Distance=-1.

ListDistance=0
DO iBGMElem = 1, nBGMElems
  ElemID = GetCNElemID(FIBGM_Element(FIBGM_offsetElem(iBGM,jBGM,kBGM)+iBGMElem))

  Distance2=(Pos3D(1)-ElemBaryNGeo(1,ElemID))*(Pos3D(1)-ElemBaryNGeo(1,ElemID)) &
           +(Pos3D(2)-ElemBaryNGeo(2,ElemID))*(Pos3D(2)-ElemBaryNGeo(2,ElemID)) &
           +(Pos3D(3)-ElemBaryNGeo(3,ElemID))*(Pos3D(3)-ElemBaryNGeo(3,ElemID))

  IF(Distance2.GT.ElemRadius2NGeo(ElemID))THEN
    Distance(iBGMElem)=-1.
  ELSE
    Distance(iBGMElem)=Distance2
  END IF
  ListDistance(iBGMElem)=ElemID
END DO ! nBGMElems

IF(ALMOSTEQUAL(MAXVAL(Distance),-1.))THEN
  RETURN
END IF

IF(nBGMElems.GT.1) CALL InsertionSort(Distance(1:nBGMElems),ListDistance(1:nBGMElems),nBGMElems)

! loop through sorted list and start by closest element
InElementCheck=.FALSE.
DO iBGMElem=1,nBGMElems
  IF (ALMOSTEQUAL(Distance(iBGMElem),-1.)) CYCLE
  ElemID=ListDistance(iBGMElem)
  IF (.NOT.DoHALO) THEN
    IF (ElemID.GT.PP_nElems) CYCLE
  END IF
  IF (TriaTracking) THEN
    CALL ParticleInsideQuad3D(Pos3D(1:3),ElemID,InElementCheck,Det)
  ELSE
    CALL GetPositionInRefElem(Pos3D(1:3),RefPos,ElemID)
    IF (MAXVAL(ABS(RefPos)).LE.1.0) InElementCheck=.TRUE.
  END IF
  IF (InElementCheck) THEN
    SinglePointToElement = ElemID
    RETURN
  END IF
END DO ! iBGMElem

END FUNCTION SinglePointToElement


END MODULE MOD_Particle_Localization
