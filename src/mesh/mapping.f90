!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_Mappings
!===================================================================================================================================
! build mappings for volume GPs to side GPs, etc
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
INTERFACE InitMappings
  MODULE PROCEDURE InitMappings
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

INTERFACE CGNS_SideToVol2
  MODULE PROCEDURE CGNS_SideToVol2
END INTERFACE

INTERFACE CGNS_VolToSide
  MODULE PROCEDURE CGNS_VolToSide
END INTERFACE

INTERFACE CGNS_VolToSide_IJK
  MODULE PROCEDURE CGNS_VolToSide_IJK
END INTERFACE

INTERFACE SideToVol
  MODULE PROCEDURE SideToVol
END INTERFACE

INTERFACE SideToAdjointLocSide
  MODULE PROCEDURE SideToAdjointLocSide
END INTERFACE

INTERFACE SideToVol2
  MODULE PROCEDURE SideToVol2
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

PUBLIC::InitMappings
PUBLIC::Flip_S2M
PUBLIC::Flip_M2S
PUBLIC::CGNS_SideToVol
PUBLIC::CGNS_SideToVol2
PUBLIC::CGNS_VolToSide
PUBLIC::CGNS_VolToSide_IJK
PUBLIC::SideToVol
PUBLIC::SideToAdjointLocSide
PUBLIC::SideToVol2
PUBLIC::VolToSide
PUBLIC::VolToVol
PUBLIC::ElemToNBElem
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMappings(N_in,V2S,V2SIJK,V2S2,CV2S,S2V,S2V2,CS2V,FS2M)
!===================================================================================================================================
! initialization of all mappings for polynomial degree N_in
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals, ONLY:abort
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: N_in
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(OUT) :: V2S(:,:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: V2SIJK(:,:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: V2S2(:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: CV2S(:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: S2V(:,:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: S2V2(:,:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: CS2V(:,:,:,:)
INTEGER,ALLOCATABLE,INTENT(OUT) :: FS2M(:,:,:,:)     !< FlipSlaveToMaster mapping
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: i,j,k,f,s, ijk(3),pq(3),p,q
!===================================================================================================================================
ALLOCATE(V2S(3,0:N_in,0:N_in,0:N_in,0:4,1:6))
ALLOCATE(V2SIJK(3,0:N_in,0:N_in,0:N_in,0:4,1:6))
ALLOCATE(V2S2(2,0:N_in,0:N_in,0:4,1:6))
ALLOCATE(CV2S(3,0:N_in,0:N_in,0:N_in,1:6))
ALLOCATE(S2V(3,0:N_in,0:N_in,0:N_in,0:4,1:6))
ALLOCATE(S2V2(2,0:N_in,0:N_in,0:4,1:6))
ALLOCATE(CS2V(2,0:N_in,0:N_in,1:6))
ALLOCATE(FS2M(2,0:N_in,0:N_in,0:4))

DO k=0,N_in; DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    DO s=1,6
      V2S(:,i,j,k,f,s) = VolToSide(i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

DO k=0,N_in; DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    DO s=1,6
      V2SIJK(:,i,j,k,f,s) = VolToSideIJK(i,j,k,s)
    END DO
  END DO
END DO; END DO; END DO

DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    DO s=1,6
      V2S2(:,i,j,f,s) = VolToSide2(i,j,f,s)
    END DO
  END DO
END DO; END DO

DO k=0,N_in; DO j=0,N_in; DO i=0,N_in
  DO s=1,6
    CV2S(:,i,j,k,s) = CGNS_VolToSide(i,j,k,s)
  END DO
END DO; END DO; END DO

DO k=0,N_in; DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    DO s=1,6
      S2V(:,i,j,k,f,s) = SideToVol(i,j,k,f,s)
    END DO
  END DO
END DO; END DO; END DO

DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    DO s=1,6
      S2V2(:,i,j,f,s) = SideToVol2(N_in,i,j,f,s)
    END DO
  END DO
END DO; END DO

DO j=0,N_in; DO i=0,N_in
  DO s=1,6
    CS2V(:,i,j,s) = CGNS_SideToVol2(N_in,i,j,s)
  END DO
END DO; END DO

! Flip_S2M
DO j=0,N_in; DO i=0,N_in
  DO f=0,4
    FS2M(:,i,j,f) = Flip_S2M(N_in,i,j,f)
  END DO
END DO; END DO

DO f = 0, 4
  DO s = 1, 6
    DO q = 0,N_in; DO p = 0,N_in
      ijk = S2V(:,0,p,q,f,s)
      pq = V2S(:,ijk(1),ijk(2),ijk(3),f,s)
      IF ((pq(1).NE.p).OR.(pq(2).NE.q)) THEN
CALL abort(&
__STAMP__&
,'SideToVol does not fit to VolToSideA')
      END IF
    END DO; END DO
  END DO ! s = 1, 6
END DO ! f = 0, 4

END SUBROUTINE InitMappings


FUNCTION Flip_S2M(N_in,p, q, flip)
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
INTEGER,INTENT(IN) :: N_in,p,q,flip
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
    Flip_S2M = (/N_in-p,     q/)
  CASE(3)
    Flip_S2M = (/N_in-q,N_in-p/)
  CASE(4)
    Flip_S2M = (/     p,N_in-q/)
END SELECT
END FUNCTION Flip_S2M


FUNCTION Flip_M2S(N_in,p, q, flip)
!===================================================================================================================================
! Transforms Coordinates from RHS of Master to RHS of Slave
!   actualy this is the same function as Flip_S2M
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in,p,q,flip
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: Flip_M2S
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Flip_M2S=Flip_S2M(N_in,p,q,flip)
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


FUNCTION CGNS_VolToSide_IJK(i,j,k, locSideID)
!===================================================================================================================================
! Transforms Volume-Coordinates into RHS of the Side (uses CGNS-Notation)
! input: i,j,k, locSideID
!   where: i,j,k = volume-indices
! output: indices in IJK of volume and  volume-index which is not used (depending on locSideID)
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
INTEGER,DIMENSION(3) :: CGNS_VolToSide_IJK
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(locSideID)
  CASE(XI_MINUS)
    CGNS_VolToSide_IJK = (/k,j,i/)
  CASE(XI_PLUS)
    CGNS_VolToSide_IJK= (/j,k,i/)
  CASE(ETA_MINUS)
    CGNS_VolToSide_IJK= (/i,k,j/)
  CASE(ETA_PLUS)
    CGNS_VolToSide_IJK = (/i,k,j/)
  CASE(ZETA_MINUS)
    CGNS_VolToSide_IJK = (/j,i,k/)
  CASE(ZETA_PLUS)
    CGNS_VolToSide_IJK = (/i,j,k/)
END SELECT
END FUNCTION CGNS_VolToSide_IJK


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


FUNCTION CGNS_SideToVol2(N_in,p, q, locSideID)
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
INTEGER,INTENT(IN) :: p,q,locSideID,N_in
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
    CGNS_SideToVol2 = (/N_in-p,q/)
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
VolToSide(1:2) = Flip_S2M(PP_N,pq(1),pq(2),flip)
VolToSide(3) = pq(3)
END FUNCTION VolToSide

FUNCTION VolToSideIJK(i,j,k,locSideID)
!===================================================================================================================================
! Transform Volume-Coordinates to RHS-Coordinates of Master. This is: VolToSide = Flip_S2M(CGNS_VolToSide(...))
! input: i,j,k, flip, locSideID
!   where: i,j,k = volume-indices
! output: indices in IJK system, but on the side
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
INTEGER,DIMENSION(3) :: VolToSideIJK
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(locSideID)
  CASE(XI_MINUS,XI_PLUS)
    VolToSideIJK = (/j,k,i/)
  CASE(ETA_MINUS,ETA_PLUS)
    VolToSideIJK = (/i,k,j/)
  CASE(ZETA_MINUS,ZETA_PLUS)
    VolToSideIJK = (/i,j,k/)
END SELECT
END FUNCTION VolToSideIJK

FUNCTION VolToSide2(ijk1,ijk2, flip, locSideID)
!===================================================================================================================================
! Transform Volume-Coordinates to RHS-Coordinates of Master. This is: VolToSide = Flip_S2M(CGNS_VolToSide(...))
! input: (ijk1,ijk2)is i,j for Zeta, ik, for eta, and jk for xi, flip, locSideID
!   where: i,j,k = volume-indices
! output: indices in Master-RHS
!===================================================================================================================================
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: ijk1,ijk2,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: VolToSide2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(3) :: pq
!===================================================================================================================================
SELECT CASE(locSideID)
CASE(XI_MINUS,XI_PLUS)
  pq = CGNS_VolToSide(0,ijk1,ijk2,locSideID)
CASE(ETA_MINUS,ETA_PLUS)
  pq = CGNS_VolToSide(ijk1,0,ijk2,locSideID)
CASE(ZETA_MINUS,ZETA_PLUS)
  pq = CGNS_VolToSide(ijk1,ijk2,0,locSideID)
END SELECT
VolToSide2(1:2) = Flip_S2M(PP_N,pq(1),pq(2),flip)
END FUNCTION VolToSide2

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
pq = Flip_M2S(PP_N,p,q,flip)
SideToVol = CGNS_SideToVol(l,pq(1),pq(2),locSideID)
END FUNCTION SideToVol


FUNCTION SideToAdjointLocSide(p, q, flip, locSideID)
!===================================================================================================================================
! Transform RHS-Coordinates of Master to Volume-Coordinates. This is: SideToVol = CGNS_SideToVol(Flip_M2S(...))
! input: l, p,q, flip, locSideID
!     where: p,q are in Master-RHS;
!            l is the xi-,eta- or zeta-index in 0:PP_N corresponding to locSideID
! output: volume-indicies
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals, ONLY:abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN) :: p,q,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: SideToAdjointLocSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!===================================================================================================================================

pq = SideToVol2(PP_N,p, q, flip, locSideID)
SELECT CASE(locSideID)
CASE(XI_MINUS)
  IF(pq(1).EQ.0)THEN
    SideToAdjointLocSide(1)=ETA_MINUS
  ELSE
    SideToAdjointLocSide(1)=ETA_PLUS
  END IF
  IF(pq(2).EQ.0)THEN
    SideToAdjointLocSide(2)=ZETA_MINUS
  ELSE
    SideToAdjointLocSide(2)=ZETA_PLUS
  END IF
CASE(XI_PLUS)
  IF(pq(1).EQ.0)THEN
    SideToAdjointLocSide(1)=ETA_MINUS
  ELSE
    SideToAdjointLocSide(1)=ETA_PLUS
  END IF
  IF(pq(2).EQ.0)THEN
    SideToAdjointLocSide(2)=ZETA_MINUS
  ELSE
    SideToAdjointLocSide(2)=ZETA_PLUS
  END IF
CASE DEFAULT
CALL abort(&
__STAMP__&
,'SideToAdjointLocSide only supports XI_PLUS,XI_MINUS')
END SELECT

END FUNCTION SideToAdjointLocSide


FUNCTION SideToVol2(N_in,p, q, flip, locSideID)
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
INTEGER,INTENT(IN) :: N_in,p,q,flip,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,DIMENSION(2) :: SideToVol2
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,DIMENSION(2) :: pq
!===================================================================================================================================
pq = Flip_M2S(N_in,p,q,flip)
SideToVol2 = CGNS_SideToVol2(N_in,pq(1),pq(2),locSideID)
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
INTEGER              :: flip, SideID
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


FUNCTION VolToVol(i,j,k,locSideID,iElem)
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
INTEGER,INTENT(IN) :: i,j,k,locSideID,iElem
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

END MODULE MOD_Mappings
