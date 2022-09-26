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

MODULE MOD_JacSurfInt
!===================================================================================================================================
! Contains the computation of the local jacobian numerical flux
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE DGJacSurfInt
  MODULE PROCEDURE DGJacSurfInt
END INTERFACE

INTERFACE DGJacSurfInt1D
  MODULE PROCEDURE DGJacSurfInt1D
END INTERFACE

INTERFACE DGJacSurfInt_Neighbor
  MODULE PROCEDURE DGJacSurfInt_Neighbor
END INTERFACE

PUBLIC::DGJacSurfInt,DGJacSurfInt1D
PUBLIC::DGJacSurfInt_Neighbor
!===================================================================================================================================

CONTAINS

SUBROUTINE DGJacSurfInt(BJ,iElem)
!===================================================================================================================================
! Contains the computation of the local jacobian of the Flux Vector Splitting scheme
! computation is done for one element!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars           ,ONLY: ElemToSide
USE MOD_Mesh_Vars           ,ONLY: BC,BoundaryType,nBCSides
USE MOD_Mesh_Vars           ,ONLY: nVecLoc, SurfLoc
USE MOD_LinearSolver_Vars   ,ONLY: nDOFElem
USE MOD_Jac_Ex_Vars         ,ONLY: LL_minus, LL_plus
USE MOD_JacExRiemann        ,ONLY: ConstructJacRiemann,ConstructJacBCRiemann,ConstructJacRiemannDielectric
USE MOD_Dielectric_vars     ,ONLY: DoDielectric,isDielectricElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: BJ(nDOFElem,nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: SideID
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,6) :: JacA
REAL,DIMENSION(0:PP_N,0:PP_N)                       :: delta
INTEGER                                             :: i,j,k,mm,nn,oo,r,s
INTEGER                                             :: vn1, vn2
INTEGER                                             :: iLocSide
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)   :: JacBC
INTEGER                                             :: BCType
!===================================================================================================================================

vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)

delta = 0.
DO i = 0,PP_N
  delta(i,i) = 1.
END DO ! i

! debug
JacA = 0.
JacBC = 0.

DO iLocSide = 1,6
  SideID = ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  IF(DoDielectric) THEN
    ! if it is a dielectric element, the preconditioner has to consider the it on the face
    IF(isDielectricElem(iElem))THEN !  master is DIELECTRIC and slave PHYSICAL
      CALL ConstructJacRiemannDielectric(nVecLoc(:,:,:,iLocSide,iElem) &
                                        ,SurfLoc(  :,:,iLocSide,iElem) &
                                        ,JacA   (:,:,:,:,iLocSide),ilocSide,iElem)
      CYCLE
    END IF
  END IF
  ! get derivative of flux maxtrix
  CALL ConstructJacRiemann(nVecLoc(:,:,:,iLocSide,iElem),SurfLoc(:,:,iLocSide,iElem),JacA(:,:,:,:,iLocSide))
  IF (SideID.LE.nBCSides) THEN
    BCType=Boundarytype(BC(SideID),BC_TYPE)
    CALL ConstructJacBCRiemann(BCType,nVecLoc(:,:,:,iLocSide,iElem),SurfLoc(:,:,iLocSide,iElem),JacBC)
    JacA(:,:,:,:,iLocSide) = JacA(:,:,:,:,iLocSide) + JacBC(:,:,:,:)
  END IF
END DO ! iLocSide

 ! loop over columns: oo,mm,nn -> s
 DO oo = 0,PP_N
   DO nn = 0,PP_N
     DO mm = 0,PP_N
       s = vn2 * oo + vn1 * nn + PP_nVar * mm
       ! loop over rows: i,j,k  -> r
       DO k = 0,PP_N
          DO j = 0,PP_N
            DO i = 0,PP_N
              r = vn2 * k + vn1 * j + PP_nVar * i
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                              + ( LL_Minus(i,mm)*JacA(:,:,j,k,XI_MINUS) &
                                                + LL_Plus (i,mm)*JacA(:,:,j,k,XI_PLUS)  )*delta(k,oo)*delta(j,nn) &
                                              + ( LL_Minus(j,nn)*JacA(:,:,i,k,ETA_MINUS) &
                                                + LL_Plus (j,nn)*JacA(:,:,i,k,ETA_PLUS)  )*delta(i,mm)*delta(k,oo) &
                                              + ( LL_Minus(k,oo)*JacA(:,:,i,j,ZETA_MINUS) &
                                                + LL_Plus (k,oo)*JacA(:,:,i,j,ZETA_PLUS)  )*delta(i,mm)*delta(j,nn)
            END DO ! i
           END DO ! j
         END DO ! k
     END DO ! nn
   END DO ! mm
 END DO ! oo

END SUBROUTINE DGJacSurfInt

SUBROUTINE DGJacSurfInt1D(dRdXi,dRdEta,dRdZeta,iElem)
!===================================================================================================================================
! Contains the computation of the local jacobian of the Flux Vector Splitting scheme
! computation is done for one element!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars           ,ONLY: ElemToSide
USE MOD_Mesh_Vars           ,ONLY: BC,BoundaryType,nBCSides
USE MOD_Mesh_Vars           ,ONLY: nVecLoc, SurfLoc
USE MOD_LinearSolver_Vars   ,ONLY: nDOFLine
USE MOD_Jac_Ex_Vars         ,ONLY: LL_minus, LL_plus
USE MOD_JacExRiemann        ,ONLY: ConstructJacRiemann,ConstructJacBCRiemann,ConstructJacRiemannDielectric
USE MOD_Dielectric_vars     ,ONLY: DoDielectric,isDielectricElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: dRdXi  (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT)                                  :: dRdEta (1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
REAL,INTENT(INOUT)                                  :: dRdZeta(1:nDOFLine,1:nDOFLine,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: SideID
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,6) :: JacA
INTEGER                                             :: i,j,k,mm,nn,oo,vni,vnj,vnk,s
INTEGER                                             :: iLocSide
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)   :: JacBC
INTEGER                                             :: BCType
!===================================================================================================================================

! debug
JacA = 0.
JacBC = 0.

DO iLocSide = 1,6
  SideID = ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  ! check for dielectric
  IF(DoDielectric) THEN
    ! if it is a dielectric element, the preconditioner has to consider the it on the face
    IF(isDielectricElem(iElem))THEN !  master is DIELECTRIC and slave PHYSICAL
      CALL ConstructJacRiemannDielectric(nVecloc(:,:,:,iLocSide,iElem),Surfloc(:,:,iLocSide,iElem) &
                                        ,JacA(:,:,:,:,iLocSide),ilocSide,iElem)
      CYCLE
    END IF
  END IF
  ! normal element
  CALL ConstructJacRiemann(nVecLoc(:,:,:,iLocSide,iElem),SurfLoc(:,:,iLocSide,iElem),JacA(:,:,:,:,iLocSide))
  IF (SideID.LE.nBCSides) THEN
   BCType=Boundarytype(BC(SideID),BC_TYPE)
   CALL ConstructJacBCRiemann(BCType,nVecLoc(:,:,:,iLocSide,iElem),SurfLoc(:,:,iLocSide,iElem),JacBC)
   JacA(:,:,:,:,iLocSide) = JacA(:,:,:,:,iLocSide) + JacBC(:,:,:,:)
  END IF
END DO ! iLocSide


DO oo = 0,PP_N
  DO nn = 0,PP_N
    DO mm = 0,PP_N
      s=mm*PP_nVar
      DO i=0,PP_N
        vni=PP_nVar*i
        dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,nn,oo) = dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,nn,oo)  &
                                              + ( LL_Minus(i,mm)*JacA(1:PP_nVar,1:PP_nVar,nn,oo,XI_MINUS) &
                                              +   LL_Plus (i,mm)*JacA(1:PP_nVar,1:PP_nVar,nn,oo,XI_PLUS)  )
      END DO ! i

      s=nn*PP_nVar
      DO j=0,PP_N
        vnj=PP_nVar*j
        dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,mm,oo) = dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,mm,oo)   &
                                            + ( LL_Minus(j,nn)*JacA(1:PP_nVar,1:PP_nVar,mm,oo,ETA_MINUS) &
                                              + LL_Plus (j,nn)*JacA(1:PP_nVar,1:PP_nVar,mm,oo,ETA_PLUS)  )
      END DO ! j

      s=oo*PP_nVar
      DO k=0,PP_N
        vnk=PP_nVar*k
        dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,mm,nn) = dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,mm,nn)   &
                                            + ( LL_Minus(k,oo)*JacA(1:PP_nVar,1:PP_nVar,mm,nn,ZETA_MINUS) &
                                              + LL_Plus (k,oo)*JacA(1:PP_nVar,1:PP_nVar,mm,nn,ZETA_PLUS)  )
      END DO ! k
    END DO ! nn
  END DO ! mm
END DO ! oo



!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        vni=PP_nVar*i
!        vnj=PP_nVar*j
!        vnk=PP_nVar*k
!        s=0
!        DO l=0,PP_N
!          dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,j,k) = dRdXi  (vni+1:vni+PP_nVar,s+1:s+PP_nVar,j,k)  &
!                                              + ( LL_Minus(i,l)*JacA(1:PP_nVar,1:PP_nVar,j,k,XI_MINUS) &
!                                                + LL_Plus (i,l)*JacA(1:PP_nVar,1:PP_nVar,j,k,XI_PLUS)  )
!
!          dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,i,k) = dRdEta (vnj+1:vnj+PP_nVar,s+1:s+PP_nVar,i,k)   &
!                                              + ( LL_Minus(j,l)*JacA(1:PP_nVar,1:PP_nVar,i,k,ETA_MINUS) &
!                                                + LL_Plus (j,l)*JacA(1:PP_nVar,1:PP_nVar,i,k,ETA_PLUS)  )
!
!          dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,i,j) = dRdZeta(vnk+1:vnk+PP_nVar,s+1:s+PP_nVar,i,j) &
!                                              + ( LL_Minus(k,l)*JacA(1:PP_nVar,1:PP_nVar,i,j,ZETA_MINUS) &
!                                                + LL_Plus (k,l)*JacA(1:PP_nVar,1:PP_nVar,i,j,ZETA_PLUS)  )
!
!          s=s+PP_nVar
!        END DO ! l
!      END DO ! i
!    END DO ! j
!  END DO ! k

END SUBROUTINE DGJacSurfInt1D

#ifdef DoNotCompileThis
SUBROUTINE DGJacSurfInt(BJ,iElem)
!===================================================================================================================================
! Debug version only for xi minus, xi plus, eta minus
! Contains the computation of the local jacobian of the Flux Vector Splitting scheme
! computation is done for one element!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: nBCSides
USE MOD_Mesh_Vars       ,ONLY: ElemToSide
USE MOD_Mesh_Vars,       ONLY: BC,BoundaryType
USE MOD_Mesh_Vars       ,ONLY: nVecLoc, SurfLoc
USE MOD_LinearSolver_Vars   ,ONLY: nDOFElem
USE MOD_Jac_Ex_Vars     ,ONLY: LL_minus, LL_plus
USE MOD_JacExRiemann    ,ONLY: ConstructJacRiemann,ConstructJacBCRiemann
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: BJ(nDOFElem,nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: SideID, flip
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)   :: JacA
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)   :: JacBC
REAL,DIMENSION(0:PP_N,0:PP_N)                       :: delta
INTEGER                                             :: i,j,k,mm,nn,oo,r,s,v,w
INTEGER                                             :: vn1, vn2
INTEGER                                             :: BCType
!===================================================================================================================================

     vn1 = PP_nVar * (PP_N + 1)
     vn2 = vn1 * (PP_N +1)


delta = 0.
DO i = 0,PP_N
  delta(i,i) = 1.
END DO ! i

! this time a different approach
! we are going side by side
! reason ElemToSide contains the correct flip

!-----------------------------------------------------------------------------------------------------------------------
!             XI MINUS
!-----------------------------------------------------------------------------------------------------------------------

! SideID = ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem)
! flip   = ElemToSide(E2S_FLIP,XI_MINUS,iELEM)
! ! get derivative of flux maxtrix
! print*,'ximinus',flip
! CALL ConstructJacRiemann(flip,nVec(:,:,:,XI_MINUS,iElem),Surf(:,:,XI_MINUS,iElem),JacA)
! IF (SideID.LE.nBCSides) THEN
!  BCType=Boundarytype(BC(SideID),BC_TYPE)
!  print*,'BC Side, BcType',BcType
!  CALL ConstructJacBCRiemann(BCType,nVec(:,:,:,XI_MINUS,iElem),Surf(:,:,XI_MINUS,iElem),JacBC)
!  JacA = JacA + JacBC
! END IF
!
! DO oo = 0,PP_N
!   DO nn = 0,PP_N
!     DO mm = 0,PP_N
!       s = vn2 * oo + vn1 * nn + PP_nVar * mm
!       DO k = 0,PP_N
!          DO j = 0,PP_N
!            DO i = 0,PP_N
!              r = vn2 * k + vn1 * j + PP_nVar * i
!              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = &
!                        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + LL_Minus(i,mm) * JacA(:,:,j,k) *delta(k,oo)*delta(j,nn)
!            END DO ! i
!           END DO ! j
!         END DO ! k
!     END DO ! nn
!   END DO ! mm
! END DO ! oo

!-----------------------------------------------------------------------------------------------------------------------
!             XI PLUS
!-----------------------------------------------------------------------------------------------------------------------

! SideID = ElemToSide(E2S_SIDE_ID,XI_PLUS,iElem)
! flip   = ElemToSide(E2S_FLIP,   XI_PLUS,iELEM)
! ! get derivative of flux maxtrix
! print*,'xiplus',flip
! CALL ConstructJacRiemann(flip,nVec(:,:,:,XI_PLUS,iElem),Surf(:,:,XI_PLUS,iElem),JacA)
! IF (SideID.LE.nBCSides) THEN
!  BCType=Boundarytype(BC(SideID),BC_TYPE)
!  print*,'BC Side, BcType',BcType
!  CALL ConstructJacBCRiemann(BCType,nVec(:,:,:,XI_PLUS,iElem),Surf(:,:,XI_PLUS,iElem),JacBC)
!  JacA = JacA + JacBC
! END IF
!
! DO oo = 0,PP_N
!   DO nn = 0,PP_N
!     DO mm = 0,PP_N
!       s = vn2 * oo + vn1 * nn + PP_nVar * mm
!       DO k = 0,PP_N
!          DO j = 0,PP_N
!            DO i = 0,PP_N
!              r = vn2 * k + vn1 * j + PP_nVar * i
!              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = &
!                        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + LL_Plus(i,mm) * JacA(:,:,j,k) *delta(k,oo)*delta(j,nn)
!            END DO ! i
!           END DO ! j
!         END DO ! k
!     END DO ! nn
!   END DO ! mm
! END DO ! oo

!-----------------------------------------------------------------------------------------------------------------------
!             Eta MINUS
!-----------------------------------------------------------------------------------------------------------------------

 SideID = ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem)
 flip   = ElemToSide(E2S_FLIP,   ETA_MINUS,iELEM)
 ! get derivative of flux maxtrix
 print*,'eta minus',flip
 CALL ConstructJacRiemann(flip,nVecLoc(:,:,:,ETA_MINUS,iElem),SurfLoc(:,:,ETA_MINUS,iElem),JacA)
 IF (SideID.LE.nBCSides) THEN
   BCType=Boundarytype(BC(SideID),BC_TYPE)
   print*,'BC Side, BcType',BcType
   CALL ConstructJacBCRiemann(BCType,nVecLoc(:,:,:,ETA_MINUS,iElem),SurfLoc(:,:,ETA_MINUS,iElem),JacBC)
   JacA = JacA + JacBC
 END IF

 DO oo = 0,PP_N
   DO nn = 0,PP_N
     DO mm = 0,PP_N
       s = vn2 * oo + vn1 * nn + PP_nVar * mm
       DO k = 0,PP_N
          DO j = 0,PP_N
            DO i = 0,PP_N
              r = vn2 * k + vn1 * j + PP_nVar * i
              BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = &
                        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) + LL_Minus(j,nn) * JacA(:,:,i,k) *delta(k,oo)*delta(i,mm)
            END DO ! i
           END DO ! j
         END DO ! k
     END DO ! nn
   END DO ! mm
 END DO ! oo

END SUBROUTINE DGJacSurfInt
#endif /* do not compile, only for debug*/


SUBROUTINE DGJacSurfInt_Neighbor(BJ,locSideID,ElemID)
!===================================================================================================================================
! Contains the computation of the local jacobian of the Flux Vector Splitting scheme
! computation is done for one element!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars              ,ONLY: ElemToSide,SideToElem
USE MOD_Mesh_Vars              ,ONLY: nBCSides
USE MOD_Mappings               ,ONLY: VolToSide,SideToVol
USE MOD_Mesh_Vars              ,ONLY: nVecLoc, SurfLoc
USE MOD_Precond_Vars           ,ONLY: neighborElemID
USE MOD_LinearSolver_Vars      ,ONLY: nDOFElem
USE MOD_Interpolation_Vars     ,ONLY: L_PlusMinus
USE MOD_CSR_Vars               ,ONLY: L_HatPlusMinus
USE MOD_JacExRiemann           ,ONLY: ConstructJacNeighborRiemann,ConstructJacBCRiemann
USE MOD_Mesh_Vars              ,ONLY: sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                  :: ElemID,locSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: BJ(nDOFElem,nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: SideID,neighbor_locSideID, flip, neighbor_flip, neighbor_ElemID
REAL,DIMENSION(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)   :: JacA
REAL,DIMENSION(0:PP_N,0:PP_N)                       :: delta
INTEGER                                             :: i,j,k,mm,nn,oo,r,s
INTEGER,DIMENSION(3)                                :: pq,mn
INTEGER                                             :: vn1, vn2
!===================================================================================================================================

vn1 = PP_nVar * (PP_N + 1)
vn2 = vn1 * (PP_N +1)

delta = 0.
DO i = 0,PP_N
  delta(i,i) = 1.
END DO ! i

! debug
JacA= 0.
BJ=0.

! get SideID
SideID = ElemToSide(E2S_SIDE_ID,locSideID,ElemID)
! my flip
flip = ElemToSide(E2S_FLIP,locSideID,ElemID)

! get derivative of flux maxtrix
IF (SideID.LE.nBCSides) RETURN
CALL ConstructJacNeighborRiemann(nVecLoc(:,:,:,locSideID,ElemID),SurfLoc(:,:,locSideID,ElemID),JacA(:,:,:,:))

! in SideTOElem, the information is stored for the master-side, the slave side is the neighbor side !!!!
IF(flip.EQ.0)THEN
    ! SideID of slave
  neighbor_locSideID=SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  neighbor_flip     =SideToElem(S2E_FLIP,SideID)
  neighbor_ElemID   =SideToElem(S2E_NB_ELEM_ID,SideID)
ELSE
  ! SideID of master
  neighbor_locSideID=SideToElem(S2E_LOC_SIDE_ID,SideID)
  neighbor_flip     =0
  neighbor_ElemID   =SideToElem(S2E_ELEM_ID,SideID)
END IF

! matrix with cell connection
neighborElemID(locSideID,ElemID)=neighbor_ElemID


#if (PP_NodeType==1) /* Gauss Points */
SELECT CASE(locSideID)
CASE(XI_MINUS,XI_PLUS)
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO mm=0,PP_N
            s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mm
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(mm,neighbor_locSideID)*JacA(:,:,j,k)

          END DO
        END DO
      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO nn=0,PP_N
            s =vn2 *mn(3) + vn1 * nn + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(nn,neighbor_locSideID)*JacA(:,:,j,k)

          END DO
        END DO
      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO oo=0,PP_N
            s =vn2 *oo + vn1 * mn(2) + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(oo,neighbor_locSideID)*JacA(:,:,j,k)

          END DO
        END DO
      END DO
    END DO
  END SELECT ! neighbor_locSideID
CASE(ETA_MINUS,ETA_PLUS)
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in j
    DO k=0,PP_N
      DO i=0,PP_N
        DO j=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO mm=0,PP_N
            s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mm
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                  -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(mm,neighbor_locSideID)*JacA(:,:,i,k)

          END DO
        END DO
      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in k
    DO k=0,PP_N
      DO i=0,PP_N
        DO j=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO nn=0,PP_N
            s =vn2 *mn(3) + vn1 * nn + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                  -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(nn,neighbor_locSideID)*JacA(:,:,i,k)

          END DO
        END DO
      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in k
    DO k=0,PP_N
      DO i=0,PP_N
        DO j=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO oo=0,PP_N
            s =vn2 *oo + vn1 * mn(2) + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(oo,neighbor_locSideID)*JacA(:,:,i,k)

          END DO
        END DO
      END DO
    END DO
  END SELECT ! neighbor_locSideID
CASE(ZETA_MINUS,ZETA_PLUS)
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        DO k=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO mm=0,PP_N
            s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mm
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(mm,neighbor_locSideID)*JacA(:,:,i,j)

          END DO
        END DO
      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        DO k=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO nn=0,PP_N
            s =vn2 *mn(3) + vn1 * nn + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(nn,neighbor_locSideID)*JacA(:,:,i,j)

          END DO
        END DO
      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        DO k=0,PP_N
          ! mapp the derivative to the correct matrix intix
          ! prolongation in k
          ! matrix index - ijk
          r = vn2 * k + vn1 * j + PP_nVar * i
          ! positionon master side
          pq = VolToSide(i,j,k, flip, locSideID)
          ! loop over second element for prolongation
          mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
          DO oo=0,PP_N
            s =vn2 *oo + vn1 * mn(2) + PP_nVar * mn(1)
            BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                                      -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(oo,neighbor_locSideID)*JacA(:,:,i,j)

          END DO
        END DO
      END DO
    END DO
  END SELECT ! neighbor_locSideID
END SELECT ! locSideID
#else
! Gauss-Lobatto
SELECT CASE(locSideID)
CASE(XI_MINUS,XI_PLUS)
  IF(locSideID.EQ.XI_MINUS) THEN
    i=0
  ELSE   IF(locSideID.EQ.XI_PLUS) THEN
    i=PP_N
  END IF
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(mn(1),neighbor_locSideID)*JacA(:,:,j,k)

      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(mn(2),neighbor_locSideID)*JacA(:,:,j,k)

      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in i
    DO k=0,PP_N
      DO j=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(i,locSideID)*L_plusminus(mn(3),neighbor_locSideID)*JacA(:,:,j,k)

      END DO
    END DO
  END SELECT ! neighbor_locSideID
CASE(ETA_MINUS,ETA_PLUS)
  IF(locSideID.EQ.ETA_MINUS) THEN
    j=0
  ELSE   IF(locSideID.EQ.ETA_PLUS) THEN
    j=PP_N
  END IF
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in j
    DO k=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(mn(1),neighbor_locSideID)*JacA(:,:,i,k)
      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in j
    DO k=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(mn(2),neighbor_locSideID)*JacA(:,:,i,k)
      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in j
    DO k=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(j,locSideID)*L_plusminus(mn(3),neighbor_locSideID)*JacA(:,:,i,k)
      END DO
    END DO
  END SELECT ! neighbor_locSideID
CASE(ZETA_MINUS,ZETA_PLUS)
  IF(locSideID.EQ.ZETA_MINUS) THEN
    k=0
  ELSE IF(locSideID.EQ.ZETA_PLUS) THEN
    k=PP_N
  END IF
  SELECT CASE(Neighbor_locSideID)
  CASE(XI_MINUS,XI_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                               -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(mn(1),neighbor_locSideID)*JacA(:,:,i,j)
      END DO
    END DO
  CASE(ETA_MINUS,ETA_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                           -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(mn(2),neighbor_locSideID)*JacA(:,:,i,j)
      END DO
    END DO
  CASE(ZETA_MINUS,ZETA_PLUS)
    ! prolongation in k
    DO j=0,PP_N
      DO i=0,PP_N
        ! mapp the derivative to the correct matrix intix
        ! prolongation in k
        ! matrix index - ijk
        r = vn2 * k + vn1 * j + PP_nVar * i
        ! positionon master side
        pq = VolToSide(i,j,k, flip, locSideID)
        ! loop over second element for prolongation
        mn =SideToVol(pq(3),pq(1),pq(2), neighbor_flip,neighbor_locSideID)
        s =vn2 *mn(3) + vn1 * mn(2) + PP_nVar * mn(1)
        BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) = BJ(r+1:r+PP_nVar,s+1:s+PP_nVar) &
                              -sJ(i,j,k,ElemID)*L_HatPlusMinus(k,locSideID)*L_plusminus(mn(3),neighbor_locSideID)*JacA(:,:,i,j)
      END DO
    END DO
  END SELECT ! neighbor_locSideID
END SELECT ! locSideID
#endif

END SUBROUTINE DGJacSurfInt_Neighbor


END MODULE MOD_JacSurfInt
