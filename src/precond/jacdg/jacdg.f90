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

MODULE MOD_JacDG
!===================================================================================================================================
! Contains the initialization of the DG global variables ! Computes the different DG spatial operators/residuals(Ut) using U
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
INTERFACE InitJacDG
  MODULE PROCEDURE InitJacDG
END INTERFACE

INTERFACE BuildJacDG
  MODULE PROCEDURE BuildJacDG
END INTERFACE

INTERFACE FinalizeJacDG
  MODULE PROCEDURE FinalizeJacDG
END INTERFACE

INTERFACE JacDG
  MODULE PROCEDURE JacDG
END INTERFACE

PUBLIC::InitJacDG,BuildJacDG,FinalizeJacDG,JacDG
!===================================================================================================================================

CONTAINS

SUBROUTINE InitJacDG()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_JacDG_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,iVar,jVar
!===================================================================================================================================
!SWRITE(UNIT_stdOut,'(A)') ' INIT CONSTANTS FOR JACAOBIAN OF DGOPERATOR ...'

ALLOCATE( JacF   (1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)  &
        , JacG   (1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)  &
        , JacH   (1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)  &
        , JacFlux(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N,1:6,1:PP_nElems)     )

ALLOCATE( deltaOld(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N) &
        , deltaNew(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N) )



deltaNew=0.
deltaOld=0.
DO jVar=1,PP_nVar
  DO iVar = 1,PP_nVar
    DO j=0,PP_N
      DO i=0,PP_N
        IF((jVar.GT.iVar).AND.(j.GT.i))THEN
          deltaOld(iVar,jVar,i,j) = 1.
        END IF
        IF((jVar.LT.iVar).AND.(j.LT.i))THEN
          deltaNew(iVar,jVar,i,j) = 1.
        END IF
      END DO ! i
    END DO ! j
  END DO ! iVar
END DO ! jVar

!SWRITE(UNIT_stdOut,'(A)')' JACOBIAN OF DGOPERATOR DONE!'

END SUBROUTINE InitJacDG


SUBROUTINE JacDG(coeff,iElem,Vin,Vout)
!===================================================================================================================================
! Matrix vector product vtilde = JacDG v
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_JacDG_Vars,         ONLY: JacF,JacG,JacH,JacFlux
USE MOD_Jac_ex_Vars,        ONLY: LL_minus, LL_plus
USE MOD_DG_Vars,            ONLY: D_hat
USE MOD_Mesh_Vars,          ONLY:sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: coeff, Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,ll
INTEGER            :: iVar, jVar
!REAL               :: delta(0:PP_N,0:PP_N)
!===================================================================================================================================

!delta=0.
!DO i=0,PP_N
!  delta(i,i)=1.
!END DO

! ! first test
! Vout=0.
! DO iElem=1,PP_nElems
!   DO k=0,PP_N
!     DO j=0, PP_N
!       DO i=0,PP_N
!         DO iVar=1,PP_nVar
!           DO oo=0,PP_N
!             DO nn=0,PP_N
!               DO mm=0,PP_N
!                 DO jVar=1,PP_nVar
!                   Vout(iVar,i,j,k,iElem) = Vout(iVar,i,j,k,iElem)   & ! physical flux
!                                                 + D_hat(i,mm)*delta(j,nn)*delta(k,oo)*JacF(iVar,jVar,i,j,k,iElem)     &
!                                                   * Vin(jVar,mm,nn,oo,iElem)                                          &
!                                                 + Delta(i,mm)*D_hat(j,nn)*delta(k,oo)*JacG(iVar,jVar,i,j,k,iElem)     &
!                                                   * Vin(jVar,mm,nn,oo,iElem)                                          &
!                                                 + Delta(i,mm)*delta(j,nn)*D_hat(k,oo)*JacH(iVar,jVar,i,j,k,iElem)     &
!                                                   *Vin(jVar,mm,nn,oo,iElem) &  ! + numerical flux
!                                          + ( LL_Minus(i,mm)*JacFLUX(iVar,jVar,j,k,XI_MINUS,iElem)                            &
!                                           + LL_Plus  (i,mm)*JacFLUX(iVar,jVar,j,k,XI_PLUS,iElem)  )*delta(k,oo)*delta(j,nn)  &
!                                              * Vin(jVar,mm,nn,oo,iElem)                                              &
!                                          + ( LL_Minus(j,nn)*JacFLUX(iVar,jVar,i,k,ETA_MINUS,iElem)                           &
!                                            + LL_Plus (j,nn)*JacFLUX(iVar,jVar,i,k,ETA_PLUS,iElem) )*delta(i,mm)*delta(k,oo)  &
!                                        * Vin(jVar,mm,nn,oo,iElem)                                              &
!                                          + ( LL_Minus(k,oo)*JacFLUX(iVar,jVar,i,j,ZETA_MINUS,iElem)                          &
!                                            + LL_Plus (k,oo)*JacFLUX(iVar,jVar,i,j,ZETA_PLUS,iElem) )*delta(i,mm)*delta(j,nn) &
!                                                    * Vin(jVar,mm,nn,oo,iElem)
!
!                 END DO ! jVar
!               END DO !mm
!             END DO ! nn
!           END DO ! oo
!         END DO ! ivar
!       END DO ! i
!     END DO ! j
!   END DO ! k
! END DO ! iElem

!! second version
!Vout=0.
!DO iElem=1,PP_nElems
!  DO k=0,PP_N
!    DO j=0, PP_N
!      DO i=0,PP_N
!        DO iVar=1,PP_nVar
!          DO oo=0,PP_N
!            DO nn=0,PP_N
!              DO mm=0,PP_N
!                DO jVar=1,PP_nVar
!                  Vout(iVar,i,j,k,iElem) = Vout(iVar,i,j,k,iElem)   &
!                           + ( (     D_hat(i,mm)*JacF   (iVar,jVar,i,j,k,iElem)                                           &
!                                + LL_Minus(i,mm)*JacFLUX(iVar,jVar,  j,k,XI_MINUS,iElem)                              &
!                                + LL_Plus (i,mm)*JacFLUX(iVar,jVar,  j,k,XI_PLUS, iElem)   )*delta(k,oo)*delta(j,nn)    &
!                             + (     D_hat(j,nn)*JacG   (iVar,jVar,i,j,k,iElem)                                           &
!                                + LL_Minus(j,nn)*JacFLUX(iVar,jVar,i,  k,ETA_MINUS,iElem)                            &
!                                + LL_Plus (j,nn)*JacFLUX(iVar,jVar,i,  k,ETA_PLUS, iElem)  )*delta(i,mm)*delta(k,oo)   &
!                             + (     D_hat(k,oo)*JacH   (iVar,jVar,i,j,k,iElem)                                           &
!                                + LL_Minus(k,oo)*JacFLUX(iVar,jVar,i,j,  ZETA_MINUS,iElem)                            &
!                                + LL_Plus (k,oo)*JacFLUX(iVar,jVar,i,j,  ZETA_PLUS, iElem) )*delta(i,mm)*delta(j,nn) ) &
!                              * Vin(jVar,mm,nn,oo,iElem)
!                END DO ! jVar
!              END DO !mm
!            END DO ! nn
!          END DO ! oo
!        END DO ! ivar
!      END DO ! i
!    END DO ! j
!  END DO ! k
!END DO ! iElem
!
!! third version
!Vout=0.
!DO iElem=1,PP_nElems
!  DO ll=0,PP_N
!    DO k=0,PP_N
!      DO j=0,PP_N
!        DO i=0,PP_N
!          DO jVar=1,PP_nVar
!            DO iVar=1,PP_nVar
!                  Vout(iVar,i,j,k,iElem) = Vout(iVar,i,j,k,iElem)                                                     &
!                             + (     D_hat(i,ll)*JacF   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(i,ll)*JacFLUX(iVar,jVar,  j,k,XI_MINUS,iElem)                              &
!                                + LL_Plus (i,ll)*JacFLUX(iVar,jVar,  j,k,XI_PLUS, iElem)   )* Vin(jVar,ll,j,k,iElem)  &
!                             + (     D_hat(j,ll)*JacG   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_MINUS,iElem)                             &
!                                + LL_Plus (j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_PLUS, iElem)  )* Vin(jVar,i,ll,k,iElem)  &
!                             + (     D_hat(k,ll)*JacH   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_MINUS,iElem)                            &
!                                + LL_Plus (k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_PLUS, iElem) )* Vin(jVar,i,j,ll,iElem)
!            END DO ! iVar
!          END DO ! jvar
!        END DO ! i
!      END DO ! j
!    END DO ! k
!  END DO ! ll
!END DO ! iElem
!
!! and the final stuff
!! I - coeff * DG(v)
!DO iVar=1,PP_nVar
!  Vout(iVar,:,:,:,:) = Vin(iVar,:,:,:,:)-coeff*sJ(:,:,:,:)*Vout(iVar,:,:,:,:)
!END DO ! iVar
!

!! fourth version
!!Vout=0.
!DO iElem=1,PP_nElems
!  DO k=0,PP_N
!    DO j=0,PP_N
!      DO i=0,PP_N
!        DO iVar=1,PP_nVar
!            Vout(iVar,i,j,k,iElem) = 0.
!            DO ll=0,PP_N
!               DO jVar=1,PP_nVar
!                  Vout(iVar,i,j,k,iElem) = Vout(iVar,i,j,k,iElem)                                                     &
!                             + (     D_hat(i,ll)*JacF   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(i,ll)*JacFLUX(iVar,jVar,  j,k,XI_MINUS,iElem)                              &
!                                + LL_Plus (i,ll)*JacFLUX(iVar,jVar,  j,k,XI_PLUS, iElem)   )* Vin(jVar,ll,j,k,iElem)  &
!                             + (     D_hat(j,ll)*JacG   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_MINUS,iElem)                             &
!                                + LL_Plus (j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_PLUS, iElem)  )* Vin(jVar,i,ll,k,iElem)  &
!                             + (     D_hat(k,ll)*JacH   (iVar,jVar,i,j,k,iElem)                                       &
!                                + LL_Minus(k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_MINUS,iElem)                            &
!                                + LL_Plus (k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_PLUS, iElem) )* Vin(jVar,i,j,ll,iElem)
!            END DO ! jvar
!          END DO ! ll
!        END DO ! iVar
!      END DO ! i
!    END DO ! j
!  END DO ! k
!END DO ! iElem

! fiths version
!DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
            Vout(iVar,i,j,k) = 0.
         END DO
      END DO
    END DO
  END DO
!END DO
!DO iElem=1,PP_nElems
  DO ll=0,PP_N
    DO k=0,PP_N
      DO j=0,PP_N
        DO i=0,PP_N
          DO jVar=1,PP_nVar
            DO iVar=1,PP_nVar
                  Vout(iVar,i,j,k) = Vout(iVar,i,j,k)                                                     &
                             + (     D_hat(i,ll)*JacF   (iVar,jVar,i,j,k,iElem)                                       &
                                + LL_Minus(i,ll)*JacFLUX(iVar,jVar,  j,k,XI_MINUS,iElem)                              &
                                + LL_Plus (i,ll)*JacFLUX(iVar,jVar,  j,k,XI_PLUS, iElem)   )* Vin(jVar,ll,j,k)  &
                             + (     D_hat(j,ll)*JacG   (iVar,jVar,i,j,k,iElem)                                       &
                                + LL_Minus(j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_MINUS,iElem)                             &
                                + LL_Plus (j,ll)*JacFLUX(iVar,jVar,i,  k,ETA_PLUS, iElem)  )* Vin(jVar,i,ll,k)  &
                             + (     D_hat(k,ll)*JacH   (iVar,jVar,i,j,k,iElem)                                       &
                                + LL_Minus(k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_MINUS,iElem)                            &
                                + LL_Plus (k,ll)*JacFLUX(iVar,jVar,i,j,  ZETA_PLUS, iElem) )* Vin(jVar,i,j,ll)
            END DO ! jvar
          END DO ! ivar
        END DO ! i
      END DO ! j
    END DO ! k
  END DO ! ll
!END DO ! iElem

! and the final stuff
! I - coeff * DGtilde(v)
! CAUTION DGtilde(v) = - sJ DG(v)
!DO iElem=1,PP_nElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        DO iVar=1,PP_nVar
          Vout(iVar,i,j,k)=Vin(iVar,i,j,k)+coeff*sJ(i,j,k,iElem)*Vout(iVar,i,j,k)
        END DO ! iVar
      END DO ! i
    END DO ! j
  END DO ! k
!END DO ! iElem


END SUBROUTINE JacDG


SUBROUTINE BuildJacDG()
!===================================================================================================================================
! create the constant physical and numerical fluxes for each element
! element local operator in the case of linear equations (Maxwell)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_JacDG_Vars,           ONLY: JacF,JacG,JacH,JacFlux
USE MOD_Mesh_Vars,            ONLY: ElemToSide
USE MOD_Mesh_Vars,            ONLY: nBCSides,BoundaryType,BC
USE MOD_Mesh_Vars,            ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars,            ONLY: nVecLoc, SurfLoc
USE MOD_JacExRiemann,         ONLY: ConstructJacRiemann,ConstructJacBCRiemann
USE MOD_Jacobian,             ONLY: EvalFluxJacobian
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: iLocSide,iElem, SideID, BCType
INTEGER                                             :: i,j,k
REAL,DIMENSION(:,:,:,:),ALLOCATABLE                 :: JacBC
REAL,DIMENSION(:,:),ALLOCATABLE                     :: fJac,gJac,hJac
!===================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' BUILDING CONSTANTS FOR JACDGOPERATOR ...'

ALLOCATE( JacBC(1:PP_nVar,1:PP_nVar,0:PP_N,0:PP_N)  &
        , fJac (1:PP_nVar,1:PP_nVar)                &
        , gJac (1:PP_nVar,1:PP_nVar)                &
        , hJac (1:PP_nVar,1:PP_nVar)                )

DO iElem=1,PP_nElems
  ! internal, physical flux
  CALL EvalFluxJacobian(fJac,gJac,hJac) !fills fJac,gJac,hJac
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        JacF(:,:,i,j,k,iElem) = fJac(:,:)*Metrics_fTilde(1,i,j,k,iElem) + &
                                gJac(:,:)*Metrics_fTilde(2,i,j,k,iElem) + &
                                hJac(:,:)*Metrics_fTilde(3,i,j,k,iElem)
        JacG(:,:,i,j,k,iElem) = fJac(:,:)*Metrics_gTilde(1,i,j,k,iElem) + &
                                gJac(:,:)*Metrics_gTilde(2,i,j,k,iElem) + &
                                hJac(:,:)*Metrics_gTilde(3,i,j,k,iElem)
        JacH(:,:,i,j,k,iElem) = fJac(:,:)*Metrics_hTilde(1,i,j,k,iElem) + &
                                gJac(:,:)*Metrics_hTilde(2,i,j,k,iElem) + &
                                hJac(:,:)*Metrics_hTilde(3,i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO ! k
  !-----------------------------------------------------------------------------------------------
  ! numerical flux || face values
  DO iLocSide = 1,6
    !IF(iLocSide.NE.ZETA_MINUS) CYCLE
    SideID = ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    CALL ConstructJacRiemann(nVecloc(:,:,:,iLocSide,iElem),Surfloc(:,:,iLocSide,iElem),JacFlux(:,:,:,:,iLocSide,iElem))
    IF (SideID.LE.nBCSides) THEN
      BCType=Boundarytype(BC(SideID),BC_TYPE)
      CALL ConstructJacBCRiemann(BCType,nVecLoc(:,:,:,iLocSide,iElem),Surfloc(:,:,iLocSide,iElem),JacBC)
      JacFlux(:,:,:,:,iLocSide,iElem) = JacFlux(:,:,:,:,iLocSide,iElem) + JacBC(:,:,:,:)
    END IF ! BC Side
  END DO ! iLocSide
END DO ! iElems

DEALLOCATE( JacBC,fJac,gJac,hJac)

SWRITE(UNIT_stdOut,'(A)')' BUILD CONSTANTS FOR JACDGOPERATOR DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE BuildJacDG


SUBROUTINE FinalizeJacDG()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_JacDG_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!SDEALLOCATE()
DEALLOCATE(JacF,JacG,JacH,JacFlux)

END SUBROUTINE FinalizeJacDG


END MODULE MOD_JacDG
