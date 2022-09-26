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

MODULE MOD_ApplyPreconditioner
!===================================================================================================================================
! Module for the Block-Jacobi Preconditioner
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE Preconditioner
  MODULE PROCEDURE Preconditioner
END INTERFACE

PUBLIC :: Preconditioner
!===================================================================================================================================

CONTAINS

SUBROUTINE Preconditioner(V,Vprecond)
!===================================================================================================================================
! select preconditioner
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars,           ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)  :: V(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL,INTENT(INOUT)  :: Vprecond(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SELECT CASE(PrecondType)
  CASE(1,2) ! Finite-Differences and analytical Preconditioner
    ! BJ per Element
    CALL ApplyPrecond_Elem(V,Vprecond)
  CASE(3)
    CALL ApplyILU(V,Vprecond)
  CASE(4)
    CALL ApplyBILU0_BCSR(V,Vprecond)
  CASE(201)
    CALL ApplyADI(V,Vprecond)
  CASE DEFAULT
    Vprecond=V
END SELECT

END SUBROUTINE Preconditioner

SUBROUTINE ApplyPrecond_Elem(v,z)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars ,ONLY:nDOFelem
USE MOD_Precond_Vars      ,ONLY:invP,PrecondMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: v(nDOFelem,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: z(nDOFelem,PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: r,s,iElem
!===================================================================================================================================

SELECT CASE(PrecondMethod)
  CASE(0) ! using a loop
    ! loop over all elements
    DO iElem=1,PP_nElems
      z(:,iElem) = 0.
      DO s=1,nDOFElem
        DO r=1,nDOFElem
          z(r,iElem) = z(r,iElem)+invP(r,s,iElem)*v(s,iElem)
        END DO ! r
      END DO ! s
    END DO ! iElem
  CASE(1) ! Matmul
    DO iElem=1,PP_nElems
      z(:,iElem) = MATMUL(invP(:,:,iElem),v(:,iElem))
    END DO !iElem
END SELECT

END SUBROUTINE  ApplyPrecond_Elem


SUBROUTINE ApplyILU(Vin,Vout)
!==================================================================================================================================
! Application of block ILU0 preconditioner
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY: nDOFelem
USE MOD_CSR_Vars               ,ONLY: DE,IL,IU
USE MOD_CSR_Vars               ,ONLY: nMTriangle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                                    :: Vin(1:nDOFElem,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                      :: Vout(1:nDOFElem,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                               :: Vcalc(1:nDOFElem,1:PP_nElems)
INTEGER                                               :: ii,k1,k2,jj,iEntry,iElem
!logical::test
!==================================================================================================================================

!test=.FALSE.
!DO iElem=1,PP_nElems
!  DO ii=1,nDOFElem
!    IF(ISNAN(Vin(ii,iElem)))THEN
!      !WRITE(*,*) 'NAN. ii,iElem',ii,iElem
!      TEST=.TRUE.
!    ENDIF
!  END DO
!END DO
!IF(TEST) THEN
!  print*,'buggy init'
!END IF

Vcalc=Vin
DO iElem=1,PP_nElems
  DO ii=1,nMTriangle,1
    iEntry=ii+1
    k1=IL(iElem)%iEntry(ii)
    k2=IL(iElem)%iEntry(ii+1)-1
    DO jj=k1,k2
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IL(iElem)%Entry(jj)*Vcalc(IL(iELEM)%jEntry(jj),iElem)
    END DO ! jj
  END DO ! ii
  ! backward elimination
  ! init backward Gauss
  Vout(nDOFElem,iElem) = Vcalc(nDOFElem,iElem)/DE(nDOFElem,iElem)
  DO ii=nMTriangle,1,-1
    k1=IU(iElem)%iEntry(ii)
    k2=IU(iElem)%iEntry(ii+1)-1
    iEntry=ii
    DO jj=k1,k2
      !Vin(iEntry)=Vin(iEntry)-UE(jj)*Vout(JU(jj))
      Vcalc(iEntry,iElem)=Vcalc(iEntry,iElem)-IU(iElem)%Entry(jj)*Vout(IU(iElem)%jEntry(jj),iElem)
    END DO ! jj
    Vout(iEntry,iElem)=Vcalc(iEntry,iElem)/DE(iEntry,iElem)
  END DO ! ii
END DO


!DO iElem=1,PP_nElems
!  DO ii=1,nMTriangle(iElem),2
!    k1=IL(iElem)%iEntry(ii)
!    k2=IL(iElem)%iEntry(ii+1)
!    iEntry=0.5*(ii+1)+1
!    DO jj=k1,k2
!      !Vin(iEntry)=Vin(iEntry)-LE(jj)*Vin(JL(jj))
!      Vin(iEntry,iElem)=Vin(iEntry,iElem)-IL(iElem)%Entry(jj)*Vin(IL(iELEM)%jEntry(jj),iElem)
!    END DO ! jj
!  END DO ! ii
!  ! backward elimination
!  ! init backward Gauss
!  Vout(nDOFElem,iElem) = Vin(nDOFElem,iElem)/DE(nDOFElem,iElem)
!  DO ii=nMTriangle(iElem)-1,1,-2
!    k1=IU(iElem)%iEntry(ii)
!    k2=IU(iElem)%iEntry(ii+1)
!    iEntry=0.5*(ii+1)
!    DO jj=k1,k2
!      !Vin(iEntry)=Vin(iEntry)-UE(jj)*Vout(JU(jj))
!      Vin(iEntry,iElem)=Vin(iEntry,iElem)-IU(iElem)%Entry(jj)*Vout(IU(iElem)%jEntry(jj),iElem)
!    END DO ! jj
!    Vout(iEntry,iElem)=Vin(iEntry,iElem)/DE(iEntry,iElem)
!  END DO ! ii
!END DO


!test=.FALSE.
!DO iElem=1,PP_nElems
!  DO ii=1,nDOFElem
!    IF(ISNAN(Vout(ii,iElem)))THEN
!      TEST=.TRUE.
!    ENDIF
!  END DO
!END DO
!IF(TEST) THEN
!  print*,'Vout NAN'
!END IF
END SUBROUTINE ApplyILU


SUBROUTINE ApplyBILU0_BCSR(Vin,Vout)
!==================================================================================================================================
! Application of block ILU0 preconditioner
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ILU_Vars               ,ONLY:BlockAA,BlockIA,BlockJA,nBDOF,BlockDiag
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,1:nBDOF,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,1:nBDOF,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: Vcalc(1:PP_nVar,1:nBDOF)
INTEGER                      :: ii,jj,k1,k2,iVar,iVar2,jrow,iElem
!==================================================================================================================================


DO iElem=1,PP_nElems
  Vout(:,:,iElem)=0.
  Vcalc(:,:) = Vin(:,:,iElem)
  ! lower sweep
  DO ii=2,nBDOF
    k1=BlockIA(ii)
    k2=BlockDiag(ii)-1
    !print*,'k1,k2',k1,k2
    DO jj=k1,k2
      jrow=BlockJA(jj)
      !print*,'jrwo',jrow
      DO iVar=1,PP_nVar
        DO iVar2=1,PP_nVar
          Vcalc(iVar,ii) = Vcalc(iVar,ii)-BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,jrow)
        END DO ! iVar2
      END DO ! ivar
    END DO ! jj
  END DO ! ii
  ! init
  jj=BlockDiag(nBDOF)
  DO iVar=1,PP_nVar
    DO iVar2=1,PP_nVar
      Vout(iVar,nBDOF,iElem)=Vout(iVar,nBDOF,iElem)+BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,nBDOF)
    END DO ! iVar2
  END DO ! iVar
  ! lower sweep
  DO ii=nBDOF-1,1,-1
    k1=BlockDiag(ii)+1
    k2=BlockIA(ii+1)-1
    !print*,'k1,k2',k1,k2
    DO jj=k1,k2
      jrow=BlockJA(jj)
      !print*,'jrwo',jrow
      DO iVar=1,PP_nVar
        DO iVar2=1,PP_nVar
          Vcalc(iVar,ii)=Vcalc(iVar,ii)-BlockAA(iVar,iVar2,jj,iElem)*Vout(iVar2,jrow,iElem)
        END DO ! iVar2
      END DO ! iVar
    END DO ! jj
    jj=BlockDiag(ii)
    jrow=BlockJA(jj)
    !print*,'jj',jj
    DO iVar=1,PP_nVar
      DO iVar2=1,PP_nVar
        Vout(iVar,ii,iElem)=Vout(iVar,ii,iElem)+BlockAA(iVar,iVar2,jj,iElem)*Vcalc(iVar2,jrow)
      END DO ! iVar2
    END DO ! iVar
  END DO ! ii
END DO ! iElem


END SUBROUTINE ApplyBILU0_BCSR


SUBROUTINE ApplyADI(Vin,Vout)
!===================================================================================================================================
! Apply BJ Preconditioner which is constructed per element
! which is infact an additive schwarz preconditioner
! the question is now, why is the multiplicative schwarz worse?
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars      ,ONLY:invXi,invEta,invZeta
USE MOD_LinearSolver_Vars ,ONLY:nDOFLine
USE MOD_LinearOperator    ,ONLY:DENSE_MATMUL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)              :: Vin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: Vout(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                         :: Vcalc (1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
!REAL                         :: Vcalc2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N)
REAL                         :: Vcalc(1:PP_nVar,0:PP_N)
INTEGER                      :: iElem,p,q!,i,j,k
!===================================================================================================================================

Vout(:,:,:,:,:)=0.
DO iElem=1,PP_nElems
  ! optimize
  DO q=0,PP_N
    DO p=0,PP_N
      ! xi direction
      CALL DENSE_MATMUL(nDOFLine,invXI(:,:,p,q,iElem),Vin(:,:,p,q,iElem),Vcalc)
      Vout(:,:,p,q,iElem)=Vout(:,:,p,q,iElem)+Vcalc
      ! eta direction
      CALL DENSE_MATMUL(nDOFLine,invEta(:,:,p,q,iElem),Vin(:,p,:,q,iElem),Vcalc)
      Vout(:,p,:,q,iElem)=Vout(:,p,:,q,iElem)+Vcalc
      ! zeta direction
      CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,p,q,iElem),Vin(:,p,q,:,iElem),Vcalc)
      Vout(:,p,q,:,iElem)=Vout(:,p,q,:,iElem)+Vcalc
    END DO ! p=0,PP_N
  END DO ! q=0,PP_N

! ! xi direction
!  DO j=0,PP_N
!    DO k=0,PP_N
!      CALL DENSE_MATMUL(nDOFLine,invXI(:,:,j,k,iElem),Vin(:,:,j,k,iElem),Vcalc(:,:,j,k))
!      !Vcalc2=MATMUL(invXI(:,:,j,k,iElem),Vin(:,:,j,k,iElem))
!      !Vout(:,:,j,k,iElem)=Vcalc2
!    END DO ! j
!  END DO ! k
! ! eta direction
!  DO k=0,PP_N
!    DO i=0,PP_N
!      CALL DENSE_MATMUL(nDOFLine,invEta(:,:,i,k,iElem),Vcalc(:,i,:,k),Vcalc2(:,i,:,k))
!      !Vcalc2=MATMUL(invEta(:,i,:,k,iElem),Vin(:,i,:,k,iElem))
!      !Vout(:,i,:,k,iElem)=Vout(:,i,:,k,iElem)+Vcalc2
!    END DO ! i
!  END DO ! k
!!  ! zeta direction
!  DO j=0,PP_N
!    DO i=0,PP_N
!      !CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vout(:,i,j,:,iElem),Vout(:,i,j,:,iElem))
!      CALL DENSE_MATMUL(nDOFLine,invZeta(:,:,i,j,iElem),Vcalc2(:,i,j,:),Vout(:,i,j,:,iElem))
!      !Vcalc2=MATMUL(invZeta(:,i,j,:,iElem),Vin(:,i,j,:,iElem))
!      !Vout(:,i,j,:,iElem)=Vout(:,i,j,:,iElem)+Vcalc2
!    END DO ! i
!  END DO ! j
!!  ! missing?
!!  !Vout=1./3.*Vout -- not working

END DO ! iElem
!Vout=1./3.*Vout

END SUBROUTINE  ApplyADI

END MODULE MOD_ApplyPreconditioner
