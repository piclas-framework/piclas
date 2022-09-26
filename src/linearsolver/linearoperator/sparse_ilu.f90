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

MODULE MOD_SparseILU
!===================================================================================================================================
! Module containing matrix vector operations
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
INTERFACE InitSparseILU
  MODULE PROCEDURE InitSparseILU
END INTERFACE

INTERFACE FinalizeSparseILU
  MODULE PROCEDURE FinalizeSparseILU
END INTERFACE

INTERFACE BuildILU0
  MODULE PROCEDURE BuildILU0
END INTERFACE

PUBLIC:: InitSparseILU,FinalizeSparseILU,BuildILU0
!===================================================================================================================================

CONTAINS

SUBROUTINE InitSparseILU()
!===================================================================================================================================
! Init Sparse ILU
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY:nDOFelem
USE MOD_CSR_Vars               ,ONLY:nUNonZeros,nLNonZeros
USE MOD_CSR_Vars               ,ONLY:DE,IL,IU
USE MOD_Precond_Vars           ,ONLY:PrecondType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

IF(PrecondType.EQ.3)THEN
  ALLOCATE(nUNonZeros(PP_nElems) &
          ,nLNonZeros(PP_nElems) )

  ALLOCATE(DE(nDOFElem,PP_nElems)&
          ,IU(PP_nElems)         &
          ,IL(PP_nElems)         )
END IF

END SUBROUTINE InitSparseILU

SUBROUTINE FinalizeSparseILU()
!===================================================================================================================================
! Init Sparse ILU
!===================================================================================================================================
! MODULES
USE MOD_CSR_Vars   ,ONLY: nUNonZeros,nLNonZeros
USE MOD_CSR_Vars   ,ONLY: DE,IU,IL
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DEALLOCATE( nUNonZeros, nLNonZeros)
DEALLOCATE( DE, IL,IU )
END SUBROUTINE FinalizeSparseILU

SUBROUTINE BuildILU0(ILU0,iElem)
!===================================================================================================================================
! Build the ILU0 per Block
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
!USE MOD_Mathtools,ONLY:INVERSE
USE MOD_LinearSolver_Vars          ,ONLY:nDOFelem
USE MOD_CSR_Vars                   ,ONLY:nUNonZeros,nLNonZeros,nMTriangle
USE MOD_CSR_VArs                   ,ONLY:DE,IL,IU
!USE MOD_CSR,                    ONLY:WriteMat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL                        :: ILU0(nDOFelem,nDOFElem)
INTEGER,INTENT(IN)          :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: ii,kk,jj,iEntry
LOGICAL                     :: first
REAL                        :: Sparsity
INTEGER                     :: lastLine,lineNonZero
LOGICAL                     :: singleValue
REAL                        :: epsZero
! check
REAL,DIMENSION(nDofElem,nDofElem)       :: PLOC!,Ainv
!===================================================================================================================================

!epsZero=0.5*EPSILON(0.0)
epsZero=EPSILON(0.0d0)


PLOC=ILU0

!----------------------------------------------------------------------------------------------------------------------------------
! ILU(0) // drop all zero values // Saad 'Iterative Methods for sparse linear systems' p.307 // p.327 in pdf
!----------------------------------------------------------------------------------------------------------------------------------


!filename='first.dat'
!OPEN (UNIT=103,FILE=filename,STATUS='UNKNOWN')
!WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
!DO ii=1,nDOFelem
!  WRITE(103,strfmt) ILU0(ii,:)
!END DO
!CLOSE(103)

!! null
!DO ii=1,nDOFElem
!  DO kk=1,nDOFElem
!    IF(ABS(ILU0(ii,kk)).LT.1e-15)THEN
!      ILU0(ii,kk) = 0.
!    END IF
!  END DO ! kk
!END DO ! ii

! Saad
DO ii=2,nDOFElem
  DO kk=1,ii-1
    !IF(ILU0(ii,kk).NE.0.)THEN
    IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
      ILU0(ii,kk) = ILU0(ii,kk)/ILU0(kk,kk)
      DO jj=kk+1,nDOFElem
        !IF(ILU0(ii,jj).NE.0.)THEN
        IF(ABS(ILU0(ii,jj)).GT.epsZero)THEN
          ILU0(ii,jj) = ILU0(ii,jj) - ILU0(ii,kk)*ILU0(kk,jj)
        END IF ! ii,jj element NZ
      END DO ! jj
    END IF ! ii,kk element NZ
  END DO ! ii
END DO ! kk

!CALL WriteMat(nDOFElem,ILU0)

!! implementation after meister, p. 217
!L=0.
!U=0.
!DO ii=1,nDOFElem ! column
!  ! row for lower
!  DO kk=ii, nDOFElem
!    IF(ABS(ILU0(kk,ii)).GT.epsZero)THEN
!      temp=0.d0
!      DO jj=1,ii-1
!        IF(ABS(ILU0(kk,jj).GT.epsZero).AND.ABS(ILU0(jj,ii).GT.epsZero))THEN
!          temp = temp - L(kk,jj) * U(jj,ii)
!        END IF
!      END DO ! jj
!      L(kk,ii) = ILU0(kk,ii) + temp
!    END IF
!  END DO ! kk
!  ! column for upper
!  DO kk=ii+1,nDOFElem
!    IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
!      temp=0.d0
!      DO jj=1,ii-1
!        IF(ABS(ILU0(ii,jj).GT.epsZero).AND.ABS(ILU0(jj,kk).GT.epsZero))THEN
!          temp = temp - L(ii,jj) * U(jj,kk)
!        END IF
!      END DO ! jj
!      U(ii,kk) = ( ILU0(ii,kk) + temp )/ L(ii,ii)
!    END IF
!  END DO ! kk
!END DO ! ii
!
!! and map all to the old matrix
!!ILU0=L+U
!DO ii=1,nDOFElem
!  DO kk=ii,nDOFElem
!    ILU0(kk,ii)=L(kk,ii)
!  END DO
!  DO kk=ii+1,nDOFElem
!    ILU0(ii,kk) = U(ii,kk)
!  END DO
!END DO

!! LU deposition
!DO kk=1,nDOFElem-1
 ! DO ii=kk+1,nDOFElem
 !   ILU0(ii,kk) = ILU0(ii,kk)/ILU0(kk,kk)
 !     DO jj=kk+1,nDOFElem
 !       ILU0(ii,jj) = ILU0(ii,jj) - ILU0(ii,kk)*ILU0(kk,jj)
 !     END DO ! jj
 ! END DO ! ii
!END DO !kk
!LU=ILU0

! nan checker
DO kk=1,nDOFElem
  DO ii=1,nDOFElem
    IF(ILU0(ii,kk).NE.ILU0(ii,kk))THEN
      CALL abort(&
      __STAMP__&
      ,' NAN! kk,ii',kk,REAL(ii))
    END IF
  END DO
END DO

!----------------------------------------------------------------------------------------------------------------------------------
! write complete LU to file
!----------------------------------------------------------------------------------------------------------------------------------
!
!! nullify
!U=0.
!L=0.
!
!! LU
!DO ii=1,nDOFElem
!  DO kk=1,nDOFElem
!    IF(kk.GE.ii)THEN
!      U(ii,kk)=ILU0(ii,KK)
!    END IF
!    IF(kk.LE.ii)THEN
!      L(ii,kk)=ILU0(ii,KK)
!    END IF
!  END DO ! kk
!END DO ! ii
!
!
!filename='LU_U.dat'
!OPEN (UNIT=103,FILE=filename,STATUS='UNKNOWN')
!WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
!DO ii=1,nDOFelem
!  WRITE(103,strfmt) U(ii,:)
!END DO
!CLOSE(103)
!filename='LU_L.dat'
!OPEN (UNIT=103,FILE=filename,STATUS='UNKNOWN')
!WRITE(strfmt,'(A1,I4,A12)')'(',nDOFelem,'(1X,E23.16))'
!DO ii=1,nDOFelem
!  WRITE(103,strfmt) L(ii,:)
!END DO
!CLOSE(103)



!----------------------------------------------------------------------------------------------------------------------------------
! Extendet CSR FORMAT
!----------------------------------------------------------------------------------------------------------------------------------
! null
!DO ii=1,nDOFElem
!  DO kk=1,nDOFElem
!    IF(ABS(ILU0(ii,kk)).LT.1e-15)THEN
!      ILU0(ii,kk) = 0.
!    END IF
!  END DO ! kk
!END DO ! ii

! get number of non-zero entries
!print*,iElem
nUNonZeros(iElem)=0
nLNonZeros(iElem)=0

DO ii=1,nDOFElem
  DO kk=1,nDOFElem
    ! upper
    IF(kk.GT.ii)THEN
      !IF(ILU0(ii,kk).NE.0.)THEN
      IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
        nUNonZeros(iElem)=nUNonZeros(iElem)+1
      END IF
    END IF ! upper
    ! lower
    IF(kk.LT.ii)THEN
      !IF(ILU0(ii,kk).NE.0.)THEN
      IF(ABS(ILU0(ii,kk)).GT.epsZERO)THEN
        nLNonZeros(iElem)=nLNonZeros(iElem)+1
      END IF
    END IF ! lower
  END DO ! kk
END DO ! ii

Sparsity= REAL(nUNonZeros(iElem))+REAL(nLNonZeros(iElem)) + REAL(nDOFElem)
Sparsity=Sparsity/REAL(nDOFElem)/REAL(nDOFElem)
!WRITE(*,'(A,I6)')       ' Element number          : ', iElem
!WRITE(*,'(A,I6)')       ' Number of non zeros in U: ', nUNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in L: ', nLNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in D: ', nDOFElem
!WRITE(*,'(A,F12.7,A)')  ' Level of Fill             ', sparsity*100.0,'%'
!print*

! Allocate
!nMTriangle(iElem)=2*(nDOFelem-1)
! simple version
nMTriangle=nDOFelem-1
SDEALLOCATE(IU(iElem)%Entry)
SDEALLOCATE(IU(iELEM)%JEntry)
SDEALLOCATE(IU(iELEM)%IEntry)
SDEALLOCATE(IL(iELEM)%Entry)
SDEALLOCATE(IL(iELEM)%JEntry)
SDEALLOCATE(IL(iELEM)%IEntry)
ALLOCATE( IU(iElem)%Entry(nUNonZeros(iElem))  &
        , IU(iELEM)%JEntry(nUNonZeros(iElem)) &
        , IU(iELEM)%IEntry(nDOFElem) &
        , IL(iELEM)%Entry(nLNonZeros(iElem))  &
        , IL(iELEM)%JEntry(nLNonZeros(iElem)) &
        , IL(iELEM)%IEntry(nDOFElem) )

! nullify
DE(:,iElem)=0.
IL(iELem)%Entry=0.
IU(iELem)%Entry=0.

! U part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nMTriangle,1
  first=.TRUE.
  lineNonZero=0
  DO kk=1,nDOFElem
    IF(kk.GT.ii)THEN
      !IF((ILU0(ii,kk)).NE.0.)THEN
      IF(ABS(ILU0(ii,kk)).GT.epsZero)THEN
        jj=jj+1
        IU(iElem)%Entry(jj)=ILU0(ii,kk)
        IU(iELem)%JENTRY(jj)=kk
        lineNonZero=lineNonZero+1
        IF(first)THEN
          IU(iElem)%IENTRY(ii)=jj
          first=.FALSE.
        END IF !first
      END IF ! zero
    END IF ! upper part
  END DO ! kk
  ! modification for zero line or single value lines
  IF(lineNonZero.EQ.0)THEN
    IF(ii.EQ.1)THEN
      IU(iElem)%iEntry(ii)=1
    ELSE
      IU(iElem)%iEntry(ii)=IU(iElem)%iEntry(ii-1)+Lastline
    END IF
  END IF
  ! last line
  LastLine=LineNonZero
END DO ! ii
 IU(iElem)%iEntry(nDOFElem)= IU(iElem)%iEntry(1)+nUNonZeros(iElem)

!jj=0
!DO ii=1,nMTriangle(iElem),2
!  first=.TRUE.
!  last=.TRUE.
!  iEntry=0.5*(ii+1)
!  DO kk=1,nDOFElem
!    IF(kk.GT.iEntry)THEN
!      IF(ILU0(iEntry,kk).NE.0)THEN
!        jj=jj+1
!        IU(iElem)%Entry(jj)=ILU0(iEntry,kk)
!        IU(iELem)%JENTRY(jj)=kk
!        IF(first)THEN
!          IU(iElem)%IENTRY(ii)=jj
!          first=.FALSE.
!        ELSE
!          IU(iELEM)%IEntry(ii+1)=jj
!          last=.FALSE.
!        END IF
!      END IF
!    END IF ! upper
!  END DO ! kk
!  IF(first)THEN
!    IU(iELEM)%iEntry(ii)  =nDOFElem
!    IU(iElem)%iEntry(ii+1)=1
!  ELSE ! only one entry .NE. zero
!    IF(last)THEN
!      IU(iElem)%iEntry(ii+1)=IU(iElem)%iEntry(ii)
!    END IF ! last
!  END IF ! first
!END DO ! ii

! L part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nMTriangle,1
  first=.TRUE.
  lineNonZero=0
  iEntry=ii+1
  DO kk=1,nDOFElem
    IF(kk.LT.iEntry)THEN
      !IF(ILU0(iEntry,kk).NE.0)THEN
      IF(ABS(ILU0(iEntry,kk)).GT.epsZero)THEN
        jj=jj+1
        IL(iElem)%Entry(jj) =ILU0(iEntry,kk)
        IL(iElem)%jEntry(jj)=kk
        lineNonZero=lineNonZero+1
        IF(first)THEN
          IL(iElem)%iEntry(ii)=jj
          first=.FALSE.
        END IF ! first
      END IF ! zero
    ELSE ! kk GE iEntry
      CYCLE
    END IF ! lower
  END DO ! kk
  ! modification for zero line or single value lines
  IF(lineNonZero.EQ.0)THEN
    IF(ii.EQ.1)THEN
      IL(iElem)%iEntry(ii)=1
    ELSE
      IL(iElem)%iEntry(ii)=IL(iElem)%iEntry(ii-1)+Lastline
    END IF
  END IF
  ! last line
  LastLine=LineNonZero
END DO ! ii
IL(iElem)%iEntry(nDOFelem)=IL(iElem)%iEntry(1)+nLNonZeros(iElem)


!jj=0
!DO ii=1,nMTriangle(iElem),2
!  first=.TRUE.
!  last=.TRUE.
!  iEntry=0.5*(ii+1)+1
!  DO kk=1,nDOFElem
!    IF(kk.LT.iEntry)THEN
!      IF(ILU0(iEntry,kk).NE.0)THEN
!        jj=jj+1
!        IL(iElem)%Entry(jj) =ILU0(iEntry,kk)
!        IL(iElem)%jEntry(jj)=kk
!        IF(first)THEN
!          IL(iElem)%iEntry(ii)=jj
!          first=.FALSE.
!        END IF
!        IF(.NOT.(first))THEN
!          IL(iElem)%iEntry(ii+1)=jj
!          last=.FALSE.
!        END IF
!      END IF
!    END IF ! upper
!  END DO ! kk
!  IF(first)THEN
!    IL(iElem)%iEntry(ii) =nDOFElem
!    IL(iElem)%iEntry(ii+1)=1
!  ELSE ! only one entry .NE. zero
!    IF(last)THEN
!      IU(iElem)%IEntry(ii+1)=IU(iElem)%IEntry(ii)
!    END IF ! last
!  END IF ! first
!END DO ! ii

! Dioganal entires
DO ii=1,nDOFElem
  DE(ii,iElem)=ILU0(ii,ii)
!  IF(ILU0(ii,ii).EQ.0)THEN
!    WRITE(*,*) "Zero element on diagonal at ",ii
!    STOP
!  END IF
!  IF(ISNAN(ILU0(ii,ii)))THEN
!    WRITE(*,*) 'Diagonal element of ',ii,' is NAN'
!    STOP
!  END IF
END DO

!----------------------------------------------------------------------------------------------------------------------------------
! verification
!----------------------------------------------------------------------------------------------------------------------------------
!
!ALLOCATE( Vin   (1:nDOFElem)&
!         ,Vresu1(1:nDOFElem)&
!         ,Vresu2(1:nDOFElem)&
!         ,B     (1:nDOFElem))
!! null
!Vin=0.
!Vresu1=0.
!Vresu2=0.
!
!! fill with bin
!DO ii=1,nDOFElem
!  !Vin(ii) =REAL(ii)
!  Vin(ii)= 4./97*REAL(ii)**2 -3*REAL(ii)
!END DO
!vin=2.34
!B=Vin
!vtild=vin
!
!DO s = 2, nDOFElem
!  DO r= 1, s-1
!    vin(s) = vin(s) - LU(s,r)*vin(r)
!  END DO ! r
!END DO ! s
!DO s=nDOFElem,1,-1
! DO r = s+1, nDOFElem
!   vin(s) = vin(s) - LU(s,r) * vresu2(r)
! END DO ! r
! vresu2(s) = vin(s) / LU(s,s)
! END DO ! s
!
!
!!! Verification LU
!!! Gauss Elimination / LU
!!DO kk=2,nDOFElem
!!  DO ii=1, kk-1
!!    Vin(kk) = Vin(kk) - ILU0(kk,ii) * Vin(ii)
!!  END DO ! ii
!!END DO ! kk
!!!Vresu2=Vin
!!!Vresu2(nDOFElem)=Vin(nDOFElem)/ILU0(nDOFElem,nDOFElem)
!!!backwards
!!DO kk=nDOFElem,1,-1
!!  DO ii=kk+1,nDOFElem
!!    Vin(kk) = Vin(kk) - ILU0(kk,ii)*Vresu2(ii)
!!  END DO
!!  Vresu2(kk) = Vin(kk) / ILU0(kk,kk)
!!END DO
!
!
!! sparse
!! new test as intent in
!DO ii=1,nMTriangle,1
!  iEntry=ii+1
!  k1=IL(iElem)%iEntry(ii)
!  k2=IL(iElem)%iEntry(ii+1)-1
!  DO jj=k1,k2
!    B(iEntry)=B(iEntry)-IL(iElem)%Entry(jj)*B(IL(iELEM)%jEntry(jj))
!  END DO ! jj
!END DO ! ii
!!DO ii=1,nMTriangle,1
!  !iEntry=ii+1
!  !k1=IL(iElem)%iEntry(ii)
!  !k2=IL(iElem)%iEntry(ii+1)-1
!  !DO jj=k1,k2
!    !B(iEntry)=B(iEntry)-IL(iElem)%Entry(jj)*B(IL(iELEM)%jEntry(jj))
!  !END DO ! jj
!!END DO ! ii
!!Vresu1=b
!! backward elimination
!! init backward Gauss
!Vresu1(nDOFElem) = B(nDOFElem)/DE(nDOFElem,iElem)
!DO ii=nMTriangle,1,-1
!  k1=IU(iElem)%iEntry(ii)
!  k2=IU(iElem)%iEntry(ii+1)-1
!  iEntry=ii
!  DO jj=k1,k2
!    !Vin(iEntry)=Vin(iEntry)-UE(jj)*Vout(JU(jj))
!    B(iEntry)=B(iEntry)-IU(iElem)%Entry(jj)*Vresu1(IU(iElem)%jEntry(jj))
!  END DO ! jj
!  Vresu1(iEntry)=B(iEntry)/DE(iEntry,iElem)
!END DO ! ii
!
!WRITE(*,'(A,E23.14)') ' Difference of both matrix vector products is', SUM(ABS(Vresu2)-ABS(Vresu1))
!
!
!Ainv(:,:)=INVERSE(Ploc)
!vResu3=MATMUL(Ainv,vtild)
!
!kk=0
!DO ii = 1, nDOFElem
!  IF(ABS(vresu3(ii)-vresu1(ii)).GT.1e-6)THEN
!    kk=kk+1
!  END IF
!END DO
!WRITE(*,*) 'Inverse M - Vsparse: ', kk
!kk=0
!DO ii = 1, nDOFElem
!  IF(ABS(vresu3(ii)-vresu2(ii)).GT.1e-6)THEN
!    kk=kk+1
!  END IF
!END DO
!WRITE(*,*) 'Inverse M - V_Gauss elimination: ', kk
!kk=0
!DO ii = 1, nDOFElem
!  IF(ABS(vresu1(ii)-vresu2(ii)).GT.1e-6)THEN
!    kk=kk+1
!  END IF
!END DO
!WRITE(*,*) 'difference between Vsparse - V_Gauss ', kk
!DO ii=1,nDOFElem
!  IF(ISNAN(Vresu1(ii)))THEN
!    WRITE(*,*) 'NAN. ii',ii
!  ENDIF
!END DO

END SUBROUTINE BuildILU0


END MODULE MOD_SparseILU
