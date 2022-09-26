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

MODULE MOD_ILU
!===================================================================================================================================
! LU SGS preconditioner
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
INTERFACE InitILU
  MODULE PROCEDURE InitILU
END INTERFACE

INTERFACE FinalizeILU
  MODULE PROCEDURE FinalizeILU
END INTERFACE

INTERFACE BuildBILU0BCSR
  MODULE PROCEDURE BuildBILU0BCSR
END INTERFACE

PUBLIC::InitILU,FinalizeILU,BuildBILU0BCSR
!===================================================================================================================================

CONTAINS

SUBROUTINE InitILU()
!===================================================================================================================================
! Init LUSGS
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_ILU_Vars
!USE MOD_ReadInTools,         ONLY:GETINT
USE MOD_LinearSolver_Vars,   ONLY:nDOFElem,nDOFLine
USE MOD_Precond_Vars,        ONLY:BlockSize,PrecondType,nBlockSize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iHelp,i,j,k,mm,nn,oo,r,s
INTEGER              :: delta(0:PP_N,0:PP_N)
LOGICAL              :: first
INTEGER,ALLOCATABLE  :: TestMat(:,:)
!CHARACTER(LEN=17)    :: strfmt
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' INIT ILU                                ...'

!REAL,ALLOCATABLE,DIMENSION(:,:,:)       :: DiagBILU0
!REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:)   :: XiBILU0,EtaBILU0,ZetaBILU0
BlockSize  = PP_nVar
nBlockSize = INT(nDOFElem/BlockSize)
IF(PrecondType.EQ.25)THEN
  ALLOCATE( DiagBILU0(PP_nVar,PP_nVar,nBlockSize,PP_nElems)   &
          , XiBILU0  (nDOFLine,nDOFLine,0:PP_N,0:PP_N,PP_nElems)  &
          , EtaBILU0 (nDOFLine,nDOFLine,0:PP_N,0:PP_N,PP_nElems)  &
          , ZetaBILU0(nDOFLine,nDOFLine,0:PP_N,0:PP_N,PP_nElems)  )
  DiagBILU0=0.
  XiBILU0=0.
  EtaBILU0=0.
  ZetaBILU0=0.
ELSE ! use block csr format
  nBlockEntries=3*(PP_N+1)**4-2*(PP_N+1)**3
  nBDOF=(PP_N+1)**3
  print*,'nBlockEntries',nBlockEntries
  print*,'nBDOF',nBDOF
  ALLOCATE( BlockAA(PP_nVar,PP_nVar,nBlockEntries,PP_nElems) &
          , BlockJA(1:nBlockEntries)                         &
          , BlockIA(1:nBDOF+1)                               &
          , BlockDiag(1:nBDOF)                               )
  ! precompute structure
  ALLOCATE(TestMat(1:nBDOF,1:nBDOF))
  delta=0
  DO i=0,PP_N
    delta(i,i)=1
  END DO
  TestMat=0
  s=1
  DO oo=0,PP_N
    DO nn=0,PP_N
      DO mm=0,PP_N
        r=1
        DO k=0,PP_N
          DO j=0,PP_N
            DO i=0,PP_N
              IF(delta(j,nn).EQ.1.AND.delta(k,oo).EQ.1)THEN
                TestMat(r,s)=1
              END IF
              IF(delta(i,mm).EQ.1.AND.delta(k,oo).EQ.1)THEN
                TestMat(r,s)=1
              END IF
              IF(delta(i,mm).EQ.1.AND.delta(j,nn).EQ.1)THEN
                TestMat(r,s)=1
              END IF
              r=r+1
            END DO !i
          END DO !j
        END DO !k
        s=s+1
      END DO ! m
    END DO ! n
  END DO ! o
  ! mapping for csr format
  ihelp=0
  DO r=1,nBDOF
    first=.TRUE.
    DO s=1,nBDOF
      IF(TestMat(r,s).EQ.1)THEN
        ihelp=ihelp+1
        !BlockEntry(:,:,ihelp,iElem)=Test(r,s)
        BlockJA(ihelp)=s
        IF(first)THEN
          BlockIA(r)=ihelp
          first=.FALSE.
        END IF ! first
      END IF ! zero
      IF(r.EQ.s)THEN
        BlockDiag(r)=ihelp
      END IF
    END DO ! s
  END DO ! r
  BlockIA(nBDOF+1)=BlockIA(1)+nBlockEntries

!  ! modification for zero line or single value lines
!  IF(MPIRoot)THEN
!    WRITE(strfmt,'(A1,I4,A12)')'(',nBDOF,'(1X,I2))'
!    WRITE(UNIT_stdOut,*)'Debug Block Jacobian to:'
!    OPEN (UNIT=103,FILE='MPIRank0_DOF_Sparsity.dat',STATUS='REPLACE')
!    DO r=1,nBDOF
!      WRITE(103,strfmt) TestMat(r,:)
!    END DO
!    CLOSE(103)
!  END IF

  DEALLOCATE(TestMat)

END IF

SWRITE(UNIT_stdOut,'(A)')' FINISHED init ILU!'

END SUBROUTINE InitILU

SUBROUTINE FinalizeILU()
!===================================================================================================================================
! Finalize LUSGS
!===================================================================================================================================
! MODULES
USE MOD_ILU_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================


END SUBROUTINE FinalizeILU


SUBROUTINE BuildBILU0BCSR(BJ,iElem)
!===================================================================================================================================
! Construct LUSGS
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mathtools         ,ONLY: INVERSE
USE MOD_LinearSolver_Vars ,ONLY: nDOFELEM
USE MOD_ILU_Vars          ,ONLY: BlockAA,BlockIA,BlockJA,nBDOF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: BJ(nDOFELem,nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: s,r,t
REAL               :: epsZero
REAL               :: dummy(PP_nVar,PP_nVar)
INTEGER            :: i,j,k1,k2,vn1,vn2
!===================================================================================================================================

epsZero=EPSILON(0.0d0)*PP_nVar*PP_nVar

! Saad
DO s=PP_nVar,nDOFElem-PP_nVar,PP_nVar
  DO r=0,s-PP_nVar,PP_nVar ! richtig? +1 fehlend?
    IF(SUM(ABS(BJ(s+1:s+PP_nVar,r+1:r+PP_nVar))).GT.epsZero)THEN
      dummy=INVERSE(BJ(r+1:r+PP_nVar,r+1:r+PP_nVar))
      BJ(s+1:s+PP_nVar,r+1:r+PP_nVar)=MATMUL(BJ(s+1:s+PP_nVar,r+1:r+PP_nVar),dummy)
      DO t=r+PP_nVar,nDOFElem-PP_nVar,PP_nVar
        IF(SUM(ABS(BJ(s+1:s+PP_nVar,t+1:t+PP_nVar))).GT.epsZero)THEN
          BJ(s+1:s+PP_nVar,t+1:t+PP_nVar) = BJ(s+1:s+PP_nVar,t+1:t+PP_nVar) &
                                          - MATMUL(BJ(s+1:s+PP_nVar,r+1:r+PP_nVar),BJ(r+1:r+PP_nVar,t+1:t+PP_nVar))
        END IF
      END DO ! t
    END IF !
  END DO ! r
END DO ! s

! store in csr format
DO i=1,nBDOF
  k1=BlockIA(i)
  k2=BlockIA(i+1)-1
  !print*,'k1,k2',k1,k2
  DO j=k1,k2
    vn1=(i-1)*PP_nVar
    vn2=(BlockJA(j)-1)*PP_nVar
    !print*,'i,j,vn1,vn2',i,j,vn1+1,vn1+PP_nVar,vn2+1,vn2+PP_nVar
    BlockAA(1:PP_nVar,1:PP_nVar,j,iElem)=BJ(vn1+1:vn1+PP_nVar,vn2+1:vn2+PP_nVar)
    IF(BlockJA(j).EQ.i)THEN
      !print*,'diag',BlockJA(j),i
      BlockAA(1:PP_nVar,1:PP_nVar,j,iElem)=INVERSE(BJ(vn1+1:vn1+PP_nVar,vn2+1:vn2+PP_nVar))
    END IF
  END DO ! j
END DO ! i

END SUBROUTINE BuildBILU0BCSR

END MODULE MOD_ILU
