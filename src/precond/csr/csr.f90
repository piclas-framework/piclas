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

MODULE MOD_CSR
!===================================================================================================================================
! Contains the initialization of the DG global variables
! Computes the different DG spatial operators/residuals(Ut) using U
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
INTERFACE InitCSR
  MODULE PROCEDURE InitCSR
END INTERFACE

!INTERFACE DGTimeDerivative_weakFormCSR
!  MODULE PROCEDURE DGTimeDerivative_weakFormCSR
!END INTERFACE

INTERFACE FinalizeCSR
  MODULE PROCEDURE FinalizeCSR
END INTERFACE

INTERFACE CSR
  MODULE PROCEDURE CSR
END INTERFACE

INTERFACE DiagCSR
  MODULE PROCEDURE DiagCSR
END INTERFACE

INTERFACE GlobalCSR
  MODULE PROCEDURE GlobalCSR
END INTERFACE

INTERFACE GlobalBCSR
  MODULE PROCEDURE GlobalBCSR
END INTERFACE

INTERFACE WriteMat
  MODULE PROCEDURE WriteMat
END INTERFACE

!PUBLIC::InitCSR,DGTimeDerivative_weakFormCSR,FinalizeCSR,CSR,DiagCSR
PUBLIC::InitCSR,FinalizeCSR,CSR,DiagCSR,GlobalCSR, GetNonZeros, WriteMat, GlobalBCSR
!===================================================================================================================================

CONTAINS

SUBROUTINE InitCSR()
!===================================================================================================================================
! Allocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CSR_Vars
!USE MOD_Mesh_Vars,          ONLY: nBCSides
USE MOD_LinearSolver_Vars,  ONLY: nDOFelem
USE MOD_Jac_ex,             ONLY: Jac_ex, Jac_Ex_Neighbor, InitJac_Ex
USE MOD_DG_Vars,            ONLY: L_HatMinus, L_HatPlus
!USE MOD_Mesh_Vars,          ONLY: ElemToSide
!USE MOD_Precond,            ONLY: BuildnVecSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: locSize
!INTEGER  :: iElem,r,SideID,ilocSide
!REAL,DIMENSION(:,:),ALLOCATABLE   :: GlobalJacobian
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT CSR...'

locSize=PP_nVar*(PP_N+1)**3
nDOFElem= locSize
! allocate the sparse stuff
ALLOCATE( nNonZeros(0:6,1:PP_nElems)            , &
          !DiagonalEntries(1:locSize,1:PP_nElems), &
          SparseMatrix(0:6,1:PP_nElems)           )


ALLOCATE(L_HatPlusMinus(0:PP_N,1:6))

L_HatPlusMinus(:,  XI_MINUS) = L_HatMinus
L_HatPlusMinus(:, ETA_MINUS) = L_HatMinus
L_HatPlusMinus(:,ZETA_MINUS) = L_HatMinus
L_HatPlusMinus(:,  XI_PLUS)  = L_HatPlus
L_HatPlusMinus(:, ETA_PLUS)  = L_HatPlus
L_HatPlusMinus(:,ZETA_PLUS)  = L_HatPlus

! machine accuracy
epsZero=EPSILON(0.0d0)

! initialize local normal vectors
!CALL BuildnVecSurf()

!CALL InitJac_Ex()

! performed in blabla
!!ALLOCATE( DebugMat(1:nDOFElem*PP_nElems,1:nDOFElem*PP_nElems)) ! zero is own element and 1:6 are the neighbor elements
!#if (PP_TimeDiscMethod==2)
!
!ALLOCATE( Ucsr (1:nDOFElem,1:PP_nElems) &
!        , Utcsr(1:nDOFElem,1:PP_nElems) )
!
!ALLOCATE( GlobalJacobian(1:nDOFElem,1:nDOFElem)) ! zero is own element and 1:6 are the neighbor elements
!!DebugMat=0.
!DO iElem=1,PP_nElems
!  ! diagonal entries
!  GlobalJacobian=0.
!  CALL Jac_ex(iElem,GlobalJacobian(:,:))
!  !DO r=1,nDOFElem
!  !  DiagonalEntries(r,iElem) = GlobalJacobian(r,r)
!  !END DO ! r
!  CALL CSR(GlobalJacobian(:,:),nDOFElem,0,iElem)
!  ! loop over all sides == Neighbours
!  DO ilocSide=1,6
!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    IF(SideID.LE.nBCSides) CYCLE
!    !SideToElem(S2E_NB_ELEM_ID,SideID)))
!    CALL Jac_Ex_Neighbor(ilocSide,iElem,GlobalJacobian(:,:))
!    CALL CSR(GlobalJacobian(:,:),nDOFElem,ilocSide,iElem)
!  END DO! ilocSide
!END DO ! iElem
!SDEALLOCATE( GlobalJacobian )
!#endif /* TimeDiscMethod==2 */

SWRITE(UNIT_stdOut,'(A)')' INIT CSR DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitCSR



SUBROUTINE CSR(Mat,nDOFElem,locSideID,ElemID)
!===================================================================================================================================
! Init Sparse ILU
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_CSR_Vars                ,ONLY: SparseMatrix,nNonZeros,epsZero
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: locSideID,ElemID,nDOFElem
REAL,INTENT(IN)                 :: Mat(nDOFElem,nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ii,kk,jj
INTEGER                         :: lastLine,lineNonZero
LOGICAL                         :: singleValue,first
!===================================================================================================================================

nNonZeros(locSideID,ElemID)=0

DO ii=1,nDOFElem
  DO kk=1,nDOFElem
    IF(ABS(Mat(ii,kk)).GT.epsZero)THEN
      nNonZeros(locSideID,ElemID)=nNonZeros(locSideID,ElemID)+1
    END IF
  END DO ! kk
END DO ! ii

!Sparsity= REAL(nUNonZeros(iElem))+REAL(nLNonZeros(iElem)) + REAL(nDOFElem)
!Sparsity=Sparsity/REAL(nDOFElem)/REAL(nDOFElem)
!WRITE(*,'(A,I6)')       ' Element number          : ', iElem
!WRITE(*,'(A,I6)')       ' Number of non zeros in U: ', nUNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in L: ', nLNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in D: ', nDOFElem
!WRITE(*,'(A,F12.7,A)')  ' Level of Fill             ', sparsity*100.0,'%'
!print*

ALLOCATE( SparseMatrix(locSideID,ElemID)%Entry(nNonZeros(locSideID,ElemID)) &
        , SparseMatrix(locSideID,ElemID)%JEntry(nNonZeros(locSideID,ElemID))  &
        , SparseMatrix(locSideID,ElemID)%IEntry(nDOFElem+1)                   )

! nullify
SparseMatrix(locSideID,ElemID)%Entry =0.
SparseMatrix(locSideID,ElemID)%JEntry  =0
SparseMatrix(locSideID,ElemID)%IEntry  =0


! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nDOFElem,1
  first=.TRUE.
  lineNonZero=0
  DO kk=1,nDOFElem
    IF(ABS(Mat(ii,kk)).GT.epsZero)THEN
      jj=jj+1
      SparseMatrix(locSideID,ElemID)%Entry(jj) =Mat(ii,kk)
      SparseMatrix(locSideID,ElemID)%JENTRY(jj)=kk
      lineNonZero=lineNonZero+1
      IF(first)THEN
        SparseMatrix(locSideID,ElemID)%IENTRY(ii)=jj
        first=.FALSE.
      END IF !first
    END IF ! zero
  END DO ! kk
  ! modification for zero line or single value lines
  IF(lineNonZero.EQ.0)THEN
    IF(ii.EQ.1)THEN
      SparseMatrix(locSideID,ElemID)%iEntry(ii)=1
    ELSE
      SparseMatrix(locSideID,ElemID)%iEntry(ii)=SparseMatrix(locSideID,ElemID)%iEntry(ii-1)+Lastline
    END IF
  END IF
  ! last line
  LastLine=LineNonZero
END DO ! ii
SparseMatrix(locSideID,ElemID)%iEntry(nDOFElem+1)= SparseMatrix(locSideID,ElemID)%iEntry(1)+nNonZeros(locSideID,ElemID)

END SUBROUTINE CSR

SUBROUTINE DiagCSR(ILU0,iElem)
!===================================================================================================================================
! Build the ILU0 per Block
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_LinearSolver_Vars      ,ONLY:nDOFelem
USE MOD_CSR_Vars               ,ONLY:nUNonZeros,nLNonZeros,nMTriangle
USE MOD_CSR_Vars               ,ONLY:DE,IL,IU
USE MOD_CSR_Vars               ,ONLY:epsZero
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
!INTEGER                     :: k1,k2
!REAL,ALLOCATABLE            :: Vin(:)
!REAL                        :: Sparsity
INTEGER                     :: lastLine,lineNonZero
LOGICAL                     :: singleValue
!REAL                        :: L(nDOFELEM,nDOFElem)
!REAL                        :: U(nDOFELEM,nDOFElem)
!===================================================================================================================================

!epsZero=0.5*EPSILON(0.0)

!----------------------------------------------------------------------------------------------------------------------------------
! Extendet CSR FORMAT
!----------------------------------------------------------------------------------------------------------------------------------

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

!Sparsity= REAL(nUNonZeros(iElem))+REAL(nLNonZeros(iElem)) + REAL(nDOFElem)
!Sparsity=Sparsity/REAL(nDOFElem)/REAL(nDOFElem)
!WRITE(*,'(A,I6)')       ' Element number          : ', iElem
!WRITE(*,'(A,I6)')       ' Number of non zeros in U: ', nUNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in L: ', nLNonZeros(iElem)
!WRITE(*,'(A,I6)')       ' Number of non zeros in D: ', nDOFElem
!WRITE(*,'(A,F12.7,A)')  ' Level of Fill             ', sparsity*100.0,'%'
!print*

! Allocate
!nMTriangle(iElem)=2*(nDOFelem-1)
! simple version
nMTriangle=nDOFelem
ALLOCATE( IU(iElem)%Entry(nUNonZeros(iElem))  &
        , IU(iELEM)%JEntry(nUNonZeros(iElem)) &
        , IU(iELEM)%IEntry(nDOFElem+1) &
        , IL(iELEM)%Entry(nLNonZeros(iElem))  &
        , IL(iELEM)%JEntry(nLNonZeros(iElem)) &
        , IL(iELEM)%IEntry(nDOFElem+1) )

! nullify
DE(:,iElem)=0.
IL(iELem)%Entry=0.
IU(iELem)%Entry=0.

! U part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nDOFElem,1
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

! L part of ILU0
! simple version
jj=0
singleValue=.FALSE.
DO ii=1,nDOFElem,1
  first=.TRUE.
  lineNonZero=0
  iEntry=ii!+1
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
      EXIT
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

! Dioganal entires
DO ii=1,nDOFElem
  DE(ii,iElem)=ILU0(ii,ii)
END DO

END SUBROUTINE DiagCSR


SUBROUTINE GlobalCSR()
!===================================================================================================================================
! Store GlobalJacobian in CSR Format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars           ,ONLY:MatPattern,MatPatternElemID,ProcJacobian
USE MOD_LinearSolver_Vars      ,ONLY:nDOFGlobal,nDOFElem
USE MOD_CSR_Vars               ,ONLY:GlobalIA,GlobalJA,GlobalAA,epsZero,nonZerosGlobal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: r,s, ii,kk
INTEGER                     :: nonZero
INTEGER                     :: iElem, iRow, conElemID, conSideID
INTEGER                     :: lineNonZero,iEntry,LastLine
LOGICAL                     :: first
!===================================================================================================================================

! get global number of non-zeros
nonZerosGlobal=0
DO iElem=1,PP_nElems
  DO iRow=1,7
    IF(MatPattern(iRow,iElem).LT.0) CYCLE
    nonZero=GetNonZeros(nDOFElem,ProcJacobian(1:nDOFElem,1:nDOFElem,MatPattern(iRow,iElem),iElem))
    nonZerosGlobal=nonZerosGlobal+nonZero
  END DO ! iRow
END DO ! iElem

!----------------------------------------------------------------------------------------------------------------------------------
! Extendet CSR FORMAT
!----------------------------------------------------------------------------------------------------------------------------------

! Allocate
ALLOCATE( GlobalAA(1:NonZerosGlobal) &
        , GlobalJA(1:NonZerosGlobal) &
        , GlobalIA(1:nDOFGlobal+1)   )

iEntry=0
DO iElem=1,PP_nElems
  DO r=1,nDOFElem
    ii=nDOFElem*(iElem-1)+r
    first=.TRUE.
    lineNonZero=0
    DO iRow=1,7
      IF(MatPatternElemID(iRow,iElem).LT.0) CYCLE
      conElemID=MatPatternElemID(iRow,iElem)
      conSideID=MatPattern      (iRow,iElem)
      DO s=1,nDOFElem
        IF(ABS(ProcJacobian(r,s,conSideID,iElem)).GT.epsZero) THEN
          iEntry=iEntry+1
          kk = nDOFElem*(conElemID-1)+s
          GlobalAA(iEntry) = ProcJacobian(r,s,conSideID,iElem)
          GlobalJA(iEntry) = kk
          IF(GlobalJA(iEntry).EQ.0) print*,'iEntry,ii,kk',iEntry,ii,kk
          lineNonZero= lineNonZero+1
          IF(first)THEN
            GlobalIA(ii)=iEntry
            first=.FALSE.
          END IF
        END IF ! zeros
      END DO ! s
    END DO ! iRow
    IF(lineNonZero.EQ.0)THEN
      IF(ii.EQ.1)THEN
        GlobalIA(ii)=1
      ELSE
        GlobalIA(ii)=GlobalIA(ii-1)+LastLine
      END IF
    END IF
    LastLine=LineNonZero
  END DO ! r
END DO ! iElem
GlobalIA(nDOFGlobal+1)=GlobalIA(1)+nonZerosGlobal

END SUBROUTINE GlobalCSR

SUBROUTINE GlobalBCSR()
!===================================================================================================================================
! Store GlobalJacobian in CSR Format
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Precond_Vars           ,ONLY:MatPattern,MatPatternElemID, ProcJacobian, BlockSize, nBlockSize
USE MOD_LinearSolver_Vars      ,ONLY:nDOFElem
USE MOD_CSR_Vars               ,ONLY:GlobalIA,GlobalJA,GlobalBAA,epsZero,nonZerosGlobal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: r,s, ii,kk
INTEGER                     :: nonZero,BlockSquare
INTEGER                     :: iElem, iRow, conElemID, conSideID
INTEGER                     :: lineNonZero,iEntry,LastLine
LOGICAL                     :: first
!===================================================================================================================================

! get global number of non-zeros
nonZerosGlobal=0
DO iElem=1,PP_nElems
  DO iRow=1,7
    IF(MatPatternElemID(iRow,iElem).LT.0) CYCLE
    nonZero=GetBlockNonZeros(BlockSize,nDOFElem,ProcJacobian(1:nDOFElem,1:nDOFElem,MatPattern(iRow,iElem),iElem))
    nonZerosGlobal=nonZerosGlobal+nonZero
  END DO ! iRow
END DO ! iElem

!----------------------------------------------------------------------------------------------------------------------------------
! Extendet CSR FORMAT
!----------------------------------------------------------------------------------------------------------------------------------

! Allocate
ALLOCATE( GlobalBAA(1:BlockSize,1:BlockSize,1:NonZerosGlobal) &
        , GlobalJA(1:NonZerosGlobal) &
        , GlobalIA(1:nBlockSize+1)   )
BlockSquare=BlockSize*BlockSize
GlobalJA=0

iEntry=0
DO iElem=1,PP_nElems
  DO r=1,nDOFElem,BlockSize
    ii=nDOFElem*(iElem-1)+r
    first=.TRUE.
    lineNonZero=0
    DO iRow=1,7
      IF(MatPatternElemID(iRow,iElem).LT.0) CYCLE
      conElemID=MatPatternElemID(iRow,iElem)
      conSideID=MatPattern      (iRow,iElem)
      DO s=1,nDOFElem,BlockSize
        IF(SUM(ABS(ProcJacobian(r:r+BlockSize-1,s:s+BlockSize-1,conSideID,iElem)))/REAL(BlockSquare).GT.epsZero) THEN
        !IF(ABS(ProcJacobian(r,s,conSideID,iElem)).GT.epsZero) THEN
          iEntry=iEntry+1
          kk = nDOFElem*(conElemID-1)+s
          GlobalBAA(1:BlockSize,1:BlockSize,iEntry) = ProcJacobian(r:r+BlockSize-1,s:s+BlockSize-1,conSideID,iElem)
          GlobalJA(iEntry) = kk
          IF(GlobalJA(iEntry).EQ.0) print*,'iEntry,ii,kk',iEntry,ii,kk
          lineNonZero= lineNonZero+1
          IF(first)THEN
            GlobalIA(ii/BlockSize+1)=iEntry
            first=.FALSE.
          END IF
        END IF ! zeros
      END DO ! s
    END DO ! iRow
    IF(lineNonZero.EQ.0)THEN
      IF(ii.EQ.1)THEN
        GlobalIA(ii/BlockSize+1)=1
      ELSE
        GlobalIA(ii/BlockSize+1)=GlobalIA(ii/BlockSize)+LastLine
        !GlobalIA(ii)=GlobalIA(ii-1)+LastLine
      END IF
    END IF
    LastLine=LineNonZero
  END DO ! r
END DO ! iElem
GlobalIA(nBlockSize+1)=GlobalIA(1)+nonZerosGlobal

END SUBROUTINE GlobalBCSR

!SUBROUTINE DGTimeDerivative_weakFormCSR(t,tStage,tDeriv)!,U,Ut)
!!==================================================================================================================================
!! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
!! U and Ut are allocated
!!==================================================================================================================================
!! MODULES
!USE MOD_Globals
!USE MOD_Preproc
!USE MOD_LinearSolver_Vars, ONLY:nDOFElem
!USE MOD_CSR_Vars,      ONLY:SparseMatrix,nNonZeros
!USE MOD_Mesh_Vars,     ONLY:nBCSides,ElemToSide,SideToElem
!USE MOD_Precond_Vars,  ONLY:neighborElemID
!#if USE_MPI
!USE MOD_MPI_Vars
!USE MOD_MPI,           ONLY:StartExchangeMPIData,FinishExchangeMPIData
!USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
!#endif
!USE MOD_DG_Vars,       ONLY:U,Ut
!USE MOD_CSR_Vars,      ONLY:Ucsr,Utcsr
!
!! old
!USE MOD_DG_Vars,       ONLY: nTotalU,U_Plus,U_Minus,Flux
!USE MOD_SurfInt,       ONLY: SurfInt
!USE MOD_VolInt,        ONLY: VolInt
!USE MOD_ProlongToFace, ONLY: ProlongToFace
!USE MOD_FillFlux,      ONLY: FillFlux,FillFlux_BC
!USE MOD_Mesh_Vars,     ONLY: sJ,Elem_xGP,nSides,nInnerSides
!USE MOD_Equation,      ONLY: CalcSource
!USE MOD_Equation_Vars, ONLY: IniExactFunc
!
!
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!REAL,INTENT(IN)                 :: t,tStage
!INTEGER,INTENT(IN)              :: tDeriv
!!REAL,INTENT(IN)                 :: U(nDOFElem,1:PP_nElems)
!!REAL,INTENT(IN)                 :: U(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!REAL,INTENT(OUT)                :: Ut(nDOFElem,1:PP_nElems)
!!REAL,INTENT(OUT)                 :: Ut(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!!REAL                             :: Ut2(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER :: iElem,ii,k1,k2,ilocSide,SideID,jj,flip,NBElemID,NBlocSideID,iSide
!INTEGER :: i ,j,k, iVar
!INTEGER :: locSideID(2),ElemID(2),ms
!===================================================================================================================================
!
!! normal DG operator, no mpi
!! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
!!CALL ProlongToFace(U,U_Minus,U_Plus,doMPISides=.FALSE.)
!!
!!! null
!!Ut=0.
!!CALL VolInt(Ut)
!!
!!! fill the all surface fluxes on this proc
!!CALL FillFlux_BC(t,tDeriv,Flux)
!!CALL FillFlux(Flux,doMPISides=.FALSE.)
!!! compute surface integral contribution and add to ut
!!CALL SurfInt(Flux,Ut,doMPISides=.FALSE.)
!!
!!DO iElem=1,PP_nElems
!!  DO k=0,PP_N
!!    DO j=0,PP_N
!!      DO i=0,PP_N
!!        DO iVar=1,PP_nVar
!!          Ut(iVar,i,j,k,iElem) = - Ut(iVar,i,j,k,iElem) * sJ(i,j,k,iElem)
!!        END DO ! iVar
!!      END DO !i
!!    END DO !j
!!  END DO !k
!!END DO ! iElem=1,nElems
!
!
!
!!Ut=0.
!!Utcsr=0.
!!DO iElem=1,PP_nElems
!!  !CALL ApplyDGCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!  CALL ApplyDGCSR(0,iElem,iElem)
!!  DO ilocSide=1,6
!!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!!    IF(SideID.LE.nBCSides) CYCLE
!!    NBElemID=neighborElemID(ilocSide,iElem)
!!  !  CALL ApplyDGCSR(U(:,:,:,:,NBElemID),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!    CALL ApplyDGCSR(ilocSide,NBElemID,iElem)
!!  END DO ! ilocSide
!!END DO ! iElem
!
!
!!!Ut=0.
!!DO iElem=1,PP_nElems
!!! ilocSide=0
!!!!  CALL ApplyDGCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!  CALL VolCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),SparseMatrix(0,iElem)%Entry,  &
!!                                                 SparseMatrix(0,iElem)%iEntry, &
!!                                                 SparseMatrix(0,iElem)%jEntry, &
!!                                                 nNonZeros(0,iElem)            )
!!  DO ilocSide=1,6
!!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!!    IF(SideID.LE.nBCSides) CYCLE
!!    NBElemID=neighborElemID(ilocSide,iElem)
!!    CALL SurfCSR(U(:,:,:,:,NBElemID),Ut(:,:,:,:,iElem),SparseMatrix(ilocSide,iElem)%Entry,  &
!!                                                       SparseMatrix(ilocSide,iElem)%iEntry, &
!!                                                       SparseMatrix(ilocSide,iElem)%jEntry, &
!!                                                       nNonZeros   (ilocSide,iElem)         )
!!   END DO ! ilocSide
!!END DO ! iElem
!
!!Ut=0.
!DO iElem=1,PP_nElems
!  ilocSide=0
!!  CALL ApplyDGCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!  CALL VolCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),SparseMatrix(0,iElem)%Entry,  &
!!                                                 SparseMatrix(0,iElem)%iEntry, &
!!                                                 SparseMatrix(0,iElem)%jEntry, &
!!                                                 nNonZeros(0,iElem)            )
!
!  CALL VolCSR(iElem,iElem,SparseMatrix(0,iElem)%Entry,  &
!                          SparseMatrix(0,iElem)%iEntry, &
!                          SparseMatrix(0,iElem)%jEntry, &
!                          nNonZeros(0,iElem)            )
!
!    !SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    !IF(SideID.LE.nBCSides) CYCLE
!    !NBElemID=neighborElemID(ilocSide,iElem)
!    !CALL ApplyDGCSR(U(:,:,:,:,NBElemID),Ut(:,:,:,:,iElem),ilocSide,iElem)
!END DO ! iElem
!
!DO iSide=nBCSides+1,nSides
!  ElemID(1)     = SideToElem(S2E_ELEM_ID,iSide)
!  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,iSide)
!  ! neighbor side
!  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,iSide)
!  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
!
!  CALL SurfCSR(ElemID(2),ElemID(1),SparseMatrix(locSideID(1),ElemID(1))%Entry,  &
!                                   SparseMatrix(locSideID(1),ElemID(1))%iEntry, &
!                                   SparseMatrix(locSideID(1),ElemID(1))%jEntry, &
!                                   nNonZeros   (locSideID(1),ElemID(1))         )
!
!  CALL SurfCSR(ElemID(1),ElemID(2),SparseMatrix(locSideID(2),ElemID(2))%Entry,  &
!                                   SparseMatrix(locSideID(2),ElemID(2))%iEntry, &
!                                   SparseMatrix(locSideID(2),ElemID(2))%jEntry, &
!                                   nNonZeros   (locSideID(2),ElemID(2))         )
!
!  !CALL SurfCSR(U(:,:,:,:,ElemID(2)),Ut(:,:,:,:,ElemID(1)),SparseMatrix(locSideID(1),ElemID(1))%Entry,  &
!  !                                                        SparseMatrix(locSideID(1),ElemID(1))%iEntry, &
!  !                                                        SparseMatrix(locSideID(1),ElemID(1))%jEntry, &
!  !                                                        nNonZeros   (locSideID(1),ElemID(1))         )
!
!  !CALL SurfCSR(U(:,:,:,:,ElemID(1)),Ut(:,:,:,:,ElemID(2)),SparseMatrix(locSideID(2),ElemID(2))%Entry,  &
!  !                                                        SparseMatrix(locSideID(2),ElemID(2))%iEntry, &
!  !                                                        SparseMatrix(locSideID(2),ElemID(2))%jEntry, &
!  !                                                        nNonZeros   (locSideID(2),ElemID(2))         )
!  !CALL ApplyDGCSR(U(:,:,:,:,ElemID(2)),Ut(:,:,:,:,ElemID(1)),locSideID(1),ElemID(1))
!  !CALL ApplyDGCSR(U(:,:,:,:,ElemID(1)),Ut(:,:,:,:,ElemID(2)),locSideID(2),ElemID(2))
!END DO ! iElem
!
!!Ut=0.
!!DO iElem=1,PP_nElems
!!  ilocSide=0
!!!  CALL ApplyDGCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!  CALL VolCSR(U(:,:,:,:,iElem),Ut(:,:,:,:,iElem),SparseMatrix(0,iElem)%Entry,  &
!!                                                 SparseMatrix(0,iElem)%iEntry, &
!!                                                 SparseMatrix(0,iElem)%jEntry, &
!!                                                 nNonZeros(0,iElem)            )
!!
!!
!!    !SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!!    !IF(SideID.LE.nBCSides) CYCLE
!!    !NBElemID=neighborElemID(ilocSide,iElem)
!!    !CALL ApplyDGCSR(U(:,:,:,:,NBElemID),Ut(:,:,:,:,iElem),ilocSide,iElem)
!!END DO ! iElem
!!
!!DO iSide=nBCSides+1,nSides
!!  ElemID(1)     = SideToElem(S2E_ELEM_ID,iSide)
!!  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,iSide)
!!  ! neighbor side
!!  ElemID(2)     = SideToElem(S2E_NB_ELEM_ID,iSide)
!!  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
!!
!!  CALL SurfCSR(U(:,:,:,:,ElemID(2)),Ut(:,:,:,:,ElemID(1)),SparseMatrix(locSideID(1),ElemID(1))%Entry,  &
!!                                                          SparseMatrix(locSideID(1),ElemID(1))%iEntry, &
!!                                                          SparseMatrix(locSideID(1),ElemID(1))%jEntry, &
!!                                                          nNonZeros   (locSideID(1),ElemID(1))         )
!!
!!  CALL SurfCSR(U(:,:,:,:,ElemID(1)),Ut(:,:,:,:,ElemID(2)),SparseMatrix(locSideID(2),ElemID(2))%Entry,  &
!!                                                          SparseMatrix(locSideID(2),ElemID(2))%iEntry, &
!!                                                          SparseMatrix(locSideID(2),ElemID(2))%jEntry, &
!!                                                          nNonZeros   (locSideID(2),ElemID(2))         )
!!  !CALL ApplyDGCSR(U(:,:,:,:,ElemID(2)),Ut(:,:,:,:,ElemID(1)),locSideID(1),ElemID(1))
!!  !CALL ApplyDGCSR(U(:,:,:,:,ElemID(1)),Ut(:,:,:,:,ElemID(2)),locSideID(2),ElemID(2))
!!END DO ! iElem
!
!!! We have to take the inverse of the Jacobians into account
!!DO iElem=1,PP_nElems
!!  DO k=0,PP_N
!!    DO j=0,PP_N
!!      DO i=0,PP_N
!!        DO iVar=1,PP_nVar
!!          Ut(iVar,i,j,k,iElem) = - Ut(iVar,i,j,k,iElem) * sJ(i,j,k,iElem)
!!        END DO ! iVar
!!      END DO !i
!!    END DO !j
!!  END DO !k
!!END DO ! iElem=1,nElems
!
!
!
!!DO iElem=1,PP_nElems
!!  DO k=0,PP_N
!!    DO j=0,PP_N
!!      DO i=0,PP_N
!!        DO iVar=1,PP_nVar
!!          IF(ABS(Ut(iVar,i,j,k,iElem)-Ut2(ivar,i,j,k,iElem)).GT.1e-5)THEN
!!            WRITE(*,'(A,I3,I3,I3,I3,I3,E24.12)') 'iVar,i,j,k,iElem,err',iVar,i,j,k,iElem,ABS(Ut(iVar,i,j,k,iElem)-&
!                                            !Ut2(ivar,i,j,k,iElem))
!!          END IF
!!        END DO ! iVar
!!      END DO !i
!!    END DO !j
!!  END DO !k
!!END DO ! iElem=1,nElems
!!
!
!END SUBROUTINE DGTimeDerivative_weakFormCSR
!
!!SUBROUTINE ApplyDGCSR(Vin,Vout,locSideID,ElemID)
!!SUBROUTINE ApplyDGCSR(locSideID,NBElemID,ElemID)
!!
!===================================================================================================================================
! Apply CSR of DG-Operator
!===================================================================================================================================
!!! MODULES
!!USE MOD_PreProc
!!USE MOD_CSR_Vars                ,ONLY: SparseMatrix,nNonZeros
!!USE MOD_CSR_Vars                ,ONLY: Ucsr,Utcsr
!!USE MOD_Implicit_Vars           ,ONLY: nDOFelem
!!! IMPLICIT VARIABLE HANDLING
!!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!!! INPUT VARIABLES
!!!REAL,INTENT(IN)                 :: Vin(1:nDOFElem)
!!REAL,INTENT(IN)                 :: AA(nNonZeros)
!!INTEGER,INTENT(IN)              :: locSideID,ElemID,NBElemID
!-----------------------------------------------------------------------------------------------------------------------------------
!!! OUTPUT VARIABLES
!!!REAL,INTENT(INOUT)                 :: Vout(1:nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
!!! LOCAL VARIABLES
!!INTEGER                         :: ii,jj,k1,k2
!===================================================================================================================================
!!
!!DO ii=1,nDOFElem
!!  k1=SparseMatrix(locSideID,ElemID)%iEntry(ii)
!!  k2=SparseMatrix(locSideID,ElemID)%iEntry(ii+1)-1
!!  !Vout(ii)=DOT_PRODUCT(AA(k1:k2),Vin(JA(k1:k2)))
!!  !Vout(ii)=DOT_PRODUCT(SparseMatrix(locSideID,ElemID)%Entry(k1:k2),Vin(SparseMatrix(locSideID,ElemID)%Jentry(k1:k2)))
!!  DO jj=k1,k2
!!    Utcsr(ii,ElemID) = Utcsr(ii,ElemID) + &
!!                        SparseMatrix(locSideID,ElemID)%Entry(jj)*Ucsr(SparseMatrix(locSideID,ElemID)%Jentry(jj),NBElemID)
!!  END DO ! jj
!!END DO ! ii
!!
!!END SUBROUTINE ApplyDGCSR
!
!!SUBROUTINE VolCSR(Vin,Vout,AA,IA,JA,nNonZeros)
!SUBROUTINE VolCSR(iElem2,iElem,AA,IA,JA,nNonZeros)
!===================================================================================================================================
!! Apply CSR of DG-Operator
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_CSR_Vars                ,ONLY: Ucsr,Utcsr
!USE MOD_LinearSolver_Vars       ,ONLY: nDOFelem
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!REAL,INTENT(IN)                 :: Vin(1:nDOFElem),AA(nNonZeros)
!REAL,INTENT(IN)                 :: AA(nNonZeros)
!INTEGER,INTENT(IN)              :: nNonZeros,IA(nDOFElem+1),JA(nNonZeros)
!INTEGER,INTENT(IN)              :: iElem,iElem2
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!REAL,INTENT(INOUT)              :: Vout(1:nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                         :: ii,jj,k1,k2
!===================================================================================================================================
!
!DO ii=1,nDOFElem
!  k1=IA(ii)
!  k2=IA(ii+1)-1
!  Utcsr(ii,iElem)=DOT_PRODUCT(AA(k1:k2),Ucsr(JA(k1:k2),iElem2))
!  !Vout(ii)=DOT_PRODUCT(AA(k1:k2),Vin(JA(k1:k2)))
!!  DO jj=k1,k2
!!    Vout(ii) = Vout(ii) + AA(jj)*Vin(JA(jj))
!!  END DO ! jj
!END DO ! ii
!
!END SUBROUTINE VolCSR
!
!!SUBROUTINE SurfCSR(Vin,Vout,AA,IA,JA,nNonZeros)
!SUBROUTINE SurfCSR(iElem2,iElem,AA,IA,JA,nNonZeros)
!===================================================================================================================================
!! Apply CSR of DG-Operator
!===================================================================================================================================
!! MODULES
!USE MOD_PreProc
!USE MOD_CSR_Vars                ,ONLY: Ucsr,Utcsr
!USE MOD_LinearSolver_Vars           ,ONLY: nDOFelem
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!!REAL,INTENT(IN)                 :: Vin(1:nDOFElem),AA(nNonZeros)
!REAL,INTENT(IN)                 :: AA(nNonZeros)
!INTEGER,INTENT(IN)              :: nNonZeros,IA(nDOFElem+1),JA(nNonZeros)
!INTEGER,INTENT(IN)              :: iElem,iElem2
!-----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!REAL,INTENT(INOUT)              :: Vout(1:nDOFElem)
!-----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                         :: ii,jj,k1,k2
!===================================================================================================================================
!
!DO ii=1,nDOFElem
!  k1=IA(ii)
!  k2=IA(ii+1)-1
!  !Vout(ii)=Vout(ii) + DOT_PRODUCT(AA(k1:k2),Vin(JA(k1:k2)))
!  Utcsr(ii,iElem)=Utcsr(ii,iElem) + DOT_PRODUCT(AA(k1:k2),Ucsr(JA(k1:k2),iElem2))
!!  DO jj=k1,k2
!!    Vout(ii) = Vout(ii) + AA(jj)*Vin(JA(jj))
!!  END DO ! jj
!END DO ! ii
!
!END SUBROUTINE SurfCSR

FUNCTION GetNonZeros( n, Mat)
!===================================================================================================================================
! Computes the number of non-zeros of a matrix with size nxn
!    input: p,q in Slave-RHS, flip;
!   output: indices in Master-RHS
!===================================================================================================================================
! MODULES
USE MOD_CSR_Vars,       ONLY: epsZero
! IMPLICIT VARIABLE HANDLING
! INPUT VARIABLES
INTEGER,INTENT(IN) :: n
REAL,INTENT(IN)    :: Mat(n,n)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: GetNonZeros
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: nonZeros, r,s
!===================================================================================================================================

nonZeros=0
DO s=1,n
  DO r=1,n
    IF(ABS(Mat(r,s)).GT.epsZero) nonZeros=nonZeros+1
  END DO ! r
END DO ! s

GetNonZeros=nonZeros

END FUNCTION GetNonZeros

FUNCTION GetBlockNonZeros( b, n, Mat)
!===================================================================================================================================
! Computes the number of non-zeros of a matrix with size nxn
!    input: p,q in Slave-RHS, flip;
!   output: indices in Master-RHS
!===================================================================================================================================
! MODULES
USE MOD_CSR_Vars,       ONLY:epsZero
!USE MOD_Precond_Vars,   ONLY:BlockSize
! IMPLICIT VARIABLE HANDLING
! INPUT VARIABLES
INTEGER,INTENT(IN) :: n
INTEGER,INTENT(IN) :: b
REAL,INTENT(IN)    :: Mat(n,n)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: GetBlockNonZeros
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: nonZeros, r,s, b2
!===================================================================================================================================

b2=b*b
nonZeros=0
DO s=1,n,b
  DO r=1,n,b
    IF(SUM(ABS(Mat(r:r+b-1,s:s+b-1)))/REAL(b2).GT.epsZero) nonZeros=nonZeros+1
  END DO ! r
END DO ! s

GetBlockNonZeros=nonZeros

END FUNCTION GetBlockNonZeros


SUBROUTINE FinalizeCSR()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
!===================================================================================================================================
! MODULES
USE MOD_CSR_Vars    ,ONLY:nNonZeros,SparseMatrix, L_HatPlusMinus!, Ucsr, Utcsr!,DiagonalEntries
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( nNonZeros      )
!SDEALLOCATE( DiagonalEntries)
SDEALLOCATE( SparseMatrix   )
SDEALLOCATE( L_HatPlusMinus )
!SDEALLOCATE( Ucsr           )
!SDEALLOCATE( Utcsr          )

END SUBROUTINE FinalizeCSR

SUBROUTINE WriteMat(n,Mat)
!===================================================================================================================================
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_CSR_Vars,             ONLY: epsZero
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: n
REAL,INTENT(IN)              :: Mat (n,n)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: ii,kk
CHARACTER(LEN=255)           :: filename
!===================================================================================================================================

WRITE(Filename,'(A,I2.2,A)')'MyRank_',myRank,'_Mat.dat'
WRITE(*,*)'Debug Jacobian to:',TRIM(Filename)
OPEN (UNIT=103,FILE=TRIM(Filename),STATUS='REPLACE')

DO ii=1,n
  DO kk=1,n
    IF(ABS(Mat(ii,kk)).GT.epsZero)THEN
      WRITE(103,'(I12,2x,I12,2x,E23.16)') ii,kk,Mat(ii,kk)
    END IF
  END DO
END DO

CLOSE(103)

END SUBROUTINE WriteMat

END MODULE MOD_CSR
