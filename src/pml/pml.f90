#include "boltzplatz.h"


MODULE MOD_PML
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitPML
  MODULE PROCEDURE InitPML
END INTERFACE
INTERFACE FinalizePML
  MODULE PROCEDURE FinalizePML
END INTERFACE
INTERFACE CalcPMLSource
  MODULE PROCEDURE CalcPMLSource
END INTERFACE
INTERFACE PMLTimeDerivative
  MODULE PROCEDURE PMLTimeDerivative
END INTERFACE

PUBLIC::InitPML,FinalizePML,CalcPMLSource,PMLTimeDerivative
!===================================================================================================================================
CONTAINS

SUBROUTINE InitPML()
!===================================================================================================================================
!  Initialize perfectly matched layer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_PML_Vars,      ONLY: PMLRamp,PMLzeta,U2,U2t,DoPML,PMLzetaGlobal
USE MOD_PML_Vars,      ONLY: nPMLElems,ElemToPML,PMLToElem,isPMLElem
USE MOD_PML_Vars,      ONLY: nPMLFaces,FaceToPML,PMLToFace,isPMLFace
USE MOD_PML_Vars,      ONLY: nPMLInterFaces,FaceToPMLInter,PMLInterToFace,isPMLInterFace
USE MOD_PML_Vars,      ONLY: PMLzeta0,PMLalpha0,xyzPhysicalMinMax,PMLzetaShape,PMLRampLength,PMLspread,PMLwriteFields, PMLzetaNorm
USE MOD_Mesh_Vars,     ONLY: Elem_xGP,BCFace_xGP,nSides!,Face_xGP  ! for PML region: xyz position of the Gauss points and Face Gauss points
USE MOD_PML_Vars,      ONLY: PMLRampLength
USE MOD_PML_Vars,      ONLY: PMLzetaEff,PMLalpha
USE MOD_PML_Vars,      ONLY: PMLnVar
USE MOD_HDF5_output,           ONLY: GatheredWriteArray,GenerateFileSkeleton,WriteAttributeToHDF5,WriteHDF5Header
USE MOD_Mesh_Vars,             ONLY: MeshFile,nGlobalElems,offsetElem
USE MOD_Output_Vars,           ONLY: ProjectName
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iPMLElem,iFace,m
REAL                :: xyzMinMax(6),xi,L,XiN,delta(3),x,y,z
REAL                :: xyzMinMaxLoc(6)
REAL                :: zetaVec,zetaVecABS
INTEGER             :: iPMLFace,iPMLInterFace,nGlobalPMLElems,nGlobalPMLFaces,nGlobalPMLInterFaces
REAL                :: fLinear,fSinus,fPolynomial
INTEGER             :: DOFcount
! for HDF5 output
INTEGER(HID_T)                 :: Dset_ID
INTEGER                        :: nVal
CHARACTER(LEN=255)             :: FileName,StrVarNames(PP_nVar)
#ifdef MPI
#endif
REAL                 :: StartT,EndT
REAL                           :: OutputTime
REAL                           :: FutureTime
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT PML...'
!===================================================================================================================================
! Readin
!===================================================================================================================================
DoPML                  = GETLOGICAL('DoPML','.FALSE.')
PMLzeta0               = GETREAL('PMLzeta0','0.')
PMLalpha0              = GETREAL('PMLalpha0','0.')
xyzPhysicalMinMax(1:6) = GETREALARRAY('xyzPhysicalMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
IF(ALMOSTEQUAL(MAXVAL(xyzPhysicalMinMax),MINVAL(xyzPhysicalMinMax)))THEN
  xyzPhysicalMinMax(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
  print*,"no PML region supplied, setting xyzPhysicalMinMax=",xyzPhysicalMinMax
END IF
PMLzetaShape           = GETINT('PMLzetaShape','0')
PMLRampLength          = GETREAL('PMLRampLength','1.')
PMLspread              = GETINT('PMLspread','0')
PMLwriteFields         = GETINT('PMLwriteFields','0')
PMLzetaNorm            = GETLOGICAL('PMLzetaNorm','.FALSE.')
! caution, in current version read in in mesh
! only for Maxwell, PP_nVar=8

IF(.NOT.DoPML) THEN
  SWRITE(UNIT_stdOut,'(A)') ' PML region deactivated. '
#if PP_nVar == 1
  CALL abort(__STAMP__,&
      'Equation system does not support a PML!',999,999.)
#endif
#if PP_nVar == 4
  CALL abort(__STAMP__,&
      'Equation system does not support a PML!',999,999.)
#endif
  PMLnVar=0
  nPMLElems=0
  RETURN
ELSE
  PMLnVar=24
END IF

CALL FindElementInRegion(isPMLElem,xyzPhysicalMinMax,ElementIsInside=.FALSE.) ! find all elements in the PML region
CALL FindInterfaces(isPMLFace,isPMLInterFace,isPMLElem)                       ! find all faces in the PML region




! Get number of PML Elems
nPMLFaces = 0
nPMLInterFaces = 0
nPMLElems = 0

DO iFace=1,nSides
  IF(isPMLFace(iFace))THEN
    nPMLFaces=nPMLFaces+1
  END IF
END DO ! iFace

DO iFace=1,nSides
  IF(isPMLInterFace(iFace))THEN
    nPMLInterFaces=nPMLInterFaces+1
  END IF
END DO ! iFace

DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    nPMLElems=nPMLElems+1
  END IF
END DO ! iElem



IF(1.EQ.2)THEN
#ifdef MPI
nGlobalPMLElems=0
nGlobalPMLFaces=0
nGlobalPMLInterfaces=0
  IF(MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              nSides              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              nSides         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  IF(.NOT.MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              nSides              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              nSides         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  !=======================================================================================================================
  CALL MPI_REDUCE(nPMLElems     ,nGlobalPMLElems     ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(nPMLFaces     ,nGlobalPMLFaces     ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(nPMLInterFaces,nGlobalPMLInterFaces,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  IF(MPIroot) print *, "============================================================================================"
  IF(MPIroot) print *, "CALL MPI_REDUCE(nPMLElems,nGlobalPMLElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)" !testing 
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  !=======================================================================================================================
  IF(MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              nSides              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              nSides         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  IF(.NOT.MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              nSides              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              nSides         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,nSides             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
#else
  nGlobalPMLElems=nPMLElems
  nGlobalPMLFaces=nPMLFaces
  nGlobalPMLInterFaces=nPMLInterFaces
  SWRITE(UNIT_stdOut,'(A,I25,A)') ' Found ', nGlobalPMLElems,' nGlobalPMLElems inside of PML region.'
  SWRITE(UNIT_stdOut,'(A,I25,A)') ' Found ', nGlobalPMLFaces,' nGlobalPMLFaces inside of PML region.'
  SWRITE(UNIT_stdOut,'(A,I25,A)') ' Found ', nGlobalPMLInterFaces,' nGlobalPMLInterFaces inside of PML region.'
#endif /*MPI*/
END IF




ALLOCATE(ElemToPML(PP_nElems)&
        ,PMLToElem(nPMLElems))
ALLOCATE(FaceToPML(nSides)&
        ,PMLToFace(nPMLFaces))
ALLOCATE(FaceToPMLInter(nSides)&
        ,PMLInterToFace(nPMLInterFaces))

ElemToPML=0
PMLToElem=0
FaceToPML=0
PMLToFace=0
FaceToPMLInter=0
PMLInterToFace=0

! Create array with mapping
iPMLElem=0
DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    iPMLElem=iPMLElem+1
    ElemToPML(iElem) = iPMLElem
    PMLToElem(iPMLElem) = iElem
  END IF
END DO
iPMLFace=0
DO iFace=1,nSides
  IF(isPMLFace(iFace))THEN
    iPMLFace=iPMLFace+1
    FaceToPML(iFace) = iPMLFace
    PMLToFace(iPMLFace) = iFace
  END IF
END DO
iPMLInterFace=0
DO iFace=1,nSides
  IF(isPMLInterFace(iFace))THEN
    iPMLInterFace=iPMLInterFace+1
    FaceToPMLInter(iFace) = iPMLInterFace
    PMLInterToFace(iPMLInterFace) = iFace
  END IF
END DO
!===================================================================================================================================
! allocate field variables
!===================================================================================================================================

!DEALLOCATE(isPMLElem) ! auskommentiert, da fuer chi scaling in PML region benoetigt
ALLOCATE(PMLRamp       (0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLzeta   (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
! CFS-PML auxiliary variables (24xNxNxNxnElemsx8)
ALLOCATE(U2       (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))        
ALLOCATE(U2t      (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLzetaEff(1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLalpha  (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
!print *, "PMLzeta size, shape and location",SIZE(PMLzeta),shape(PMLzeta),loc(PMLzeta)
!print *, "U2 size, shape and location     ",SIZE(U2),shape(U2),loc(U2)
! zero and unity
! the local DG solution
ALLOCATE(PMLzetaGlobal(3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
PMLzetaGlobal=0.
PMLzeta=0.
U2 =0.
U2t=0.
PMLRamp=1.

PMLzetaEff=0.
PMLalpha=PMLalpha0
DOFcount=0
!print *, "PMLzeta size, shape and location",SIZE(PMLzeta),shape(PMLzeta),loc(PMLzeta)
!print *, "U2 size, shape and location     ",SIZE(U2),shape(U2),loc(U2)
!===================================================================================================================================
!determine PMLzeta values for each interpolation point
!===================================================================================================================================
xyzMinMaxloc(:) = (/MINVAL(BCFace_xGP(1,:,:,:)),MAXVAL(BCFace_xGP(1,:,:,:)),MINVAL(BCFace_xGP(2,:,:,:)),&
                    MAXVAL(BCFace_xGP(2,:,:,:)),MINVAL(BCFace_xGP(3,:,:,:)),MAXVAL(BCFace_xGP(3,:,:,:))/)

#ifdef MPI
   CALL MPI_ALLREDUCE(xyzMinMaxloc(1),xyzMinMax(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(xyzMinMaxloc(2),xyzMinMax(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(xyzMinMaxloc(3),xyzMinMax(3), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(xyzMinMaxloc(4),xyzMinMax(4), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(xyzMinMaxloc(5),xyzMinMax(5), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(xyzMinMaxloc(6),xyzMinMax(6), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
#else
   xyzMinMax=xyzMinMaxloc
#endif /*MPI*/

print*,"xyzMinMax - X",xyzMinMax(1:2)
print*,"xyzMinMax - Y",xyzMinMax(3:4)
print*,"xyzMinMax - Z",xyzMinMax(5:6)
!read*
print*,"xyzPhysicalMinMax - X",xyzPhysicalMinMax(1:2)
print*,"xyzPhysicalMinMax - Y",xyzPhysicalMinMax(3:4)
print*,"xyzPhysicalMinMax - Z",xyzPhysicalMinMax(5:6)
!read*


!print *, "xyzMinMax",xyzMinMax
SELECT CASE (PMLzetaShape)
CASE(0) ! Constant Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1) .OR. Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN        
          PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0
        END IF
        ! y-PML region
        IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3) .OR. Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN        
          PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0
        END IF
        ! z-PML region
        IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5) .OR. Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN        
          PMLzeta(3,i,j,k,ElemToPML(iElem)) = PMLzeta0
         ! 1.) DEBUGPML: changed from 3 to 1+2 (i.e. that for z-dir PMLs the x- and y- values are set)
         ! because generally the Gyrotron only has z-PML
         ! 3.) 2nd test: apply zeta to all 3 directions
         !PMLzeta(1,i,j,k,ElemToPML(iElem)) = PMLzeta0
         !PMLzeta(2,i,j,k,ElemToPML(iElem)) = PMLzeta0


        END IF
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(1) ! Linear Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
!        IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN        
!          XiN = (ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1)))/(ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1)))
!          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
!          XiN = (ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2)))/(ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2)))
!          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        END IF
!        ! y-PML region
!        IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN        
!          XiN = (ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3)))/(ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3)))
!          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN    
!          XiN = (ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4)))/(ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4)))
!          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        END IF
!        ! z-PML region
!        IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN        
!          XiN = (ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5)))/(ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5)))
!          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN    
!          XiN = (ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6)))/(ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6)))
!          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
!        END IF



        DO m=1,3 ! m=x,y,z
          IF (Elem_xGP(m,i,j,k,iElem) .LT. xyzPhysicalMinMax(2*m-1)) THEN ! 1, 3, 5 
            XiN = ABS (( Elem_xGP(m,i,j,k,iElem)-xyzPhysicalMinMax(2*m-1)   ) / &
                       (        xyzMinMax(2*m-1)-xyzPhysicalMinMax(2*m-1)   ))
            PMLzeta(m,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
          ELSEIF (Elem_xGP(m,i,j,k,iElem) .GT. xyzPhysicalMinMax(2*m)) THEN ! 2, 4, 6
            XiN = ABS (( Elem_xGP(m,i,j,k,iElem)-xyzPhysicalMinMax(2*m)   ) / &
                       (        xyzMinMax(2*m)  -xyzPhysicalMinMax(2*m)   ))
            PMLzeta(m,i,j,k,ElemToPML(iElem))   = PMLzeta0*fLinear(XiN)
          END IF
        END DO
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(2) ! Sinusoidal  Distribution of the Damping Coefficient
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN
          xi                                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1))
          L                                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength) 
        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
          xi                                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2))
          L                                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength) 
        END IF
        ! y-PML region
        IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN
          xi                                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3))
          L                                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength) 
        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN
          xi                                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4))
          L                                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength) 
        END IF
        ! z-PML region
        IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN
          xi                                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5))
          L                                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength) 
        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN
          xi                                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6))
          L                                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fSinus(xi/L,PMLRampLength)
        ENDIF 
  END DO; END DO; END DO; END DO !iElem,k,i,j
CASE(3) ! polynomial
        ! f''(0) = 0, f'(1) = 0
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        ! x-PML region
        IF (Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(1)) THEN
          xi                                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(1))
          L                                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
          XiN                                 = xi/L
          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
        ELSEIF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(2)) THEN
          xi                                  = ABS(Elem_xGP(1,i,j,k,iElem))-ABS(xyzPhysicalMinMax(2))
          L                                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
          XiN                                 = xi/L
          PMLzeta(1,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
        END IF
        ! y-PML region
        IF (Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(3)) THEN
          xi                                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(3))
          L                                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
          XiN                                 = xi/L
          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
        ELSEIF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(4)) THEN
          xi                                  = ABS(Elem_xGP(2,i,j,k,iElem))-ABS(xyzPhysicalMinMax(4))
          L                                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
          XiN                                 = xi/L
          PMLzeta(2,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
        END IF
        ! z-PML region
        IF (Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(5)) THEN
          xi                                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(5))
          L                                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
          XiN                                 = xi/L
          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
        ELSEIF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(6)) THEN
          xi                                  = ABS(Elem_xGP(3,i,j,k,iElem))-ABS(xyzPhysicalMinMax(6))
          L                                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
          XiN                                 = xi/L
          PMLzeta(3,i,j,k,ElemToPML(iElem))   = PMLzeta0*fPolynomial(Xin)
      ENDIF 
  END DO; END DO; END DO; END DO !iElem,k,i,j

CASE DEFAULT
!  CALL abort(__STAMP__,'Shape function for damping coefficient in PML region not specified!',999,999.)
END SELECT ! PMLzetaShape



!===================================================================================================================================
! Modification to zeta values
!===================================================================================================================================
! cutout
! 2.) DEBUGPML: exclude the cylinder from the PML, where particles travel -> r < 4mm
!DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!  IF (isPMLElem(iElem)) THEN
!    IF (SQRT(Elem_xGP(1,i,j,k,iElem)**2+Elem_xGP(2,i,j,k,iElem)**2)<4E-3)THEN
!      !PMLzeta(1,i,j,k,ElemToPML(iElem)) = 0.0
!      !PMLzeta(2,i,j,k,ElemToPML(iElem)) = 0.0
!      !PMLzeta(3,i,j,k,ElemToPML(iElem)) = 0.0
!      ! 3.) set the whole cell to zero
!      PMLzeta(:,:,:,:,ElemToPML(iElem)) = 0.0
!      DOFcount=DOFcount+1
!    END IF
!  END IF
!END DO; END DO; END DO; END DO !iElem,k,i,j
!IF(DOFcount.GT.0)THEN
!  print*,"CutOut cylidner: DOFcount = ",DOFcount
!END IF

!Spreading: Set All PMLzeta values for a direction to PMLzeta
IF (PMLspread.EQ.1) THEN
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        IF (PMLzeta(1,i,j,k,iPMLElem) .GT. 0.0) PMLzeta(:,i,j,k,iPMLElem)=PMLzeta(1,i,j,k,iPMLElem)
        IF (PMLzeta(2,i,j,k,iPMLElem) .GT. 0.0) PMLzeta(:,i,j,k,iPMLElem)=PMLzeta(2,i,j,k,iPMLElem)
        IF (PMLzeta(3,i,j,k,iPMLElem) .GT. 0.0) PMLzeta(:,i,j,k,iPMLElem)=PMLzeta(3,i,j,k,iPMLElem) 
  END DO; END DO; END DO; END DO !iPMLElem,k,i,j
END IF

!===================================================================================================================================
! Modification to zeta values
!===================================================================================================================================
!PMLzetaNorm=.TRUE.
! Normalizing: recalculate zeta if multiple direction
IF (PMLzetaNorm) THEN
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        zetaVecABS=SQRT(PMLzeta(1,i,j,k,iPMLElem)**2 &
                       +PMLzeta(2,i,j,k,iPMLElem)**2 &
                       +PMLzeta(3,i,j,k,iPMLElem)**2 )
        zetaVec=MAX(PMLzeta(1,i,j,k,iPMLElem),0.)
        zetaVec=MAX(PMLzeta(2,i,j,k,iPMLElem),zetaVec)
        zetaVec=MAX(PMLzeta(3,i,j,k,iPMLElem),zetaVec)
        PMLzeta(:,i,j,k,iPMLElem) = PMLzeta(:,i,j,k,iPMLElem)/zetaVecABS*zetaVec
  END DO; END DO; END DO; END DO !iPMLElem,k,i,j
END IF

!===================================================================================================================================
! CFS-PML formulation: calculate zeta  using the complex frequency shift
!===================================================================================================================================
!ACHTUNG: da am interface pmlzeta=0 wird somit durch null geteilt
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,3 ! m=x,y,z
    PMLzetaEff(m,i,j,k,iPMLElem) = ( PMLalpha(m,i,j,k,iPMLElem)+PMLzeta(m,i,j,k,iPMLElem) )
  END DO
END DO; END DO; END DO; END DO !iPMLElem,k,i,j
DEALLOCATE(PMLalpha)

!===================================================================================================================================
! create global zeta field for parallel output of zeta distribution
!===================================================================================================================================
IF(PMLzeta0.GT.0)THEN
  DO iElem=1,PP_nElems
    IF(isPMLElem(iElem))THEN
      PMLzetaGlobal(:,:,:,:,iElem)=PMLzeta(:,:,:,:,ElemToPML(iElem))/PMLzeta0
  !print*,"PMLzetaGlobal",PMLzetaGlobal
    END IF
  END DO!iElem
END IF

!===================================================================================================================================
! determine Elem_xGP distance to PML interface for PMLRamp
!===================================================================================================================================
!DO iPMLElem=1,nPMLElems; DO p=0,PP_N; DO q=0,PP_N
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! x-PML region
  !x = Face_xGP(1,p,q,PMLToFace(iPMLFace))
  !y = Face_xGP(2,p,q,PMLToFace(iPMLFace))
  !z = Face_xGP(3,p,q,PMLToFace(iPMLFace))
  x = Elem_xGP(1,i,j,k,PMLToElem(iPMLElem))
  y = Elem_xGP(2,i,j,k,PMLToElem(iPMLElem))
  z = Elem_xGP(3,i,j,k,PMLToElem(iPMLElem))
  delta=0.

  ! x-PML region
  IF (x .LT. xyzPhysicalMinMax(1)) THEN
    xi                  = ABS(x)-ABS(xyzPhysicalMinMax(1))
    L                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
  ELSEIF (x .GT. xyzPhysicalMinMax(2)) THEN
    xi                  = ABS(x)-ABS(xyzPhysicalMinMax(2))
    L                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
  ELSE
    xi=0
    L=1
  END IF
  delta(1)=MAXVAL((/0.,xi/L/))
  ! y-PML region
  IF (y .LT. xyzPhysicalMinMax(3)) THEN
    xi                  = ABS(y)-ABS(xyzPhysicalMinMax(3))
    L                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
  ELSEIF (y .GT. xyzPhysicalMinMax(4)) THEN
    xi                  = ABS(y)-ABS(xyzPhysicalMinMax(4))
    L                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
  ELSE
    xi=0
    L=1
  END IF
  delta(2)=MAXVAL((/0.,xi/L/))
  ! x-PML region
  IF (z .LT. xyzPhysicalMinMax(5)) THEN
    xi                  = ABS(z)-ABS(xyzPhysicalMinMax(5))
    L                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
  ELSEIF (z .GT. xyzPhysicalMinMax(6)) THEN
    xi                  = ABS(z)-ABS(xyzPhysicalMinMax(6))
    L                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
  ELSE
    xi=0
    L=1
  END IF
  delta(3)=MAXVAL((/0.,xi/L/))
  ! set the ramp value from 1 down to 0
  !PMLRamp(p,q,iPMLFace)=1.-( MAXVAL(delta)-SIN(2*ACOS(-1.)*MAXVAL(delta))/(2*ACOS(-1.)) )
  PMLRamp(i,j,k,iPMLElem) = 1. - fLinear(MAXVAL(delta))

  ! set the ramp value from 1 down to 0.82 (measured power loss)
  ! add ramp from 0 to 0.82 (power drain 30GHz Gyrotron over 2mm PML)
  !PMLRamp(i,j,k,iPMLElem) = PMLRamp(i,j,k,iPMLElem) + 0.82*fLinear(MAXVAL(delta))
!END DO; END DO; END DO !iFace,p,q
END DO; END DO; END DO; END DO !iPMLElem,k,i,j

!===================================================================================================================================
! write PMLRamp field
!===================================================================================================================================
IF (PMLwriteFields.EQ.1) THEN 
  OPEN(unit=110,file='PMLRamp.dat',status='unknown')
  !DO iPMLFace=1,nPMLFaces; DO p=0,PP_N; DO q=0,PP_N
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    write(110,'(4(ES15.7,2X))')&
                       Elem_xGP(1,i,j,k,PMLToElem(iPMLElem)), &
                       Elem_xGP(2,i,j,k,PMLToElem(iPMLElem)), &
                       Elem_xGP(3,i,j,k,PMLToElem(iPMLElem)), &
                       PMLRamp(i,j,k,iPMLElem)
  !END DO; END DO; END DO !iFace,p,q
  END DO; END DO; END DO; END DO !iPMLElem,k,i,j
  CLOSE(110)
END IF
!===================================================================================================================================
! write PMLInterFace
!===================================================================================================================================
!IF (PMLwriteFields.EQ.1) THEN 
!  OPEN(unit=110,file='PMLInterFace.dat',status='unknown')
!  DO iPMLInterFace=1,nPMLInterFaces; DO p=0,PP_N; DO q=0,PP_N
!    write(110,'(4(ES15.7,2X))')&
!                       Face_xGP(1,p,q,PMLInterToFace(iPMLInterFace)), &
!                       Face_xGP(2,p,q,PMLInterToFace(iPMLInterFace)), &
!                       Face_xGP(3,p,q,PMLInterToFace(iPMLInterFace)), &
!                       PMLRamp(p,q,FaceToPML(PMLInterToFace(iPMLInterFace)))
!  END DO; END DO; END DO !iFace,p,q
!  CLOSE(110)
!END IF
!===================================================================================================================================
! write PMLzeta field
!===================================================================================================================================
IF (PMLwriteFields.EQ.1) THEN
  OPEN(unit=110,file='PMLzeta.dat',status='unknown')
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        write(110,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)')&
                   Elem_xGP(1,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   Elem_xGP(2,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   Elem_xGP(3,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   PMLzeta(1,i,j,k,iPMLElem),'  ',             &
                   PMLzeta(2,i,j,k,iPMLElem),'  ',             &
                   PMLzeta(3,i,j,k,iPMLElem)
  END DO; END DO; END DO; END DO !iPMLElem,k,j,i
  CLOSE(110)
END IF
!===================================================================================================================================
! write PMLzetaEff field
!===================================================================================================================================
IF (PMLwriteFields.EQ.1) THEN
  OPEN(unit=110,file='PMLzetaEff.dat',status='unknown')
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
        write(110,'(ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7,A,ES15.7)')&
                   Elem_xGP(1,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   Elem_xGP(2,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   Elem_xGP(3,i,j,k,PMLToElem(iPMLElem)),'  ', &
                   PMLzetaEff(1,i,j,k,iPMLElem),'  ',             &
                   PMLzetaEff(2,i,j,k,iPMLElem),'  ',             &
                   PMLzetaEff(3,i,j,k,iPMLElem)
  END DO; END DO; END DO; END DO !iPMLElem,k,j,i
  CLOSE(110)
END IF

return
!===================================================================================================================================
! write PMLzetaGlobal field to HDF5 file
!===================================================================================================================================
OutputTime=0.0
FutureTime=0.0
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PMLZetaGlobal TO HDF5 FILE...'
END IF
StartT=BOLTZPLATZTIME()
StrVarNames(1)='PMLzetaGlobalX'
StrVarNames(2)='PMLzetaGlobalY'
StrVarNames(3)='PMLzetaGlobalZ'


! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_PMLZetaGlobal',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('PMLZetaGlobal',3,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif


!#ifdef MPI
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.)
!#else
!  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
!#endif

! Write DG solution ----------------------------------------------------------------------------------------------------------------
nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
!print*,"nGlobalElems=",nGlobalElems
!read*

! muss DG solution heißen, damit es vom Paraview State Reader erkannt wird
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/3,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/3,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=PMLzetaGlobal)



! Close the dataset and property list.
CALL H5DCLOSE_F(Dset_id, iError)

! Close the file.
CALL CloseDataFile()

IF(MPIROOT)THEN
  EndT=BOLTZPLATZTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF

END SUBROUTINE InitPML



SUBROUTINE CalcPMLSource()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY: Ut
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLToElem
USE MOD_PML_Vars,      ONLY: PMLzeta,U2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iPMLElem,m
!===================================================================================================================================
! sources for the standard variables
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,8
    Ut(m,i,j,k,PMLToElem(iPMLElem)) = Ut(m,i,j,k,PMLToElem(iPMLElem))  &
                                     -PMLzeta(1,i,j,k,iPMLElem)*U2(m*3-2,i,j,k,iPMLElem) &   ! = 1,4,7,10,13,16,19,22
                                     -PMLzeta(2,i,j,k,iPMLElem)*U2(m*3-1,i,j,k,iPMLElem) &   ! = 2,5,8,11,12,17,20,23
                                     -PMLzeta(3,i,j,k,iPMLElem)*U2(m*3  ,i,j,k,iPMLElem)     ! = 3,6,9,12,15,18,21,24
  END DO
END DO; END DO; END DO !nPMLElems,k,j,i
END DO
END SUBROUTINE CalcPMLSource


SUBROUTINE PMLTimeDerivative()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars,      ONLY: U2,U2t
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLToElem,PMLnVar
USE MOD_Mesh_Vars,     ONLY: sJ
USE MOD_PML_Vars,      ONLY: PMLzetaEff
USE MOD_Equation_Vars, ONLY: fDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iPMLElem,iPMLVar
!===================================================================================================================================
! We have to take the inverse of the Jacobians into account
! the '-' sign is due to the movement of the term to the right-hand-side of the equation
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO iPMLVar=1,PMLnVar
    U2t(iPMLVar,i,j,k,iPMLElem) = - U2t(iPMLVar,i,j,k,iPMLElem) * sJ(i,j,k,PMLToElem(iPMLElem))
  END DO
END DO; END DO; END DO; END DO !nPMLElems,k,j,i


! Add Source Terms
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  U2t(1 : 3,i,j,k,iPMLElem) = U2t(1 : 3,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(1 : 3,i,j,k,iPMLElem)
  U2t(4 : 6,i,j,k,iPMLElem) = U2t(4 : 6,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(4 : 6,i,j,k,iPMLElem)
  U2t(7 : 9,i,j,k,iPMLElem) = U2t(7 : 9,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(7 : 9,i,j,k,iPMLElem)
  U2t(10:12,i,j,k,iPMLElem) = U2t(10:12,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(10:12,i,j,k,iPMLElem)
  U2t(13:15,i,j,k,iPMLElem) = U2t(13:15,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(13:15,i,j,k,iPMLElem)
  U2t(16:18,i,j,k,iPMLElem) = U2t(16:18,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(16:18,i,j,k,iPMLElem)
  U2t(19:21,i,j,k,iPMLElem) = U2t(19:21,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(19:21,i,j,k,iPMLElem)
  U2t(22:24,i,j,k,iPMLElem) = U2t(22:24,i,j,k,iPMLElem) - PMLzetaEff(1:3,i,j,k,iPMLElem) * U2(22:24,i,j,k,iPMLElem)
END DO; END DO; END DO; END DO !nPMLElems,k,j,i


! 1.) DEBUGPML: apply the damping factor also to PML source terms
! copied from: U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
!U2 = U2 * fDamping

! 2.) DEBUGPML: apply the damping factor only to PML variables for Phi_E and Phi_B
!               to prevent charge-related instabilities (accumulation of divergence compensation over time)
U2(19:24,:,:,:,:) = fDamping* U2(19:24,:,:,:,:) 

END SUBROUTINE PMLTimeDerivative


SUBROUTINE FindElementInRegion(isElem,region,ElementIsInside)
!===================================================================================================================================
! Check 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,       ONLY: abort
USE MOD_Mesh_Vars,     ONLY: Elem_xGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: ElementIsInside
REAL,INTENT(IN)    :: region(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isElem(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,m
!===================================================================================================================================
ALLOCATE(isElem(1:PP_nElems))
isElem=.FALSE.
! PML elements are inside of the xyPMLMinMax region
print*,"Checking region:", region
IF(ElementIsInside)THEN
  isElem(:)=.TRUE.
  print*,"for elemens inside"
ELSE
  isElem(:)=.FALSE.
  print*,"for elemens outside"
END IF

DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,3 ! m=x,y,z
    IF ( (Elem_xGP(m,i,j,k,iElem) .LT. region(2*m-1)) .OR. & ! 1,3,5
         (Elem_xGP(m,i,j,k,iElem) .GT. region(2*m)) ) THEN   ! 2,4,6 ! element is outside
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outisde the region
    END IF
  END DO
END DO; END DO; END DO; END DO !iElem,k,j,i

IF(ElementIsInside)THEN
  print*,"No. of elements INSIDE region: ",COUNT(isElem)
ELSE
  print*,"No. of elements OUTSIDE region: ",COUNT(isElem)
END IF

END SUBROUTINE  FindElementInRegion


SUBROUTINE FindInterfaces(isFace,isInterFace,isElem)
!===================================================================================================================================
! Check if a face is in a special region (e.g. PML) and/or connects a special region (e.g. PML) to the physical region
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,       ONLY: abort,myrank,MPI_COMM_WORLD
USE MOD_Mesh_Vars,     ONLY: SideToElem,ElemToSide,nSides,nBCSides
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: isElem(1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isFace(:)
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isInterFace(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                 :: Plus(:,:,:,:),Minus(:,:,:,:)
INTEGER                          :: iElem,ilocSide,iSide,SideID
INTEGER                          :: ElemID(2)
LOGICAL                          :: printInfo
!===================================================================================================================================
ALLOCATE(isFace(1:nSides))
ALLOCATE(isInterFace(1:nSides))
isFace=.FALSE.
isInterFace=.FALSE.
printInfo=.FALSE.

! Check each element for being part of the, e.g. PML, region: set each of the 6 sides to be .TRUE.
DO iElem=1,PP_nElems
  IF(.NOT.isElem(iElem))CYCLE
  DO ilocSide =1,6
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    isFace(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,PP_nElems

#ifdef MPI
! ---------------------------------------------
! For MPI sides send the info to all other procs
! use numbering:    Plus+Minus  = 1: isFace
!                                 2: isInterFace
!                                 0: normal face in physical region
ALLOCATE(Plus(1,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Minus(1,0:PP_N,0:PP_N,1:nSides))
Plus=0.
Minus=0.
DO iSide=1,nBCSides ! 1.) do BC sides - when they are in a special region (e.g. PML) they need special treatment
  IF(isFace(iSide))THEN
    Plus( 1,:,:,iSide)=1. ! set my side true
    Minus(1,:,:,iSide)=1. ! set my neighbor side true
  ELSE
    Minus(1,:,:,iSide)=5. ! for BC sides set value that this side cannot be an InterFace!
  END IF
END DO
DO iSide=nBCSides+1,nSides ! 2.) do inner sides
  IF(isFace(iSide)) THEN
    Plus(1,:,:,iSide)=1.  ! set my side true
    Minus(1,:,:,iSide)=1. ! set my neighbor side true
  END IF
END DO
IF(printInfo)THEN
print*,"vorher (Send/Receive)"
print*,"Plus ",NINT(Plus(1,1,1,:))
print*,"Minus",NINT(Minus(1,1,1,:))
read*
END IF

! send my info to neighbor 
CALL StartReceiveMPIData(1,Minus(1,0:PP_N,0:PP_N,SideID_plus_lower:SideID_plus_upper) &
                                                ,SideID_plus_lower,SideID_plus_upper,RecRequest_U,SendID=2) ! Receive MINE
CALL StartSendMPIData(1,Plus(1,0:PP_N,0:PP_N,SideID_plus_lower:SideID_plus_upper) &
                                            ,SideID_plus_lower,SideID_plus_upper,SendRequest_U,SendID=2) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE

! add Minus to Plus and send
Plus=Plus+Minus
IF(printInfo)THEN
print*,"nachher (Send/Receive)"
print*,"Plus ",NINT(Plus(1,1,1,:))
read*
END IF

CALL StartReceiveMPIData(1,Plus,1,nSides,RecRequest_Flux ,SendID=1) ! Receive MINE
CALL StartSendMPIData(   1,Plus,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR

DO iSide=1,nSides ! get MPI Interfaces
  IF(Plus(1,0,0,iSide).EQ.1)THEN ! check MPI sides: if one side is .TRUE. and the other is .FALSE. -> InterFace
    isFace(iSide)=.TRUE. ! when my side is not PML but neighbor is PML
    isInterFace(iSide)=.TRUE.
  END IF
END DO

DEALLOCATE(Plus,Minus)
#endif /*MPI*/
! ---------------------------------------------

! and fill non-mpi sides: get local Interfaces
DO iSide=1,nSides
  IF(isInterFace(iSide)) CYCLE ! MPI sides are already finished
  IF(isFace(iSide))THEN
    ElemID(1)=SideToElem(S2E_ELEM_ID,iSide)    ! SideID
    ElemID(2)=SideToElem(S2E_NB_ELEM_ID,iSide) ! SideID
!print*,"ElemID(:)",ElemID
    IF(ElemID(1).LT.1) CYCLE
    IF(ElemID(2).LT.1) CYCLE ! neighbor is boundary side (either BC or MPI side)
    IF((isElem(ElemID(1)).AND. .NOT.isElem(ElemID(2))) .OR. &
  (.NOT.isElem(ElemID(1)).AND.      isElem(ElemID(2))) )THEN
      isInterFace(iSide)=.TRUE.
   END IF
  END IF
END DO

! test
DO iSide=nSides,nSides
  print*,"myrank=",myrank,"--------isInterFace(",iSide,")=",isInterFace(iSide),"  of total= ",COUNT(isInterFace),&
          "    isFace(",iSide,")=",     isFace(iSide),"  of total= ",COUNT(isFace)
END DO

!read*
!stop
!print*,"should be total 16 (PML-interfaces) and total 100 (PML-faces)"
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
!stop
END SUBROUTINE  FindInterFaces


SUBROUTINE FinalizePML()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_PML_Vars,            ONLY: PMLzeta,U2,U2t
USE MOD_PML_Vars,            ONLY: ElemToPML,PMLToElem,DoPML,isPMLElem,isPMLFace,PMLToFace,FaceToPML
USE MOD_PML_Vars,            ONLY: PMLRamp
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!RETURN
IF(.NOT.DoPML) RETURN
SDEALLOCATE(PMLzeta)
SDEALLOCATE(U2)
SDEALLOCATE(U2t)
SDEALLOCATE(PMLToElem)
SDEALLOCATE(ElemToPML)
SDEALLOCATE(PMLToFace)
SDEALLOCATE(FaceToPML)
SDEALLOCATE(PMLRamp)
SDEALLOCATE(isPMLElem)
SDEALLOCATE(isPMLFace)
END SUBROUTINE FinalizePML


END MODULE MOD_PML
























!===================================================================================================================================
! local SUBROUTINES and FUNCTIONS


REAL FUNCTION fLinear(x)
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_PML_Vars,            ONLY: PMLRampLength 
! IMPLICIT VARIABLE HANDLING 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.PMLRampLength) THEN
  x_temp = x/PMLRampLength
  fLinear = x_temp
ELSE
  fLinear = 1.
END IF
END FUNCTION fLinear


REAL FUNCTION fSinus(x,RampLength)
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
!USE MOD_PML_Vars,            ONLY: PMLRampLength
! IMPLICIT VARIABLE HANDLING 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
REAL, INTENT(IN) :: RampLength
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.RampLength) THEN
  x_temp = x/RampLength
  fSinus = x_temp-SIN(2*ACOS(-1.)*x_temp)/(2*ACOS(-1.))
ELSE
  fSinus = 1.
END IF
END FUNCTION fSinus



REAL FUNCTION fPolynomial(x)
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_PML_Vars,            ONLY: PMLRampLength
! IMPLICIT VARIABLE HANDLING 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.PMLRampLength) THEN
  x_temp = x/PMLRampLength
  fPolynomial = -3*x_temp**4+4*x_temp**3
ELSE
  fPolynomial = 1.
END IF
END FUNCTION fPolynomial


