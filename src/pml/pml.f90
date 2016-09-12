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
INTERFACE PMLsetZero
  MODULE PROCEDURE PMLsetZero
END INTERFACE

PUBLIC::InitPML,FinalizePML,CalcPMLSource,PMLTimeDerivative,PMLsetZero
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
USE MOD_Mesh_Vars,     ONLY: Elem_xGP,BCFace_xGP,Face_xGP  ! for PML region: xyz position of the Gauss points and Face Gauss points
USE MOD_Mesh_Vars,     ONLY: nSides
USE MOD_PML_Vars,      ONLY: PMLRampLength
USE MOD_PML_Vars,      ONLY: PMLzetaEff,PMLalpha
USE MOD_PML_Vars,      ONLY: nPMLVars
#ifdef PARTICLES
USE MOD_Interpolation_Vars,    ONLY:InterpolationInitIsDone
#endif
USE MOD_HDF5_output,           ONLY: GatheredWriteArray,GenerateFileSkeleton,WriteAttributeToHDF5,WriteHDF5Header
USE MOD_Mesh_Vars,             ONLY: MeshFile,nGlobalElems,offsetElem,SideToElem,ElemToSide
USE MOD_Output_Vars,           ONLY: ProjectName
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
#endif
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem,iProbe,iPMLElem,iFace,PP_nFaces,p,q,m
REAL                :: xyzMinMax(6),xi,L,XiN,delta(3),x,y,z
REAL                :: xyzMinMaxLoc(6)
REAL                :: zetaVec,zetaVecABS
INTEGER             :: multVec
INTEGER             :: iPMLFace,iPMLInterFace,nGlobalPMLElems,nGlobalPMLFaces,nGlobalPMLInterFaces
REAL                :: fLinear,fSinus,fPolynomial
REAL                :: dxMax,dxMin,dyMax,dyMin,dzMax,dzMin
REAL                :: FaceTolerance
INTEGER             :: DOFcount,ElemID1,ElemID2,iSide,ilocSide,SideID
! for HDF5 output
INTEGER(HID_T)                 :: Dset_ID
INTEGER                        :: nVal
CHARACTER(LEN=255)             :: FileName,FileString,MeshFile255,StrVarNames(PP_nVar),Statedummy
#ifdef MPI
REAL,ALLOCATABLE            :: PMLPlus(:,:,:,:)
REAL,ALLOCATABLE            :: PMLMinus(:,:,:,:)
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
! get information of PML size
PMLzeta0               = GETREAL('PMLzeta0','0.')
PMLalpha0              = GETREAL('PMLalpha0','0.')
xyzPhysicalMinMax(1:6) = GETREALARRAY('xyzPhysicalMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
IF(ALL(xyzPhysicalMinMax.EQ.0.))THEN
  xyzPhysicalMinMax(1:6)=HUGE(1.0) ! min
END IF
PMLzetaShape           = GETINT('PMLzetaShape','0')
PMLRampLength          = GETREAL('PMLRampLength','1.')
PMLspread              = GETINT('PMLspread','0')
PMLwriteFields         = GETINT('PMLwriteFields','0')
PMLzetaNorm            = GETLOGICAL('PMLzetaNorm','.FALSE.')
! caution, in current version read in in mesh
! only for Maxwell, PP_nVar=8

ALLOCATE(isPMLElem(1:PP_nElems))
isPMLElem=.FALSE.


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
  nPMLVars=0
  nPMLElems=0
  RETURN
ELSE
  nPMLVars=24
END IF

!===================================================================================================================================
! check if Element is PMLElem
!===================================================================================================================================
PP_nFaces=SIZE(Face_xGP(1,1,1,:))

ALLOCATE(isPMLFace(1:PP_nFaces))
ALLOCATE(isPMLInterFace(1:PP_nFaces))
isPMLFace=.FALSE.
isPMLInterFace=.FALSE.


isPMLElem=.TRUE.
DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! x-PML region
  IF (Elem_xGP(1,i,j,k,iElem) .GT. xyzPhysicalMinMax(1) .AND. Elem_xGP(1,i,j,k,iElem) .LT. xyzPhysicalMinMax(2)) THEN        
    isPMLElem(iElem) = .FALSE.
  END IF
  ! y-PML region
  IF (Elem_xGP(2,i,j,k,iElem) .GT. xyzPhysicalMinMax(3) .AND. Elem_xGP(2,i,j,k,iElem) .LT. xyzPhysicalMinMax(4)) THEN        
    isPMLElem(iElem) = .FALSE.
  END IF
  ! z-PML region
  IF (Elem_xGP(3,i,j,k,iElem) .GT. xyzPhysicalMinMax(5) .AND. Elem_xGP(3,i,j,k,iElem) .LT. xyzPhysicalMinMax(6)) THEN        
    isPMLElem(iElem) = .FALSE.
  END IF
END DO; END DO; END DO; END DO !iElem,k,i,j

DO iElem=1,PP_nElems
  IF(.NOT.isPMLElem(iElem))CYCLE
  DO ilocSide =1,6
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    isPMLFace(SideID)=.TRUE.
  END DO ! ilocSide=1,6
END DO ! iElem=1,PP_nElems

#ifdef MPI
ALLOCATE(PMLPlus(1,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(PMLMinus(1,0:PP_N,0:PP_N,1:nSides))
PMLPlus=0.
PMLMinus=0.
DO iSide=1,nSides
  IF(isPMLFace(iSide)) PMLPlus(1,:,:,iSide)=1.
END DO

! send my info to neighbor 
CALL StartReceiveMPIData(1,PMLMinus(1,0:PP_N,0:PP_N,SideID_plus_lower:SideID_plus_upper) &
                          ,SideID_plus_lower,SideID_plus_upper,RecRequest_U,SendID=2) ! Receive MINE
CALL StartSendMPIData(1,PMLPlus(1,0:PP_N,0:PP_N,SideID_plus_lower:SideID_plus_upper) &
                       ,SideID_plus_lower,SideID_plus_upper,SendRequest_U,SendID=2) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE

! add PMLMinus to PMLPlus and send
PMLPlus=PMLPlus+PMLMinus

CALL StartReceiveMPIData(1,PMLPlus,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
CALL StartSendMPIData(1,PMLPlus,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR

DO iSide=1,nSides
  IF(PMLPLUS(1,0,0,iSide).EQ.1)THEN
    iSPMLInterface(iSide)=.TRUE.
  END IF
END DO

DEALLOCATE(PMLPlus,PMLMinus)
#endif /*MPI*/

! and fill non-mpi sides
DO iSide=1,nSides
  IF(iSPMLInterface(iSide)) CYCLE
  IF(isPMLFace(iSide))THEN
    ElemID1=SideToElem(S2E_ELEM_ID,SideID)
    ElemID2=SideToElem(S2E_NB_ELEM_ID,SideID)
    IF(ElemID1.LT.1) CYCLE
    IF(ElemID2.LT.1) CYCLE
    IF((isPMLElem(ElemID1).AND. .NOT.isPMLElem(ElemID2)) .OR. &
       (.NOT.isPMLElem(ElemID1).AND.  isPMLElem(ElemID2)) )THEN
    !IF(isPMLElem(ElemID1).NOT..isPMLElem(ElemID2)) THEN
      iSPMLInterface(iSide)=.TRUE.
   END IF
  END IF
END DO


! Get number of PML Elems
nPMLFaces = 0
nPMLInterFaces = 0
nPMLElems = 0

DO iFace=1,PP_nFaces
  IF(isPMLFace(iFace))THEN
    nPMLFaces=nPMLFaces+1
  END IF
END DO ! iFace

DO iFace=1,PP_nFaces
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
    write(*,'(A92)') "      myrank              PP_nFaces              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              PP_nFaces         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  IF(.NOT.MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              PP_nFaces              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              PP_nFaces         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLInterFaces  ,nGlobalPMLInterFaces
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
    write(*,'(A92)') "      myrank              PP_nFaces              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              PP_nFaces         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLInterFaces  ,nGlobalPMLInterFaces
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
  IF(.NOT.MPIroot)THEN
    write(*,'(A92)') "      myrank              PP_nElems              nPMLElems        nGlobalPMLElems"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nElems             ,nPMLElems       ,nGlobalPMLElems
    write(*,'(A92)') "      myrank              PP_nFaces              nPMLFaces        nGlobalPMLFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLFaces       ,nGlobalPMLFaces
    write(*,'(A92)') "      myrank              PP_nFaces         nPMLInterFaces   nGlobalPMLInterFaces"
    write(*,'(I23,I23,I23,I23)')       myrank             ,PP_nFaces             ,nPMLInterFaces  ,nGlobalPMLInterFaces
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
ALLOCATE(FaceToPML(PP_nFaces)&
        ,PMLToFace(nPMLFaces))
ALLOCATE(FaceToPMLInter(PP_nFaces)&
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
DO iFace=1,PP_nFaces
  IF(isPMLFace(iFace))THEN
    iPMLFace=iPMLFace+1
    FaceToPML(iFace) = iPMLFace
    PMLToFace(iPMLFace) = iFace
  END IF
END DO
iPMLInterFace=0
DO iFace=1,PP_nFaces
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
ALLOCATE(U2       (1:24,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))        
ALLOCATE(U2t      (1:24,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
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



        DO m=1,3
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
  ! x-direction
  !IF (PMLzeta(1,i,j,k,iPMLElem).EQ.0.) THEN
    !PMLzetaEff(1,i,j,k,iPMLElem) = 0.
  !ELSE
    PMLzetaEff(1,i,j,k,iPMLElem) = ( PMLalpha(1,i,j,k,iPMLElem)+PMLzeta(1,i,j,k,iPMLElem) )
  !END IF
  ! y-direction
  !IF (PMLzeta(2,i,j,k,iPMLElem).EQ.0.) THEN
    !PMLzetaEff(2,i,j,k,iPMLElem) = 0.
  !ELSE
    PMLzetaEff(2,i,j,k,iPMLElem) = ( PMLalpha(2,i,j,k,iPMLElem)+PMLzeta(2,i,j,k,iPMLElem) )
  !END IF
  ! z-direction
  !IF (PMLzeta(3,i,j,k,iPMLElem).EQ.0.) THEN
    !PMLzetaEff(3,i,j,k,iPMLElem) = 0.
  !ELSE
    PMLzetaEff(3,i,j,k,iPMLElem) = ( PMLalpha(3,i,j,k,iPMLElem)+PMLzeta(3,i,j,k,iPMLElem) )
  !END IF
  !check NaN
  !if (ISNAN(PMLzetaEff(1,i,j,k,iPMLElem))) THEN
    !CALL abort(__STAMP__,&
    !'"PMLzetaEff(1,i,j,k,iPMLElem)" is a NaN',999,999.)
  !END IF
  !if (ISNAN(PMLzetaEff(2,i,j,k,iPMLElem))) THEN
    !CALL abort(__STAMP__,&
    !'"PMLzetaEff(2,i,j,k,iPMLElem)" is a NaN',999,999.)
  !END IF
  !if (ISNAN(PMLzetaEff(3,i,j,k,iPMLElem))) THEN
    !CALL abort(__STAMP__,&
    !'"PMLzetaEff(3,i,j,k,iPMLElem)" is a NaN',999,999.)
  !END IF
END DO; END DO; END DO; END DO !iPMLElem,k,i,j
DEALLOCATE(PMLalpha)




!===================================================================================================================================
! create global zeta field for parallel output of zeta distribution
!===================================================================================================================================
DO iElem=1,PP_nElems
  IF(isPMLElem(iElem))THEN
    PMLzetaGlobal(:,:,:,:,iElem)=PMLzeta(:,:,:,:,ElemToPML(iElem))/PMLzeta0
!print*,"PMLzetaGlobal",PMLzetaGlobal
  END IF
END DO!iElem



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
USE MOD_Globals,       ONLY: abort
USE MOD_Equation_Vars, ONLY: c_corr,c
USE MOD_DG_Vars,       ONLY: Ut,U
USE MOD_PML_Vars,      ONLY: nPMLElems,ElemToPML,PMLToElem
USE MOD_PML_Vars,      ONLY: nPMLFaces,FaceToPML,PMLToFace
USE MOD_PML_Vars,      ONLY: PMLzeta,U2
USE MOD_DG_Vars,       ONLY: D !for
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
!print*,"==========================================================================================================================="
!print *," CALL CalcPMLSource()"
!print*,"==========================================================================================================================="
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! sources for the standard variables
  DO m=1,8
    !print*,"m",m            ! = 1,2,3, 4, 5, 6, 7, 8
    !print*,"m*3-2",m*3-2    ! = 1,4,7,10,13,16,19,22
    !print*,"m*3-1",m*3-1    ! = 1,5,8,11,12,17,20,23
    !print*,"m*3",m*3        ! = 3,6,9,12,15,18,21,24
    Ut(m,i,j,k,PMLToElem(iPMLElem)) = Ut(m,i,j,k,PMLToElem(iPMLElem))  &
                                     -PMLzeta(1,i,j,k,iPMLElem)*U2(m*3-2,i,j,k,iPMLElem) &   ! = 1,4,7,10,13,16,19,22
                                     -PMLzeta(2,i,j,k,iPMLElem)*U2(m*3-1,i,j,k,iPMLElem) &   ! = 2,5,8,11,12,17,20,23
                                     -PMLzeta(3,i,j,k,iPMLElem)*U2(m*3  ,i,j,k,iPMLElem)     ! = 3,6,9,12,15,18,21,24
                                    !-DOT_PRODUCT(PMLzeta(:,i,j,k,iPMLElem),U2(m*3-2:m*3,i,j,k,iPMLElem))

  END DO


END DO; END DO; END DO !nPMLElems,k,j,i
    !print*,"Element:                         ",iPMLElem
    !print*,"MAXVAL(PMLzeta(:,:,:,:,iPMLElem))",MAXVAL(PMLzeta(1,:,:,:,iPMLElem)),MAXVAL(PMLzeta(2,:,:,:,iPMLElem)),&
!MAXVAL(PMLzeta(3,:,:,:,iPMLElem))
    !print*,"MINVAL(PMLzeta(:,:,:,:,iPMLElem))",MINVAL(PMLzeta(1,:,:,:,iPMLElem)),MINVAL(PMLzeta(2,:,:,:,iPMLElem)),&
!MINVAL(PMLzeta(3,:,:,:,iPMLElem))
    !print*,"MAXVAL(U2(:,:,:,:,iPMLElem))     ",MAXVAL(U2(:,:,:,:,iPMLElem))
END DO
END SUBROUTINE CalcPMLSource


SUBROUTINE PMLTimeDerivative()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY: abort
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY: U,Ut
USE MOD_PML_Vars,      ONLY: PMLzeta,U2,U2t,nPMLVars
USE MOD_PML_Vars,      ONLY: nPMLElems,ElemToPML,PMLToElem
USE MOD_PML_Vars,      ONLY: nPMLFaces,FaceToPML,PMLToFace
USE MOD_Mesh_Vars,     ONLY: Elem_xGP,sJ
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
  DO iPMLVar=1,24
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
U2(19:24,:,:,:,:) = fDamping* U2(19:24,:,:,:,:) 



END SUBROUTINE PMLTimeDerivative


SUBROUTINE PMLsetZero()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY: abort
USE MOD_PreProc
USE MOD_PML_Vars,      ONLY: nPMLElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPMLElem
!===================================================================================================================================
! null here to increase time for communication
DO iPMLElem=1,nPMLElems
  !FPt(1:3,:,:,:,iPMLElem) = 0.0 
  !FQt(1:3,:,:,:,iPMLElem) = 0.0 
  !FRt(1:3,:,:,:,iPMLElem) = 0.0 
  !FLt(1:3,:,:,:,iPMLElem) = 0.0 
  !FMt(1:3,:,:,:,iPMLElem) = 0.0 
  !FNt(1:3,:,:,:,iPMLElem) = 0.0 
  !FSt(1:3,:,:,:,iPMLElem) = 0.0 
  !FTt(1:3,:,:,:,iPMLElem) = 0.0 
END DO ! iPMLElem=1,nPMLElems
END SUBROUTINE PMLsetZero

SUBROUTINE Print_U(t)
!===============================================================================================================================
!  
!===============================================================================================================================
! MODULES
!USE MOD_PML,              ONLY: Print_U
USE MOD_PreProc
USE MOD_DG_Vars,          ONLY: U,Ut
USE MOD_MESH_VARS,         ONLY: sJ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: I
REAL,PARAMETER                :: infinity=HUGE(1.)
!===================================================================================================================================
print*,"MAXVAL(Ut(:,:,:,:,:))              ",MAXVAL(Ut(:,:,:,:,:))
do I=1,PP_nElems
!print*,"U = ",Ut(:,1,1,1,I)
IF (ISNAN(MAXVAL(Ut(:,1,1,1,I)))) THEN
  print*,"======================== NaN TERMINATION ========================"
  print*,"t=",t
  print*,"I=",I
  print*,"MAXVAL(Ut(:,1,1,1,I))",MAXVAL(Ut(:,1,1,1,I))
  stop
end if
IF (MAXVAL(Ut(:,1,1,1,I)).GE.infinity) THEN
  print*,"====================== INFINITY TERMINATION ======================"
  print*,"t=                                ",t
  print*,"I=1,PP_nElems                     ",I," (of ",PP_nElems,")"
  print*,"MAXVAL(sJ(:,:,:,I))              ",MAXVAL(sJ(:,:,:,I))
  print*,"MINVAL(sJ(:,:,:,I))              ",MINVAL(sJ(:,:,:,I))
  print*,"MAXVAL(Ut(:,1,1,1,I))              ",MAXVAL(Ut(:,1,1,1,I))
  print*,"(MAXVAL(Ut(:,1,1,1,I)).GE.infinity)",(MAXVAL(Ut(:,1,1,1,I)).GE.infinity)
  print*,"ISNAN(MAXVAL(Ut(:,1,1,1,I)))       ",ISNAN(MAXVAL(Ut(:,1,1,1,I)))
  stop
end if
end do
END SUBROUTINE Print_U



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
!print *, "PMLzeta size, shape and location",SIZE(PMLzeta),shape(PMLzeta),loc(PMLzeta)
!print *, "U2 size, shape and location     ",SIZE(U2),shape(U2),loc(U2)
IF(.NOT.DoPML) RETURN
!print*,"ALLOCATED(PMLzeta)",ALLOCATED(PMLzeta)
SDEALLOCATE(PMLzeta)
!print*,"ALLOCATED(U2)",ALLOCATED(U2)
SDEALLOCATE(U2)
!print*,"ALLOCATED(U2t)",ALLOCATED(U2t)
SDEALLOCATE(U2t) !momentan verursacht dies noch einen crash (resolve it !!)
!print*,"4"
SDEALLOCATE(PMLToElem)
!print*,"5"
SDEALLOCATE(ElemToPML)
!print*,"6"
SDEALLOCATE(PMLToFace)
!print*,"7"
SDEALLOCATE(FaceToPML)
SDEALLOCATE(PMLRamp)
SDEALLOCATE(isPMLElem)
!print*,"10"
SDEALLOCATE(isPMLFace)
!print*,"11"
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


