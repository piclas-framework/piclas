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
USE MOD_PML_Vars
USE MOD_HDF5_output,   ONLY: GatheredWriteArray,GenerateFileSkeleton,WriteAttributeToHDF5,WriteHDF5Header
USE MOD_HDF5_output,   ONLY: WritePMLzetaGlobalToHDF5
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER             :: i,j,k,iElem,iPMLElem,m
!REAL                :: xyzMinMax(6)!,xi,L,XiN,delta(3),x,y,z
!REAL                :: xyzMinMaxLoc(6)
!REAL                :: zetaVec,zetaVecABS
!INTEGER             :: nGlobalPMLElems,nGlobalPMLFaces,nGlobalPMLInterFaces
!INTEGER             :: DOFcount
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
xyzPMLMinMax(1:6)      = GETREALARRAY('xyzPMLMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
IF(ALMOSTEQUAL(MAXVAL(xyzPhysicalMinMax),MINVAL(xyzPhysicalMinMax)))THEN
  xyzPhysicalMinMax(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
  IF(ALMOSTEQUAL(MAXVAL(xyzPMLMinMax),MINVAL(xyzPMLMinMax)))THEN
    xyzPMLMinMax(1:4)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
    usePMLMinMax=.FALSE.
    SWRITE(UNIT_stdOut,'(A)')"no PML region supplied, setting xyzPhysicalMinMax =",xyzPhysicalMinMax
    SWRITE(UNIT_stdOut,'(A)')"no PML region supplied, setting xyzPMLMinMax      =",xyzPMLMinMax
  ELSE
    usePMLMinMax=.FALSE.
  END IF
END IF
PMLzetaShape           = GETINT('PMLzetaShape','0')
PMLRampLength          = GETREAL('PMLRampLength','1.')
PMLspread              = GETINT('PMLspread','0')
PMLwriteFields         = GETINT('PMLwriteFields','0')
PMLzetaNorm            = GETLOGICAL('PMLzetaNorm','.FALSE.')

PMLprintInfo           = GETINT('PMLprintInfo','0') ! 0=only root prints PML info, 1=all procs print PML info
IF(PMLprintInfo.EQ.0)THEN
  PMLprintInfoProcs=1 ! only root prints infos
ELSE
  PMLprintInfoProcs=nProcessors ! all procs print their infos
END IF
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

! find all elements in the PML region. Here: find all elements located outside of 'xyzPhysicalMinMax' 
CALL FindElementInRegion(isPMLElem,xyzPhysicalMinMax,ElementIsInside=.FALSE.)

! find all faces in the PML region
CALL FindInterfaces(isPMLFace,isPMLInterFace,isPMLElem)

! Get number of PML Elems, Faces and Interfaces. Create Mappngs PML <-> physical region
CALL CountAndCreateMappings('PML',&
                            isPMLElem,isPMLFace,isPMLInterFace,&
                            nPMLElems,nPMLFaces, nPMLInterFaces,&
                            ElemToPML,PMLToElem,&
                            FaceToPML,PMLToFace,&
                            FaceToPMLInter,PMLInterToFace)

! nPMLElems is determined, now allocate the PML field correnspondingly
ALLOCATE(U2       (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))        
ALLOCATE(U2t      (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
U2=0.
U2t=0.

! Set the damping profile function zeta=f(x,y) in the PML region
CALL SetPMLdampingProfile()

! create a HDF5 file containing the PMLzetaGlobal field
CALL WritePMLzetaGlobalToHDF5()

PMLInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT PML DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
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
USE MOD_Globals,       ONLY:abort,myrank,MPI_COMM_WORLD,UNIT_stdOut
USE MOD_Mesh_Vars,     ONLY:Elem_xGP
USE MOD_PML_Vars,      ONLY:PMLprintInfoProcs
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
#ifdef MPI
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
#endif /*MPI*/
  WRITE(UNIT_stdOut,'(A,6E15.6)')"Checking region:", region
#ifdef MPI
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
#endif /*MPI*/
IF(ElementIsInside)THEN
  isElem(:)=.TRUE.  !print*,"for elemens inside"
ELSE
  isElem(:)=.FALSE. !print*,"for elemens outside"
END IF

DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,3 ! m=x,y,z
    IF ( (Elem_xGP(m,i,j,k,iElem) .LT. region(2*m-1)) .OR. & ! 1,3,5
         (Elem_xGP(m,i,j,k,iElem) .GT. region(2*m)) ) THEN   ! 2,4,6 ! element is outside
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outisde the region
    END IF
  END DO
END DO; END DO; END DO; END DO !iElem,k,j,i


#ifdef MPI
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
#endif /*MPI*/
    IF(ElementIsInside)THEN
      WRITE(UNIT_stdOut,'(A,I12)')"No. of elements INSIDE region: ",COUNT(isElem)
    ELSE
      WRITE(UNIT_stdOut,'(A,I12)')"No. of elements OUTSIDE region: ",COUNT(isElem)
    END IF
#ifdef MPI
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
#endif /*MPI*/

END SUBROUTINE  FindElementInRegion


SUBROUTINE FindInterfaces(isFace,isInterFace,isElem)
!===================================================================================================================================
! Check if a face is in a special region (e.g. PML) and/or connects a special region (e.g. PML) to the physical region
!===================================================================================================================================
! MODULES
USE MOD_PreProc
!USE MOD_Globals,       ONLY: abort,myrank,MPI_COMM_WORLD!,nProcessors
USE MOD_Globals!,       ONLY: UNIT_stdOut,iError
USE MOD_Mesh_Vars,     ONLY: nSides,nBCSides
USE MOD_PML_vars,      ONLY: PMLprintInfoProcs
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
!USE MOD_Mesh_Vars,     ONLY:SideID_plus_upper,SideID_plus_lower
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: isElem(1:PP_nElems)
!CHARACTER(LEN=*),INTENT(IN)       :: TypeName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isFace(:)
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isInterFace(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides):: isFace_Plus,isFace_Minus,isFace_combined
INTEGER                                 :: iSide
INTEGER                                 :: I!,SideCounterUnity,SideCounterUnityGlobal
LOGICAL                                 :: printInfo
!===================================================================================================================================
ALLOCATE(isFace(1:nSides))
ALLOCATE(isInterFace(1:nSides))
isFace=.FALSE.
isInterFace=.FALSE.
printInfo=.FALSE.

! Check each element for being part of the, e.g. PML, region: set each of the 6 sides to be .TRUE.
!DO iElem=1,PP_nElems
  !IF(.NOT.isElem(iElem))CYCLE
  !DO ilocSide =1,6
    !SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    !isFace(SideID)=.TRUE.
  !END DO ! ilocSide=1,6
!END DO ! iElem=1,PP_nElems

! ---------------------------------------------
! For MPI sides send the info to all other procs
isFace_Plus=0.
isFace_Minus=0.
isFace_combined=0.
CALL ProlongToFace_PMLInfo(isElem,isFace_Minus,isFace_Plus,doMPISides=.FALSE.)
#ifdef MPI
CALL ProlongToFace_PMLInfo(isElem,isFace_Minus,isFace_Plus,doMPISides=.TRUE.)

! send my info to neighbor 
CALL StartReceiveMPIData(1,isFace_Plus,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Plus,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

CALL StartReceiveMPIData(1,isFace_Minus,1,nSides ,RecRequest_U,SendID=1) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Minus,1,nSides,SendRequest_U,SendID=1) ! Send YOUR

CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE -receive YOUR
CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U,SendID=1) !Send MINE -receive YOUR
#endif /*MPI*/

! add isFace_Minus to isFace_Plus and send
isFace_combined=2*isFace_Plus+isFace_Minus
! use numbering:    2*isFace_Plus+isFace_Minus  = 1: Minus side is PML
!                                                 2: Plus  side is PML
!                                                 3: both sides are PML sides
!                                                 0: normal face in physical region

DO iSide=1,nSides
  IF(isFace_combined(1,0,0,iSide).GT.0)THEN
    isFace(iSide)=.TRUE. ! mixed or pure PML face:  when my side is not PML but neighbor is PML
    IF((isFace_combined(1,0,0,iSide).EQ.1).OR.&
       (isFace_combined(1,0,0,iSide).EQ.2))THEN
        isInterFace(iSide)=.TRUE. ! set all mixed faces as InterFaces, exclude BCs later on
    END IF
  END IF
END DO
isInterFace(1:nBCSides)=.FALSE. ! BC sides cannot be interfaces!

!SideCounterUnity=0
!DO iSide=1,nSides
  !IF(isFace_combined(1,0,0,iSide).GT.0)THEN
    !SideCounterUnity=SideCounterUnity+1
  !END IF
!END DO
!print*,SideCounterUnity
!CALL MPI_REDUCE(SideCounterUnity,SideCounterUnityGlobal,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
!CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
!DO I=0,PMLprintInfoProcs-1
  !IF(I.EQ.myrank)THEN
    !WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
    !' myrank=',myrank,' Found ',SideCounterUnityGlobal,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
  !END IF
  !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
!END DO

! test
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
    DO iSide=nSides,nSides
      WRITE(UNIT_stdOut,'(A8,I5,A15,I5,A2,L5,A15,I5,A8,I5,A2,L5,A12,I5)')&
              "myrank=",myrank,&
       ": isInterFace(",iSide,")=",isInterFace(iSide),&
       " of total= ",COUNT(isInterFace),&
       " isFace(",iSide,")=",isFace(iSide),"  of total= ",COUNT(isFace)
    END DO
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO

!print*,"should be total 16 (PML-interfaces) and total 100 (PML-faces)"
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END SUBROUTINE  FindInterFaces


SUBROUTINE CountAndCreateMappings(TypeName,&
                                  isElem,isFace,isInterFace,&
                                  nElems,nFaces, nInterFaces,&
                                  ElemToX,XToElem,&
                                  FaceToX,XToFace,&
                                  FaceToXInter,XInterToFace)
!===================================================================================================================================
! 1.) Count the number of Elements, Faces and Interfaces of the PML/BGK/... region
! 2.) Create mappings from general element to PML/BGK/... element and vice versa
!                                  face    to PML/BGK/... face or interface and vice vesa
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,     ONLY: nSides
USE MOD_PML_vars,      ONLY: PMLprintInfoProcs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                :: isElem(:),isFace(:),isInterFace(:)
CHARACTER(LEN=*),INTENT(IN)       :: TypeName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)             :: nFaces,nInterFaces,nElems
INTEGER,ALLOCATABLE,INTENT(INOUT) :: ElemToX(:),XToElem(:),FaceToX(:),XToFace(:),FaceToXInter(:),XInterToFace(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iElem,iFace,nGlobalElems,nGlobalFaces,nGlobalInterFaces
INTEGER                           :: iXElem,iXFace,iXInterFace
INTEGER                           :: I
!===================================================================================================================================
! Get number of Elems
nFaces = 0
nInterFaces = 0
nElems = 0
DO iFace=1,nSides
  IF(isFace(iFace))THEN
    nFaces=nFaces+1
  END IF
END DO ! iFace
DO iFace=1,nSides
  IF(isInterFace(iFace))THEN
    nInterFaces=nInterFaces+1
  END IF
END DO ! iFace
DO iElem=1,PP_nElems
  IF(isElem(iElem))THEN
    nElems=nElems+1
  END IF
END DO ! iElem
!IF(1.EQ.2)THEN
!===================================================================================================================================
! print face number infos
!===================================================================================================================================
#ifdef MPI
nGlobalElems=0     
nGlobalFaces=0     
nGlobalInterfaces=0
IF(0.EQ.myrank) WRITE(UNIT_stdOut,'(A)') "========================================================================================="
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
    !write(*,'(A8,I5,A11,I5,A11,I5,A17,I5)')&
    !" myrank=",myrank," PP_nElems=",PP_nElems," nElems=",nElems," nGlobalElems=",nGlobalElems
    !write(*,'(A8,I5,A11,I5,A11,I5,A17,I5)')&
    !" myrank=",myrank," nSides=",nSides," nFaces=",nFaces," nGlobalFaces=",nGlobalFaces
    !write(*,'(A8,I5,A11,I5,A11,I5,A17,I5)')&
    !" myrank=",myrank," nSides=",nSides," nFaces=",nInterFaces," nGlobalFaces=",nGlobalInterFaces
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalInterFaces,' nGlobal',TRIM(TypeName),"-InterFaces inside of ",TRIM(TypeName),'-region.'
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
!=======================================================================================================================
! only send info to root
CALL MPI_REDUCE(nElems     ,nGlobalElems     ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
CALL MPI_REDUCE(nFaces     ,nGlobalFaces     ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
CALL MPI_REDUCE(nInterFaces,nGlobalInterFaces,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
IF(MPIroot) WRITE(UNIT_stdOut,'(A)') "============================================================================================"
IF(MPIroot) WRITE(UNIT_stdOut,'(A)') "CALL MPI_REDUCE(nElems,nGlobalElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)" !testing 
IF(MPIroot) WRITE(UNIT_stdOut,'(A)') "============================================================================================"
CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
!=======================================================================================================================
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalInterFaces,' nGlobal',TRIM(TypeName),"-InterFaces inside of ",TRIM(TypeName),'-region.'
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
IF(0.EQ.myrank) WRITE(UNIT_stdOut,'(A)') "========================================================================================="
#else
nGlobalElems=nElems
nGlobalFaces=nFaces
nGlobalInterFaces=nInterFaces
WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalInterFaces,' nGlobal',TRIM(TypeName),"-InterFaces inside of ",TRIM(TypeName),'-region.'
#endif /*MPI*/

!===================================================================================================================================
! create  mappings: element<->pml-element
!                      face<->pml-face
!                      face<->interface
!===================================================================================================================================
ALLOCATE(ElemToX(PP_nElems)&
        ,XToElem(nElems))
ALLOCATE(FaceToX(nSides)&
        ,XToFace(nFaces))
ALLOCATE(FaceToXInter(nSides)&
        ,XInterToFace(nInterFaces))
ElemToX=0
XToElem=0
FaceToX=0
XToFace=0
FaceToXInter=0
XInterToFace=0
! Create array with mapping
iXElem=0
DO iElem=1,PP_nElems
  IF(isElem(iElem))THEN
    iXElem=iXElem+1
    ElemToX(iElem) = iXElem
    XToElem(iXElem) = iElem
  END IF
END DO
iXFace=0
DO iFace=1,nSides
  IF(isFace(iFace))THEN
    iXFace=iXFace+1
    FaceToX(iFace) = iXFace
    XToFace(iXFace) = iFace
  END IF
END DO
iXInterFace=0
DO iFace=1,nSides
  IF(isInterFace(iFace))THEN
    iXInterFace=iXInterFace+1
    FaceToXInter(iFace) = iXInterFace
    XInterToFace(iXInterFace) = iFace
  END IF
END DO

END SUBROUTINE CountAndCreateMappings


SUBROUTINE SetPMLdampingProfile()
!===================================================================================================================================
! Determine the local PML damping factor in x,y and z-direction using a constant/linear/polynomial/... function
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,     ONLY: Elem_xGP,Face_xGP,nBCSides!,Face_xGP
USE MOD_PML_Vars,      ONLY: PMLzeta,PMLzetaEff,PMLalpha,usePMLMinMax,xyzPMLzetaShapeBase!,xyPMLMinMax,PMLRamp
USE MOD_PML_Vars,      ONLY: nPMLElems,PMLToElem,PMLprintInfoProcs
USE MOD_PML_Vars,      ONLY: PMLzeta0,PMLalpha0,xyzPhysicalMinMax,PMLzetaShape!,PMLRampLength!,PMLwriteFields
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iPMLElem
REAL                :: xyzMinMax(6),xi,XiN!,delta(2),x,y
REAL                :: xyzMinMaxloc(6)
REAL                :: function_type
INTEGER             :: DOFcount,iDir
REAL                :: L_vec(6)
!===================================================================================================================================
!ALLOCATE(PMLRamp          (0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLzeta      (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLzetaEff   (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PMLalpha     (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
PMLzeta=0.
!PMLRamp=1. ! goes from 1 to 0
PMLzetaEff=0.
PMLalpha=PMLalpha0 ! currently only constant a alpha distribution in the PML region is used
DOFcount=0
! get processor local bounding box of faces for damping value ramp
xyzMinMaxloc(:) = (/MINVAL(Face_xGP(1,:,:,1:nBCSides)),MAXVAL(Face_xGP(1,:,:,1:nBCSides)),&
                    MINVAL(Face_xGP(2,:,:,1:nBCSides)),MAXVAL(Face_xGP(2,:,:,1:nBCSides)),&
                    MINVAL(Face_xGP(3,:,:,1:nBCSides)),MAXVAL(Face_xGP(3,:,:,1:nBCSides))/)
! get global bounding box of faces for damping value ramp
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

#ifdef MPI
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
#endif /*MPI*/
  print*,"xyzMinMax - X",xyzMinMax(1:2)
  print*,"xyzMinMax - Y",xyzMinMax(3:4)
  print*,"xyzMinMax - Z",xyzMinMax(5:6)
#ifdef MPI
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
#endif /*MPI*/

#ifdef MPI
DO I=0,PMLprintInfoProcs-1
  IF(I.EQ.myrank)THEN
#endif /*MPI*/
  print*,"xyzPhysicalMinMax - X",xyzPhysicalMinMax(1:2)
  print*,"xyzPhysicalMinMax - Y",xyzPhysicalMinMax(3:4)
  print*,"xyzPhysicalMinMax - Z",xyzPhysicalMinMax(5:6)
#ifdef MPI
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
#endif /*MPI*/

! ----------------------------------------------------------------------------------------------------------------------------------
!determine PMLzeta values for each interpolation point according to ramping function (const., linear, sinusoidal, polynomial)
IF(usePMLMinMax)THEN ! use xyPMLMinMax -> define the PML region
  DO I=1,6
    IF( ALMOSTEQUAL(ABS(xyzMinMax(I)),ABS(xyzPMLzetaShapeBase(I))) )THEN
      L_vec(I)=HUGE(1.)
    ELSE 
      L_vec(I)=ABS(xyzMinMax(I)-xyzPMLzetaShapeBase(I))
    END IF
    print*,"L_vec   =", NINT(L_vec(I))
  END DO
  DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO iDir=1,3 !1=x, 2=y, 3=z
      IF       (Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)) .LT.xyzPMLzetaShapeBase(2*iDir-1)) THEN ! 1,3,5
        xi =ABS(Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem))    -xyzPMLzetaShapeBase(2*iDir-1))      ! 1,3,5
        PMLzeta(iDir,i,j,k,iPMLElem) = PMLzeta0*function_type(xi/L_vec(2*iDir-1),PMLzetashape)   ! 1,3,5
      ELSEIF   (Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)) .GT. xyzPMLzetaShapeBase(2*iDir )) THEN ! 2,4,6
        xi =ABS(Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem))     -xyzPMLzetaShapeBase(2*iDir ))      ! 2,4,6
        PMLzeta(iDir,i,j,k,iPMLElem) = PMLzeta0*function_type(xi/L_vec(2*iDir ),PMLzetashape)    ! 2,4,6
      END IF
    END DO
  END DO; END DO; END DO; END DO !iPMLElem,k,j,i
! ----------------------------------------------------------------------------------------------------------------------------------
ELSE ! use xyzPhysicalMinMax -> define the physical region
  !SELECT CASE (PMLzetaShape)
  !CASE(0) !Constant Distribution of the Damping Coefficient
    !DO iPMLElem=1,nPMLElems
          !!IF (isPMLElem(iElem)) THEN
            !PMLzeta( 1:3,:,:,:,iPMLElem) = PMLzeta0
            !PMLalpha(1:3,:,:,:,iPMLElem) = PMLalpha0
          !!END IF
    !END DO
  !CASE(1,2,3) ! Linear/Sinusoidal/POlynomial Distribution of the Damping Coefficient
    DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      DO iDir=1,3 !1=x, 2=y, 3=z
        IF          (Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)) .LT. xyzPhysicalMinMax(2*iDir-1)) THEN ! region is in lower part of the domain
          XiN = (ABS(Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)))-ABS(xyzPhysicalMinMax(2*iDir-1)))/&
                                                        (ABS(xyzMinMax(2*iDir-1))-ABS(xyzPhysicalMinMax(2*iDir-1)))
                      PMLzeta(iDir,i,j,k,iPMLElem)   = PMLzeta0*function_type(XiN,PMLzetaShape)
        ELSEIF      (Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)) .GT. xyzPhysicalMinMax(2*iDir)) THEN ! region is in upper part of the domain
          XiN = (ABS(Elem_xGP(iDir,i,j,k,PMLToElem(iPMLElem)))-ABS(xyzPhysicalMinMax(2*iDir)))/&
                                                        (ABS(xyzMinMax(2*iDir))-ABS(xyzPhysicalMinMax(2*iDir)))
                      PMLzeta(iDir,i,j,k,iPMLElem)   = PMLzeta0*function_type(XiN,PMLzetaShape)
        END IF
      END DO
    END DO; END DO; END DO; END DO !iElem,k,j,i
  !CASE DEFAULT
    !CALL abort(&
    !__STAMP__&
    !,'Shape function for damping coefficient in PML region not specified!',999,999.)
  !END SELECT ! PMLzetaShape





!    FIX this   ! determine Elem_xGP distance to PML interface for PMLRamp
!    FIX this   DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N
!    FIX this     ! x-PML region
!    FIX this     x = Elem_xGP(1,j,k,PMLToElem(iPMLElem))
!    FIX this     y = Elem_xGP(2,j,k,PMLToElem(iPMLElem))
!    FIX this     delta=0.
!    FIX this     ! x-PML region --------------------------------------------------------------
!    FIX this     IF (x .LT. xyzPhysicalMinMax(1)) THEN
!    FIX this       xi                  = ABS(x)          -ABS(xyzPhysicalMinMax(1))
!    FIX this       L                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
!    FIX this     ELSEIF (x .GT. xyzPhysicalMinMax(2)) THEN
!    FIX this       xi                  = ABS(x)          -ABS(xyzPhysicalMinMax(2))
!    FIX this       L                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
!    FIX this     ELSE
!    FIX this       xi=0
!    FIX this       L=1
!    FIX this     END IF
!    FIX this     delta(1)=MAXVAL((/0.,xi/L/))
!    FIX this     ! y-PML region --------------------------------------------------------------
!    FIX this     IF (y .LT. xyzPhysicalMinMax(3)) THEN
!    FIX this       xi                  = ABS(y)          -ABS(xyzPhysicalMinMax(3))
!    FIX this       L                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
!    FIX this     ELSEIF (y .GT. xyzPhysicalMinMax(4)) THEN
!    FIX this       xi                  = ABS(y)          -ABS(xyzPhysicalMinMax(4))
!    FIX this       L                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
!    FIX this     ELSE
!    FIX this       xi=0
!    FIX this       L=1
!    FIX this     END IF
!    FIX this     delta(2)=MAXVAL((/0.,xi/L/))
!    FIX this     ! set the ramp value from 1 down to 0: use the larged value of "delta"
!    FIX this     PMLRamp(j,k,iPMLElem) = 1. - function_type(MAXVAL(delta),PMLzetaShape)
!    FIX this   END DO; END DO; END DO !iPMLElem,k,j
END IF ! usePMLMinMax
! ----------------------------------------------------------------------------------------------------------------------------------
! CFS-PML formulation: calculate zeta eff using the complex frequency shift PMLalpha
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  PMLzetaEff(:,i,j,k,iPMLElem) = ( PMLalpha(:,i,j,k,iPMLElem)+PMLzeta(:,i,j,k,iPMLElem) )
END DO; END DO; END DO; END DO !iPMLElem,k,j,i
DEALLOCATE(PMLalpha)












! OLD!!!!!!!!!!!!!!!!!!
!===================================================================================================================================
! Modification to zeta values
!===================================================================================================================================
!PMLzetaNorm=.TRUE.
! Normalizing: recalculate zeta if multiple direction
!       IF (PMLzetaNorm) THEN
!         DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!               zetaVecABS=SQRT(PMLzeta(1,i,j,k,iPMLElem)**2 &
!                              +PMLzeta(2,i,j,k,iPMLElem)**2 &
!                              +PMLzeta(3,i,j,k,iPMLElem)**2 )
!               zetaVec=MAX(PMLzeta(1,i,j,k,iPMLElem),0.)
!               zetaVec=MAX(PMLzeta(2,i,j,k,iPMLElem),zetaVec)
!               zetaVec=MAX(PMLzeta(3,i,j,k,iPMLElem),zetaVec)
!               PMLzeta(:,i,j,k,iPMLElem) = PMLzeta(:,i,j,k,iPMLElem)/zetaVecABS*zetaVec
!         END DO; END DO; END DO; END DO !iPMLElem,k,i,j
!       END IF



!===================================================================================================================================
! determine Elem_xGP distance to PML interface for PMLRamp
!===================================================================================================================================
!         !DO iPMLElem=1,nPMLElems; DO p=0,PP_N; DO q=0,PP_N
!         DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!           ! x-PML region
!           !x = Face_xGP(1,p,q,PMLToFace(iPMLFace))
!           !y = Face_xGP(2,p,q,PMLToFace(iPMLFace))
!           !z = Face_xGP(3,p,q,PMLToFace(iPMLFace))
!           x = Elem_xGP(1,i,j,k,PMLToElem(iPMLElem))
!           y = Elem_xGP(2,i,j,k,PMLToElem(iPMLElem))
!           z = Elem_xGP(3,i,j,k,PMLToElem(iPMLElem))
!           delta=0.
!         
!           ! x-PML region
!           IF (x .LT. xyzPhysicalMinMax(1)) THEN
!             xi                  = ABS(x)-ABS(xyzPhysicalMinMax(1))
!             L                   = ABS(xyzMinMax(1))-ABS(xyzPhysicalMinMax(1))
!           ELSEIF (x .GT. xyzPhysicalMinMax(2)) THEN
!             xi                  = ABS(x)-ABS(xyzPhysicalMinMax(2))
!             L                   = ABS(xyzMinMax(2))-ABS(xyzPhysicalMinMax(2))
!           ELSE
!             xi=0
!             L=1
!           END IF
!           delta(1)=MAXVAL((/0.,xi/L/))
!           ! y-PML region
!           IF (y .LT. xyzPhysicalMinMax(3)) THEN
!             xi                  = ABS(y)-ABS(xyzPhysicalMinMax(3))
!             L                   = ABS(xyzMinMax(3))-ABS(xyzPhysicalMinMax(3))
!           ELSEIF (y .GT. xyzPhysicalMinMax(4)) THEN
!             xi                  = ABS(y)-ABS(xyzPhysicalMinMax(4))
!             L                   = ABS(xyzMinMax(4))-ABS(xyzPhysicalMinMax(4))
!           ELSE
!             xi=0
!             L=1
!           END IF
!           delta(2)=MAXVAL((/0.,xi/L/))
!           ! x-PML region
!           IF (z .LT. xyzPhysicalMinMax(5)) THEN
!             xi                  = ABS(z)-ABS(xyzPhysicalMinMax(5))
!             L                   = ABS(xyzMinMax(5))-ABS(xyzPhysicalMinMax(5))
!           ELSEIF (z .GT. xyzPhysicalMinMax(6)) THEN
!             xi                  = ABS(z)-ABS(xyzPhysicalMinMax(6))
!             L                   = ABS(xyzMinMax(6))-ABS(xyzPhysicalMinMax(6))
!           ELSE
!             xi=0
!             L=1
!           END IF
!           delta(3)=MAXVAL((/0.,xi/L/))
!           ! set the ramp value from 1 down to 0
!           !PMLRamp(p,q,iPMLFace)=1.-( MAXVAL(delta)-SIN(2*ACOS(-1.)*MAXVAL(delta))/(2*ACOS(-1.)) )
!           PMLRamp(i,j,k,iPMLElem) = 1. - fLinear(MAXVAL(delta))
!         
!           ! set the ramp value from 1 down to 0.82 (measured power loss)
!           ! add ramp from 0 to 0.82 (power drain 30GHz Gyrotron over 2mm PML)
!           !PMLRamp(i,j,k,iPMLElem) = PMLRamp(i,j,k,iPMLElem) + 0.82*fLinear(MAXVAL(delta))
!         !END DO; END DO; END DO !iFace,p,q
!         END DO; END DO; END DO; END DO !iPMLElem,k,i,j

END SUBROUTINE SetPMLdampingProfile


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


SUBROUTINE ProlongToFace_PMLInfo(isElem,isFace_Minus,isFace_Plus,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
!USE MOD_Interpolation_Vars, ONLY: L_Minus,L_Plus
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
!USE MOD_Mesh_Vars,          ONLY: SideID_minus_lower,SideID_minus_upper
!USE MOD_Mesh_Vars,          ONLY: SideID_plus_lower,SideID_plus_upper
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE 
LOGICAL,INTENT(IN)              :: isElem(1:PP_nElems) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: isFace_Minus(1,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT)              :: isFace_Plus( 1,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: i,ElemID(2),SideID,flip(2),LocSideID(2),firstSideID,lastSideID
!===================================================================================================================================
IF(doMPISides)THEN
  ! only YOUR MPI Sides are filled
  firstSideID = nBCSides+nInnerSides+nMPISides_MINE+1
  lastSideID  = firstSideID-1+nMPISides_YOUR 
  flip(1)     = -1
ELSE
  ! BCSides, InnerSides and MINE MPISides are filled
  firstSideID = 1
  lastSideID  = nBCSides+nInnerSides+nMPISides_MINE
  flip(1)     = 0
END IF
DO SideID=firstSideID,lastSideID
  ! master side, flip=0
  ElemID(1)    = SideToElem(S2E_ELEM_ID,SideID)  
  locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! neighbor side !ElemID,locSideID and flip =-1 if not existing
  ElemID(2)    = SideToElem(S2E_NB_ELEM_ID,SideID)
  locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  flip(2)      = SideToElem(S2E_FLIP,SideID)
  DO i=1,2 !first maste then slave side
    SELECT CASE(Flip(i))
      CASE(0) ! master side
        isFace_Minus(:,:,:,SideID)=MERGE(1,0,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
      CASE(1:4) ! slave side
        isFace_Plus( :,:,:,SideID)=MERGE(1,0,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
    END SELECT
  END DO !i=1,2, masterside & slave side 
END DO !SideID
END SUBROUTINE ProlongToFace_PMLInfo
END MODULE MOD_PML
























!===================================================================================================================================
! local SUBROUTINES and FUNCTIONS


REAL FUNCTION function_type(x,PMLzetaShape)
!===================================================================================================================================
! switch between different types of ramping functions for the calculation of the local zeta damping value field 
!===================================================================================================================================
! MODULES
USE MOD_Globals,       ONLY: abort
! IMPLICIT VARIABLE HANDLING 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,    INTENT(IN) :: x
INTEGER, INTENT(IN) :: PMLzetaShape ! linear, polynomial, const., sinusoidal ramping function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: fLinear,fSinus,fPolynomial
!===================================================================================================================================
SELECT CASE (PMLzetaShape)
CASE(0) !Constant Distribution of the Damping Coefficient
  function_type=1.
CASE(1) ! Linear Distribution of the Damping Coefficient
  function_type=fLinear(x)
CASE(2) ! Sinusoidal  Distribution of the Damping Coefficient
  function_type=fSinus(x)
CASE(3) ! polynomial
  function_type=fPolynomial(x)
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,'Shape function for damping coefficient in PML region not specified!',999,999.)
END SELECT ! PMLzetaShape

END FUNCTION function_type


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


REAL FUNCTION fSinus(x)
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


