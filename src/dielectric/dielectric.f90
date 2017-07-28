#include "boltzplatz.h"


MODULE MOD_Dielectric
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
INTERFACE InitDielectric
  MODULE PROCEDURE InitDielectric
END INTERFACE
INTERFACE FinalizeDielectric
  MODULE PROCEDURE FinalizeDielectric
END INTERFACE

PUBLIC::InitDielectric,FinalizeDielectric
!===================================================================================================================================
CONTAINS

SUBROUTINE InitDielectric()
!===================================================================================================================================
!  Initialize perfectly matched layer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: eps0,mu0
USE MOD_ReadInTools
USE MOD_Dielectric_Vars
!USE MOD_HDF5_output,   ONLY: GatheredWriteArray,GenerateFileSkeleton,WriteAttributeToHDF5,WriteHDF5Header
!USE MOD_HDF5_output,   ONLY: WritePMLzetaGlobalToHDF5
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: I
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Dielectric...'
!===================================================================================================================================
! Readin
!===================================================================================================================================
DoDielectric                     = GETLOGICAL('DoDielectric','.FALSE.')
DielectricEps0                   = GETREAL('DielectricEps0','0.')
DielectricMu0                    = GETREAL('DielectricMu0','0.')
xyzPhysicalMinMaxDielectric(1:6) = GETREALARRAY('xyzPhysicalMinMaxDielectric',6,'0.0,0.0,0.0,0.0,0.0,0.0')
xyzDielectricMinMax(1:6)         = GETREALARRAY('xyzDielectricMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
! use xyzPhysicalMinMax before xyzDielectricMinMax: 1.) check for xyzPhysicalMinMax 2.) check for xyzDielectricMinMax
IF(ALMOSTEQUAL(MAXVAL(xyzPhysicalMinMax),MINVAL(xyzPhysicalMinMax)))THEN ! if still the initialized values
  xyzPhysicalMinMax(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
  IF(ALMOSTEQUAL(MAXVAL(xyzDielectricMinMax),MINVAL(xyzDielectricMinMax)))THEN ! if still the initialized values
    xyzDielectricMinMax(1:4)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
    useDielectricMinMax=.FALSE. ! ! xyzPhysicalMinMax and xyzDielectricMinMax are undefined -> use HUGE for both
    SWRITE(UNIT_stdOut,'(A)')"no Dielectric region supplied, setting xyzPhysicalMinMax(1:6): Setting [+-HUGE]"
    SWRITE(UNIT_stdOut,'(A)')"no Dielectric region supplied, setting xyzDielectricMinMax(1:6)     : Setting [+-HUGE]"
  ELSE
    SWRITE(UNIT_stdOut,'(A)')"Dielectric region supplied via xyzDielectricMinMax(1:6)"
    useDielectricMinMax=.TRUE. ! xyzPhysicalMinMax is undefined but xyzDielectricMinMax is not
  END IF
ELSE
  SWRITE(UNIT_stdOut,'(A)')"Dielectric region supplied via xyzPhysicalMinMax(1:6)"
END IF
! display ranges of Dielectric region depending on useDielectricMinMax
SWRITE(UNIT_stdOut,'(A,L)') 'useDielectricMinMax=',useDielectricMinMax
IF(.NOT.useDielectricMinMax)THEN
  SWRITE(UNIT_stdOut,'(A)') '  Ranges for xyzPhysicalMinMax(1:6) are'
  SWRITE(UNIT_stdOut,'(A)') '       [        x-dir         ] [        y-dir         ] [         z-dir        ]'
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MIN'
  DO I=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzPhysicalMinMax(2*I-1)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MAX'
  DO I=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzPhysicalMinMax(2*I)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
ELSE
  SWRITE(UNIT_stdOut,'(A)') 'Ranges for xyzDielectricMinMax(1:6) are'
  SWRITE(UNIT_stdOut,'(A)') '       [        x-dir         ] [        y-dir         ] [         z-dir        ]'
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MIN'
  DO I=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzDielectricMinMax(2*I-1)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MAX'
  DO I=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzDielectricMinMax(2*I)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
END IF

!DielectriczetaShape           = GETINT('DielectriczetaShape','0')
!DielectricRampLength          = GETREAL('DielectricRampLength','1.')
!Dielectricspread              = GETINT('Dielectricspread','0')
!DielectricwriteFields         = GETINT('DielectricwriteFields','0')
!DielectriczetaNorm            = GETLOGICAL('DielectriczetaNorm','.FALSE.')

DielectricprintInfo           = GETINT('DielectricprintInfo','0') ! 0=only root prints Dielectric info
!                                                                 ! 1=all ranks print Dielectric info
IF(DielectricprintInfo.EQ.0)THEN
  DielectricprintInfoProcs=1 ! only root prints infos
ELSE
  DielectricprintInfoProcs=nProcessors ! all procs print their infos
END IF
! caution, in current version read in in mesh
! only for Maxwell, PP_nVar=8

IF(.NOT.DoDielectric) THEN
  SWRITE(UNIT_stdOut,'(A)') ' Dielectric region deactivated. '
  !DielectricnVar=0
  nDielectricElems=0
  RETURN
ELSE
  !DielectricnVar=24
END IF

! find all elements in the Dielectric region. Here: find all elements located outside of 'xyzPhysicalMinMax' 
CALL FindElementInRegion(isDielectricElem,xyzPhysicalMinMax,ElementIsInside=.FALSE.)

! find all faces in the Dielectric region
CALL FindInterfaces(isDielectricFace,isDielectricInterFace,isDielectricElem)

! Get number of Dielectric Elems, Faces and Interfaces. Create Mappngs Dielectric <-> physical region
CALL CountAndCreateMappings('Dielectric',&
                            isDielectricElem,isDielectricFace,isDielectricInterFace,&
                            nDielectricElems,nDielectricFaces, nDielectricInterFaces,&
                            ElemToDielectric,DielectricToElem,&
                            FaceToDielectric,DielectricToFace,&
                            FaceToDielectricInter,DielectricInterToFace)

! nDielectricElems is determined, now allocate the Dielectric field correnspondingly
!ALLOCATE(U2       (1:DielectricnVar,0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))        
!ALLOCATE(U2t      (1:DielectricnVar,0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))
!U2=0.
!U2t=0.

! Set the damping profile function zeta=f(x,y) in the Dielectric region
CALL SetDielectricdampingProfile()

! create a HDF5 file containing the DielectriczetaGlobal field
CALL WriteDielectriczetaGlobalToHDF5()

DielectricInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Dielectric DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDielectric


SUBROUTINE FindElementInRegion(isElem,region,ElementIsInside)
!===================================================================================================================================
! Check 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,       ONLY:abort,myrank,UNIT_stdOut
#ifdef MPI
USE MOD_Globals,       ONLY:MPI_COMM_WORLD
#endif /*MPI*/
USE MOD_Mesh_Vars,     ONLY:Elem_xGP
USE MOD_Dielectric_Vars,      ONLY:DielectricprintInfoProcs
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
! Dielectric elements are inside of the xyDielectricMinMax region
#ifdef MPI
DO I=0,DielectricprintInfoProcs-1
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
DO I=0,DielectricprintInfoProcs-1
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
! Check if a face is in a special region (e.g. Dielectric) and/or connects a special region (e.g. Dielectric) to the physical region
!===================================================================================================================================
! MODULES
USE MOD_PreProc
!USE MOD_Globals,       ONLY: abort,myrank,MPI_COMM_WORLD!,nProcessors
USE MOD_Globals!,       ONLY: UNIT_stdOut,iError
USE MOD_Mesh_Vars,     ONLY: nSides,nBCSides
USE MOD_Dielectric_vars,      ONLY: DielectricprintInfoProcs
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

! Check each element for being part of the, e.g. Dielectric, region: set each of the 6 sides to be .TRUE.
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
CALL ProlongToFace_DielectricInfo(isElem,isFace_Minus,isFace_Plus,doMPISides=.FALSE.)
#ifdef MPI
CALL ProlongToFace_DielectricInfo(isElem,isFace_Minus,isFace_Plus,doMPISides=.TRUE.)

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
! use numbering:    2*isFace_Plus+isFace_Minus  = 1: Minus side is Dielectric
!                                                 2: Plus  side is Dielectric
!                                                 3: both sides are Dielectric sides
!                                                 0: normal face in physical region

DO iSide=1,nSides
  IF(isFace_combined(1,0,0,iSide).GT.0)THEN
    isFace(iSide)=.TRUE. ! mixed or pure Dielectric face:  when my side is not Dielectric but neighbor is Dielectric
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
!DO I=0,DielectricprintInfoProcs-1
  !IF(I.EQ.myrank)THEN
    !WRITE(UNIT_stdOut,'(A8,I5,A,I10,A10,A10,A22,A10,A8)') &
    !' myrank=',myrank,' Found ',SideCounterUnityGlobal,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
  !END IF
  !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
!END DO

! test
DO I=0,DielectricprintInfoProcs-1
  IF(I.EQ.myrank)THEN
    DO iSide=nSides,nSides
      WRITE(UNIT_stdOut,'(A8,I5,A15,I5,A2,L5,A15,I5,A8,I5,A2,L5,A12,I5)')&
              "myrank=",myrank,&
       ": isInterFace(",iSide,")=",isInterFace(iSide),&
       " of total= ",COUNT(isInterFace),&
       " isFace(",iSide,")=",isFace(iSide),"  of total= ",COUNT(isFace)
    END DO
  END IF
#ifdef MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
#endif /*MPI*/
END DO

!print*,"should be total 16 (Dielectric-interfaces) and total 100 (Dielectric-faces)"
!CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END SUBROUTINE  FindInterFaces


SUBROUTINE CountAndCreateMappings(TypeName,&
                                  isElem,isFace,isInterFace,&
                                  nElems,nFaces, nInterFaces,&
                                  ElemToX,XToElem,&
                                  FaceToX,XToFace,&
                                  FaceToXInter,XInterToFace)
!===================================================================================================================================
! 1.) Count the number of Elements, Faces and Interfaces of the Dielectric/BGK/... region
! 2.) Create mappings from general element to Dielectric/BGK/... element and vice versa
!                                  face    to Dielectric/BGK/... face or interface and vice vesa
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,     ONLY: nSides
USE MOD_Dielectric_vars,      ONLY: DielectricprintInfoProcs
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
DO I=0,DielectricprintInfoProcs-1
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
DO I=0,DielectricprintInfoProcs-1
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


SUBROUTINE FinalizeDielectric()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_Dielectric_Vars,            ONLY: Dielectriczeta,U2,U2t
USE MOD_Dielectric_Vars,            ONLY: ElemToDielectric,DielectricToElem,DoDielectric,isDielectricElem,isDielectricFace,DielectricToFace,FaceToDielectric
USE MOD_Dielectric_Vars,            ONLY: DielectricRamp
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
IF(.NOT.DoDielectric) RETURN
SDEALLOCATE(Dielectriczeta)
SDEALLOCATE(U2)
SDEALLOCATE(U2t)
SDEALLOCATE(DielectricToElem)
SDEALLOCATE(ElemToDielectric)
SDEALLOCATE(DielectricToFace)
SDEALLOCATE(FaceToDielectric)
SDEALLOCATE(DielectricRamp)
SDEALLOCATE(isDielectricElem)
SDEALLOCATE(isDielectricFace)
END SUBROUTINE FinalizeDielectric


SUBROUTINE ProlongToFace_DielectricInfo(isElem,isFace_Minus,isFace_Plus,doMPISides)
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
END SUBROUTINE ProlongToFace_DielectricInfo
END MODULE MOD_Dielectric
