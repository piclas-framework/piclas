#include "boltzplatz.h"


MODULE MOD_Dielectric
!===================================================================================================================================
! Dielectric material handling in Maxwell's equations (HDG dielectric is done elsewhere)
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
USE MOD_ReadInTools
USE MOD_Dielectric_Vars
USE MOD_HDF5_output,     ONLY: WriteDielectricGlobalToHDF5
USE MOD_Equation_Vars,   ONLY: c_corr,c
USE MOD_Mesh_Vars,       ONLY: Elem_xGP
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: i
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Dielectric...'
!===================================================================================================================================
! Readin
!===================================================================================================================================
DoDielectric                     = GETLOGICAL('DoDielectric','.FALSE.')
DielectricEpsR                   = GETREAL('DielectricEpsR','1.')
DielectricMuR                    = GETREAL('DielectricMuR','1.')
DielectricTestCase               = GETSTR('DielectricTestCase','default')
DielectricRmax                   = GETREAL('DielectricRmax','1.')
IF((DielectricEpsR.LT.0.0).OR.(DielectricMuR.LT.0.0))THEN
  CALL abort(&
  __STAMP__&
  ,'Dielectric: MuR or EpsR cannot be negative.')
END IF
DielectricEpsR_inv               = 1./(DielectricEpsR)                   ! 1./EpsR
!DielectricConstant_inv           = 1./(DielectricEpsR*DielectricMuR)     !             1./(EpsR*MuR)
DielectricConstant_RootInv       = 1./sqrt(DielectricEpsR*DielectricMuR) !         1./sqrt(EpsR*MuR)
eta_c_dielectric                 = (c_corr-DielectricConstant_RootInv)*c ! ( chi - 1./sqrt(EpsR*MuR) ) * c
c_dielectric                     = c*DielectricConstant_RootInv          !          c/sqrt(EpsR*MuR)
c2_dielectric                    = c*c/(DielectricEpsR*DielectricMuR)            !           c**2/(EpsR*MuR)

xyzPhysicalMinMaxDielectric(1:6) = GETREALARRAY('xyzPhysicalMinMaxDielectric',6,'0.0,0.0,0.0,0.0,0.0,0.0')
xyzDielectricMinMax(1:6)         = GETREALARRAY('xyzDielectricMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
! use xyzPhysicalMinMaxDielectric before xyzDielectricMinMax: 1.) check for xyzPhysicalMinMaxDielectric 2.) check for xyzDielectricMinMax
IF(ALMOSTEQUAL(MAXVAL(xyzPhysicalMinMaxDielectric),MINVAL(xyzPhysicalMinMaxDielectric)))THEN ! if still the initialized values
  xyzPhysicalMinMaxDielectric(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
  IF(ALMOSTEQUAL(MAXVAL(xyzDielectricMinMax),MINVAL(xyzDielectricMinMax)))THEN ! if still the initialized values
    xyzDielectricMinMax(1:4)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
    useDielectricMinMax=.FALSE. ! ! xyzPhysicalMinMaxDielectric and xyzDielectricMinMax are undefined -> use HUGE for both
    SWRITE(UNIT_stdOut,'(A)')"no Dielectric region supplied, setting xyzPhysicalMinMaxDielectric(1:6): Setting [+-HUGE]"
    SWRITE(UNIT_stdOut,'(A)')"no Dielectric region supplied, setting xyzDielectricMinMax(1:6)     : Setting [+-HUGE]"
  ELSE
    SWRITE(UNIT_stdOut,'(A)')"Dielectric region supplied via xyzDielectricMinMax(1:6)"
    useDielectricMinMax=.TRUE. ! xyzPhysicalMinMaxDielectric is undefined but xyzDielectricMinMax is not
  END IF
ELSE
  SWRITE(UNIT_stdOut,'(A)')"Dielectric region supplied via xyzPhysicalMinMaxDielectric(1:6)"
END IF
! display ranges of Dielectric region depending on useDielectricMinMax
SWRITE(UNIT_stdOut,'(A,L)') 'useDielectricMinMax=',useDielectricMinMax
IF(.NOT.useDielectricMinMax)THEN
  SWRITE(UNIT_stdOut,'(A)') '  Ranges for xyzPhysicalMinMaxDielectric(1:6) are'
  SWRITE(UNIT_stdOut,'(A)') '       [        x-dir         ] [        y-dir         ] [         z-dir        ]'
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MIN'
  DO i=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzPhysicalMinMaxDielectric(2*i-1)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MAX'
  DO i=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzPhysicalMinMaxDielectric(2*i)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
ELSE
  SWRITE(UNIT_stdOut,'(A)') 'Ranges for xyzDielectricMinMax(1:6) are'
  SWRITE(UNIT_stdOut,'(A)') '       [        x-dir         ] [        y-dir         ] [         z-dir        ]'
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MIN'
  DO i=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzDielectricMinMax(2*i-1)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MAX'
  DO i=1,3
    SWRITE(UNIT_stdOut,OUTPUTFORMAT,ADVANCE='NO')  xyzDielectricMinMax(2*i)
  END DO
  SWRITE(UNIT_stdOut,'(A)') ''
END IF

DielectricprintInfo           = GETINT('DielectricprintInfo','0') ! 0=only root prints Dielectric info
!                                                                 ! 1=all ranks print Dielectric info
IF(DielectricprintInfo.EQ.0)THEN
  DielectricprintInfoProcs=1 ! only root prints infos
ELSE
  DielectricprintInfoProcs=nProcessors ! all procs print their infos
END IF

IF(.NOT.DoDielectric) THEN
  SWRITE(UNIT_stdOut,'(A)') ' Dielectric region deactivated. '
  nDielectricElems=0
  RETURN
END IF

! find all elements in the Dielectric region. Here: find all elements located outside of 'xyzPhysicalMinMaxDielectric' 
IF(useDielectricMinMax)THEN
  CALL FindElementInRegion(isDielectricElem,xyzDielectricMinMax,ElementIsInside=.TRUE.)
ELSE
  CALL FindElementInRegion(isDielectricElem,xyzPhysicalMinMaxDielectric,ElementIsInside=.FALSE.)
END IF

! find all faces in the Dielectric region
CALL FindInterfaces(isDielectricFace,isDielectricInterFace,isDielectricElem)

! Get number of Dielectric Elems, Faces and Interfaces. Create Mappngs Dielectric <-> physical region
CALL CountAndCreateMappings('Dielectric',&
                            isDielectricElem,isDielectricFace,isDielectricInterFace,&
                            nDielectricElems,nDielectricFaces, nDielectricInterFaces,&
                            ElemToDielectric,DielectricToElem,&
                            FaceToDielectric,DielectricToFace,&
                            FaceToDielectricInter,DielectricInterToFace)

! Set the dielectric profile function EpsR,MuR=f(x,y,z) in the Dielectric region
CALL SetDielectricVolumeProfile()

! Determine dielectric Values on faces and communicate them
CALL SetDielectricFaceProfile()

! create a HDF5 file containing the DielectriczetaGlobal field
CALL WriteDielectricGlobalToHDF5()

DielectricInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Dielectric DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDielectric


SUBROUTINE FindElementInRegion(isElem,region,ElementIsInside)
!===================================================================================================================================
! As soon as only one DOF is not inside/outside of the region, the complete element is excluded 
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
USE MOD_Globals
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides
USE MOD_Dielectric_vars, ONLY: DielectricprintInfoProcs
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,             ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: isElem(1:PP_nElems) ! True/False element: special region
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isFace(:)           ! True/False face: special region <-> special region
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isInterFace(:)      ! True/False face: special region <-> physical region (or vice versa)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides):: isFace_Slave,isFace_Master,isFace_combined ! the dimension is only used because of
                                                                                    ! the prolong to face routine and MPI logic
INTEGER                                 :: iSide
INTEGER                                 :: I
!===================================================================================================================================
ALLOCATE(isFace(1:nSides))
ALLOCATE(isInterFace(1:nSides))
isFace=.FALSE.
isInterFace=.FALSE.

! ---------------------------------------------
! For MPI sides send the info to all other procs
isFace_Slave=0.
isFace_Master=0.
isFace_combined=0.
CALL ProlongToFace_DielectricInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.FALSE.)
#ifdef MPI
CALL ProlongToFace_DielectricInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.TRUE.)

! send my info to neighbor 
CALL StartReceiveMPIData(1,isFace_Slave,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Slave,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

CALL StartReceiveMPIData(1,isFace_Master,1,nSides ,RecRequest_U,SendID=1) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Master,1,nSides,SendRequest_U,SendID=1) ! Send YOUR

CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE -receive YOUR
CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U,SendID=1) !Send MINE -receive YOUR
#endif /*MPI*/

! add isFace_Master to isFace_Slave and send
isFace_combined=2*isFace_Slave+isFace_Master
! use numbering:  2*isFace_Slave+isFace_Master  = 1: Master side is Dielectric
!                                                 2: Slave  side is Dielectric
!                                                 3: both sides are Dielectric sides
!                                                 0: normal face in physical region (no dielectric involved)

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


SUBROUTINE SetDielectricVolumeProfile()
!===================================================================================================================================
! Determine the local Dielectric damping factor in x,y and z-direction using a constant/linear/polynomial/... function
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,            ONLY: Elem_xGP
USE MOD_Dielectric_Vars,      ONLY: DielectricEps,DielectricMu,DielectricConstant_inv
USE MOD_Dielectric_Vars,      ONLY: nDielectricElems,DielectricToElem
USE MOD_Dielectric_Vars,      ONLY: DielectricRmax,DielectricEpsR,DielectricMuR,DielectricTestCase
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iDielectricElem
REAL                :: r
!===================================================================================================================================
ALLOCATE(         DielectricEps(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))
ALLOCATE(          DielectricMu(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))
ALLOCATE(DielectricConstant_inv(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))
DielectricEps=0.
DielectricMu=0.
DielectricConstant_inv=0.
IF(TRIM(DielectricTestCase).EQ.'FishEyeLens')THEN
  ! use function with radial dependence: EpsR=n0^2 / (1 + (r/r_max)^2)^2
  DO iDielectricElem=1,nDielectricElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r = SQRT(Elem_xGP(1,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(2,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(3,i,j,k,DielectricToElem(iDielectricElem))**2  )
    DielectricEps(i,j,k,iDielectricElem) = 4./((1+(r/DielectricRmax)**2)**2)
  END DO; END DO; END DO; END DO !iDielectricElem,k,j,i
  DielectricMu(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems) = DielectricMuR
ELSE
  DielectricEps(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems) = DielectricEpsR
  DielectricMu( 0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems) = DielectricMuR
END IF
DielectricConstant_inv(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems) = 1./& ! 1./(EpsR*MuR)
                                                                 (DielectricEps(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems)*&
                                                                  DielectricMu( 0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems))
END SUBROUTINE SetDielectricVolumeProfile


SUBROUTINE SetDielectricFaceProfile()
!===================================================================================================================================
! set the dielectric factor 1./SQRT(EpsR*MuR) for each face DOF in the array "Dielectric_Master"
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars, ONLY:DielectricConstant_inv,dielectric_Master,dielectric_Slave,isDielectricElem,ElemToDielectric
USE MOD_Mesh_Vars,       ONLY:nSides
USE MOD_ProlongToFace,   ONLY:ProlongToFace
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,             ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: Dielectric_dummy_Master ! deallocate at the end of this routine
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: Dielectric_dummy_Slave  ! deallocate at the end of this routine
REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: Dielectric_dummy_elem   ! deallocate at the end of this routine
INTEGER                               :: iElem
!===================================================================================================================================
ALLOCATE(Dielectric_dummy_elem(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
Dielectric_dummy_elem=0. ! default is an invalid number
ALLOCATE(Dielectric_dummy_Master(PP_nVar,0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Dielectric_dummy_Slave( PP_nVar,0:PP_N,0:PP_N,1:nSides))
Dielectric_dummy_Master=0.
Dielectric_dummy_Slave =0.

! fill dummy values for non-Dielectric sides
DO iElem=1,PP_nElems
  IF(isDielectricElem(iElem))THEN
    ! set only the first dimension to 1./SQRT(EpsR*MuR) (the rest are dummies)
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=SQRT(DielectricConstant_inv(0:PP_N,0:PP_N,0:PP_N,ElemToDielectric(iElem)))
  ELSE
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=1.0
  END IF
END DO

CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.FALSE.)
#ifdef MPI
CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.TRUE.)

! send Slave Dielectric info (real array with dimension (N+1)*(N+1)) to Master procs
CALL StartReceiveMPIData(1,Dielectric_dummy_Slave(1,0:PP_N,0:PP_N,1:nSides) ,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,Dielectric_dummy_Slave(1,0:PP_N,0:PP_N,1:nSides) ,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

! send Master Dielectric info (real array with dimension (N+1)*(N+1)) to Slave procs
CALL StartReceiveMPIData(1,Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides),1,nSides ,RecRequest_U ,SendID=1) ! Receive YOUR
CALL StartSendMPIData(   1,Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides),1,nSides,SendRequest_U ,SendID=1) ! Send MINE

CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE - receive YOUR
CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U ,SendID=1) !Send YOUR - receive MINE 
#endif /*MPI*/

ALLOCATE(Dielectric_Master(0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Dielectric_Slave( 0:PP_N,0:PP_N,1:nSides))
Dielectric_Master=Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides)
Dielectric_Slave =Dielectric_dummy_Slave( 1,0:PP_N,0:PP_N,1:nSides)

SDEALLOCATE(Dielectric_dummy_Slave)
SDEALLOCATE(Dielectric_dummy_Master)
SDEALLOCATE(Dielectric_dummy_elem)
END SUBROUTINE SetDielectricFaceProfile


SUBROUTINE FinalizeDielectric()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_Dielectric_Vars,            ONLY: DoDielectric,DielectricEps,DielectricMu
USE MOD_Dielectric_Vars,            ONLY: ElemToDielectric,DielectricToElem,isDielectricElem
USE MOD_Dielectric_Vars,            ONLY: FaceToDielectric,DielectricToFace,isDielectricFace
!USE MOD_Dielectric_Vars,            ONLY: DielectricRamp
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
SDEALLOCATE(DielectricEps)
SDEALLOCATE(DielectricMu)
SDEALLOCATE(DielectricToElem)
SDEALLOCATE(ElemToDielectric)
SDEALLOCATE(DielectricToFace)
SDEALLOCATE(FaceToDielectric)
SDEALLOCATE(isDielectricElem)
SDEALLOCATE(isDielectricFace)
END SUBROUTINE FinalizeDielectric


SUBROUTINE ProlongToFace_DielectricInfo(isElem,isFace_Master,isFace_Slave,doMPISides)
!===================================================================================================================================
! Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
! integration points, using fast 1D Interpolation and store in global side structure
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE 
LOGICAL,INTENT(IN)              :: isElem(1:PP_nElems) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: isFace_Master(1,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT)              :: isFace_Slave( 1,0:PP_N,0:PP_N,1:nSides)
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
        isFace_Master(:,:,:,SideID)=MERGE(1,0,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
      CASE(1:4) ! slave side
        isFace_Slave( :,:,:,SideID)=MERGE(1,0,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
    END SELECT
  END DO !i=1,2, masterside & slave side 
END DO !SideID
END SUBROUTINE ProlongToFace_DielectricInfo


END MODULE MOD_Dielectric
