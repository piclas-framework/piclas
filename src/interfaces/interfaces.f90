#include "boltzplatz.h"

MODULE MOD_Interfaces
!===================================================================================================================================
!> Contains the routines to
!> - identify interfaces between special regions, e.g., PML <-> physical region
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitInterfaces
  MODULE PROCEDURE InitInterfaces
END INTERFACE
INTERFACE FindInterfacesInRegion
  MODULE PROCEDURE FindInterfacesInRegion
END INTERFACE
INTERFACE FindElementInRegion
  MODULE PROCEDURE FindElementInRegion
END INTERFACE
INTERFACE CountAndCreateMappings
  MODULE PROCEDURE CountAndCreateMappings
END INTERFACE
INTERFACE FinalizeInterfaces
  MODULE PROCEDURE FinalizeInterfaces
END INTERFACE

PUBLIC::InitInterfaces
PUBLIC::FindElementInRegion
PUBLIC::FindInterfacesInRegion
PUBLIC::CountAndCreateMappings
PUBLIC::FinalizeInterfaces
!===================================================================================================================================
CONTAINS


SUBROUTINE InitInterfaces
!===================================================================================================================================
!> Check every face and set the correct identifier for selecting the corresponding Riemann solver
!> possible connections are (Master <-> Slave direction is important):
!>   - vaccuum <-> vacuum          : RIEMANN_VACUUM         = 0
!>   - PML <-> vacuum              : RIEMANN_PML            = 1
!>   - PML <-> PML                 : RIEMANN_PML            = 1
!>   - dielectric <-> dielectric   : RIEMANN_DIELECTRIC     = 2
!>   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC = 3
!>   - vacuum      -> dielectri    : RIEMANN_VAC2DIELECTRIC = 4
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,       ONLY:nSides
#ifndef PP_HDG
USE MOD_PML_vars,        ONLY:DoPML,isPMLFace,isPMLInterFace
#endif /*NOT HDG*/
USE MOD_Dielectric_vars, ONLY:DoDielectric,isDielectricFace,isDielectricInterFace,isDielectricElem
USE MOD_Interfaces_Vars, ONLY:InterfaceRiemann,InterfacesInitIsDone
USE MOD_Globals,         ONLY:abort,myrank,UNIT_stdOut,mpiroot,iError
USE MOD_Mesh_Vars,       ONLY:SideToElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,ElemID
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT INTERFACES...' 
ALLOCATE(InterfaceRiemann(1:nSides))
DO SideID=1,nSides
  InterfaceRiemann(SideID)=-1 ! set default to invalid number: check later
  ! 0.) Sanity: It is forbidden to connect a PML to a dielectric region because it is not implemented!
#ifndef PP_HDG ! pure Maxwell simulations
  IF(DoPML.AND.DoDielectric)THEN
    IF(isPMLFace(SideID).AND.isDielectricFace(SideID))THEN
      CALL abort(&
      __STAMP__&
      ,'It is forbidden to connect a PML to a dielectric region! (Not implemented)')
    END IF
  END IF

  ! 1.) Check Perfectly Matched Layer
  IF(DoPML) THEN
    IF (isPMLFace(SideID))THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
                                 !     (flux-splitting!)
      InterfaceRiemann(SideID)=RIEMANN_PML
      CYCLE ! don't check the following if the flux has already been calculated here -> continue with next side
    END IF
  END IF ! DoPML
#endif /*NOT HDG*/

  ! 2.) Check Dielectric Medium
  IF(DoDielectric) THEN
    IF (isDielectricFace(SideID))THEN ! 1.) RiemannDielectric
      IF(isDielectricInterFace(SideID))THEN
        ! a) physical <-> dielectric region: for Riemann solver, select A+ and A- as functions of f(Eps0,Mu0) or f(EpsR,MuR)
        ElemID = SideToElem(S2E_ELEM_ID,SideID) ! get master element ID for checking if it is in a physical or dielectric region
        IF(isDielectricElem(ElemID))THEN !  master is DIELECTRIC and slave PHYSICAL
          InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC2VAC ! A+(Eps0,Mu0) and A-(EpsR,MuR)
        ELSE ! master is PHYSICAL and slave DIELECTRIC
          InterfaceRiemann(SideID)=RIEMANN_VAC2DIELECTRIC ! A+(EpsR,MuR) and A-(Eps0,Mu0)
        END IF
      ELSE ! b) dielectric region <-> dielectric region
        InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC
      END IF
    ELSE ! c) no Dielectric, standard flux
      InterfaceRiemann(SideID)=RIEMANN_VACUUM
    END IF ! IF(isDielectricFace(SideID))
  ELSE ! d) no Dielectric, standard flux
    InterfaceRiemann(SideID)=RIEMANN_VACUUM
  END IF ! DoDielectric

  IF(InterfaceRiemann(SideID).EQ.-1)THEN
    CALL abort(&
    __STAMP__&
    ,'Interface for Riemann solver not correctly determined (vacuum, dielectric, PML ...)')
  END IF
END DO ! SideID

InterfacesInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT INTERFACES DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitInterfaces


SUBROUTINE FindElementInRegion(isElem,region,ElementIsInside,DisplayInfoProcs)
!===================================================================================================================================
!> Determine whether an element resides within or outside of a special region (e.g. PML or dielectric region)
!> Note: As soon as only one DOF is not inside/outside of the region, the complete element is excluded 
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals,              ONLY:abort,myrank,UNIT_stdOut,mpiroot,iError
#ifdef MPI
USE MOD_Globals,              ONLY:MPI_COMM_WORLD
#endif /*MPI*/
USE MOD_Mesh_Vars,            ONLY:Elem_xGP
USE MOD_Dielectric_Vars,      ONLY:DielectricprintInfoProcs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: ElementIsInside
REAL,INTENT(IN)                  :: region(1:6)
INTEGER,INTENT(IN),OPTIONAL      :: DisplayInfoProcs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isElem(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,m,NProcs
!===================================================================================================================================
! set logical vector for each element
ALLOCATE(isElem(1:PP_nElems))
isElem=.FALSE.
IF(ElementIsInside)THEN
  SWRITE(UNIT_stdOut,'(A,6E15.6)')"Checking for elements INSIDE region:", region
  isElem(:)=.TRUE.
ELSE
  SWRITE(UNIT_stdOut,'(A,6E15.6)')"Checking for elements OUTSIDE region:", region
  isElem(:)=.FALSE.
END IF

DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,3 ! m=x,y,z
    IF ( (Elem_xGP(m,i,j,k,iElem) .LT. region(2*m-1)) .OR. & ! 1,3,5
         (Elem_xGP(m,i,j,k,iElem) .GT. region(2*m)) ) THEN   ! 2,4,6 ! element is outside
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outisde the region
    END IF
  END DO
END DO; END DO; END DO; END DO !iElem,k,j,i

! Display debugging output by each rank
IF(PRESENT(DisplayInfoProcs))THEN
  NProcs=DisplayInfoProcs
ELSE
  NProcs=0
END IF
#ifdef MPI
DO I=0,NProcs-1
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


SUBROUTINE FindInterfacesInRegion(isFace,isInterFace,isElem,DisplayInfoProcs)
!===================================================================================================================================
!> Check if a face is in a special region (e.g. Dielectric) and/or connects a special region (e.g. Dielectric) to the physical
!> region. This is used, e.g., for dielectric or PML regions
!> indentifies the following connections and stores them in "isFace_combined"
!> use numbering: 1: Master side is special (e.g. dielectric)
!>                2: Slave  side is special (e.g. dielectric)
!>                3: both sides are special (e.g. dielectric) sides
!>                0: normal face in physical region (no special region involved)
!> ToDo: adjust (reduce) the dimension of "isFace_combined", which currently corresponds to the dimensions of the MPI exchange
!> routines -> need for template fortran files
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,       ONLY: nSides,nBCSides
#ifdef MPI
USE MOD_MPI_Vars
USE MOD_MPI,             ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)               :: isElem(1:PP_nElems) ! True/False element: special region
INTEGER,INTENT(IN),OPTIONAL      :: DisplayInfoProcs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isFace(:)           ! True/False face: special region <-> special region
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isInterFace(:)      ! True/False face: special region <-> physical region (or vice versa)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides):: isFace_Slave,isFace_Master,isFace_combined ! the dimension is only used because of
                                                                                      ! the prolong to face routine and MPI logic
INTEGER                                 :: iSide,I,NProcs
!===================================================================================================================================
IF(PRESENT(DisplayInfoProcs))THEN
  NProcs=DisplayInfoProcs
ELSE
  NProcs=0
END IF
ALLOCATE(isFace(1:nSides))
ALLOCATE(isInterFace(1:nSides))
isFace=.FALSE.
isInterFace=.FALSE.

! For MPI sides send the info to all other procs
isFace_Slave=0.
isFace_Master=0.
isFace_combined=0.
CALL ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.FALSE.)
#ifdef MPI
CALL ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.TRUE.)

! send Slave special region info (real with [0=no special region] or [1=special region] as (N+1)*(N+1) array) to Master procs
CALL StartReceiveMPIData(1,isFace_Slave,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Slave,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

! send Master special region info (real with [0=no special region] or [1=special region] as (N+1)*(N+1) array) to Slave procs
CALL StartReceiveMPIData(1,isFace_Master,1,nSides ,RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIData(   1,isFace_Master,1,nSides,SendRequest_U,SendID=1) ! Send MINE

CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE -receive YOUR
CALL FinishExchangeMPIData(SendRequest_U ,RecRequest_U ,SendID=1) !Send YOUR -receive MINE
#endif /*MPI*/

! add isFace_Master to isFace_Slave and send
isFace_combined=2*isFace_Slave+isFace_Master
! use numbering:  2*isFace_Slave+isFace_Master  = 1: Master side is special (e.g. dielectric)
!                                                 2: Slave  side is special (e.g. dielectric)
!                                                 3: both sides are special (e.g. dielectric) sides
!                                                 0: normal face in physical region (no special region involved)

DO iSide=1,nSides
  IF(isFace_combined(1,0,0,iSide).GT.0)THEN
    isFace(iSide)=.TRUE. ! mixed or pure special region face: when my side is not special but neighbor is special
    IF((isFace_combined(1,0,0,iSide).EQ.1).OR.&
       (isFace_combined(1,0,0,iSide).EQ.2))THEN
        isInterFace(iSide)=.TRUE. ! set all mixed faces as InterFaces, exclude BCs later on
    END IF
  END IF
END DO
isInterFace(1:nBCSides)=.FALSE. ! BC sides cannot be interfaces!

! Debugging output (optional)
DO I=0,NProcs-1
  IF(I.EQ.myrank)THEN
    DO iSide=nSides,nSides
      WRITE(UNIT_stdOut,'(A8,I10,A15,I10,A2,L5,A15,I10,A8,I10,A2,L5,A12,I10)')&
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

END SUBROUTINE FindInterfacesInRegion


SUBROUTINE ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides)
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
END SUBROUTINE ProlongToFace_ElementInfo

SUBROUTINE CountAndCreateMappings(TypeName,&
                                  isElem,isFace,isInterFace,&
                                  nElems,nFaces, nInterFaces,&
                                  ElemToX,XToElem,&
                                  FaceToX,XToFace,&
                                  FaceToXInter,XInterToFace,&
                                  DisplayInfoProcs)
!===================================================================================================================================
!> 1.) Count the number of Elements, Faces and Interfaces of the PML/Dielectric/... region
!> 2.) Create mappings from general element to PML/Dielectric/... element and vice versa
!>                                  face    to PML/Dielectric/... face or interface and vice vesa
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,     ONLY: nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                :: isElem(:),isFace(:),isInterFace(:)
CHARACTER(LEN=*),INTENT(IN)       :: TypeName
INTEGER,INTENT(IN),OPTIONAL       :: DisplayInfoProcs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)             :: nFaces,nInterFaces,nElems
INTEGER,ALLOCATABLE,INTENT(INOUT) :: ElemToX(:),XToElem(:),FaceToX(:),XToFace(:),FaceToXInter(:),XInterToFace(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iElem,iFace,nGlobalElems,nGlobalFaces,nGlobalInterFaces
INTEGER                           :: iXElem,iXFace,iXInterFace
INTEGER                           :: I,NProcs
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
IF(PRESENT(DisplayInfoProcs))THEN
  NProcs=DisplayInfoProcs
ELSE
  NProcs=0
END IF
#ifdef MPI
nGlobalElems=0
nGlobalFaces=0
nGlobalInterfaces=0
IF(0.EQ.myrank) WRITE(UNIT_stdOut,'(A)') "========================================================================================="
DO I=0,NProcs-1
  IF(I.EQ.myrank)THEN
    !write(*,'(A8,I10,A11,I10,A11,I10,A17,I10)')&
    !" myrank=",myrank," PP_nElems=",PP_nElems," nElems=",nElems," nGlobalElems=",nGlobalElems
    !write(*,'(A8,I10,A11,I10,A11,I10,A17,I10)')&
    !" myrank=",myrank," nSides=",nSides," nFaces=",nFaces," nGlobalFaces=",nGlobalFaces
    !write(*,'(A8,I10,A11,I10,A11,I10,A17,I10)')&
    !" myrank=",myrank," nSides=",nSides," nFaces=",nInterFaces," nGlobalFaces=",nGlobalInterFaces
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
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
DO I=0,NProcs-1
  IF(I.EQ.myrank)THEN
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
    WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
          ' myrank=',myrank,' Found ', nGlobalInterFaces,' nGlobal',TRIM(TypeName),"-InterFaces inside of ",TRIM(TypeName),'-region.'
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
END DO
IF(0.EQ.myrank) WRITE(UNIT_stdOut,'(A)') "========================================================================================="
#else
nGlobalElems=nElems
nGlobalFaces=nFaces
nGlobalInterFaces=nInterFaces
WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalElems     ,' nGlobal',TRIM(TypeName),"-Elems inside of      ",TRIM(TypeName),'-region.'
WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalFaces     ,' nGlobal',TRIM(TypeName),"-Faces inside of      ",TRIM(TypeName),'-region.'
WRITE(UNIT_stdOut,'(A8,I10,A,I10,A10,A10,A22,A10,A8)') &
      ' myrank=',myrank,' Found ', nGlobalInterFaces,' nGlobal',TRIM(TypeName),"-InterFaces inside of ",TRIM(TypeName),'-region.'
#endif /*MPI*/

!===================================================================================================================================
! create  mappings: element<->pml-element
!                      face<->pml-face
!                      face<->interface
!===================================================================================================================================
ALLOCATE(ElemToX(PP_nElems))
ALLOCATE(XToElem(nElems))
ALLOCATE(FaceToX(nSides))
ALLOCATE(XToFace(nFaces))
ALLOCATE(FaceToXInter(nSides))
ALLOCATE(XInterToFace(nInterFaces))
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


SUBROUTINE FinalizeInterfaces()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_Interfaces_Vars, ONLY:InterfaceRiemann,InterfacesInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
InterfacesInitIsDone = .FALSE.
SDEALLOCATE(InterfaceRiemann)
END SUBROUTINE FinalizeInterfaces


END MODULE MOD_Interfaces
