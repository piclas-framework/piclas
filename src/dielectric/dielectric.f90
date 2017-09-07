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
USE MOD_Interfaces,      ONLY: FindInterfacesInRegion,FindElementInRegion,CountAndCreateMappings
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
  CALL FindElementInRegion(isDielectricElem,xyzDielectricMinMax,&
                           ElementIsInside=.TRUE.,DisplayInfoProcs=DielectricprintInfoProcs)
ELSE
  CALL FindElementInRegion(isDielectricElem,xyzPhysicalMinMaxDielectric,&
                           ElementIsInside=.FALSE.,DisplayInfoProcs=DielectricprintInfoProcs)
END IF

! find all faces in the Dielectric region
CALL FindInterfacesInRegion(isDielectricFace,isDielectricInterFace,isDielectricElem,DielectricprintInfoProcs)

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
! LOCAL VARIABLE,Dielectric_dummy_Master2S
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,1:nSides)           :: Dielectric_dummy_Master 
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,1:nSides)           :: Dielectric_dummy_Slave  
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) :: Dielectric_dummy_elem   
#ifdef MPI
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Master2
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Slave2 
#endif /*MPI*/
INTEGER                                                  :: iElem,I,J,iSide
!===================================================================================================================================
Dielectric_dummy_elem    = 0. ! default is an invalid number
Dielectric_dummy_Master  = 0.
Dielectric_dummy_Slave   = 0.

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

! re-map data (for MPI communication)
Dielectric_dummy_Master2 = 0.
Dielectric_dummy_Slave2  = 0.
DO I=0,PP_N
  DO J=0,PP_N
    DO iSide=1,nSides
      Dielectric_dummy_Master2(1,I,J,iSide)=Dielectric_dummy_Master(1,I,J,iSide)
      Dielectric_dummy_Slave2 (1,I,J,iSide)=Dielectric_dummy_Slave( 1,I,J,iSide)
    END DO
  END DO
END DO

! send Slave Dielectric info (real array with dimension (N+1)*(N+1)) to Master procs
CALL StartReceiveMPIData(1,Dielectric_dummy_Slave2 ,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,Dielectric_dummy_Slave2 ,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

! send Master Dielectric info (real array with dimension (N+1)*(N+1)) to Slave procs
CALL StartReceiveMPIData(1,Dielectric_dummy_Master2,1,nSides ,RecRequest_U ,SendID=1) ! Receive YOUR
CALL StartSendMPIData(   1,Dielectric_dummy_Master2,1,nSides,SendRequest_U ,SendID=1) ! Send MINE

CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE - receive YOUR
CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U ,SendID=1) !Send YOUR - receive MINE 
#endif /*MPI*/

ALLOCATE(Dielectric_Master(0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Dielectric_Slave( 0:PP_N,0:PP_N,1:nSides))


#ifdef MPI
Dielectric_Master=Dielectric_dummy_Master2(1,0:PP_N,0:PP_N,1:nSides)
Dielectric_Slave =Dielectric_dummy_Slave2( 1,0:PP_N,0:PP_N,1:nSides)
#else
Dielectric_Master=Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides)
Dielectric_Slave =Dielectric_dummy_Slave( 1,0:PP_N,0:PP_N,1:nSides)
#endif /*MPI*/

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





END MODULE MOD_Dielectric
