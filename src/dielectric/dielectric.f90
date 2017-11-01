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
#ifndef PP_HDG
USE MOD_HDF5_output,     ONLY: WriteDielectricGlobalToHDF5
#endif /*PP_HDG*/
USE MOD_Equation_Vars,   ONLY: c_corr,c
USE MOD_Interfaces,      ONLY: FindInterfacesInRegion,FindElementInRegion,CountAndCreateMappings,DisplayRanges,SelectMinMaxRegion
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
IF(.NOT.DoDielectric) THEN
  SWRITE(UNIT_stdOut,'(A)') ' Dielectric region deactivated. '
  nDielectricElems=0
  RETURN
END IF
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
! determine Dielectric elements
xyzPhysicalMinMaxDielectric(1:6) = GETREALARRAY('xyzPhysicalMinMaxDielectric',6,'0.0,0.0,0.0,0.0,0.0,0.0')
xyzDielectricMinMax(1:6)         = GETREALARRAY('xyzDielectricMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
! use xyzPhysicalMinMaxDielectric before xyzDielectricMinMax: 
! 1.) check for xyzPhysicalMinMaxDielectric 
! 2.) check for xyzDielectricMinMax
CALL SelectMinMaxRegion('Dielectric',useDielectricMinMax,&
                        'xyzPhysicalMinMaxDielectric',xyzPhysicalMinMaxDielectric,&
                        'xyzDielectricMinMax',xyzDielectricMinMax)

! display ranges of Dielectric region depending on useDielectricMinMax
!CALL DisplayRanges('useDielectricMinMax',useDielectricMinMax,&
                   !'xyzDielectricMinMax',xyzDielectricMinMax(1:6),&
           !'xyzPhysicalMinMaxDielectric',xyzPhysicalMinMaxDielectric(1:6))

! find all elements in the Dielectric region
IF(useDielectricMinMax)THEN ! find all elements located inside of 'xyzMinMax'
  CALL FindElementInRegion(isDielectricElem,xyzDielectricMinMax,&
                           ElementIsInside=.TRUE.,DisplayInfo=.TRUE.) ! pure Maxwell simulations
ELSE ! find all elements located outside of 'xyzPhysicalMinMaxDielectric'
  CALL FindElementInRegion(isDielectricElem,xyzPhysicalMinMaxDielectric,&
                           ElementIsInside=.FALSE.,DisplayInfo=.TRUE.)
END IF

! find all faces in the Dielectric region
CALL FindInterfacesInRegion(isDielectricFace,isDielectricInterFace,isDielectricElem)

! Get number of Dielectric Elems, Faces and Interfaces. Create Mappngs Dielectric <-> physical region
CALL CountAndCreateMappings('Dielectric',&
                            isDielectricElem,isDielectricFace,isDielectricInterFace,&
                            nDielectricElems,nDielectricFaces, nDielectricInterFaces,&
                            ElemToDielectric,DielectricToElem,& ! these two are allocated
                            FaceToDielectric,DielectricToFace,& ! these two are allocated
                            FaceToDielectricInter,DielectricInterToFace,& ! these two are allocated
                            DisplayInfo=.TRUE.)
#ifndef PP_HDG
! Set the dielectric profile function EpsR,MuR=f(x,y,z) in the Dielectric region: only for Maxwell
CALL SetDielectricVolumeProfile()

! Determine dielectric Values on faces and communicate them: only for Maxwell
CALL SetDielectricFaceProfile()

! create a HDF5 file containing the DielectriczetaGlobal field: only for Maxwell
CALL WriteDielectricGlobalToHDF5()
#endif /*PP_HDG*/

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
! Fish eye lens: half sphere filled with gradually changing dielectric medium
IF(TRIM(DielectricTestCase).EQ.'FishEyeLens')THEN
  ! use function with radial dependence: EpsR=n0^2 / (1 + (r/r_max)^2)^2
  DO iDielectricElem=1,nDielectricElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r = SQRT(Elem_xGP(1,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(2,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(3,i,j,k,DielectricToElem(iDielectricElem))**2  )
    DielectricEps(i,j,k,iDielectricElem) = 4./((1+(r/DielectricRmax)**2)**2)
  END DO; END DO; END DO; END DO !iDielectricElem,k,j,i
  DielectricMu(0:PP_N,0:PP_N,0:PP_N,1:nDielectricElems) = DielectricMuR
! Sphere filled with constant dielectric medium
ELSEIF(TRIM(DielectricTestCase).EQ.'Sphere')THEN
  ! use function with radial dependence: EpsR=n0^2 / (1 + (r/r_max)^2)^2
  DO iDielectricElem=1,nDielectricElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r = SQRT(Elem_xGP(1,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(2,i,j,k,DielectricToElem(iDielectricElem))**2+&
             Elem_xGP(3,i,j,k,DielectricToElem(iDielectricElem))**2  )
    DielectricEps(i,j,k,iDielectricElem) = DielectricRmax
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


#ifndef PP_HDG
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
#endif /*PP_HDG*/


SUBROUTINE FinalizeDielectric()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
USE MOD_Dielectric_Vars,            ONLY: DoDielectric,DielectricEps,DielectricMu
USE MOD_Dielectric_Vars,            ONLY: ElemToDielectric,DielectricToElem,isDielectricElem
USE MOD_Dielectric_Vars,            ONLY: FaceToDielectric,DielectricToFace,isDielectricFace
USE MOD_Dielectric_Vars,            ONLY: FaceToDielectricInter,DielectricInterToFace,isDielectricInterFace
USE MOD_Dielectric_Vars,            ONLY: DielectricConstant_inv,Dielectric_Master,Dielectric_Slave
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
SDEALLOCATE(FaceToDielectricInter)
SDEALLOCATE(DielectricInterToFace)
SDEALLOCATE(FaceToDielectric)
SDEALLOCATE(isDielectricElem)
SDEALLOCATE(isDielectricFace)
SDEALLOCATE(isDielectricInterFace)
SDEALLOCATE(DielectricConstant_inv)
SDEALLOCATE(Dielectric_Master)
SDEALLOCATE(Dielectric_Slave)
END SUBROUTINE FinalizeDielectric





END MODULE MOD_Dielectric
