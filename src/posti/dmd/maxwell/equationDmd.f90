#include "piclas.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_EquationDMD
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitEquationDMD
  MODULE PROCEDURE InitEquationDMD
END INTERFACE

INTERFACE CalcEquationDMD
  MODULE PROCEDURE CalcEquationDMD
END INTERFACE

INTERFACE FinalizeEquationDMD
  MODULE PROCEDURE FinalizeEquationDMD
END INTERFACE

PUBLIC::InitEquationDMD,CalcEquationDMD,FinalizeEquationDMD
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize the visualization and map the variable names to classify these in conservative and derived quantities.
!===================================================================================================================================
SUBROUTINE InitEquationDMD()
! MODULES
USE MOD_Globals
USE MOD_DMD_Vars          ,ONLY:nVar_State,VarNames_State
USE MOD_EquationDMD_Vars
!USE MOD_EOS               ,ONLY: InitEOS
!USE MOD_EOS_Posti_Vars    ,ONLY: nVarDepEOS,DepTableEOS,DepNames
USE MOD_Readintools       ,ONLY: CountOption,GETSTR
USE MOD_StringTools       ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: iVar,countCons
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT EquationDMD ...'

!check if Varnames in HDF5 file are conservatives
! TODO: also generate mapping in case conservatives are there but not in the correct order
countCons=0
!WRITE (*,*) "nVar_State =", nVar_State
IF(nVar_State.GE.PP_nVar)THEN
  DO iVar=1,PP_nVar
    !IF (STRICMP(VarNames_State(iVar), DepNames(iVar))) THEN
      countCons=countCons+1
    !END IF
  END DO
END IF
IF(countCons.NE.PP_nVar) THEN
  CALL CollectiveStop(__STAMP__,'Not all necessary variables are present in HDF5 files')
END IF
!nVarDep=nVarDepEOS
nVarDep=nVar_State
ALLOCATE(VarNamesAll(nVarDep))
!VarNamesAll=DepNames
ALLOCATE(DepTable(nVarDep,0:nVarDep))
!DepTable=DepTableEOS

! generate mappings
!CALL Build_mapCalc_mapVisu()

! initialize EOS
!CALL InitEOS()


WRITE(UNIT_stdOut,'(A)')' INIT EquationDMD DONE!'
WRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitEquationDMD



!===================================================================================================================================
!> This routine computes the state on the visualization grid 
!===================================================================================================================================
SUBROUTINE CalcEquationDMD(DMDData,DMDData_out)
! MODULES
USE MOD_Globals
USE MOD_DMD_Vars            ,ONLY: nVar_State,VarNames_State,N_State,N_StateZ,nDoFs,nElems_State,nVarDMD,use2D
USE MOD_EquationDMD_Vars
!USE MOD_EOS_Posti           ,ONLY: CalcQuantities
USE MOD_StringTools         ,ONLY: STRICMP
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)    :: DMDData(nVar_State,N_State+1,N_State+1,N_StateZ+1,nElems_State)
REAL,INTENT(OUT)   :: DMDData_out(nDoFs*nVarDMD)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: maskCalc(nVarDep),nVal(4)
INTEGER            :: iVarOut,iVarIn,iVar
REAL,ALLOCATABLE   :: UCalc(:,:,:,:,:)
REAL,ALLOCATABLE   :: DMDData_tmp(:,:,:,:,:)
REAL,ALLOCATABLE   :: DMDData_tmp2(:,:)
REAL,ALLOCATABLE   :: DMDData_tmp3(:,:)
!INTEGER            :: OutputVarsIndex(MAXVAL(mapVisu))
!===================================================================================================================================
!WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)', ADVANCE='NO')" CONVERT DERIVED QUANTITIES..."
! CALCULATE DERIVED QUATITIES -----------------------------------------------------------------------------------------------------!
maskCalc=1
nVal=(/N_State+1,N_State+1,N_StateZ+1,nElems_State/)
nVarCalc=nVarDep
nVarVisuTotal=nVarDep
! Copy existing variables from solution array
! Attention: nVarCalc must be last dimension (needed for CalcQuantities from flexilib!)
IF(.NOT. use2D) THEN
  ALLOCATE(UCalc(N_State+1,N_State+1,N_State+1,nElems_State,nVarCalc))
  ALLOCATE(DMDData_tmp(nVarCalc,N_State+1,N_State+1,N_State+1,nElems_State))
  ALLOCATE(DMDData_tmp2(nElems_State*(N_State+1)**3,nVarCalc))
  ALLOCATE(DMDData_tmp3(nElems_State*(N_State+1)**3,nVarVisuTotal))
ELSE
  ALLOCATE(UCalc(N_State+1,N_State+1,1        ,nElems_State,nVarCalc))
  ALLOCATE(DMDData_tmp(nVarCalc,N_State+1,N_State+1,1,nElems_State))
  ALLOCATE(DMDData_tmp2(nElems_State*(N_State+1)**2,nVarCalc))
  ALLOCATE(DMDData_tmp3(nElems_State*(N_State+1)**2,nVarVisuTotal))
END IF

DO iVarOut=1,nVarDep ! iterate over all out variables
  !WRITE (*,*) "iVarOut,mapCalc(iVarOut) =", iVarOut,mapCalc(iVarOut)
  !IF (mapCalc(iVarOut).LT.1) CYCLE ! check if variable must be calculated
  DO iVarIn=1,nVar_State ! iterate over all in variables
    !IF( STRICMP(VarNamesAll(iVarOut),VarNames_State(iVarIn))) THEN
      !UCalc(:,:,:,:,mapCalc(iVarOut))=DMDData(iVarIn,:,:,:,:)
      UCalc(:,:,:,:,iVarOut)=DMDData(iVarIn,:,:,:,:)
      !maskCalc(iVarOut)=0 ! remove variable from maskCalc, since they now got copied and must not be calculated.
    !END IF
  END DO
eND DO

! calculate all quantities
!CALL CalcQuantities(nVarCalc,nVal,(/1/),mapCalc,UCalc,maskCalc)

! fill output array
DMDData_tmp2=RESHAPE(UCalc,(/nDofs,nVarCalc/))

DO iVar=1,nVarDep
  !IF (mapVisu(iVar) .GT. 0 ) THEN
    !DMDData_tmp3(:,mapVisu(ivar))=DMDData_tmp2(:,mapCalc(iVar))
    DMDData_tmp3(:,ivar)=DMDData_tmp2(:,iVar)
  !END IF
END Do

DMDData_out(:)=RESHAPE(TRANSPOSE(DMDData_tmp3(:,:)), (/nDoFs*nVarDMD/))

DEALLOCATE(UCalc)
DEALLOCATE(DMDData_tmp)

WRITE(UNIT_stdOut,'(A)')" DONE!"
END SUBROUTINE CalcEquationDMD


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE FinalizeEquationDMD()
! MODULES
USE MOD_Globals
USE MOD_EquationDMD_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
DEALLOCATE(TransMap,is2D) 
WRITE(UNIT_stdOut,'(A)') '  EquationDMD FINALIZED'
END SUBROUTINE FinalizeEquationDMD

END MODULE MOD_EquationDMD
