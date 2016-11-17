#include "boltzplatz.h"

MODULE MOD_Output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE Visualize
  MODULE PROCEDURE Visualize
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,Visualize,FinalizeOutput
!===================================================================================================================================

CONTAINS

SUBROUTINE InitOutput()
!===================================================================================================================================
! Initialize all output (and analyze) variables.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools,ONLY:GetStr,GetLogical,GETINT
USE MOD_Output_Vars,ONLY:NVisu,ProjectName,OutputInitIsDone,OutputFormat
USE MOD_Interpolation_Vars,ONLY:xGP,wBary,InterpolationInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
CHARACTER(LEN=255)             :: LogFile
LOGICAL                        :: WriteErrorFiles
CHARACTER(LEN=40)              :: DefStr
!===================================================================================================================================
IF ((.NOT.InterpolationInitIsDone).OR.OutputInitIsDone) THEN
  CALL abort(&
      __STAMP__&
      ,'InitOutput not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

WRITE(DefStr,'(i4)') PP_N
NVisu=GETINT('NVisu',DefStr)
CALL initOutputBasis(PP_N,NVisu,xGP,wBary)
! Name for all output files
ProjectName=GETSTR('ProjectName')
Logging    =GETLOGICAL('Logging','.FALSE.')

WriteErrorFiles=GETLOGICAL('WriteErrorFiles','.TRUE.')
IF(WriteErrorFiles)THEN
  ! Open file for error output
  WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'
  OPEN(UNIT=UNIT_errOut,  &
       FILE=ErrorFileName,&
       STATUS='REPLACE',  &
       ACTION='WRITE',    &
       IOSTAT=OpenStat)
ELSE
  UNIT_errOut=UNIT_stdOut
END IF

OutputFormat = GETINT('OutputFormat','1')
! Open file for logging
IF(Logging)THEN
  WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
  OPEN(UNIT=UNIT_logOut,  &
       FILE=LogFile,      &
       STATUS='UNKNOWN',  &
       ACTION='WRITE',    &
       POSITION='APPEND', &
       IOSTAT=OpenStat)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging

OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput



SUBROUTINE InitOutputBasis(N_in,NVisu_in,xGP,wBary)
!===================================================================================================================================
! Initialize all output variables.
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars, ONLY:Vdm_GaussN_NVisu
USE MOD_Basis,       ONLY:InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                  :: N_in,NVisu_in
REAL,INTENT(IN),DIMENSION(0:N_in)   :: xGP,wBary
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:NVisu_in)  :: XiVisu
INTEGER                     :: i
!===================================================================================================================================
!equidistant visu points
DO i=0,NVisu_in
  XiVisu(i) = 2./REAL(NVisu_in) * REAL(i) - 1. 
END DO
! Gauss/Gl -> Visu : computation -> visualization
ALLOCATE(Vdm_GaussN_NVisu(0:NVisu_in,0:N_in))
CALL InitializeVandermonde(N_in,NVisu_in,wBary,xGP,XiVisu,Vdm_GaussN_NVisu)
END SUBROUTINE InitOutputBasis


#ifdef PARTICLES
SUBROUTINE Visualize(OutputTime)
#else
SUBROUTINE Visualize()
#endif /*PARTICLES/*
!===================================================================================================================================
! Simple visualization of conservative variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Output_Vars,ONLY:ProjectName,OutputFormat
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars,ONLY:Elem_xGP
USE MOD_Output_Vars,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasis,ONLY:ChangeBasis3D
USE MOD_Tecplot,ONLY:WriteDataToTecplotBinary
#ifdef PARTICLES
USE MOD_OutPutVTK,ONLY:WriteDataToVTK
USE MOD_Particle_Output_Vars, ONLY: WriteFieldsToVTK
#endif /*PARTICLES*/
!USE MOD_Eval_xyz,ONLY:eval_xyz
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
#ifdef PARTICLES
REAL,INTENT(IN)               :: OutputTime
#endif /*PARTICLES/*
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef PARTICLES
INTEGER                       :: iElem
REAL                          :: Coords_NVisu(3,0:NVisu,0:NVisu,0:NVisu,1:PP_nElems)
REAL                          :: U_NVisu(PP_nVar,0:NVisu,0:NVisu,0:NVisu,1:PP_nElems)
CHARACTER(LEN=255)            :: FileString
#endif /*PARTICLES/*
CHARACTER(LEN=32),ALLOCATABLE :: VarNames(:) ! Output variable names
!===================================================================================================================================
IF(outputFormat.LE.0) RETURN
! Specify output names
#if PP_nVar == 5
ALLOCATE(VarNames(PP_nVar))
VarNames(1)       ='Density'
VarNames(2)       ='MomentumX'
VarNames(3)       ='MomentumY'
VarNames(4)       ='MomentumZ'
VarNames(PP_nVar) ='EnergyStagnationDensity'
#endif
#if PP_nVar == 1
ALLOCATE(VarNames(PP_nVar))
VarNames(PP_nVar) ='Solution'
#endif
#if PP_nVar == 8
ALLOCATE(VarNames(PP_nVar))
VarNames(1)='ElectricFieldX'
VarNames(2)='ElectricFieldY'
VarNames(3)='ElectricFieldZ'
VarNames(4)='MagneticFieldX'
VarNames(5)='MagneticFieldY'
VarNames(6)='MagneticFieldZ'
VarNames(7)='Phi'       
VarNames(8)='Psi'       
#endif
#if PP_nVar == 4
ALLOCATE(VarNames(PP_nVar))
VarNames(1)='ElectricFieldX'
VarNames(2)='ElectricFieldY'
VarNames(3)='ElectricFieldZ'    
VarNames(4)='Psi'       
#endif

#ifdef PARTICLES
FileString=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Solution',OutputTime))
IF(WriteFieldsToVTK) THEN
  FileString=TRIM(INTSTAMP(TRIM(FileString),myRank))//'.vtk'
ELSE
  FileString=TRIM(INTSTAMP(TRIM(FileString),myRank))//'.dat'
END IF
DO iElem=1,PP_nElems
  ! Create coordinates of visualization points
  CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  ! Interpolate solution onto visu grid
  CALL ChangeBasis3D(PP_nVar,PP_N,NVisu,Vdm_GaussN_NVisu,U(1:PP_nVar,:,:,:,iElem),U_NVisu(1:PP_nVar,:,:,:,iElem))
END DO !iElem
! Visualize data
IF(WriteFieldsToVTK) THEN
  CALL WriteDataToVTK(NVisu,PP_nElems,PP_nVar,VarNames,Coords_NVisu(1:3,:,:,:,:),U_NVisu,FileString)
ELSE
  CALL WriteDataToTecplotBinary(NVisu,PP_nElems,PP_nVar,0,VarNames,Coords_NVisu(1:3,:,:,:,:),U_NVisu,TRIM(FileString))
END IF

IF(1.EQ.2)THEN
  WRITE(*,*) (OutputTime)
END IF
#endif /*PARTICLES*/

! test eval_xyz 
!DO iElem=1,PP_nElems
!  eval_vec = (/0.3,0.4,0.5/)
!  !eval_vec = (/30,40,50/)
!  CALL eval_xyz(eval_vec,6,PP_N,U(1:6,:,:,:,iElem),U_eval,iElem)
!END DO !iElem

END SUBROUTINE Visualize



SUBROUTINE FinalizeOutput()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:Vdm_GaussN_NVisu,OutputInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(Vdm_GaussN_NVisu)
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
