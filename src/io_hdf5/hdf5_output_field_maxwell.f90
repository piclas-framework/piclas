!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_HDF5_Output_Fields_Maxwell
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
USE MOD_io_HDF5
USE MOD_HDF5_output
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
#if !(USE_HDG) && !(USE_FV)
PUBLIC :: WritePMLzetaGlobalToHDF5
#endif /*!(USE_HDG) && !(USE_FV)*/
#if (PP_nVar==8)
PUBLIC :: WritePMLDataToHDF5
#endif /*(PP_nVar==8)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
!===================================================================================================================================

CONTAINS

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
#if !(USE_HDG) && !(USE_FV)
SUBROUTINE WritePMLzetaGlobalToHDF5()
!===================================================================================================================================
! write PMLzetaGlobal field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_PML_Vars             ,ONLY: PMLzeta0,PML,ElemToPML,isPMLElem
USE MOD_Mesh_Vars            ,ONLY: MeshFile,offsetElem,nElems
USE MOD_io_HDF5
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_ChangeBasis          ,ONLY: ChangeBasis3D
USE MOD_DG_vars              ,ONLY: N_DG_Mapping,nDofsMapping
USE MOD_HDF5_Output_ElemData ,ONLY: WriteAdditionalElemData
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: nVarOut=3
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName
REAL                :: StartT,EndT
REAL                :: OutputTime
INTEGER             :: iElem,Nloc
! p-adaption output
REAL,ALLOCATABLE    :: U_N_2D_local(:,:)
INTEGER             :: i,j,k,iDOF,nDOFOutput,offsetDOF
!===================================================================================================================================
#if USE_LOADBALANCE
IF(PerformLoadBalance) RETURN
#endif /*USE_LOADBALANCE*/
! ---------------------------------------------------------
! Prepare U_N_2D_local array for output as DG_Solution
! ---------------------------------------------------------
! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF
! Allocate local 2D array: create global Eps field for parallel output of Eps distribution
ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))
ALLOCATE(StrVarNames(1:nVarOut))
StrVarNames(1)='PMLzetaGlobalX'
StrVarNames(2)='PMLzetaGlobalY'
StrVarNames(3)='PMLzetaGlobalZ'

! Write into 2D array
iDOF = 0
IF(.NOT.ALMOSTZERO(PMLzeta0))THEN
  DO iElem=1,nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    IF(isPMLElem(iElem))THEN
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        iDOF = iDOF + 1
        U_N_2D_local(1:nVarOut,iDOF) = PML(ElemToPML(iElem))%zeta(1:nVarOut,i,j,k)/PMLzeta0
      END DO; END DO; END DO
    ELSE
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        iDOF = iDOF + 1
        U_N_2D_local(1:nVarOut,iDOF) = 0.0
      END DO; END DO; END DO
    END IF ! isPMLElem(iElem)
  END DO ! iElem=1,nElems
ELSE
  U_N_2D_local = 0.0
END IF

SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE PMLZetaGlobal TO HDF5 FILE...'
GETTIME(StartT)
OutputTime=0.0
! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
CALL GenerateFileSkeleton('PMLZetaGlobal',nVarOut,StrVarNames,TRIM(MeshFile),OutputTime,FileNameOut=FileName)
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesPMLzetaGlobal',nVarOut,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF ! MPIRoot

! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
CALL WriteAdditionalElemData(FileName,ElementOutNloc)

! ---------------------------------------------------------
! Output of DG_Solution
! ---------------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE(nVarOut         => INT(nVarOut,IK)           ,&
          nDofsMapping    => INT(nDofsMapping,IK)      ,&
          nDOFOutput      => INT(nDOFOutput,IK)        ,&
          offsetDOF       => INT(offsetDOF,IK)         )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName = 'DG_Solution' , rank = 2                , &
                        nValGlobal  = (/nVarOut     , nDofsMapping/)          , &
                        nVal        = (/nVarOut     , nDOFOutput/)            , &
                        offset      = (/0_IK        , offsetDOF/)             , &
                        collective  = .TRUE.        , RealArray = U_N_2D_local)
END ASSOCIATE
SDEALLOCATE(U_N_2D_local)
GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WritePMLzetaGlobalToHDF5
#endif /*!(USE_HDG) && !(USE_FV)*/


#if (PP_nVar==8)
SUBROUTINE WritePMLDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars          ,ONLY: offsetElem,nElems
USE MOD_PML_Vars           ,ONLY: DoPML,PMLnVar,isPMLElem
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping,nDofsMapping
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: iElem,Nloc
! p-adaption output
REAL,ALLOCATABLE               :: U_N_2D_local(:,:)
INTEGER                        :: i,j,k,iDOF,nDOFOutput,offsetDOF
!===================================================================================================================================

IF(.NOT.DoPML) RETURN

! ---------------------------------------------------------
! Prepare U_N_2D_local array for output as DG_Solution
! ---------------------------------------------------------
! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF
! Allocate local 2D array: create global Eps field for parallel output of Eps distribution
ALLOCATE(U_N_2D_local(1:PMLnVar,1:nDOFOutput))
ALLOCATE(StrVarNames(PMLnVar))
StrVarNames( 1)='PML-ElectricFieldX-P1'
StrVarNames( 2)='PML-ElectricFieldX-P2'
StrVarNames( 3)='PML-ElectricFieldX-P3'
StrVarNames( 4)='PML-ElectricFieldY-P4'
StrVarNames( 5)='PML-ElectricFieldY-P5'
StrVarNames( 6)='PML-ElectricFieldY-P6'
StrVarNames( 7)='PML-ElectricFieldZ-P7'
StrVarNames( 8)='PML-ElectricFieldZ-P8'
StrVarNames( 9)='PML-ElectricFieldZ-P9'
StrVarNames(10)='PML-MagneticFieldX-P10'
StrVarNames(11)='PML-MagneticFieldX-P11'
StrVarNames(12)='PML-MagneticFieldX-P12'
StrVarNames(13)='PML-MagneticFieldY-P13'
StrVarNames(14)='PML-MagneticFieldY-P14'
StrVarNames(15)='PML-MagneticFieldY-P15'
StrVarNames(16)='PML-MagneticFieldZ-P16'
StrVarNames(17)='PML-MagneticFieldZ-P17'
StrVarNames(18)='PML-MagneticFieldZ-P18'
StrVarNames(19)='PML-PhiB-P19'
StrVarNames(20)='PML-PhiB-P20'
StrVarNames(21)='PML-PhiB-P21'
StrVarNames(22)='PML-PsiE-P22'
StrVarNames(23)='PML-PsiE-P23'
StrVarNames(24)='PML-PsiE-P24'

DO iElem=1,nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  IF(isPMLElem(iElem))THEN
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:PMLnVar,iDOF) = U_N(iElem)%U2(1:PMLnVar,i,j,k)
    END DO; END DO; END DO
  ELSE
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:PMLnVar,iDOF) = 0.0
    END DO; END DO; END DO
  END IF ! isPMLElem(iElem)
END DO ! iElem=1,nElems

IF(MPIRoot) THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesPML',PMLnVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

! ---------------------------------------------------------
! Output of PML_Solution
! ---------------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE(PMLnVar         => INT(PMLnVar,IK)           ,&
          nDofsMapping    => INT(nDofsMapping,IK)      ,&
          nDOFOutput      => INT(nDOFOutput,IK)        ,&
          offsetDOF       => INT(offsetDOF,IK)         )
CALL GatheredWriteArray(FileName, create = .FALSE.                            , &
                        DataSetName = 'PML_Solution' , rank = 2                , &
                        nValGlobal  = (/PMLnVar      , nDofsMapping/)          , &
                        nVal        = (/PMLnVar      , nDOFOutput/)            , &
                        offset      = (/0_IK         , offsetDOF/)             , &
                        collective  = .TRUE.         , RealArray = U_N_2D_local)
END ASSOCIATE
SDEALLOCATE(U_N_2D_local)
DEALLOCATE(StrVarNames)


END SUBROUTINE WritePMLDataToHDF5
#endif /*(PP_nVar==8)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/

END MODULE MOD_HDF5_Output_Fields_Maxwell
