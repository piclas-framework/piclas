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

MODULE MOD_DVM_Boundary_Analyze
!===================================================================================================================================
!> Particle boundary sampling: calculation and output of heat flux, forces, and impact properties at boundaries
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::WriteDVMSurfToHDF5, InitDVMBoundaryAnalyze, FinalizeDVMBoundaryAnalyze
!===================================================================================================================================

CONTAINS

SUBROUTINE InitDVMBoundaryAnalyze()
!===================================================================================================================================
! Initialize DVM surface output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Equation_Vars_FV         ,ONLY: DVMSurfaceValues, nVarDVMSurf
USE MOD_Mesh_Vars                ,ONLY: nBCSides,nBCs,BoundaryName
USE MOD_Mesh_Vars_FV             ,ONLY: nGlobalBCSides,offsetBCSide,FVSurfBCName,nFVSurfBC,BoundaryType_FV
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                :: iBC
CHARACTER(LEN=255),ALLOCATABLE         :: BCName(:)
#if USE_MPI
INTEGER                                :: sendbuf,recvbuf
#endif /*USE_MPI*/
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') ' INIT FV SURFACE SIDES ...'

ALLOCATE(DVMSurfaceValues(1:nVarDVMSurf,1:nBCSides))

! Find global number of surf sides and offset for each proc
#if USE_MPI
sendbuf = nBCSides
recvbuf = 0
CALL MPI_EXSCAN(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_PICLAS,iError)
offsetBCSide = recvbuf
! last proc knows total number of BC sides
sendbuf = offsetBCSide + nBCSides
CALL MPI_BCAST(sendbuf,1,MPI_INTEGER,nProcessors-1,MPI_COMM_PICLAS,iError)
nGlobalBCSides = sendbuf
#else
offsetBCSide  = 0
nGlobalBCSides = nBCSides
#endif /*USE_MPI*/

! create boundary name mapping for surfaces SurfaceBC number mapping
nFVSurfBC = 0
ALLOCATE(BCName(1:nBCs))
DO iBC=1,nBCs
  BCName=''
END DO
DO iBC=1,nBCs
  IF (BoundaryType_FV(iBC,BC_TYPE).EQ.1) CYCLE !periodic BC
  nFVSurfBC         = nFVSurfBC + 1
  BCName(nFVSurfBC) = BoundaryName(iBC)
END DO

IF (nFVSurfBC.GE.1) THEN
ALLOCATE(FVSurfBCName(1:nFVSurfBC))
  DO iBC=1,nFVSurfBC
    FVSurfBCName(iBC) = BCName(iBC)
  END DO
END IF
DEALLOCATE(BCName)

LBWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')'DONE'

END SUBROUTINE InitDVMBoundaryAnalyze

SUBROUTINE WriteDVMSurfToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
!> write the suerface values to a HDF5 state file
!> now with general communicator as surface communicator only built for particles
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: ProjectName
USE MOD_HDF5_Output             ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_IO_HDF5
USE MOD_Mesh_Vars               ,ONLY: nBCSides
USE MOD_Mesh_Vars_FV            ,ONLY: nGlobalBCSides,offsetBCSide,FVSurfBCName,nFVSurfBC
USE MOD_Equation_Vars_FV        ,ONLY: DVMSurfaceValues,nVarDVMSurf
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
REAL,INTENT(IN)                      :: OutputTime
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: FileName,FileString,Statedummy
CHARACTER(LEN=255)                  :: H5_Name
CHARACTER(LEN=255)                  :: NodeTypeTemp
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: nVarCount
REAL                                :: StartT,EndT
!===================================================================================================================================

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iERROR)
#endif /*USE_MPI*/

! Return if no sampling sides
IF (nGlobalBCSides      .EQ.0) RETURN

IF (myRank.EQ.0) THEN
    WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE DVMSurfSTATE TO HDF5 FILE...'
    StartT=LOCALTIME()
END IF
FileName   = TIMESTAMP(TRIM(ProjectName)//'_DVMSurfState',OutputTime)
FileString = TRIM(FileName)//'.h5'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (myRank.EQ.0) THEN
#endif
    CALL OpenDataFile(FileString,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    Statedummy = 'DVMSurfState'

    ! Write file header
    CALL WriteHDF5Header(Statedummy,File_ID)
    CALL WriteAttributeToHDF5(File_ID,'DSMC_nSurfSample',1,IntegerScalar=1) ! dummy for piclas2vtk
    CALL WriteAttributeToHDF5(File_ID,'MeshFile'        ,1,StrScalar=(/TRIM(MeshFileName)/))
    CALL WriteAttributeToHDF5(File_ID,'Time'            ,1,RealScalar=OutputTime)
    CALL WriteAttributeToHDF5(File_ID,'BC_Surf'         ,nFVSurfBC,StrArray=FVSurfBCName)
    NodeTypeTemp='VISU'
    CALL WriteAttributeToHDF5(File_ID,'NodeType'        ,1,StrScalar=(/NodeTypeTemp/))

    ALLOCATE(Str2DVarNames(1:nVarDVMSurf))
    Str2DVarNames(:) = ''
    nVarCount        = 1

    ! fill varnames for total values (add new variables to MACROSURF_NVARS)
    CALL AddVarName(Str2DVarNames,nVarDVMSurf,nVarCount,'ForcePerAreaX')
    CALL AddVarName(Str2DVarNames,nVarDVMSurf,nVarCount,'ForcePerAreaY')
    CALL AddVarName(Str2DVarNames,nVarDVMSurf,nVarCount,'ForcePerAreaZ')
    CALL AddVarName(Str2DVarNames,nVarDVMSurf,nVarCount,'HeatFlux')
    ! CALL AddVarName(Str2DVarNames,nVarDVMSurf,nVarCount,'iBC')

    CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVarDVMSurf,StrArray=Str2DVarNames)

    CALL CloseDataFile()
    DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF

CALL MPI_BARRIER(MPI_COMM_PICLAS,iERROR)

CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
#else
CALL OpenDataFile(FileString,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

WRITE(H5_Name,'(A)') 'SurfaceData'
! WARNING: Only the sampling leaders write the data to .h5
ASSOCIATE (&
        nGlobalSides      => INT(nGlobalBCSides,IK), &
        nLocalSides       => INT(nBCSides,IK)      , &
        offsetSurfSide    => INT(offsetBCSide,IK)  , &
        SurfOutputSize    => INT(nVarDVMSurf,IK))
    CALL WriteArrayToHDF5(DataSetName=H5_Name            , rank=2                   , &
                        nValGlobal =(/SurfOutputSize, nGlobalSides/)    , &
                        nVal       =(/SurfOutputSize, nLocalSides/)     , &
                        offset     =(/0_IK,  offsetSurfSide/)           , &
                        collective =.FALSE.                                         , &
                        RealArray  = DVMSurfaceValues(1:SurfOutputSize,1:nLocalSides))
END ASSOCIATE

CALL CloseDataFile()

IF (myRank.EQ.0) THEN
    EndT=LOCALTIME()
    CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE., rank=myRank)
END IF

END SUBROUTINE WriteDVMSurfToHDF5

SUBROUTINE AddVarName(StrArray,ArrayDim,idx,VarName)
!----------------------------------------------------------------------------------------------------------------------------------!
! description
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: ArrayDim
CHARACTER(LEN=*),INTENT(INOUT) :: StrArray(ArrayDim)
INTEGER,INTENT(INOUT)          :: idx
CHARACTER(LEN=*),INTENT(IN)    :: VarName
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
StrArray(idx)=TRIM(VarName)
idx=idx+1
END SUBROUTINE AddVarName

SUBROUTINE FinalizeDVMBoundaryAnalyze()
!===================================================================================================================================
! Finalize DVM surface output
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars_FV         ,ONLY: DVMSurfaceValues
USE MOD_Mesh_Vars_FV             ,ONLY: FVSurfBCName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===================================================================================================================================

SDEALLOCATE(DVMSurfaceValues)
SDEALLOCATE(FVSurfBCName)

END SUBROUTINE FinalizeDVMBoundaryAnalyze

END MODULE MOD_DVM_Boundary_Analyze