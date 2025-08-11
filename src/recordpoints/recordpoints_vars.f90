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

!==================================================================================================================================
!> Variables needed for the evaluation of the record points
!==================================================================================================================================
MODULE MOD_RecordPoints_Vars
! MODULES
#if USE_MPI
USE mpi_f08
#endif /*USE_MPI*/

IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255) :: RPDefFile                        !< File containing element-local parametric recordpoint coordinates and structure
LOGICAL            :: RecordPointsInitIsDone = .FALSE. !< mark wheter recordpoints init routine is finished
LOGICAL            :: RP_inUse  = .FALSE.              !< mark whether recordpoints should be evaluated during computation
LOGICAL            :: RP_onProc = .FALSE.              !< marks wheter current proc has RPs
LOGICAL            :: RP_fileExists = .FALSE.          !< flag if RP file for analyze level has been created
INTEGER            :: RP_Buffersize                    !< no. of time samples (size of RP_Data)
INTEGER            :: RP_MaxBuffersize                 !< max. allowed no. of time samples
INTEGER            :: RP_SamplingOffset                !< sampling rate (each .. iterations)
INTEGER            :: nRP                              !< no. of RP on proc
INTEGER            :: nGlobalRP                        !< total no. of RP
INTEGER            :: offsetRP                         !< offset for each proc in global RP list
INTEGER            :: iSample=0                        !< no of samples in array
INTEGER            :: nSamples=0                       !< total no. samples in case of multiple io steps
INTEGER            :: chunkSamples=0                   !< time samples per chunk for IO (first iSample in file)
INTEGER,ALLOCATABLE:: RP_ElemID(:)                     !< mapping from RP->Elem (nRP)
REAL,ALLOCATABLE   :: L_xi_RP(:,:)                     !< Lagrange basis evaluated at RP coords (xi-dir)
REAL,ALLOCATABLE   :: L_eta_RP(:,:)                    !< Lagrange basis evaluated at RP coords (eta-dir)
REAL,ALLOCATABLE   :: L_zeta_RP(:,:)                   !< Lagrange basis evaluated at RP coords (zeta-dir)
REAL,ALLOCATABLE   :: RP_Data(:,:,:)                   !< solution evaluated at RPs (nvar,nRP,nSamples)
REAL,ALLOCATABLE   :: lastSample(:,:)                  !< solution evaluated at RPs (nvar,nRP,nSamples)
CHARACTER(LEN=255) :: StrVarNames(PP_nVar)             !< RP variables names for output

!----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for RPs
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER            :: myRPrank
#if USE_MPI
TYPE(MPI_Comm)     :: RP_COMM=MPI_COMM_NULL   !< MPI RP communicator
#endif /*USE_MPI*/
INTEGER            :: nRP_Procs

END MODULE MOD_recordPoints_Vars