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

MODULE MOD_MPI_Output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
#if USE_MPI
USE mpi_f08
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if USE_MPI && defined(MEASURE_MPI_WAIT)
PUBLIC::OutputMPIW8Time
#endif /*USE_MPI && defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

CONTAINS

#if USE_MPI && defined(MEASURE_MPI_WAIT)
!===================================================================================================================================
!> Root writes measured MPI_WAIT() times to disk
!===================================================================================================================================
SUBROUTINE OutputMPIW8Time()
! MODULES
USE MOD_Globals
USE MOD_MPI_Vars          ,ONLY: MPIW8TimeGlobal      , MPIW8TimeProc      , MPIW8TimeField      , MPIW8Time      , MPIW8TimeBaS
USE MOD_MPI_Vars          ,ONLY: MPIW8CountGlobal , MPIW8CountProc , MPIW8CountField , MPIW8Count , MPIW8CountBaS
USE MOD_MPI_Vars          ,ONLY: MPIW8TimeSim,MPIW8TimeMM, MPIW8CountMM
USE MOD_StringTools       ,ONLY: INTTOSTR
#if defined(PARTICLES)
USE MOD_Particle_MPI_Vars ,ONLY: MPIW8TimePart,MPIW8CountPart
#endif /*defined(PARTICLES)*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                                :: WriteHeader
INTEGER                                :: ioUnit,i
CHARACTER(LEN=150)                     :: formatStr
CHARACTER(LEN=22),PARAMETER            :: outfile    ='MPIW8Time.csv'
CHARACTER(LEN=22),PARAMETER            :: outfilePerc='MPIW8TimePercent.csv'
CHARACTER(LEN=22),PARAMETER            :: outfileProc='MPIW8TimeProc'
CHARACTER(LEN=30)                      :: outfileProc_loc
CHARACTER(LEN=10)                      :: hilf
INTEGER,PARAMETER                      :: nTotalVars =2*MPIW8SIZE+2
CHARACTER(LEN=255),DIMENSION(nTotalVars) :: StrVarNames(nTotalVars)=(/ CHARACTER(LEN=255) :: &
    'nProcessors'                 , &
    'WallTimeSim'                 , &
    'Barrier-and-Sync'            , &
    'Barrier-and-Sync-Counter'    , &
    'RAM-Measure-Reduce'          , &
    'RAM-Measure-Reduce-Counter'    &
#if USE_HDG
   ,'HDG-SendLambda'              , & ! (1)
    'HDG-SendLambda-Counter'      , & ! (1)
    'HDG-ReceiveLambda'           , & ! (2)
    'HDG-ReceiveLambda-Counter'   , & ! (2)
    'HDG-Broadcast'               , & ! (3)
    'HDG-Broadcast-Counter'       , & ! (3)
    'HDG-Allreduce'               , & ! (4)
    'HDG-Allreduce-Counter'         & ! (4)
#else
   ,'DGSEM-Send'                  , &     ! (1)
    'DGSEM-Send-Counter'          , &     ! (1)
    'DGSEM-Receive'               , &     ! (2)
    'DGSEM-Receive-Counter'         &     ! (2)
#endif /*USE_HDG*/
#if defined(PARTICLES)
   ,'SendNbrOfParticles'          , & ! (1)
    'SendNbrOfParticles-Counter'  , & ! (1)
    'RecvNbrOfParticles'          , & ! (2)
    'RecvNbrOfParticles-Counter'  , & ! (2)
    'SendParticles'               , & ! (3)
    'SendParticles-Counter'       , & ! (3)
    'RecvParticles'               , & ! (4)
    'RecvParticles-Counter'       , & ! (4)
    'EmissionParticles'           , & ! (5)
    'EmissionParticles-Counter'   , & ! (5)
    'PIC-depo-Wait'               , & ! (6)
    'PIC-depo-Wait-Counter'         & ! (6)
#endif /*defined(PARTICLES)*/
    /)
! CHARACTER(LEN=255),DIMENSION(nTotalVars) :: StrVarNamesProc(nTotalVars)=(/ CHARACTER(LEN=255) :: &
!     'Rank'                  &
! #if defined(PARTICLES)
!    ,'SendNbrOfParticles'  , &
!     'RecvNbrOfParticles'  , &
!     'SendParticles'       , &
!     'RecvParticles'         &
! #endif /*defined(PARTICLES)*/
!     /)
CHARACTER(LEN=255)         :: tmpStr(nTotalVars)
CHARACTER(LEN=1000)        :: tmpStr2
CHARACTER(LEN=1),PARAMETER :: delimiter=","
REAL                       :: TotalSimTime,MPIW8TimeSimeGlobal,TotalCounter
!===================================================================================================================================
MPIW8Time(                1:1)                              = MPIW8TimeBaS
MPIW8Count(               1:1)                              = MPIW8CountBaS
MPIW8Time(                2:2)                              = MPIW8TimeMM
MPIW8Count(               2:2)                              = MPIW8CountMM
MPIW8Time(                3:MPIW8SIZEFIELD+2)               = MPIW8TimeField
MPIW8Count(               3:MPIW8SIZEFIELD+2)               = MPIW8CountField
#if defined(PARTICLES)
MPIW8Time( MPIW8SIZEFIELD+3:MPIW8SIZEFIELD+MPIW8SIZEPART+2) = MPIW8TimePart
MPIW8Count(MPIW8SIZEFIELD+3:MPIW8SIZEFIELD+MPIW8SIZEPART+2) = MPIW8CountPart
#endif /*defined(PARTICLES)*/

! Collect and output measured MPI_WAIT() times
IF(MPIroot)THEN
  ALLOCATE(MPIW8TimeProc(MPIW8SIZE*nProcessors))
  ALLOCATE(MPIW8CountProc(MPIW8SIZE*nProcessors))
  CALL MPI_REDUCE(MPIW8TimeSim , MPIW8TimeSimeGlobal , 1         , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(MPIW8Time    , MPIW8TimeGlobal     , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_REDUCE(MPIW8Count   , MPIW8CountGlobal    , MPIW8SIZE , MPI_INTEGER8         , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)

  CALL MPI_GATHER(MPIW8Time  , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPIW8TimeProc  , MPIW8SIZE , MPI_DOUBLE_PRECISION , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_GATHER(MPIW8Count , MPIW8SIZE , MPI_INTEGER8         , MPIW8CountProc , MPIW8SIZE , MPI_INTEGER8         , 0 , MPI_COMM_PICLAS , iError)
ELSE
  CALL MPI_REDUCE(MPIW8TimeSim , 0 , 1         , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IError)
  CALL MPI_REDUCE(MPIW8Time    , 0 , MPIW8SIZE , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IError)
  CALL MPI_REDUCE(MPIW8Count   , 0 , MPIW8SIZE , MPI_INTEGER8         , MPI_SUM , 0 , MPI_COMM_PICLAS , IError)

  CALL MPI_GATHER(MPIW8Time  , MPIW8SIZE , MPI_DOUBLE_PRECISION , 0              , 0         , MPI_DOUBLE_PRECISION , 0 , MPI_COMM_PICLAS , iError)
  CALL MPI_GATHER(MPIW8Count , MPIW8SIZE , MPI_INTEGER8         , 0              , 0         , MPI_INTEGER8         , 0 , MPI_COMM_PICLAS , iError)
END IF

! --------------------------------------------------
! Only MPI root outputs the data to file
! --------------------------------------------------
IF(.NOT.MPIRoot)RETURN

! --------------------------------------------------
! Output file 1 of 3: MPIW8Time.csv
! --------------------------------------------------
! Either create new file or add info to existing file
WriteHeader = .TRUE.
IF(FILEEXISTS(outfile)) WriteHeader = .FALSE.

IF(WriteHeader)THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=' '
  DO i=1,nTotalVars
    WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO i=1,nTotalVars
    IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF ! WriteHeader

IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
  WRITE(tmpStr2,formatStr)                                  &
      " ",REAL(nProcessors)                                ,& !     'nProcessors'
      delimiter,MPIW8TimeSimeGlobal                        ,& !     'WallTimeSim'
      delimiter,MPIW8TimeGlobal(1)                         ,& !     'Barrier-and-Sync'
      delimiter,REAL(MPIW8CountGlobal(1))                  ,& !     'Barrier-and-Sync-Counter'
      delimiter,MPIW8TimeGlobal(2)                         ,& !     'RAM-Measure-Reduce'
      delimiter,REAL(MPIW8CountGlobal(2))                  ,& !     'RAM-Measure-Reduce-Counter'
      delimiter,MPIW8TimeGlobal(3)                         ,& ! (1) 'HDG-SendLambda'    or 'DGSEM-Send'
      delimiter,REAL(MPIW8CountGlobal(3))                  ,& ! (1) 'HDG-SendLambda-Counter'    or 'DGSEM-Send-Counter'
      delimiter,MPIW8TimeGlobal(4)                         ,& ! (2) 'HDG-ReceiveLambda' or 'DGSEM-Receive'
      delimiter,REAL(MPIW8CountGlobal(4))                   & ! (2) 'HDG-ReceiveLambda-Counter' or 'DGSEM-Receive-Counter'
#if USE_HDG
     ,delimiter,MPIW8TimeGlobal(5)                         ,& ! (3) 'HDG-Broadcast'
      delimiter,REAL(MPIW8CountGlobal(5))                  ,& ! (3) 'HDG-Broadcast-Counter'
      delimiter,MPIW8TimeGlobal(6)                         ,& ! (4) 'HDG-Allreduce'
      delimiter,REAL(MPIW8CountGlobal(6))                   & ! (4) 'HDG-Allreduce-Counter'
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+1)  ,& ! (1) 'SendNbrOfParticles'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+1)) ,& ! (1) 'SendNbrOfParticles-Counter'
      delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+2)  ,& ! (2) 'RecvNbrOfParticles'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+2)) ,& ! (2) 'RecvNbrOfParticles-Counter'
      delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+3)  ,& ! (3) 'SendParticles'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+3)) ,& ! (3) 'SendParticles-Counter'
      delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+4)  ,& ! (4) 'RecvParticles'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+4)) ,& ! (4) 'RecvParticles-Counter'
      delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+5)  ,& ! (5) 'EmissionParticles'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+5)) ,& ! (5) 'EmissionParticles-Counter'
      delimiter,MPIW8TimeGlobal(      MPIW8SIZEFIELD+2+6)  ,& ! (6) 'PIC-depo-Wait'
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+6))  & ! (6) 'PIC-depo-Wait-Counter'
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  WRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write MPI_WAIT wall time info!"
END IF

! --------------------------------------------------
! Output file 2 of 3: MPIW8Time.csv
! --------------------------------------------------
! Either create new file or add info to existing file
WriteHeader = .TRUE.
IF(FILEEXISTS(outfilePerc)) WriteHeader = .FALSE.

IF(WriteHeader)THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfilePerc),STATUS="UNKNOWN")
  tmpStr=""
  DO i=1,nTotalVars
    WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO i=1,nTotalVars
    IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF ! WriteHeader

IF(FILEEXISTS(outfilePerc))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfilePerc),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
  !TotalSimTime = MPIW8TimeSim*REAL(nProcessors)
  TotalSimTime = MPIW8TimeSimeGlobal*0.01         ! Convert to [%]
  TotalCounter = REAL(SUM(MPIW8CountGlobal))*0.01 ! Convert to [%]
  WRITE(tmpStr2,formatStr)                              &
      " ",REAL(nProcessors)                            ,&
      delimiter,100.                                   ,& ! MPIW8TimeSim*nProcessors / TotalSimTime
      delimiter,      MPIW8TimeGlobal(1) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(1))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(2) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(2))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(3) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(3))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(4) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(4))/TotalCounter  &
#if USE_HDG
     ,delimiter,      MPIW8TimeGlobal(5) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(5))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(6) /TotalSimTime ,&
      delimiter,REAL(MPIW8CountGlobal(6))/TotalCounter  &
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+1) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+1))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+2) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+2))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+3) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+3))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+4) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+4))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+5) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+5))/TotalCounter ,&
      delimiter,      MPIW8TimeGlobal(MPIW8SIZEFIELD+2+6) /TotalSimTime        ,&
      delimiter,REAL(MPIW8CountGlobal(MPIW8SIZEFIELD+2+6))/TotalCounter  &
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
  CLOSE(ioUnit)
ELSE
  WRITE(UNIT_StdOut,'(A)')TRIM(outfilePerc)//" does not exist. Cannot write MPI_WAIT wall time info!"
END IF

! --------------------------------------------------
! Output file 3 of 3: MPIW8TimeProc-XXX.csv
! --------------------------------------------------
! Cannot append to proc file, iterate name  -000.csv, -001.csv, -002.csv
WRITE(UNIT=hilf,FMT='(A1,I3.3,A4)') '-',0,'.csv'
outfileProc_loc=TRIM(outfileProc)//'-'//TRIM(ADJUSTL(INTTOSTR(nProcessors)))
DO WHILE(FILEEXISTS(TRIM(outfileProc_loc)//TRIM(hilf)))
  i = i + 1
  WRITE(UNIT=hilf,FMT='(A1,I3.3,A4)') '-',i,'.csv'
END DO
outfileProc_loc=TRIM(outfileProc_loc)//TRIM(hilf)

! Write the file header
OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfileProc_loc),STATUS="UNKNOWN")
tmpStr=' '
DO i=1,nTotalVars
  WRITE(tmpStr(i),'(A)')delimiter//'"'//TRIM(StrVarNames(i))//'"'
END DO
WRITE(formatStr,'(A1)')'('
DO i=1,nTotalVars
  IF(i.EQ.nTotalVars)THEN ! skip writing "," and the end of the line
    WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i))
  ELSE
    WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(i)),','
  END IF
END DO

WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

! Output the processor wait times
WRITE(formatStr,'(A2,I2,A14,A1)')'(',nTotalVars,CSVFORMAT,')'
DO i = 0,nProcessors-1
  WRITE(tmpStr2,formatStr)                                            &
            " ",REAL(i)                                              ,&
      delimiter,MPIW8TimeSim                                         ,&
      delimiter,      MPIW8TimeProc(i*MPIW8SIZE+1)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+1))                  ,&
      delimiter,      MPIW8TimeProc(i*MPIW8SIZE+2)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+2))                  ,&
      delimiter,      MPIW8TimeProc(i*MPIW8SIZE+2)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+2))                  ,&
      delimiter,      MPIW8TimeProc(i*MPIW8SIZE+4)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+4))                   &
#if USE_HDG
     ,delimiter,      MPIW8TimeProc(i*MPIW8SIZE+5)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+5))                  ,&
      delimiter,      MPIW8TimeProc(i*MPIW8SIZE+6)                   ,&
      delimiter,REAL(MPIW8CountProc(i*MPIW8SIZE+6))                   &
#endif /*USE_HDG*/
#if defined(PARTICLES)
     ,delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+1)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+1)) ,&
      delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+2)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+2)) ,&
      delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+3)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+3)) ,&
      delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+4)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+4)) ,&
      delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+5)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+5)) ,&
      delimiter,      MPIW8TimeProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+6)  ,&
      delimiter,REAL(MPIW8CountProc(MPIW8SIZEFIELD+2+i*MPIW8SIZE+6))  &
#endif /*defined(PARTICLES)*/
  ; ! this is required for terminating the "&" when particles=off
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
END DO
CLOSE(ioUnit)

DEALLOCATE(MPIW8TimeProc, MPIW8CountProc)

END SUBROUTINE OutputMPIW8Time
#endif /*USE_MPI && defined(MEASURE_MPI_WAIT)*/

END MODULE MOD_MPI_Output
