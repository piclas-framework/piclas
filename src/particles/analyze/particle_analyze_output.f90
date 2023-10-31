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

MODULE MOD_Particle_Analyze_Output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
#ifdef PARTICLES
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: WriteParticleTrackingData
PUBLIC :: DisplayCoupledPowerPart
!===================================================================================================================================

CONTAINS


!----------------------------------------------------------------------------------------------------------------------------------!
!> Write particle info to ParticlePosition.csv file
!> time, pos, velocity, gamma, element
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE WriteParticleTrackingData(time,iter)
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals               ,ONLY: FILEEXISTS,unit_stdout,DOTPRODUCT,abort,MPI_COMM_PICLAS,MPIRoot
#if USE_MPI
USE MOD_Globals               ,ONLY: iError
#endif /*USE_MPI*/
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_Particle_Vars         ,ONLY: PartState, PDM, PEM
USE MOD_Particle_Analyze_Vars ,ONLY: printDiff,printDiffVec,printDiffTime
USE MOD_part_tools            ,ONLY: UpdateNextFreePosition
USE MOD_Globals_Vars          ,ONLY: c2_inv
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                  :: time
INTEGER(KIND=8),INTENT(IN)       :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=20),PARAMETER              :: outfile='ParticlePosition.csv'
INTEGER                                  :: ioUnit,I
CHARACTER(LEN=150)                       :: formatStr
INTEGER,PARAMETER                        :: nOutputVar=10
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    '001-time',     &
    'PartNum',  &
    'PartPosX', &
    'PartPosY', &
    'PartPosZ', &
    'PartVelX', &
    'PartVelY', &
    'PartVelZ', &
    'gamma',    &
    'Element'/)
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: tmpStr ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=1000)                      :: tmpStr2
CHARACTER(LEN=1),PARAMETER               :: delimiter=","
LOGICAL                                  :: FileExist,CreateFile
REAL                                     :: diffPos,diffVelo
INTEGER                                  :: iPartState
REAL                                     :: gamma1
!===================================================================================================================================
! only the root shall write this file
!IF(.NOT.MPIRoot)RETURN

! check if file is to be created
CreateFile=.TRUE.
IF(iter.GT.0)CreateFile=.FALSE.                             ! don't create new file if this is not the first iteration
IF((DoRestart).AND.(FILEEXISTS(outfile)))CreateFile=.FALSE. ! don't create new file if this is a restart and the file already exists
!                                                           ! assume continued simulation and old load balance data is still needed

! check if new file with header is to be created
INQUIRE(FILE = outfile, EXIST=FileExist)
IF(.NOT.FileExist)CreateFile=.TRUE.                         ! if no file exists, create one

! create file with header
IF(MPIRoot.AND.CreateFile) THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),STATUS="UNKNOWN")
  tmpStr=""
  DO I=1,nOutputVar
    WRITE(tmpStr(I),'(A)')delimiter//'"'//TRIM(StrVarNames(I))//'"'
  END DO
  WRITE(formatStr,'(A1)')'('
  DO I=1,nOutputVar
    IF(I.EQ.nOutputVar)THEN ! skip writing "," and the end of the line
      WRITE(formatStr,'(A,A1,I2)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I))
    ELSE
      WRITE(formatStr,'(A,A1,I2,A1)')TRIM(formatStr),'A',LEN_TRIM(tmpStr(I)),','
    END IF
  END DO

  WRITE(formatStr,'(A,A1)')TRIM(formatStr),')' ! finish the format
  WRITE(tmpStr2,formatStr)tmpStr               ! use the format and write the header names to a temporary string
  tmpStr2(1:1) = " "                           ! remove possible delimiter at the beginning (e.g. a comma)
  WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2))    ! clip away the front and rear white spaces of the temporary string

  CLOSE(ioUnit)
END IF

#if USE_MPI
! Barrier is required is root creates file and other processor prints to this file
IF(CreateFile) CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

CALL UpdateNextFreePosition()

! Print info to file
IF(FILEEXISTS(outfile))THEN
  OPEN(NEWUNIT=ioUnit,FILE=TRIM(outfile),POSITION="APPEND",STATUS="OLD")
  WRITE(formatStr,'(A2,I2,A14,A1)')'(',nOutputVar,CSVFORMAT,')'
  DO i=1,PDM%ParticleVecLength

    ! Calculate Lorentz factor
    gamma1 = DOTPRODUCT(PartState(4:6,i))*c2_inv
    ! Sanity check: Lorentz factor must be below 1.0
    IF(gamma1.GE.1.0)THEN
      gamma1=-1.0
    ELSE
      gamma1=1.0/SQRT(1.-gamma1)
    END IF

    IF (PDM%ParticleInside(i)) THEN
      WRITE(tmpStr2,formatStr)&
          " ",time, &                 ! time
          delimiter,REAL(i), &        ! PartNum
          delimiter,PartState(1,i), & ! PartPosX
          delimiter,PartState(2,i), & ! PartPosY
          delimiter,PartState(3,i), & ! PartPosZ
          delimiter,PartState(4,i), & ! PartVelX
          delimiter,PartState(5,i), & ! PartVelY
          delimiter,PartState(6,i), & ! PartVelZ
          delimiter,        gamma1, & ! gamma
          delimiter,REAL(PEM%GlobalElemID(i))                                                  ! Element
      WRITE(ioUnit,'(A)')TRIM(ADJUSTL(tmpStr2)) ! clip away the front and rear white spaces of the data line
    END IF
  END DO
  CLOSE(ioUnit)
ELSE
  SWRITE(UNIT_StdOut,'(A)')TRIM(outfile)//" does not exist. Cannot write particle tracking info!"
END IF

! printDiff
IF (printDiff) THEN
  diffPos=0.
  diffVelo=0.
  IF (time.GE.printDiffTime) THEN
    printDiff=.FALSE.
    DO iPartState=1,3
      diffPos=diffPos+(printDiffVec(iPartState)-PartState(iPartState,1))**2
      diffVelo=diffVelo+(printDiffVec(iPartState+3)-PartState(iPartState+3,1))**2
    END DO
    WRITE(*,'(A,e24.14,1X,e24.14)') 'L2-norm from printDiffVec: ',SQRT(diffPos),SQRT(diffVelo)
  END IF
END IF
END SUBROUTINE WriteParticleTrackingData


SUBROUTINE DisplayCoupledPowerPart()
!===================================================================================================================================
!> Print accumulated power transferred to particles to std out (power for each species and total power over all particles)
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: Time
USE MOD_Restart_Vars          ,ONLY: RestartTime
USE MOD_Globals               ,ONLY: abort,mpiroot
USE MOD_Particle_Analyze_Vars ,ONLY: PCouplSpec
USE MOD_Particle_Vars         ,ONLY: nSpecies,Species
USE MOD_Mesh_Vars             ,ONLY: nElems, offSetElem
#if USE_MPI
USE MOD_Globals
#endif /*USE_MPI*/
USE MOD_Globals               ,ONLY: UNIT_StdOut
USE MOD_Particle_Mesh_Vars    ,ONLY: ElemVolume_Shared
USE MOD_Mesh_Tools            ,ONLY: GetCNElemID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL          :: timediff,PTotal(1:nSpecies),SumPTotal(1:nSpecies)
INTEGER       :: iSpec,iElem
CHARACTER(5)  :: hilf
!===================================================================================================================================

IF(ABS(Time-RestartTime).LE.0.0) RETURN

timediff = 1.0 / (Time-RestartTime)
! Sanity check: integrate cell-averaged PCoupl
PTotal = 0.

DO iSpec = 1, nSpecies
  IF(ABS(Species(iSpec)%ChargeIC).GT.0.0)THEN
    DO iElem = 1, nElems
      PTotal(iSpec) = PTotal(iSpec) + PCouplSpec(iSpec)%DensityAvgElem(iElem) * ElemVolume_Shared(GetCNElemID(iElem+offSetElem))
    END DO ! iElem = 1, nElems
  END IF ! ABS(Species(iSpec)%ChargeIC).GT.0.0)
END DO ! iSpec = 1, nSpecies

! Sum the power
#if USE_MPI
CALL MPI_REDUCE(PTotal(1:nSpecies),SumPTotal(1:nSpecies),nSpecies,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS,iError)
#else
SumPTotal(1:nSpecies)=PTotal(1:nSpecies)
#endif
IF(mpiroot)THEN
  SumPTotal = SumPTotal * timediff
  WRITE(UNIT_StdOut,*) "    Averaged coupled power per species [W]"
  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilf,FMT='(I0)') iSpec
    WRITE (UNIT_StdOut,*) "    "//hilf//" : ",SumPTotal(iSpec)
  END DO ! iSpec = 1, nSpecies
  WRITE (UNIT_StdOut,*) "    Total : ",SUM(SumPTotal)
END IF ! mpiroot

END SUBROUTINE DisplayCoupledPowerPart



#endif /*PARTICLES*/
END MODULE MOD_Particle_Analyze_Output
