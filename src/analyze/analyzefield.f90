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

MODULE MOD_AnalyzeField
!===================================================================================================================================
! Contains the Poynting Vector Integral part for the power analysis of the field vector
!===================================================================================================================================
USE MOD_Globals, ONLY:UNIT_stdout
USE MOD_PreProc
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC:: AnalyzeField
!===================================================================================================================================

CONTAINS

SUBROUTINE AnalyzeField(Time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars          ,ONLY: DoFieldAnalyze,CalcEpot,WEl
USE MOD_Analyze_Vars          ,ONLY: CalcBoundaryFieldOutput,BFO
USE MOD_Mesh_Vars             ,ONLY: BoundaryName
#if (PP_nVar>=6)
USE MOD_Analyze_Vars          ,ONLY: CalcPoyntingInt,nPoyntingIntPlanes,PosPoyntingInt
USE MOD_AnalyzeField_Poynting ,ONLY: CalcPoyntingIntegral
#endif /*PP_nVar>=6*/
#if (PP_nVar==8)
USE MOD_Analyze_Vars          ,ONLY: WMag,WPhi,WPsi
#endif /*PP_nVar=8*/
USE MOD_Particle_Analyze_Vars ,ONLY: IsRestart
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_Dielectric_Vars       ,ONLY: DoDielectric
#if USE_HDG
USE MOD_AnalyzeField_HDG      ,ONLY: CalculateAverageElectricPotential,CalculateElectricDisplacementCurrentSurface
USE MOD_HDG_Vars              ,ONLY: HDGNorm,iterationTotal,RunTimeTotal,UseFPC,FPC,UseEPC,EPC
USE MOD_Analyze_Vars          ,ONLY: AverageElectricPotential,CalcAverageElectricPotential,EDC,CalcElectricTimeDerivative
USE MOD_TimeDisc_Vars         ,ONLY: dt
#endif /*USE_HDG*/
#ifdef PARTICLES
USE MOD_PICInterpolation_Vars ,ONLY: DoInterpolation
#endif /*PARTICLES*/
USE MOD_AnalyzeField_Tools    ,ONLY: CalculateBoundaryFieldOutput
#if !(USE_FV) || (USE_HDG)
USE MOD_AnalyzeField_Tools    ,ONLY: CalcPotentialEnergy,CalcPotentialEnergy_Dielectric
#endif /*no FV alone - !(USE_FV) || (USE_HDG)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL            :: isOpen
CHARACTER(LEN=350) :: outfile
INTEGER            :: unit_index, nOutputVarTotal
#if (PP_nVar>=6)
INTEGER            :: iPlane
REAL               :: PoyntingIntegral(1:nPoyntingIntPlanes)
#endif /*PP_nVar>=6*/
CHARACTER(LEN=1000) :: formatStr
#if (PP_nVar==8)
INTEGER,PARAMETER  :: helpInt=4
#else
INTEGER,PARAMETER  :: helpInt=0
#endif /*PP_nVar=8*/
#if USE_HDG
INTEGER,PARAMETER  :: helpInt2=4
INTEGER            :: iEDCBC,iUniqueFPCBC,iUniqueEPCBC
#else
INTEGER,PARAMETER  :: helpInt2=0
#endif /*USE_HDG*/
INTEGER,PARAMETER  :: nOutputVar=2+helpInt+helpInt2
CHARACTER(LEN=255),DIMENSION(nOutputVar) :: StrVarNames(nOutputVar)=(/ CHARACTER(LEN=255) :: &
    'time', &
    'E-El'  &
#if (PP_nVar==8)
   ,'E-Mag', &
    'E-phi', &
    'E-psi', &
    'E-pot'  &
#endif /*PP_nVar=8*/
#if USE_HDG
   ,'HDG-#iterations', &
    'HDG-RunTime', &
    'HDG-RunTimePerIteration', &
    'HDG-Norm' &
#endif /*USE_HDG*/
    /)
CHARACTER(LEN=500),ALLOCATABLE :: tmpStr(:) ! needed because PerformAnalyze is called multiple times at the beginning
CHARACTER(LEN=5000)            :: tmpStr2
CHARACTER(LEN=1),PARAMETER     :: delimiter=","
INTEGER                        :: I,iBoundary
CHARACTER(LEN=255) :: StrVarNameTmp
REAL               :: BoundaryFieldOutput(1:PP_nVar)
!===================================================================================================================================
IF ( DoRestart ) THEN
  isRestart = .true.
END IF
IF (.NOT.DoFieldAnalyze) RETURN
unit_index = 537
#if USE_MPI
IF(MPIROOT)THEN
#endif /*USE_MPI*/
  INQUIRE(UNIT   = unit_index , OPENED = isOpen)
  IF (.NOT.isOpen) THEN
    outfile = 'FieldAnalyze.csv'
    IF (isRestart .and. FILEEXISTS(outfile)) THEN
       OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
       !CALL FLUSH (unit_index)
    ELSE
      OPEN(unit_index,file=TRIM(outfile))
      !CALL FLUSH (unit_index)
      !--- insert header

      ! Set the header line content
#if (PP_nVar>=6)
      nPoyntingIntPlanes=MERGE(nPoyntingIntPlanes,0,CalcPoyntingInt) ! set to zero if the flag is false (otherwise not initialized)
      nOutputVarTotal = nOutputVar + nPoyntingIntPlanes
#else
      nOutputVarTotal = nOutputVar
#endif /*PP_nVar>=6*/
#if USE_HDG
      ! Add averaged electric field (integrated and averaged along y-z-direction)
      IF(CalcAverageElectricPotential) nOutputVarTotal = nOutputVarTotal + 1
      !-- Electric displacement current
      IF(CalcElectricTimeDerivative) nOutputVarTotal = nOutputVarTotal + EDC%NBoundaries
      !-- Floating boundary condition
      IF(UseFPC) nOutputVarTotal = nOutputVarTotal + 2*FPC%nUniqueFPCBounds ! Charge and Voltage on each FPC
      !-- Electric potential condition
      IF(UseEPC) nOutputVarTotal = nOutputVarTotal + 2*EPC%nUniqueEPCBounds ! Current and Voltage on each EPC
#endif /*USE_HDG*/
#if (PP_nVar==8)
      IF(.NOT.CalcEpot) nOutputVarTotal = nOutputVarTotal - 5
#else
      IF(.NOT.CalcEpot) nOutputVarTotal = nOutputVarTotal - 1
#endif /*PP_nVar=8*/
      IF(CalcBoundaryFieldOutput) nOutputVarTotal = nOutputVarTotal + BFO%NFieldBoundaries
      ALLOCATE(tmpStr(1:nOutputVarTotal))
      tmpStr=""

      nOutputVarTotal = 0
      DO I=1,nOutputVar
        ! When NOT CalcEpot, skip entries 2,...,6
#if (PP_nVar==8)
        IF((.NOT.CalcEpot).AND.((1.LT.I).AND.(I.LE.6))) CYCLE
#else
        IF((.NOT.CalcEpot).AND.(I.EQ.2)) CYCLE
#endif /*PP_nVar=8*/
        nOutputVarTotal = nOutputVarTotal + 1
        WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNames(I))//'"'
      END DO

#if (PP_nVar>=6)
      ! Add Poynting vector integral (integrated energy density through plane)
      IF(CalcPoyntingInt)THEN
        DO iPlane=1,nPoyntingIntPlanes
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3,A1,E23.16E3,A1)') 'Plane-Pos-',iPlane,'(',PosPoyntingInt(iPlane),')'
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
        END DO
      END IF
#endif /*PP_nVar>=6*/

#if USE_HDG
      ! Add averaged electric field (integrated and averaged along y-z-direction)
      IF(CalcAverageElectricPotential)THEN
        nOutputVarTotal = nOutputVarTotal + 1
        WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-AverageElectricPotential"'
      END IF ! CalcAverageElectricPotential

      !-- Electric displacement current
      IF(CalcElectricTimeDerivative)THEN
        DO iEDCBC = 1, EDC%NBoundaries
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3,A)') 'ElecDisplCurrent-',iEDCBC,'-'//TRIM(BoundaryName(EDC%FieldBoundaries(iEDCBC)))
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
        END DO ! iEDCBC = 1, EDC%NBoundaries
      END IF

      !-- Floating boundary condition
      IF(UseFPC)THEN
        DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3)') 'FPC-Charge-BCState-',FPC%BCState(iUniqueFPCBC)
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3)') 'FPC-Voltage-BCState-',FPC%BCState(iUniqueFPCBC)
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
        END DO ! iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      END IF

      !-- Electric potential condition
      IF(UseEPC)THEN
        DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3)') 'EPC-Current-BCState-',EPC%BCState(iUniqueEPCBC)
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3)') 'EPC-Voltage-BCState-',EPC%BCState(iUniqueEPCBC)
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'"'
        END DO ! iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      END IF
#endif /*USE_HDG*/

      ! Add BoundaryFieldOutput for each boundary that is required
      IF(CalcBoundaryFieldOutput)THEN
        DO iBoundary=1,BFO%NFieldBoundaries
          nOutputVarTotal = nOutputVarTotal + 1
          WRITE(StrVarNameTmp,'(A,I0.3)') 'BFO-boundary-',iBoundary
          WRITE(tmpStr(nOutputVarTotal),'(A,I0.3,A)')delimiter//'"',nOutputVarTotal,'-'//TRIM(StrVarNameTmp)//'-'//&
              TRIM(BoundaryName(BFO%FieldBoundaries(iBoundary)))//'"'
        END DO
      END IF

      ! Set the format
      WRITE(formatStr,'(A1)')'('
      DO I=1,nOutputVarTotal
        WRITE(formatStr,'(A,A1,A1,I2)')TRIM(formatStr),',','A',LEN_TRIM(tmpStr(I))
      END DO
      formatStr(2:2) = ' ' ! remove comma
      WRITE(formatStr,'(A,A1)')TRIM(formatStr),')'      ! finish the format

      WRITE(tmpStr2,formatStr)tmpStr(1:nOutputVarTotal) ! use the format and write the header names to a temporary string
      tmpStr2(1:1) = " "                                ! remove possible delimiter at the beginning (e.g. a comma)
      WRITE(unit_index,'(A)')TRIM(ADJUSTL(tmpStr2))     ! clip away the front and rear white spaces of the temporary string
    END IF
  END IF
#if USE_MPI
END IF
#endif /*USE_MPI*/

IF(CalcEpot)THEN
  ! energy of
  ! 1) electric field
  ! 2) magnetic field
  ! 3) divergence correction magnetic
  ! 4) divergence correction electric + charge
#if (PP_nVar==8)
  IF(DoDielectric)THEN
    CALL CalcPotentialEnergy_Dielectric(WEl,WMag,Wphi,Wpsi)
  ELSE
    CALL CalcPotentialEnergy(WEl,WMag,Wphi,Wpsi)
  END IF ! DoDielectric
#else
  IF(DoDielectric)THEN
    CALL CalcPotentialEnergy_Dielectric(WEl)
  ELSE
    CALL CalcPotentialEnergy(WEl)
  END IF ! DoDielectric
#endif /*PP_nVar=8*/
END IF ! CalcEpot
#if (PP_nVar>=6)
IF(CalcPoyntingInt) CALL CalcPoyntingIntegral(PoyntingIntegral)
#endif /*PP_nVar>=6*/

#ifdef PARTICLES
IF(.NOT.DoInterpolation)THEN
#endif /*PARTICLES*/
#if USE_HDG
  !1.2 Calculate external E-field
  IF(CalcAverageElectricPotential) CALL CalculateAverageElectricPotential()
#endif /*USE_HDG*/
#ifdef PARTICLES
END IF ! .NOT.DoInterpolation
#endif /*PARTICLES*/
#if USE_HDG
IF(CalcElectricTimeDerivative) CALL CalculateElectricDisplacementCurrentSurface()
#endif /*USE_HDG*/

IF(MPIROOT)THEN
  WRITE(unit_index,'(E23.16E3)',ADVANCE='NO') Time
  IF (CalcEpot) THEN
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', WEl
#if (PP_nVar==8)
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', WMag
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', WPhi
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', WPsi
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', WEl + WMag + WPhi + WPsi
#endif /*PP_nVar=8*/
  END IF
#if USE_HDG
  WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', REAL(iterationTotal)
  WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', RunTimeTotal
  IF(iterationTotal.GT.0)THEN
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', RunTimeTotal/REAL(iterationTotal)
  ELSE
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', 0.
  END IF ! iterationTotal.GT.0
  WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',', HDGNorm
#endif /*USE_HDG*/
#if (PP_nVar>=6)
  ! Add Poynting vector integral (integrated energy density through plane)
  IF(CalcPoyntingInt)THEN
    DO iPlane=1,nPoyntingIntPlanes
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',PoyntingIntegral(iPlane)
    END DO
  END IF
#endif /*PP_nVar>=6*/
#if USE_HDG
  ! Add averaged electric field (integrated and averaged along y-z-direction)
  IF(CalcAverageElectricPotential)THEN
    WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',AverageElectricPotential
  END IF ! CalcAverageElectricPotential

  !-- Electric displacement current
  IF(CalcElectricTimeDerivative)THEN
    DO iEDCBC = 1, EDC%NBoundaries
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',EDC%Current(iEDCBC)
    END DO ! iEDCBC = 1, EDC%NBoundaries
  END IF

  !-- Floating boundary condition
  IF(UseFPC)THEN
    DO iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',FPC%Charge(iUniqueFPCBC)
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',FPC%Voltage(iUniqueFPCBC)
    END DO !iUniqueFPCBC = 1, FPC%nUniqueFPCBounds
  END IF

  !-- Electric potential condition
  IF(UseEPC)THEN
    DO iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',-EPC%Charge(iUniqueEPCBC)/dt
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',EPC%Voltage(iUniqueEPCBC)
    END DO !iUniqueEPCBC = 1, EPC%nUniqueEPCBounds
  END IF
#endif /*USE_HDG*/
  ! ! Add BoundaryFieldOutput for each boundary that is required
  IF(CalcBoundaryFieldOutput)THEN
    DO iBoundary=1,BFO%NFieldBoundaries
      CALL CalculateBoundaryFieldOutput(BFO%FieldBoundaries(iBoundary),Time,BoundaryFieldOutput)
      WRITE(unit_index,CSVFORMAT,ADVANCE='NO') ',',BoundaryFieldOutput
    END DO
  END IF ! CalcBoundaryFieldOutput
  write(unit_index,'(A)') '' ! write 'newline' to file to finish the current line
END IF

END SUBROUTINE AnalyzeField


END MODULE MOD_AnalyzeField