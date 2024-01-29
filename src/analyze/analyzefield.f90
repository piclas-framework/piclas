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

#if (PP_nVar>=6)
INTERFACE CalcPoyntingIntegral
  MODULE PROCEDURE CalcPoyntingIntegral
END INTERFACE

PUBLIC:: CalcPoyntingIntegral,GetPoyntingIntPlane,FinalizePoyntingInt
#endif

INTERFACE CalcPotentialEnergy
  MODULE PROCEDURE CalcPotentialEnergy
END INTERFACE

INTERFACE CalcPotentialEnergy_Dielectric
  MODULE PROCEDURE CalcPotentialEnergy_Dielectric
END INTERFACE

PUBLIC:: CalcPotentialEnergy,CalcPotentialEnergy_Dielectric
PUBLIC:: AnalyzeField
#if USE_HDG
PUBLIC:: GetAverageElectricPotentialPlane,CalculateAverageElectricPotential,FinalizeAverageElectricPotential
#endif /*USE_HDG*/
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
#endif /*PP_nVar>=6*/
#if (PP_nVar==8)
USE MOD_Analyze_Vars          ,ONLY: WMag,WPhi,WPsi
#endif /*PP_nVar=8*/
USE MOD_Particle_Analyze_Vars ,ONLY: IsRestart
USE MOD_Restart_Vars          ,ONLY: DoRestart
USE MOD_Dielectric_Vars       ,ONLY: DoDielectric
#if USE_HDG
USE MOD_HDG_Vars              ,ONLY: HDGNorm,iterationTotal,RunTimeTotal,UseFPC,FPC,UseEPC,EPC
USE MOD_Analyze_Vars          ,ONLY: AverageElectricPotential,CalcAverageElectricPotential,EDC,CalcElectricTimeDerivative
USE MOD_TimeDisc_Vars         ,ONLY: dt
#endif /*USE_HDG*/
#ifdef PARTICLES
USE MOD_PICInterpolation_Vars ,ONLY: DoInterpolation
#endif /*PARTICLES*/
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
  END IF
#else
  IF(DoDielectric)THEN
    CALL CalcPotentialEnergy_Dielectric(WEl)
  ELSE
    CALL CalcPotentialEnergy(WEl)
  END IF
#endif /*PP_nVar=8*/
END IF
#if (PP_nVar>=6)
IF(CalcPoyntingInt) CALL CalcPoyntingIntegral(PoyntingIntegral,doProlong=.TRUE.)
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

#if (PP_nVar>=6)
SUBROUTINE CalcPoyntingIntegral(PoyntingIntegral,doProlong)
!===================================================================================================================================
! Calculation of Poynting Integral with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is missing ... ups
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems, SurfElem, NormVec
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_Analyze_Vars       ,ONLY: nPoyntingIntPlanes,S,isPoyntingIntSide,SideIDToPoyntingSide,PoyntingMainDir
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus,wGPSurf
USE MOD_DG_Vars            ,ONLY: U,U_master
USE MOD_Globals_Vars       ,ONLY: smu0
USE MOD_Dielectric_Vars    ,ONLY: isDielectricFace,PoyntingUseMuR_Inv,Dielectric_MuR_Master_inv,DoDielectric
#if USE_MPI
USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN),OPTIONAL :: doProlong
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: PoyntingIntegral(1:nPoyntingIntPlanes)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem, SideID,ilocSide,iPoyntingSide
INTEGER          :: p,q,l
REAL             :: Uface(PP_nVar,0:PP_N,0:PP_N)
REAL             :: SIP(0:PP_N,0:PP_N)
#if USE_MPI
REAL             :: SumSabs(nPoyntingIntPlanes)
#endif
LOGICAL          :: Prolong=.TRUE.
!REAL             :: sresvac
!===================================================================================================================================

IF(PRESENT(doProlong))THEN
  Prolong=doProlong
ELSE
  Prolong=.TRUE.
ENDIF

S    = 0.
PoyntingIntegral = 0.

iPoyntingSide = 0 ! only if all Poynting vectors are desired
DO iELEM = 1, nElems
  Do ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isPoyntingIntSide(SideID)) CYCLE
      IF(Prolong)THEN
#if (PP_NodeType==1) /* for Gauss-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,l,p,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,0,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,l,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,p,q,l,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! qfirst stuff
        CASE(XI_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,0,p,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,l,p,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,0,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+U(:,p,l,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,0,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,q,l,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        END SELECT
#else /* for Gauss-Lobatto-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          Uface(:,:,:)=U(:,:,0,:,iElem)
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
          Uface(:,:,:)=U(:,PP_N,:,:,iElem)
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,PP_N,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,PP_N,iElem)
            END DO ! p
          END DO ! q
        END SELECT
#endif
        ELSE ! no prolonge to face
          Uface=U_master(:,:,:,SideID)
        END IF ! Prolong
        ! calculate Poynting vector
        iPoyntingSide = iPoyntingSide + 1

        ! check if dielectric regions are involved
        IF(DoDielectric)THEN
          IF(PoyntingUseMuR_Inv.AND.isDielectricFace(SideID))THEN
            CALL PoyntingVectorDielectric(Uface(:,:,:),S(:,:,:,iPoyntingSide),Dielectric_MuR_Master_inv(0:PP_N,0:PP_N,SideID))
          ELSE
            CALL PoyntingVector(Uface(:,:,:),S(:,:,:,iPoyntingSide))
          END IF
        ELSE
          CALL PoyntingVector(Uface(:,:,:),S(:,:,:,iPoyntingSide))
        END IF

        IF ( NormVec(PoyntingMainDir,0,0,SideID) .GT. 0.0 ) THEN
          SIP(:,:) =   S(1,:,:,iPoyntingSide) * NormVec(1,:,:,SideID) &
                     + S(2,:,:,iPoyntingSide) * NormVec(2,:,:,SideID) &
                     + S(3,:,:,iPoyntingSide) * NormVec(3,:,:,SideID)
        ELSE ! NormVec(PoyntingMainDir,:,:,iPoyningSide) < 0
          SIP(:,:) = - S(1,:,:,iPoyntingSide) * NormVec(1,:,:,SideID) &
                     - S(2,:,:,iPoyntingSide) * NormVec(2,:,:,SideID) &
                     - S(3,:,:,iPoyntingSide) * NormVec(3,:,:,SideID)
        END IF ! NormVec(PoyntingMainDir,:,:,iPoyntingSide)
        ! multiplied by surface element and  Gauss Points
        SIP(:,:) = SIP(:,:) * SurfElem(:,:,SideID) * wGPSurf(:,:)

        ! total flux through each plane
        PoyntingIntegral(SideIDToPoyntingSide(SideID)) = PoyntingIntegral(SideIDToPoyntingSide(SideID)) + smu0* SUM(SIP(:,:))
    END IF ! flip =0
  END DO ! iSides
END DO ! iElems

#if USE_MPI
  CALL MPI_REDUCE(PoyntingIntegral(:) , SumSabs(:) , nPoyntingIntPlanes , MPI_DOUBLE_PRECISION ,MPI_SUM, 0, MPI_COMM_PICLAS,IERROR)
  PoyntingIntegral(:) = SumSabs(:)
#endif /*USE_MPI*/

END SUBROUTINE CalcPoyntingIntegral
#endif


#if (PP_nVar>=6)
PPURE SUBROUTINE PoyntingVector(Uface_in,Sloc)
!===================================================================================================================================
!> Calculate the Poynting Vector on a certain face for vacuum properties
!>
!> ATTENTION: permeability is not applied here due to performance gain
!> Definition: S = E x H = 1/mu0 * ( E x H )
!> Here      : S = E x B (i.e. mu0 is applied later)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: Uface_in(PP_nVar,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: Sloc(1:3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
!===================================================================================================================================

! calculate the Poynting vector at each node, additionally the abs of the Poynting vector only based on E
DO p = 0,PP_N
  DO q = 0,PP_N
    Sloc(1,p,q)  =  Uface_in(2,p,q)*Uface_in(6,p,q) - Uface_in(3,p,q)*Uface_in(5,p,q)
    Sloc(2,p,q)  = -Uface_in(1,p,q)*Uface_in(6,p,q) + Uface_in(3,p,q)*Uface_in(4,p,q)
    Sloc(3,p,q)  =  Uface_in(1,p,q)*Uface_in(5,p,q) - Uface_in(2,p,q)*Uface_in(4,p,q)
  END DO ! q - PP_N
END DO  ! p - PP_N

END SUBROUTINE PoyntingVector


PPURE SUBROUTINE PoyntingVectorDielectric(Uface_in,Sloc,mu_r_inv)
!===================================================================================================================================
!> Calculate the Poynting Vector on a certain face for dielectric properties (consider mu_r here, but not mu0)
!>
!> ATTENTION: permeability is not applied here due to performance gain
!> Definition: S = E x H = 1/(mu_r*mu_0) * ( E x H )
!> Here      : S = 1/mu_r * E x B (i.e. mu0 is applied later)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: Uface_in(PP_nVar,0:PP_N,0:PP_N)
REAL,INTENT(IN)       :: mu_r_inv(0:PP_N,0:PP_N)         ! 1/mu_r for every face DOF (may vary on face depending on position)
!                                                        ! (isotropic property for permittivity)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: Sloc(1:3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
!===================================================================================================================================

! calculate the Poynting vector at each node, additionally the abs of the Poynting vector only based on E
DO p = 0,PP_N
  DO q = 0,PP_N
    Sloc(1,p,q)  = (  Uface_in(2,p,q)*Uface_in(6,p,q) - Uface_in(3,p,q)*Uface_in(5,p,q) ) * mu_r_inv(p,q)
    Sloc(2,p,q)  = ( -Uface_in(1,p,q)*Uface_in(6,p,q) + Uface_in(3,p,q)*Uface_in(4,p,q) ) * mu_r_inv(p,q)
    Sloc(3,p,q)  = (  Uface_in(1,p,q)*Uface_in(5,p,q) - Uface_in(2,p,q)*Uface_in(4,p,q) ) * mu_r_inv(p,q)
  END DO ! q - PP_N
END DO  ! p - PP_N

END SUBROUTINE PoyntingVectorDielectric


SUBROUTINE GetPoyntingIntPlane()
!===================================================================================================================================
!> Initializes Poynting vector integral variables and check every side: set "isPoyntingIntSide(SideID) = .TRUE." if a side coincides
!> with a defined Poynting vector integral plane.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars       ,ONLY: nSides,nElems,Face_xGP
USE MOD_Mesh_Vars       ,ONLY: ElemToSide,normvec
USE MOD_Analyze_Vars    ,ONLY: PoyntingIntCoordErr,nPoyntingIntPlanes,PosPoyntingInt,S,STEM
USE MOD_Analyze_Vars    ,ONLY: isPoyntingIntSide,SideIDToPoyntingSide,PoyntingMainDir
USE MOD_ReadInTools     ,ONLY: GETINT,GETREAL
USE MOD_Dielectric_Vars ,ONLY: DoDielectric,nDielectricElems,DielectricMu,ElemToDielectric,isDielectricInterFace
USE MOD_Dielectric_Vars ,ONLY: isDielectricFace,PoyntingUseMuR_Inv
USE MOD_Globals         ,ONLY: abort
#if USE_MPI
USE MOD_Globals
#else
USE MOD_Globals         ,ONLY: CollectiveStop
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iSide, iPlane, SideID
INTEGER,ALLOCATABLE :: nFaces(:)
REAL                :: diff
INTEGER             :: p,q
CHARACTER(LEN=32)   :: index_plane
INTEGER,ALLOCATABLE :: sumFaces(:)
INTEGER             :: sumAllfaces
LOGICAL             :: CheckDielectricSides
INTEGER             :: PoyntingNormalDir1,PoyntingNormalDir2
INTEGER             :: nPoyntingIntSides    !< Sides for the calculation of the Poynting vector integral
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' GET PLANES TO CALCULATE POYNTING VECTOR INTEGRAL ...'

! Initialize number of Poynting plane sides zero and set all sides to false
nPoyntingIntSides=0
ALLOCATE(isPoyntingIntSide(1:nSides))
isPoyntingIntSide = .FALSE.

! Get the number of Poynting planes and coordinates
nPoyntingIntPlanes = GETINT('PoyntingVecInt-Planes')
PoyntingMainDir = GETINT('PoyntingMainDir') ! default "3" is z-direction
SELECT CASE (PoyntingMainDir)
  CASE (1) ! poynting vector integral in x-direction
    PoyntingNormalDir1=2
    PoyntingNormalDir2=3
  CASE (2) ! poynting vector integral in y-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=3
  CASE (3) ! poynting vector integral in z-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=2
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Poynting vector itnegral currently only in x,y,z!')
END SELECT
ALLOCATE(PosPoyntingInt(nPoyntingIntPlanes))
ALLOCATE(SideIDToPoyntingSide(nSides))
ALLOCATE(nFaces(nPoyntingIntPlanes))
SideIDToPoyntingSide = -1
nFaces(:) = 0

! Get z-coordinates and factors for every Poynting plane
DO iPlane=1,nPoyntingIntPlanes
 WRITE(UNIT=index_plane,FMT='(I2.2)') iPlane
 SELECT CASE (PoyntingMainDir)
    CASE (1)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-x-coord')
    CASE (2)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-y-coord')
    CASE (3)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-z-coord')
  END SELECT
END DO
PoyntingIntCoordErr=GETREAL('Plane-Tolerance')

! Dielectric Sides:
! 1.) check if a dielectric region (only permeability, NOT permittivity is important) coincides with a Poynting vector
!     integral plane. Dielectric interfaces with mu_r .NE. 1.0 cannot compute a Poynting vector because of the jump in material
!     parameter of mu_r
CheckDielectricSides=.FALSE.
IF(DoDielectric)THEN
  IF(ANY(ABS(DielectricMu(:,:,:,1:nDielectricElems)-1.0).GT.0.0))THEN
    CheckDielectricSides=.TRUE.
  END IF
END IF

! 2.) for dielectric sides (NOT interface sides between dielectric and some other region), determine mu_r on face for Poynting vector
PoyntingUseMuR_Inv=.FALSE.

! Loop over all planes
DO iPlane = 1, nPoyntingIntPlanes
  ! Loop over all elements
  DO iElem=1,nElems
    ! Loop over all local sides
    DO iSide=1,6
      IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
        SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
        ! First search only planes with normal vector parallel to direction of "MainDir"
        IF((     NormVec(PoyntingNormalDir1,0,0,SideID)  < PoyntingIntCoordErr) .AND. &
           (     NormVec(PoyntingNormalDir2,0,0,SideID)  < PoyntingIntCoordErr) .AND. &
           ( ABS(NormVec(PoyntingMainDir   ,0,0,SideID)) > PoyntingIntCoordErr))THEN
        ! Loop over all Points on Face
          DO q=0,PP_N
            DO p=0,PP_N
              diff = ABS(Face_xGP(PoyntingMainDir,p,q,SideID) - PosPoyntingInt(iPlane))
              IF (diff < PoyntingIntCoordErr) THEN
                IF (.NOT.isPoyntingIntSide(SideID)) THEN
                  nPoyntingIntSides = nPoyntingIntSides +1
                  SideIDToPoyntingSide(SideID) = iPlane
                  isPoyntingIntSide(SideID) = .TRUE.
                  nFaces(iPlane) = nFaces(iPlane) + 1

                  ! Dielectric sides
                  IF(CheckDielectricSides)THEN
                    ! 1.) Check for illegal sides in dielectrics: mu_r != 1.0 on dielectric interface
                    IF(isDielectricInterFace(SideID))THEN
                      IF(ANY(ABS(DielectricMu(:,:,:,ElemToDielectric(iElem))-1.0).GT.0.0))THEN
                        ! If the Poynting vector integral SideID additionally is a dielectric interface between a dielectric region
                        ! with a permittivity and vacuum, then mu_r might be unequal to 1.0 on the interface and the calculation of
                        ! the Poynting vector is not implemented for this case
                        IPWRITE(UNIT_stdOut,*) " "
                        IPWRITE(UNIT_stdOut,*) "Found illegal Poyting plane side. SideID= ",SideID,&
                            " z-coordinate= ",PosPoyntingInt(iPlane)
                        CALL abort(__STAMP__&
                            ,'GetPoyntingIntPlane: Found SideID for Poynting vector integral which is attached to an element'//&
                            ' within which the dielectric permittivity mu_r is not euqal to 1.0 everywhere. The value could be'//&
                            ' unequal to 1.0 on the interface and this is not implemented. TODO: determine mu_r on interface,'//&
                            ' communicate it via MPI (do not forget Mortar sides) and calculate the Poynting vector on that'//&
                            ' interface via some method.')
                      END IF
                    END IF

                    ! 2.) Check for legal sides in dielectrics: mu_r != 1.0 within dielectric region
                    IF(isDielectricFace(SideID))THEN
                      !IPWRITE(UNIT_stdOut,*) "found dielectric face: ",SideID,"z= ",PosPoyntingInt(iPlane)
                      PoyntingUseMuR_Inv=.TRUE.
                    END IF
                  END IF

                END IF
              END IF ! diff < eps
            END DO !p
          END DO !q
        END IF ! n parallel gyrotron axis
      END IF ! flip = 0 master side
    END DO ! iSides
  END DO !iElem=1,nElems
END DO ! iPlanes

! Dielectric sides:
#if USE_MPI
! Send info to ALL MPI ranks:
! TODO: If 1/mu_r is never needed on master AND slave procs, this routine can be adjusted so that only master procs determine the
! prolonged values of mu_r and no MPI information has to be sent. The master side cannot currently be outside of the dielectric
! region (e.g. in vacuum) because that is not allowed. If this would be allowed that MPI rank would need the information of the
! prolonged dielectric material properties from the slave side
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,PoyntingUseMuR_Inv,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#endif
! Determine mu_r on faces within a dielectric region for calculating the Poynting vector and communicate the
! prolonged values via MPI
#if (PP_nVar>=6)
IF(PoyntingUseMuR_Inv) CALL SetDielectricFaceProfileForPoynting()
#endif /*(PP_nVar>=6)*/

ALLOCATE(sumFaces(nPoyntingIntPlanes))
#if USE_MPI
sumFaces=0
sumAllFaces=0
  CALL MPI_REDUCE(nFaces , sumFaces , nPoyntingIntPlanes , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  !nFaces(:) = sumFaces(:)
  CALL MPI_REDUCE(nPoyntingIntSides , sumAllFaces , 1 , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  !nPoyntingIntSides = sumAllFaces
#else
sumFaces=nFaces
sumAllFaces=nPoyntingIntSides
#endif /*USE_MPI*/

DO iPlane= 1, nPoyntingIntPlanes
  LBWRITE(UNIT_stdOut,'(A,I2,A,I10,A)') 'Processed plane no.: ',iPlane,'. Found ',sumFaces(iPlane),' surfaces.'
END DO
LBWRITE(UNIT_stdOut,'(A,I10,A)') 'A total of',sumAllFaces, &
                        ' surfaces for the poynting vector integral calculation are found.'

ALLOCATE(S    (1:3,0:PP_N,0:PP_N,1:nPoyntingIntSides) , &
         STEM     (0:PP_N,0:PP_N,1:nPoyntingIntSides)  )

LBWRITE(UNIT_stdOut,'(A)') ' ... POYNTING VECTOR INTEGRAL INITIALIZATION DONE.'

END SUBROUTINE GetPoyntingIntPlane


SUBROUTINE FinalizePoyntingInt()
!===================================================================================================================================
! Finalize Poynting Integral
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars ,ONLY:PosPoyntingInt,S,STEM,isPoyntingIntSide,SideIDToPoyntingSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! DEALLOCATE ALL
SDEALLOCATE(isPoyntingIntSide)
SDEALLOCATE(PosPoyntingInt)
SDEALLOCATE(SideIDToPoyntingSide)
SDEALLOCATE(S)
SDEALLOCATE(STEM)

END SUBROUTINE FinalizePoyntingInt
#endif /*PP_nVar>=6*/


#if (PP_nVar==8)
SUBROUTINE CalcPotentialEnergy(WEl, WMag, Wphi, Wpsi)
#else
SUBROUTINE CalcPotentialEnergy(WEl)
#endif /*PP_nVar=8*/
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,          ONLY : nElems, sJ
USE MOD_Interpolation_Vars, ONLY : wGP
#if (PP_nVar==8)
USE MOD_Globals_Vars       ,ONLY: smu0
#endif /*PP_nVar=8*/
USE MOD_Globals_Vars       ,ONLY: eps0
#if !(USE_HDG)
USE MOD_DG_Vars,            ONLY : U
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars      ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B,E
#endif /*PP_nVar==1*/
#else
USE MOD_PML_Vars           ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl
#if (PP_nVar==8)
REAL,INTENT(OUT)                :: WMag,Wpsi,Wphi
#endif /*PP_nVar=8*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL              :: WEl_tmp, WMag_tmp, E_abs
#if !(USE_HDG)
REAL              :: B_abs , Phi_abs, Psi_abs
#endif
#if USE_MPI
REAL              :: RD
#endif
#if (PP_nVar==8)
REAL              :: Wphi_tmp, Wpsi_tmp
#endif /*PP_nVar=8*/
!===================================================================================================================================

Wel=0.
#if (PP_nVar==8)
WMag=0.
Wphi=0.
Wpsi=0.
#endif /*PP_nVar=8*/
DO iElem=1,nElems
#if !(USE_HDG)
  IF(DoPML)THEN
    IF(isPMLElem(iElem))CYCLE
  END IF
#endif
  !--- Calculate and save volume of element iElem
  WEl_tmp=0.
  WMag_tmp=0.
#if (PP_nVar==8)
    Wphi_tmp = 0.
    Wpsi_tmp = 0.
#endif /*PP_nVar=8*/
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
! in electromagnetische felder by henke 2011 - springer
! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
    E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#elif PP_nVar==3
    B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
    E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
    B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
    E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
#endif /*USE_HDG*/

#if (PP_nVar==8)
    B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
    Phi_abs = U(7,i,j,k,iElem)*U(7,i,j,k,iElem)
    Psi_abs = U(8,i,j,k,iElem)*U(8,i,j,k,iElem)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
    WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs
#if (PP_nVar==8)
    WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
    Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Phi_abs
    Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Psi_abs
#endif /*PP_nVar=8*/
  END DO; END DO; END DO
  WEl = WEl + WEl_tmp
#if (PP_nVar==8)
  WMag = WMag + WMag_tmp
  Wphi = Wphi + Wphi_tmp
  Wpsi = Wpsi + Wpsi_tmp
#endif /*PP_nVar=8*/
END DO

WEl = WEl * eps0 * 0.5
#if (PP_nVar==8)
WMag = WMag * smu0 * 0.5
! caution: change of coefficients for divergence energies
Wphi = Wphi * eps0*0.5
Wpsi = Wpsi * smu0*0.5
#endif /*PP_nVar=8*/

#if USE_MPI
! todo: only one reduce with array
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,WEl  , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(MPI_IN_PLACE,WMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,Wphi , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,Wpsi , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
ELSE
  CALL MPI_REDUCE(WEl         ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(WMag        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(Wphi        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  CALL MPI_REDUCE(Wpsi        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
END IF
#endif /*USE_MPI*/

END SUBROUTINE CalcPotentialEnergy


#if (PP_nVar==8)
SUBROUTINE CalcPotentialEnergy_Dielectric(WEl, WMag, Wphi, Wpsi)
#else
SUBROUTINE CalcPotentialEnergy_Dielectric(WEl)
#endif /*PP_nVar=8*/
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
#if USE_HDG
#if PP_nVar==3 || PP_nVar==4
USE MOD_Dielectric_vars    ,ONLY: DielectricMu
#endif /*PP_nVar==3 or 4*/
#endif /*USE_HDG or PP_nVar==8*/
USE MOD_Mesh_Vars          ,ONLY: nElems, sJ
USE MOD_Interpolation_Vars ,ONLY: wGP
#if (PP_nVar==8)
USE MOD_Dielectric_vars    ,ONLY: DielectricMu
USE MOD_Globals_Vars       ,ONLY: smu0
#endif /*PP_nVar=8*/
USE MOD_Globals_Vars       ,ONLY: eps0
USE MOD_Dielectric_vars    ,ONLY: isDielectricElem,DielectricEps,ElemToDielectric
#if !(USE_HDG)
USE MOD_DG_Vars            ,ONLY: U
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars      ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars      ,ONLY: B
#else
USE MOD_Equation_Vars      ,ONLY: B,E
#endif /*PP_nVar==1*/
#else
USE MOD_PML_Vars           ,ONLY: DoPML,isPMLElem
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl
#if (PP_nVar==8)
REAL,INTENT(OUT)                :: WMag,Wpsi,Wphi
#endif /*PP_nVar=8*/
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL              :: WEl_tmp, WMag_tmp, E_abs
#if !(USE_HDG)
REAL              :: B_abs , Phi_abs, Psi_abs
#endif
#if USE_MPI
REAL              :: RD
#endif
#if (PP_nVar==8)
REAL              :: Wphi_tmp, Wpsi_tmp
#endif /*PP_nVar=8*/
!===================================================================================================================================

Wel=0.
#if (PP_nVar==8)
WMag=0.
Wphi=0.
Wpsi=0.
#endif /*PP_nVar=8*/

DO iElem=1,nElems
#if !(USE_HDG)
  IF(DoPML)THEN
    IF(isPMLElem(iElem))CYCLE
  END IF
#endif
  !--- Calculate and save volume of element iElem
  WEl_tmp=0.
  WMag_tmp=0.
  J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)





  IF(isDielectricElem(iElem))THEN
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! in electromagnetische felder by henke 2011 - springer
      ! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#elif PP_nVar==3
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
      E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
#endif /*USE_HDG*/

#if (PP_nVar==8)
      B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
      Phi_abs = U(7,i,j,k,iElem)*U(7,i,j,k,iElem)
      Psi_abs = U(8,i,j,k,iElem)*U(8,i,j,k,iElem)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs / DielectricMu( i,j,k,ElemToDielectric(iElem))
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs / DielectricMu( i,j,k,ElemToDielectric(iElem))
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
      WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs * DielectricEps(i,j,k,ElemToDielectric(iElem))
#if (PP_nVar==8)
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs / DielectricMu(i,j,k,ElemToDielectric(iElem))
      Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Phi_abs * DielectricEps(i,j,k,ElemToDielectric(iElem))
      Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Psi_abs / DielectricMu(i,j,k,ElemToDielectric(iElem))
#endif /*PP_nVar=8*/
    END DO; END DO; END DO
  ELSE
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! in electromagnetische felder by henke 2011 - springer
      ! WMag = 1/(2mu) * int_V B^2 dV
#if USE_HDG
#if PP_nVar==1
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#elif PP_nVar==3
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
      E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
#endif /*USE_HDG*/

#if (PP_nVar==8)
      B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
      Phi_abs = U(7,i,j,k,iElem)*U(7,i,j,k,iElem)
      Psi_abs = U(8,i,j,k,iElem)*U(8,i,j,k,iElem)
#endif /*PP_nVar=8*/
#if USE_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#endif /*PP_nVar==3*/
#endif /*USE_HDG*/
      WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs
#if (PP_nVar==8)
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
      Wphi_tmp = Wphi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Phi_abs
      Wpsi_tmp = Wpsi_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * Psi_abs
#endif /*PP_nVar=8*/
    END DO; END DO; END DO
  END IF


    WEl = WEl + WEl_tmp
#if (PP_nVar==8)
    WMag = WMag + WMag_tmp
#endif /*PP_nVar=8*/




END DO

WEl = WEl * eps0 * 0.5
#if (PP_nVar==8)
WMag = WMag * smu0 * 0.5
! caution: change of coefficients for divergence energies
Wphi = Wphi * eps0*0.5
Wpsi = Wpsi * smu0*0.5
#endif /*PP_nVar=8*/

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,WEl  , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(MPI_IN_PLACE,WMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
ELSE
  CALL MPI_REDUCE(WEl         ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#if (PP_nVar==8)
  CALL MPI_REDUCE(WMag        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
#endif /*PP_nVar=8*/
END IF
#endif /*USE_MPI*/

END SUBROUTINE CalcPotentialEnergy_Dielectric


#if (PP_nVar>=6)
SUBROUTINE SetDielectricFaceProfileForPoynting()
!===================================================================================================================================
!> THIS ROUTINE IS ONLY CALLED FOR THE POYNTING VECTOR INTEGRAL CALCULATION ON INITIALIZATION
!>
!> Set the dielectric factor 1./MuR for each face DOF in the array "Dielectric_MuR_Master_inv" (needed for S = E X H calculation).
!> Only the array "Dielectric_MuR_Master_inv" is used in the Riemann solver, as only the master calculates the flux array
!> (maybe slave information is used in the future)
!>
!> Note:
!> for MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!> threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!> is not allowed)
!> re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!> information)
!> This could be overcome by using template subroutines .t90 (see FlexiOS)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: Dielectric_MuR_Master_inv,Dielectric_MuR_Slave_inv
USE MOD_Dielectric_Vars ,ONLY: isDielectricElem,ElemToDielectric,DielectricMu
USE MOD_Mesh_Vars       ,ONLY: nSides
USE MOD_ProlongToFace   ,ONLY: ProlongToFace
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI             ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
USE MOD_FillMortar      ,ONLY: U_Mortar
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
#if USE_MPI
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Master2
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Slave2
INTEGER                                                  :: I,J,iSide
#endif /*USE_MPI*/
INTEGER                                                  :: iElem
!===================================================================================================================================
! General workflow:
! 1.  Initialize dummy arrays for Elem/Face
! 2.  Fill dummy element values for non-Dielectric sides
! 3.  Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar)
! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!     is not allowed)
!     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!     information)
! 5.  Send/Receive MPI data
! 6.  Allocate the actually needed arrays containing the dielectric material information on the sides
! 7.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)

! 1.  Initialize dummy arrays for Elem/Face
Dielectric_dummy_elem    = -1.
Dielectric_dummy_Master  = -1.
Dielectric_dummy_Slave   = -1.

! 2.  Fill dummy element values for non-Dielectric sides
DO iElem=1,PP_nElems
  IF(isDielectricElem(iElem))THEN
    ! set only the first dimension to 1./MuR (the rest are dummies)
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=1.0 / DielectricMu(0:PP_N,0:PP_N,0:PP_N,ElemToDielectric(iElem))
  ELSE
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=1.0
  END IF
END DO

!3.   Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar)
CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.FALSE.)
CALL U_Mortar(Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.FALSE.)
#if USE_MPI
  CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.TRUE.)
  CALL U_Mortar(Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.TRUE.)

  ! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
  !     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
  !     is not allowed)
  !     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
  !     information)
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

  ! 5.  Send Slave Dielectric info (real array with dimension (N+1)*(N+1)) to Master procs
  CALL StartReceiveMPIData(1,Dielectric_dummy_Slave2 ,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
  CALL StartSendMPIData(   1,Dielectric_dummy_Slave2 ,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

  ! Send Master Dielectric info (real array with dimension (N+1)*(N+1)) to Slave procs
  CALL StartReceiveMPIData(1,Dielectric_dummy_Master2,1,nSides ,RecRequest_U ,SendID=1) ! Receive YOUR
  CALL StartSendMPIData(   1,Dielectric_dummy_Master2,1,nSides,SendRequest_U ,SendID=1) ! Send MINE

  CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE - receive YOUR
  CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U ,SendID=1) !Send YOUR - receive MINE
#endif /*USE_MPI*/

! 6.  Allocate the actually needed arrays containing the dielectric material information on the sides
ALLOCATE(Dielectric_MuR_Master_inv(0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Dielectric_MuR_Slave_inv( 0:PP_N,0:PP_N,1:nSides))


! 7.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
#if USE_MPI
  Dielectric_MuR_Master_inv=Dielectric_dummy_Master2(1,0:PP_N,0:PP_N,1:nSides)
  Dielectric_MuR_Slave_inv =Dielectric_dummy_Slave2( 1,0:PP_N,0:PP_N,1:nSides)
#else
  Dielectric_MuR_Master_inv=Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides)
  Dielectric_MuR_Slave_inv =Dielectric_dummy_Slave( 1,0:PP_N,0:PP_N,1:nSides)
#endif /*USE_MPI*/

! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)
IF(MINVAL(Dielectric_MuR_Master_inv).LE.0.0)THEN
  CALL abort(&
  __STAMP__&
  ,'Dielectric material values for Riemann solver not correctly determined. MINVAL(Dielectric_MuR_Master_inv)=',&
  RealInfoOpt=MINVAL(Dielectric_MuR_Master_inv))
END IF
END SUBROUTINE SetDielectricFaceProfileForPoynting
#endif /*(PP_nVar>=6)*/


#if USE_HDG
SUBROUTINE CalculateAverageElectricPotential()
!===================================================================================================================================
! Calculation of the average electric potential with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is
! missing ... ups
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems, SurfElem
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_Analyze_Vars       ,ONLY: isAverageElecPotSide,AverageElectricPotential,AverageElectricPotentialFaces
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus,wGPSurf
USE MOD_DG_Vars            ,ONLY: U,U_master
#if USE_MPI
USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem, SideID,ilocSide
INTEGER          :: p,q,l
REAL             :: Uface(PP_nVar,0:PP_N,0:PP_N)
!REAL             :: SIP(0:PP_N,0:PP_N)
LOGICAL,PARAMETER:: Prolong=.TRUE.
REAL             :: AverageElectricPotentialProc
REAL             :: area_loc,integral_loc
!===================================================================================================================================

!IF(PRESENT(doProlong))THEN
  !Prolong=doProlong
!ELSE
  !Prolong=.TRUE.
!ENDIF

AverageElectricPotentialProc = 0.

DO iELEM = 1, nElems
  Do ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isAverageElecPotSide(SideID)) CYCLE
      IF(Prolong)THEN
#if (PP_NodeType==1) /* for Gauss-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,l,p,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,0,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,l,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,p,q,l,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! qfirst stuff
        CASE(XI_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,0,p,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,l,p,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,0,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+U(:,p,l,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,0,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,q,l,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        END SELECT
#else /* for Gauss-Lobatto-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          Uface(:,:,:)=U(:,:,0,:,iElem)
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
          Uface(:,:,:)=U(:,PP_N,:,:,iElem)
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,PP_N,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,PP_N,iElem)
            END DO ! p
          END DO ! q
        END SELECT
#endif
        ELSE ! no prolonge to face
          Uface=U_master(:,:,:,SideID)
        END IF ! Prolong

        ! multiplied by surface element and Gauss Points
        area_loc     = SUM(SurfElem(:,:,SideID) * wGPSurf(:,:))
        integral_loc = SUM(Uface(1,:,:) * SurfElem(:,:,SideID) * wGPSurf(:,:))
        AverageElectricPotentialProc = AverageElectricPotentialProc + integral_loc/area_loc
    END IF ! flip =0
  END DO ! iSides
END DO ! iElems
!AverageElectricPotentialProc = AverageElectricPotentialProc / (1e-4 * 1.28e-2)

#if USE_MPI
  CALL MPI_ALLREDUCE(AverageElectricPotentialProc , AverageElectricPotential , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
  AverageElectricPotential = AverageElectricPotentialProc
#endif /*USE_MPI*/

! Get global average value
AverageElectricPotential = AverageElectricPotential / AverageElectricPotentialFaces

END SUBROUTINE CalculateAverageElectricPotential


SUBROUTINE GetAverageElectricPotentialPlane()
!===================================================================================================================================
!> Initializes Poynting vector integral variables and check every side: set "isPoyntingIntSide(SideID) = .TRUE." if a side coincides
!> with a defined Poynting vector integral plane.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars    ,ONLY: nSides,nElems,Face_xGP
USE MOD_Mesh_Vars    ,ONLY: ElemToSide,normvec
USE MOD_Analyze_Vars ,ONLY: isAverageElecPotSide,AverageElectricPotentialCoordErr,PosAverageElectricPotential
USE MOD_ReadInTools  ,ONLY: GETINT,GETREAL
USE MOD_Globals      ,ONLY: abort
#if USE_MPI
USE MOD_Globals
#endif
USE MOD_Analyze_Vars ,ONLY: AverageElectricPotential,AverageElectricPotentialFaces
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iSide, SideID
INTEGER             :: nAverageElecPotSides
REAL                :: diff
INTEGER             :: p,q
!===================================================================================================================================

AverageElectricPotentialCoordErr = GETREAL('AvgPotential-Plane-Tolerance')
SWRITE(UNIT_stdOut,'(A)') ' GET PLANE TO CALCULATE AVERAGE ELECTRIC POTENTIAL ...'

! Initialize number of Poynting plane sides zero and set all sides to false
ALLOCATE(isAverageElecPotSide(1:nSides))
isAverageElecPotSide = .FALSE.

! Counter
AverageElectricPotential=0. ! initialize
nAverageElecPotSides = 0

! Loop over all elements
DO iElem=1,nElems
  ! Loop over all local sides
  DO iSide=1,6
    IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
      ASSOCIATE( AverageElecPotNormalDir1 => 2 ,&
                 AverageElecPotNormalDir2 => 3 ,&
                 AverageElecPotMainDir    => 1 )
        ! First search only planes with normal vector parallel to direction of "MainDir"
        IF((     NormVec(AverageElecPotNormalDir1,0,0,SideID)  < AverageElectricPotentialCoordErr) .AND. &
           (     NormVec(AverageElecPotNormalDir2,0,0,SideID)  < AverageElectricPotentialCoordErr) .AND. &
           ( ABS(NormVec(AverageElecPotMainDir   ,0,0,SideID)) > AverageElectricPotentialCoordErr))THEN
        ! Loop over all Points on Face
          DO q=0,PP_N
            DO p=0,PP_N
              diff = ABS(Face_xGP(AverageElecPotMainDir,p,q,SideID) - PosAverageElectricPotential)
              IF (diff < AverageElectricPotentialCoordErr) THEN
                IF (.NOT.isAverageElecPotSide(SideID)) THEN
                  nAverageElecPotSides = nAverageElecPotSides +1
                  isAverageElecPotSide(SideID) = .TRUE.
                END IF
              END IF ! diff < eps
            END DO !p
          END DO !q
        END IF
      END ASSOCIATE
    END IF ! flip = 0 master side
  END DO ! iSides
END DO !iElem=1,nElems

#if USE_MPI
CALL MPI_ALLREDUCE(nAverageElecPotSides , AverageElectricPotentialFaces , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
AverageElectricPotentialFaces=nAverageElecPotSides
#endif /*USE_MPI*/
SWRITE(UNIT_stdOut,'(A,I10,A)') ' A total of',AverageElectricPotentialFaces,&
    ' surfaces for the average electric potential calculation are found.'
SWRITE(UNIT_stdOut,'(A)') ' ... AVERAGE ELECTRIC POTENTIAL INITIALIZATION DONE.'
#if USE_MPI
IF(MPIRoot)THEN
#endif /*USE_MPI*/
  IF(AverageElectricPotentialFaces.EQ.0)THEN
    SWRITE(UNIT_stdOut,*) 'ERROR with: PosAverageElectricPotential = ',PosAverageElectricPotential
    CALL abort(__STAMP__&
    ,'Found zero faces for averaging the electric potential. Please make sure \nthat the x-coordinate coincides with element'//&
    ' interfaces. Planes cutting through elements in currently not implemented.')
  END IF ! AverageElectricPotentialFaces.EQ.0
#if USE_MPI
END IF ! MPIRoot
#endif /*USE_MPI*/

END SUBROUTINE GetAverageElectricPotentialPlane


SUBROUTINE FinalizeAverageElectricPotential()
!===================================================================================================================================
! Finalize Poynting Integral
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars ,ONLY:isAverageElecPotSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! DEALLOCATE ALL
SDEALLOCATE(isAverageElecPotSide)
END SUBROUTINE FinalizeAverageElectricPotential


SUBROUTINE CalculateElectricDisplacementCurrentSurface()
!===================================================================================================================================
!> Calculation of the average electric potential with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is
!> missing ... ups
!>
!> 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!      interpolate the vector field Et = (/Etx, Ety, Etz/) to the boundary face
!> 2.) Apply the normal vector: Uface(1,:,:)=DOT_PRODUCT(Uface(1:3,:,:),NormVec(1:3,:,:,SideID))
!      Store result of dot product in first array index
!> 3.) Get BC index and EDC index and the mapping of the SideID boundary to the EDC boundary ID and store the integrated current
!> 4.) Communicate the integrated current values on each boundary to the MPI root process (the root outputs the values to .csv)
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: SurfElem,SideToElem,nBCSides,NormVec,BC
USE MOD_Analyze_Vars       ,ONLY: EDC
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus,wGPSurf
USE MOD_Equation_Vars      ,ONLY: Et
#if USE_MPI
USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ElemID,SideID,ilocSide
INTEGER          :: p,q,l
REAL             :: Uface(1:3,0:PP_N,0:PP_N)
INTEGER          :: iBC,iEDCBC
!REAL             :: SIP(0:PP_N,0:PP_N)
!REAL             :: AverageElectricPotentialProc
!REAL             :: area_loc,integral_loc
!===================================================================================================================================
! Nullify
EDC%Current = 0.

! 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!     interpolate the vector field Et = (/Etx, Ety, Etz/) to the boundary face
DO SideID=1,nBCSides
  ElemID   = SideToElem(S2E_ELEM_ID,SideID)
  ilocSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
#if (PP_NodeType==1) /* for Gauss-points*/
  SELECT CASE(ilocSide)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Et(:,0,p,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+Et(:,l,p,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Et(:,p,0,q,ElemID)*L_Minus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Et(:,p,l,q,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Et(:,p,q,0,ElemID)*L_Minus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,q,p)=Uface(:,q,p)+Et(:,p,q,l,ElemID)*L_Minus(l)
        END DO ! l
      END DO ! p
    END DO ! qfirst stuff
  CASE(XI_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Et(:,0,p,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Et(:,l,p,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,PP_N-p,q)=Et(:,p,0,q,ElemID)*L_Plus(0)
        DO l=1,PP_N
          ! switch to right hand system
          Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+Et(:,p,l,q,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Et(:,p,q,0,ElemID)*L_Plus(0)
        DO l=1,PP_N
          Uface(:,p,q)=Uface(:,p,q)+Et(:,p,q,l,ElemID)*L_Plus(l)
        END DO ! l
      END DO ! p
    END DO ! q
  END SELECT
#else /* for Gauss-Lobatto-points*/
  SELECT CASE(ilocSide)
  CASE(XI_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Et(:,0,p,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ETA_MINUS)
    Uface(:,:,:)=Et(:,:,0,:,ElemID)
  CASE(ZETA_MINUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,q,p)=Et(:,p,q,0,ElemID)
      END DO ! p
    END DO ! q
  CASE(XI_PLUS)
    Uface(:,:,:)=Et(:,PP_N,:,:,ElemID)
  CASE(ETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,PP_N-p,q)=Et(:,p,PP_N,q,ElemID)
      END DO ! p
    END DO ! q
  CASE(ZETA_PLUS)
    DO q=0,PP_N
      DO p=0,PP_N
        Uface(:,p,q)=Et(:,p,q,PP_N,ElemID)
      END DO ! p
    END DO ! q
  END SELECT
#endif


  ! 2.) Apply the normal vector: Uface(1,:,:)=DOT_PRODUCT(Uface(1:3,:,:),NormVec(1:3,:,:,SideID))
  !     Store result of dot product in first array index
  Uface(1,:,:) =   Uface(1,:,:) * NormVec(1,:,:,SideID) &
                 + Uface(2,:,:) * NormVec(2,:,:,SideID) &
                 + Uface(3,:,:) * NormVec(3,:,:,SideID)

  ! 3.) Get BC index and EDC index and the mapping of the SideID boundary to the EDC boundary ID and store the integrated current
  iBC    = BC(SideID)
  iEDCBC = EDC%BCIDToEDCBCID(iBC)
  EDC%Current(iEDCBC) = EDC%Current(iEDCBC) + SUM(Uface(1,:,:) * SurfElem(:,:,SideID) * wGPSurf(:,:))

END DO

#if USE_MPI
! 4.) Communicate the integrated current values on each boundary to the MPI root process (the root outputs the values to .csv)
DO iEDCBC = 1, EDC%NBoundaries
  IF(EDC%COMM(iEDCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
    ASSOCIATE( Current => EDC%Current(iEDCBC), COMM => EDC%COMM(iEDCBC)%UNICATOR)
      IF(MPIroot)THEN
        CALL MPI_REDUCE(MPI_IN_PLACE , Current , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , COMM , IERROR)
      ELSE
        CALL MPI_REDUCE(Current      , 0       , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , COMM , IERROR)
      END IF ! myLeaderGroupRank.EQ.0
    END ASSOCIATE
  END IF ! EDC%COMM(iEDCBC)%UNICATOR.NE.MPI_COMM_NULL
END DO ! iEDCBC = 1, EDC%NBoundaries
#endif /*USE_MPI*/

END SUBROUTINE CalculateElectricDisplacementCurrentSurface
#endif /*USE_HDG*/

!===================================================================================================================================
!> Determine the field boundary output (BFO) values for the given iBC
!===================================================================================================================================
SUBROUTINE CalculateBoundaryFieldOutput(iBC,Time,BoundaryFieldOutput)
! MODULES
#if USE_HDG
USE MOD_Mesh_Vars ,ONLY: BoundaryType
USE MOD_Equation  ,ONLY: ExactFunc
#else
USE MOD_Globals   ,ONLY: abort
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: iBC
REAL,INTENT(IN)    :: Time
REAL,INTENT(OUT)   :: BoundaryFieldOutput(1:PP_nVar)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_HDG
INTEGER           :: BCType,BCState
#else
INTEGER           :: dummy
#endif /*USE_HDG*/
!===================================================================================================================================
BoundaryFieldOutput=0.!Initialize
#if USE_HDG
#if (PP_nVar==1)
BCType =BoundaryType(iBC,BC_TYPE)
BCState=BoundaryType(iBC,BC_STATE)
ASSOCIATE( x => (/0., 0., 0./) )
  SELECT CASE(BCType)
  CASE(2) ! exact BC = Dirichlet BC !! ExactFunc via BCState (time is optional)
    CALL ExactFunc(BCState , x , BoundaryFieldOutput , t=Time)
  CASE(4) ! exact BC = Dirichlet BC !! Zero potential
    BoundaryFieldOutput = 0.
  CASE(5,51,52,60) ! exact BC = Dirichlet BC !! ExactFunc via RefState (time is optional)
    CALL ExactFunc(  -1    , x , BoundaryFieldOutput , t=Time  , iRefState=BCState)
  CASE(6) ! exact BC = Dirichlet BC !! ExactFunc via RefState (Time is optional)
    CALL ExactFunc(  -2    , x , BoundaryFieldOutput , t=time  , iRefState=BCState)
  CASE(50) ! exact BC = Dirichlet BC !! ExactFunc via bias voltage DC
    CALL ExactFunc(  -5    , x , BoundaryFieldOutput )
  END SELECT ! BCType
END ASSOCIATE
#else
CALL abort(__STAMP__,'CalculateBoundaryFieldOutput is not implemented for PP_nVar>1')
#endif /*PP_nVar==1*/
#else
CALL abort(__STAMP__,'CalculateBoundaryFieldOutput is not implemented for other equation systems yet (only HDG)')
! Suppress warnings
dummy=iBC
dummy=INT(Time)
#endif /*USE_HDG*/

END SUBROUTINE CalculateBoundaryFieldOutput

END MODULE MOD_AnalyzeField
