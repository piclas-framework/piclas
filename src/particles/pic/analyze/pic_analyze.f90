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

MODULE MOD_PIC_Analyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CalcDepositedCharge
  MODULE PROCEDURE CalcDepositedCharge
END INTERFACE

#if USE_HDG
INTERFACE CalculateBRElectronsPerCell
  MODULE PROCEDURE CalculateBRElectronsPerCell
END INTERFACE

PUBLIC:: CalculateBRElectronsPerCell
#endif /*USE_HDG*/

PUBLIC:: VerifyDepositedCharge, CalcDepositedCharge
!===================================================================================================================================

CONTAINS

SUBROUTINE VerifyDepositedCharge()
!===================================================================================================================================
! Calculate the deposited charges
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: nElems, N_VolMesh
USE MOD_Interpolation_Vars    ,ONLY: NAnalyze,N_InterAnalyze,wAnalyze
USE MOD_Particle_Vars         ,ONLY: PDM, Species, PartSpecies ,PartMPF,usevMPF
USE MOD_Particle_Analyze_Vars ,ONLY: ChargeCalcDone
USE MOD_PICDepo_Vars          ,ONLY: sfDepo3D,dimFactorSF,VerifyChargeStr,DepositionType
#if defined(IMPA)
USE MOD_LinearSolver_Vars     ,ONLY: ImplicitSource
#else
USE MOD_PICDepo_Vars          ,ONLY: PS_N
#endif
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_DG_Vars               ,ONLY: N_DG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: k,l,m,i,iElem,Nloc
REAL              :: ChargeNumerical, ChargeLoc, ChargeAnalytical
REAL              :: Coords_NAnalyze(3,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL              :: U_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL              :: J_NAnalyze(1,0:NAnalyze,0:NAnalyze,0:NAnalyze)
REAL              :: IntegrationWeight
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' PERFORMING CHARGE DEPOSITION PLAUSIBILITY CHECK...'

ChargeNumerical=0. ! Nullify
DO iElem=1,nElems
  Nloc = N_DG(iElem)
  ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
  CALL ChangeBasis3D(3,Nloc,NAnalyze,N_InterAnalyze(Nloc)%Vdm_GaussN_NAnalyze,N_VolMesh(iElem)%Elem_xGP(1:3,:,:,:),Coords_NAnalyze(1:3,:,:,:))
  ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the Jacobian ;-)
  CALL ChangeBasis3D(1,Nloc,NAnalyze,N_InterAnalyze(Nloc)%Vdm_GaussN_NAnalyze,1./N_VolMesh(iElem)%sJ(:,:,:),J_NAnalyze(1:1,:,:,:))
  ! Interpolate the solution to the analyze grid
#if defined(IMPA)
  CALL ChangeBasis3D(1,PP_N,NAnalyze,N_InterAnalyze(PP_N)%Vdm_GaussN_NAnalyze,ImplicitSource(4,:,:,:,iElem),U_NAnalyze(1,:,:,:))
#else
  CALL ChangeBasis3D(1,Nloc,NAnalyze,N_InterAnalyze(Nloc)%Vdm_GaussN_NAnalyze,PS_N(iElem)%PartSource(4,0:Nloc,0:Nloc,0:Nloc),U_NAnalyze(1,:,:,:))
#endif
  ChargeLoc=0. ! Nullify
  DO m=0,NAnalyze
    DO l=0,NAnalyze
      DO k=0,NAnalyze
        IntegrationWeight = wAnalyze(k)*wAnalyze(l)*wAnalyze(m)*J_NAnalyze(1,k,l,m)
        ChargeLoc = ChargeLoc + IntegrationWeight*U_NAnalyze(1,k,l,m)
      END DO ! k
    END DO ! l
  END DO ! m
  ChargeNumerical = ChargeNumerical + ChargeLoc
END DO

ChargeAnalytical=0.
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    IF(usevMPF)THEN
      ChargeAnalytical = ChargeAnalytical + Species(PartSpecies(i))%ChargeIC * PartMPF(i)
    ELSE
      ChargeAnalytical = ChargeAnalytical + Species(PartSpecies(i))%ChargeIC * Species(PartSpecies(i))%MacroParticleFactor
    END IF
  END IF
END DO

! Collect info on MPI root process
#if USE_MPI
   IF(MPIRoot) THEN
     CALL MPI_REDUCE(MPI_IN_PLACE     , ChargeAnalytical , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
     CALL MPI_REDUCE(MPI_IN_PLACE     , ChargeNumerical  , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
   ELSE
     CALL MPI_REDUCE(ChargeAnalytical , 0                , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
     CALL MPI_REDUCE(ChargeNumerical  , 0                , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , IERROR)
   END IF
#endif

! Output info to std.out
IF(MPIRoot)THEN
  ! Check if the charge is to be distributed over a line (1D) or area (2D)
  IF(StringBeginsWith(DepositionType,'shape_function').AND.(.NOT.sfDepo3D)) ChargeAnalytical = ChargeAnalytical * dimFactorSF
  WRITE(*,'(A,ES25.14E3,A)') 'On the grid deposited charge (numerical) : ', ChargeNumerical ,' [C] '//TRIM(VerifyChargeStr)
  WRITE(*,'(A,ES25.14E3,A)') 'Charge over by the particles (analytical): ', ChargeAnalytical,' [C] '//TRIM(VerifyChargeStr)
  WRITE(*,'(A,ES25.14E3,A)') 'Absolute deposition error                : ', ABS(ChargeAnalytical-ChargeNumerical),' [C]'
  IF(ABS(ChargeAnalytical).GT.0.0)THEN
    WRITE(*,'(A,ES25.14E3,A)') 'Relative deposition error in percent     : ', &
        ABS(ChargeAnalytical-ChargeNumerical)/ChargeAnalytical*100,' [%]'
  ELSE
    WRITE(*,'(A)') 'Relative deposition error in percent     :                       100 [%]'
  END IF
  WRITE(UNIT_stdOut,'(A)')' CHARGE DEPOSITION PLAUSIBILITY CHECK DONE!'
  WRITE(UNIT_StdOut,'(132("-"))')
END IF ! MPIRoot
ChargeCalcDone = .TRUE.

END SUBROUTINE VerifyDepositedCharge


SUBROUTINE CalcDepositedCharge()
!===================================================================================================================================
! Calculation of deposited charge and compute the absolute and relative error
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: N_VolMesh
USE MOD_Particle_Vars         ,ONLY: PDM, Species, PartSpecies, usevmpf, PartMPF
USE MOD_Interpolation_Vars    ,ONLY: N_Inter
USE MOD_Particle_Analyze_Vars ,ONLY: PartCharge
USE MOD_TimeDisc_Vars         ,ONLY: iter
#if defined(IMPA)
USE MOD_LinearSolver_Vars     ,ONLY: ImplicitSource
#else
USE MOD_PICDepo_Vars          ,ONLY: PS_N
#endif
USE MOD_DG_Vars               ,ONLY: N_DG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,iPart,iElem,Nloc
REAL              :: Charge(2)
#if USE_MPI
REAL              :: RECBR(2)
#endif /*USE_MPI*/
!===================================================================================================================================


! compute local charge
Charge=0.
PartCharge=0.
IF(iter.EQ.0) RETURN
DO iElem=1,PP_nElems
  Nloc = N_DG(iElem)
  ! compute the deposited charge
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    ASSOCIATE( wGPijk => N_Inter(Nloc)%wGP(i)*N_Inter(Nloc)%wGP(j)*N_Inter(Nloc)%wGP(k) ,&
               J_Nloc => 1./N_VolMesh(iElem)%sJ(i,j,k) )
#if defined(IMPA)
#if USE_HDG
      Charge(1) = Charge(1)+ wGPijk * ImplicitSource(1,i,j,k,iElem) * J_Nloc
#else /* DG */
      Charge(1) = Charge(1)+ wGPijk * ImplicitSource(4,i,j,k,iElem) * J_Nloc
#endif
#else
      Charge(1) = Charge(1)+ wGPijk * PS_N(iElem)%PartSource(4,i,j,k) * J_Nloc
#endif
    END ASSOCIATE
  END DO; END DO; END DO
END DO

! charge of all particles inside of domain
DO iPart=1,PDM%ParticleVecLength
  IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  IF(usevMPF)THEN
    Charge(2) = Charge(2) + Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
  ELSE
    Charge(2) = Charge(2) + Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF
END DO

! MPI Communication
#if USE_MPI
IF (MPIRoot) THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Charge , 2 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
ELSE ! no Root
  CALL MPI_REDUCE(Charge,RECBR  ,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_PICLAS, IERROR)
END IF
#endif

IF (MPIRoot) THEN
  PartCharge(1)=Charge(1)
  ! absolute error
  PartCharge(2)=ABS(Charge(2)-Charge(1))
  ! relative error
  IF(ALMOSTZERO(Charge(2)))THEN
    PartCharge(3)=0.
  ELSE
    PartCharge(3)=ABS(Charge(2)-Charge(1))/Charge(2)
  END IF
END IF

END SUBROUTINE CalcDepositedCharge


#if USE_HDG
SUBROUTINE CalculateBRElectronsPerCell(iElem,RegionID,ElectronNumberCell)
!===================================================================================================================================
! calcs integrated (physical) number of BR electrons in cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ElementaryCharge
USE MOD_Preproc
USE MOD_Mesh_Vars          ,ONLY: N_VolMesh
USE MOD_Interpolation_Vars ,ONLY: N_Inter
#if PP_nVar==1
USE MOD_DG_Vars            ,ONLY: U_N
#endif
USE MOD_HDG_Vars           ,ONLY: RegionElectronRef
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: iElem, RegionID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: ElectronNumberCell
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,k,Nloc
REAL              :: source_e
!===================================================================================================================================
ElectronNumberCell=0.
Nloc = N_DG(iElem)
DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
#if PP_nVar==1
  source_e = U_N(iElem)%U(1,i,j,k)-RegionElectronRef(2,RegionID)
#else
  CALL abort(__STAMP__,' CalculateBRElectronsPerCell only implemented for electrostatic HDG!')
#endif
  IF (source_e .LT. 0.) THEN
    source_e = RegionElectronRef(1,RegionID) &         !--- boltzmann relation (electrons as isothermal fluid!)
    * EXP( (source_e) / RegionElectronRef(3,RegionID) )
  ELSE
    source_e = RegionElectronRef(1,RegionID) &         !--- linearized boltzmann relation at positive exponent
    * (1. + ((source_e) / RegionElectronRef(3,RegionID)) )
  END IF
  ElectronNumberCell = ElectronNumberCell + N_Inter(Nloc)%wGP(i)*N_Inter(Nloc)%wGP(j)*N_Inter(Nloc)%wGP(k) * source_e * &
      (1./N_VolMesh(iElem)%sJ(i,j,k))
END DO; END DO; END DO
ElectronNumberCell=ElectronNumberCell/ElementaryCharge

END SUBROUTINE CalculateBRElectronsPerCell
#endif /*USE_HDG*/


END MODULE MOD_PIC_Analyze
