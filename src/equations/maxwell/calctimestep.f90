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

MODULE MOD_CalcTimeStep
!===================================================================================================================================
! Low-Storage Runge-Kutta integration of degree 3 for one step.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!===================================================================================================================================

CONTAINS

FUNCTION CALCTIMESTEP()
!===================================================================================================================================
! Calculate the time step for the current update of U for the Euler-Equations
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars     ,ONLY: N_VolMesh
USE MOD_Equation_Vars ,ONLY: c_corr
USE MOD_Globals_Vars  ,ONLY: c
USE MOD_TimeDisc_Vars ,ONLY: CFLScale
USE MOD_DG_Vars       ,ONLY: N_DG
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                         :: CalcTimeStep
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem,Nloc
REAL                         :: Max_Lambda1,Max_Lambda2,Max_Lambda3
REAL                         :: TimeStepConv,locTimeStepConv
!===================================================================================================================================
locTimeStepConv=HUGE(1.)
DO iElem=1,PP_nElems
  Nloc = N_DG(iElem)
  Max_Lambda1=0.
  Max_Lambda2=0.
  Max_Lambda3=0.
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        Max_Lambda1=MAX(Max_Lambda1,N_VolMesh(iElem)%sJ(i,j,k)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k)*N_VolMesh(iElem)%Metrics_fTilde(:,i,j,k)))))
        Max_Lambda2=MAX(Max_Lambda2,N_VolMesh(iElem)%sJ(i,j,k)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k)*N_VolMesh(iElem)%Metrics_gTilde(:,i,j,k)))))
        Max_Lambda3=MAX(Max_Lambda3,N_VolMesh(iElem)%sJ(i,j,k)*(MAX(1.,c_corr)*c &
                        *SQRT(SUM(N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k)*N_VolMesh(iElem)%Metrics_hTilde(:,i,j,k)))))
      END DO ! i
    END DO ! j
  END DO ! k

! VERSION 3: ---------------------------------
!   locTimeStepConv=MIN(locTimeStepConv,CFLScale*2./(Max_Lambda1))
! --------------------------------------------
! VERSION 2: quadratic superposition
!  locTimeStepConv=MIN(locTimeStepConv,CFLScale*2./SQRT(Max_Lambda1**2+Max_Lambda2**2+Max_Lambda3**2))
! --------------------------------------------
! VERSION 1: linear superposition
  locTimeStepConv=MIN(locTimeStepConv,CFLScale*2./(Max_Lambda1+Max_Lambda2+Max_Lambda3))
! --------------------------------------------
  IF(locTimeStepConv.NE.locTimeStepConv)THEN
    ERRWRITE(*,'(A,I0,A,I0,A,I0,A,I0)')'Convective timestep NaN on proc ',myRank,' at global position (iElem): ',iElem
    ERRWRITE(*,*)'dt_conv=',locTimeStepConv
    CALL abort(__STAMP__,'Convective timestep NaN!',999,999.)
  END IF
END DO ! iElem
#if USE_MPI
CALL MPI_ALLREDUCE(locTimeStepConv,TimeStepConv,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#else
TimeStepConv=locTimeStepConv
#endif /*USE_MPI*/
CalcTimeStep=TimeStepConv

END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
