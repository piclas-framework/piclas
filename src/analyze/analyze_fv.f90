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

MODULE MOD_Analyze_FV
!===================================================================================================================================
! Contains DG analyze
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

!===================================================================================================================================
#if USE_FV
PUBLIC:: CalcError_FV
#endif /*USE_FV*/
!===================================================================================================================================

CONTAINS

#if USE_FV
SUBROUTINE CalcError_FV(time,L_2_Error,L_Inf_Error)
!===================================================================================================================================
! Calculates L_infinfity and L_2 norms of finite volume state variables (cell-integrated values)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FV_Vars            ,ONLY: U_FV
USE MOD_Equation_FV        ,ONLY: ExactFunc_FV
USE MOD_Equation_Vars_FV   ,ONLY: IniExactFunc_FV
USE MOD_Mesh_Vars_FV       ,ONLY: Elem_xGP_FV
USE MOD_Particle_Mesh_Vars ,ONLY: MeshVolume
USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared
USE MOD_Analyze_Vars       ,ONLY: OutputErrorNormsToH5
#ifdef discrete_velocity
USE MOD_TimeDisc_Vars      ,ONLY: dt
USE MOD_DistFunc,           ONLY: MacroValuesFromDistribution
#endif /*discrete_velocity*/
#ifdef PARTICLES
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_Particle_Mesh_Vars ,ONLY: offsetComputeNodeElem
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)               :: time
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#ifdef discrete_velocity
REAL,INTENT(OUT)              :: L_2_Error(14)   !< L2 error of the solution
REAL,INTENT(OUT)              :: L_Inf_Error(14) !< LInf error of the solution
#else
REAL,INTENT(OUT)              :: L_2_Error(PP_nVar_FV)   !< L2 error of the solution
REAL,INTENT(OUT)              :: L_Inf_Error(PP_nVar_FV) !< LInf error of the solution
#endif /*discrete_velocity*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
REAL                          :: U_exact(1:PP_nVar_FV)
INTEGER                       :: offsetElemCNProc,CNElemID
#ifdef discrete_velocity
REAL                          :: MacroVal(14), MacroVal_exact(14), tau, real_dt
#endif /*discrete_velocity*/
!===================================================================================================================================
IF (OutputErrorNormsToH5) CALL abort(__STAMP__,'OutputErrorNormsToH5 not implemented for FV')
L_Inf_Error(:)=-1.E10
L_2_Error(:)=0.

#ifdef discrete_velocity
IF (time.EQ.0.) THEN
  real_dt = 0.
ELSE
  real_dt = dt
END IF
#endif /*discrete_velocity*/

DO iElem=1,PP_nElems
#if USE_MPI && defined(PARTICLES)
! J_N is only built for local DG elements. Therefore, array is only filled for elements on the same compute node
  offsetElemCNProc = offsetElem - offsetComputeNodeElem
#else
  offsetElemCNProc = 0
#endif  /*USE_MPI && defined(PARTICLES)*/
  CNElemID=iElem+offsetElemCNProc
  CALL ExactFunc_FV(IniExactFunc_FV,time,Elem_xGP_FV(1:3,0,0,0,iElem),U_exact(1:PP_nVar_FV))
#ifdef discrete_velocity
  ! DVM: calculate errors for the macroscopic values
  CALL MacroValuesFromDistribution(MacroVal,U_FV(:,0,0,0,iElem),real_dt,tau,1)
  CALL MacroValuesFromDistribution(MacroVal_exact,U_exact(:),real_dt,tau,1)
  L_Inf_Error = MAX(L_Inf_Error,abs(MacroVal(1:14) - MacroVal_exact(1:14)))
  ! To sum over the elements, We compute here the square of the L_2 error
  L_2_Error = L_2_Error+(MacroVal(1:14) - MacroVal_exact(1:14))*&
                        (MacroVal(1:14) - MacroVal_exact(1:14))*ElemVolume_Shared(CNElemID)
#else
  L_Inf_Error = MAX(L_Inf_Error,abs(U_FV(:,0,0,0,iElem) - U_exact(1:PP_nVar_FV)))
  ! To sum over the elements, We compute here the square of the L_2 error
  L_2_Error = L_2_Error+(U_FV(:,0,0,0,iElem) - U_exact(1:PP_nVar_FV))*&
                        (U_FV(:,0,0,0,iElem) - U_exact(1:PP_nVar_FV))*ElemVolume_Shared(CNElemID)
#endif /*discrete_velocity*/
END DO ! iElem=1,PP_nElems
#if USE_MPI
#ifdef discrete_velocity
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error   , 14 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error , 14 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error   , 0            , 14 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(L_Inf_Error , 0            , 14 , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#else
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE , L_2_Error   , PP_nVar_FV , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(MPI_IN_PLACE , L_Inf_Error , PP_nVar_FV , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
  ELSE
    CALL MPI_REDUCE(L_2_Error   , 0            , PP_nVar_FV , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , MPI_COMM_PICLAS , iError)
    CALL MPI_REDUCE(L_Inf_Error , 0            , PP_nVar_FV , MPI_DOUBLE_PRECISION , MPI_MAX , 0 , MPI_COMM_PICLAS , iError)
    ! in this case the receive value is not relevant.
  END IF
#endif /*discrete_velocity*/
#endif /*USE_MPI*/

! We normalize the L_2 Error with the Volume of the domain and take into account that we have to use the square root
L_2_Error = SQRT(L_2_Error/MeshVolume)

END SUBROUTINE CalcError_FV
#endif /*USE_FV*/


END MODULE MOD_Analyze_FV
