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

MODULE MOD_Analyze_DVM

INTERFACE AnalyzeDVM
  MODULE PROCEDURE AnalyzeDVM
END INTERFACE

PUBLIC :: AnalyzeDVM

CONTAINS

SUBROUTINE AnalyzeDVM()
!===================================================================================================================================
! Macro values for ElemData output
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars
USE MOD_TimeDisc_Vars         ,ONLY: dt, dt_Min
USE MOD_DistFunc              ,ONLY: MacroValuesFromDistribution
USE MOD_Mesh_Vars             ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tau=1., MacroVal(8)=0.
INTEGER                         :: iElem
!===================================================================================================================================
print*, 'analyzedvm'
DO iElem=1,nElems
    !CALL MacroValuesFromDistribution(MacroVal,U_FV(:,0,0,0,iElem),dt,tau,1)
    DVM_ElemData1(iElem)=MacroVal(1)
    DVM_ElemData2(iElem)=MacroVal(2)
    DVM_ElemData3(iElem)=MacroVal(3)
    DVM_ElemData4(iElem)=MacroVal(4)
    DVM_ElemData5(iElem)=MacroVal(5)
    DVM_ElemData6(iElem)=MacroVal(6)
    DVM_ElemData7(iElem)=MacroVal(7)
    DVM_ElemData8(iElem)=MacroVal(8)
    DVM_ElemData9(iElem)=dt_Min(DT_MIN)/tau !relaxation factor for regular timesteps, higher than for the actual analyze timestep
END DO

END SUBROUTINE

END MODULE MOD_Analyze_DVM