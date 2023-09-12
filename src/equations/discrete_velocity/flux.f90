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

MODULE MOD_Flux
!===================================================================================================================================
! Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if !(USE_FV)
INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D
END INTERFACE

INTERFACE EvalFlux3DDielectric
  MODULE PROCEDURE EvalFlux3DDielectric
END INTERFACE

PUBLIC::EvalFlux3D,EvalFlux3DDielectric
#endif /*!(USE_FV)*/
!===================================================================================================================================

CONTAINS

#if !(USE_FV)
SUBROUTINE EvalFlux3D(iElem,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_TimeDisc_Vars,ONLY : dt
USE MOD_FV_Vars      ,ONLY: U
USE MOD_Equation_Vars,ONLY: DVMnVelos, DVMVelos, DVMBGKModel, DVMMethod, DVMDim
USE MOD_DistFunc     ,ONLY: MacroValuesFromDistribution
USE MOD_DistFunc     ,ONLY: MaxwellDistribution, MaxwellDistributionCons, ShakhovDistribution, ESBGKDistribution
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: MacroVal(14), tau, fTarget(PP_nVar), UTemp(PP_nVar), gamma
INTEGER             :: i,j,k, kVel, jVel, iVel, upos
!==================================================================================================================================
DO k=0,PP_N
  DO j=0,PP_N
    DO i=0,PP_N
      UTemp=U(:,i,j,k,iElem)
      CALL MacroValuesFromDistribution(MacroVal,UTemp,dt/2.,tau,1)
      ! wrong because no reconstruction/relaxation was done before ---> DVM only works with FV
      SELECT CASE (DVMBGKModel)
        CASE(1)
          CALL ESBGKDistribution(MacroVal,fTarget)
        CASE(2)
          CALL ShakhovDistribution(MacroVal,fTarget)
        CASE(3)
          CALL MaxwellDistribution(MacroVal,fTarget)
        CASE(4)
          CALL MaxwellDistributionCons(MacroVal,fTarget)
        CASE DEFAULT
          CALL abort(__STAMP__,'DVM BGK Model not implemented')
      END SELECT
      IF (dt.GT.0.) THEN
        SELECT CASE (DVMMethod)
        CASE(1)
          gamma = 2.*tau*(1.-EXP(-dt/2./tau))/dt
        CASE(2)
          gamma = 2.*tau/(2.*tau+dt)
        END SELECT
        UTemp = gamma*UTemp + (1.-gamma)*fTarget
      ELSE
        UTemp = 0.
      END IF
      DO kVel=1, DVMnVelos(3);   DO jVel=1, DVMnVelos(2);   DO iVel=1, DVMnVelos(1)
        upos= iVel+(jVel-1)*DVMnVelos(1)+(kVel-1)*DVMnVelos(1)*DVMnVelos(2)
        f(upos,i,j,k) = DVMVelos(iVel,1) * UTemp(upos)
        g(upos,i,j,k) = DVMVelos(jVel,2) * UTemp(upos)
        h(upos,i,j,k) = DVMVelos(kVel,3) * UTemp(upos)
        IF (DVMDim.LT.3) THEN
          f(PP_nVar/2+upos,i,j,k) = DVMVelos(iVel,1) * UTemp(PP_nVar/2+upos)
          g(PP_nVar/2+upos,i,j,k) = DVMVelos(jVel,2) * UTemp(PP_nVar/2+upos)
          h(PP_nVar/2+upos,i,j,k) = DVMVelos(kVel,3) * UTemp(PP_nVar/2+upos)
        END IF
      END DO; END DO; END DO;
    END DO ! i
  END DO ! j
END DO ! k

END SUBROUTINE EvalFlux3D

SUBROUTINE EvalFlux3DDielectric(iElem,f,g,h)
!===================================================================================================================================
! dummy routine for call in dg/volint.f90
!===================================================================================================================================
! MODULES
USE MOD_PreProc

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(8,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h    ! Cartesian fluxes (iVar,i,j,k)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

END SUBROUTINE EvalFlux3DDielectric
#endif /*!(USE_FV)*/

END MODULE MOD_Flux