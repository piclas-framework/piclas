!==================================================================================================================================
! Copyright (c) 2015 - 2019 Wladimir Reschke
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_MacroBody
!===================================================================================================================================
!> Main Module for macroscopic bodies inside particle domain
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

PUBLIC :: MacroBody_main
!PUBLIC :: UpdateMacroBodyVars
!===================================================================================================================================

CONTAINS

SUBROUTINE MacroBody_main()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_TimeDisc_Vars ,ONLY: dt
USE MOD_Particle_Vars ,ONLY: MacroPartFluxesEnabled, MacroPart, UseMacroPart, nMacroParticle,MacroPartAccelerationEnabled
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                    :: iMP
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.UseMacroPart) RETURN

MacroPart(:)%center(1) = MacroPart(:)%center(1) + MacroPart(:)%velocity(1)*dt
MacroPart(:)%center(2) = MacroPart(:)%center(2) + MacroPart(:)%velocity(2)*dt
MacroPart(:)%center(3) = MacroPart(:)%center(3) + MacroPart(:)%velocity(3)*dt
IF(MacroPartAccelerationEnabled) THEN
  DO iMP=1,nMacroParticle
    MacroPart(iMP)%velocity(1:6) = MacroPart(iMP)%velocity(1:6) + MacroPart(iMP)%RHS(1:6)
    MacroPart(iMP)%RHS(1:6)=0.
  END DO
END IF
IF(MacroPartFluxesEnabled) THEN
  DO iMP=1,nMacroParticle
    MacroPart(iMP)%radius = MacroPart(iMP)%radius + MacroPart(iMP)%RHS(7)
    MacroPart(iMP)%temp   = MacroPart(iMP)%temp   + MacroPart(iMP)%RHS(8)
    MacroPart(iMP)%mass   = MacroPart(iMP)%mass   + MacroPart(iMP)%RHS(9)
    MacroPart(iMP)%RHS(7:9)=0.
  END DO
END IF

END SUBROUTINE MacroBody_main


!SUBROUTINE UpdateMacroBodyVars()
!!===================================================================================================================================
!!>
!!===================================================================================================================================
!! MODULES                                                                                                                          !
!USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
!USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues, KeepWallParticles, Species, nSpecies
!USE MOD_DSMC_Vars              ,ONLY: DSMC
!USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall, PartBound, SurfCOMM
!USE MOD_TimeDisc_Vars          ,ONLY: tend,time
!#if USE_LOADBALANCE
!USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
!#endif /*USE_LOADBALANCE*/
!USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo, SurfModel
!USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbProb, CalcDesorbProb
!USE MOD_SurfaceModel_Tools     ,ONLY: SMCR_AdjustMapNum, IsReactiveSurface, SurfaceHasModelNum
!#if USE_MPI
!USE MOD_SurfaceModel_MPI       ,ONLY: ExchangeSurfaceHaloToOrigin, ExchangeSurfaceOriginToHalo, ExchangeSurfDistInfo
!USE MOD_SurfaceModel_MPI       ,ONLY: MapHaloInnerToOriginInnerSurf
!#endif /*USE_MPI*/
!!----------------------------------------------------------------------------------------------------------------------------------!
!IMPLICIT NONE
!! INPUT / OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------!
!! LOCAL VARIABLES
!INTEGER                          :: iSpec, iSurfSide, p, q, new_adsorbates, numSites
!REAL                             :: maxPart
!REAL                             :: coverage_tmp, coverage_corrected
!#if USE_LOADBALANCE
!REAL                             :: tLBStart
!#endif /*USE_LOADBALANCE*/
!#if (PP_TimeDiscMethod==42)
!REAL                             :: desorbnum_covreduce
!#endif
!!----------------------------------------------------------------------------------------------------------------------------------!
!IF (.NOT.UseMacroPart) RETURN

!END SUBROUTINE UpdateMacroBodyVars


END MODULE MOD_MacroBody
