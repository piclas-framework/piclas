!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_FPFlow_Init
!===================================================================================================================================
! Initialization of DSMC
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
PUBLIC :: InitFPFlow, DefineParametersFPFlow, FinalizeFPFlow
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FP-Flow
!==================================================================================================================================
SUBROUTINE DefineParametersFPFlow()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("FP-Flow")

CALL prms%CreateIntOption(    'Particles-FP-CollModel',       'Select the Fokker-Planck method:\n'//&
                                                              '1: Cubic (only atomic species)\n'//&
                                                              '2: Ellipsoidal Statistical (ESFP)')
CALL prms%CreateIntOption(    'Particles-ESFP-Model',         'Select the ESFP model:\n'//&
                                                              '1: Exact\n'//&
                                                              '2: Approximative', '1')
CALL prms%CreateLogicalOption('Particles-FP-DoVibRelaxation', 'Enable modelling of vibrational excitation','.TRUE.')
CALL prms%CreateLogicalOption('Particles-FP-UseQuantVibEn',   'Enable quantized modelling of vibrational energy levels','.TRUE.')
CALL prms%CreateLogicalOption('Particles-FP-DoCellAdaptation','Enables octree cell refinement until the given number of '//&
                                                              'particles is reached. Equal refinement in all three '//&
                                                              'directions (x,y,z)','.FALSE.')
CALL prms%CreateIntOption(    'Particles-FP-MinPartsPerCell', 'Define minimum number of particles per cell for octree '//&
                                                              'cell refinement')
CALL prms%CreateLogicalOption('Particles-CoupledFPDSMC',      'Perform a coupled DSMC-FP simulation with a given number density'//&
                                                              'as a switch parameter','.FALSE.')
CALL prms%CreateRealOption(   'Particles-FP-DSMC-SwitchDens', 'Number density [1/m3] above which the FP method is used, below'//&
                                                              'which DSMC is performed.','0.0')

END SUBROUTINE DefineParametersFPFlow

SUBROUTINE InitFPFlow()
!===================================================================================================================================
! Initialization of the Fokker-Planck method
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars             ,ONLY: NGeo, nElems
USE MOD_Globals_Vars          ,ONLY: PI, BoltzmannConst
USE MOD_ReadInTools
USE MOD_DSMC_Vars             ,ONLY: DSMC, CollInf
USE MOD_DSMC_ParticlePairing  ,ONLY: DSMC_init_octree
USE MOD_PARTICLE_Vars         ,ONLY: nSpecies, Species, DoVirtualCellMerge
USE MOD_FPFlow_Vars
USE MOD_BGK_Vars              ,ONLY: DoBGKCellAdaptation, BGKMinPartPerCell
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER               :: iSpec, iSpec2
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') ' INIT FP Solver...'

ALLOCATE(SpecFP(nSpecies))

DO iSpec = 1, nSpecies
  ALLOCATE(SpecFP(iSpec)%CollFreqPreFactor(nSpecies))
  DO iSpec2=1, nSpecies
    SpecFP(iSpec)%CollFreqPreFactor(iSpec2)= 0.5*(CollInf%dref(iSpec,iSpec) + CollInf%dref(iSpec2,iSpec2))**2.0 &
        * SQRT(2.*PI*BoltzmannConst*CollInf%Tref(iSpec,iSpec)*(Species(iSpec)%MassIC + Species(iSpec2)%MassIC) &
        /(Species(iSpec)%MassIC * Species(iSpec2)%MassIC))/CollInf%Tref(iSpec,iSpec)**(-CollInf%omega(iSpec,iSpec) +0.5)
  END DO
END DO

FPCollModel = GETINT('Particles-FP-CollModel')
ESFPModel = GETINT('Particles-ESFP-Model')
DoBGKCellAdaptation = GETLOGICAL('Particles-FP-DoCellAdaptation')
IF(DoBGKCellAdaptation) THEN
  BGKMinPartPerCell = GETINT('Particles-FP-MinPartsPerCell')
  IF(.NOT.DSMC%UseOctree) THEN
    DSMC%UseOctree = .TRUE.
    IF(NGeo.GT.PP_N) CALL abort(&
__STAMP__&
,' Set PP_N to NGeo, otherwise the volume is not computed correctly.')
    CALL DSMC_init_octree()
  END IF
END IF

IF(DSMC%CalcQualityFactors) THEN
  ALLOCATE(FP_QualityFacSamp(1:6,nElems))
  FP_QualityFacSamp(1:6,1:nElems) = 0.0
END IF

FPDoVibRelaxation = GETLOGICAL('Particles-FP-DoVibRelaxation')
FPUseQuantVibEn = GETLOGICAL('Particles-FP-UseQuantVibEn')
CoupledFPDSMC = GETLOGICAL('Particles-CoupledFPDSMC')
IF(CoupledFPDSMC) THEN
  IF (DoVirtualCellMerge) THEN  
    CALL abort(__STAMP__,' Virtual cell merge not implemented for coupled DSMC-FP simulations!')
  END IF
  FPDSMCSwitchDens = GETREAL('Particles-FP-DSMC-SwitchDens')
END IF

FPInitDone = .TRUE.
LBWRITE(UNIT_stdOut,'(A)') ' INIT FP-FLOW DONE!'

END SUBROUTINE InitFPFlow

SUBROUTINE FinalizeFPFlow()
!===================================================================================================================================
!> Deallocating FP-specific variables, for a complete finalization FinalizeDSMC is also required
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_FPFlow_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(SpecFP)
SDEALLOCATE(FP_QualityFacSamp)

END SUBROUTINE FinalizeFPFlow

END MODULE MOD_FPFlow_Init
