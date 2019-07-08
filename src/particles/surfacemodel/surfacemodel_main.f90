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

MODULE MOD_SurfaceModel
!===================================================================================================================================
!> Main Module for surface model
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
PUBLIC :: SurfaceModel_main
PUBLIC :: UpdateSurfModelVars
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceModel_main()
!===================================================================================================================================
!> Main Routine treating all surface and calculating desorbing / evaporating number of particles
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_SMCR                   ,ONLY: SMCR_PartDesorb, SMCR_Diffusion
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.ANY(PartBound%Reactive)) RETURN
! evaporating particles for surface models using mean values on surfaces
CALL CalcEvapPartNum()
! desorbing particles for surface model 3
CALL SMCR_Diffusion()
CALL SMCR_PartDesorb()

END SUBROUTINE SurfaceModel_main


SUBROUTINE UpdateSurfModelVars()
!===================================================================================================================================
!> Update and sample surface values for adsorption, desorption and reactions on surfaces (heterogeneous reactions)
!> Only procs with surface enter function, all other exit routine
!> 1. Communicate number of particles that were absorbed by halo-sides of neighbour procs (SumAdsorbPart and SumERDesorbed),
!>    so that own proc has number of total adsorption particles for each species on own surfaces
!> 2. After communication go through all own sides (no mpi halo sides) and adjust coverage resulting from changes through
!>    adsorption and reactions for all surfacemodels
!> 2.1  For surfacemodel=3 calculate number of adsorbate change (surfacempf!=gasmpf) and if changed Call AdjustReconstructMapNum
!>      to adjust number of adsorbates on reconstructed Monte Carlo surface
!> 2.2  Sample macroscopic surface coverage values
!> 3. Reset/Adjust surface model sum counters
!> 4. Calculate global mean probabilities for surface models
!> 5. Send Coverages and Distribution info of own sides to halo sides of other procs (other procs have own sides as halo sides
!>    and need coverage info for adsorption calculation (mpi routines know which sides to communicate)
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues, KeepWallParticles, Species, nSpecies
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall, PartBound
USE MOD_TimeDisc_Vars          ,ONLY: tend,time
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools      ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo, SurfModel
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbProb, CalcDesorbProb
USE MOD_SurfaceModel_Tools     ,ONLY: SMCR_AdjustMapNum, IsReactiveSurface, SurfaceHasModelNum
#ifdef MPI
USE MOD_SurfaceModel_MPI       ,ONLY: ExchangeAdsorbNum, ExchangeCoverageInfo, ExchangeSurfDistInfo
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSpec, iSurfSide, p, q, new_adsorbates, numSites
REAL                             :: maxPart
REAL                             :: coverage_tmp, coverage_corrected
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
#if (PP_TimeDiscMethod==42)
REAL                             :: desorbnum_covreduce
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.SurfMesh%SurfOnProc) RETURN
IF (.NOT.ANY(Partbound%Reactive)) RETURN

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
IF (.NOT.KeepWallParticles) THEN
#ifdef MPI
!----- 1.
  CALL ExchangeAdsorbNum()
#if USE_LOADBALANCE
  CALL LBSplitTime(LB_SURFCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*MPI*/

  ! adjust coverages of all species on surfaces
  DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
    IF (DSMC%ReservoirRateStatistic) SurfModel%Info(iSpec)%WallSpecNumCount = 0
#endif
!----- 2.
    DO iSurfSide = 1,SurfMesh%nSides
      IF (.NOT.IsReactiveSurface(iSurfSide)) CYCLE
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
#if (PP_TimeDiscMethod==42)
          ! write number of adsorbates into info array
          IF (DSMC%ReservoirRateStatistic) THEN
            maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
            SurfModel%Info(iSpec)%WallSpecNumCount = SurfModel%Info(iSpec)%WallSpecNumCount &
                                                          + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                                                          * maxPart/Species(iSpec)%MacroParticleFactor)
          END IF
#endif
          SELECT CASE (SurfaceHasModelNum(iSurfSide))
          CASE(1)
            maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + ( SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                * Species(iSpec)%MacroParticleFactor / maxPart
          CASE(2)
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec)
          CASE(3)
#if (PP_TimeDiscMethod==42)
            IF (Adsorption%CoverageReduction) THEN
              ! calculate number of desorbed particles for each species
              desorbnum_covreduce = (REAL(Adsorption%CovReductionStep(iSpec))/REAL(SurfDistInfo(p,q,iSurfSide)%nSites(3))) &
                  * REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) &
                  * SurfMesh%SurfaceArea(p,q,iSurfSide),8)) / Species(iSpec)%MacroParticleFactor
              coverage_tmp = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  - REAL(desorbnum_covreduce) * Species(iSpec)%MacroParticleFactor &
                  / REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
              IF(coverage_tmp.GT.0.) THEN
                Adsorption%Coverage(p,q,iSurfSide,iSpec) = coverage_tmp
                new_adsorbates = -1*Adsorption%CovReductionStep(iSpec)
                CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec)
              END IF
            END IF
#endif
            maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
            coverage_tmp = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                + REAL(( SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                * Species(iSpec)%MacroParticleFactor) / maxPart
            IF (coverage_tmp.LT.0.) THEN ! can only happen for (ER + desorption) or (ER + MPF_surf<MPF_gas) --> SumAdsorbPart<0
              coverage_corrected = 0. !Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  !+ REAL(( - (SurfModel%SumDesorbPart(p,q,iSurfSide,iSpec) - SurfModel%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                  !* Species(iSpec)%MacroParticleFactor) / maxPart
            ELSE
              coverage_corrected = coverage_tmp
            END IF
!----- 2.1
            IF (SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec).NE.0) THEN
              ! calculate number of adsorbed particles on reconstructed surface for each species
              numSites = SurfDistInfo(p,q,iSurfSide)%nSites(3) !number of simulated surface atoms
              SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                    + (REAL(SurfModel%SumAdsorbPart(p,q,iSurfSide,iSpec)) * Species(iSpec)%MacroParticleFactor &
                    / maxPart) * REAL(numSites)
            END IF ! SumAdsorbPart!=0
            ! convert to integer adsorbates
            new_adsorbates = INT(SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec))
            IF (new_adsorbates.NE.0) THEN
              ! Adjust tracking of adsorbing background particles
              CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec,SampleFlag=.TRUE.)
              SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                                                                - new_adsorbates
            END IF
            Adsorption%Coverage(p,q,iSurfSide,iSpec) = coverage_corrected
          END SELECT
!----- 2.2
          IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
              .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
            SampWall(iSurfSide)%SurfModelState(5+iSpec,p,q) = SampWall(iSurfSide)%SurfModelState(5+iSpec,p,q) &
                                                        + Adsorption%Coverage(p,q,iSurfSide,iSpec)
          END IF
        END DO
      END DO
    END DO
  END DO
!----- 3.
  ! SumEvapPart is nullified in particle emission (surface flux) after inserting particles at corresponding surfaces
  SurfModel%SumEvapPart(:,:,:,:) = SurfModel%SumEvapPart(:,:,:,:) + SurfModel%SumERDesorbed(:,:,1:SurfMesh%nSides,:)
  SurfModel%SumERDesorbed(:,:,:,:) = 0
  SurfModel%SumDesorbPart(:,:,:,:) = 0
  SurfModel%SumAdsorbPart(:,:,:,:) = 0
  SurfModel%SumReactPart(:,:,:,:) = 0
END IF
Adsorption%NumCovSamples = Adsorption%NumCovSamples + 1
!----- 4.
CALL CalcDesorbProb()
CALL CalcAdsorbProb()
#if USE_LOADBALANCE
CALL LBPauseTime(LB_SURF,tLBStart)
#endif /*USE_LOADBALANCE*/

!----- 5.
#ifdef MPI
! communicate coverage and probabilities to halo sides of neighbour procs
CALL ExchangeCoverageInfo()
! communicate distribution to halo-sides of neighbour procs
CALL ExchangeSurfDistInfo()
#endif

END SUBROUTINE UpdateSurfModelVars


SUBROUTINE CalcEvapPartNum()
!===================================================================================================================================
!> calculation of number of evaporating/desorbing particles when mean surface densities are used (mean probabilities)
!===================================================================================================================================
USE MOD_Globals_Vars           ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, surfmodel, SpecSurf
USE MOD_SurfaceModel_Tools     ,ONLY: SurfaceHasModelNum
USE MOD_Particle_Boundary_Tools,ONLY: TSURUTACONDENSCOEFF
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_TimeDisc_Vars          ,ONLY: dt
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSurfSide, iSpec, p, q, NPois, PartEvapInfo
REAL                             :: WallPartNum, PartEvap, RanNum, Tpois, LiquidSurfTemp
REAL                             :: A, B, C, pressureVapor
INTEGER                          :: iPart, PartEvap_temp
REAl                             :: veloPart, sigma
#if USE_LOADBALANCE
INTEGER                          :: ElemID
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.SurfMesh%SurfOnProc) RETURN

DO iSpec = 1,nSpecies
  DO iSurfSide = 1,SurfMesh%nSides
#if USE_LOADBALANCE
    IF(PerformLBSample) THEN
      ElemID = PartSideToElem(S2E_ELEM_ID,SurfMesh%SurfIDToSideID(iSurfSide))
      nSurfacePartsPerElem(ElemID) = nSurfacePartsPerElem(ElemID) + 1
    END IF
#endif /*USE_LOADBALANCE*/
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        SELECT CASE (SurfaceHasModelNum(iSurfSide))
        CASE(1)
          IF (Adsorption%Coverage(p,q,iSurfSide,iSpec).GT.0) THEN
            WallPartNum = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                      / Species(iSpec)%MacroParticleFactor
            PartEvap = WallPartNum * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
          ELSE
            PartEvap = 0.
          END IF
        CASE(101,102)
          IF (Adsorption%SurfaceSpec(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )),iSpec)) THEN
            LiquidSurfTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))
            ! Antoine parameters defined in ini file are chosen so pressure given is in bar
            A = SpecSurf(iSpec)%ParamAntoine(1)
            B = SpecSurf(iSpec)%ParamAntoine(2)
            C = SpecSurf(iSpec)%ParamAntoine(3)
            ! Use Antoine Eq. to calculate pressure vapor
            pressureVapor = 10 ** (A- B/(C+LiquidSurfTemp)) * 1e5 !transformation bar -> Pa
            ! Use Hertz-Knudsen equation to calculate number of evaporating liquid particles from surface
            PartEvap = pressureVapor / ( 2*PI*Species(iSpec)%MassIC*BoltzmannConst*LiquidSurfTemp)**0.5 &
                     * SurfMesh%SurfaceArea(p,q,iSurfSide) / Species(iSpec)%MacroParticleFactor * dt &
                     * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
            IF (SurfaceHasModelNum(iSurfSide).EQ.102) THEN
              sigma=SQRT(BoltzmannConst*LiquidSurfTemp/Species(iSpec)%MassIC)
              PartEvap_temp = PartEvap
              PartEvap = 0
              DO iPart = 1,PartEvap_temp
                CALL RANDOM_NUMBER(RanNum)
                veloPart = SQRT(-2*LOG(RanNum))*sigma
                CALL RANDOM_NUMBER(RanNum)
                IF ( (TSURUTACONDENSCOEFF(iSpec,veloPart,LiquidSurfTemp).GE.RanNum) ) PartEvap=PartEvap+1
              END DO
            END IF
            WallPartNum = PartEvap
          ELSE
            WallPartNum = 0
          END IF
        CASE DEFAULT
          WallPartNum = 0
        END SELECT

        IF (WallPartNum.GT.0) THEN
          IF (PartEvap.GT.WallPartNum .OR. PartEvap.LT.0.) THEN
            PartEvapInfo = INT(WallPartNum)
          ELSE
            CALL RANDOM_NUMBER(RanNum)
            IF (EXP(-PartEvap).LE.TINY(PartEvap) &
#if (PP_TimeDiscMethod==42)
                .OR. Adsorption%TPD &
#endif
                ) THEN
              PartEvapInfo = INT(PartEvap + RanNum)
            ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
              Npois=0
              Tpois=1.0
              DO
                Tpois=RanNum*Tpois
                IF (Tpois.LT.TINY(Tpois)) THEN
                  PartEvapInfo = INT(PartEvap + RanNum)
                  EXIT
                END IF
                IF (Tpois.GT.EXP(-PartEvap)) THEN
                  Npois=Npois+1
                  CALL RANDOM_NUMBER(RanNum)
                ELSE
                  PartEvapInfo = Npois
                  EXIT
                END IF
              END DO
            END IF
          END IF !PartEvap.GT.WallPartNum
          SurfModel%SumEvapPart(p,q,iSurfSide,iSpec) = SurfModel%SumEvapPart(p,q,iSurfSide,iSpec) + INT(PartEvapInfo)
#if (PP_TimeDiscMethod==42)
          SurfModel%Info(iSpec)%NumOfDes = SurfModel%Info(iSpec)%NumOfDes + INT(PartEvapInfo)
#endif
        END IF !WallPartNum.GT.0
      END DO
    END DO
  END DO
END DO

END SUBROUTINE CalcEvapPartNum


END MODULE MOD_SurfaceModel
