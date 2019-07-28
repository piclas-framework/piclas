!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
PUBLIC :: Evaporation
PUBLIC :: SurfaceModel_main
PUBLIC :: UpdateSurfModelVars
!===================================================================================================================================

CONTAINS

SUBROUTINE SurfaceModel_main()
!===================================================================================================================================
!> Main Routine treating all surface and calculating desorbing / evaporating number of particles
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars ,ONLY: PartSurfaceModel, SolidSimFlag, LiquidSimFlag
USE MOD_SMCR          ,ONLY: SMCR_PartDesorb
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
!first do solid surfaces
IF (SolidSimFlag) THEN
  SELECT CASE (PartSurfaceModel)
  CASE (1)
    CALL Calc_DesorbPartNum()
  CASE (3)
    CALL SMCR_PartDesorb()
    !CALL AnalyzePartitionTemp()
  END SELECT
END IF

!then do liquid surfaces
IF (LiquidSimFlag) CALL Evaporation()

END SUBROUTINE SurfaceModel_main

SUBROUTINE UpdateSurfModelVars()
!===================================================================================================================================
!> Update and sample surface values for adsorption, desorption and reactions on surfaces (heterogeneous reactions)
!> Only procs with surface enter function, all other exit routine
!> 1. Communicate number of particles that were adsorbed on halo-sides of neighbour procs,
!>    so that own proc has number of total adsorption particles for each species on own surfaces
!> 2. After communication go through all own sides (no mpi halo sides) and adjust coverage resulting from changes through
!>    adsorption and reactions for all surfacemodels
!> 2.1  For surfacemodel=3 calculate number of adsorbate change (surfacempf!=gasmpf) and if changed Call AdjustReconstructMapNum
!>      to adjust number of adsorbates on reconstructed Monte Carlo surface
!> 3. Sample macroscopic surface coverage values
!> 4. Reinitialized surface reaction counters
!> 5. Calculated global adsorption probabilities for surfacemodel=1,2
!> 6. Send Coverages and Distribution info of own sides to halo sides of other procs (other procs have own sides as halo sides
!>    and need coverage info for adsorption calculation
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Particle_Vars          ,ONLY: WriteMacroSurfaceValues, KeepWallParticles, Species, nSpecies, PartSurfaceModel
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, SampWall, PartBound
USE MOD_TimeDisc_Vars          ,ONLY: tend,time
USE MOD_Mesh_Vars              ,ONLY: BC
#if USE_LOADBALANCE
USE MOD_LoadBalance_tools      ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption, SurfDistInfo
USE MOD_SurfaceModel_Tools     ,ONLY: CalcAdsorbProb, CalcDesorbProb
USE MOD_SurfaceModel_Tools     ,ONLY: SMCR_AdjustMapNum
#ifdef MPI
USE MOD_SurfaceModel_MPI       ,ONLY: ExchangeAdsorbNum, ExchangeCoverageInfo, ExchangeSurfDistInfo
#endif /*MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSpec, iSurfSide, p, q, new_adsorbates, numSites, SideID, PartboundID
REAL                             :: maxPart
REAL                             :: coverage_tmp, coverage_corrected
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!

IF (PartSurfaceModel.GT.0) THEN
  IF (.NOT.SurfMesh%SurfOnProc) RETURN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  IF (.NOT.KeepWallParticles) THEN
#ifdef MPI
! 1. communicate number of particles that were adsorbed on halo-sides of neighbour procs
    CALL ExchangeAdsorbNum()
#if USE_LOADBALANCE
    CALL LBSplitTime(LB_SURFCOMM,tLBStart)
#endif /*USE_LOADBALANCE*/
#endif /*MPI*/
    ! adjust coverages of all species on surfaces
    DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
      IF (DSMC%ReservoirRateStatistic) Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = 0
#endif
! 2. go through all sides and adjust coverages
      DO iSurfSide = 1,SurfMesh%nSides
        SideID = Adsorption%SurfSideToGlobSideMap(iSurfSide)
        PartboundID = PartBound%MapToPartBC(BC(SideID))
        IF (.NOT.PartBound%SolidReactive(PartboundID)) CYCLE
        DO q = 1,nSurfSample
          DO p = 1,nSurfSample
#if (PP_TimeDiscMethod==42)
            ! write number of adsorbates into info array
            IF (DSMC%ReservoirRateStatistic) THEN
              maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
              Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount = Adsorption%AdsorpInfo(iSpec)%WallSpecNumCount &
                                                            + INT( Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                                                            * maxPart/Species(iSpec)%MacroParticleFactor)
            END IF
#endif
            SELECT CASE (PartSurfaceModel)
            CASE(1)
              maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
              Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  + ( Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                  - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                  * Species(iSpec)%MacroParticleFactor / maxPart
            CASE(2)
              Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  + Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec)
            CASE(3)
              maxPart = REAL(INT(Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide),8))
              coverage_tmp = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  + REAL(( Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec) &
                  - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                  * Species(iSpec)%MacroParticleFactor) / maxPart
              IF (coverage_tmp.LT.0.) THEN ! can only happen for (ER + desorption) or (ER + MPF_surf<MPF_gas) --> SumAdsorbPart<0
                coverage_corrected = 0. !Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                    !+ REAL(( - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) ) &
                    !* Species(iSpec)%MacroParticleFactor) / maxPart
              ELSE
                coverage_corrected = coverage_tmp
              END IF
! 2.1  adjust number of mapped adsorbates on reconstructed surface if SumAdsorbPart > 0
              IF (Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec).NE.0) THEN
                ! calculate number of adsorbed particles on reconstructed surface for each species
                numSites = SurfDistInfo(p,q,iSurfSide)%nSites(3) !number of simulated surface atoms
                SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                      + (REAL(Adsorption%SumAdsorbPart(p,q,iSurfSide,iSpec)) * Species(iSpec)%MacroParticleFactor &
                      / maxPart) * REAL(numSites)
                ! convert to integer adsorbates
                new_adsorbates = INT(SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec))
                IF (new_adsorbates.NE.0) THEN
                  ! Adjust tracking of adsorbing background particles
                  CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec)
                  SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                                                                    - new_adsorbates
                END IF
              ELSE ! additional check if coverage changed but adsorbnum_tmp is not zero while sumadsorbpart is
                ! convert to integer adsorbates
                new_adsorbates = INT(SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec))
                IF (new_adsorbates.NE.0) THEN
                  ! Adjust tracking of adsorbing background particles
                  CALL SMCR_AdjustMapNum(p,q,iSurfSide,new_adsorbates,iSpec)
                  SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) = SurfDistInfo(p,q,iSurfSide)%adsorbnum_tmp(iSpec) &
                                                                    - new_adsorbates
                END IF
              END IF ! SumAdsorbPart!=0
              Adsorption%Coverage(p,q,iSurfSide,iSpec) = coverage_corrected
            END SELECT
! 3. sample adsorption coverage
            IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd))&
                .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
              SampWall(iSurfSide)%Adsorption(1+iSpec,p,q) = SampWall(iSurfSide)%Adsorption(1+iSpec,p,q) &
                                                          + Adsorption%Coverage(p,q,iSurfSide,iSpec)
            END IF
          END DO
        END DO
      END DO
    END DO
! 4. Reinitialize surface reaction counters
    IF (PartSurfaceModel.GE.2) THEN
      Adsorption%SumDesorbPart(:,:,:,:) = Adsorption%SumERDesorbed(:,:,1:SurfMesh%nSides,:)
      Adsorption%SumERDesorbed(:,:,:,:) = 0
    ELSE
      Adsorption%SumDesorbPart(:,:,:,:) = 0
    END IF
    Adsorption%SumAdsorbPart(:,:,:,:) = 0
    Adsorption%SumReactPart(:,:,:,:) = 0
  END IF
! 5. calculate probabiities
  IF (PartSurfaceModel.EQ.1) THEN
    CALL CalcAdsorbProb()
    IF (KeepWallParticles) CALL CalcDesorbprob()
  ELSE IF (PartSurfaceModel.EQ.2) THEN
    CALL CalcDesorbProb()
    CALL CalcAdsorbProb()
  END IF
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_SURF,tLBStart)
#endif /*USE_LOADBALANCE*/

! 6. communicate surface state to halo sides of neighbours
#ifdef MPI
  ! communicate coverage and probabilities to halo sides of neighbour procs
  IF (PartSurfaceModel.EQ.2) CALL ExchangeCoverageInfo()
  ! communicate distribution to halo-sides of neighbour procs
  IF (PartSurfaceModel.EQ.3) CALL ExchangeSurfDistInfo()
#endif

END IF

END SUBROUTINE UpdateSurfModelVars


SUBROUTINE Calc_DesorbPartNum()
!===================================================================================================================================
!> calculation of desorbing particle number when mean desorption probabilities are used (surfacemodel 1)
!===================================================================================================================================
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species, PartSurfaceModel
USE MOD_SurfaceModel_Vars      ,ONLY: Adsorption
USE MOD_SurfaceModel_Tools     ,ONLY: CalcDesorbProb
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                          :: iSurfSide, iSpec, p, q, NPois, WallPartNum
REAL                             :: PartAds, PartDes, RanNum, Tpois
#if USE_LOADBALANCE
INTEGER                          :: globSide, ElemID
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IF (.NOT.SurfMesh%SurfOnProc) RETURN
#if (PP_TimeDiscMethod==42)
Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
Adsorption%AdsorpInfo(:)%NumOfDes = 0
#endif
CALL CalcDesorbProb()

DO iSpec = 1,nSpecies
  DO iSurfSide = 1,SurfMesh%nSides
#if USE_LOADBALANCE
    IF(PerformLBSample) THEN
      globSide = Adsorption%SurfSideToGlobSideMap(iSurfSide)
      ElemID = PartSideToElem(S2E_ELEM_ID,globSide)
      nSurfacePartsPerElem(ElemID) = nSurfacePartsPerElem(ElemID) + 1
    END IF
#endif /*USE_LOADBALANCE*/
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        IF (Adsorption%Coverage(p,q,iSurfSide,iSpec).GT.0) THEN
          WallPartNum = INT(Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                    * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                    / Species(iSpec)%MacroParticleFactor)
          IF (WallPartNum .GT. 0) THEN
          PartAds = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                      / Species(iSpec)%MacroParticleFactor

          IF (PartSurfaceModel.EQ.1) THEN
            PartDes = PartAds * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
          END IF
            IF (PartDes.GT.PartAds) THEN
              Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) + INT(PartAds)
            ELSE
              CALL RANDOM_NUMBER(RanNum)
              IF (EXP(-PartDes).LE.TINY(PartDes) &
#if (PP_TimeDiscMethod==42)
              .OR. Adsorption%TPD &
#endif
              ) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                +INT(PartDes + RanNum)
              ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
                Npois=0
                Tpois=1.0
                DO
                  Tpois=RanNum*Tpois
                  IF (Tpois.LT.TINY(Tpois)) THEN
                    Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                    +INT(PartDes + RanNum)
                    EXIT
                  END IF
                  IF (Tpois.GT.EXP(-PartDes)) THEN
                    Npois=Npois+1
                    CALL RANDOM_NUMBER(RanNum)
                  ELSE
                    Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)+Npois
                    EXIT
                  END IF
                END DO
              END IF
            END IF !PartDes.GT.WallPartNum
            IF (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).GT.WallPartNum  &
            .OR.Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).LT.0) THEN
              Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = WallPartNum
            END IF
          ELSE !not PartAds.GT.0
            Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = 0
          END IF !PartAds.GT.0
#if (PP_TimeDiscMethod==42)
          Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                                + Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)
#endif
        END IF
      END DO
    END DO
  END DO
END DO

END SUBROUTINE Calc_DesorbPartNum


SUBROUTINE Evaporation()
!===================================================================================================================================
!> calculation of evaporating particle number when particles are deleted at condensation and inserted at evaporation
!===================================================================================================================================
USE MOD_Particle_Vars          ,ONLY: nSpecies, Species
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, PI
USE MOD_SurfaceModel_Vars      ,ONLY: Liquid
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample, SurfMesh, PartBound
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_TimeDisc_Vars          ,ONLY: dt
#if USE_LOADBALANCE
USE MOD_Particle_Mesh_Vars     ,ONLY: PartSideToElem
USE MOD_LoadBalance_Vars       ,ONLY: nSurfacePartsPerElem, PerformLBSample
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
INTEGER                          :: iSurfSide, iSpec, p, q, Npois
REAL                             :: PartEvap, RanNum, Tpois
REAL                             :: LiquidSurfTemp
REAL                             :: pressure_vapor, A, B, C
#if USE_LOADBALANCE
INTEGER                          :: globSide, ElemID
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
IF (.NOT.SurfMesh%SurfOnProc) RETURN
#if (PP_TimeDiscMethod==42)
Liquid%Info(:)%NumOfDes = 0
#endif
Liquid%SumEvapPart(:,:,:,:) = 0
DO iSpec = 1,nSpecies
  DO iSurfSide = 1,SurfMesh%nSides
    IF (PartBound%SolidState(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))) CYCLE
    IF (PartBound%LiquidSpec(PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) ))).NE.iSpec) CYCLE
#if USE_LOADBALANCE
    IF(PerformLBSample) THEN
      globSide = SurfMesh%SurfIDToSideID(iSurfSide)
      ElemID = PartSideToElem(S2E_ELEM_ID,globSide)
      nSurfacePartsPerElem(ElemID) = nSurfacePartsPerElem(ElemID) + 1
    END IF
#endif /*USE_LOADBALANCE*/
    LiquidSurfTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(SurfMesh%SurfIDToSideID(iSurfSide))))
    DO q = 1,nSurfSample
      DO p = 1,nSurfSample
        ! Antoine parameters defined in ini file are chosen so pressure given is in bar
        A = PartBound%ParamAntoine(1,PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))
        B = PartBound%ParamAntoine(2,PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))
        C = PartBound%ParamAntoine(3,PartBound%MapToPartBC(BC( SurfMesh%SurfIDToSideID(iSurfSide) )))
        ! Use Antoine Eq. to calculate pressure vapor
        pressure_vapor = 10 ** (A- B/(C+LiquidSurfTemp)) * 1e5 !transformation bar -> Pa
        ! Use Hertz-Knudsen equation to calculate number of evaporating liquid particles from surface
        PartEvap = pressure_vapor / ( 2*PI*Species(iSpec)%MassIC*BoltzmannConst*LiquidSurfTemp)**0.5 &
                 * SurfMesh%SurfaceArea(p,q,iSurfSide) / Species(iSpec)%MacroParticleFactor * dt
        CALL RANDOM_NUMBER(RanNum)
        IF (EXP(-PartEvap).LE.TINY(PartEvap)) THEN
          Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)&
          +INT(PartEvap + RanNum)
        ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
          Npois=0
          Tpois=1.0
          DO
            Tpois=RanNum*Tpois
            IF (Tpois.LT.TINY(Tpois)) THEN
              Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)&
              +INT(PartEvap + RanNum)
              EXIT
            END IF
            IF (Tpois.GT.EXP(-PartEvap)) THEN
              Npois=Npois+1
              CALL RANDOM_NUMBER(RanNum)
            ELSE
              Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)+Npois
              EXIT
            END IF
          END DO
        END IF
#if (PP_TimeDiscMethod==42)
        Liquid%Info(iSpec)%NumOfDes = Liquid%Info(iSpec)%NumOfDes + Liquid%SumEvapPart(p,q,iSurfSide,iSpec)
#endif
      END DO
    END DO
  END DO
END DO

END SUBROUTINE Evaporation


END MODULE MOD_SurfaceModel
