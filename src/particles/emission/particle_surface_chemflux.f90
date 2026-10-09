!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_Particle_SurfChemFlux
!===================================================================================================================================
!> Module for particle insertion through the surface flux
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: PureSurfChemistry, ParticleSurfDiffusion
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> Particle insertion by pure surface reactions (independent of gas-collisions on the surface)
!> 1.) Determine the surface parameters
!> 2.) Calculate the number of newly created products and update the surface properties
!>  a) Langmuir-Hinshelwood reaction with instantaneous desorption (Arrhenius model)
!>  b) Langmuir-Hinshelwood reaction (Arrhenius model)
!>  c) Thermal desorption (Polanyi-Wigner equation)
!> 3.) Insert the product species into the gas phase
!===================================================================================================================================
SUBROUTINE PureSurfChemistry()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Part_Tools              ,ONLY: CalcRadWeightMPF, CalcVarWeightMPF, GetNextFreePosition
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Pure   ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, GlobalSide2SurfSide, SurfSideArea, SurfTotalSideOnNode
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, SurfChemReac, ChemWallProp
USE MOD_DSMC_Vars               ,ONLY: DSMC,SamplingActive
#if USE_MPI
USE MOD_MPI_Shared_vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                     :: iSpec , iSF, iSide, SideID
INTEGER                     :: BCSideID, ElemID, iLocSide, globElemId, CNElemID
REAL                        :: ReacHeat, DesHeat, RanNum, Area, AdsDens
REAL                        :: DesCount, TotalCoverage
REAL                        :: nu, E_act, Coverage, Rate, DissOrder, AdCount
REAL                        :: BetaCoeff, WallTemp, SurfMol, SurfMolDens
INTEGER                     :: iReac, iProd, ReactantCount, BoundID
INTEGER                     :: iVal, iReactant, iValReac, SurfSideID, iBias
INTEGER                     :: SubP, SubQ
REAL                        :: DesCount_Spec(1:nSpecies,1:nSurfSample,1:nSurfSample)
INTEGER                     :: SurfReacBias(SurfChem%NumOfReact)
INTEGER                     :: iShuffle, jShuffle, iTmp
LOGICAL                     :: DoReac, DoSampleHeat
!===================================================================================================================================
IF(.NOT.SurfTotalSideOnNode) RETURN

DoSampleHeat = (DSMC%CalcSurfaceVal.AND.(SamplingActive.OR.WriteMacroSurfaceValues))

DO iSF = 1, SurfChem%CatBoundNum
  BoundID = SurfChem%SurfacefluxBC(iSF)
  ! No pure surface reactions on this boundary
  IF(.NOT.ANY(SurfChem%PSMap(BoundID)%PureSurfReac)) CYCLE

  DO iSide = 1, BCdata_auxSF(BoundID)%SideNumber
    BCSideID   = BCdata_auxSF(BoundID)%SideList(iSide)
    ElemID     = SideToElem(S2E_ELEM_ID,BCSideID)
    iLocSide   = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
    globElemId = ElemID + offSetElem
    CNElemID   = GetCNElemID(globElemId)
    SideID     = GetGlobalNonUniqueSideID(globElemId,iLocSide)
    SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

    IF (SurfSideID.LT.1) CALL abort(__STAMP__,'Chemical Surface Flux is not allowed on non-sampling sides!')

    WallTemp      = PartBound%WallTemp(BoundID) ! Boundary temperature
    DesCount_Spec = 0.0

    ! Number of surface molecules PER AREA - independent of the sub-surface element
    IF(PartBound%LatticeVec(BoundID).GT.0.) THEN
      SurfMolDens = PartBound%MolPerUnitCell(BoundID)/(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
    ELSE
      ! Alternative: average number of surface molecules per area for a monolayer
      SurfMolDens = 10.**19
    END IF

    DO SubQ = 1, nSurfSample
      DO SubP = 1, nSurfSample

        Area    = SurfSideArea(SubP,SubQ,SurfSideID)
        SurfMol = SurfMolDens * Area

        DO iShuffle = 1, SurfChem%NumOfReact
          SurfReacBias(iShuffle) = iShuffle
        END DO
        DO iShuffle = SurfChem%NumOfReact, 2, -1
          CALL RANDOM_NUMBER(RanNum)
          jShuffle = 1 + INT(RanNum*REAL(iShuffle))
          IF (jShuffle.GT.iShuffle) jShuffle = iShuffle
          iTmp                   = SurfReacBias(iShuffle)
          SurfReacBias(iShuffle) = SurfReacBias(jShuffle)
          SurfReacBias(jShuffle) = iTmp
        END DO

        ! Loop over the different types of pure surface reactions
        DO iBias = 1, SurfChem%NumOfReact
          iReac = SurfReacBias(iBias)

          IF (.NOT.SurfChem%PSMap(BoundID)%PureSurfReac(iReac)) CYCLE

          TotalCoverage = SUM(ChemWallProp(1:nSpecies,SubP,SubQ,SurfSideID))
          ! Check for the correct coverage range of the reaction
          IF (SurfChemReac(iReac)%CovDep) THEN
            IF (TotalCoverage.LT.SurfChemReac(iReac)%MinCovTotal.OR.TotalCoverage.GT.SurfChemReac(iReac)%MaxCovTotal) CYCLE

            DoReac = .TRUE.
            DO iSpec=1, nSpecies
              IF (ChemWallProp(iSpec,SubP,SubQ,SurfSideID).LT.SurfChemReac(iReac)%MinCov(iSpec)) DoReac = .FALSE.
              IF (ChemWallProp(iSpec,SubP,SubQ,SurfSideID).GT.SurfChemReac(iReac)%MaxCov(iSpec)) DoReac = .FALSE.
            END DO
            IF (.NOT.DoReac) CYCLE
          END IF
          ! Check if the reaction only takes place in a given temperature range
          IF (SurfChemReac(iReac)%TempDep) THEN
            IF (WallTemp.LT.SurfChemReac(iReac)%MinTemp.OR.WallTemp.GT.SurfChemReac(iReac)%MaxTemp) CYCLE
          END IF

          ReactantCount = 1
          Coverage      = 0.
          AdCount       = 0.

          ! 2.) Calculate the number of newly created products and update the surface properties
          SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

          ! 2a) Langmuir-Hinshelwood reaction with instantaneous desorption (Arrhenius model)
          CASE('LHD')
            AdsDens  = 1.
            Coverage = 1.
            ! Product of the reactant coverage values
            DO iVal=1,SIZE(SurfChemReac(iReac)%Reactants(:))
              IF(SurfChemReac(iReac)%Reactants(iVal).GT.0) THEN
                iReactant = SurfChemReac(iReac)%Reactants(iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iReactant)
                IF(iReactant.NE.SurfChem%SurfSpecies) THEN ! Coverage set to 1 for surface species
                  ! Higher coverage due to bulk species still enters into the equation as 1
                  IF (ChemWallProp(iReactant,SubP,SubQ,SurfSideID).GT.1.) THEN
                    Coverage = Coverage * 1.
                    AdsDens  = AdsDens * 1. * SurfMolDens
                  ELSE
                    Coverage = Coverage * ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                    ! Particle density of the adsorbate on the surface
                    AdsDens  = AdsDens * ChemWallProp(iReactant,SubP,SubQ,SurfSideID) * SurfMolDens
                  END IF
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            BetaCoeff = SurfChemReac(iReac)%HeatAccommodation ! Incomplete energy accomodation
            ReacHeat  = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst
            nu        = SurfChemReac(iReac)%Prefactor
            E_act     = SurfChemReac(iReac)%ArrheniusEnergy
            ! Randomization of the reaction process
            CALL RANDOM_NUMBER(RanNum)

            ! Calculate the rate in dependence of the temperature and coverage
            Rate     = nu * AdsDens * exp(-E_act/WallTemp) ! Energy in K
            DesCount = Rate * dt * Area * (-LOG(RanNum))/ ReactantCount

            IF (DesCount.GT.0.) THEN
              DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                IF (SurfChemReac(iReac)%Reactants(iValReac).EQ.0) CYCLE
                iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                IF (iReactant.EQ.SurfChem%SurfSpecies) CYCLE
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iReactant)
                ! Number of available adsorbates
                AdCount = ChemWallProp(iReactant,SubP,SubQ,SurfSideID) * SurfMol
                ! Check if enough adsorbate reactants are available
                IF(DesCount.GT.AdCount/ReactantCount) DesCount = AdCount/ReactantCount
              END DO ! iValReac

              !2nd pass: consume the adsorbates. A species listed twice is visited twice and
              !           therefore decremented twice, which is the correct stoichiometry.
              DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                IF (SurfChemReac(iReac)%Reactants(iValReac).EQ.0) CYCLE
                iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                IF (iReactant.EQ.SurfChem%SurfSpecies) CYCLE
                ChemWallProp(iReactant,SubP,SubQ,SurfSideID) =  MAX(0.,ChemWallProp(iReactant,SubP,SubQ,SurfSideID) &
                                                             - DesCount/SurfMol)
              END DO ! iValReac

              ! Assign the desorption count to the product
              DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
                IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                  iProd = SurfChemReac(iReac)%Products(iVal)
                  DesCount_Spec(iProd,SubP,SubQ) = DesCount_Spec(iProd,SubP,SubQ) + DesCount
                END IF
              END DO ! iVal

              ! Update the catalytic heat flux onto the surface element
              IF (DoSampleHeat)ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) = ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) &
                                                            + DesCount*ReacHeat*BetaCoeff
            END IF

          ! b) Langmuir-Hinshelwood reaction (Arrhenius model)
          CASE('LH')
            Coverage = 1.
            AdsDens  = 1.
            ! Product of the reactant coverage values
            DO iVal=1,SIZE(SurfChemReac(iReac)%Reactants(:))
              IF(SurfChemReac(iReac)%Reactants(iVal).NE.0) THEN
                iReactant = SurfChemReac(iReac)%Reactants(iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iReactant)
                IF(iReactant.NE.SurfChem%SurfSpecies) THEN ! Coverage set to 1 for surface species
                  ! Higher coverage due to bulk species still enters into the equation as 1
                  IF (ChemWallProp(iReactant,SubP,SubQ,SurfSideID).GT.1.) THEN
                    Coverage = Coverage * 1.
                    AdsDens  = AdsDens * 1. * SurfMolDens
                  ELSE
                    Coverage = Coverage * ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                    ! Particle density of the adsorbate on the surface
                    AdsDens  = AdsDens * ChemWallProp(iReactant,SubP,SubQ,SurfSideID) * SurfMolDens
                  END IF
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            ! Complete accommodation due to the intermediate desorption step
            BetaCoeff = SurfChemReac(iReac)%HeatAccommodation
            ReacHeat  = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst
            nu        = SurfChemReac(iReac)%Prefactor
            E_act     = SurfChemReac(iReac)%ArrheniusEnergy

            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * AdsDens * exp(-E_act/WallTemp)! Energy in K

            ! Randomize the reaction process
            CALL RANDOM_NUMBER(RanNum)
            ! Reaction product number, here as a COVERAGE fraction rather than an absolute count
            DesCount = Rate * dt * Area * (-LOG(RanNum))/(SurfMol*ReactantCount)

            IF (DesCount.GT.0.) THEN
              ! Check if enough adsorbate particles are available for the process
              DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                IF (SurfChemReac(iReac)%Reactants(iValReac).EQ.0) CYCLE
                iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                IF (iReactant.EQ.SurfChem%SurfSpecies) CYCLE
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iReactant)
                Coverage = ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                IF(DesCount.GT.Coverage/ReactantCount) DesCount = Coverage/ReactantCount
              END DO ! iValReac

              ! Check if the maximum coverage of the products is reached after the reaction
              DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
                IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                  iProd = SurfChemReac(iReac)%Products(iVal)
                  ! Test for the maximum of the product coverage and set the desorpton count to the difference
                  IF((ChemWallProp(iProd,SubP,SubQ,SurfSideID)+DesCount).GT.PartBound%MaxCoverage(BoundID,iProd)) THEN
                    DesCount = PartBound%MaxCoverage(BoundID,iProd) - ChemWallProp(iProd,SubP,SubQ,SurfSideID)
                  END IF
                END IF ! iVal in Products
              END DO ! iVal

              ! Update the reactant coverage
              DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                IF (SurfChemReac(iReac)%Reactants(iValReac).EQ.0) CYCLE
                iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                IF (iReactant.EQ.SurfChem%SurfSpecies) CYCLE
                ChemWallProp(iReactant,SubP,SubQ,SurfSideID) =  MAX(0.,ChemWallProp(iReactant,SubP,SubQ,SurfSideID) - DesCount)
              END DO ! iValReac

              ! Update the product coverage
              DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
                IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                  iProd = SurfChemReac(iReac)%Products(iVal)
                  ChemWallProp(iProd,SubP,SubQ,SurfSideID) = ChemWallProp(iProd,SubP,SubQ,SurfSideID) + DesCount
                END IF ! iVal in Products
              END DO ! iVal

              ! Update the surface heatflux (DesCount is a coverage fraction here, hence the SurfMol factor)
              IF (DoSampleHeat)ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) = ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) &
                                                            + DesCount*ReacHeat*SurfMol*BetaCoeff
            END IF

          ! c) Thermal desorption (Polanyi-Wigner equation)
          CASE('D')
            ! Number of adsorbed particles on the subside.
            DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
              IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                  ! Higher coverage due to bulk species still enters into the equation as 1
                  IF (ChemWallProp(iReactant,SubP,SubQ,SurfSideID).GT.1.) THEN
                    Coverage = 1.
                  ELSE
                    Coverage = ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                  END IF
                ELSE
                  Coverage = 1.
                END IF
                AdCount = Coverage * SurfMol
              END IF
            END DO

            ! Absolute particle density of the surface element
            AdsDens = Coverage * SurfMolDens

            ! Calculate the desorption energy in dependence of the coverage [J]
            BetaCoeff = SurfChemReac(iReac)%HeatAccommodation
            DesHeat   = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst

            ! Define the variables
            DissOrder = SurfChemReac(iReac)%DissOrder
            nu        = SurfChemReac(iReac)%Prefactor
            ! Calculate the desorption prefactor in dependence of coverage and temperature of the boundary
            IF(nu.EQ.0.) THEN
              nu = 10.**(SurfChemReac(iReac)%C_a + SurfChemReac(iReac)%C_b * Coverage)
              IF (DissOrder.EQ.2) THEN
                ! Convert the prefactor to coverage values for the associative desorption
                nu = 10.**(SurfChemReac(iReac)%C_a + SurfChemReac(iReac)%C_b * Coverage) *10.**(15)
              END IF
            END IF

            ! Randomization of the desorption process
            CALL RANDOM_NUMBER(RanNum)

            E_act = SurfChemReac(iReac)%E_initial + Coverage * SurfChemReac(iReac)%W_interact
            ! Reaction rate according to the Polanyi-Wigner equation
            Rate     = nu * AdsDens**DissOrder * exp(-E_act/WallTemp)  ! Energy in K
            DesCount = Rate * dt * Area * (-LOG(RanNum)) /DissOrder

            IF (DesCount.GT.0.) THEN
              ! Upper bound for the desorption number
              IF(DesCount.GE.(AdCount/DissOrder)) DesCount = AdCount/DissOrder

              DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                  iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                  IF (iReactant.NE.SurfChem%SurfSpecies) THEN
                    ChemWallProp(iReactant,SubP,SubQ,SurfSideID) =  MAX(0.,ChemWallProp(iReactant,SubP,SubQ,SurfSideID) &
                                                                 - DissOrder*DesCount/SurfMol)
                  END IF
                END IF
              END DO

              ! Assign the desorption count to the product
              DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
                IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                  iProd = SurfChemReac(iReac)%Products(iVal)
                  DesCount_Spec(iProd,SubP,SubQ) = DesCount_Spec(iProd,SubP,SubQ) + DesCount
                END IF
              END DO ! iVal

              ! Update the surface heatflux (desorption consumes heat)
              IF (DoSampleHeat)ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) = ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) &
                                                            - DesCount*DesHeat*BetaCoeff
            END IF

          CASE DEFAULT
            ! Reaction types not handled here ('A','ER','P') are treated during the particle-wall interaction
          END SELECT
        END DO !iBias

      END DO ! SubP
    END DO ! SubQ

    ! Insert the calculated number of product species into the gas phase, once for the whole side
    CALL ParticleSurfChemFlux(DesCount_Spec,BoundID,iSide,iSF)

  END DO !iSide
END DO !iSF

#if USE_MPI
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE PureSurfChemistry

!===================================================================================================================================
!> Particle insertion by pure surface reactions (independent of gas-collisions on the surface)
!> Insertion of the determined number of product species for each of the surface reactions
!===================================================================================================================================
SUBROUTINE ParticleSurfChemFlux(DesCount_Spec,BoundID,iSide,iSF)
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars
USE MOD_Part_Tools              ,ONLY: CalcRadWeightMPF, CalcVarWeightMPF, GetNextFreePosition
USE MOD_DSMC_Vars               ,ONLY: CollisMode, DoRadialWeighting, DoLinearWeighting, DoCellLocalWeighting, DSMC
USE MOD_DSMC_Vars               ,ONLY: SamplingActive
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemMidPoint_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Pure   ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize, SurfMeshSubSideData
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide, nSurfSample, SurfSideSamplingMidPoints
USE MOD_Particle_SurfFlux       ,ONLY: CalcPartPosTriaSurface, DefineSideDirectVec2D
USE MOD_SurfaceModel_Vars       ,ONLY: ChemDesorpWall
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr
USE MOD_Particle_Boundary_Tools ,ONLY: CalcWallSample
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
#if USE_MPI
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: nSurfacefluxPerElem
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: DesCount_Spec(1:nSpecies,1:nSurfSample,1:nSurfSample)
INTEGER, INTENT(IN)              :: BoundID, iSide, iSF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                     :: iSpec, SideID, NbrOfParticle, PartID, SurfSideID
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, iPart, iPartTotal
INTEGER                     :: PartsEmitted, Node1, Node2, globElemId, CNElemID
REAL                        :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3), RVec(2), minPos(2)
REAL                        :: SurfElemMPF
INTEGER                     :: SubP, SubQ, PartInsSide, iTry, p, q, pFound, qFound, iSampleAcc, jSampleAcc
INTEGER,PARAMETER           :: MaxTries = 1000
REAL                        :: AreaTot, AreaCum, RanNum, NewPos(1:3), distance, distanceMin
LOGICAL                     :: PosAccepted
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

BCSideID   = BCdata_auxSF(BoundID)%SideList(iSide)
ElemID     = SideToElem(S2E_ELEM_ID,BCSideID)
iLocSide   = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
globElemId = ElemID + offSetElem
CNElemID   = GetCNElemID(globElemId)
SideID     = GetGlobalNonUniqueSideID(globElemId,iLocSide)
SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

! Geometry of the side, identical for all species and sub-surface elements
xyzNod(1:3) = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%xyzNod(1:3)
IF(Symmetry%Axisymmetric) CALL DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)

AreaTot = 0.
DO jSample = 1,SurfFluxSideSize(2)
  DO iSample = 1,SurfFluxSideSize(1)
    AreaTot = AreaTot + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
  END DO
END DO

! Current boundary condition
PartsEmitted  = 0
NbrOfParticle = 0
iPartTotal    = 0

! Insert the product species into the gas phase by a surface flux from the boundary
DO iSpec = 1, nSpecies
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

  ! MPF of the element
  IF (DoRadialWeighting) THEN
    SurfElemMPF = CalcRadWeightMPF(ElemMidPoint_Shared(2,CNElemID),iSpec,ElemID)
  ELSE IF (DoLinearWeighting.OR.DoCellLocalWeighting) THEN
    SurfElemMPF = CalcVarWeightMPF(ElemMidPoint_Shared(:,CNElemID), ElemID, iSpec)
  ELSE
    SurfElemMPF = Species(iSpec)%MacroParticleFactor
  END IF

  DO SubQ = 1, nSurfSample
    DO SubP = 1, nSurfSample

      ChemDesorpWall(iSpec,SubP,SubQ,SurfSideID) = ChemDesorpWall(iSpec,SubP,SubQ,SurfSideID) &
                                                 + DesCount_Spec(iSpec,SubP,SubQ)

      PartInsSide = INT(ChemDesorpWall(iSpec,SubP,SubQ,SurfSideID)/SurfElemMPF)
      IF (PartInsSide.LT.1) CYCLE

      ChemDesorpWall(iSpec,SubP,SubQ,SurfSideID) = ChemDesorpWall(iSpec,SubP,SubQ,SurfSideID) &
                                                 - REAL(PartInsSide)*SurfElemMPF

      !-- Fill Particle Informations (PartState, Partelem, etc.)
      PartID = 1
      DO iPart = 1,PartInsSide
        IF ((iPart.EQ.1).OR.PDM%ParticleInside(PartID)) PartID = GetNextFreePosition(iPartTotal+1)
        PosAccepted = .FALSE.
        iSampleAcc  = 1
        jSampleAcc  = 1
        DO iTry = 1,MaxTries
          ! Draw the sub-side triangle proportional to its area
          iSample = 1
          jSample = 1
          IF (AreaTot.GT.0.) THEN
            CALL RANDOM_NUMBER(RanNum)
            RanNum  = RanNum*AreaTot
            AreaCum = 0.
            TriaLoop: DO q = 1,SurfFluxSideSize(2)
              DO p = 1,SurfFluxSideSize(1)
                AreaCum = AreaCum + SurfMeshSubSideData(p,q,BCSideID)%area
                iSample = p
                jSample = q
                IF (RanNum.LT.AreaCum) EXIT TriaLoop
              END DO
            END DO TriaLoop
          END IF

          ! Position within that triangle
          Node1 = jSample+1
          Node2 = jSample+2
          Vector1       = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
          Vector2       = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
          midpoint(1:3) = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
          ndist(1:3)    = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)

          IF(Symmetry%Axisymmetric) THEN
            NewPos(1:3) = CalcPartPosAxisym(minPos, RVec)
          ELSE
            NewPos(1:3) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
          END IF

          iSampleAcc = iSample
          jSampleAcc = jSample

          ! Without a subdivision every position is valid
          IF (nSurfSample.EQ.1) THEN
            PosAccepted = .TRUE.
            EXIT
          END IF

          ! Same nearest-midpoint criterion that SurfaceModelling uses to assign TrackInfo%p/q
          distanceMin = HUGE(1.)
          pFound      = 1
          qFound      = 1
          DO q = 1,nSurfSample
            DO p = 1,nSurfSample
              distance = VECNORM3D(NewPos(1:3) - SurfSideSamplingMidPoints(1:3,p,q,SurfSideID))
              IF (distance.LT.distanceMin) THEN
                distanceMin = distance
                pFound      = p
                qFound      = q
              END IF
            END DO
          END DO

          IF ((pFound.EQ.SubP).AND.(qFound.EQ.SubQ)) THEN
            PosAccepted = .TRUE.
            EXIT
          END IF
        END DO ! iTry

        ! Fallback: place the particle on the midpoint of the sub-element. With an acceptance rate of about
        ! 1/nSurfSample**2 this is unreachable in practice and only guards against a degenerate geometry.
        IF (.NOT.PosAccepted) NewPos(1:3) = SurfSideSamplingMidPoints(1:3,SubP,SubQ,SurfSideID)

        PartState(1:3,PartID)   = NewPos(1:3)
        PartSpecies(PartID)     = iSpec
        LastPartPos(1:3,PartID) = NewPos(1:3)
        IF(CollisMode.GT.1) CALL DSMC_SetInternalEnr(iSpec, BoundID, PartID, 3)
        PDM%ParticleInside(PartID)   = .TRUE.
        PDM%dtFracPush(PartID)       = .TRUE.
        PDM%IsNewPart(PartID)        = .TRUE.
        PEM%GlobalElemID(PartID)     = globElemId
        PEM%LastGlobalElemID(PartID) = globElemId
        iPartTotal    = iPartTotal + 1
        NbrOfParticle = NbrOfParticle + 1
        IF(usevMPF) PartMPF(PartID) = SurfElemMPF
        ! The velocity has to use the triangle the position was accepted in. It is set before the sampling and the
        ! particle balance below, both of which evaluate PartState(4:6) and the particle weight.
        CALL SetChemFluxVelocities(PartID,iSpec,iSF,iSampleAcc,jSampleAcc,BCSideID)
        ! Sampling of the newly created particle
        IF(DSMC%CalcSurfaceVal.AND.(SamplingActive.OR.WriteMacroSurfaceValues)) THEN
          ! CalcWallSample takes the sub-surface indices from TrackInfo, which is only set during tracking
          TrackInfo%p = SubP
          TrackInfo%q = SubQ
          CALL CalcWallSample(PartID,SurfSideID,'new',PartPosImpact_opt=NewPos)
        END IF
        IF(CalcPartBalance) THEN
          ! Compute number of input particles and energy
          nPartIn(iSpec)    = nPartIn(iSpec) + 1
          PartEkinIn(iSpec) = PartEkinIn(iSpec)+CalcEkinPart(PartID)
        END IF ! CalcPartBalance
      END DO ! iPart

      PartsEmitted = PartsEmitted + PartInsSide
#if USE_LOADBALANCE
      !used for calculating LoadBalance of tCurrent(LB_SURFFLUX)
      nSurfacefluxPerElem(ElemID) = nSurfacefluxPerElem(ElemID) + PartInsSide
#endif /*USE_LOADBALANCE*/
    END DO ! SubP
  END DO ! SubQ
#if USE_LOADBALANCE
  CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/

  IF (NbrOfParticle.NE.iPartTotal) CALL abort(__STAMP__, 'ERROR in ParticleSurfChemFlux: NbrOfParticle.NE.iPartTotal')
  IF (iPartTotal.GT.0) THEN
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
    PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))
  END IF
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/

  IF (NbrOfParticle.NE.PartsEmitted) THEN
    ! should be equal for including the following lines in tSurfaceFlux
    CALL abort(__STAMP__,'ERROR in ParticleSurfChemFlux: NbrOfParticle.NE.PartsEmitted')
  END IF
END DO ! iSpec

END SUBROUTINE ParticleSurfChemFlux

!===================================================================================================================================
!> (Instantaneous) Diffusion of particles along the surface corresponding to an averaging over the surface elements
!===================================================================================================================================
SUBROUTINE ParticleSurfDiffusion()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode, SurfSide2GlobalSide, PartBound, nPartBound, nSurfSample
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, ChemWallProp
#if USE_MPI
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
#else
USE MOD_Particle_Boundary_Vars  ,ONLY: nGlobalSurfSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                     :: iSurfSide, GlobalSideID, iPartBound, iBin, nSurfSideLoc
INTEGER                     :: SubP, SubQ
REAL                        :: CovSum(1:nSpecies,1:nPartBound)
REAL                        :: nSidesBin(1:nPartBound)   !< REAL so that it can be reduced together with CovSum
!===================================================================================================================================
IF(.NOT.SurfTotalSideOnNode) RETURN
IF(.NOT.(SurfChem%TotDiffusion.OR.SurfChem%Diffusion)) RETURN

#if USE_MPI
! Only the compute-node roots take part: they own the shared array and they are the ranks that form
! MPI_COMM_LEADERS_SURF (same pattern as ExchangeChemSurfData)
IF (myComputeNodeRank.EQ.0) THEN
  nSurfSideLoc = nComputeNodeSurfSides
#else
  nSurfSideLoc = nGlobalSurfSides
#endif /*USE_MPI*/

  DO SubQ = 1, nSurfSample
    DO SubP = 1, nSurfSample

      CovSum    = 0.
      nSidesBin = 0.

      DO iSurfSide = 1, nSurfSideLoc
        GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
        iPartBound   = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))
        IF (.NOT.SurfChem%BoundIsChemSurf(iPartBound)) CYCLE

        ! Average over all catalytic boundaries -> single bin, otherwise one bin per boundary
        iBin = MERGE(1, iPartBound, SurfChem%TotDiffusion)

        CovSum(1:nSpecies,iBin) = CovSum(1:nSpecies,iBin) + ChemWallProp(1:nSpecies,SubP,SubQ,iSurfSide)
        nSidesBin(iBin)         = nSidesBin(iBin) + 1.
      END DO

#if USE_MPI
      IF (MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,CovSum   ,nSpecies*nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,IERROR)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,nSidesBin,         nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SURF,IERROR)
      END IF
#endif /*USE_MPI*/

      ! --- Mean value per bin
      DO iBin = 1, nPartBound
        IF (nSidesBin(iBin).LE.0.) CYCLE
        CovSum(1:nSpecies,iBin) = CovSum(1:nSpecies,iBin)/nSidesBin(iBin)
      END DO

      ! --- 2nd pass: redistribute the coverage equally over the sides of the bin
      DO iSurfSide = 1, nSurfSideLoc
        GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
        iPartBound   = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))
        IF (.NOT.SurfChem%BoundIsChemSurf(iPartBound)) CYCLE

        iBin = MERGE(1, iPartBound, SurfChem%TotDiffusion)
        IF (nSidesBin(iBin).LE.0.) CYCLE

        ChemWallProp(1:nSpecies,SubP,SubQ,iSurfSide) = CovSum(1:nSpecies,iBin)
      END DO

    END DO ! SubP
  END DO ! SubQ

#if USE_MPI
END IF ! myComputeNodeRank.EQ.0

CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE ParticleSurfDiffusion

!===================================================================================================================================
!>
!===================================================================================================================================
FUNCTION CalcPartPosAxisym(minPos,RVec)
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)            :: minPos(2), RVec(2)
REAL                        :: CalcPartPosAxisym(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal1, Particle_pos(3)
!===================================================================================================================================
IF ((.NOT.(ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))))) THEN
  CALL RANDOM_NUMBER(RandVal1)
  Particle_pos(2) = minPos(2) + RandVal1 * RVec(2)
  ! x-position depending on the y-location
  Particle_pos(1) = minPos(1) + (Particle_pos(2)-minPos(2)) * RVec(1) / RVec(2)
  Particle_pos(3) = 0.
ELSE
  CALL RANDOM_NUMBER(RandVal1)
  IF (ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
    ! y_min = y_max, faces parallel to x-direction, constant distribution
    Particle_pos(1:2) = minPos(1:2) + RVec(1:2) * RandVal1
  ELSE
  ! No LinearWeighting, regular linear distribution of particle positions
    Particle_pos(1:2) = minPos(1:2) + RVec(1:2) &
        * ( SQRT(RandVal1*((minPos(2) + RVec(2))**2-minPos(2)**2)+minPos(2)**2) - minPos(2) ) / (RVec(2))
  END IF
  Particle_pos(3) = 0.
END IF

CalcPartPosAxisym = Particle_pos

END FUNCTION CalcPartPosAxisym


!===================================================================================================================================
!> Chemistry SurfaceFlux: Simplified version of SetSurfacefluxVelocities under the assumption of velocity magnitude = 0
!===================================================================================================================================
SUBROUTINE SetChemFluxVelocities(PartID,iSpec,iSF,iSample,jSample,BCSideID)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars              ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound
USE MOD_Particle_Surfaces_Vars    ,ONLY: SurfMeshSubSideData
USE MOD_Part_Tools                ,ONLY: InRotRefFrameCheck
USE MOD_SurfaceModel_Vars         ,ONLY: SurfChem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: PartID
INTEGER,INTENT(IN)               :: iSpec,iSF,iSample,jSample,BCSideID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: Vec3D(3), vec_nIn(1:3), vec_t1(1:3), vec_t2(1:3)
REAL                             :: RandVal1,RandVal2(2),Velo1,Velo2,Velosq,Temp
!===================================================================================================================================

Temp = PartBound%WallTemp(SurfChem%SurfacefluxBC(iSF))

vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
CALL RANDOM_NUMBER(RandVal1)
!-- 1.: sample normal directions and build complete velo-vector
Vec3D(1:3) = vec_nIn(1:3) * SQRT(2.*BoltzmannConst*Temp/Species(iSpec)%MassIC)*SQRT(-LOG(RandVal1))
Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal2)
  Velo1 = 2.*RandVal2(1) - 1.
  Velo2 = 2.*RandVal2(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
Velo1 = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
Velo2 = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
Vec3D(1:3) = Vec3D(1:3) + vec_t1(1:3) * (Velo1*SQRT(BoltzmannConst*Temp/Species(iSpec)%MassIC))
Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3) * (Velo2*SQRT(BoltzmannConst*Temp/Species(iSpec)%MassIC))
PartState(4:6,PartID) = Vec3D(1:3)

IF(UseRotRefFrame) THEN
  ! Detect if particle is within a RotRefDomain
  InRotRefFrame(PartID) = InRotRefFrameCheck(PartID)
  ! Initialize velocity in the rotational frame of reference
  IF(InRotRefFrame(PartID)) THEN
    PartVeloRotRef(1:3,PartID) = PartState(4:6,PartID) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,PartID))
  END IF
END IF

END SUBROUTINE SetChemFluxVelocities


END MODULE MOD_Particle_SurfChemFlux
