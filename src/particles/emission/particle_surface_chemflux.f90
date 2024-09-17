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
PUBLIC :: ParticleSurfChemFlux, ParticleSurfDiffusion, RemoveBias
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
SUBROUTINE ParticleSurfChemFlux()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Part_Tools              ,ONLY: CalcRadWeightMPF, GetNextFreePosition
USE MOD_DSMC_Vars               ,ONLY: CollisMode, RadialWeighting
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemMidPoint_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, GlobalSide2SurfSide, SurfSideArea, SurfTotalSideOnNode
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, SurfChemReac, ChemWallProp, ChemDesorpWall
USE MOD_Particle_SurfFlux       ,ONLY: CalcPartPosTriaSurface, DefineSideDirectVec2D
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr
#if USE_MPI
USE MOD_MPI_Shared_vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: nSurfacefluxPerElem
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
INTEGER                     :: iSpec , iSF, iSide, SideID, NbrOfParticle, PartID
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSubSide, iPart, iPartTotal
INTEGER                     :: PartsEmitted, Node1, Node2, globElemId, CNElemID
REAL                        :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3), RVec(2), minPos(2)
REAL                        :: ReacHeat, DesHeat, RanNum, Area, AdsDens
REAL                        :: DesCount, SurfElemMPF
REAL                        :: nu, E_act, Coverage, Rate, DissOrder, AdCount
REAL                        :: BetaCoeff
REAL                        :: WallTemp
REAL                        :: SurfMol, SurfMolDens
INTEGER                     :: iReac, ReactantCount, BoundID
INTEGER                     :: iVal, iReactant, iValReac, SurfSideID, iBias
INTEGER                     :: SubP, SubQ
INTEGER                     :: SurfReacBias(SurfChem%NumOfReact)
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

IF(.NOT.SurfTotalSideOnNode) RETURN

! 1.) Determine the surface parameters
! TODO: TrackInfo is only available during tracking
! SubP = TrackInfo%p
! SubQ = TrackInfo%q
SubP = 1
SubQ = 1

DO iSF = 1, SurfChem%CatBoundNum
  BoundID = SurfChem%SurfacefluxBC(iSF)
  IF(ANY(SurfChem%PSMap(BoundID)%PureSurfReac)) THEN

    DO iSide = 1, BCdata_auxSF(BoundID)%SideNumber
      BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      globElemId = ElemID + offSetElem
      CNElemID = GetCNElemID(globElemId)
      SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      IF (SurfSideID.LT.1) CALL abort(__STAMP__,'Chemical Surface Flux is not allowed on non-sampling sides!')

      WallTemp = PartBound%WallTemp(BoundID) ! Boundary temperature

      Area = SurfSideArea(SubP, SubQ,SurfSideID)
      IF(PartBound%LatticeVec(BoundID).GT.0.) THEN
      ! Absolute number of surface molecules in dependence of the occupancy of the unit cell
        SurfMol = PartBound%MolPerUnitCell(BoundID) * Area &
                  /(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
      ! Number of surface molecules per m^3 surface
        SurfMolDens = PartBound%MolPerUnitCell(BoundID) /(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
      ELSE
      ! Alternative calculation by the average number of surface molecules per area for a monolayer
        SurfMol = 10.**19 * Area
        ! Number of surface molecules per m^3 surface
        SurfMolDens = 10.**19
      END IF

      ! Randomize the order in which the reactions are called to remove biases
      SurfReacBias = RemoveBias(SurfChem%NumOfReact)

      ! Loop over the different types of pure surface reactions
      DO iBias = 1, SurfChem%NumOfReact
        iReac = SurfReacBias(iBias)
        IF (SurfChem%PSMap(BoundID)%PureSurfReac(iReac)) THEN

          ! 2.) Calculate the number of newly created products and update the surface properties
          SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

          ! 2a) Langmuir-Hinshelwood reaction with instantaneous desorption (Arrhenius model)
          CASE('LHD')
            AdsDens = 1.
            Coverage = 1.
            ! Product of the reactant coverage values
            DO iVal=1,SIZE(SurfChemReac(iReac)%Reactants(:))
              IF(SurfChemReac(iReac)%Reactants(iVal).GT.0) THEN
                iSpec = SurfChemReac(iReac)%Reactants(iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iSpec)
                IF(iSpec.NE.SurfChem%SurfSpecies) THEN ! Coverage set to 1 for surface species
                  Coverage = Coverage * ChemWallProp(iSpec,SubP,SubQ,SurfSideID)
                  ! Particle density of the adsorbate on the surface
                  AdsDens = AdsDens * ChemWallProp(iSpec,SubP,SubQ,SurfSideID) * SurfMolDens
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            BetaCoeff = SurfChemReac(iReac)%HeatAccommodation !! Incomplete energy accomodation
            ReacHeat = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst

            nu = SurfChemReac(iReac)%Prefactor
            E_act =  SurfChemReac(iReac)%ArrheniusEnergy

            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * AdsDens * exp(-E_act/WallTemp) ! Energy in K

            DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
              IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                iSpec = SurfChemReac(iReac)%Products(iVal)

                ! Randomize the reaction process
                CALL RANDOM_NUMBER(RanNum)

                ! Number of products to be inserted into the gas phase
                ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) =  Rate * dt * Area * (-LOG(RanNum)) / ReactantCount + &
                                                                    ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID)

                DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                  IF (SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                    iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                    IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                      Coverage = ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                    ELSE
                      Coverage = 1.
                    END IF

                    ! Number of available adsorbates
                    AdCount = Coverage * SurfMol

                    ! Check if enough adsorbate reactants are available
                    IF(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) .GT. AdCount/ReactantCount) THEN
                      ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) = AdCount/ReactantCount
                    END IF

                  END IF
                END DO

                ! Update the surface coverage values and the heat flux
                IF(INT(ChemDesorpWall(iSpec,SubP, SubQ, SurfSideID),8).GE.1) THEN
                  ChemWallProp(nSpecies+1, SubP, SubQ, SurfSideID) = ChemWallProp(nSpecies+1, SubP, SubQ, SurfSideID) &
                                                 + INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8) * ReacHeat * BetaCoeff
                  DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                    IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                      iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                      IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                        ChemWallProp(iReactant, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant, SubP, SubQ, SurfSideID) &
                                                                  - INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8)/SurfMol
                      END IF
                    END IF
                  END DO ! iValReac
                END IF !ChemDesorpWall.GE.1
              END IF ! iVal in Products
            END DO ! iVal

          ! b) Langmuir-Hinshelwood reaction (Arrhenius model)
          CASE('LH')
            Coverage = 1.
            AdsDens = 1.
            ! Product of the reactant coverage values
            DO iVal=1,SIZE(SurfChemReac(iReac)%Reactants(:))
              IF(SurfChemReac(iReac)%Reactants(iVal).GT.0) THEN
                iSpec = SurfChemReac(iReac)%Reactants(iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac(iReac)%Reactants(:).EQ.iSpec)
                IF(iSpec.NE.SurfChem%SurfSpecies) THEN ! Coverage set to 1 for surface species
                  Coverage = Coverage * ChemWallProp(iSpec,SubP,SubQ,SurfSideID)
                  ! Calculate the adsorbate particle density per unit element
                  AdsDens = AdsDens * ChemWallProp(iSpec,SubP,SubQ,SurfSideID) * SurfMolDens
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            ! Complete accommodation due to the intermediate desorption step
            ReacHeat = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst

            nu = SurfChemReac(iReac)%Prefactor
            E_act =  SurfChemReac(iReac)%ArrheniusEnergy
            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * AdsDens * exp(-E_act/WallTemp)! Energy in K

            ! Randomize the reaction process
            CALL RANDOM_NUMBER(RanNum)

            ! Reaction product number
            DesCount =  Rate * dt * Area * (-LOG(RanNum))/(SurfMol*ReactantCount)

            DO iVal=1,SIZE(SurfChemReac(iReac)%Products(:))
              IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                iSpec = SurfChemReac(iReac)%Products(iVal)

                DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                  IF (SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                    iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                    IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                      Coverage = ChemWallProp(iReactant,SubP,SubQ,SurfSideID)
                    ELSE
                      Coverage = 1.
                    END IF

                    ! Check if enough adsorbate reactants are available
                    IF(DesCount .GT. Coverage/ReactantCount) THEN
                      DesCount = Coverage/ReactantCount
                    END IF

                  END IF
                END DO

                ! Consider reactions between the same species
                DesCount = DesCount / ReactantCount

                ! Update the surface coverage values and the heat flux
                ChemWallProp(iSpec, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec, SubP, SubQ, SurfSideID) + DesCount
                ! Test for the maximum of the product coverage
                IF(ChemWallProp(iSpec, SubP, SubQ, SurfSideID).GT.PartBound%MaxCoverage(BoundID, iSpec)) THEN
                  ChemWallProp(iSpec, SubP, SubQ, SurfSideID) = PartBound%MaxCoverage(BoundID, iSpec)
                END IF
                ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) = ChemWallProp(nSpecies+1,SubP,SubQ,SurfSideID) + DesCount*ReacHeat*SurfMol
                DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                  IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                    iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                    IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                      ChemWallProp(iReactant, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,SubP,SubQ,SurfSideID) - DesCount
                    END IF
                  END IF
                END DO ! iValReac
              END IF ! iVal in Products
            END DO ! iVal


          ! c) Thermal desorption (Polanyi-Wigner equation)
          CASE('D')
            DO iVal=1, SIZE(SurfChemReac(iReac)%Products(:))
              IF (SurfChemReac(iReac)%Products(iVal).NE.0) THEN
                iSpec = SurfChemReac(iReac)%Products(iVal)

                ! Number of adsorbed particles on the subside
                IF(ANY(SurfChemReac(iReac)%Reactants(:).NE.0)) THEN
                  DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                    IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                      iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                      IF(iReactant.NE.SurfChem%SurfSpecies) THEN
                        Coverage = ChemWallProp(iReactant,SubP, SubQ, SurfSideID)
                      ELSE
                        ! Coverage set to 1 for surface species
                        Coverage = 1.
                      END IF
                      AdCount = Coverage * SurfMol
                    END IF
                  END DO
                ELSE
                  Coverage = ChemWallProp(iSpec,SubP, SubQ, SurfSideID)
                  AdCount = Coverage * SurfMol
                END IF
                ! Absolute particle density of the surface element
                AdsDens = Coverage * SurfMolDens

                ! Calculate the desorption energy in dependence of the coverage [J]
                DesHeat = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst

                ! Define the variables
                DissOrder = SurfChemReac(iReac)%DissOrder
                nu = SurfChemReac(iReac)%Prefactor

                ! Calculate the desorption prefactor in dependence of coverage and temperature of the boundary
                IF(nu.EQ.0.) THEN
                  nu = 10.**(SurfChemReac(iReac)%C_a + SurfChemReac(iReac)%C_b * Coverage)
                  IF (DissOrder.EQ.2) THEN
                    ! Convert the prefactor to coverage values for the associative desorption
                    nu = 10.**(SurfChemReac(iReac)%C_a + SurfChemReac(iReac)%C_b * Coverage) *10.**(15) !!
                  END IF
                END IF

                E_act = SurfChemReac(iReac)%E_initial + Coverage * SurfChemReac(iReac)%W_interact
                ! Reaction rate according to the Polanyi-Wigner equation
                Rate = nu * AdsDens**DissOrder * exp(-E_act/WallTemp)  ! Energy in K
                ! Randomize the desorption process
                CALL RANDOM_NUMBER(RanNum)

                ! Determine the desorption probability
                ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) = (Rate * dt * Area*(-LOG(RanNum)))/DissOrder + &
                                                                  ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID)

                IF(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID).GE.(AdCount/DissOrder)) THEN
                  ! Upper bound for the desorption number
                  ChemDesorpWall(iSpec,  SubP, SubQ, SurfSideID) = AdCount/DissOrder
                END IF

                ! Update the adsorbtion and desorption count together with the heat flux
                IF(INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID)/Species(iSpec)%MacroParticleFactor).GE.1) THEN
                  ChemWallProp(nSpecies+1, SubP, SubQ, SurfSideID) = ChemWallProp(nSpecies+1, SubP, SubQ, SurfSideID) &
                                                              - INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8) * DesHeat
                  IF(ANY(SurfChemReac(iReac)%Reactants(:).NE.0)) THEN
                    DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
                      IF(SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
                        iReactant = SurfChemReac(iReac)%Reactants(iValReac)
                        IF (iReactant.NE.SurfChem%SurfSpecies) THEN
                          ChemWallProp(iReactant,SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,SubP, SubQ, SurfSideID) &
                                                        - DissOrder*INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8)/SurfMol
                        END IF
                      END IF
                    END DO
                  ELSE
                    ChemWallProp(iSpec,SubP,SubQ,SurfSideID) = ChemWallProp(iSpec,SubP,SubQ,SurfSideID) &
                                                        - DissOrder*INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8)/SurfMol
                  END IF
                END IF !ChemDesorbWall .GE. 1
              END IF ! Products .NE. 1
            END DO !iSpec

          CASE DEFAULT
          END SELECT

        END IF !iReac.EQ.PureSurfReac

      END DO !iBias

      ! Current boundary condition
      PartsEmitted = 0
      NbrOfParticle = 0
      iPartTotal = 0

      ! 3.) Insert the product species into the gas phase by a surface flux from the boundary
      DO iSpec = 1, nSpecies
#if USE_LOADBALANCE
        CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
        IF (RadialWeighting%DoRadialWeighting) THEN
          SurfElemMPF = CalcRadWeightMPF(ElemMidPoint_Shared(2,CNElemID), iSpec, ElemID)
        ELSE
          SurfElemMPF = Species(iSpec)%MacroParticleFactor
        END IF
        IF (INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID)/SurfElemMPF).GE.1) THEN

          ! Define the necessary variables
          xyzNod(1:3) = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%xyzNod(1:3)

          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            Node1 = jSample+1
            Node2 = jSample+2
            Vector1 = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
            Vector2 = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
            midpoint(1:3) = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
            ndist(1:3) = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)

            ! REQUIRED LATER FOR THE POSITION START
            IF(Symmetry%Axisymmetric) CALL DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)

            PartInsSubSide = INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID)/SurfElemMPF)

            ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) = ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID) &
                                                            - INT(ChemDesorpWall(iSpec, SubP, SubQ, SurfSideID),8)
            NbrOfParticle = NbrOfParticle + PartInsSubSide

            !-- Fill Particle Informations (PartState, Partelem, etc.)
            PartID = 1
            DO iPart=1,PartInsSubSide
              IF ((iPart.EQ.1).OR.PDM%ParticleInside(PartID)) PartID = GetNextFreePosition(iPartTotal+1)
              IF(Symmetry%Axisymmetric) THEN
                PartState(1:3,PartID) = CalcPartPosAxisym(minPos, RVec)
              ELSE
                PartState(1:3,PartID) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
              END IF
              PartSpecies(PartID) = iSpec
              LastPartPos(1:3,PartID)=PartState(1:3,PartID)
              IF(CollisMode.GT.1) CALL DSMC_SetInternalEnr(iSpec, BoundID, PartID, 3)
              PDM%ParticleInside(PartID) = .TRUE.
              PDM%dtFracPush(PartID) = .TRUE.
              PDM%IsNewPart(PartID) = .TRUE.
              PEM%GlobalElemID(PartID) = globElemId
              PEM%LastGlobalElemID(PartID) = globElemId
              iPartTotal = iPartTotal + 1
              IF(usevMPF)THEN
                PartMPF(PartID) = SurfElemMPF
              END IF ! usevMPF
              IF(CalcPartBalance) THEN
                ! Compute number of input particles and energy
                  nPartIn(iSpec) = nPartIn(iSpec) + 1
                  PartEkinIn(iSpec) = PartEkinIn(iSpec)+CalcEkinPart(PartID)
              END IF ! CalcPartBalance
              CALL SetChemFluxVelocities(PartID,iSpec,iSF,iSample,jSample,BCSideID)
            END DO
            PartsEmitted = PartsEmitted + PartInsSubSide
#if USE_LOADBALANCE
            !used for calculating LoadBalance of tCurrent(LB_SURFFLUX)
            nSurfacefluxPerElem(ElemID)=nSurfacefluxPerElem(ElemID)+PartInsSubSide
#endif /*USE_LOADBALANCE*/
          END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
#if USE_LOADBALANCE
      CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
        END IF ! iSide
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


    END DO !iSide

  ELSE
    CYCLE
  END IF !ANY PureSurfReac
END DO !iSF

#if USE_MPI
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE ParticleSurfChemFlux

!===================================================================================================================================
!> Bias treatment for multiple reactions on the same surface element by randomization of the reaction order
!===================================================================================================================================
FUNCTION RemoveBias(SurfNumOfReac)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: SurfNumOfReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                     :: RemoveBias(SurfNumOfReac)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i, j, k, m
INTEGER                     :: temp
REAL                        :: RanNum
!===================================================================================================================================

RemoveBias = [(i,i=1,SurfNumOfReac)]
! Shuffle
m = SurfNumOfReac
DO k = 1, 2
  DO i = 1, m
    CALL RANDOM_NUMBER(RanNum)
    j = 1 + FLOOR(m*RanNum)
    temp = RemoveBias(j)
    RemoveBias(j) = RemoveBias(i)
    RemoveBias(i) = temp
  END DO
END DO

END FUNCTION RemoveBias

!===================================================================================================================================
!> (Instantaneous) Diffusion of particles along the surface corresponding to an averaging over the surface elements
!===================================================================================================================================
SUBROUTINE ParticleSurfDiffusion()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode, GlobalSide2SurfSide
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, ChemWallProp
#if USE_MPI
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfTotalSides
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, nComputeNodeProcessors
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
INTEGER                     :: firstSide, lastSide, SideNumber
INTEGER                     :: iSpec , iSF, iSide, BoundID, SideID
INTEGER                     :: BCSideID, ElemID, iLocSide
INTEGER                     :: globElemId
INTEGER                     :: SurfSideID
INTEGER                     :: SubP, SubQ
REAL                        :: Coverage_Sum
!===================================================================================================================================
IF(.NOT.SurfTotalSideOnNode) RETURN

! TODO: TrackInfo is only available during tracking
! SubP = TrackInfo%p
! SubQ = TrackInfo%q
SubP = 1
SubQ = 1

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nGlobalSurfSides
#endif /*USE_MPI*/

SideNumber = lastSide - firstSide + 1

! Average/diffusion over all catalytic boundaries
IF(SurfChem%TotDiffusion) THEN
  DO iSpec = 1, nSpecies
    ChemWallProp(iSpec,SubP,SubQ,:) = SUM(ChemWallProp(iSpec,SubP,SubQ,:))/SideNumber
  END DO

! Diffusion over a single reactive boundary
ELSE IF(SurfChem%Diffusion) THEN
  DO iSF = 1, SurfChem%CatBoundNum
    BoundID = SurfChem%SurfacefluxBC(iSF)
    SideNumber = BCdata_auxSF(BoundID)%SideNumber

    ! Determine the sum of the coverage on all individual subsides
    DO iSpec = 1, nSpecies
      Coverage_Sum = 0.0
      DO iSide = 1, SideNumber
        BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        globElemId = ElemID + offSetElem
        SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
        SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

        Coverage_Sum = Coverage_Sum + ChemWallProp(iSpec,SubP,SubQ,SurfSideID)

      END DO

      ! Redistribute the coverage equally over all subsides
      DO iSide = 1, SideNumber
        BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        globElemId = ElemID + offSetElem
        SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
        SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

        ChemWallProp(iSpec,SubP,SubQ,SurfSideID) = Coverage_Sum/SideNumber

      END DO

    END DO !iSpec
  END DO !iSF
END IF !Diffusion

#if USE_MPI
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
  ! No VarWeighting, regular linear distribution of particle positions
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
  PDM%InRotRefFrame(PartID) = InRotRefFrameCheck(PartID)
  ! Initialize velocity in the rotational frame of reference
  IF(PDM%InRotRefFrame(PartID)) THEN
    PartVeloRotRef(1:3,PartID) = PartState(4:6,PartID) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,PartID))
  END IF
END IF

END SUBROUTINE SetChemFluxVelocities


END MODULE MOD_Particle_SurfChemFlux