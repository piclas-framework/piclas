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
!> Particle insertion by pure surface reactions
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
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, RadialWeighting
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemMidPoint_Shared
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Part_Tools              ,ONLY: GetParticleWeight
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass, SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Boundary_Vars
USE MOD_SurfaceModel_Vars      ,ONLY: ChemWallProp_Shared_Win,SurfChemReac, ChemWallProp, ChemDesorpWall
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
USE MOD_Particle_SurfFlux       ,ONLY: SetSurfacefluxVelocities, CalcPartPosTriaSurface, DefineSideDirectVec2D
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr
#if USE_MPI
USE MOD_MPI_Shared_vars         ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
#endif
!#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
!#endif /*IMPA*/
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
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, SideID, NbrOfParticle, ParticleIndexNbr
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSubSide, iPart, iPartTotal
INTEGER                     :: PartsEmitted, Node1, Node2, globElemId, CNElemID
REAL                        :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3), RVec(2), minPos(2)
REAL                        :: ReacHeat, DesHeat, RanNum, Area, AdsDens
REAL                        :: DesCount, SurfElemMPF
REAL                        :: nu, E_act, Coverage, Rate, DissOrder, AdCount
REAL                        :: BetaCoeff
REAL                        :: WallTemp
REAL                        :: SurfMol, SurfMolAbs
REAL                        :: MPF
INTEGER                     :: iReac, ReactantCount, BoundID, nSF
INTEGER                     :: iVal, iReactant, iValReac, SurfSideID, iBias
INTEGER                     :: SubP, SubQ
INTEGER                     :: SurfReacBias(SurfChemReac%NumOfReact)
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
! 1.) Determine the surface parameters
nSF = SurfChemReac%CatBoundNum
SubP = TrackInfo%p
SubQ = TrackInfo%q

DO iSF = 1, nSF
  BoundID = SurfChemReac%Surfaceflux(iSF)%BC
  IF(ANY(SurfChemReac%PSMap(BoundID)%PureSurfReac)) THEN

    DO iSide = 1, BCdata_auxSF(BoundID)%SideNumber
      BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      globElemId = ElemID + offSetElem
      CNElemID = GetCNElemID(globElemId)
      SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      IF (RadialWeighting%DoRadialWeighting) THEN
        SurfElemMPF = CalcRadWeightMPF(ElemMidPoint_Shared(2,CNElemID), iSpec, ElemID)
      ELSE
        SurfElemMPF = Species(1)%MacroParticleFactor
      END IF

      IF (SurfSideID.LT.1) CALL abort(__STAMP__,'Chemical Surface Flux is not allowed on non-sampling sides!')

      WallTemp = PartBound%WallTemp(BoundID) ! Boundary temperature

      Area = SurfSideArea_Shared(SubP, SubQ,SurfSideID)
      IF(PartBound%LatticeVec(BoundID).GT.0.) THEN
      ! Number of surface molecules in dependence of the occupancy of the unit cell
        SurfMol = PartBound%MolPerUnitCell(BoundID) * Area &
                  /(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
        SurfMolAbs = PartBound%MolPerUnitCell(BoundID) /(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
      ELSE
      ! Alternative calculation by the average number of surface molecules per area for a monolayer
        SurfMol = 10.**19 * Area
        SurfMolAbs = 10.**19
      END IF

      ! Randomize the order in which the reactions are called to remove biases
      SurfReacBias = RemoveBias(SurfChemReac%NumOfReact)

      ! Loop over the different types of pure surface reactions
      DO iBias = 1, SurfChemReac%NumOfReact
        iReac = SurfReacBias(iBias)
        IF (SurfChemReac%PSMap(BoundID)%PureSurfReac(iReac)) THEN

          ! 2.) Calculate the number of newly created products and update the surface properties
          SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))

          ! 2a) Langmuir-Hinshelwood reaction with instantaneous desorption (Arrhenius model)
          CASE('LHD')
            AdsDens = 1.
            Coverage = 1.
            ! Product of the reactant coverage values
            DO iVal=1,SIZE(SurfChemReac%Reactants(iReac,:))
              IF(SurfChemReac%Reactants(iReac,iVal).GT.0) THEN
                iSpec = SurfChemReac%Reactants(iReac,iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac%Reactants(iReac,:).EQ.iSpec)
                IF(iSpec.NE.SurfChemReac%SurfSpecies) THEN
                  Coverage = Coverage * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID)
                  ! Get the particle density per surface area from the coverage
                  AdsDens = AdsDens * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) * SurfMolAbs
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            BetaCoeff = SurfChemReac%HeatAccommodation(iReac)
            ReacHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

            nu = SurfChemReac%Prefactor(iReac)
            E_act =  SurfChemReac%ArrheniusEnergy(iReac)

            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * AdsDens * exp(-E_act/WallTemp) ! Energy in K

            DO iVal=1,SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)

                ! Randomize the reaction process
                CALL RANDOM_NUMBER(RanNum)

                ! Number of products to be inserted into the gas phase
                ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) =  Rate * dt * Area * (-LOG(RanNum)) / ReactantCount + &
                                                                    ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID)

                DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                  IF (SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                    iReactant = SurfChemReac%Reactants(iReac,iValReac)
                    IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
                      Coverage = ChemWallProp(iReactant,1,SubP,SubQ,SurfSideID)
                    ELSE
                      Coverage = 1.
                    END IF

                    AdCount = Coverage * SurfMol

                    ! Check if enough adsorbate reactants are available
                    IF(ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) .GT. AdCount/ReactantCount) THEN
                      ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) = AdCount/ReactantCount
                    END IF

                  END IF
                END DO

                ! Update the surface coverage values and the heat flux
                IF(INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8).GE.1) THEN
                  ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) &
                                                 + INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8) * ReacHeat * BetaCoeff
                  DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                    IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                      iReactant = SurfChemReac%Reactants(iReac,iValReac)
                      IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
                        ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) &
                                                                  - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol
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
            DO iVal=1,SIZE(SurfChemReac%Reactants(iReac,:))
              IF(SurfChemReac%Reactants(iReac,iVal).GT.0) THEN
                iSpec = SurfChemReac%Reactants(iReac,iVal)
                ! Test for multiples of the same reactant
                ReactantCount = COUNT(SurfChemReac%Reactants(iReac,:).EQ.iSpec)
                IF(iSpec.NE.SurfChemReac%SurfSpecies) THEN
                  Coverage = Coverage * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID)
                  ! Calculate the absolute particle density per unit element
                  AdsDens = AdsDens * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) * SurfMolAbs
                END IF
              END IF
            END DO

            ! Determine the reaction energy in dependence of the surface coverage [J]
            ! Complete accomodation due to the intermediate desorption step
            ReacHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

            nu = SurfChemReac%Prefactor(iReac)
            E_act =  SurfChemReac%ArrheniusEnergy(iReac)
            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * AdsDens * exp(-E_act/WallTemp)! Energy in K

            ! Randomize the reaction process
            CALL RANDOM_NUMBER(RanNum)

            ! Reaction product number
            DesCount =  Rate * dt * Area * (-LOG(RanNum))/(SurfMol*ReactantCount)

            DO iVal=1,SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)

                DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                  IF (SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                    iReactant = SurfChemReac%Reactants(iReac,iValReac)
                    IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
                      Coverage = ChemWallProp(iReactant,1,SubP,SubQ,SurfSideID)
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
                ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID) + DesCount
                ! Test for the maximum of the product coverage
                IF(ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID).GT.PartBound%MaxCoverage(BoundID, iSpec)) THEN
                  ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID) = PartBound%MaxCoverage(BoundID, iSpec)
                END IF
                ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) + DesCount*ReacHeat*SurfMol
                DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                  IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                    iReactant = SurfChemReac%Reactants(iReac,iValReac)
                    IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
                      ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1,SubP,SubQ,SurfSideID) - DesCount
                    END IF
                  END IF
                END DO ! iValReac
              END IF ! iVal in Products
            END DO ! iVal


          ! c) Thermal desorption (Polanyi-Wigner equation)
          CASE('D')
            DO iVal=1, SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)

                ! Number of adsorbed particles on the subside
                IF(ANY(SurfChemReac%Reactants(iReac,:).NE.0)) THEN
                  DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                    IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                      iReactant = SurfChemReac%Reactants(iReac,iValReac)
                      IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
                        Coverage = ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID)
                      ELSE
                        Coverage = 1.
                      END IF
                      AdCount = Coverage * SurfMol
                    END IF
                  END DO
                ELSE
                  Coverage = ChemWallProp(iSpec,1,SubP, SubQ, SurfSideID)
                  AdCount = Coverage * SurfMol
                END IF
                ! Absolute particle density of the surface element
                AdsDens = Coverage * SurfMolAbs

                ! Calculate the desorption energy in dependence of the coverage [J]
                DesHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

                ! Define the variables
                DissOrder = SurfChemReac%DissOrder(iReac)
                nu = SurfChemReac%Prefactor(iReac)

                ! Calculate the desorption prefactor in dependence of coverage and temperature of the boundary
                IF(nu.EQ.0.) THEN
                  nu = 10.**(SurfChemReac%C_a(iReac) + SurfChemReac%C_b(iReac) * Coverage)
                  IF (DissOrder.EQ.2) THEN
                    ! Convert the prefactor to coverage values for the associative desorption
                    nu = 10.**(SurfChemReac%C_a(iReac) + SurfChemReac%C_b(iReac) * Coverage) *10.**(15)
                  END IF
                END IF

                E_act = SurfChemReac%E_initial(iReac) + Coverage * SurfChemReac%W_interact(iReac)
                Rate = nu * AdsDens**DissOrder * exp(-E_act/WallTemp)  ! Energy in K
                ! Randomize the desorption process
                CALL RANDOM_NUMBER(RanNum)

                ! Determine the desorption probability
                ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) = (Rate * dt * Area*(-LOG(RanNum)))/DissOrder + &
                                                                  ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)

                IF(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID).GE.(AdCount/DissOrder)) THEN
                  ! Upper bound for the desorption number
                  ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) = AdCount/DissOrder
                END IF

                ! Update the adsorbtion and desorption count together with the heat flux
                IF(INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/Species(iSpec)%MacroParticleFactor).GE.1) THEN
                  ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) &
                                                              - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8) * DesHeat
                  IF(ANY(SurfChemReac%Reactants(iReac,:).NE.0)) THEN
                    DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                      IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                        iReactant = SurfChemReac%Reactants(iReac,iValReac)
                        IF (iReactant.NE.SurfChemReac%SurfSpecies) THEN
                          ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID) &
                                                        - DissOrder*INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol
                        END IF
                      END IF
                    END DO
                  ELSE
                    ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) = ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) &
                                                        - DissOrder*INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol
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

      ! 3.) Insert the product species into the gas phase
      DO iSpec = 1, nSpecies

#if USE_LOADBALANCE
        CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

        IF (INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/SurfElemMPF).GE.1) THEN

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

            PartInsSubSide = INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/SurfElemMPF)

            ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) = ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) &
                                                            - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)
            NbrOfParticle = NbrOfParticle + PartInsSubSide

            !-- Fill Particle Informations (PartState, Partelem, etc.)
            ParticleIndexNbr = 1
            DO iPart=1,PartInsSubSide
              IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
                ParticleIndexNbr = PDM%nextFreePosition(iPartTotal + 1 + PDM%CurrentNextFreePosition)
              END IF
              IF (ParticleIndexNbr .ne. 0) THEN
                IF(Symmetry%Axisymmetric) THEN
                  PartState(1:3,ParticleIndexNbr) = CalcPartPosAxisym(minPos, RVec)
                ELSE
                  PartState(1:3,ParticleIndexNbr) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
                END IF
                LastPartPos(1:3,ParticleIndexNbr)=PartState(1:3,ParticleIndexNbr)
                PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
                PDM%dtFracPush(ParticleIndexNbr) = .TRUE.
                PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
                PEM%GlobalElemID(ParticleIndexNbr) = globElemId
                PEM%LastGlobalElemID(ParticleIndexNbr) = globElemId
                iPartTotal = iPartTotal + 1
                IF(usevMPF)THEN
                  PartMPF(ParticleIndexNbr) = SurfElemMPF
                END IF ! usevMPF
              ELSE
                CALL abort(__STAMP__,'ERROR in ParticleSurfChemFlux: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
              END IF
            END DO

            CALL SetSurfacefluxVelocities(2,iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,NbrOfParticle,PartInsSubSide)
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

        ! Set the particle properties
        CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)

        IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(iSpec,-1,NbrOfParticle)

        IF (useDSMC.AND.(CollisMode.GT.1)) THEN
          DO iPart = 1,NbrOfParticle
            PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
            IF (PositionNbr .NE. 0) CALL DSMC_SetInternalEnr(iSpec, BoundID, PositionNbr, 3, iReac)
          END DO
        END IF

        IF(CalcPartBalance) THEN
        ! Compute number of input particles and energy
          nPartIn(iSpec)=nPartIn(iSpec) + NBrofParticle

          DO iPart=1,NbrOfparticle
            PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
            IF (PositionNbr .ne. 0) PartEkinIn(PartSpecies(PositionNbr))= &
                                    PartEkinIn(PartSpecies(PositionNbr))+CalcEkinPart(PositionNbr)
          END DO ! iPart
        END IF ! CalcPartBalance

        PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
        PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
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
!> Bias treatment for multiple reactions on the same surface element
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
!> (Instantaneous) Diffusion of particles along the surface
!===================================================================================================================================
SUBROUTINE ParticleSurfDiffusion()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, nComputeNodeProcessors
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Boundary_Vars
USE MOD_SurfaceModel
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars      ,ONLY: ChemWallProp_Shared_Win,SurfChemReac, ChemWallProp
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
!#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackInfo
!#endif /*IMPA*/
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
INTEGER                     :: CatBoundNum
INTEGER                     :: SurfSideID
INTEGER                     :: SubP, SubQ
REAL                        :: Coverage_Sum
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
CatBoundNum = SurfChemReac%CatBoundNum

SubP = TrackInfo%p
SubQ = TrackInfo%q

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nSurfTotalSides
#endif /*USE_MPI*/

SideNumber = lastSide - firstSide + 1

! Average/diffusion over all catalytic boundaries
IF(SurfChemReac%TotDiffusion) THEN
  DO iSpec = 1, nSpecies
    ChemWallProp(iSpec,1,SubP,SubQ,:) = SUM(ChemWallProp(iSpec,1,SubP,SubQ,:))/SideNumber
  END DO

! Diffusion over a single reactive boundary
ELSE IF(SurfChemReac%Diffusion) THEN
  DO iSF = 1, CatBoundNum
    BoundID = SurfChemReac%Surfaceflux(iSF)%BC
    SideNumber = BCdata_auxSF(BoundID)%SideNumber

    ! Determine the sum of the coverage on all indivual subsides
    DO iSpec = 1, nSpecies
      Coverage_Sum = 0.0
      DO iSide = 1, SideNumber
        BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        globElemId = ElemID + offSetElem
        SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
        SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

        Coverage_Sum = Coverage_Sum + ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID)

      END DO

      ! Redistribute the coverage equally over all subsides
      DO iSide = 1, SideNumber
        BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        globElemId = ElemID + offSetElem
        SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
        SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

        ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) = Coverage_Sum/SideNumber

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

END MODULE MOD_Particle_SurfChemFlux