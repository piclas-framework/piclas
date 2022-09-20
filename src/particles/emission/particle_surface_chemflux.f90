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
!> Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
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
USE MOD_Part_Tools              ,ONLY: GetParticleWeight
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass, SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Boundary_Vars 
USE MOD_SurfaceModel_Vars      ,ONLY: ChemWallProp_Shared_Win,SurfChemReac, ChemWallProp, ChemDesorpWall, ChemCountReacWall
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
!#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod, TrackInfo
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
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, SideID, NbrOfParticle, ParticleIndexNbr
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSubSide, iPart, iPartTotal
INTEGER                     :: PartsEmitted, Node1, Node2, globElemId
REAL                        :: Particle_pos(3), xyzNod(3), Vector1(3), Vector2(3)
REAL                        :: ndist(3), midpoint(3)
REAL,ALLOCATABLE            :: particle_positions(:)
REAL                        :: ReacHeat, DesHeat
REAL                        :: DesCount
REAL                        :: nu, E_act, Coverage, Rate, DissOrder, AdCount
REAL                        :: BetaCoeff
REAL                        :: WallTemp
REAL                        :: SurfMol
REAL                        :: MPF
INTEGER                     :: SurfNumOfReac, iReac, ReactantCount, BoundID, nSF
INTEGER                     :: iVal, iReactant, iValReac, SurfSideID, iBias
INTEGER                     :: SubP, SubQ
INTEGER, ALLOCATABLE        :: SurfReacBias(:)

#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================  
SurfNumOfReac = SurfChemReac%NumOfReact
nSF = SurfChemReac%CatBoundNum

SubP = TrackInfo%p
SubQ = TrackInfo%q

ALLOCATE(SurfReacBias(SurfNumOfReac))

DO iSF = 1, nSF
  BoundID = SurfChemReac%Surfaceflux(iSF)%BC
  IF(ANY(SurfChemReac%PSMap(BoundID)%PureSurfReac)) THEN

    DO iSide = 1, BCdata_auxSF(BoundID)%SideNumber
      BCSideID=BCdata_auxSF(BoundID)%SideList(iSide)
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      globElemId = ElemID + offSetElem
      SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)

      IF (SurfSideID.LT.1) THEN
        CALL abort(&
        __STAMP__&
        ,'Chemical Surface Flux is not allowed on non-sampling sides!')
      END IF  

      ! Boundary temperature
      WallTemp = PartBound%WallTemp(BoundID)

      ! Calculate the number of molecules on the surface
      IF(PartBound%LatticeVec(BoundID).GT.0.) THEN
      ! Surface molecules in dependence of the occupancy of the unit cell
        SurfMol = PartBound%MolPerUnitCell(BoundID) * SurfSideArea_Shared(SubP, SubQ,SurfSideID) &
                  /(PartBound%LatticeVec(BoundID)*PartBound%LatticeVec(BoundID))
      ELSE
      ! Alternative calculation by average number of molecules per area for a monolayer
        SurfMol = 10.**19 * SurfSideArea_Shared(SubP, SubQ,SurfSideID)
      END IF

      ! Randomize the order in which the reactions are called
      CALL RemoveBias(SurfNumOfReac, SurfReacBias)

      ! Loop over the reaction types 
      DO iBias = 1, SurfNumOfReac
        iReac = SurfReacBias(iBias)
        IF (SurfChemReac%PSMap(BoundID)%PureSurfReac(iReac)) THEN

          SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))

          ! LH-reaction with immediate desorption
          CASE('LHD')
            ! Define the variables
            Coverage = 1.
            ! prefactor for the rate equation including the dependency on all coverage values
            DO iVal=1,SIZE(SurfChemReac%Reactants(iReac,:))

              IF(SurfChemReac%Reactants(iReac,iVal).GT.0) THEN
                iSpec = SurfChemReac%Reactants(iReac,iVal)
                MPF = Species(iSpec)%MacroParticleFactor
                Coverage = Coverage * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID)
              END IF
            END DO

            ! Determine the reaction energy [K]
            BetaCoeff = SurfChemReac%HeatAccomodation(iReac)
            ReacHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

            nu = SurfChemReac%Prefactor(iReac)
            E_act =  SurfChemReac%ArrheniusEnergy(iReac)       

            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * Coverage * exp(-E_act/WallTemp) ! Energy in K

            DO iVal=1,SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)

                ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) =  Rate * dt * SurfMol + ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID)

                DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                  IF (SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                    iReactant = SurfChemReac%Reactants(iReac,iValReac)
                    ! Test for multiples of the same reactant
                    ReactantCount = COUNT(SurfChemReac%Reactants(iReac,:) == iReactant)

                    Coverage = ChemWallProp(iReactant,1,SubP,SubQ,SurfSideID)
                    AdCount = Coverage * SurfMol

                    ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) = AdCount/ReactantCount

                    ! Test if enough reactants exist on the surface
                    IF(ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) .GT. AdCount/ReactantCount) THEN
                      ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) = AdCount/ReactantCount
                    END IF

                  END IF
                END DO

                ! Update the surface coverage values 
                IF(INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/MPF).GE.1) THEN
                  ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) + INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8) * ReacHeat * BetaCoeff
                  DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:)) 
                    IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN 
                      iReactant = SurfChemReac%Reactants(iReac,iValReac)
                      ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol  
                    END IF
                  END DO
                  ! Count the number of surface reactions
                  ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) = ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) + INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/MPF)  
                END IF
              END IF ! iVal in Products 
            END DO ! iVal            

          ! LH-reaction
          CASE('LH')
            ! Define the variables
            Coverage = 1.
            ! prefactor for the rate equation including the dependency on all coverage values
            DO iVal=1,SIZE(SurfChemReac%Reactants(iReac,:))

              IF(SurfChemReac%Reactants(iReac,iVal).GT.0) THEN
                iSpec = SurfChemReac%Reactants(iReac,iVal)
                MPF = Species(iSpec)%MacroParticleFactor
                Coverage = Coverage * ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID)
              END IF
            END DO

            ! Determine the reaction energy [K]
            BetaCoeff = SurfChemReac%HeatAccomodation(iReac)
            ReacHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

            nu = SurfChemReac%Prefactor(iReac)    
            E_act =  SurfChemReac%ArrheniusEnergy(iReac)     
            ! Calculate the rate in dependence of the temperature and coverage
            Rate = nu * Coverage * exp(-E_act/WallTemp) ! Energy in K

            DO iVal=1,SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)

                DesCount =  Rate * dt * SurfMol 

                DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
                  IF (SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                    iReactant = SurfChemReac%Reactants(iReac,iValReac)
                    ! Check for multiples of the same reactant
                    ReactantCount = COUNT(SurfChemReac%Reactants(iReac,:) == iReactant)

                    Coverage = ChemWallProp(iReactant,1,SubP,SubQ,SurfSideID)
                    AdCount = Coverage * SurfMol
                    
                    !Check if enough reactants exist on the surface
                    IF(DesCount .GT. AdCount/ReactantCount) THEN
                      DesCount = AdCount/ReactantCount
                    END IF

                  END IF
                END DO

                ! Update the surface coverage values
                IF(INT(DesCount/MPF).GE.1) THEN
                  ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,1, SubP, SubQ, SurfSideID) + INT(DesCount,8)/SurfMol 
                  ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) + INT(DesCount,8) * ReacHeat * BetaCoeff
                  DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:)) 
                    IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                      iReactant = SurfChemReac%Reactants(iReac,iValReac)
                      ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1, SubP, SubQ, SurfSideID) - INT(DesCount,8)/SurfMol  
                    END IF
                  END DO
                  ! Count the number of surface reactions
                  ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) = ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) + INT(DesCount/MPF)     
                END IF
              END IF ! iVal in Products
            END DO ! iVal    


          ! Desorption
          CASE('D')
            DO iVal=1, SIZE(SurfChemReac%Products(iReac,:))
              IF (SurfChemReac%Products(iReac,iVal).NE.0) THEN
                iSpec = SurfChemReac%Products(iReac,iVal)
                MPF = Species(iSpec)%MacroParticleFactor
                
                IF(ANY(SurfChemReac%Reactants(iReac,:).NE.0)) THEN
                  DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:)) 
                    IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                      iReactant = SurfChemReac%Reactants(iReac,iValReac)
                      Coverage = ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID)
                      AdCount = Coverage * SurfMol
                    END IF
                  END DO
                ELSE 
                  Coverage = ChemWallProp(iSpec,1,SubP, SubQ, SurfSideID)
                  AdCount = Coverage * SurfMol
                END IF        

                ! Calculate the desorption energy in dependence of the coverage [K]
                DesHeat = (SurfChemReac%EReact(iReac) - Coverage*SurfChemReac%EScale(iReac)) * BoltzmannConst

                ! Define the variables
                DissOrder = SurfChemReac%DissOrder(iReac)
                nu = SurfChemReac%Prefactor(iReac)
                E_act = SurfChemReac%ArrheniusEnergy(iReac)        
                Rate = SurfChemReac%Rate(iReac)

                ! Calculate the desorption parameters in dependence of coverage and temperature of the boundary
                IF(nu.EQ.0.) THEN
                  nu = 10.**(SurfChemReac%C_a(iReac) + SurfChemReac%C_b(iReac) * Coverage)     
                  IF (DissOrder.EQ.2) THEN
                    ! Convert the prefactor to coverage values for the associative desorption
                    nu = 10.**(SurfChemReac%C_a(iReac) + SurfChemReac%C_b(iReac) * Coverage) *10.**(15)
                  END IF
                END IF
                  
                E_act = SurfChemReac%E_initial(iReac) + Coverage * SurfChemReac%W_interact(iReac)
                Rate = nu * Coverage**DissOrder * exp(-E_act/WallTemp)  ! Energy in K

                ! Determine the desorption probability
                ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) =  Rate * dt * SurfMol + ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/DissOrder

                IF(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID).GE.(AdCount/DissOrder)) THEN
                  ! Upper bound for the desorption number
                  ChemDesorpWall(iSpec, 1, SubP, SubQ, SurfSideID) = AdCount/DissOrder
                END IF

                ! Update the adsorbtion and desorption count
                IF(INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/MPF).GE.1) THEN
                  ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) = ChemWallProp(iSpec,2, SubP, SubQ, SurfSideID) - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8) * DesHeat             
                  IF(ANY(SurfChemReac%Reactants(iReac,:).NE.0)) THEN
                    DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:)) 
                      IF(SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
                        iReactant = SurfChemReac%Reactants(iReac,iValReac)
                        ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID) = ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID) - DissOrder*INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol
                      END IF
                    END DO
                  ELSE 
                    ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) = ChemWallProp(iSpec,1,SubP,SubQ,SurfSideID) - DissOrder*INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)/SurfMol
                  END IF

                  ! Count the number of surface reactions
                    ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) = ChemCountReacWall(iReac, 1, SubP, SubQ, SurfSideID) + INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/MPF)     
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

      DO iSpec = 1, nSpecies

        ! Loop over the side numbers
        IF (INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/Species(iSpec)%MacroParticleFactor).GE.1) THEN

          ! Random distribution of the particles on the surface
          ! Define the necessary variables
          xyzNod(1:3) = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%xyzNod(1:3)

          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            ! Number of additional particles to be inserted 
            Node1 = jSample+1    
            Node2 = jSample+2     !        
            Vector1 = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
            Vector2 = BCdata_auxSF(BoundID)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
            midpoint(1:3) = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
            ndist(1:3) = BCdata_auxSF(BoundID)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)

            !PartInsSubSide = INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID))
            PartInsSubSide = INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID)/Species(iSpec)%MacroParticleFactor)

            ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) = ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID) - INT(ChemDesorpWall(iSpec,1, SubP, SubQ, SurfSideID),8)
            NbrOfParticle = NbrOfParticle + PartInsSubSide
            ALLOCATE(particle_positions(1:PartInsSubSide*3))

            iPart = 1
            !-- Set Positions
            DO WHILE (iPart .LE. PartInsSubSide)
              Particle_pos(1:3) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
              particle_positions(iPart*3-2) = Particle_pos(1)
              particle_positions(iPart*3-1) = Particle_pos(2)
              particle_positions(iPart*3  ) = Particle_pos(3)
              iPart=iPart+1
            END DO

            !-- Fill Particle Informations (PartState, Partelem, etc.)
            ParticleIndexNbr = 1
            DO iPart=1,PartInsSubSide
              IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
                ParticleIndexNbr = PDM%nextFreePosition(iPartTotal + 1 + PDM%CurrentNextFreePosition)
              END IF
              IF (ParticleIndexNbr .ne. 0) THEN
                PartState(1:3,ParticleIndexNbr) = particle_positions(3*(iPart-1)+1:3*(iPart-1)+3)
                LastPartPos(1:3,ParticleIndexNbr)=PartState(1:3,ParticleIndexNbr)
                PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
                PDM%dtFracPush(ParticleIndexNbr) = .TRUE.
                PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
                PEM%GlobalElemID(ParticleIndexNbr) = globElemId
                PEM%LastGlobalElemID(ParticleIndexNbr) = globElemId !needed when ParticlePush is not executed, e.g. "delay"
                iPartTotal = iPartTotal + 1
                IF (RadialWeighting%DoRadialWeighting) THEN
                  PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(2,ParticleIndexNbr), iSpec,ParticleIndexNbr)
                END IF
              ELSE
                CALL abort(&
            __STAMP__&
            ,'ERROR in ParticleSurfaceflux: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
              END IF
            END DO
            DEALLOCATE(particle_positions)
            
            CALL SetSurfacefluxVelocities(iSpec,iReac,iSF,iSample,jSample,iSide,BCSideID,SideID,NbrOfParticle,PartInsSubSide)
                  
            PartsEmitted = PartsEmitted + PartInsSubSide
          END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
        END IF ! iSide
        IF (NbrOfParticle.NE.iPartTotal) CALL abort(__STAMP__, 'Error 2 in ParticleSurfaceflux!')

        CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)

        IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(iSpec,-1,NbrOfParticle)
        ! define molecule stuff

        IF (useDSMC.AND.(CollisMode.GT.1)) CALL SetInnerEnergies(iSpec, BoundID, NbrOfParticle,iReac)

        IF(CalcPartBalance) THEN
        ! Compute number of input particles and energy
          nPartIn(iSpec)=nPartIn(iSpec) + NBrofParticle

          DO iPart=1,NbrOfparticle
            PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
            IF (PositionNbr .ne. 0) PartEkinIn(PartSpecies(PositionNbr))= &
                                    PartEkinIn(PartSpecies(PositionNbr))+CalcEkinPart(PositionNbr)
          END DO ! iPart
        END IF ! CalcPartBalance
        ! instead of an UpdateNextfreePosition we update the particleVecLength only - enough ?!?
        PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
        PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
        ! Sample Energies on Surfaces when particles are emitted from them
        IF (NbrOfParticle.NE.PartsEmitted) THEN
          ! should be equal for including the following lines in tSurfaceFlux
          CALL abort(&
        __STAMP__&
        ,'ERROR in ParticleSurfaceflux: NbrOfParticle.NE.PartsEmitted')
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
SUBROUTINE RemoveBias(SurfNumOfReac, SurfReacBias)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: SurfNumOfReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, ALLOCATABLE, INTENT(OUT) :: SurfReacBias(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i, j, k, m
INTEGER                     :: temp
REAL                        :: RanNum
!===================================================================================================================================

SurfReacBias = [(i,i=1,SurfNumOfReac)]

m = SurfNumOfReac
DO k = 1, 2
    DO i = 1, m
      CALL RANDOM_NUMBER(RanNum)
      j = 1 + FLOOR(m*RanNum)
      temp = SurfReacBias(j)
      SurfReacBias(j) = SurfReacBias(i)
      SurfReacBias(i) = temp
    END DO
END DO

END SUBROUTINE RemoveBias

!===================================================================================================================================
!> Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
!===================================================================================================================================
SUBROUTINE ParticleSurfDiffusion()
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank, nComputeNodeProcessors
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Surfaces_Vars
USE MOD_Particle_Boundary_Vars 
USE MOD_SurfaceModel
USE MOD_SurfaceModel_Chemistry
USE MOD_SurfaceModel_Vars      ,ONLY: ChemWallProp_Shared_Win,SurfChemReac, ChemWallProp
USE MOD_MPI_Shared_vars        ,ONLY: MPI_COMM_SHARED
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
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
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, BoundID, SideID
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample
INTEGER                     :: Node1, Node2, globElemId
INTEGER                     :: CatBoundNum   
INTEGER                     :: SurfSideID
INTEGER                     :: SubP, SubQ
REAL                        :: Coverage_Sum

#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================  
CatBoundNum = SurfChemReac%CatBoundNum

! SubP = TrackInfo%p
! SubQ = TrackInfo%q
SubP = 1
SubQ = 1

#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nSurfTotalSides
#endif /*USE_MPI*/

SideNumber = lastSide - firstSide + 1

IF(SurfChemReac%TotDiffusion) THEN

  DO iSpec = 1, nSpecies
    ChemWallProp(iSpec,1,SubP,SubQ,:) = SUM(ChemWallProp(iSpec,1,SubP,SubQ,:))/SideNumber
  END DO

ELSE IF(SurfChemReac%Diffusion) THEN
  DO iSF = 1, CatBoundNum
    BoundID = SurfChemReac%Surfaceflux(iSF)%BC
    SideNumber = BCdata_auxSF(BoundID)%SideNumber
    
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
SUBROUTINE SetInnerEnergies(iSpec, BoundID, NbrOfParticle,iReac)
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
USE MOD_Particle_Vars           ,ONLY: PDM
USE MOD_DSMC_PolyAtomicModel ,ONLY: DSMC_SetInternalEnr
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                        :: iSpec, BoundID, NbrOfParticle,iReac
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iPart, PositionNbr
!===================================================================================================================================
iPart = 1

DO WHILE (iPart .le. NbrOfParticle)
  PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)

  IF (PositionNbr .ne. 0) THEN
    CALL DSMC_SetInternalEnr(iSpec, BoundID, PositionNbr,3,iReac)
  END IF
  iPart = iPart + 1
END DO
END SUBROUTINE SetInnerEnergies


!===================================================================================================================================
!> Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
FUNCTION CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
! MODULES
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)            :: xyzNod(3), Vector1(3), Vector2(3), ndist(3), midpoint(3)
REAL                        :: CalcPartPosTriaSurface(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal2(2), PartDistance
REAL, PARAMETER             :: eps_nontria=1.0E-6
!===================================================================================================================================
  CALL RANDOM_NUMBER(RandVal2)
  IF (TrackingMethod.NE.TRIATRACKING) THEN !prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
    RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2) !shift randVal off from 0 and 1
    DO WHILE (ABS(RandVal2(1)+RandVal2(2)-1.0).LT.eps_nontria) !sum must not be 1, since this corresponds to third egde
      CALL RANDOM_NUMBER(RandVal2)
      RandVal2 = RandVal2 + eps_nontria*(1. - 2.*RandVal2)
    END DO
  END IF
  CalcPartPosTriaSurface = xyzNod + Vector1 * RandVal2(1)
  CalcPartPosTriaSurface = CalcPartPosTriaSurface + Vector2 * RandVal2(2)
  PartDistance = ndist(1)*(CalcPartPosTriaSurface(1)-midpoint(1)) & !Distance from v1-v2
              + ndist(2)*(CalcPartPosTriaSurface(2)-midpoint(2)) &
              + ndist(3)*(CalcPartPosTriaSurface(3)-midpoint(3))
  IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
    CalcPartPosTriaSurface(1:3) = 2.*midpoint(1:3)-CalcPartPosTriaSurface(1:3)
  END IF

END FUNCTION CalcPartPosTriaSurface

!===================================================================================================================================
!> Determine the particle velocity of each inserted particle
!===================================================================================================================================
SUBROUTINE SetSurfacefluxVelocities(iSpec,iReac,iSF,iSample,jSample,iSide,BCSideID,SideID,NbrOfParticle,PartIns)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_Particle_Surfaces_Vars, ONLY : SurfMeshSubSideData, TriaSurfaceFlux
USE MOD_Particle_Surfaces,      ONLY : CalcNormAndTangBezier
USE MOD_SurfaceModel_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iSpec,iReac,iSF,iSample,jSample,iSide,BCSideID,SideID,NbrOfParticle,PartIns
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,PositionNbr,envelope,BoundID
REAL                             :: Vec3D(3), vec_nIn(1:3), vec_t1(1:3), vec_t2(1:3)
REAL                             :: a,zstar,RandVal1,RandVal2(2),RandVal3(3),u,RandN,RandN_save,Velo1,Velo2,Velosq,T,beta,z
LOGICAL                          :: RandN_in_Mem
REAL                             :: projFak                          ! VeloVecIC projected to inwards normal of tria
REAL                             :: Velo_t1                          ! Velo comp. of first orth. vector in tria
REAL                             :: Velo_t2                          ! Velo comp. of second orth. vector in tria
REAL                             :: VeloIC
!===================================================================================================================================
IF(PartIns.LT.1) RETURN

RandN_in_Mem=.FALSE.
envelope=-1
BoundID = SurfChemReac%Surfaceflux(iSF)%BC

IF (.NOT.SurfChemReac%Surfaceflux(iSF)%VeloIsNormal) THEN
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
  vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
END IF !.NOT.VeloIsNormal

VeloIC = SurfChemReac%Surfaceflux(iSF)%VeloIC
T = SurfChemReac%Surfaceflux(iSF)%MWTemperatureIC
a = SurfChemReac%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn(iSpec)
projFak = SurfChemReac%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak
Velo_t1 = SurfChemReac%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1
Velo_t2 = SurfChemReac%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2

! Set velocities
  !-- determine envelope for most efficient ARM [Garcia and Wagner 2006, JCP217-2]
  IF (ALMOSTZERO(VeloIC*projFak)) THEN
    ! Rayleigh distri
    envelope = 0
  ELSE IF (-0.4.LT.a .AND. a.LT.1.3) THEN
    ! low speed flow
    IF (a.LE.0.) THEN
      envelope = 1
    ELSE
      envelope = 3
    END IF !choose envelope based on flow direction
  ELSE
    ! high speed / general flow
    IF (a.LT.0.) THEN
      envelope = 2
    ELSE
      envelope = 4
    END IF !choose envelope based on flow direction
  END IF !low speed / high speed / rayleigh flow

  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      !-- 0a.: In case of side-normal velocities: calc n-/t-vectors at particle position, xi was saved in PartState(4:5)
      IF (SurfChemReac%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
        vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
        vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
        vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
      END IF !VeloIsNormal

      SELECT CASE(envelope)
      CASE(0)
        CALL RANDOM_NUMBER(RandVal1)
        zstar = -SQRT(-LOG(RandVal1))
      CASE(1)
        DO
          CALL RANDOM_NUMBER(RandVal2)
          zstar = -SQRT(a*a-LOG(RandVal2(1)))
          IF ( -(a-zstar)/zstar .GT. RandVal2(2)) THEN
            EXIT
          END IF
        END DO
      CASE(2)
        z = 0.5*(a-SQRT(a*a+2.))
        beta  = a-(1.0-a)*(a-z)
        DO
          CALL RANDOM_NUMBER(RandVal3)
          IF (EXP(-(beta*beta))/(EXP(-(beta*beta))+2.0*(a-z)*(a-beta)*EXP(-(z*z))).GT.RandVal3(1)) THEN
            zstar=-SQRT(beta*beta-LOG(RandVal3(2)))
            IF ( -(a-zstar)/zstar .GT. RandVal3(3)) THEN
              EXIT
            END IF
          ELSE
            zstar=beta+(a-beta)*RandVal3(2)
            IF ( (a-zstar)/(a-z)*EXP(z*z-(zstar*zstar)) .GT. RandVal3(3)) THEN
              EXIT
            END IF
          END IF
        END DO
      CASE(3)
        DO
          CALL RANDOM_NUMBER(RandVal3)
          u = RandVal3(1)
          IF ( a*SQRT(PI)/(a*SQRT(PI)+1+a*a) .GT. u) THEN

            zstar = -1./SQRT(2.)*ABS(RandN)
            EXIT
          ELSE IF ( (a*SQRT(PI)+1.)/(a*SQRT(PI)+1+a*a) .GT. u) THEN
            zstar = -SQRT(-LOG(RandVal3(2)))
            EXIT
          ELSE
            zstar = (1.0-SQRT(RandVal3(2)))*a
            IF (EXP(-(zstar*zstar)).GT.RandVal3(3)) THEN
              EXIT
            END IF
          END IF
        END DO
      CASE(4)
        DO
          CALL RANDOM_NUMBER(RandVal3)
          IF (1.0/(2.0*a*SQRT(PI)+1.0).GT.RandVal3(1)) THEN
            zstar=-SQRT(-LOG(RandVal3(2)))
          ELSE

              IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
                RandN = RandN_save
                RandN_in_Mem=.FALSE.
              ELSE
                Velosq = 2
                DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
                  CALL RANDOM_NUMBER(RandVal2)
                  Velo1 = 2.*RandVal2(1) - 1.
                  Velo2 = 2.*RandVal2(2) - 1.
                  Velosq = Velo1**2 + Velo2**2
                END DO
                RandN = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
                RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
                RandN_in_Mem=.TRUE.
              END IF

            zstar = 1./SQRT(2.)*RandN
          END IF
          IF ( (a-zstar)/a .GT. RandVal3(3)) THEN
            EXIT
          END IF
        END DO
      CASE DEFAULT
        CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong envelope in SetSurfacefluxVelocities!')
      END SELECT

      !-- 2.: sample normal directions and build complete velo-vector
      Vec3D(1:3) = vec_nIn(1:3) * SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC)*(a-zstar)
!      IF (.NOT.DoZigguratSampling) THEN !polar method
        Velosq = 2
        DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
          CALL RANDOM_NUMBER(RandVal2)
          Velo1 = 2.*RandVal2(1) - 1.
          Velo2 = 2.*RandVal2(2) - 1.
          Velosq = Velo1**2 + Velo2**2
        END DO
        Velo1 = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
        Velo2 = Velo2*SQRT(-2*LOG(Velosq)/Velosq)

      Vec3D(1:3) = Vec3D(1:3) + vec_t1(1:3) * ( Velo_t1+Velo1*SQRT(BoltzmannConst*T/Species(iSpec)%MassIC) )
      Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3) * ( Velo_t2+Velo2*SQRT(BoltzmannConst*T/Species(iSpec)%MassIC) )
      PartState(4:6,PositionNbr) = Vec3D(1:3)
    ELSE !PositionNbr .EQ. 0
      CALL abort(__STAMP__,'PositionNbr .EQ. 0!')
    END IF !PositionNbr .NE. 0
  END DO !i = ...NbrOfParticle 

END SUBROUTINE SetSurfacefluxVelocities

END MODULE MOD_Particle_SurfChemFlux
