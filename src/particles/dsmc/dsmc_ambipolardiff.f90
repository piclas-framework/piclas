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

MODULE MOD_DSMC_AmbipolarDiffusion
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
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
PUBLIC :: InitializeVariablesAmbipolarDiff, AD_SetInitElectronVelo, AD_InsertParticles, AD_DeleteParticles
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeVariablesAmbipolarDiff()
!===================================================================================================================================
!> Ambipolar Diffusion: Electrons are attached to and move with the ions, but still have their own velocity vector for collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars        ,ONLY: ElementaryCharge
USE MOD_Particle_Vars       ,ONLY: nSpecies,Species, PDM
USE MOD_DSMC_Vars           ,ONLY: useDSMC, DSMC, AmbipolElecVelo
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSpec
!===================================================================================================================================
IF(useDSMC) THEN
  DSMC%DoAmbipolarDiff = GETLOGICAL('Particles-DSMC-AmbipolarDiffusion')
  IF (DSMC%DoAmbipolarDiff) THEN
    DSMC%AmbiDiffElecSpec = 0
    DO iSpec = 1, nSpecies
      IF (Species(iSpec)%ChargeIC.GT.0.0) CYCLE
      IF(NINT(Species(iSpec)%ChargeIC/(-ElementaryCharge)).EQ.1) DSMC%AmbiDiffElecSpec=iSpec
    END DO
    IF(DSMC%AmbiDiffElecSpec.EQ.0) THEN
      CALL abort(__STAMP__&
          ,'ERROR: No electron species found for ambipolar diffusion: ' &
          ,IntInfoOpt=DSMC%AmbiDiffElecSpec)
    END IF
    IF(.NOT.ALLOCATED(AmbipolElecVelo)) ALLOCATE(AmbipolElecVelo(PDM%maxParticleNumber))
  END IF
END IF

END SUBROUTINE InitializeVariablesAmbipolarDiff


SUBROUTINE AD_SetInitElectronVelo(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
!> Initialize the electron velocity vector during the initial particle insertion
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_part_tools              ,ONLY: BuildTransGaussNums
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: FractNbr,iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)           :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! LOCAL VARIABLES
INTEGER                         :: i, PositionNbr
CHARACTER(30)                   :: velocityDistribution
REAL                            :: VeloIC, VeloVecIC(3), maxwellfac, VeloVecNorm
REAL                            :: iRanPart(3, NbrOfParticle), Vec3D(3)
!===================================================================================================================================
IF(NbrOfParticle.LT.1) RETURN
IF(Species(FractNbr)%ChargeIC.LE.0.0) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PDM%maxParticleNumber!')
END IF

velocityDistribution=Species(FractNbr)%Init(iInit)%velocityDistribution
VeloIC=Species(FractNbr)%Init(iInit)%VeloIC
VeloVecIC=Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
VeloVecNorm = VECNORM(VeloVecIC(1:3))
IF (VeloVecNorm.GT.0.0) THEN
  VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
END IF

SELECT CASE(TRIM(velocityDistribution))
CASE('constant')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloVecIC(1:3) * VeloIC 
    END IF
  END DO
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_maxwell_lpn(DSMC%AmbiDiffElecSpec, Vec3D, Temperature=Species(FractNbr)%Init(iInit)%MWTemperatureIC)
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + Vec3D(1:3)
    END IF
  END DO
CASE('maxwell')
  CALL BuildTransGaussNums(NbrOfParticle, iRanPart)
  maxwellfac = SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(DSMC%AmbiDiffElecSpec)%MassIC)
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + iRanPart(1:3,i)*maxwellfac
    END IF
  END DO
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'Velo-Distri not implemented for ambipolar diffusion!')
END SELECT

END SUBROUTINE AD_SetInitElectronVelo


SUBROUTINE AD_InsertParticles(iPartIndx_Node, nPart, iPartIndx_NodeTotalAmbi, TotalPartNum)
!===================================================================================================================================
!> Creating electrons for each actual ion simulation particle, using the stored velocity vector for the electron
!===================================================================================================================================
! MODULES
USE MOD_Globals                
USE MOD_DSMC_Vars               ,ONLY: BGGas, CollisMode, DSMC, PartStateIntEn, AmbipolElecVelo, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: DSMCSumOfFormedParticles, newAmbiParts, iPartIndx_NodeNewAmbi, DSMC_RHS
USE MOD_PARTICLE_Vars           ,ONLY: PDM, PartSpecies, PartState, PEM, Species, VarTimeStep, PartMPF, Symmetry
USE MOD_Particle_Tracking       ,ONLY: ParticleInsideCheck
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)             :: nPart, TotalPartNum
INTEGER,INTENT(INOUT)             :: iPartIndx_Node(1:nPart)
INTEGER,INTENT(INOUT),ALLOCATABLE :: iPartIndx_NodeTotalAmbi(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr, iLoop, nNewElectrons, IonIndX(nPart), iElem, PartNum
REAL              :: MaxPos(3), MinPos(3), Vec3D(3), RandomPos(3)
LOGICAL           :: InsideFlag
!===================================================================================================================================

MaxPos = -HUGE(MaxPos)
MinPos = HUGE(MinPos)
iNewPart=0
PositionNbr = 0
nNewElectrons = 0

DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  MaxPos(1) = MAX(MaxPos(1),PartState(1,iPart))
  MaxPos(2) = MAX(MaxPos(2),PartState(2,iPart))
  MaxPos(3) = MAX(MaxPos(3),PartState(3,iPart))
  MinPos(1) = MIN(MinPos(1),PartState(1,iPart))
  MinPos(2) = MIN(MinPos(2),PartState(2,iPart))
  MinPos(3) = MIN(MinPos(3),PartState(3,iPart))
  IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
    nNewElectrons = nNewElectrons + 1
    IonIndX(nNewElectrons) = iPart
  END IF
END DO
ALLOCATE(iPartIndx_NodeTotalAmbi(nPart + nNewElectrons))
iPartIndx_NodeTotalAmbi(1:nPart) = iPartIndx_Node(1:nPart)
TotalPartNum = nPart

DO iLoop = 1, nNewElectrons
  DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
  PositionNbr = PDM%nextFreePosition(DSMCSumOfFormedParticles+PDM%CurrentNextFreePosition)
  IF (PositionNbr.EQ.0) THEN
    CALL Abort(&
__STAMP__&
,'ERROR in Ambipolar Diffusion: MaxParticleNumber too small!')
  END IF
  InsideFlag=.FALSE.
  iElem = PEM%GlobalElemID(iPartIndx_Node(1))
  DO WHILE(.NOT.InsideFlag)
    CALL RANDOM_NUMBER(Vec3D(1:3))
    RandomPos(1:3) = MinPos(1:3)+Vec3D(1:3)*(MaxPos(1:3)-MinPos(1:3))
    IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
    IF(Symmetry%Order.LE.1) RandomPos(2) = 0.
    InsideFlag = ParticleInsideCheck(RandomPos,iPart,iElem)
  END DO
  PartState(1:3,PositionNbr) = RandomPos(1:3)
  PartSpecies(PositionNbr) = DSMC%AmbiDiffElecSpec
  PEM%GlobalElemID(PositionNbr) = iElem
  PDM%ParticleInside(PositionNbr) = .true.
  PartState(4:6,PositionNbr) = AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(1:3)
  DSMC_RHS(1:3,PositionNbr) = 0.0
  DEALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
    PartStateIntEn( 1,PositionNbr) = 0.
    PartStateIntEn( 2,PositionNbr) = 0.
    IF (DSMC%ElectronicModel.GT.0)   PartStateIntEn( 3,PositionNbr) = 0.
  END IF
  IF(RadialWeighting%DoRadialWeighting) PartMPF(PositionNbr) = PartMPF(IonIndX(iLoop))
  IF(VarTimeStep%UseVariableTimeStep) VarTimeStep%ParticleTimeStep(PositionNbr) = VarTimeStep%ParticleTimeStep(IonIndX(iLoop))
  iPartIndx_NodeTotalAmbi(nPart+iLoop) = PositionNbr
END DO
! Output variable
TotalPartNum = nPart + nNewElectrons
! Variable used for allocation
PartNum = TotalPartNum
! In case of the background gas, additional particles will be added
IF(BGGas%NumberOfSpecies.GT.0) PartNum = PartNum + TotalPartNum

newAmbiParts = 0
IF (ALLOCATED(iPartIndx_NodeNewAmbi)) DEALLOCATE(iPartIndx_NodeNewAmbi)
ALLOCATE(iPartIndx_NodeNewAmbi(PartNum))

END SUBROUTINE AD_InsertParticles


SUBROUTINE AD_DeleteParticles(iPartIndx_Node, nPart_opt)
!===================================================================================================================================
!> Deletes all electron created for the collision process and saves their velocity vector to the ion they are attached to
!> 1) Counting the ions/electrons within particles that already existed within the cell (but may have changed their species)
!> 2) Counting the ions/electrons within particles that were newly created during the time step
!> 3) Assinging each ion an electron, saving its velocity vector and deleting it
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PARTICLE_Vars       ,ONLY: PDM, PartSpecies, Species, PartState, PEM
USE MOD_Particle_Analyze    ,ONLY: PARTISELECTRON
USE MOD_DSMC_Vars           ,ONLY: AmbipolElecVelo, newAmbiParts, iPartIndx_NodeNewAmbi, DSMC_RHS
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: iPartIndx_Node(:), nPart_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart, iLoop, nElectron, nIon, nPart
INTEGER, ALLOCATABLE        :: ElecIndx(:), IonIndX(:)
!===================================================================================================================================
nElectron =0; nIon = 0

IF(PRESENT(nPart_opt)) THEN
  ALLOCATE(ElecIndx(2*nPart_opt), IonIndX(2*nPart_opt))
  nPart = nPart_opt
ELSE
  ALLOCATE(ElecIndx(newAmbiParts), IonIndX(newAmbiParts))
  nPart = 0
END IF

! 1) Treating the particles that already existed within the cell (but may have changed their species)
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    IF(PARTISELECTRON(iPart)) THEN
      nElectron = nElectron + 1
      ElecIndx(nElectron) = iPart
    END IF
    IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
      nIon = nIon + 1
      IonIndX(nIon) = iPart
    END IF
  END IF
END DO
! 2) Treating particles that were newly created during the time step
DO iLoop = 1, newAmbiParts
  iPart = iPartIndx_NodeNewAmbi(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    IF(PARTISELECTRON(iPart)) THEN
      nElectron = nElectron + 1
      ElecIndx(nElectron) = iPart
    END IF
    IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
      nIon = nIon + 1
      IonIndX(nIon) = iPart
    END IF
  END IF
END DO
! Sanity check: number of electrons and ions has to be identical
IF(nIon.NE.nElectron) THEN
  IPWRITE(*,*) 'Initial number of particles:', nPart, 'Number of new particles', newAmbiParts
  DO iLoop = 1, nPart
    iPart = iPartIndx_Node(iLoop)
    IPWRITE(*,*) 'Index, Element ID, Particle Inside?:', iPart, PDM%ParticleInside(iPart), PEM%GlobalElemID(iPart)
    IPWRITE(*,*) 'Species, Charge:', PartSpecies(iPart), Species(PartSpecies(iPart))%ChargeIC
  END DO
  IPWRITE(*,*) 'List of new particles:'
  DO iLoop = 1, newAmbiParts
    iPart = iPartIndx_NodeNewAmbi(iLoop)
    IPWRITE(*,*) 'Index, Element ID, Particle Inside?:', iPart, PDM%ParticleInside(iPart), PEM%GlobalElemID(iPart)
    IPWRITE(*,*) 'Species, Charge:', PartSpecies(iPart), Species(PartSpecies(iPart))%ChargeIC
  END DO
  CALL abort(__STAMP__,&
      'ERROR: Number of electrons and ions is not equal for ambipolar diffusion: ',IntInfoOpt=nIon-nElectron)
END IF

! 3) Assinging each ion an electron, saving its velocity vector and deleting it
DO iLoop = 1, nElectron
  IF (ALLOCATED(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)) DEALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)
  ALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(3))
  AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(1:3) = PartState(4:6,ElecIndx(iLoop)) + DSMC_RHS(1:3,ElecIndx(iLoop))
  PDM%ParticleInside(ElecIndx(iLoop)) = .FALSE.
END DO

DEALLOCATE(iPartIndx_NodeNewAmbi)

END SUBROUTINE AD_DeleteParticles

END MODULE MOD_DSMC_AmbipolarDiffusion
