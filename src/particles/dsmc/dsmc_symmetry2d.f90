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

MODULE MOD_DSMC_Symmetry2D
!===================================================================================================================================
!> Routines for 2D (planar/axisymmetric) simulations
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
PUBLIC :: DSMC_2D_InitVolumes, DSMC_2D_InitRadialWeighting, DSMC_2D_RadialWeighting, DSMC_2D_SetInClones, DSMC_2D_CalcSymmetryArea
PUBLIC :: CalcRadWeightMPF
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_2D_InitVolumes()
!===================================================================================================================================
!> Routine determines a symmetry side and calculates the 2D (area faces in symmetry plane) and axisymmetric volumes (cells are
!> revolved around the symmetry axis). The symmetry side will be used later on to determine in which two directions the quadtree
!> shall refine the mesh, skipping the z-dimension to avoid an unnecessary refinement.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: Pi
USE MOD_PreProc                 ,ONLY: PP_N
USE MOD_Mesh_Vars               ,ONLY: nElems, nBCSides, BC, SideToElem, SurfElem
USE MOD_Interpolation_Vars      ,ONLY: wGP
USE MOD_Particle_Vars           ,ONLY: Symmetry2DAxisymmetric
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_DSMC_Vars               ,ONLY: SymmetrySide
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: SideID, iLocSide, i, j, ElemID, iNode
REAL                            :: radius
LOGICAL                         :: SymmetryBCExists
!===================================================================================================================================

SymmetryBCExists = .FALSE.
ALLOCATE(SymmetrySide(1:nElems,1:2))                ! 1: GlobalSide, 2: LocalSide
SymmetrySide = -1

IF(.NOT.ALMOSTEQUALRELATIVE(GEO%zmaxglob,ABS(GEO%zminglob),1e-5)) THEN
  SWRITE(*,*) 'Maximum dimension in z:', GEO%zmaxglob
  SWRITE(*,*) 'Minimum dimension in z:', GEO%zminglob
  SWRITE(*,*) 'Deviation', (ABS(GEO%zmaxglob)-ABS(GEO%zminglob))/ABS(GEO%zminglob), ' > 1e-5'
  CALL abort(__STAMP__&
    ,'ERROR: Please orient your mesh with one cell in z-direction around 0, |z_min| = z_max !')
END IF

DO SideID=1,nBCSides
  ! Get the SideID of the symmetry sides defined as symmetry BC ('symmetric')
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(SideID))).EQ.PartBound%SymmetryBC) THEN
    ElemID = SideToElem(1,SideID)
    IF (ElemID.LT.1) THEN
      ElemID = SideToElem(2,SideID)
      iLocSide = SideToElem(4,SideID)
    ELSE
      iLocSide = SideToElem(3,SideID)
    END IF
    ! Exclude the symmetry axis (y=0)
    IF(MAXVAL(GEO%NodeCoords(2,GEO%ElemSideNodeID(:,iLocSide,ElemID))).GT.0.0) THEN
      ! The z-plane with the positive z component is chosen
      IF(MINVAL(GEO%NodeCoords(3,GEO%ElemSideNodeID(:,iLocSide,ElemID))).GT.(GEO%zmaxglob+GEO%zminglob)/2.) THEN
        IF(SymmetrySide(ElemID,1).GT.0) THEN
          CALL abort(__STAMP__&
            ,'ERROR: PICLas could not determine a unique symmetry surface for 2D/axisymmetric calculation!'//&
            ' Please orient your mesh with x as the symmetry axis and positive y as the second/radial direction!')
        END IF
        SymmetrySide(ElemID,1) = SideID
        SymmetrySide(ElemID,2) = iLocSide
        ! The volume calculated at this point (final volume for the 2D case) corresponds to the cell face area (z-dimension=1) in
        ! the xy-plane.
        GEO%Volume(ElemID) = 0.0
        DO j=0,PP_N; DO i=0,PP_N
          GEO%Volume(ElemID) = GEO%Volume(ElemID) + wGP(i)*wGP(j)*SurfElem(i,j,SideID)
        END DO; END DO
        ! Characteristic length is compared to the mean free path as the condition to refine the mesh. For the 2D/axisymmetric case
        ! the third dimension is not considered as particle interaction occurs in the xy-plane, effectively reducing the refinement
        ! requirement.
        GEO%CharLength(ElemID) = SQRT(GEO%Volume(ElemID))
        ! Axisymmetric case: The volume is multiplied by the circumference to get the volume of the ring. The cell face in the
        ! xy-plane is rotated around the x-axis. The radius is the middle point of the cell face.
        IF (Symmetry2DAxisymmetric) THEN
          radius = 0.
          DO iNode = 1, 4
            radius = radius + GEO%NodeCoords(2,GEO%ElemSideNodeID(iNode,iLocSide,ElemID))
          END DO
          radius = radius / 4.
          GEO%Volume(ElemID) = GEO%Volume(ElemID) * 2. * Pi * radius
        END IF
        SymmetryBCExists = .TRUE.
      END IF    ! y-coord greater 0.0
    END IF      ! Greater z-coord
  END IF
END DO

IF(.NOT.SymmetryBCExists) THEN
  CALL abort(__STAMP__&
    ,'At least one symmetric BC (in the xy-plane) has to be defined for 2D simulations')
END IF

! LocalVolume & MeshVolume: Recalculate the volume of the mesh of a single process and the total mesh volume
GEO%LocalVolume = SUM(GEO%Volume)
#ifdef MPI
CALL MPI_ALLREDUCE(GEO%LocalVolume,GEO%MeshVolume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
#else
GEO%MeshVolume=GEO%LocalVolume
#endif /*MPI*/

END SUBROUTINE DSMC_2D_InitVolumes

SUBROUTINE DSMC_2D_InitRadialWeighting()
!===================================================================================================================================
!> Read-in and initialize the variables required for the cloning procedures. Two modes with a delayed clone insertion are available:
!> 1: Insert the clones after the delay in the same chronological order as they were created
!> 2: Choose a random list of particles to insert after the delay buffer is full with clones
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Restart_Vars            ,ONLY: DoRestart
USE MOD_PARTICLE_Vars           ,ONLY: PDM
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting, ClonedParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Linear increasing weighting factor in the radial direction up to the domain boundary
RadialWeighting%PartScaleFactor = GETREAL('Particles-RadialWeighting-PartScaleFactor')
IF(RadialWeighting%PartScaleFactor.LT.1.) THEN
  CALL Abort(&
      __STAMP__,&
    'ERROR in 2D axisymmetric simulation: PartScaleFactor has to be greater than 1!',RealInfoOpt=RadialWeighting%PartScaleFactor)
END IF
RadialWeighting%CloneMode = GETINT('Particles-RadialWeighting-CloneMode')
RadialWeighting%CloneInputDelay = GETINT('Particles-RadialWeighting-CloneDelay')
! Cell local radial weighting (all particles have the same weighting factor within a cell)
RadialWeighting%CellLocalWeighting = GETLOGICAL('Particles-RadialWeighting-CellLocalWeighting')

RadialWeighting%NextClone = 0

SELECT CASE(RadialWeighting%CloneMode)
  CASE(1)
    IF(RadialWeighting%CloneInputDelay.LT.1) THEN
      CALL Abort(&
          __STAMP__,&
        'ERROR in 2D axisymmetric simulation: Clone delay should be greater than 0')
    END IF
    ALLOCATE(RadialWeighting%ClonePartNum(0:(RadialWeighting%CloneInputDelay-1)))
    ALLOCATE(ClonedParticles(1:INT(PDM%maxParticleNumber/RadialWeighting%CloneInputDelay),0:(RadialWeighting%CloneInputDelay-1)))
    RadialWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) RadialWeighting%CloneDelayDiff = 1
  CASE(2)
    IF(RadialWeighting%CloneInputDelay.LT.2) THEN
      CALL Abort(&
          __STAMP__,&
        'ERROR in 2D axisymmetric simulation: Clone delay should be greater than 1')
    END IF
    ALLOCATE(RadialWeighting%ClonePartNum(0:RadialWeighting%CloneInputDelay))
    ALLOCATE(ClonedParticles(1:INT(PDM%maxParticleNumber/RadialWeighting%CloneInputDelay),0:RadialWeighting%CloneInputDelay))
    RadialWeighting%ClonePartNum = 0
    IF(.NOT.DoRestart) RadialWeighting%CloneDelayDiff = 0
  CASE DEFAULT
    CALL Abort(&
        __STAMP__,&
      'ERROR in Radial Weighting of 2D/Axisymmetric: The selected cloning mode is not available! Choose between 1 and 2.'//&
        ' CloneMode=1: Delayed insertion of clones; CloneMode=2: Delayed randomized insertion of clones')
END SELECT

END SUBROUTINE DSMC_2D_InitRadialWeighting


SUBROUTINE DSMC_2D_RadialWeighting(iPart,iElem)
!===================================================================================================================================
!> Routine for the treatment of particles with enabled radial weighting (weighting factor is increasing linearly with increasing y)
!> 1.) Determine the new particle weight and decide whether to clone or to delete the particle
!> 2a.) Particle cloning, if the local weighting factor is smaller than the previous (particle travelling downwards)
!> 2b.) Particle deletion, if the local weighting factor is greater than the previous (particle travelling upwards)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting, DSMC, PartStateIntEn, useDSMC, CollisMode
USE MOD_DSMC_Vars               ,ONLY: ClonedParticles, VibQuantsPar, SpecDSMC, PolyatomMolDSMC, CollInf
USE MOD_Particle_Vars           ,ONLY: PartMPF, PDM, PartSpecies, PartState, Species, LastPartPos
USE MOD_TimeDisc_Vars           ,ONLY: iter
USE Ziggurat
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iPart, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPolyatMole, cloneIndex, DelayCounter
REAL                            :: DeleteProb, iRan, NewMPF, CloneProb, OldMPF
LOGICAL                         :: DoCloning
!===================================================================================================================================
DoCloning = .FALSE.
DeleteProb = 0.

IF (.NOT.(PartMPF(iPart).GT.Species(PartSpecies(iPart))%MacroParticleFactor)) RETURN

! 1.) Determine the new particle weight and decide whether to clone or to delete the particle
NewMPF = CalcRadWeightMPF(PartState(iPart,2), PartSpecies(iPart),iPart)
OldMPF = PartMPF(iPart)
CloneProb = (OldMPF/NewMPF)-INT(OldMPF/NewMPF)
CALL RANDOM_NUMBER(iRan)
IF((CloneProb.GT.iRan).AND.(NewMPF.LT.OldMPF)) THEN
  DoCloning = .TRUE.
  IF(INT(OldMPF/NewMPF).GT.1) THEN
    IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
    CALL Abort(&
        __STAMP__,&
      'ERROR in 2D axisymmetric simulation: More than one clone per particle is not allowed! Reduce the time step or'//&
        ' the radial weighting factor! Cloning probability is:',RealInfoOpt=CloneProb)
  END IF
END IF
PartMPF(iPart) = NewMPF

IF(DoCloning) THEN
  ! 2a.) Particle cloning, if the local weighting factor is smaller than the previous (particle travelling downwards)
  ! Get the list number to store the clones, depending on the chosen clone mode
  SELECT CASE(RadialWeighting%CloneMode)
  CASE(1)
  ! ######## Clone Delay ###################################################################################################
  ! Insertion of the clones after a defined delay, all clones are collected in a single list and inserted before boundary
  ! treatment in the next time step at their original positions
    DelayCounter = MOD((INT(iter,4)+RadialWeighting%CloneDelayDiff-1),RadialWeighting%CloneInputDelay)
  CASE(2)
  ! ######## Clone Random Delay #############################################################################################
  ! A list, which is RadialWeighting%CloneInputDelay + 1 long, is filled with clones to be inserted. After the list
  ! is full, NextClone gives the empty particle list, whose clones were inserted during the last SetInClones step
    IF((INT(iter,4)+RadialWeighting%CloneDelayDiff).LE.RadialWeighting%CloneInputDelay) THEN
      DelayCounter = INT(iter,4)+RadialWeighting%CloneDelayDiff
    ELSE
      DelayCounter = RadialWeighting%NextClone
    END IF
  END SELECT
  ! Storing the particle information
  RadialWeighting%ClonePartNum(DelayCounter) = RadialWeighting%ClonePartNum(DelayCounter) + 1
  cloneIndex = RadialWeighting%ClonePartNum(DelayCounter)
  ClonedParticles(cloneIndex,DelayCounter)%PartState(1:6)= PartState(iPart,1:6)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(1:2) = PartStateIntEn(iPart,1:2)
    IF(DSMC%ElectronicModel) ClonedParticles(cloneIndex,DelayCounter)%PartStateIntEn(3) =   PartStateIntEn(iPart,3)
    IF(SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
      IF(ALLOCATED(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)) &
        DEALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants)
      ALLOCATE(ClonedParticles(cloneIndex,DelayCounter)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      ClonedParticles(cloneIndex,DelayCounter)%VibQuants(:) = VibQuantsPar(iPart)%Quants(:)
    END IF
  END IF
  ClonedParticles(cloneIndex,DelayCounter)%Species = PartSpecies(iPart)
  ClonedParticles(cloneIndex,DelayCounter)%Element = iElem
  ClonedParticles(cloneIndex,DelayCounter)%LastPartPos(1:3) = LastPartPos(iPart,1:3)
  ClonedParticles(cloneIndex,DelayCounter)%WeightingFactor = PartMPF(iPart)
ELSE
! ######## Particle Delete #######################################################################################################
! 2b.) Particle deletion, if the local weighting factor is greater than the previous (particle travelling upwards)
  IF(NewMPF.GT.OldMPF) THEN
    DeleteProb = 1. - CloneProb
    IF (DeleteProb.GT.0.5) THEN
      IPWRITE(*,*) 'New weighting factor:', NewMPF, 'Old weighting factor:', OldMPF
      CALL abort(__STAMP__,&
        'ERROR in Radial Weighting of 2D/Axisymmetric: The deletion probability is higher than 0.5! Reduce the time step or'//&
        ' the radial weighting factor! Deletion probability is:',RealInfoOpt=DeleteProb)
    END IF
    CALL RANDOM_NUMBER(iRan)
    IF(DeleteProb.GT.iRan) THEN
      PDM%ParticleInside(iPart) = .FALSE.
      IF (CollInf%ProhibitDoubleColl) CollInf%OldCollPartner(iPart) = 0
    END IF
  END IF
END IF

END SUBROUTINE DSMC_2D_RadialWeighting


SUBROUTINE DSMC_2D_SetInClones()
!===================================================================================================================================
!> Insertion of cloned particles during the previous time steps. Clones insertion is delayed by at least one time step to avoid the
!> avalanche phenomenon (identical particles travelling on the same path, not colliding due to zero relative velocity).
!> 1.) Chose which list to insert depending on the clone mode
!> 2.) Insert the clones at the position they were created
!> 3.) Reset the list
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars             ,ONLY: ClonedParticles, PartStateIntEn, useDSMC, CollisMode, DSMC, RadialWeighting
  USE MOD_DSMC_Vars             ,ONLY: VibQuantsPar, SpecDSMC, PolyatomMolDSMC, SamplingActive
  USE MOD_Particle_Vars         ,ONLY: PDM, PEM, PartSpecies, PartState, LastPartPos, PartMPF, WriteMacroVolumeValues, VarTimeStep
  USE MOD_Particle_VarTimeStep  ,ONLY: CalcVarTimeStep
  USE MOD_TimeDisc_Vars         ,ONLY: iter
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart, PositionNbr, iPolyatMole, DelayCounter
REAL                            :: iRan
!===================================================================================================================================

! 1.) Chose which list to insert depending on the clone mode
SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  ! During the first iterations the delay counter refers to the empty clone array (which is filled during the following tracking)
  ! Afterwards, the MODULUS counts up from zero to CloneInputDelay-1
  DelayCounter = MOD((INT(iter,4)+RadialWeighting%CloneDelayDiff-1),RadialWeighting%CloneInputDelay)
CASE(2)
  ! During the first iterations, check if number of iterations is less than the input delay and leave routine. Afterwards, a
  ! random clone list from the previous time steps is chosen.
  IF((INT(iter,4)+RadialWeighting%CloneDelayDiff).GT.RadialWeighting%CloneInputDelay) THEN
    CALL RANDOM_NUMBER(iRan)
    ! Choosing random clone between 0 and CloneInputDelay
    DelayCounter = INT((RadialWeighting%CloneInputDelay+1)*iRan)
    DO WHILE (DelayCounter.EQ.RadialWeighting%NextClone)
      CALL RANDOM_NUMBER(iRan)
      DelayCounter = INT((RadialWeighting%CloneInputDelay+1)*iRan)
    END DO
    ! Save the chosen list as the next available list to store clones in the next time step
    RadialWeighting%NextClone = DelayCounter
  ELSE
    RETURN
  END IF
END SELECT

IF(RadialWeighting%ClonePartNum(DelayCounter).EQ.0) RETURN

! 2.) Insert the clones at the position they were created
DO iPart = 1, RadialWeighting%ClonePartNum(DelayCounter)
  PDM%ParticleVecLength = PDM%ParticleVecLength + 1
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1
  PositionNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
  IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
    CALL Abort(&
       __STAMP__,&
      'ERROR in 2D axisymmetric simulation: New Particle Number greater max Part Num!')
  END IF
  ! Copy particle parameters
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PDM%IsNewPart(PositionNbr) = .TRUE.
  PDM%dtFracPush(PositionNbr) = .FALSE.
  PartState(PositionNbr,1:5) = ClonedParticles(iPart,DelayCounter)%PartState(1:5)
  ! Creating a relative velocity in the z-direction
  PartState(PositionNbr,6) = - ClonedParticles(iPart,DelayCounter)%PartState(6)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    PartStateIntEn(PositionNbr,1:2) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(1:2)
    IF(DSMC%ElectronicModel) PartStateIntEn(PositionNbr,3) = ClonedParticles(iPart,DelayCounter)%PartStateIntEn(3)
    IF(SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(ClonedParticles(iPart,DelayCounter)%Species)%SpecToPolyArray
      IF(ALLOCATED(VibQuantsPar(PositionNbr)%Quants)) DEALLOCATE(VibQuantsPar(PositionNbr)%Quants)
      ALLOCATE(VibQuantsPar(PositionNbr)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
      VibQuantsPar(PositionNbr)%Quants(:) = ClonedParticles(iPart,DelayCounter)%VibQuants(:)
    END IF
  END IF
  PartSpecies(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Species
  PEM%Element(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Element
  PEM%lastElement(PositionNbr) = ClonedParticles(iPart,DelayCounter)%Element
  LastPartPos(PositionNbr,1:3) = ClonedParticles(iPart,DelayCounter)%LastPartPos(1:3)
  PartMPF(PositionNbr) =  ClonedParticles(iPart,DelayCounter)%WeightingFactor
  IF (VarTimeStep%UseVariableTimeStep) THEN
    VarTimeStep%ParticleTimeStep(PositionNbr) = CalcVarTimeStep(PartState(PositionNbr,1),PartState(PositionNbr,2),&
                                                                PEM%Element(PositionNbr))
  END IF
  ! Counting the number of clones per cell
  IF(SamplingActive.OR.WriteMacroVolumeValues) THEN
    IF(DSMC%CalcQualityFactors) DSMC%QualityFacSamp(PEM%Element(PositionNbr),5) = &
                                            DSMC%QualityFacSamp(PEM%Element(PositionNbr),5) + 1
  END IF
END DO

! 3.) Reset the list
RadialWeighting%ClonePartNum(DelayCounter) = 0

END SUBROUTINE DSMC_2D_SetInClones


REAL FUNCTION DSMC_2D_CalcSymmetryArea(iLocSide,iElem, ymin, ymax)
!===================================================================================================================================
!> Calculates the actual area of an element for 2D simulations (plane/axisymmetric) regardless of the mesh dimension in z
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars          ,ONLY: Pi
USE MOD_Particle_Vars         ,ONLY: Symmetry2DAxisymmetric
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: iElem,iLocSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, OPTIONAL, INTENT(OUT)   :: ymax,ymin
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iNode
REAL                          :: P(1:2,1:4), Pmin(2), Pmax(2), Length, MidPoint
!===================================================================================================================================

Pmin = HUGE(Pmin)
Pmax = -HUGE(Pmax)

DO iNode = 1,4
  P(1:2,iNode) = GEO%NodeCoords(1:2,GEO%ElemSideNodeID(iNode,iLocSide,iElem))
END DO

Pmax(1) = MAXVAL(P(1,:))
Pmax(2) = MAXVAL(P(2,:))
Pmin(1) = MINVAL(P(1,:))
Pmin(2) = MINVAL(P(2,:))

IF (PRESENT(ymax).AND.PRESENT(ymin)) THEN
  ymin = Pmin(2)
  ymax = Pmax(2)
END IF

Length = SQRT((Pmax(1)-Pmin(1))**2 + (Pmax(2)-Pmin(2))**2)

MidPoint = (Pmax(2)+Pmin(2)) / 2.
IF(Symmetry2DAxisymmetric) THEN
  DSMC_2D_CalcSymmetryArea = Length * MidPoint * Pi * 2.
  ! Area of the cells on the rotational symmetry axis is set to one
  IF(.NOT.(DSMC_2D_CalcSymmetryArea.GT.0.0)) DSMC_2D_CalcSymmetryArea = 1.
ELSE
  DSMC_2D_CalcSymmetryArea = Length
END IF
RETURN

END FUNCTION DSMC_2D_CalcSymmetryArea


REAL FUNCTION CalcRadWeightMPF(yPos, iSpec, iPart)
!===================================================================================================================================
!> Determines the weighting factor when using an additional radial weighting for axisymmetric simulations. Linear increase from the
!> rotational axis (y=0) to the outer domain boundary (y=ymax).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
USE MOD_Particle_Vars           ,ONLY: Species, PEM
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                :: yPos
INTEGER, INTENT(IN)             :: iSpec
INTEGER, OPTIONAL,INTENT(IN)    :: iPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                 :: yPosIn
!===================================================================================================================================

IF(RadialWeighting%CellLocalWeighting.AND.PRESENT(iPart)) THEN
  yPosIn = GEO%ElemMidPoint(2,PEM%Element(iPart))
ELSE
  yPosIn = yPos
END IF

CalcRadWeightMPF = (1. + yPosIn/GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))*Species(iSpec)%MacroParticleFactor

RETURN

END FUNCTION CalcRadWeightMPF

END MODULE MOD_DSMC_Symmetry2D
