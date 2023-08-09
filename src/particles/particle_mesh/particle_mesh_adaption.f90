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

MODULE MOD_Cell_Adaption
!===================================================================================================================================
!> Module containing routines for the recursive octree cell refinement algorithm for the BGK and FP-Flow particle methods.
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
PUBLIC :: DefineParametersAdaptMesh, Init_MeshAdaption, MeshAdaption, CalcGradients
!===================================================================================================================================

CONTAINS

SUBROUTINE DefineParametersAdaptMesh()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE

CALL prms%SetSection("Mesh Adaption")
CALL prms%CreateIntOption( 'Part-MeshAdapt-MinPartNum', 'Minimum particle number per subcell for the mesh adaption', '6')
CALL prms%CreateRealOption('Part-MeshAdapt-RefineLimitGradient', 'Gradient above which the subcell adaption is called', '0.05')                       

END SUBROUTINE DefineParametersAdaptMesh

SUBROUTINE Init_MeshAdaption()
!===================================================================================================================================
! Read-in of octree variables and building of the octree for a node depth of 2 during the initialization
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars     ,ONLY: AdaptMesh, MinPartCell, AdaptOrder, RefineFactorGrad
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem
!===================================================================================================================================
! READ-IN of the minimum particle number per subcell
MinPartCell = GETINT('Part-MeshAdapt-MinPartNum')
IF (MinPartCell.LT.5) CALL abort(__STAMP__,'ERROR: Given minimum particle number per cell is less than 5!')

ALLOCATE(AdaptMesh(nElems))
RefineFactorGrad = GETREAL('Part-MeshAdapt-RefineLimitGradient')

! Initialization of the relative orientation of the basis vectors in the physical and reference space
DO iElem=1, nElems
  AdaptMesh(iElem)%CellOrientation(1:3) = (/1,2,3/)
  CALL DefineElementOrientation(iElem)
END DO

END SUBROUTINE Init_MeshAdaption


SUBROUTINE MeshAdaption(iElem) 
!===================================================================================================================================
!> Adaption routine for an anisotropic mesh adaption, based on flow gradients
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_BGK_Vars
USE MOD_FPFlow_Vars
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_TimeDisc_Vars          ,ONLY: TEnd, Time
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars     
USE MOD_Particle_Vars          ,ONLY: PEM, WriteMacroVolumeValues, DoVirtualCellMerge, VirtMergedCells
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)      :: iElem 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE     :: PartIndx(:)
INTEGER                  :: nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem
INTEGER                  :: PartNumAdapt, SubcellID, GlobalElemID, CNElemID, nPart, iLoop, iPart
REAL, ALLOCATABLE        :: MappedPart_Subcell(:,:,:), SubcellVolTemp(:)
REAL                     :: MaxGradient
LOGICAL                  :: DoMergedCell
!===================================================================================================================================
DoMergedCell = .FALSE.

CNElemID = GetCNElemID(iElem+offSetElem)

! Calculation of the quality factors
IF(DSMC%CalcQualityFactors) THEN
  IF(BGKInitDone) THEN
    BGK_MeanRelaxFactorCounter = 0; BGK_MeanRelaxFactor = 0.; BGK_MaxRelaxFactor = 0.; BGK_MaxRotRelaxFactor = 0.
    BGK_PrandtlNumber=0.; BGK_ExpectedPrandtlNumber=0.; BGK_Viscosity=0.; BGK_ThermalConductivity=0.
  END IF
  IF(FPInitDone) THEN
    FP_MeanRelaxFactorCounter = 0; FP_MeanRelaxFactor = 0.; FP_MaxRelaxFactor = 0.; FP_MaxRotRelaxFactor = 0.; FP_PrandtlNumber = 0.
  END IF
END IF

! Skip cell if number of particles is less than 2, create particle list (iPartIndx_Node) and sum-up bulk velocity
nPart = PEM%pNumber(iElem)
IF (DoVirtualCellMerge) THEN
  IF(VirtMergedCells(iElem)%isMerged) RETURN
  IF(VirtMergedCells(iElem)%NumOfMergedCells.GT.0) THEN
    nPartMerged = nPart
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      nPartMerged = nPartMerged + PEM%pNumber(VirtMergedCells(iElem)%MergedCellID(iMergeElem))
    END DO
    IF (nPartMerged.LE.1) RETURN
    ALLOCATE(PartIndx(nPartMerged))
    iPart = PEM%pStart(iElem)
    iLoopLoc = 0
    DO iLoop = 1, nPart
      iLoopLoc = iLoopLoc + 1
      PartIndx(iLoopLoc) = iPart
      iPart = PEM%pNext(iPart)
    END DO
    DO iMergeElem = 1, VirtMergedCells(iElem)%NumOfMergedCells
      locElem = VirtMergedCells(iElem)%MergedCellID(iMergeElem)
      nPartLoc = PEM%pNumber(locElem)
      iPart = PEM%pStart(locElem)
      DO iLoop = 1, nPartLoc
        iLoopLoc = iLoopLoc + 1
        PartIndx(iLoopLoc) = iPart
        iPart = PEM%pNext(iPart)
      END DO
    END DO
    DoMergedCell = .TRUE.
  END IF  
ELSE IF ((nPart.EQ.0).OR.(nPart.EQ.1)) THEN
  RETURN
END IF

! Call to the collision routines for the merged cells
IF (DoMergedCell) THEN
#if (PP_TimeDiscMethod==300)  
  CALL FP_CollisionOperator(PartIndx, nPartMerged, VirtMergedCells(iElem)%MergedVolume)
#else
  ! IF (BGKMovingAverage) THEN
  !   CALL BGK_CollisionOperator(PartIndx, nPartMerged,VirtMergedCells(iElem)%MergedVolume, &
  !           ElemNodeAveraging(iElem)%Root%AverageValues(:))
  ! ELSE
    CALL BGK_CollisionOperator(PartIndx, nPartMerged, VirtMergedCells(iElem)%MergedVolume)
 ! END IF
#endif
  DEALLOCATE(PartIndx)
END IF

! Calculate the relative gradient of density, velocity, temperature and pressure to check for an adaption
! Calculate the temperature gradient in all three directions to determine the order of refinement directions
CALL CalcGradients(iElem,MaxGradient)

AdaptMesh(iElem)%SplitOrder = INT(nPart/MinPartCell)
AdaptMesh(iElem)%SplitOrder = INT(AdaptMesh(iElem)%SplitOrder)**(1./3.)

! Adaption requirements: enough particles for at least two subcells and gradient larger than predefined limit
IF ((nPart.GE.(2.*MinPartCell)).AND.(MaxGradient.GE.RefineFactorGrad).AND.(.NOT.DoMergedCell)) THEN
  ! Calculate the splitting order (1-8, 2-27, 3-64...)
  AdaptMesh(iElem)%SplitOrder = INT(nPart/MinPartCell)
  AdaptMesh(iElem)%SplitOrder = INT(AdaptMesh(iElem)%SplitOrder)**(1./3.)

  IF (AdaptMesh(iElem)%SplitOrder.GT.0) THEN
    CALL SubcellAssingement(iElem) ! Split the element into the defined number of subcells
  END IF

ELSE IF (.NOT.DoMergedCell) THEN! Normal routine without merged or split cells

  iPart = PEM%pStart(iElem)
  ALLOCATE(PartIndx(nPart))
  DO iLoop = 1, nPart
    PartIndx(iLoop) = iPart
    iPart = PEM%pNext(iPart)
  END DO

  ! No mesh adaption, direct call to the collision operators
#if (PP_TimeDiscMethod==300)
  CALL FP_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
#else
  IF (BGKMovingAverage) THEN
    CALL BGK_CollisionOperator(PartIndx, nPart,ElemVolume_Shared(CNElemID), ElemNodeAveraging(iElem)%Root%AverageValues(:))
  ELSE
    CALL BGK_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
  END IF
#endif
  DEALLOCATE(PartIndx)
END IF

! Sampling of the quality factors
IF(DSMC%CalcQualityFactors) THEN
  IF((Time.GE.(1-DSMC%TimeFracSamp)*TEnd).OR.WriteMacroVolumeValues) THEN
    IF(BGKInitDone) THEN
      BGK_QualityFacSamp(1,iElem) = BGK_QualityFacSamp(1,iElem) + BGK_MeanRelaxFactor
      BGK_QualityFacSamp(2,iElem) = BGK_QualityFacSamp(2,iElem) + REAL(BGK_MeanRelaxFactorCounter)
      BGK_QualityFacSamp(3,iElem) = BGK_QualityFacSamp(3,iElem) + BGK_MaxRelaxFactor
      BGK_QualityFacSamp(4,iElem) = BGK_QualityFacSamp(4,iElem) + 1.
      BGK_QualityFacSamp(5,iElem) = BGK_QualityFacSamp(5,iElem) + BGK_MaxRotRelaxFactor
      BGK_QualityFacSamp(6,iElem) = BGK_QualityFacSamp(6,iElem) + BGK_PrandtlNumber
      BGK_QualityFacSamp(7,iElem) = BGK_QualityFacSamp(7,iElem) + BGK_ExpectedPrandtlNumber
      BGK_QualityFacSamp(8,iElem) = BGK_QualityFacSamp(8,iElem) + BGK_Viscosity
      BGK_QualityFacSamp(9,iElem) = BGK_QualityFacSamp(9,iElem) + BGK_ThermalConductivity
    END IF
    IF(FPInitDone) THEN
      FP_QualityFacSamp(1,iElem) = FP_QualityFacSamp(1,iElem) + FP_MeanRelaxFactor
      FP_QualityFacSamp(2,iElem) = FP_QualityFacSamp(2,iElem) + REAL(FP_MeanRelaxFactorCounter)
      FP_QualityFacSamp(3,iElem) = FP_QualityFacSamp(3,iElem) + FP_MaxRelaxFactor
      FP_QualityFacSamp(4,iElem) = FP_QualityFacSamp(4,iElem) + 1.
      FP_QualityFacSamp(5,iElem) = FP_QualityFacSamp(5,iElem) + FP_MaxRotRelaxFactor
      FP_QualityFacSamp(6,iElem) = FP_QualityFacSamp(6,iElem) + FP_PrandtlNumber
    END IF
  END IF
END IF

END SUBROUTINE MeshAdaption


SUBROUTINE SubcellAssingement(iElem)  
!============================================================================================================================
!============================================================================================================================
! use MODULES                
USE MOD_BGK_Vars
USE MOD_FPFlow_Vars
USE MOD_Preproc
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars     
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_DSMC_Vars              ,ONLY: DSMC                   
USE MOD_Particle_Vars          ,ONLY: PEM, PartState, LastPartPos, PartPosRef
USE MOD_Mesh_Vars              ,ONLY: nElems, offSetElem, sJ
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Interpolation_Vars     ,ONLY: xGP, wBary
USE MOD_Basis                  ,ONLY: InitializeVandermonde
!-----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: iCell, nCells, iPart, nPart, iLoop, iDir, iVal, nSplit, CellCoord(3), GlobalElemID
INTEGER                       :: j, k, l, nOrder
REAL, ALLOCATABLE             :: Subcell_xGP(:),  Subcell_Vdm(:,:), DetJac(:,:,:,:), LocalVdm(:,:)
REAL                          :: Subcell_wGP, PartPosMapped(3), DetLocal(1,0:PP_N,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! Number of subcells from the split order
nCells = (AdaptMesh(iElem)%SplitOrder+1)**3

! Particle number in the element
nPart = PEM%pNumber(iElem)
iPart = PEM%pStart(iElem) 

! Allocate the particle number per subcell
ALLOCATE(SubPartNum(nCells)) 
SubPartNum = 0

! Allocate the particle indices for the subcell
ALLOCATE(SubPartIndx(nCells,nPart)) 
SubPartIndx = 0

! Allocate the subvolumes for each cell
ALLOCATE(SubVolume(nCells)) 

nOrder = AdaptMesh(iElem)%SplitOrder
! Allocate: Local determinant and vandermonde
ALLOCATE(DetJac(1,0:nOrder,0:nOrder,0:nOrder))
ALLOCATE(LocalVdm(0:nOrder,0:nOrder))

! Initalize interpolation points and subcell vandermonde
ALLOCATE(Subcell_xGP(0:nOrder))
ALLOCATE(Subcell_Vdm(0:nOrder, 0:nOrder))

! Coordinates of the interpolation points -> midpoints of the subcells
DO iVal = 0, nOrder
  Subcell_xGP(iVal) = (2.*REAL(iVal) + 1.)/(REAL(nOrder+1)) - 1.   
END DO
! Weights of the interpolation points, total volume = 8
Subcell_wGP = 2./REAL(1.0+nOrder)

! Local determinant from sJ (global value?)
DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
  DetLocal(1,j,k,l)=1./sJ(j,k,l,iElem)
END DO; END DO; END DO

! Call to the functions
CALL InitializeVandermonde(PP_N,nOrder,wBary,xGP,Subcell_xGP,Subcell_Vdm)
CALL ChangeBasis3D(1,PP_N,nOrder,Subcell_Vdm,DetLocal(:,:,:,:),DetJac(:,:,:,:))

! Build the connection of the subcells in the form of the Peano curve
CALL BuildPeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir)

DO iCell=1, nCells
  j = PeanoCurve(iCell,1)
  k = PeanoCurve(iCell,2)
  l = PeanoCurve(iCell,3)
  ! Calculate the subcell volume from the weights and the Jacobian
  Subvolume(iCell) = DetJac(1,j,k,l)*Subcell_wGP**3
END DO

DEALLOCATE(Subcell_xGP)
DEALLOCATE(Subcell_Vdm)
DEALLOCATE(DetJac)
DEALLOCATE(LocalVdm)

! Build the connection of the subcells in the form of the Peano curve
CALL BuildPeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir)

GlobalElemID = iElem + offsetElem
DO iLoop = 1, nPart
  ! Map Particle to -1|1 space 
  IF (TrackingMethod.EQ.REFMAPPING) THEN ! Refmapping
    PartPosMapped(1:3) = PartPosRef(1:3,iPart)
  ELSE
#if PP_TimeDiscMethod==300
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosMapped(1:3,iPart),GlobalElemID)
#else
    ! Attention: LastPartPos is the reference position here
    PartPosMapped(1:3)=LastPartPos(1:3,iPart)
#endif 
  END IF
  
  ! Determine the subcell coordinates from the position in the reference element
  DO iDir = 1, 3
    IF (PartPosMapped(iDir).GE.1.) THEN
      PartPosMapped(iDir) = 0.99999
    ELSE IF (PartPosMapped(iDir).LE.-1.) THEN
      PartPosMapped(iDir) = -0.99999
    END IF
    CellCoord(iDir) = INT((PartPosMapped(iDir)+1.)/2. * REAL(nOrder+1))
  END DO

  ! Get the subcell ID for the Peano-curve connection
  CALL SubcellID_PeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir, CellCoord, iCell)

  SubPartNum(iCell) = SubPartNum(iCell) + 1
  SubPartIndx(iCell,SubPartNum(iCell)) = iPart
  iPart = PEM%pNext(iPart)
END DO

 ! If the minimum particle number is not given in all subcells, merge subcells together
IF (MINVAL(SubPartNum,MASK = SubPartNum .GT.0).LT.MinPartCell) THEN
  CALL CombineSubcells(iElem)
END IF

! Collision routines for the subcells, considered are only the particles in the respective subcell
DO iCell = 1, nCells
#if (PP_TimeDiscMethod==300)
  CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), SubVolume(iCell))
#else
  CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), SubVolume(iCell))
  ! IF (BGKMovingAverage) THEN
  !   CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), SubVolume(iCell), &
  !   Averaging%SubNode(iLoop)%AverageValues(:)) !! Averaging%SubNode
  ! ELSE
  !   CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), SubVolume(iCell))
  ! END IF
#endif
END DO

DEALLOCATE(SubPartNum)
DEALLOCATE(SubPartIndx) 
DEALLOCATE(SubVolume)

RETURN
END SUBROUTINE SubcellAssingement
  
  
SUBROUTINE CombineSubcells(iElem) 
!============================================================================================================================
!============================================================================================================================
! MODULES                                   
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PEM
USE MOD_Particle_Mesh_Vars     
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: iLoop, nCells, tempCell, iCell, tempCellLast, nAdapt, nSplit, nPart
!-----------------------------------------------------------------------------------------------------------------------------------
nPart = PEM%pNumber(iElem)
nSplit = AdaptMesh(iElem)%SplitOrder     ! 1 -> 8, 2 -> 27, 3 -> 64, ...
nCells = (nSplit+1)**3  

tempCell = 1
tempCellLast = 1

DO iCell = 2, nCells
  IF (SubPartNum(tempCell).LT.MinPartCell) THEN
    ! Move the particle indices
    DO iLoop=1, SubPartNum(iCell) 
      SubPartIndx(tempCell,SubPartNum(tempCell)+iLoop) = SubPartIndx(iCell,iLoop)
    END DO
    ! Move the particle number 
    SubPartNum(tempCell) = SubPartNum(tempCell) + SubPartNum(iCell)
    SubPartNum(iCell) = 0 
    ! Move the subcell volume
    SubVolume(tempCell) = SubVolume(tempCell) + SubVolume(iCell)
    SubVolume(iCell) = 0.
  ELSE
    IF (iCell.EQ.nCells) THEN ! Seperate check for the last subcell
      IF (SubPartNum(iCell).LT.MinPartCell) THEN
        ! Move particle indices
        DO iLoop=1, SubPartNum(iCell) 
          SubPartIndx(tempCell,SubPartNum(tempCell)+iLoop) = SubPartIndx(iCell,iLoop)
        END DO
        ! Move the particle number
        SubPartNum(tempCell) = SubPartNum(tempCell) + SubPartNum(iCell)
        SubPartNum(iCell) = 0
        ! Move the subcell volume
        SubVolume(tempCell) = Subvolume(tempCell) + SubVolume(iCell)
        SubVolume(iCell) = 0.
      END IF
    ELSE
      tempCellLast = tempCell
      tempCell = iCell
    END IF
  END IF
END DO

! Merge remaining fragments which do not meet the minimum particle numbers
IF ((SubPartNum(tempCell).LT.MinPartCell).AND.(tempCell.NE.tempCellLast)) THEN
  ! Move the particle indices
  DO iLoop=1, SubPartNum(tempCell) 
    SubPartIndx(tempCellLast,SubPartNum(tempCellLast)+iLoop) = SubPartIndx(tempCell,iLoop)
  END DO
  ! Move the particle number
  SubPartNum(tempCellLast) = SubPartNum(tempCellLast) + SubPartNum(tempCell)
  SubPartNum(tempCell) = 0 
  ! Move the subcell volume
  SubVolume(tempCellLast) = SubVolume(tempCellLast) + SubVolume(tempCell)
  SubVolume(tempCell) = 0.
END IF

RETURN
END SUBROUTINE CombineSubcells


SUBROUTINE CalcGradients(ElemID,MaxGradient)
!===================================================================================================================================
!> Calculate the density, velocity, pressure and temperature gradients for the subcell adaption process
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: PEM, nSpecies, PartSpecies, Species, usevMPF, PartState
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared, ElemToElemMapping, ElemToElemInfo, ElemMidPoint_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: AdaptMesh
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars              ,ONLY: DSMC, RadialWeighting, VarWeighting
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)              :: MaxGradient
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: CNElemID, CnNbElem, GlobalElemID, GlobNbElem, iVal, LocNbElem, nNbElems
INTEGER                       :: adaptDir, iPart, iLoop, iSpec, iCoord
REAL                          :: RefDens, RefVelo, RefTemp, RefPres, NbDens, NbVelo, NbTemp, NbPres, RefNonAvTemp(3), NbNonAvTemp(3)
REAL                          :: Velo(nSpecies,3), Velo2(nSpecies,3), totalWeight, Temp(3), SpecPartNum(nSpecies)
REAL                          :: Gradient(4), TempGradient(3), MeanDistance, NbDistance
!-----------------------------------------------------------------------------------------------------------------------------------
SpecPartNum = 0.0; Velo = 0.0; Velo2 = 0.0; RefVelo = 0.0; Temp = 0.0; RefTemp = 0.0; MaxGradient  = 0.0; MeanDistance = 0.0
Gradient = 0.0; TempGradient = 0.0; NbVelo = 0.0; NbTemp = 0.0; RefNonAvTemp = 0.0; NbNonAvTemp = 0.0 

GlobalElemID = ElemID+offSetElem
CNElemID     = GetCNElemID(GlobalElemID)

! Determine the number of different species and the total particle number in the element by a loop over all particles 
iPart = PEM%pStart(ElemID)
DO iLoop = 1, PEM%pNumber(ElemID)
  SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

  ! Particle velocities (squared) weighted by the particle number per element
  Velo(PartSpecies(iPart),1:3) = Velo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
  Velo2(PartSpecies(iPart),1:3) = Velo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)
  iPart = PEM%pNext(iPart)
END DO

! Total particle number in the element
totalWeight = SUM(SpecPartNum)

! Number density in the element = Particle number/volume
IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
  RefDens = totalWeight / ElemVolume_Shared(CNElemID)
ELSE
  RefDens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CNElemID)
END IF

! Calculation of the total mean velocity under consideration of all directions
DO iLoop = 1, 3
  IF (totalWeight.GT.0.) THEN
    RefVelo = RefVelo + (SUM(Velo(:,iLoop))/totalWeight)**2
  END IF
END DO
RefVelo = SQRT(RefVelo)

! Calculation of the total mean translational temperature
DO iSpec = 1, nSpecies
  IF (SpecPartNum(iSpec).GE.1.0) THEN 
    Temp  = Species(iSpec)%MassIC/BoltzmannConst * ((Velo2(iSpec,:)/SpecPartNum(iSpec)) - (Velo(iSpec,:)/SpecPartNum(iSpec))**2)
    RefTemp = RefTemp + ((Temp(1) + Temp(2) + Temp(3)) / 3.)*SpecPartNum(iSpec)

    ! Translational temperature in x/y/z direction
    RefNonAvTemp = RefNonAvTemp + Temp*SpecPartNum(iSpec)
  END IF
END DO  
IF (RefTemp.GT.0.) THEN
  IF (totalWeight.GT.0.) THEN
    RefTemp = RefTemp/totalweight
  ELSE
    RefTemp = 0.
  END IF
ELSE 
  RefTemp = 0.
END IF

! Temperature values in x/y/z-direction
IF (totalWeight.GT.0.) THEN
  RefNonAvTemp = RefNonAvTemp / totalWeight
ELSE
  RefNonAvTemp = 0.
END IF

! Pressure calculation (ideal gas law)
RefPres = BoltzmannConst * RefDens * RefTemp

! Relative gradient between the cell and all neighboring cells for the particle density, velocity, translational temperatures
nNbElems = ElemToElemMapping(2,CNElemID)
! Loop over all neighbouring elements
DO iVal = 1, nNbElems 
  CnNbElem = ElemToElemInfo(ElemToElemMapping(1,CNElemID)+iVal)
  GlobNbElem = GetGlobalElemID(CnNbElem)
  LocNbElem = GlobNbElem-offSetElem

  SpecPartNum = 0.0; Velo = 0.0; Velo2 = 0.0; Temp = 0.0; NbVelo = 0.0; NbTemp = 0.0; NbNonAvTemp = 0.0 

  IF ((LocNBElem.LT.1).OR.(LocNBElem.GT.nElems)) CYCLE

  ! Distance between the two midpoints of the cell
  NbDistance = 0.0
  DO iCoord=1, 3
    NbDistance = NbDistance + (ElemMidPoint_Shared(iCoord,CNElemID)-ElemMidPoint_Shared(iCoord,CnNbElem))**2
  END DO
  NbDistance = SQRT(NbDistance)
  MeanDistance = MeanDistance + NbDistance

  ! Determine the number of different species and the total particle number in the element by a loop over all particles 
  iPart = PEM%pStart(LocNbElem)
  DO iLoop = 1, PEM%pNumber(LocNbElem)
    SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)

    ! Particle velocities (squared) weighted by the particle number per element
    Velo(PartSpecies(iPart),1:3) = Velo(PartSpecies(iPart),1:3) + PartState(4:6,iPart) * GetParticleWeight(iPart)
    Velo2(PartSpecies(iPart),1:3) = Velo2(PartSpecies(iPart),1:3) + PartState(4:6,iPart)**2 * GetParticleWeight(iPart)
    iPart = PEM%pNext(iPart)
  END DO

  ! Total particle number in the element
  totalWeight = SUM(SpecPartNum)

  ! Number density in the element = Particle number/volume
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting.OR.VarWeighting%DoVariableWeighting) THEN
    NbDens = totalWeight / ElemVolume_Shared(CnNbElem)
  ELSE
    NbDens = totalWeight * Species(1)%MacroParticleFactor / ElemVolume_Shared(CnNbElem)
  END IF

  ! Calculation of the total mean velocity under consideration of all directions
  DO iLoop = 1, 3
    IF (totalWeight.GT.0.) THEN
      NbVelo = NbVelo + (SUM(Velo(:,iLoop))/totalWeight)**2
    END IF
  END DO
  NbVelo = SQRT(NbVelo)

  ! Calculation of the total mean translational temperature
  DO iSpec = 1, nSpecies
    IF (SpecPartNum(iSpec).GE.1.0) THEN 
      Temp  = Species(iSpec)%MassIC/BoltzmannConst * ((Velo2(iSpec,:)/SpecPartNum(iSpec)) - (Velo(iSpec,:)/SpecPartNum(iSpec))**2)
      NbTemp = NbTemp + ((Temp(1) + Temp(2) + Temp(3)) / 3.)*SpecPartNum(iSpec)

      ! Translational temperature in x/y/z direction
      NbNonAvTemp = NbNonAvTemp + Temp*SpecPartNum(iSpec)
    END IF
  END DO  
  IF (NbTemp.GT.0.) THEN
    IF (totalWeight.GT.0.) THEN
      NbTemp = NbTemp/totalweight
    ELSE 
      NbTemp = 0.
    END IF
  ELSE 
    NbTemp = 0.
  END IF

  ! Temperature values in x/y/z-direction
  IF (totalWeight.GT.0.) THEN 
    NbNonAvTemp = NbNonAvTemp / totalWeight
  ELSE
    NbNonAvTemp = 0.
  END IF

  ! Pressure calculation (ideal gas law)
  NbPres = BoltzmannConst * NbDens * NbTemp

  Gradient(1) = Gradient(1) + ABS((RefDens-NbDens)/NbDistance)
  Gradient(2) = Gradient(2) + ABS((RefVelo-NbVelo)/NbDistance)
  Gradient(3) = Gradient(3) + ABS((RefTemp-NbTemp)/NbDistance)
  Gradient(4) = Gradient(4) + ABS((RefPres-NbPres)/NbDistance)
  TempGradient = TempGradient + ABS((RefNonAvTemp-NbNonAvTemp)/NbDistance) 
END DO ! iNbElem

MeanDistance = MeanDistance/nNbElems

Gradient(1) = Gradient(1)*MeanDistance/(nNbElems*RefDens)
Gradient(2) = Gradient(2)*MeanDistance/(nNbElems*RefVelo)
Gradient(3) = Gradient(3)*MeanDistance/(nNbElems*RefTemp)
Gradient(4) = Gradient(4)*MeanDistance/(nNbElems*RefPres)
TempGradient = TempGradient*MeanDistance/(nNbElems*RefNonAvTemp)

DO iVal = 1, 3
  IF (ISNAN(TempGradient(iVal))) THEN
    TempGradient(iVal) = 0.
  END IF
END DO

! Determine the maximum gradient value
MaxGradient = MAXVAL(Gradient)

IF (ANY(TempGradient(:).EQ.0.)) THEN
  AdaptMesh(ElemID)%RefineDir(1:3) = (/1,2,3/)
ELSE
  adaptDir = MAXLOC(TempGradient,1)
  AdaptMesh(ElemID)%RefineDir(1) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
  TempGradient(adaptDir) = 0.0
  adaptDir = MAXLOC(TempGradient,1)
  AdaptMesh(ElemID)%RefineDir(2) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
  TempGradient(adaptDir) = 0.0 
  adaptDir = MAXLOC(TempGradient,1)
  AdaptMesh(ElemID)%RefineDir(3) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
END IF

END SUBROUTINE CalcGradients


SUBROUTINE DefineElementOrientation(iElem)
!===================================================================================================================================
! Element Orientation between the reference and physical space
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Mesh_Vars          ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
! Output variable: Cell orientation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL         :: NodeCoords(3,8), Vec(3), ResVec1(3,3), ResVec2(3,3), DiffVec(3,3), VecMag(3), MinVec(3)
INTEGER      :: iNode, iVal, jVal, dir1, dir2, angle(3,3), CNElemID
!===================================================================================================================================
CNElemID     = GetCNElemID(iElem+offSetElem) 

DO iNode = 1, 8
  NodeCoords(1:3, iNode) = NodeCoords_Shared(1:3,ElemNodeID_Shared(iNode,CNElemID))
END DO
Vec(1:3) = (/-0.99,0.,0./) ! Border vectors of the reference space [-1,1] ! Determination of Vec and P necessary
ResVec1(1:3,1) = MapToGeo(Vec,NodeCoords) ! Change of the border vectors to 1. and -1. ??
Vec(1:3) = (/0.99,0.,0./)
ResVec2(1:3,1) = MapToGeo(Vec,NodeCoords)
Vec(1:3) = (/0.,-0.99,0./)
ResVec1(1:3,2) = MapToGeo(Vec,NodeCoords)
Vec(1:3) = (/0.,0.99,0./)
ResVec2(1:3,2) = MapToGeo(Vec,NodeCoords)
Vec(1:3) = (/0.,0.,-0.99/)
ResVec1(1:3,3) = MapToGeo(Vec,NodeCoords)
Vec(1:3) = (/0.,0.,0.99/)
ResVec2(1:3,3) = MapToGeo(Vec,NodeCoords)

DO iVal=1,3; DO jVal=1,3
  DiffVec(iVal,jVal) = ResVec1(iVal,jVal)-ResVec2(iVal,jVal) ! Difference vector
END DO; END DO
DO iVal=1,3 ! Magnitude of the three new basis vectors
  VecMag(iVal) = SQRT(DiffVec(1,iVal)*DiffVec(1,iVal)+DiffVec(2,iVal)*DiffVec(2,iVal)+DiffVec(3,iVal)*DiffVec(3,iVal))
END DO
DO iVal=1,3; DO jVal=1,3
  angle(iVal,jVal) = ACOS(DiffVec(iVal,jVal)/ VecMag(jVal)) ! Angle between the basis vectors
  IF (angle(iVal,jVal).GT.PI/2.) THEN 
    angle(iVal,jVal) = PI - angle(iVal,jVal)
  END IF
END DO; END DO

DO iVal=1,3
  MinVec(iVal) = MINVAL(angle(:,iVal)) ! Minimum value of the angle for all three vectors
END DO

dir1 = MINLOC(MinVec,1) ! Vector with the minimum angle
dir2 = MINLOC(angle(:,dir1),1)  !x/y/z component of the minimum angle vector
AdaptMesh(iElem)%CellOrientation(dir2) = dir1 ! Map from x/y/z to the basis vector with the smallest angle
DO iVal=1,3 ! Set the first values to a high number to find the next values in line
  angle(iVal,dir1) = 1000
  angle(dir2,iVal) = 1000 
END DO
DO iVal=1,3 ! Repeat for the next vectors
  MinVec(iVal) = MINVAL(angle(:,iVal))
END DO
dir1 = MINLOC(MinVec,1) 
dir2 = MINLOC(angle(:,dir1),1)  
AdaptMesh(iElem)%CellOrientation(dir2) = dir1
DO iVal=1,3
  angle(iVal,dir1) = 1000
  angle(dir2,iVal) = 1000
END DO
DO iVal=1,3
  MinVec(iVal) = MINVAL(angle(:,iVal))
END DO
dir1 = MINLOC(MinVec,1)
dir2 = MINLOC(angle(:,dir1),1)  
AdaptMesh(iElem)%CellOrientation(dir2) = dir1

END SUBROUTINE DefineElementOrientation

  
FUNCTION MapToGeo(xi,P)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN)          :: xi(3)      ! Basis vectors 
  REAL,INTENT(IN)          :: P(3,8)     ! Coordinates of the element nodes
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                     :: MapToGeo(3)
!===================================================================================================================================

MapToGeo = 0.125  *(P(:,1)*(1-xi(1)) * (1-xi(2)) * (1-xi(3))  &
                  + P(:,2)*(1+xi(1)) * (1-xi(2)) * (1-xi(3))  &
                  + P(:,3)*(1+xi(1)) * (1+xi(2)) * (1-xi(3))  &
                  + P(:,4)*(1-xi(1)) * (1+xi(2)) * (1-xi(3))  &
                  + P(:,5)*(1-xi(1)) * (1-xi(2)) * (1+xi(3))  &
                  + P(:,6)*(1+xi(1)) * (1-xi(2)) * (1+xi(3))  &
                  + P(:,7)*(1+xi(1)) * (1+xi(2)) * (1+xi(3))  &
                  + P(:,8)*(1-xi(1)) * (1+xi(2)) * (1+xi(3))) 

END FUNCTION MapToGeo 


SUBROUTINE BuildPeanoCurve(Order, dir)
!===================================================================================================================================
! Build the Peano curve in dependence of the order of discretization
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Mesh_Vars     ,ONLY: AdaptMesh, PeanoCurve
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES          
  INTEGER, INTENT(IN) :: Order, dir(3)      
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER             :: nCells, iCell
!===================================================================================================================================
nCells = (Order+1)**3

IF (.NOT.ALLOCATED(PeanoCurve)) THEN
  ALLOCATE(PeanoCurve(nCells,3))
ELSE 
  DEALLOCATE(PeanoCurve)
  ALLOCATE(PeanoCurve(nCells,3))
END IF  

DO iCell = 1, nCells
  ! Direction which should be refined last (highest number of subcells)
  PeanoCurve(iCell,dir(1)) = INT((iCell-1)*(Order+1)/nCells) 

  ! Direction which is merged first (lowest number of subcells)
  IF (MOD(INT((iCell-1)/(Order+1)) ,2).EQ.0) THEN ! Increasing values starting from nOrder+1
    PeanoCurve(iCell,dir(3)) = MOD((iCell-1), (Order+1))
  ELSE ! Decreasing values starting from 0
    PeanoCurve(iCell,dir(3)) = Order - MOD((iCell-1),(Order+1))
  END IF

  ! Remaining direction
  IF (MOD(INT((iCell-1)*(Order+1)/nCells) ,2).EQ.0) THEN ! Increasing values in steps of nOrder+1
    PeanoCurve(iCell,dir(2)) = MOD((INT((iCell-1)/(Order+1))), (Order+1)) 
  ELSE ! Decreasing values starting in steps of nOrder+1
    PeanoCurve(iCell,dir(2)) = Order - MOD((INT((iCell-1)/(Order+1))),(Order+1))
  END IF

END DO

END SUBROUTINE BuildPeanoCurve


SUBROUTINE SubcellID_PeanoCurve(Order, dir, Coord, SubCellID)
!===================================================================================================================================
! Determine the subcell-ID based on the Peano curve coordinates
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Particle_Mesh_Vars     ,ONLY: AdaptMesh, PeanoCurve
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES          
  INTEGER, INTENT(IN) :: Order, dir(3), Coord(3)      
  INTEGER, INTENT(OUT):: SubcellID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER             :: nCells, iCell, MinCellVal1, MinCellVal2, iCoord
!===================================================================================================================================
nCells = (Order+1)**3

! First lower bound for the cell ID from the first refinement direction
MinCellVal1 = Coord(dir(1))*nCells/(Order+1) + 1

! Second lower bound for the cell ID from the second refinement direction
IF (MOD(Coord(dir(1)) ,2).EQ.0) THEN ! Increasing order
  MinCellVal2 = MinCellVal1 + (Order+1)*Coord(dir(2))
ELSE ! Decreasing order
  MinCellVal2 = MinCellVal1 + (Order+1)*(Order-Coord(dir(2)))
END IF

! Exact coordinate from the last refinement direction
IF (MOD((MinCellVal2-1)/(Order+1),2).EQ.0) THEN ! Increasing Order
  SubcellID = MinCellVal2 + Coord(dir(3))
ELSE ! Decreasing Order
  SubcellID = MinCellVal2 + Order - Coord(dir(3))
END IF

! Test for correct assingement
DO iCoord =1, 3
  IF (Coord(iCoord).NE.PeanoCurve(SubcellID,iCoord)) THEN
    CALL abort(__STAMP__,'Wrong assingement of the subcells to the Peano curve.')
  END IF
END DO
    
END SUBROUTINE SubcellID_PeanoCurve
  

END MODULE MOD_Cell_Adaption
  