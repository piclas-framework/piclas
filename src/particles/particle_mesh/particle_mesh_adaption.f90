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
PUBLIC :: DefineParametersAdaptMesh, Init_MeshAdaption,MeshAdaption, CalcGradients
!===================================================================================================================================

CONTAINS

SUBROUTINE DefineParametersAdaptMesh()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE

CALL prms%SetSection("Mesh Adaption")
CALL prms%CreateLogicalOption(  'Part-DoSubcellAdaption', 'Anisotropic mesh adaption', '.FALSE.')
CALL prms%CreateIntOption( 'Part-MeshAdapt-MinPartNum', 'Minimum particle number per subcell for the mesh adaption', '6')
CALL prms%CreateIntOption( 'Part-MeshAdapt-MaxPartNum', 'Maximum particle number per subcell for the mesh adaption', '60')
CALL prms%CreateIntOption( 'Part-MeshAdapt-IterationNum', 'Iteration number for the adaption process', '1000')
CALL prms%CreateRealOption('Part-MeshAdapt-RefineLimitGradient', 'Gradient above which the subcell adaption is called', '0.05')

END SUBROUTINE DefineParametersAdaptMesh

SUBROUTINE Init_MeshAdaption()
!===================================================================================================================================
! Read-in of the variables for the subcell adaption
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars
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
! Read-in of the minimum particle number per subcell
MinPartCell = GETINT('Part-MeshAdapt-MinPartNum')
IF (MinPartCell.LT.5) CALL abort(__STAMP__,'ERROR: Given minimum particle number per cell is less than 5!')

MaxPartCell = GETINT('Part-MeshAdapt-MaxPartNum')
IterAdapt   = GETINT('Part-MeshAdapt-IterationNum')

ALLOCATE(AdaptMesh(nElems))
RefineFactorGrad = GETREAL('Part-MeshAdapt-RefineLimitGradient')

! Initialization of the relative orientation of the basis vectors in the physical and reference space
DO iElem=1, nElems
  AdaptMesh(iElem)%CellOrientation(1:3) = (/1,2,3/)
  CALL DefineElementOrientation(iElem)

  AdaptMesh(iElem)%SplitOrder = 0
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
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_Particle_Mesh_Vars
USE MOD_DSMC_Analyze           ,ONLY: CalcMeanFreePath
USE MOD_Particle_Vars          ,ONLY: PEM, WriteMacroVolumeValues, DoVirtualCellMerge, VirtMergedCells, Symmetry
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision
USE MOD_Particle_Vars          ,ONLY: nSpecies, PartSpecies
USE MOD_TimeDisc_Vars
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
INTEGER                  :: nPartMerged, nPartLoc, locElem, iLoopLoc, iMergeElem, CNElemID, nPart, iLoop, iPart
REAL                     :: MaxGradient, SpecPartNum(nSpecies)
LOGICAL                  :: DoMergedCell, DoAdaptCell
!===================================================================================================================================
DoMergedCell = .FALSE.
SpecPartNum = 0.

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

! Skip cell if number of particles is less than 2 or if a cell merge is performed
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
#elif (PP_TimeDiscMethod==400)
  IF (BGKMovingAverage) THEN
    CALL BGK_CollisionOperator(PartIndx, nPartMerged,VirtMergedCells(iElem)%MergedVolume, &
            ElemNodeAveraging(iElem)%Root%AverageValues(:))
  ELSE
    CALL BGK_CollisionOperator(PartIndx, nPartMerged, VirtMergedCells(iElem)%MergedVolume)
 END IF
#else
  CALL PerformPairingAndCollision(PartIndx, nPartMerged, iElem, VirtMergedCells(iElem)%MergedVolume)
#endif
  DEALLOCATE(PartIndx)
END IF

! Update the Subcell-Adaption process only after a certain number of iterations
IF ((MOD(iter,IterAdapt).EQ.0.).OR.(iter.EQ.1)) THEN

  ! Calculate the relative gradient of density, velocity, temperature and pressure to check for an adaption for BGK or FP
  ! Calculate the temperature gradient in all three directions to determine the order of refinement directions
  CALL CalcGradients(iElem,MaxGradient)

  iPart = PEM%pStart(iElem)
  DO iLoop = 1, nPart
    ! Determination of the particle number per species for the calculation of the reference diameter for the mixture
    SpecPartNum(PartSpecies(iPart)) = SpecPartNum(PartSpecies(iPart)) + GetParticleWeight(iPart)
    iPart = PEM%pNext(iPart)
  END DO

  ! Use the mean free path to decide if a subcell adaption is necessary for DSMC
  DSMC%MeanFreePath = CalcMeanFreePath(SpecPartNum, SUM(SpecPartNum), ElemVolume_Shared(CNElemID))

! Check for the refinement criterium based on the timedisc method
#if (PP_TimeDiscMethod==300)
  IF (MaxGradient.GE.RefineFactorGrad) THEN
    DoAdaptCell = .TRUE.
  ELSE
    DoAdaptCell = .FALSE.
  END IF
#elif (PP_TimeDiscMethod==400)
  IF (MaxGradient.GE.RefineFactorGrad) THEN
    DoAdaptCell = .TRUE.
  ELSE
    DoAdaptCell = .FALSE.
  END IF
#else
  IF ((DSMC%MeanFreePath.LT.ElemCharLength_Shared(CNElemID)).OR.(nPart.GT.MaxPartCell)) THEN
    DoAdaptCell = .TRUE.
  ELSE
    DoAdaptCell = .FALSE.
  END IF
#endif

  ! Perform the adaption if the cell has enough particles
  IF ((nPart.GE.(2.*MinPartCell)).AND.DoAdaptCell.AND.(.NOT.DoMergedCell)) THEN
    ! Calculate the splitting order (1-8, 2-27, 3-64...)
    AdaptMesh(iElem)%SplitOrder = INT(nPart/MinPartCell)
    IF (Symmetry%Order.EQ.2 ) THEN
      AdaptMesh(iElem)%SplitOrder = INT(AdaptMesh(iElem)%SplitOrder)**(1./2.)
    ELSE
      AdaptMesh(iElem)%SplitOrder = INT(AdaptMesh(iElem)%SplitOrder)**(1./3.)
    END IF

    IF (Symmetry%Order.EQ.2) THEN
      CALL SubcellAdaption2D(iElem) ! Split the element into the defined number of subcells for the 2D case
    ELSE
      CALL SubcellAdaption(iElem) ! Split the element into the defined number of subcells
    END IF
  ELSE IF (.NOT.DoMergedCell) THEN! Normal routine without merged or split cells

    AdaptMesh(iElem)%SplitOrder = 0
    iPart = PEM%pStart(iElem)
    ALLOCATE(PartIndx(nPart))
    DO iLoop = 1, nPart
      PartIndx(iLoop) = iPart
      iPart = PEM%pNext(iPart)
    END DO

    ! No mesh adaption, direct call to the collision operators
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
#elif (PP_TimeDiscMethod==400)
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(PartIndx, nPart,ElemVolume_Shared(CNElemID), ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
    END IF
#else
    CALL PerformPairingAndCollision(PartIndx, nPart, iElem , ElemVolume_Shared(CNElemID))
#endif
    DEALLOCATE(PartIndx)
  END IF
ELSE
  IF (AdaptMesh(iElem)%SplitOrder.GT.0) THEN
    IF (Symmetry%Order.EQ.2) THEN
      CALL MergedCellAssingement2D(iElem) ! Use the predefined subcell assingement
    ELSE
      CALL MergedCellAssingement(iElem) ! Use the predefined subcell assingement
    END IF
  ELSE
    iPart = PEM%pStart(iElem)
    ALLOCATE(PartIndx(nPart))
    DO iLoop = 1, nPart
      PartIndx(iLoop) = iPart
      iPart = PEM%pNext(iPart)
    END DO

    ! No mesh adaption, direct call to the collision operators
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
#elif (PP_TimeDiscMethod==400)
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(PartIndx, nPart,ElemVolume_Shared(CNElemID), ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(PartIndx, nPart, ElemVolume_Shared(CNElemID))
    END IF
#else
    CALL PerformPairingAndCollision(PartIndx, nPart, iElem , ElemVolume_Shared(CNElemID))
#endif
    DEALLOCATE(PartIndx)
  END IF
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


SUBROUTINE SubcellAdaption(iElem)
!============================================================================================================================
! Split the cell into a certain number of subcells for the 3D case
! Calculate the subcell volume based on the Vandermonde and the Jacobian
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
USE MOD_Mesh_Vars              ,ONLY: offSetElem, sJ
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision
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
INTEGER                       :: iCell, nCells, iPart, nPart, iLoop, iDir, iVal, CellCoord(3), GlobalElemID
INTEGER                       :: j, k, l, nOrder
REAL, ALLOCATABLE             :: Subcell_xGP(:),  Subcell_Vdm(:,:), DetJac(:,:,:,:), LocalVdm(:,:)
REAL                          :: Subcell_wGP, PartPosMapped(3), DetLocal(1,0:PP_N,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! Number of subcells from the split order
nCells = (AdaptMesh(iElem)%SplitOrder+1)**3
MeshAdapt(1,iElem) = REAL(nCells)

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
IF (.NOT.ALLOCATED(AdaptMesh(iElem)%Subvolume)) THEN
  ALLOCATE(AdaptMesh(iElem)%SubVolume(nCells))
ELSE
  DEALLOCATE(AdaptMesh(iElem)%Subvolume)
  ALLOCATE(AdaptMesh(iElem)%Subvolume(nCells))
END IF

! Allocate the merged subcell mapping for each cell
IF (.NOT.ALLOCATED(AdaptMesh(iElem)%SubcellMap)) THEN
  ALLOCATE(AdaptMesh(iElem)%SubcellMap(nCells))
  AdaptMesh(iElem)%SubcellMap = 0
ELSE
  DEALLOCATE(AdaptMesh(iElem)%SubcellMap)
  ALLOCATE(AdaptMesh(iElem)%SubcellMap(nCells))
  AdaptMesh(iElem)%SubcellMap = 0
END IF

nOrder = AdaptMesh(iElem)%SplitOrder
! Allocate: Local determinant and vandermonde
ALLOCATE(DetJac(1,0:nOrder,0:nOrder,0:nOrder))
ALLOCATE(LocalVdm(0:nOrder,0:PP_N))

! Initalize interpolation points and subcell vandermonde
ALLOCATE(Subcell_xGP(0:nOrder))
ALLOCATE(Subcell_Vdm(0:nOrder, 0:nOrder))

! Coordinates of the interpolation points -> midpoints of the subcells
DO iVal = 0, nOrder
  Subcell_xGP(iVal) = (2.*REAL(iVal) + 1.)/(REAL(nOrder+1)) - 1.
END DO
! Weights of the interpolation points, total volume = 8
Subcell_wGP = 2./REAL(1.0+nOrder)

! Local determinant from sJ
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
  AdaptMesh(iElem)%Subvolume(iCell) = DetJac(1,j,k,l)*Subcell_wGP**3
END DO

DEALLOCATE(Subcell_xGP)
DEALLOCATE(Subcell_Vdm)
DEALLOCATE(DetJac)
DEALLOCATE(LocalVdm)

GlobalElemID = iElem + offsetElem
DO iLoop = 1, nPart
  ! Map Particle to -1|1 space
  IF (TrackingMethod.EQ.REFMAPPING) THEN ! Refmapping
    PartPosMapped(1:3) = PartPosRef(1:3,iPart)
  ELSE
#if PP_TimeDiscMethod==300
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosMapped(1:3),GlobalElemID)
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
IF (MINVAL(SubPartNum).LT.MinPartCell) THEN
  CALL CombineSubcells(iElem)
ELSE
  ! Collision routines for the subcells, considered are only the particles in the respective subcell
  DO iCell = 1, nCells
    AdaptMesh(iElem)%SubcellMap(iCell) = iCell
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
#elif (PP_TimeDiscMethod==400)
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell), &
      ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
    END IF
#else
    CALL PerformPairingAndCollision(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), iElem, &
     AdaptMesh(iElem)%SubVolume(iCell))
#endif
  END DO
END IF

DEALLOCATE(SubPartNum)
DEALLOCATE(SubPartIndx)

CALL CalcSubcellNum_X(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Y(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Z(iElem,nOrder,AdaptMesh(iElem)%RefineDir)

RETURN
END SUBROUTINE SubcellAdaption


SUBROUTINE SubcellAdaption2D(iElem)
!============================================================================================================================
! Split the cell into a certain number of subcells for the 2D case
! Calculate the subcell volume based on the Vandermonde and the Jacobian
!============================================================================================================================
! use MODULES
USE MOD_Globals_Vars
USE MOD_BGK_Vars
USE MOD_FPFlow_Vars
USE MOD_Preproc
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars
USE MOD_ChangeBasis            ,ONLY: ChangeBasis2D
USE MOD_DSMC_Vars              ,ONLY: DSMC, SymmetrySide
USE MOD_Particle_Vars          ,ONLY: PEM, PartState, LastPartPos, PartPosRef
USE MOD_Mesh_Vars              ,ONLY: offSetElem, SurfElem, sJ, Face_xGP
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision, MapToGeo2D, GeoCoordToMap2D
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
INTEGER                       :: iCell, nCells, iPart, nPart, iLoop, iDir, iVal, CellCoord(3), GlobalElemID, SideID
INTEGER                       :: j, k, l, nOrder, nCells2D
REAL, ALLOCATABLE             :: Subcell_xGP(:),  Subcell_Vdm(:,:), DetJac(:,:,:), LocalVdm(:,:), CellFace_xGP(:,:,:)
REAL                          :: Subcell_wGP, PartPosMapped(3), DetLocal(1,0:PP_N,0:PP_N), Pos(2), CellFaceLocal(2,0:PP_N,0:PP_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! Number of subcells from the split order
nCells2D = (AdaptMesh(iElem)%SplitOrder+1)**2
nCells = (AdaptMesh(iElem)%SplitOrder+1)**3
MeshAdapt(1,iElem) = REAL(nCells2D)
MeshAdapt(4,iElem) = 1.
! Particle number in the element
nPart = PEM%pNumber(iElem)
iPart = PEM%pStart(iElem)

GlobalElemID = iElem+offSetElem

! Allocate the particle number per subcell
ALLOCATE(SubPartNum(nCells2D))
SubPartNum = 0

! Allocate the particle indices for the subcell
ALLOCATE(SubPartIndx(nCells2D,nPart))
SubPartIndx = 0

! Allocate the subvolumes for each cell
IF (.NOT.ALLOCATED(AdaptMesh(iElem)%Subvolume)) THEN
  ALLOCATE(AdaptMesh(iElem)%SubVolume(nCells2D))
ELSE
  DEALLOCATE(AdaptMesh(iElem)%Subvolume)
  ALLOCATE(AdaptMesh(iElem)%Subvolume(nCells2D))
END IF

! Allocate the merged subcell mapping for each cell
IF (.NOT.ALLOCATED(AdaptMesh(iElem)%SubcellMap)) THEN
  ALLOCATE(AdaptMesh(iElem)%SubcellMap(nCells2D))
  AdaptMesh(iElem)%SubcellMap = 0
ELSE
  DEALLOCATE(AdaptMesh(iElem)%SubcellMap)
  ALLOCATE(AdaptMesh(iElem)%SubcellMap(nCells2D))
  AdaptMesh(iElem)%SubcellMap = 0
END IF

nOrder = AdaptMesh(iElem)%SplitOrder
! Allocate: Local determinant and vandermonde
ALLOCATE(DetJac(1,0:nOrder,0:nOrder))
ALLOCATE(LocalVdm(0:nOrder,0:PP_N))
ALLOCATE(CellFace_xGP(2,0:nOrder,0:nOrder))

! Initalize interpolation points and subcell vandermonde
ALLOCATE(Subcell_xGP(0:nOrder))
ALLOCATE(Subcell_Vdm(0:nOrder, 0:nOrder))

! Coordinates of the interpolation points -> midpoints of the subcells
DO iVal = 0, nOrder
  Subcell_xGP(iVal) = (2.*REAL(iVal) + 1.)/(REAL(nOrder+1)) - 1.
END DO
! Weights of the interpolation points, total volume = 8
Subcell_wGP = 2./REAL(1.0+nOrder)

SideID = SymmetrySide(iElem,1)

DO j=0, PP_N; DO k=0, PP_N
  DetLocal(1,j,k)=SurfElem(j,k,SideID)
END DO; END DO

DO j=0, PP_N; DO k=0, PP_N
  CellFaceLocal(1:2,j,k) = Face_xGP(1:2,j,k,SideID)
END DO; END DO

! Call to the functions
CALL InitializeVandermonde(PP_N,nOrder,wBary,xGP,Subcell_xGP,Subcell_Vdm)
CALL ChangeBasis2D(1,PP_N,nOrder,Subcell_Vdm,DetLocal(:,:,:),DetJac(:,:,:))
CALL ChangeBasis2D(2, PP_N, nOrder, Subcell_Vdm ,CellFaceLocal(:,:,:),CellFace_xGP(:,:,:))

! Build the connection of the subcells in the form of the Peano curve
CALL BuildPeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir)

DO iCell=1, nCells2D
  j = PeanoCurve(iCell,1)
  k = PeanoCurve(iCell,2)
  l = PeanoCurve(iCell,3)

  Pos(1) = Subcell_xGP(j)
  Pos(2) = Subcell_xGP(k)
  Pos(1:2) = MapToGeo2D(Pos, iElem)
  ! Calculate the subcell volume from the weights and the Jacobian
  AdaptMesh(iElem)%Subvolume(iCell) = DetJac(1,j,k)*Subcell_wGP**2 * 2. * Pi * Pos(2)
END DO

DEALLOCATE(Subcell_xGP)
DEALLOCATE(Subcell_Vdm)
DEALLOCATE(DetJac)
DEALLOCATE(LocalVdm)
DEALLOCATE(CellFace_xGP)

GlobalElemID = iElem + offsetElem
DO iLoop = 1, nPart
  ! Map Particle to -1|1 space
  IF (TrackingMethod.EQ.REFMAPPING) THEN ! Refmapping
    PartPosMapped(1:3) = PartPosRef(1:3,iPart)
  ELSE
#if PP_TimeDiscMethod==300
  CALL GeoCoordToMap2D(PartState(1:2,iPart), PartPosMapped(1:2), iElem)
  PartPosMapped(3) = 0.
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
  CellCoord(AdaptMesh(iElem)%RefineDir(1)) = 0
  CALL SubcellID_PeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir, CellCoord, iCell)

  SubPartNum(iCell) = SubPartNum(iCell) + 1
  SubPartIndx(iCell,SubPartNum(iCell)) = iPart
  iPart = PEM%pNext(iPart)
END DO

 ! If the minimum particle number is not given in all subcells, merge subcells together
IF (MINVAL(SubPartNum).LT.MinPartCell) THEN
  CALL CombineSubcells(iElem)
ELSE
  ! Collision routines for the subcells, considered are only the particles in the respective subcell
  DO iCell = 1, nCells2D
    AdaptMesh(iElem)%SubcellMap(iCell) = iCell
#if (PP_TimeDiscMethod==300)
    CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
#elif (PP_TimeDiscMethod==400)
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell), &
      ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
    END IF
#else
    CALL PerformPairingAndCollision(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), iElem, &
    AdaptMesh(iElem)%Subvolume(iCell))
#endif
  END DO
END IF

DEALLOCATE(SubPartNum)
DEALLOCATE(SubPartIndx)

CALL CalcSubcellNum_X(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Y(iElem,nOrder,AdaptMesh(iElem)%RefineDir)

RETURN
END SUBROUTINE SubcellAdaption2D

SUBROUTINE MergedCellAssingement(iElem)
!============================================================================================================================
! Assign the subcells to a previously split-and-merged subcell for the 3D case
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_BGK_Vars
USE MOD_FPFlow_Vars
USE MOD_Preproc
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars
USE MOD_DSMC_Vars              ,ONLY: DSMC
USE MOD_Particle_Vars          ,ONLY: PEM, PartState, LastPartPos, PartPosRef
USE MOD_Mesh_Vars              ,ONLY: offSetElem
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
!-----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: iCell, iCellMerge, nCells, iPart, nPart, iLoop, iDir, CellCoord(3), GlobalElemID
INTEGER                       :: nOrder
REAL                          :: PartPosMapped(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! Number of subcells from the split order
nCells = (AdaptMesh(iElem)%SplitOrder+1)**3
MeshAdapt(1,iElem) = REAL(nCells)

! Particle number in the element
nPart = PEM%pNumber(iElem)
iPart = PEM%pStart(iElem)

! Allocate the particle number per subcell
ALLOCATE(SubPartNum(nCells))
SubPartNum = 0

! Allocate the particle indices for the subcell
ALLOCATE(SubPartIndx(nCells,nPart))
SubPartIndx = 0

nOrder = AdaptMesh(iElem)%SplitOrder

! Build the connection of the subcells in the form of the Peano curve
CALL BuildPeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir)

GlobalElemID = iElem + offsetElem
DO iLoop = 1, nPart
  ! Map Particle to -1|1 space
  IF (TrackingMethod.EQ.REFMAPPING) THEN ! Refmapping
    PartPosMapped(1:3) = PartPosRef(1:3,iPart)
  ELSE
#if PP_TimeDiscMethod==300
    CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosMapped(1:3),GlobalElemID)
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

  ! Get the ID of the merged cells
  iCellMerge = AdaptMesh(iElem)%SubcellMap(iCell)
  IF (iCellMerge.EQ.0) THEN
    CALL abort(__STAMP__,'Wrong assingement in the mapping for the merged subcells.')
  END IF

  SubPartNum(iCellMerge) = SubPartNum(iCellMerge) + 1
  SubPartIndx(iCellMerge,SubPartNum(iCellMerge)) = iPart
  iPart = PEM%pNext(iPart)
END DO

! Call to the collision operator
DO iCell = 1, nCells
  IF (AdaptMesh(iElem)%SubVolume(iCell).GT.0.) THEN
    IF (SubPartNum(iCell).GT.2) THEN
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
#elif (PP_TimeDiscMethod==400)
      IF (BGKMovingAverage) THEN
        CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell), &
        ElemNodeAveraging(iElem)%Root%AverageValues(:))
      ELSE
        CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
      END IF
#else
      CALL PerformPairingAndCollision(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), iElem, AdaptMesh(iElem)%SubVolume(iCell))
#endif
    END IF
  END IF
END DO

DEALLOCATE(SubPartNum)
DEALLOCATE(SubPartIndx)

CALL CalcSubcellNum_X(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Y(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Z(iElem,nOrder,AdaptMesh(iElem)%RefineDir)

RETURN
END SUBROUTINE MergedCellAssingement


SUBROUTINE MergedCellAssingement2D(iElem)
!============================================================================================================================
! Assign the subcells to a previously split-and-merged subcell for the 2D case
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_Globals_Vars
USE MOD_BGK_Vars
USE MOD_FPFlow_Vars
USE MOD_Preproc
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Particle_Mesh_Vars
USE MOD_Particle_Vars          ,ONLY: PEM, PartState, LastPartPos, PartPosRef
USE MOD_Mesh_Vars              ,ONLY: offSetElem, SurfElem
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision, GeoCoordToMap2D
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
!-----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: iCell, nCells, iPart, nPart, iLoop, iDir, CellCoord(3), GlobalElemID, iCellMerge
INTEGER                       :: nOrder, nCells2D
REAL                          :: PartPosMapped(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! Number of subcells from the split order
nCells2D = (AdaptMesh(iElem)%SplitOrder+1)**2
nCells = (AdaptMesh(iElem)%SplitOrder+1)**3
MeshAdapt(1,iElem) = REAL(nCells2D)
MeshAdapt(4,iElem) = 1.

! Particle number in the element
nPart = PEM%pNumber(iElem)
iPart = PEM%pStart(iElem)

GlobalElemID = iElem+offSetElem

! Allocate the particle number per subcell
ALLOCATE(SubPartNum(nCells2D))
SubPartNum = 0

! Allocate the particle indices for the subcell
ALLOCATE(SubPartIndx(nCells2D,nPart))
SubPartIndx = 0

nOrder = AdaptMesh(iElem)%SplitOrder

! Build the connection of the subcells in the form of the Peano curve
CALL BuildPeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir)

GlobalElemID = iElem + offsetElem
DO iLoop = 1, nPart
  ! Map Particle to -1|1 space
  IF (TrackingMethod.EQ.REFMAPPING) THEN ! Refmapping
    PartPosMapped(1:3) = PartPosRef(1:3,iPart)
  ELSE
#if PP_TimeDiscMethod==300
  CALL GeoCoordToMap2D(PartState(1:2,iPart), PartPosMapped(1:2), iElem)
  PartPosMapped(3) = 0.
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
  CellCoord(AdaptMesh(iElem)%RefineDir(1)) = 0
  CALL SubcellID_PeanoCurve(nOrder,AdaptMesh(iElem)%RefineDir, CellCoord, iCell)

  ! Map the subcell to the merged cell
  iCellMerge = AdaptMesh(iElem)%SubcellMap(iCell)
  IF (iCellMerge.EQ.0) THEN
    CALL abort(__STAMP__,'Wrong assingement in the mapping for the merged subcells.')
  END IF

  SubPartNum(iCellMerge) = SubPartNum(iCellMerge) + 1
  SubPartIndx(iCellMerge,SubPartNum(iCellMerge)) = iPart
  iPart = PEM%pNext(iPart)
END DO

! Call to the collision routines
DO iCell = 1, nCells2D
  IF (AdaptMesh(iElem)%SubVolume(iCell).GT.0.) THEN
    IF (SubPartNum(iCell).GT.2) THEN
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
#elif (PP_TimeDiscMethod==400)
      IF (BGKMovingAverage) THEN
        CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell), &
        ElemNodeAveraging(iElem)%Root%AverageValues(:))
      ELSE
        CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%SubVolume(iCell))
      END IF
#else
      CALL PerformPairingAndCollision(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), iElem, AdaptMesh(iElem)%Subvolume(iCell))
#endif
    END IF
  END IF
END DO

DEALLOCATE(SubPartNum)
DEALLOCATE(SubPartIndx)

CALL CalcSubcellNum_X(iElem,nOrder,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Y(iElem,nOrder,AdaptMesh(iElem)%RefineDir)

RETURN
END SUBROUTINE MergedCellAssingement2D


SUBROUTINE CombineSubcells(iElem)
!============================================================================================================================
! Merge the subcells along a Peano curve to reach the minimum particle number per cell
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: PEM, Symmetry
USE MOD_Particle_Mesh_Vars
USE MOD_BGK_CollOperator       ,ONLY: BGK_CollisionOperator
USE MOD_FP_CollOperator        ,ONLY: FP_CollisionOperator
USE MOD_DSMC_ParticlePairing   ,ONLY: PerformPairingAndCollision
USE MOD_BGK_Vars               ,ONLY: BGKMovingAverage,ElemNodeAveraging
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                       :: iLoop, nCells, tempCell, iCell, tempCellLast, nSplit, nPart
!-----------------------------------------------------------------------------------------------------------------------------------
nPart = PEM%pNumber(iElem)
nSplit = AdaptMesh(iElem)%SplitOrder     ! 1 -> 8, 2 -> 27, 3 -> 64, ...
! Calculate the number of subcells for the 2D or the 3D case
IF (Symmetry%Order.EQ.2) THEN
  nCells = (nSplit+1)**2
ELSE
  nCells = (nSplit+1)**3
END IF

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
    AdaptMesh(iElem)%Subvolume(tempCell) = AdaptMesh(iElem)%Subvolume(tempCell) + AdaptMesh(iElem)%Subvolume(iCell)
    AdaptMesh(iElem)%Subvolume(iCell) = 0.
    ! Store the merged values
    AdaptMesh(iElem)%SubcellMap(iCell) = tempCell
    AdaptMesh(iElem)%SubcellMap(tempCell) = tempCell
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
        AdaptMesh(iElem)%Subvolume(tempCell) = AdaptMesh(iElem)%Subvolume(tempCell) + AdaptMesh(iElem)%Subvolume(iCell)
        AdaptMesh(iElem)%Subvolume(iCell) = 0.
        ! Store the merged values
        AdaptMesh(iElem)%SubcellMap(iCell) = tempCell
        AdaptMesh(iElem)%SubcellMap(tempCell) = tempCell
      ELSE
        AdaptMesh(iElem)%SubcellMap(iCell) = iCell
      END IF
    ELSE
      tempCellLast = tempCell
      tempCell = iCell
      AdaptMesh(iElem)%SubcellMap(tempCell) = tempCell
      AdaptMesh(iElem)%SubcellMap(tempCellLast) = tempCellLast
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
  AdaptMesh(iElem)%Subvolume(tempCellLast) = AdaptMesh(iElem)%Subvolume(tempCellLast) + AdaptMesh(iElem)%Subvolume(tempCell)
  AdaptMesh(iElem)%Subvolume(tempCell) = 0.
  ! Store the merged values
  AdaptMesh(iElem)%SubcellMap(tempCell) = tempCellLast
  AdaptMesh(iElem)%SubcellMap(tempCellLast) = tempCellLast
END IF

MeshAdapt(1,iElem) = 0.
DO iCell = 1, nCells
  IF (SubPartNum(iCell).GT.0) THEN
    MeshAdapt(1,iElem) = MeshAdapt(1,iElem) + 1.
    ! Collision routines for the subcells, considered are only the particles in the respective subcell
#if (PP_TimeDiscMethod==300)
      CALL FP_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%Subvolume(iCell))
#elif (PP_TimeDiscMethod==400)
    IF (BGKMovingAverage) THEN
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%Subvolume(iCell), &
      ElemNodeAveraging(iElem)%Root%AverageValues(:))
    ELSE
      CALL BGK_CollisionOperator(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), AdaptMesh(iElem)%Subvolume(iCell))
    END IF
#else
    CALL PerformPairingAndCollision(SubPartIndx(iCell,1:SubPartNum(iCell)), SubPartNum(iCell), iElem, AdaptMesh(iElem)%Subvolume(iCell))
#endif
  END IF
END DO

CALL CalcSubcellNum_X(iElem,nSplit,AdaptMesh(iElem)%RefineDir)
CALL CalcSubcellNum_Y(iElem,nSplit,AdaptMesh(iElem)%RefineDir)
IF (Symmetry%Order.NE.2) THEN
  CALL CalcSubcellNum_Z(iElem,nSplit,AdaptMesh(iElem)%RefineDir)
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
USE MOD_Particle_Vars          ,ONLY: PEM, nSpecies, PartSpecies, Species, usevMPF, PartState, Symmetry
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared, ElemToElemMapping, ElemToElemInfo, ElemMidPoint_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: AdaptMesh
USE MOD_Mesh_Vars              ,ONLY: offsetElem, nElems
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_part_tools             ,ONLY: GetParticleWeight
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting, VarWeighting
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

  ! Gradient weighted by the distance to the neighbour elements
  Gradient(1) = Gradient(1) + ABS((RefDens-NbDens)/NbDistance)
  Gradient(2) = Gradient(2) + ABS((RefVelo-NbVelo)/NbDistance)
  Gradient(3) = Gradient(3) + ABS((RefTemp-NbTemp)/NbDistance)
  Gradient(4) = Gradient(4) + ABS((RefPres-NbPres)/NbDistance)
  TempGradient = TempGradient + ABS((RefNonAvTemp-NbNonAvTemp)/NbDistance)
END DO ! iNbElem

MeanDistance = MeanDistance/nNbElems

IF (RefDens.GT.0.) THEN
  Gradient(1) = Gradient(1)*MeanDistance/(nNbElems*RefDens)
ELSE
  Gradient(1) = 0.
END IF
IF (RefVelo.GT.0.) THEN
  Gradient(2) = Gradient(2)*MeanDistance/(nNbElems*RefVelo)
ELSE
  Gradient(2) = 0.
END IF
IF (RefTemp.GT.0.) THEN
  Gradient(3) = Gradient(3)*MeanDistance/(nNbElems*RefTemp)
ELSE
  Gradient(3) = 0.
END IF
Gradient(4) = Gradient(4)*MeanDistance/(nNbElems*RefPres)
IF (RefPres.GT.0.) THEN
  Gradient(4) = Gradient(4)*MeanDistance/(nNbElems*RefPres)
ELSE
  Gradient(4) = 0.
END IF
DO iVal = 1, 3
  IF (RefNonAvTemp(iVal).GT.0.) THEN
    TempGradient(iVal) = TempGradient(iVal)*MeanDistance/(nNbElems*RefNonAvTemp(iVal))
  ELSE
    TempGradient(iVal) = 0.
  END IF
END DO

DO iVal = 1, 3
  IF (ISNAN(TempGradient(iVal))) THEN
    TempGradient(iVal) = 0.
  END IF
END DO

! Determine the maximum gradient value
MaxGradient = MAXVAL(Gradient)

! Determine the order of the Peano curve merge from the size of the directional temperature gradient
IF (Symmetry%Order.EQ.2) THEN
  IF (ANY(TempGradient(:).EQ.0.)) THEN
    AdaptMesh(ElemID)%RefineDir(1:3) = (/3,2,1/) ! z-direction always last
  ELSE
    adaptdir = 3
    AdaptMesh(ElemID)%RefineDir(1) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
    TempGradient(adaptDir) = -1.
    adaptDir = MAXLOC(TempGradient,1)
    AdaptMesh(ElemID)%RefineDir(2) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
    TempGradient(adaptDir) = -1.
    adaptDir = MAXLOC(TempGradient,1)
    AdaptMesh(ElemID)%RefineDir(3) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
  END IF
ELSE
  IF (ANY(TempGradient(:).EQ.0.)) THEN
    AdaptMesh(ElemID)%RefineDir(1:3) = (/1,2,3/)
  ELSE
    adaptDir = MAXLOC(TempGradient,1)
    AdaptMesh(ElemID)%RefineDir(1) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
    TempGradient(adaptDir) = -1.
    adaptDir = MAXLOC(TempGradient,1)
    AdaptMesh(ElemID)%RefineDir(2) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
    TempGradient(adaptDir) = -1.
    adaptDir = MAXLOC(TempGradient,1)
    AdaptMesh(ElemID)%RefineDir(3) = AdaptMesh(ElemID)%CellOrientation(adaptDir)
  END IF
END IF

END SUBROUTINE CalcGradients


SUBROUTINE DefineElementOrientation(iElem)
!===================================================================================================================================
! Element Orientation between the reference -1|1 space and the physical space
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Mesh_Vars          ,ONLY: offsetElem
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


SUBROUTINE CalcSubcellNum_X(iElem,Order, dir)
!===================================================================================================================================
! Determine the number of subcells in x-direction
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Mesh_Vars     ,ONLY: PeanoCurve, AdaptMesh, MeshAdapt
  USE MOD_Particle_Vars          ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: iElem,Order, dir(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: nCells, TestCellID, TestCoord(3), xVal, yVal, zVal
  INTEGER              :: SubcellID, LastID, z_Order
  INTEGER, ALLOCATABLE :: CellCount(:,:)
!===================================================================================================================================
! TODO : first determine the element orientation
  ! Total and directional subcell number
nCells = (Order+1)**3

ALLOCATE(CellCount(0:Order,0:Order))
CellCount = 0

! Set the cell number in z-direction to 1
IF (Symmetry%Order.EQ.2) THEN
  z_Order = 0
ELSE
  z_Order = Order
END IF

! Loop over the cells in y- and z-direction
DO yVal = 0, Order; DO zVal = 0, z_Order
  TestCoord(AdaptMesh(iElem)%CellOrientation(2)) = yVal
  TestCoord(AdaptMesh(iElem)%CellOrientation(3)) = zVal
  LastID = 0
  ! Loop over the cells in x-direction
  DO xVal = 0, Order
    TestCoord(AdaptMesh(iElem)%CellOrientation(1)) = xVal
    ! Determine the subcell ID for the test coordinate
    CALL SubcellID_PeanoCurve(Order, dir, TestCoord, TestCellID)
    ! Check if the cell is merged
    SubcellID = AdaptMesh(iElem)%SubcellMap(TestCellID)
    IF (SubcellID.NE.LastID) THEN
      ! Set the new cell ID
      LastID = SubcellID
      CellCount(yVal,zVal) = CellCount(yVal,zVal) + 1
    END IF
  END DO
END DO; END DO

MeshAdapt(2,iElem) = SUM(CellCount(:,:))/((Order+1)*(z_Order+1))

DEALLOCATE(CellCount)

END SUBROUTINE CalcSubcellNum_X

SUBROUTINE CalcSubcellNum_Y(iElem,Order, dir)
!===================================================================================================================================
! Determine the number of subcells in y-direction
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Mesh_Vars     ,ONLY: PeanoCurve, AdaptMesh, MeshAdapt
  USE MOD_Particle_Vars          ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: iElem,Order, dir(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: nCells, TestCellID, TestCoord(3), xVal, yVal, zVal
  INTEGER              :: SubcellID, LastID, z_Order
  INTEGER, ALLOCATABLE :: CellCount(:,:)
!===================================================================================================================================
! TODO : first determine the element orientation
  ! Total and directional subcell number
nCells = (Order+1)**3

ALLOCATE(CellCount(0:Order,0:Order))
CellCount = 0

! Set the cell number in z-direction to 1
IF (Symmetry%Order.EQ.2) THEN
  z_Order = 0
ELSE
  z_Order = Order
END IF

! Loop over the cells in x- and z-direction
DO xVal = 0, Order; DO zVal = 0, z_Order
  TestCoord(AdaptMesh(iElem)%CellOrientation(1)) = xVal
  TestCoord(AdaptMesh(iElem)%CellOrientation(3)) = zVal
  LastID = 0
  ! Loop over the cells in y-direction
  DO yVal = 0, Order
    TestCoord(AdaptMesh(iElem)%CellOrientation(2)) = yVal
    ! Determine the subcell ID for the test coordinate
    CALL SubcellID_PeanoCurve(Order, dir, TestCoord, TestCellID)
    ! Check if the cell is merged
    SubcellID = AdaptMesh(iElem)%SubcellMap(TestCellID)
    IF (SubcellID.NE.LastID) THEN
      ! Set the new cell ID
      LastID = SubcellID
      CellCount(xVal,zVal) = CellCount(xVal,zVal) + 1
    END IF
  END DO
END DO; END DO

  MeshAdapt(3,iElem) = SUM(CellCount)/((Order+1)*(z_Order+1))

DEALLOCATE(CellCount)

END SUBROUTINE CalcSubcellNum_Y

SUBROUTINE CalcSubcellNum_Z(iElem,Order, dir)
!===================================================================================================================================
! Determine the number of subcells in z-direction
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Mesh_Vars     ,ONLY: PeanoCurve, AdaptMesh, MeshAdapt
  USE MOD_Particle_Vars          ,ONLY: Symmetry
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN) :: iElem,Order, dir(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER              :: nCells, TestCellID, TestCoord(3), xVal, yVal, zVal
  INTEGER              :: SubcellID, LastID
  INTEGER, ALLOCATABLE :: CellCount(:,:)
!===================================================================================================================================
! TODO : first determine the element orientation
  ! Total and directional subcell number
nCells = (Order+1)**3

ALLOCATE(CellCount(0:Order,0:Order))
CellCount = 0

! Loop over the cells in x- and z-direction
DO xVal = 0, Order; DO yVal = 0, Order
  TestCoord(AdaptMesh(iElem)%CellOrientation(1)) = xVal
  TestCoord(AdaptMesh(iElem)%CellOrientation(2)) = yVal
  LastID = 0
  ! Loop over the cells in y-direction
  DO zVal = 0, Order
    TestCoord(AdaptMesh(iElem)%CellOrientation(3)) = zVal
    ! Determine the subcell ID for the test coordinate
    CALL SubcellID_PeanoCurve(Order, dir, TestCoord, TestCellID)
    ! Check if the cell is merged
    SubcellID = AdaptMesh(iElem)%SubcellMap(TestCellID)
    IF (SubcellID.NE.LastID) THEN
      ! Set the new cell ID
      LastID = SubcellID
      CellCount(xVal,yVal) = CellCount(xVal,yVal) + 1
    END IF
  END DO
END DO; END DO

  MeshAdapt(4,iElem) = SUM(CellCount)/((Order+1)*(Order+1))

DEALLOCATE(CellCount)

END SUBROUTINE CalcSubcellNum_Z

SUBROUTINE BuildPeanoCurve(Order, dir)
!===================================================================================================================================
! Build the Peano curve in dependence of the order of discretization
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Mesh_Vars     ,ONLY: PeanoCurve
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
! Subcell number
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
  USE MOD_Particle_Mesh_Vars     ,ONLY: PeanoCurve
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
  INTEGER             :: nCells, MinCellVal1, MinCellVal2, iCoord
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
    print*, 'Coordinates of the subcell: ', Coord
    print*, 'Coordinates of the subcell generated by the Peano curve: ', PeanoCurve(SubCellID,:)
    CALL abort(__STAMP__,'Wrong assingement of the subcells to the Peano curve.')
  END IF
END DO

END SUBROUTINE SubcellID_PeanoCurve


END MODULE MOD_Cell_Adaption