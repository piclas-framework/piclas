!==================================================================================================================================
! Copyright (c) 2021 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
 
MODULE MOD_Particle_SurfChemFlux_Init
!===================================================================================================================================
!> Module for the initialization of particle emission through the chemistry surface flux
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: InitializeParticleSurfChemFlux
!===================================================================================================================================
CONTAINS

SUBROUTINE InitializeParticleSurfChemFlux()
!===================================================================================================================================
! Init Particle Inserting via Surface Flux
!===================================================================================================================================
! Modules
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars              ,ONLY: nBCSides, SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_SurfaceModel_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER               :: iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,currentBC
INTEGER               :: MaxSurfacefluxBCs,nDataBC
INTEGER               :: iReac, SurfNumReac
REAL                  :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
TYPE(tBCdata_auxSCFRadWeight),ALLOCATABLE          :: BCdata_auxSCFTemp(:)
!===================================================================================================================================
! ALLOCATE(SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:nBCSides),SurfMeshSideAreas(1:nBCSides))
! SurfMeshSideAreas=0.
! ! global calculations for sampling the faces for area and vector calculations (checks the integration with CODE_ANALYZE)
! CALL BCSurfMeshSideAreasandNormals()

!MaxSurfacefluxBCs=0
!nDataBC=0 

!read/prepare parameters and determine nec. BCs
!CALL ReadInAndPrepareSurfChemFlux(MaxSurfacefluxBCs, nDataBC)

!CALL CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSCFTemp)

SurfNumReac = SurfChemReac%NumOfReact
!initialize Surfaceflux-specific data
DO iReac=1,SurfNumReac
  DO iSF=1,SurfChemReac%NumOfBounds(iReac)
    currentBC = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC
    IF (BCdata_auxSCF(currentBC)%SideNumber.GT.0) THEN
      
      ! Loop over sides on the surface flux
      DO iSide=1,BCdata_auxSCF(currentBC)%SideNumber
        BCSideID=BCdata_auxSCF(currentBC)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
        SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
        ! Calculate the total area of the surface flux
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            tmp_SubSideAreas(iSample,jSample)=SurfMeshSubSideData(iSample,jSample,BCSideID)%area
          END DO; END DO
        ! Initialize surface flux
        CALL InitSurfFlux(iReac, iSF, iSide, tmp_SubSideAreas, BCdata_auxSCFTemp)
        
      END DO ! iSide
    ELSE IF (BCdata_auxSCF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(__STAMP__&
        ,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
  END DO !iSF
END DO !iSpec

! Deallocate auxiliary variable container (no pointers used inside container)
SDEALLOCATE(BCdata_auxSCFTemp)

END SUBROUTINE InitializeParticleSurfChemFlux


SUBROUTINE ReadInAndPrepareSurfChemFlux(nDataBC)
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_SurfaceModel_Tools     ,ONLY: GetWallTemperature
USE MOD_SurfaceModel_Vars
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT) :: nDataBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSF
INTEGER               :: iReac, SurfNumReac
!===================================================================================================================================
SurfNumReac = SurfChemReac%NumOfReact
ALLOCATE(SurfChemReac%SFMap(SurfChemReac%NumOfReact))
DO iReac=1,SurfNumReac
  IF (SurfChemReac%NumOfBounds(iReac).EQ.0) THEN
    CYCLE
  ELSE
    ALLOCATE(SurfChemReac%SFMap(iReac)%Surfaceflux(1:SurfChemReac%NumOfBounds(iReac)))
    ! Initialize Surfaceflux to BC mapping
    SurfChemReac%SFMap(iReac)%Surfaceflux(:)%BC=-1
    DO iSF=1,SurfChemReac%NumOfBounds(iReac)
      SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC = SurfChemReac%BoundMap(iReac)%Boundaries(iSF)
    END DO
  END IF

  DO iSF=1,SurfChemReac%NumOfBounds(iReac)
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%InsertedParticle = 0
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%InsertedParticleSurplus = 0
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total = 0
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total_allProcsTotal = 0

    ! get surfaceflux data
    IF (SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC.LT.1 .OR. SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC.GT.nPartBound) THEN
      CALL abort(&
__STAMP__&
, 'SurfacefluxBCs must be between 1 and nPartBound!')
    ELSE IF (BCdata_auxSF(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)%SideNumber.EQ. -1) THEN !not set yet
      BCdata_auxSF(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)%SideNumber=0
      nDataBC=nDataBC+1
    END IF

    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%velocityDistribution = 'maxwell_lpn'
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC = 0.
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal = .FALSE.
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%MWTemperatureIC = PartBound%WallTemp(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC)
  
  END DO !iSF
END DO ! iReac

END SUBROUTINE ReadInAndPrepareSurfChemFlux


SUBROUTINE BCSurfMeshSideAreasandNormals()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_SurfaceModel_Vars
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, SideToElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: totalArea
INTEGER                 :: BCSideID, ElemID, iLocSide, SideID, jSample, iSample
REAL                    :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_nOut(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t1(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                    :: tmp_Vec_t2(3,SurfFluxSideSize(1),SurfFluxSideSize(2))
!===================================================================================================================================
totalArea=0.
DO BCSideID=1,nBCSides
  ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
  iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
  IF (SurfFluxSideSize(1).NE.1 .OR. SurfFluxSideSize(2).NE.2) CALL abort(&
__STAMP__&
, 'SurfFluxSideSize must be 1,2 for TriaSurfaceFlux!')
   DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    CALL CalcNormAndTangTriangle(SideID=SideID,nVec=tmp_Vec_nOut(:,iSample,jSample) &
      ,tang1=tmp_Vec_t1(:,iSample,jSample) &
      ,tang2=tmp_Vec_t2(:,iSample,jSample) &
      ,area=tmp_SubSideAreas(iSample,jSample) &
      ,TriNum=jSample)
    SurfMeshSideAreas(BCSideID)=SurfMeshSideAreas(BCSideID)+tmp_SubSideAreas(iSample,jSample)
  END DO; END DO
  totalArea=totalArea+SurfMeshSideAreas(BCSideID)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn=-tmp_Vec_nOut(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1=tmp_Vec_t1(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2=tmp_Vec_t2(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%area=tmp_SubSideAreas(iSample,jSample)
  END DO; END DO
END DO

END SUBROUTINE BCSurfMeshSideAreasandNormals

SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSCFTemp)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSCF is created. Furthermore, the side areas are corrected for the 1D/2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_SurfaceModel_Vars
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, BC, SideToElem
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO, ElemMidPoint_Shared, SideInfo_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Vars          ,ONLY: Symmetry
USE MOD_DSMC_Symmetry          ,ONLY: DSMC_1D_CalcSymmetryArea, DSMC_2D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                           :: nDataBC
TYPE(tBCdata_auxSCFRadWeight), ALLOCATABLE, INTENT(INOUT)        :: BCdata_auxSCFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: TmpMapToBC(1:nDataBC), TmpSideStart(1:nDataBC), TmpSideNumber(1:nDataBC), TmpSideEnd(1:nDataBC)
! PartBC, Start of Linked List for Sides in SurfacefluxBC, Number of Particles in Sides in SurfacefluxBC, End of Linked List for Sides in SurfacefluxBC
INTEGER               :: TmpSideNext(1:nBCSides) !Next: Sides of diff. BCs ar not overlapping!
INTEGER               :: countDataBC,iBC,BCSideID,currentBC,iSF,iCount,iLocSide,SideID
INTEGER               :: ElemID,CNElemID,GlobalElemID
INTEGER               :: iSample,jSample,iSub
INTEGER               :: iReac, SurfNumReac
REAL                  :: ymax,ymin,yMaxTemp,yMinTemp
!===================================================================================================================================
!create Side lists for applicable BCs
!temporary (linked) lists
TmpMapToBC = 0; TmpSideStart = 0; TmpSideNumber = 0; TmpSideEnd = 0; TmpSideNext = 0
countDataBC=0
DO iBC=1,nPartBound
  IF (BCdata_auxSCF(iBC)%SideNumber.EQ. -1) CYCLE !not set for SFs
  countDataBC=countDataBC+1
  TmpMapToBC(countDataBC)=iBC
END DO

DO BCSideID=1,nBCSides
  currentBC=0
  DO iBC=1,countDataBC
    IF (PartBound%MapToPartBC(BC(BCSideID)) .EQ. TmpMapToBC(iBC)) currentBC=iBC
  END DO
  IF (currentBC.EQ.0) CYCLE
  IF (TmpSideNumber(currentBC).EQ.0) THEN
    TmpSideStart(currentBC) = BCSideID ! Start of Linked List for Sides
  ELSE
    TmpSideNext(TmpSideEnd(currentBC)) = BCSideID ! Next Side
  END IF
  !-- prepare for next entry in list
  TmpSideEnd(currentBC) = BCSideID
  TmpSideNumber(currentBC) = TmpSideNumber(currentBC) + 1  ! Number of Sides
END DO ! BCSideID

!save sequential lists in BCdata_auxSCF
DO iBC=1,countDataBC
  BCdata_auxSCF(TmpMapToBC(iBC))%SideNumber=TmpSideNumber(iBC)
  IF (TmpSideNumber(iBC).EQ.0) CYCLE
  ALLOCATE(BCdata_auxSCF(TmpMapToBC(iBC))%SideList(1:TmpSideNumber(iBC)))
  ALLOCATE(BCdata_auxSCF(TmpMapToBC(iBC))%TriaSwapGeo(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
  ALLOCATE(BCdata_auxSCF(TmpMapToBC(iBC))%TriaSideGeo(1:TmpSideNumber(iBC)))

  SurfNumReac = SurfChemReac%NumOfReact
  DO iReac=1,SurfNumReac
    DO iSF=1,SurfChemReac%NumOfBounds(iReac)
      IF (TmpMapToBC(iBC).EQ.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC) THEN !only surfacefluxes with iBC
        ALLOCATE(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
      END IF
    END DO
  END DO
  BCSideID=TmpSideStart(iBC)
  iCount=0

  DO !follow BCSideID list seq. with iCount
    iCount=iCount+1
    BCdata_auxSCF(TmpMapToBC(iBC))%SideList(iCount)=BCSideID
      ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
      iLocSide = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
      SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
      !----- symmetry specific area calculation start
      IF(Symmetry%Order.EQ.2) THEN
        GlobalElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
        CNElemID     = GetCNElemID(GlobalElemID)
        iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
        IF(Symmetry%Axisymmetric) THEN
          ! Calculate the correct area for the axisymmetric (ring area) and 2D (length) and get ymin and ymax for element
          SurfMeshSubSideData(1,1,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,CNElemID, ymin, ymax)
          SurfMeshSubSideData(1,2,BCSideID)%area = 0.0
          ! Determination of the mean radial weighting factor for calculation of the number of particles to be inserted !!! necessary?
          IF (RadialWeighting%DoRadialWeighting) THEN
            IF((ymax - ymin).GT.0.0) THEN
              ! Surfaces that are NOT parallel to the YZ-plane
              IF(RadialWeighting%CellLocalWeighting) THEN
                ! Cell local weighting
                BCdata_auxSCFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = (1. + ElemMidPoint_Shared(2,CNElemID) &
                                                                        / GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))
              ELSE
                BCdata_auxSCFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1.
                DO iSub = 1, RadialWeighting%nSubSides
                  yMinTemp = ymin + (iSub-1) * (ymax - ymin) / RadialWeighting%nSubSides
                  yMaxTemp = ymin + iSub * (ymax - ymin) / RadialWeighting%nSubSides
                  BCdata_auxSCFTemp(TmpMapToBC(iBC))%SubSideWeight(iCount,iSub) = 1.          &
                    + (yMaxTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.) &
                    -  yMinTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.))/(yMaxTemp - yMinTemp)
                END DO
                BCdata_auxSCFTemp(TmpMapToBC(iBC))%SubSideArea(iCount,:) = DSMC_2D_CalcSymmetryAreaSubSides(iLocSide,CNElemID)
              END IF
            ELSE ! surfaces parallel to the x-axis (ymax = ymin)
              BCdata_auxSCFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1. &
                                                                        + ymax/(GEO%ymaxglob)*(RadialWeighting%PartScaleFactor-1.)
            END IF
          END IF
        ELSE
          SurfMeshSubSideData(1,1:2,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,CNElemID) / 2.
        END IF
      ELSE IF(Symmetry%Order.EQ.1) THEN
        SurfMeshSubSideData(1,1:2,BCSideID)%area = DSMC_1D_CalcSymmetryArea(iLocSide,ElemID) / 2.
      END IF

      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        CALL CalcNormAndTangTriangle(SideID=SideID &
          ,midpoint=BCdata_auxSCF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%midpoint &
          ,ndist=BCdata_auxSCF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%ndist &
          ,xyzNod=BCdata_auxSCF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%xyzNod &
          ,Vectors=BCdata_auxSCF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%Vectors &
          ,TriNum=jSample)
      END DO; END DO

    !-- BC-list specific data
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      BCdata_auxSCF(TmpMapToBC(iBC))%LocalArea = BCdata_auxSCF(TmpMapToBC(iBC))%LocalArea &
        + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END DO; END DO

    !-- next Side
    IF (BCSideID .EQ. TmpSideEnd(iBC)) THEN
      IF (TmpSideNumber(iBC).NE.iCount) THEN
        CALL abort(&
__STAMP__&
,'Someting is wrong with TmpSideNumber of iBC',iBC,999.)
      ELSE
        EXIT
      END IF
    END IF
    BCSideID=TmpSideNext(BCSideID)
  END DO ! BCSideID (iCount)
END DO !iBC

END SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux

SUBROUTINE InitSurfFlux(iReac, iSF, iSide, tmp_SubSideAreas, BCdata_auxSCFTemp)
!===================================================================================================================================
!> Initialize surface flux variables in SurfFluxSubSideData type
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, PI
USE MOD_SurfaceModel_Vars
USE MOD_Particle_Vars           ,ONLY: Species,nSpecies
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptive
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iReac, iSF, iSide
REAL, INTENT(IN)      :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
TYPE(tBCdata_auxSCFRadWeight), ALLOCATABLE, INTENT(IN)        :: BCdata_auxSCFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: jSample, iSample, iSub, currentBC, BCSideID
INTEGER               :: iSpec
REAL                  :: vec_nIn(3), nVFR, vec_t1(3), vec_t2(3), projFak, v_thermal, a, vSF
!===================================================================================================================================
currentBC = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%BC
BCSideID=BCdata_auxSCF(currentBC)%SideList(iSide)
DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
  vec_nIn = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn
  vec_t1 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1
  vec_t2 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2
  IF (.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal) THEN
    projFak = DOT_PRODUCT(vec_nIn,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !VeloVecIC projected to inwards normal
  ELSE
    projFak = 1.
  END IF

  DO iSpec=1,nSpecies  
    IF(ANY(SurfChemReac%Products(iReac,:).EQ.iSpec)) THEN
      v_thermal = SQRT(2.*BoltzmannConst*SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) !thermal speed
    ELSE 
      v_thermal = 0.
    END IF
  END DO

  a = 0 !dummy for projected speed ratio in constant v-distri
  !-- compute total volume flow rate through surface
  SELECT CASE(TRIM(SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%velocityDistribution))
  !CASE('constant')
   ! vSF = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak !Velo proj. to inwards normal
   ! nVFR = MAX(tmp_SubSideAreas(iSample,jSample) * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
  CASE('maxwell','maxwell_lpn')
    IF(v_thermal.NE.0.) THEN
      a = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
    ELSE 
      a = 0.
    END IF
    vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
    nVFR = tmp_SubSideAreas(iSample,jSample) * vSF !VFR projected to inwards normal of sub-side
    IF(RadialWeighting%DoRadialWeighting) THEN
      nVFR = nVFR / BCdata_auxSCFTemp(currentBC)%WeightingFactor(iSide)
      DO iSub = 1, RadialWeighting%nSubSides
        IF(ABS(BCdata_auxSCFTemp(currentBC)%SubSideWeight(iSide,iSub)).GT.0.)THEN
          SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) = BCdata_auxSCFTemp(currentBC)%SubSideArea(iSide,iSub) &
                                                             * vSF / BCdata_auxSCFTemp(currentBC)%SubSideWeight(iSide,iSub)
        END IF
      END DO
    END IF
  CASE DEFAULT
    CALL abort(__STAMP__,&
      'ERROR in SurfaceFlux: Wrong velocity distribution!')
  END SELECT
  SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VFR_total + nVFR
  !-- store SF-specific SubSide data in SurfFluxSubSideData (incl. projected velos)
  SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR = nVFR
  SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak = projFak
  SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn = a
  IF (.NOT.SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIsNormal) THEN
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 &
      = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t1,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !v in t1-dir
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 &
      = SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t2,SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%VeloVecIC) !v in t2-dir
  ELSE
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 = 0. !v in t1-dir
    SurfChemReac%SFMap(iReac)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 = 0. !v in t2-dir
  END IF! .NOT.VeloIsNormal
END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)

END SUBROUTINE InitSurfFlux

END MODULE MOD_Particle_SurfChemFlux_Init
