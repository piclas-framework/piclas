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

MODULE MOD_surface_flux
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

!===================================================================================================================================
PUBLIC         :: InitializeParticleSurfaceflux,ParticleSurfaceflux
!===================================================================================================================================
CONTAINS

SUBROUTINE ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Vars          ,ONLY: AdaptiveWeightFac,nSpecies, Species, VarTimeStep, DoPoissonRounding, DoTimeDepInflow
USE MOD_Particle_Vars          ,ONLY: nMacroRestartFiles, Symmetry2D, UseAdaptive, UseCircularInflow
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, TriaSurfaceFlux
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Mesh_Vars              ,ONLY: NGeo, nElems
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(INOUT) :: MaxSurfacefluxBCs, nDataBC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(42)         :: hilf, hilf2, hilf3
INTEGER               :: iSpec, iElem, iSF
!===================================================================================================================================
AdaptiveWeightFac = GETREAL('Part-AdaptiveWeightingFactor','0.001')

DO iSpec=1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  Species(iSpec)%nSurfacefluxBCs = GETINT('Part-Species'//TRIM(hilf)//'-nSurfacefluxBCs','0')
  IF (useDSMC.AND.(BGGas%NumberOfSpecies.GT.0)) THEN
    IF (BGGas%BackgroundSpecies(iSpec)) THEN
      IF (Species(iSpec)%nSurfacefluxBCs.GT.0) CALL abort(&
  __STAMP__&
  , 'SurfaceFlux is not implemented for the BGG-species!')
    END IF
  END IF
  IF (Species(iSpec)%nSurfacefluxBCs.EQ.0) THEN
    CYCLE
  ELSE
    ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
    ! Initialize Surfaceflux to BC mapping
    Species(iSpec)%Surfaceflux(:)%BC=-1
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      WRITE(UNIT=hilf2,FMT='(I0)') iSF
      hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
      Species(iSpec)%Surfaceflux(iSF)%BC = GETINT('Part-Species'//TRIM(hilf2)//'-BC','0')
    END DO
  END IF

  MaxSurfacefluxBCs=MAX(MaxSurfacefluxBCs,Species(iSpec)%nSurfacefluxBCs)
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    WRITE(UNIT=hilf2,FMT='(I0)') iSF
    hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = 0
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total = 0
    Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal = 0
    ! get surfaceflux data
    IF (Species(iSpec)%Surfaceflux(iSF)%BC.LT.1 .OR. Species(iSpec)%Surfaceflux(iSF)%BC.GT.nPartBound) THEN
      CALL abort(&
__STAMP__&
, 'SurfacefluxBCs must be between 1 and nPartBound!')
    ELSE IF (BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber.EQ. -1) THEN !not set yet
      BCdata_auxSF(Species(iSpec)%Surfaceflux(iSF)%BC)%SideNumber=0
      nDataBC=nDataBC+1
    END IF
    Species(iSpec)%Surfaceflux(iSF)%velocityDistribution  = &
        TRIM(GETSTR('Part-Species'//TRIM(hilf2)//'-velocityDistribution','constant'))
    IF (TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'constant' .AND. &
        TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'liquid' .AND. &
        TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'maxwell' .AND. &
        TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'maxwell_lpn') THEN
      CALL abort(&
__STAMP__&
, 'Only constant or maxwell-like velodistri implemented for surfaceflux!')
    END IF
    Species(iSpec)%Surfaceflux(iSF)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC')
    Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-VeloIsNormal')
    IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow=.FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%VeloVecIC          =GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3)
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow=GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CircularInflow')
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        UseCircularInflow=.TRUE.
        Species(iSpec)%Surfaceflux(iSF)%dir(1)       = GETINT('Part-Species'//TRIM(hilf2)//'-axialDir')
        IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.1) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=2
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=3
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.2) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=3
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=1
        ELSE IF (Species(iSpec)%Surfaceflux(iSF)%dir(1).EQ.3) THEN
          Species(iSpec)%Surfaceflux(iSF)%dir(2)=1
          Species(iSpec)%Surfaceflux(iSF)%dir(3)=2
        ELSE
          CALL abort(__STAMP__&
            ,'ERROR in init: axialDir for SFradial must be between 1 and 3!')
        END IF
        IF ( Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(2)).NE.0. .OR. &
             Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(Species(iSpec)%Surfaceflux(iSF)%dir(3)).NE.0. ) THEN
          CALL abort(__STAMP__&
            ,'ERROR in init: axialDir for SFradial do not correspond to VeloVecIC!')
        END IF
        Species(iSpec)%Surfaceflux(iSF)%origin       = GETREALARRAY('Part-Species'//TRIM(hilf2)//'-origin',2,'0. , 0.')
        WRITE(UNIT=hilf3,FMT='(E16.8)') HUGE(Species(iSpec)%Surfaceflux(iSF)%rmax)
        Species(iSpec)%Surfaceflux(iSF)%rmax     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmax',TRIM(hilf3))
        Species(iSpec)%Surfaceflux(iSF)%rmin     = GETREAL('Part-Species'//TRIM(hilf2)//'-rmin','0.')
      END IF
    END IF !.NOT.VeloIsNormal
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      !--- normalize VeloVecIC
      IF (nMacroRestartFiles.EQ.0) THEN
        IF (.NOT. ALL(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(:).eq.0.)) THEN
          Species(iSpec)%Surfaceflux(iSF)%VeloVecIC = Species(iSpec)%Surfaceflux(iSF)%VeloVecIC &
            /SQRT(DOT_PRODUCT(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC))
        END IF
      END IF
    END IF
    Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC       = GETREAL('Part-Species'//TRIM(hilf2)//'-MWTemperatureIC','0.')
    Species(iSpec)%Surfaceflux(iSF)%PartDensity           = GETREAL('Part-Species'//TRIM(hilf2)//'-PartDensity','0.')
    Species(iSpec)%Surfaceflux(iSF)%ReduceNoise           = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-ReduceNoise','.FALSE.')
    IF (DoPoissonRounding .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(*,*)'WARNING: Poisson sampling not possible for noise reduction of surfacefluxes:'
      SWRITE(*,*)'switching now to Random rounding...'
      DoPoissonRounding   = .FALSE.
    END IF
    IF (DoTimeDepInflow .AND. Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      SWRITE(*,*)'WARNING: Time-dependent inflow is not possible for noise reduction of surfacefluxes:'
      SWRITE(*,*)'switching now to Random rounding...'
      DoTimeDepInflow   = .FALSE.
    END IF
    IF (TriaSurfaceFlux) THEN
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject=.FALSE.
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%AcceptReject          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-AcceptReject','.TRUE.')
    END IF
    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.GT.1) THEN
      SWRITE(*,*)'WARNING: BezierSampleN > 0 may not be necessary as ARM is used for SurfaceFlux!'
    ELSE IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%AcceptReject .AND. BezierSampleN.LE.NGeo .AND. .NOT.TriaSurfaceFlux) THEN
      SWRITE(*,*)'WARNING: The choosen small BezierSampleN (def.: NGeo) might result in inhom. SurfFluxes without ARM!'
    END IF
    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
      WRITE( hilf3, '(I0.2)') NGeo*NGeo*NGeo !1 for linear elements, this is an arbitray estimation for higher N!
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = GETINT('Part-Species'//TRIM(hilf2)//'-ARM_DmaxSampleN',hilf3)
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = 0
    END IF
    ! ================================= ADAPTIVE BC READ IN START =================================================================!
    Species(iSpec)%Surfaceflux(iSF)%Adaptive         = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-Adaptive','.FALSE.')
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      DoPoissonRounding = .TRUE.
      UseAdaptive  = .TRUE.
      IF(DoRefMapping) THEN
        CALL abort(__STAMP__&
            ,'ERROR: Adaptive surface flux boundary conditions are not implemented with DoRefMapping!')
      END IF
      IF(Symmetry2D.OR.VarTimeStep%UseVariableTimeStep) THEN
        CALL abort(__STAMP__&
            ,'ERROR: Adaptive surface flux boundary conditions are not implemented with 2D/axisymmetric or variable time step!')
      END IF
      ! Total area of surface flux
      IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = Pi*(Species(iSpec)%Surfaceflux(iSF)%rmax &
          *Species(iSpec)%Surfaceflux(iSF)%rmax - Species(iSpec)%Surfaceflux(iSF)%rmin*Species(iSpec)%Surfaceflux(iSF)%rmin)
      ELSE
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = 0.
      END IF
      Species(iSpec)%Surfaceflux(iSF)%AdaptiveType         = GETINT('Part-Species'//TRIM(hilf2)//'-Adaptive-Type')
      SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
      ! Farbar2014 - Case 1: Inlet Type 1, constant pressure and temperature
      !              Case 2: Outlet Type 1, constant pressure
      !              Case 3: Inlet Type 2, constant mass flow and temperature
      ! Lei2017    - Case 4: cf_3, constant mass flow and temperature N through mass flow and particles out
      CASE(1,2)
        Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure  = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Pressure')
        Species(iSpec)%Surfaceflux(iSF)%PartDensity       = Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure &
                                                            / (BoltzmannConst * Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC)
      CASE(3,4)
        Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow     = GETREAL('Part-Species'//TRIM(hilf2)//'-Adaptive-Massflow')
        IF(Species(iSpec)%Surfaceflux(iSF)%VeloIC.LE.0.0) THEN
          CALL abort(__STAMP__&
            ,'ERROR in init of adaptive inlet: positive initial guess of velocity for Type 3/Type 4 condition required!')
        END IF
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1:3,1:nElems))
        ! Initializing the array with the given velocity vector and magnitude. It is used as a fallback, when the sampled velocity
        ! in the cell is zero (e.g. when starting a simulation with zero particles)
        DO iElem = 1, nElems
          Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1:3,iElem) = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
                                                                                * Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(1:3)
        END DO
        Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut = 0
      END SELECT
      IF (PartBound%TargetBoundCond(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.PartBound%ReflectiveBC) THEN ! iSF on reflective BC
        IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow.AND.(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.NE.4)) THEN
          CALL abort(__STAMP__&
            ,'ERROR in adaptive surface flux: using a reflective BC without circularInflow is only allowed for Type 4!')
        END IF
      END IF
    ELSE
      Species(iSpec)%Surfaceflux(iSF)%AdaptiveType = 0
    END IF
    ! ================================= ADAPTIVE BC READ IN END ===================================================================!
  END DO !iSF
END DO ! iSpec

END SUBROUTINE ReadInAndPrepareSurfaceFlux


SUBROUTINE BCSurfMeshSideAreasandNormals()
!===================================================================================================================================
! Initialize the variables first
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, SurfMeshSubSideData, BezierSampleN, SurfMeshSideAreas, TriaSurfaceFlux
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, SideToElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas, CalcNormAndTangTriangle
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
  ElemID = SideToElem(1,BCSideID)
  IF (ElemID.LT.1) THEN !not sure if necessary
    ElemID   = SideToElem(2,BCSideID)
    iLocSide = SideToElem(4,BCSideID)
  ELSE
    iLocSide = SideToElem(3,BCSideID)
  END IF
  SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
  IF (TriaSurfaceFlux) THEN
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
  ELSE
    IF (ANY(SurfFluxSideSize.NE.BezierSampleN)) CALL abort(&
__STAMP__&
, 'SurfFluxSideSize must be BezierSampleN,BezierSampleN for .NOT.TriaSurfaceFlux!')
    CALL GetBezierSampledAreas(SideID=SideID &
      ,BezierSampleN=BezierSampleN &
      ,SurfMeshSubSideAreas=tmp_SubSideAreas &
      ,SurfMeshSideArea_opt=SurfMeshSideAreas(BCSideID) &
      ,SurfMeshSubSideVec_nOut_opt=tmp_Vec_nOut &
      ,SurfMeshSubSideVec_t1_opt=tmp_Vec_t1 &
      ,SurfMeshSubSideVec_t2_opt=tmp_Vec_t2)
  END IF
  totalArea=totalArea+SurfMeshSideAreas(BCSideID)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn=-tmp_Vec_nOut(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1=tmp_Vec_t1(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2=tmp_Vec_t2(:,iSample,jSample)
    SurfMeshSubSideData(iSample,jSample,BCSideID)%area=tmp_SubSideAreas(iSample,jSample)
  END DO; END DO
END DO
#ifdef CODE_ANALYZE
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
IPWRITE(*,*)"totalArea       = ",totalArea
IPWRITE(*,*)"totalArea/(pi) = ",totalArea/(ACOS(-1.))
IPWRITE(*,*)" ===== TOTAL AREA (all BCsides) ====="
#endif /*CODE_ANALYZE*/
END SUBROUTINE BCSurfMeshSideAreasandNormals

SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSFTemp)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, SurfMeshSubSideData, SurfFluxSideSize, TriaSurfaceFlux, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SideType
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
USE MOD_Mesh_Vars              ,ONLY: nBCSides, offsetElem, BC, SideToElem
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO, ElemMidPoint_Shared, SideInfo_Shared
USE MOD_Particle_Vars          ,ONLY: UseCircularInflow, Species, DoSurfaceFlux, nSpecies, Symmetry2D, Symmetry2DAxisymmetric
USE MOD_DSMC_Symmetry2D        ,ONLY: DSMC_2D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_Particle_Surfaces      ,ONLY: CalcNormAndTangTriangle
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                           :: nDataBC
TYPE(tBCdata_auxSFRadWeight), ALLOCATABLE, INTENT(INOUT)        :: BCdata_auxSFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: TmpMapToBC(1:nDataBC), TmpSideStart(1:nDataBC), TmpSideNumber(1:nDataBC), TmpSideEnd(1:nDataBC)
! PartBC, Start of Linked List for Sides in SurfacefluxBC, Number of Particles in Sides in SurfacefluxBC, End of Linked List for Sides in SurfacefluxBC
INTEGER               :: TmpSideNext(1:nBCSides) !Next: Sides of diff. BCs ar not overlapping!
INTEGER               :: countDataBC, iBC, BCSideID, currentBC, iSF, ElemID, iCount, iLocSide, SideID, iPartBound
INTEGER               :: iSample, jSample, iSpec, iSub
REAL, ALLOCATABLE     :: areasLoc(:),areasGlob(:)
LOGICAL               :: OutputSurfaceFluxLinked
REAL                  :: ymax, ymin, yMaxTemp, yMinTemp
!===================================================================================================================================
!-- 2.: create Side lists for applicable BCs
!--- 2a: temporary (linked) lists
OutputSurfaceFluxLinked=GETLOGICAL('OutputSurfaceFluxLinked','.FALSE.')
TmpMapToBC = 0; TmpSideStart = 0; TmpSideNumber = 0; TmpSideEnd = 0; TmpSideNext = 0
countDataBC=0
DO iBC=1,nPartBound
  IF (BCdata_auxSF(iBC)%SideNumber.EQ. -1) CYCLE !not set for SFs or CollectCharges
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

IF(RadialWeighting%DoRadialWeighting) THEN
  ALLOCATE(BCdata_auxSFTemp(1:nPartBound))
END IF

!--- 2b: save sequential lists in BCdata_auxSF
DO iBC=1,countDataBC
  BCdata_auxSF(TmpMapToBC(iBC))%SideNumber=TmpSideNumber(iBC)
  IF (TmpSideNumber(iBC).EQ.0) CYCLE
  ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%SideList(1:TmpSideNumber(iBC)))
  IF (TriaSurfaceFlux) THEN
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
    ALLOCATE(BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(1:TmpSideNumber(iBC)))
  END IF
  IF(RadialWeighting%DoRadialWeighting) THEN
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(1:TmpSideNumber(iBC)))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor = 1.
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight = 0.
    ALLOCATE(BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
    BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea = 0.
  END IF
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      IF (TmpMapToBC(iBC).EQ.Species(iSpec)%Surfaceflux(iSF)%BC) THEN !only surfacefluxes with iBC
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)))
        IF (UseCircularInflow .AND. (iSF .LE. Species(iSpec)%nSurfacefluxBCs)) THEN
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(1:TmpSideNumber(iBC)) )
        END IF
        IF(RadialWeighting%DoRadialWeighting) THEN
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%nVFRSub(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
        END IF
      END IF
    END DO
  END DO
  BCSideID=TmpSideStart(iBC)
  iCount=0
  DO !follow BCSideID list seq. with iCount
    iCount=iCount+1
    BCdata_auxSF(TmpMapToBC(iBC))%SideList(iCount)=BCSideID
    IF (TriaSurfaceFlux) THEN
      ElemID = SideToElem(1,BCSideID)
      IF (ElemID.LT.1) THEN !not sure if necessary
        ElemID = SideToElem(2,BCSideID)
        iLocSide = SideToElem(4,BCSideID)
      ELSE
        iLocSide = SideToElem(3,BCSideID)
      END IF
      SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
      !----- symmetry specific area calculation start
      IF(Symmetry2D) THEN
        ElemID = SideInfo_Shared(SIDE_ELEMID,SideID)
        iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
        IF(Symmetry2DAxisymmetric) THEN
          ! Calculate the correct area for the axisymmetric (ring area) and 2D (length) and get ymin and ymax for element
          SurfMeshSubSideData(1,1,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,ElemID, ymin, ymax)
          SurfMeshSubSideData(1,2,BCSideID)%area = 0.0
          ! Determination of the mean radial weighting factor for calculation of the number of particles to be inserted
          IF (RadialWeighting%DoRadialWeighting) THEN
            IF((ymax - ymin).GT.0.0) THEN
              ! Surfaces that are NOT parallel to the YZ-plane
              IF(RadialWeighting%CellLocalWeighting) THEN
                ! Cell local weighting
                BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = (1. + ElemMidPoint_Shared(2,ElemID+offsetElem)&
                                                                        / GEO%ymaxglob*(RadialWeighting%PartScaleFactor-1.))
              ELSE
                BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1.
                DO iSub = 1, RadialWeighting%nSubSides
                  yMinTemp = ymin + (iSub-1) * (ymax - ymin) / RadialWeighting%nSubSides
                  yMaxTemp = ymin + iSub * (ymax - ymin) / RadialWeighting%nSubSides
                  BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideWeight(iCount,iSub) = 1.          &
                    + (yMaxTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.) &
                    -  yMinTemp**2/(GEO%ymaxglob*2.)*(RadialWeighting%PartScaleFactor-1.))/(yMaxTemp - yMinTemp)
                END DO
                BCdata_auxSFTemp(TmpMapToBC(iBC))%SubSideArea(iCount,:) = DSMC_2D_CalcSymmetryAreaSubSides(iLocSide,ElemID)
              END IF
            ELSE ! surfaces parallel to the x-axis (ymax = ymin)
              BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = 1. &
                                                                        + ymax/(GEO%ymaxglob)*(RadialWeighting%PartScaleFactor-1.)
            END IF
          END IF
        ELSE
          SurfMeshSubSideData(1,1:2,BCSideID)%area = DSMC_2D_CalcSymmetryArea(iLocSide,ElemID) / 2.
        END IF
      END IF
      !----- symmetry specific area calculation end
      IF (.NOT.TriaTracking) THEN !check that all sides are planar if TriaSurfaceFlux is used for tracing or refmapping
        IF (SideType(SideID).NE.PLANAR_RECT .AND. SideType(SideID).NE.PLANAR_NONRECT) CALL abort(&
__STAMP__&
,'every surfaceflux-sides must be planar if TriaSurfaceFlux is used for tracing or refmapping!!!')
      END IF !.NOT.TriaTracking

      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        CALL CalcNormAndTangTriangle(SideID=SideID &
          ,midpoint=BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%midpoint &
          ,ndist=BCdata_auxSF(TmpMapToBC(iBC))%TriaSwapGeo(iSample,jSample,iCount)%ndist &
          ,xyzNod=BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%xyzNod &
          ,Vectors=BCdata_auxSF(TmpMapToBC(iBC))%TriaSideGeo(iCount)%Vectors &
          ,TriNum=jSample)
      END DO; END DO
    END IF !TriaSurfaceFlux

    !-- BC-list specific data
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      BCdata_auxSF(TmpMapToBC(iBC))%LocalArea = BCdata_auxSF(TmpMapToBC(iBC))%LocalArea &
        + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END DO; END DO

    !-- next Side
    IF (BCSideID .EQ. TmpSideEnd(iBC)) THEN
      IF (TmpSideNumber(iBC).NE.iCount) THEN
        CALL abort(&
__STAMP__&
,'Someting is wrong with TmpSideNumber of iBC',iBC,999.)
      ELSE
        IF(OutputSurfaceFluxLinked)THEN
          IPWRITE(*,'(I4,I7,A53,I0)') iCount,' Sides have been found for Surfaceflux-linked PartBC ',TmpMapToBC(iBC)
        END IF
        DoSurfaceFlux=.TRUE.
        EXIT
      END IF
    END IF
    BCSideID=TmpSideNext(BCSideID)
  END DO ! BCSideID (iCount)
END DO !iBC
!-- communicate areas
#if USE_MPI
   ALLOCATE( areasLoc(1:nPartBound) , areasGlob(1:nPartBound) )
   areasLoc=0.
   areasGlob=0.
   DO iPartBound=1,nPartBound
     areasLoc(iPartBound)=BCdata_auxSF(iPartBound)%LocalArea
   END DO
   CALL MPI_ALLREDUCE(areasLoc,areasGlob,nPartBound,MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
#endif
   DO iPartBound=1,nPartBound
#if USE_MPI
     BCdata_auxSF(iPartBound)%GlobalArea=areasGlob(iPartBound)
#else
     BCdata_auxSF(iPartBound)%GlobalArea=BCdata_auxSF(iPartBound)%LocalArea
#endif
!     IPWRITE(*,'(I4,A,I4,2(x,E16.8))') 'areas:-', &
!       iPartBound,BCdata_auxSF(iPartBound)%GlobalArea,BCdata_auxSF(iPartBound)%LocalArea
   END DO
#if USE_MPI
   DEALLOCATE(areasLoc,areasGlob)
#endif

END SUBROUTINE CreateSideListAndFinalizeAreasSurfFlux

SUBROUTINE AllocateAdaptiveBCSampling(AdaptiveInitDone)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars               ,ONLY: offsetElem, nElems
USE MOD_Particle_Vars           ,ONLY: nSpecies, Adaptive_MacroVal
USE MOD_Restart_Vars            ,ONLY: DoRestart,RestartFile
USE MOD_HDF5_INPUT              ,ONLY: ReadArray, DatasetExists
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(OUT)              :: AdaptiveInitDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: AdaptiveDataExists
REAL,ALLOCATABLE                  :: ElemData_HDF5(:,:,:)
INTEGER                           :: iElem
!===================================================================================================================================
ALLOCATE(Adaptive_MacroVal(1:13,1:nElems,1:nSpecies))
Adaptive_MacroVal(:,:,:)=0
! If restart is done, check if adptiveinfo exists in state, read it in and write to adaptive_macrovalues
AdaptiveInitDone = .FALSE.
IF (DoRestart) THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  ! read local ParticleInfo from HDF5
  CALL DatasetExists(File_ID,'AdaptiveInfo',AdaptiveDataExists)
  IF(AdaptiveDataExists)THEN
    AdaptiveInitDone = .TRUE.
    ALLOCATE(ElemData_HDF5(1:10,1:nSpecies,1:nElems))
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nSpecies   => INT(nSpecies,IK) ,&
          offsetElem => INT(offsetElem,IK),&
          nElems     => INT(nElems,IK)    )
      CALL ReadArray('AdaptiveInfo',3,(/10_IK, nSpecies, nElems/),offsetElem,3,RealArray=ElemData_HDF5(:,:,:))
    END ASSOCIATE
    DO iElem = 1,nElems
      Adaptive_MacroVal(DSMC_VELOX,iElem,:)   = ElemData_HDF5(1,:,iElem)
      Adaptive_MacroVal(DSMC_VELOY,iElem,:)   = ElemData_HDF5(2,:,iElem)
      Adaptive_MacroVal(DSMC_VELOZ,iElem,:)   = ElemData_HDF5(3,:,iElem)
      Adaptive_MacroVal(DSMC_TEMPX,iElem,:)   = ElemData_HDF5(4,:,iElem)
      Adaptive_MacroVal(DSMC_TEMPY,iElem,:)   = ElemData_HDF5(5,:,iElem)
      Adaptive_MacroVal(DSMC_TEMPZ,iElem,:)   = ElemData_HDF5(6,:,iElem)
      Adaptive_MacroVal(DSMC_NUMDENS,iElem,:) = ElemData_HDF5(7,:,iElem)
      ! Porous BC parameter (11: Pumping capacity [m3/s], 12: Static pressure [Pa], 13: Integral pressure difference [Pa])
      Adaptive_MacroVal(11:13,iElem,:)        = ElemData_HDF5(8:10,:,iElem)
    END DO
    SDEALLOCATE(ElemData_HDF5)
  END IF
  CALL CloseDataFile()
END IF

END SUBROUTINE AllocateAdaptiveBCSampling

SUBROUTINE DefineCircInflowRejectType(iSpec, iSF, iSide)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars                 ,ONLY: offsetElem, SideToElem
USE MOD_Particle_Surfaces         ,ONLY: GetSideBoundingBox
USE MOD_Particle_Surfaces_Vars    ,ONLY: BCdata_auxSF
USE MOD_Particle_Vars             ,ONLY: Species
USE MOD_Particle_Mesh_Tools       ,ONLY: GetGlobalNonUniqueSideID, GetSideBoundingBoxTria
USE MOD_Particle_Tracking_Vars    ,ONLY: TriaTracking
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, iSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: BoundingBox(1:3,1:8), origin(2), Vector1(3), Vector2(3), Vector3(3), xyzNod(3), VecBoundingBox(3)
REAL                  :: corner(3), corners(2,4), radiusCorner(2,4), rmax, rmin, point(2), vec(2)
LOGICAL               :: r0inside, intersecExists(2,2)
INTEGER               :: dir(3), iNode, currentBC, BCSideID, ElemID, iLocSide, SideID
!===================================================================================================================================
!-- check where the sides are located relative to rmax (based on corner nodes of bounding box)
!- RejectType=0 : complete side is inside valid bounds
!- RejectType=1 : complete side is outside of valid bounds
!- RejectType=2 : side is partly inside valid bounds
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
ElemID = SideToElem(1,BCSideID)
IF (ElemID.LT.1) THEN !not sure if necessary
  ElemID = SideToElem(2,BCSideID)
  iLocSide = SideToElem(4,BCSideID)
ELSE
  iLocSide = SideToElem(3,BCSideID)
END IF
SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
IF  (TriaTracking) THEN
  CALL GetSideBoundingBoxTria(SideID,BoundingBox)
ELSE
  CALL GetSideBoundingBox(SideID,BoundingBox)
END IF
intersecExists=.FALSE.
!atan2Shift=0.
r0inside=.FALSE.
dir=Species(iSpec)%Surfaceflux(iSF)%dir
origin=Species(iSpec)%Surfaceflux(iSF)%origin
Vector1(:)=0.
Vector2(:)=0.
Vector3(:)=0.
xyzNod(1)=MINVAL(BoundingBox(1,:))
xyzNod(2)=MINVAL(BoundingBox(2,:))
xyzNod(3)=MINVAL(BoundingBox(3,:))
VecBoundingBox(1) = MAXVAL(BoundingBox(1,:)) -MINVAL(BoundingBox(1,:))
VecBoundingBox(2) = MAXVAL(BoundingBox(2,:)) -MINVAL(BoundingBox(2,:))
VecBoundingBox(3) = MAXVAL(BoundingBox(3,:)) -MINVAL(BoundingBox(3,:))
Vector1(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(2)) = VecBoundingBox(dir(2))
Vector2(dir(3)) = VecBoundingBox(dir(3))
Vector3(dir(3)) = VecBoundingBox(dir(3))

!-- determine rmax (and corners)
DO iNode=1,4
  SELECT CASE(iNode)
  CASE(1)
    corner = xyzNod
  CASE(2)
    corner = xyzNod + Vector1
  CASE(3)
    corner = xyzNod + Vector2
  CASE(4)
    corner = xyzNod + Vector3
  END SELECT
  corner(dir(2)) = corner(dir(2)) - origin(1)
  corner(dir(3)) = corner(dir(3)) - origin(2)
  corners(1:2,iNode)=(/corner(dir(2)),corner(dir(3))/) !coordinates of orth. dirs
  radiusCorner(1,iNode)=SQRT(corner(dir(2))**2+corner(dir(3))**2)
END DO !iNode
rmax=MAXVAL(radiusCorner(1,1:4))

!-- determine rmin
DO iNode=1,4
  SELECT CASE(iNode)
  CASE(1)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
    vec=(/Vector1(dir(2)),Vector1(dir(3))/)
  CASE(2)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)-origin
    vec=(/Vector3(dir(2)),Vector3(dir(3))/)
  CASE(3)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
    vec=(/-Vector1(dir(2)),-Vector1(dir(3))/)
  CASE(4)
    point=(/xyzNod(dir(2)),xyzNod(dir(3))/)+(/Vector2(dir(2)),Vector2(dir(3))/)-origin
    vec=(/-Vector3(dir(2)),-Vector3(dir(3))/)
  END SELECT
  vec=point + MIN(MAX(-DOT_PRODUCT(point,vec)/DOT_PRODUCT(vec,vec),0.),1.)*vec
  radiusCorner(2,iNode)=SQRT(DOT_PRODUCT(vec,vec)) !rmin
END DO !iNode

!-- determine if r0 is inside of bounding box
IF ((origin(1) .GE. MINVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(2),:))) .AND. &
   (origin(1) .LE. MAXVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(2),:))) .AND. &
   (origin(2) .GE. MINVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(3),:))) .AND. &
   (origin(2) .LE. MAXVAL(BoundingBox(Species(iSpec)%Surfaceflux(iSF)%dir(3),:))) ) THEN
   r0inside = .TRUE.
END IF
IF (r0inside) THEN
  rmin = 0.
ELSE
  rmin=MINVAL(radiusCorner(2,1:4))
END IF
! define rejecttype
IF ( (rmin .GT. Species(iSpec)%Surfaceflux(iSF)%rmax) .OR. (rmax .LT. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=1
#ifdef CODE_ANALYZE
  CountCircInflowType(2,iSF,iSpec)=CountCircInflowType(2,iSF,iSpec)+1
#endif
ELSE IF ( (rmax .LE. Species(iSpec)%Surfaceflux(iSF)%rmax) .AND. (rmin .GE. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=0
#ifdef CODE_ANALYZE
  CountCircInflowType(1,iSF,iSpec)=CountCircInflowType(1,iSF,iSpec)+1
#endif
ELSE
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=2
#ifdef CODE_ANALYZE
  CountCircInflowType(3,iSF,iSpec)=CountCircInflowType(3,iSF,iSpec)+1
#endif
END IF !  (rmin > Surfaceflux-rmax) .OR. (rmax < Surfaceflux-rmin)
IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
  IF((Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).NE.1)   &
      .AND.(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)) CALL CircularInflow_Area(iSpec,iSF,iSide,BCSideID)
END IF
END SUBROUTINE DefineCircInflowRejectType

SUBROUTINE InitNonAdaptiveSurfFlux(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, PI
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, SurfMeshSubSideData, tBCdata_auxSFRadWeight, BCdata_auxSF
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, iSide
REAL, INTENT(IN)      :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
TYPE(tBCdata_auxSFRadWeight), ALLOCATABLE, INTENT(IN)        :: BCdata_auxSFTemp(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: jSample, iSample, iSub, currentBC, BCSideID
REAL                  :: vec_nIn(3), nVFR, vec_t1(3), vec_t2(3), projFak, v_thermal, a, vSF
!===================================================================================================================================
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
  vec_nIn = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn
  vec_t1 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1
  vec_t2 = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    projFak = DOT_PRODUCT(vec_nIn,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !VeloVecIC projected to inwards normal
  ELSE
    projFak = 1.
  END IF
  v_thermal = SQRT(2.*BoltzmannConst*Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) !thermal speed
  a = 0 !dummy for projected speed ratio in constant v-distri
  !-- compute total volume flow rate through surface
  SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
  CASE('constant')
    vSF = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak !Velo proj. to inwards normal
    nVFR = MAX(tmp_SubSideAreas(iSample,jSample) * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
  CASE('liquid')
    vSF = v_thermal / (2.0*SQRT(PI)) !mean flux velocity through normal sub-face
    nVFR = tmp_SubSideAreas(iSample,jSample) * vSF !VFR projected to inwards normal of sub-side
  CASE('maxwell','maxwell_lpn')
    IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
      CALL abort(&
__STAMP__&
,'Something is wrong with the Surfaceflux parameters!')
    END IF
    a = Species(iSpec)%Surfaceflux(iSF)%VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
    vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
    nVFR = tmp_SubSideAreas(iSample,jSample) * vSF !VFR projected to inwards normal of sub-side
    IF(RadialWeighting%DoRadialWeighting) THEN
      nVFR = nVFR / BCdata_auxSFTemp(currentBC)%WeightingFactor(iSide)
      DO iSub = 1, RadialWeighting%nSubSides
        IF(ABS(BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)).GT.0.)THEN
          Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub) = BCdata_auxSFTemp(currentBC)%SubSideArea(iSide,iSub) &
                                                             * vSF / BCdata_auxSFTemp(currentBC)%SubSideWeight(iSide,iSub)
        END IF
      END DO
    END IF
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'wrong velo-distri for Surfaceflux!')
  END SELECT
  IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN !check rmax-rejection
    IF (Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) THEN ! complete side is outside of valid bounds
      nVFR = 0.
    END IF
  END IF
  Species(iSpec)%Surfaceflux(iSF)%VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total + nVFR
  !-- store SF-specific SubSide data in SurfFluxSubSideData (incl. projected velos)
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR = nVFR
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak = projFak
  Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn = a
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 &
      = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t1,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !v in t1-dir
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 &
      = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
      * DOT_PRODUCT(vec_t2,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC) !v in t2-dir
  ELSE
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1 = 0. !v in t1-dir
    Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2 = 0. !v in t2-dir
  END IF! .NOT.VeloIsNormal
END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)

END SUBROUTINE InitNonAdaptiveSurfFlux

SUBROUTINE InitAdaptiveSurfFlux(iSpec, iSF, ElemID, AdaptiveInitDone)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species, Adaptive_MacroVal
USE MOD_Particle_Boundary_Vars ,ONLY: nPorousBC
USE MOD_Restart_Vars           ,ONLY: DoMacroscopicRestart, MacroRestartValues
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF, ElemID
LOGICAL, INTENT(IN)   :: AdaptiveInitDone
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF (.NOT.AdaptiveInitDone) THEN
  ! initialize velocity, trans_temperature and density of macrovalues
  IF (DoMacroscopicRestart) THEN
    Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec) = MacroRestartValues(ElemID,iSpec,DSMC_VELOX)
    Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec) = MacroRestartValues(ElemID,iSpec,DSMC_VELOY)
    Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec) = MacroRestartValues(ElemID,iSpec,DSMC_VELOZ)
    Adaptive_MacroVal(DSMC_TEMPX,ElemID,iSpec) = MAX(0.,MacroRestartValues(ElemID,iSpec,DSMC_TEMPX))
    Adaptive_MacroVal(DSMC_TEMPY,ElemID,iSpec) = MAX(0.,MacroRestartValues(ElemID,iSpec,DSMC_TEMPY))
    Adaptive_MacroVal(DSMC_TEMPZ,ElemID,iSpec) = MAX(0.,MacroRestartValues(ElemID,iSpec,DSMC_TEMPZ))
    Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec) = MacroRestartValues(ElemID,iSpec,DSMC_NUMDENS)
    IF(nPorousBC.GT.0) THEN
      CALL abort(&
__STAMP__&
,'Macroscopic restart with porous BC and without state file including adaptive BC info not implemented!')
    END IF
  ELSE
    Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
        * Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(1)
    Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
        * Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(2)
    Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%VeloIC &
        * Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(3)
    Adaptive_MacroVal(DSMC_TEMPX,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
    Adaptive_MacroVal(DSMC_TEMPY,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
    Adaptive_MacroVal(DSMC_TEMPZ,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
    Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%PartDensity
  END IF
END IF
END SUBROUTINE InitAdaptiveSurfFlux

SUBROUTINE InitReduceNoiseSF(iSpec, iSF)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: iSpec, iSF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: iProc
!===================================================================================================================================
IF(MPIroot)THEN
  ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(0:nProcessors-1))
  Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs=0.
ELSE
  ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(1)) !dummy for debug
END IF !MPIroot
#if USE_MPI
CALL MPI_GATHER(Species(iSpec)%Surfaceflux(iSF)%VFR_total,1,MPI_DOUBLE_PRECISION &
  ,Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs,1,MPI_DOUBLE_PRECISION,0,PartMPI%COMM,iError)
IF(MPIroot)THEN
  DO iProc=0,PartMPI%nProcs-1
    Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal &
      + Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(iProc)
  END DO
END IF
#else  /*USE_MPI*/
Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs=Species(iSpec)%Surfaceflux(iSF)%VFR_total
Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal=Species(iSpec)%Surfaceflux(iSF)%VFR_total
#endif  /*USE_MPI*/
END SUBROUTINE InitReduceNoiseSF

SUBROUTINE InitializeParticleSurfaceflux()
!===================================================================================================================================
! Init Particle Inserting via Surface Flux
!===================================================================================================================================
! Modules
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Mesh_Vars              ,ONLY: nBCSides, SideToElem, offsetElem, NGeo
USE MOD_Particle_Boundary_Vars ,ONLY: nPorousBC
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas, tBCdata_auxSFRadWeight
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux, DoPoissonRounding, nDataBC_CollectCharges, DoTimeDepInflow
USE MOD_Particle_Vars          ,ONLY: UseAdaptive, UseCircularInflow, DoForceFreeSurfaceFlux
USE MOD_Restart_Vars           ,ONLY: DoRestart, RestartTime
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Vars          ,ONLY: CountCircInflowType
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER               :: iSpec,iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,currentBC
INTEGER               :: iCopy1, iCopy2, iCopy3, MaxSurfacefluxBCs,nDataBC
REAL                  :: tmp_SubSideDmax(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: tmp_BezierControlPoints2D(2,0:NGeo,0:NGeo,SurfFluxSideSize(1),SurfFluxSideSize(2))
REAL                  :: VFR_total
LOGICAL               :: AdaptiveInitDone
TYPE(tBCdata_auxSFRadWeight),ALLOCATABLE          :: BCdata_auxSFTemp(:)
#if USE_MPI
REAL                  :: totalAreaSF_global
#endif
!===================================================================================================================================
ALLOCATE(SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:nBCSides),SurfMeshSideAreas(1:nBCSides))
SurfMeshSideAreas=0.
! global calculations for sampling the faces for area and vector calculations (checks the integration with CODE_ANALYZE)
CALL BCSurfMeshSideAreasandNormals()

UseCircularInflow=.FALSE.
UseAdaptive=.FALSE.
MaxSurfacefluxBCs=0
nDataBC=nDataBC_CollectCharges !sides may be also used for collectcharges of floating potential!!!
DoSurfaceFlux=.FALSE.

!-- 1.: read/prepare parameters and determine nec. BCs
CALL ReadInAndPrepareSurfaceFlux(MaxSurfacefluxBCs, nDataBC)

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoPoissonRounding,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoTimeDepInflow,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
#endif /*USE_MPI*/

CALL CreateSideListAndFinalizeAreasSurfFlux(nDataBC, BCdata_auxSFTemp)

#ifdef CODE_ANALYZE
IF (UseCircularInflow) THEN
  ALLOCATE(CountCircInflowType(1:3,1:MaxSurfacefluxBCs,1:nSpecies))
  CountCircInflowType = 0
END IF
#endif

!-- 3.: initialize Surfaceflux-specific data
! Allocate sampling of near adaptive boundary element values
IF(UseAdaptive.OR.(nPorousBC.GT.0)) CALL AllocateAdaptiveBCSampling(AdaptiveInitDone)

DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    !--- 3a: SF-specific data of Sides
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC !go through sides if present in proc...
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) THEN
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2), &
                  1:BCdata_auxSF(currentBC)%SideNumber))
        Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight = 0.0
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2), &
                1:BCdata_auxSF(currentBC)%SideNumber))
        Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide = 0.0
      END IF
    END IF
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
      DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
        BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
        ElemID = SideToElem(1,BCSideID)
        IF (ElemID.LT.1) THEN !not sure if necessary
          ElemID = SideToElem(2,BCSideID)
          iLocSide = SideToElem(4,BCSideID)
        ELSE
          iLocSide = SideToElem(3,BCSideID)
        END IF
        SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
        IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
          CALL GetBezierSampledAreas(SideID=SideID &
            ,BezierSampleN=BezierSampleN &
            ,BezierSurfFluxProjection_opt=.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal &
            ,SurfMeshSubSideAreas=tmp_SubSideAreas &  !SubSide-areas proj. to inwards normals
            ,DmaxSampleN_opt=Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN &
            ,Dmax_opt=tmp_SubSideDmax &
            ,BezierControlPoints2D_opt=tmp_BezierControlPoints2D)
        ELSE IF (.NOT.TriaSurfaceFlux) THEN
          CALL GetBezierSampledAreas(SideID=SideID &
            ,BezierSampleN=BezierSampleN &
            ,BezierSurfFluxProjection_opt=.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal &
            ,SurfMeshSubSideAreas=tmp_SubSideAreas)  !SubSide-areas proj. to inwards normals
        ELSE !TriaSurfaceFlux
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            tmp_SubSideAreas(iSample,jSample)=SurfMeshSubSideData(iSample,jSample,BCSideID)%area
            IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
              IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
                Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = Species(iSpec)%Surfaceflux(iSF)%totalAreaSF &
                                                              + SurfMeshSubSideData(iSample,jSample,BCSideID)%area
              END IF
            END IF
          END DO; END DO
        END IF
        !Circular Inflow Init
        IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) CALL DefineCircInflowRejectType(iSpec, iSF, iSide)
        !Init non-adaptive SF
        IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
          CALL InitNonAdaptiveSurfFlux(iSpec, iSF, iSide, tmp_SubSideAreas, BCdata_auxSFTemp)
        END IF
        !Init Stuff for acceptance-rejection sampling on SF
        IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax = tmp_SubSideDmax(iSample,jSample)
            IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
              ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                          ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo))
              DO iCopy1=0,NGeo; DO iCopy2=0,NGeo; DO iCopy3=1,2
                Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
                                                                   ,iSide)%BezierControlPoints2D(iCopy3,iCopy2,iCopy1) &
                  = tmp_BezierControlPoints2D(iCopy3,iCopy2,iCopy1,iSample,jSample)
              END DO; END DO; END DO
            END IF !.NOT.VeloIsNormal
          END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
        END IF
        !Init adaptive SF
        IF (Species(iSpec)%Surfaceflux(iSF)%Adaptive) CALL InitAdaptiveSurfFlux(iSpec, iSF, ElemID, AdaptiveInitDone)
      END DO ! iSide

    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
    !--- 3b: ReduceNoise initialization
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) CALL InitReduceNoiseSF(iSpec, iSF)
    !--- Finalize adaptive SF
#if USE_MPI
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        totalAreaSF_global = 0.0
        CALL MPI_ALLREDUCE(Species(iSpec)%Surfaceflux(iSF)%totalAreaSF,totalAreaSF_global,1, &
                            MPI_DOUBLE_PRECISION,MPI_SUM,PartMPI%COMM,IERROR)
        Species(iSpec)%Surfaceflux(iSF)%totalAreaSF = totalAreaSF_global
      END IF
    END IF
#endif

#ifdef CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(x,I0),A,3(x,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ' &
                                            , CountCircInflowType(1,iSF,iSpec),CountCircInflowType(2, iSF,iSpec) &
                                            , CountCircInflowType(3, iSF,iSpec)
    END IF
#endif /*CODE_ANALYZE*/
  END DO !iSF
END DO !iSpec

#ifdef CODE_ANALYZE
SDEALLOCATE(CountCircInflowType)
#endif
! Deallocate auxiliary variable container (no pointers used inside container)
SDEALLOCATE(BCdata_auxSFTemp)

! Setting variables required after a restart
IF(DoRestart) THEN
  DO iSpec=1,nSpecies
    DO iSF = 1, Species(iSpec)%NumberOfInits
      Species(iSpec)%Init(iSF)%InsertedParticle = INT(Species(iSpec)%Init(iSF)%ParticleEmission * RestartTime,8)
    END DO
    DO iSF = 1, Species(iSpec)%nSurfacefluxBCs
      IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal !proc global total (for non-root: dummy!!!)
      ELSE
        VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total               !proc local total
      END IF
      Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity * RestartTime &
        / Species(iSpec)%MacroParticleFactor * VFR_total,8)
    END DO
  END DO
END IF

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoSurfaceFlux,1,MPI_LOGICAL,MPI_LOR,PartMPI%COMM,iError) !set T if at least 1 proc have SFs
#endif  /*USE_MPI*/
IF (.NOT.DoSurfaceFlux) THEN !-- no SFs defined
  SWRITE(*,*) 'WARNING: No Sides for SurfacefluxBCs found! DoSurfaceFlux is now disabled!'
END IF
DoForceFreeSurfaceFlux = GETLOGICAL('DoForceFreeSurfaceFlux','.FALSE.')

END SUBROUTINE InitializeParticleSurfaceflux

SUBROUTINE CalcPartInsSubSidesStandardCase(iSpec, iSF, PartInsSubSides)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_TimeDisc_Vars           ,ONLY: dt, RKdtFrac, RKdtFracTotal, Time
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, BCdata_auxSF
USE MOD_Part_Emission_Tools     ,ONLY: IntegerDivide, SamplePoissonDistri
#if USE_MPI
USE MOD_Particle_MPI_Vars       ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iSpec, iSF
INTEGER, INTENT(OUT), ALLOCATABLE   :: PartInsSubSides(:,:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)        :: inserted_Particle_iter,inserted_Particle_time,inserted_Particle_diff
INTEGER                :: currentBC, PartInsSF, IntSample
REAL                   :: VFR_total, PartIns, RandVal1
INTEGER, ALLOCATABLE   :: PartInsProc(:)
!===================================================================================================================================
  !--- Noise reduction (both ReduceNoise=T (with comm.) and F (proc local), but not for DoPoissonRounding)
  currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
  IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
    !-- calc global to-be-inserted number of parts and distribute to procs (root)
    ALLOCATE(PartInsProc(0:nProcessors-1))
    PartInsProc=0
  END IF !ReduceNoise
  IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%ReduceNoise .OR. MPIroot) THEN !ReduceNoise: root only
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
      VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcsTotal !proc global total
    ELSE
      VFR_total = Species(iSpec)%Surfaceflux(iSF)%VFR_total               !proc local total
    END IF
    PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
      * dt*RKdtFrac * VFR_total
    inserted_Particle_iter = INT(PartIns,8)
    PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
      * (Time + dt*RKdtFracTotal) * VFR_total
    !-- random-round the inserted_Particle_time for preventing periodicity
    IF (inserted_Particle_iter.GE.1) THEN
      CALL RANDOM_NUMBER(RandVal1)
      inserted_Particle_time = INT(PartIns+RandVal1,8)
    ELSE IF (inserted_Particle_iter.GE.0) THEN !needed, since InsertedParticleSurplus can increase
                                               !and _iter>1 needs to be possible for preventing periodicity
      IF (ALMOSTEQUAL(PartIns,0.)) THEN !dummy for procs without SFs (needed for mpi-comm, are cycled later)
        inserted_Particle_time = INT(PartIns,8)
      ELSE !poisson-distri of PartIns-INT(PartIns)
        CALL SamplePoissonDistri( PartIns-INT(PartIns) , IntSample )
        inserted_Particle_time = INT(INT(PartIns)+IntSample,8) !INT(PartIns) + POISDISTRI( PartIns-INT(PartIns) )
      END IF
    ELSE !dummy for procs without SFs (needed for mpi-comm, are cycled later)
      inserted_Particle_time = INT(PartIns,8)
    END IF
    !-- evaluate inserted_Particle_time and inserted_Particle_iter
    inserted_Particle_diff = inserted_Particle_time - Species(iSpec)%Surfaceflux(iSF)%InsertedParticle &
      - inserted_Particle_iter - Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
    PartInsSF = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
    Species(iSpec)%Surfaceflux(iSF)%InsertedParticle = Species(iSpec)%Surfaceflux(iSF)%InsertedParticle + INT(PartInsSF,8)
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
#if USE_MPI
      CALL IntegerDivide(PartInsSF,nProcessors,Species(iSpec)%Surfaceflux(iSF)%VFR_total_allProcs(0:nProcessors-1) &
        ,PartInsProc(0:nProcessors-1))
#else  /*USE_MPI*/
      PartInsProc=PartInsSF
#endif  /*USE_MPI*/
    END IF !ReduceNoise
  END IF !ReduceNoise, MPIroot
#if USE_MPI
  IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN !scatter PartInsProc into PartInsSF of procs
    CALL MPI_SCATTER(PartInsProc(0:nProcessors-1),1,MPI_INTEGER,PartInsSF,1,MPI_INTEGER,0,PartMPI%COMM,IERROR)
  END IF !ReduceNoise
#endif  /*USE_MPI*/
  !-- calc global to-be-inserted number of parts and distribute to SubSides (proc local)
  SDEALLOCATE(PartInsSubSides)
  ALLOCATE(PartInsSubSides(SurfFluxSideSize(1),SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber))
  PartInsSubSides=0
  IF (BCdata_auxSF(currentBC)%SideNumber.LT.1) THEN
    IF (PartInsSF.NE.0) CALL abort(&
  __STAMP__&
  ,'ERROR in ParticleSurfaceflux: Someting is wrong with PartInsSF of BC ',currentBC)
  ELSE
    CALL IntegerDivide(PartInsSF,BCdata_auxSF(currentBC)%SideNumber*SurfFluxSideSize(1)*SurfFluxSideSize(2) &
      ,Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2) &
                                                          ,1:BCdata_auxSF(currentBC)%SideNumber)%nVFR &
      ,PartInsSubSides(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber) )
  END IF
END SUBROUTINE CalcPartInsSubSidesStandardCase

SUBROUTINE CalcExtraPartsReactiveBC(iSpec, iSF)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: iSpec, iSF
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
      CALL abort(&
__STAMP__&
,'Reactive Boundaries not implemented in this PICLas Version!')
!  IF (SurfMesh%SideIDToSurfID(SideID).GT.0) THEN
!    ! sumEvapPart for triatracking is only allocated over (1,1,nsurfsides,nspecies)
!    IF (.NOT.TriaSurfaceFlux .OR. (iSample.EQ.1 .AND. jSample.EQ.1)) THEN
!      ExtraParts = SurfModel%SumEvapPart(iSample,jSample,SurfMesh%SideIDToSurfID(SideID),iSpec)
!      SurfModel%SumEvapPart(iSample,jSample,SurfMesh%SideIDToSurfID(SideID),iSpec) = 0
!    END IF
!    IF (TriaSurfaceFlux) THEN
!      IF (iSample.EQ.1 .AND. jSample.EQ.1) THEN !first tria
!        AreasTria(1)=SurfMeshSubSideData(1,1,BCSideID)%area
!        AreasTria(2)=SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),BCSideID)%area
!        ExtraPartsTria(:) = 0
!        CALL IntegerDivide(ExtraParts, 2, AreasTria, ExtraPartsTria)
!        ExtraParts = ExtraPartsTria(1)
!      ELSE !second tria
!        ExtraParts = ExtraPartsTria(2)
!      END IF
!    END IF !TriaSurfaceFlux
!    SurfModel%Info(iSpec)%NumOfDes=SurfModel%Info(iSpec)%NumOfDes+ExtraParts
!  END IF !SurfMesh%SideIDToSurfID(SideID).GT.0
END SUBROUTINE CalcExtraPartsReactiveBC

SUBROUTINE DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars        ,ONLY: NodeCoords_Shared, ElemSideNodeID_Shared, SideInfo_Shared
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                 :: SideID
REAL, INTENT(IN)                    :: xyzNod(3)
REAL, INTENT(OUT)                   :: minPos(2), RVec(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iLocSide, globElemId, Node1, Node2, minVec
REAL                    :: Vector1(3), Vector2(3), Vector2D(2)
!===================================================================================================================================
iLocSide = SideInfo_Shared(SIDE_LOCALID,SideID)
globElemId = SideInfo_Shared(SIDE_ELEMID,SideID)
!-- compute parallelogram of triangle (only simple 2 value adds/subs, other from init)
Node1 = 2     ! normal = cross product of 1-2 and 1-3 for first triangle
Node2 = 4     !          and 1-3 and 1-4 for second triangle
Vector1(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node1,iLocSide,globElemId)+1) - xyzNod(1:3)
Vector2(1:3) = NodeCoords_Shared(1:3,ElemSideNodeID_Shared(Node2,iLocSide,globElemId)+1) - xyzNod(1:3)
IF (ABS(Vector1(3)).GT.ABS(Vector2(3))) THEN
  Vector2D(1:2) = Vector2(1:2)
ELSE
  Vector2D(1:2) = Vector1(1:2)
END IF
minVec = MINLOC((/xyzNod(2), xyzNod(2)+Vector2D(2)/),1)
SELECT CASE(minVec)
CASE(1)
  minPos(1:2) = xyzNod(1:2)
  RVec(1:2) =  Vector2D(1:2)
CASE(2)
  minPos(1:2) = xyzNod(1:2) + Vector2D(1:2)
  RVec(1:2) = - Vector2D(1:2)
END SELECT
END SUBROUTINE DefineSideDirectVec2D

SUBROUTINE CreateLinkedListReactiveBC(iSpec, iSF, SideID, iSide, iSample, jSample, ParticleIndexNbr, currentSurfFluxPart)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars                     ,ONLY: DSMC
USE MOD_TimeDisc_Vars                 ,ONLY: Time, TEnd
USE MOD_Particle_Boundary_Vars        ,ONLY: SurfMesh, PartBound, nSurfSample
USE MOD_Particle_Tracking_Vars        ,ONLY: TriaTracking
USE MOD_Particle_Vars                 ,ONLY: Species, tSurfFluxLink, WriteMacroSurfaceValues
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSurfFluxLink),POINTER, INTENT(INOUT) :: currentSurfFluxPart
INTEGER, INTENT(IN)                        :: iSpec, iSF, SideID, iSide, iSample, jSample, ParticleIndexNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: currentBC
!===================================================================================================================================
! check if surfaceflux is used for surface sampling (neccessary for desorption and evaporation)
! create linked list of surfaceflux-particle-info for sampling case
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
 IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)) &
      .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
    IF (PartBound%TargetBoundCond(CurrentBC).EQ.PartBound%ReflectiveBC) THEN
      ! first check if linked list is initialized and initialize if neccessary
      IF (.NOT. ASSOCIATED(currentSurfFluxPart)) THEN
        ALLOCATE(currentSurfFluxPart)
        IF (.NOT. ASSOCIATED(Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart)) THEN
          Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart => currentSurfFluxPart
          Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart  => currentSurfFluxPart
        END IF
      ! check if surfaceflux has already list (happens if second etc. surfaceflux is considered)
      ! create linke to next surfflux-part from current list
      ELSE IF (.NOT. ASSOCIATED(Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart)) THEN
        IF (.NOT. ASSOCIATED(currentSurfFluxPart%next)) THEN
          ALLOCATE(currentSurfFluxPart%next)
        END IF
        currentSurfFluxPart => currentSurfFluxPart%next
        Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart => currentSurfFluxPart
        Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart  => currentSurfFluxPart
      ! surfaceflux has already list but new particle is being inserted
      ! create linke to next surfflux-part from current list
      ELSE
        IF (.NOT. ASSOCIATED(currentSurfFluxPart%next)) THEN
          ALLOCATE(currentSurfFluxPart%next)
        END IF
        currentSurfFluxPart => currentSurfFluxPart%next
        Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart  => currentSurfFluxPart
      END IF
      ! save index and sideinfo for current to be inserted particle
      currentSurfFluxPart%PartIdx = ParticleIndexNbr
      IF (.NOT.TriaTracking .AND. (nSurfSample.GT.1)) THEN
        IF (.NOT. ALLOCATED(currentSurfFluxPart%SideInfo)) ALLOCATE(currentSurfFluxPart%SideInfo(1:3))
        currentSurfFluxPart%SideInfo(1) = iSide
        currentSurfFluxPart%SideInfo(2) = iSample
        currentSurfFluxPart%SideInfo(3) = jSample
      ELSE
        IF (.NOT. ALLOCATED(currentSurfFluxPart%SideInfo)) ALLOCATE(currentSurfFluxPart%SideInfo(1))
        currentSurfFluxPart%SideInfo(1) = SurfMesh%SideIDToSurfID(SideID)
      END IF
    END IF ! sampling is on (CalcSurfaceVal)
  END IF ! wallmodel or liquidsim
END SUBROUTINE CreateLinkedListReactiveBC

SUBROUTINE SamplingForReactiveBC(iSpec, iSF, PartsEmitted, currentSurfFluxPart)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Boundary_Vars        ,ONLY: PartBound, nSurfSample
USE MOD_Particle_Tracking_Vars        ,ONLY: TriaTracking
USE MOD_Particle_Vars                 ,ONLY: Species, tSurfFluxLink, PartState, PEM
USE MOD_Particle_Boundary_Tools       ,ONLY: CalcWallSample
USE MOD_DSMC_Vars                     ,ONLY: CollisMode, PartStateIntEn
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers            ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSurfFluxLink),POINTER, INTENT(INOUT) :: currentSurfFluxPart
INTEGER, INTENT(IN)                        :: iSpec, iSF, PartsEmitted
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: currentBC, PartID, SurfSideID, p, q
REAL                        :: VelXold, VelYold, VelZold, VeloReal, TransArray(1:6),IntArray(1:6)
REAL                        :: EtraOld, EtraWall, EtraNew, ErotOld, ErotWall, ErotNew, EvibOld, EvibWall, EVibNew
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
IF ((PartBound%TargetBoundCond(CurrentBC).EQ.PartBound%ReflectiveBC) .AND. (PartsEmitted.GT.0)) THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  ! check if surfaceflux is used for surface sampling (neccessary for desorption and evaporation)
  ! only allocated if sampling and surface model enabled
  currentSurfFluxPart => Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart
  DO WHILE(ASSOCIATED(currentSurfFluxPart))
    PartID     = currentSurfFluxPart%PartIdx
    SurfSideID = currentSurfFluxPart%SideInfo(1)
    IF (TriaTracking.OR.(nSurfSample.EQ.1)) THEN
      p = 1
      q = 1
    ELSE
      p = currentSurfFluxPart%SideInfo(2)
      q = currentSurfFluxPart%SideInfo(3)
    END IF
    ! set velocities and translational energies
    VelXold  = PartBound%WallVelo(1,CurrentBC)
    VelYold  = PartBound%WallVelo(2,CurrentBC)
    VelZold  = PartBound%WallVelo(3,CurrentBC)
    EtraOld = 0.0
    EtraWall = EtraOld
    VeloReal = VECNORM(PartState(4:6,PartID))
    EtraNew = 0.5 * Species(iSpec)%MassIC * VeloReal**2
    ! fill Transarray
    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    ! must be old_velocity-new_velocity
    TransArray(4) = VelXold-PartState(4,PartID)
    TransArray(5) = VelYold-PartState(5,PartID)
    TransArray(6) = VelZold-PartState(6,PartID)
    IF (CollisMode.GT.1) THEN
      ! set rotational energies
      ErotWall = 0
      ErotOld  = ErotWall
      ErotNew  = PartStateIntEn(2,PartID)
      ! fill rotational internal array
      IntArray(1) = ErotOld
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      ! set vibrational energies
      EvibWall = 0 ! calculated and added in particle desorption calculation
      EvibOld  = EvibWall ! calculated and added in particle desorption calculation
      EvibNew  = PartStateIntEn(1,PartID)
      ! fill vibrational internal array
      IntArray(4) = EvibOld
      IntArray(5) = EvibWall
      IntArray(6) = EvibNew
    ELSE
      IntArray(:) = 0.
    END IF
    ! sample values
    CALL CalcWallSample(PartID,SurfSideID,p,q,TransArray,IntArray,.False.,emission_opt=.TRUE.)
    currentSurfFluxPart => currentSurfFluxPart%next
#if USE_LOADBALANCE
    CALL LBElemSplitTime(PEM%LocalElemID(PartID),tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (ASSOCIATED(currentSurfFluxPart,Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart%next)) THEN
      currentSurfFluxPart => Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart
      EXIT
    END IF
  END DO
END IF ! reflective bc
IF (ASSOCIATED(Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart)) THEN
  Species(iSpec)%Surfaceflux(iSF)%firstSurfFluxPart  => NULL()
  Species(iSpec)%Surfaceflux(iSF)%lastSurfFluxPart  => NULL()
END IF

END SUBROUTINE SamplingForReactiveBC


#ifdef CODE_ANALYZE
SUBROUTINE AnalyzePartPos(ParticleIndexNbr)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars                 ,ONLY: LastPartPos, PDM, PartDtFrac,PartIsImplicit
USE MOD_Particle_Mesh_Vars            ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                        :: ParticleIndexNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(   (LastPartPos(1,ParticleIndexNbr).GT.GEO%xmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,ParticleIndexNbr),GEO%xmaxglob) &
  .OR.(LastPartPos(1,ParticleIndexNbr).LT.GEO%xminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(1,ParticleIndexNbr),GEO%xminglob) &
  .OR.(LastPartPos(2,ParticleIndexNbr).GT.GEO%ymaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,ParticleIndexNbr),GEO%ymaxglob) &
  .OR.(LastPartPos(2,ParticleIndexNbr).LT.GEO%yminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(2,ParticleIndexNbr),GEO%yminglob) &
  .OR.(LastPartPos(3,ParticleIndexNbr).GT.GEO%zmaxglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,ParticleIndexNbr),GEO%zmaxglob) &
  .OR.(LastPartPos(3,ParticleIndexNbr).LT.GEO%zminglob).AND. .NOT.ALMOSTEQUAL(LastPartPos(3,ParticleIndexNbr),GEO%zminglob) ) THEN
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' ParticleInside ',PDM%ParticleInside(ParticleIndexNbr)
#ifdef IMPA
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PartIsImplicit ', PartIsImplicit(ParticleIndexNbr)
  IPWRITE(UNIt_stdOut,'(I0,A18,ES25.14)')                       ' PartDtFrac ', PartDtFrac(ParticleIndexNbr)
#endif /*IMPA*/
  IPWRITE(UNIt_stdOut,'(I0,A18,L)')                            ' PDM%IsNewPart ', PDM%IsNewPart(ParticleIndexNbr)
  IPWRITE(UNIt_stdOut,'(I0,A18,x,A18,x,A18)')                  '    min ', ' value ', ' max '
  IPWRITE(UNIt_stdOut,'(I0,A2,x,ES25.14,x,ES25.14,x,ES25.14)') ' x', GEO%xminglob, LastPartPos(1,ParticleIndexNbr) &
                                                                , GEO%xmaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,ES25.14,x,ES25.14,x,ES25.14)') ' y', GEO%yminglob, LastPartPos(2,ParticleIndexNbr) &
                                                                , GEO%ymaxglob
  IPWRITE(UNIt_stdOut,'(I0,A2,x,ES25.14,x,ES25.14,x,ES25.14)') ' z', GEO%zminglob, LastPartPos(3,ParticleIndexNbr) &
                                                                , GEO%zmaxglob
  CALL abort(&
     __STAMP__ &
#if  defined(IMPA) || defined(ROS)
     ,' LastPartPos outside of mesh. iPart=, iStage',ParticleIndexNbr,REAL(iStage))
#else
     ,' LastPartPos outside of mesh. iPart=',ParticleIndexNbr)
#endif
END IF
END SUBROUTINE AnalyzePartPos
#endif /*CODE_ANALYZE*/

SUBROUTINE SetInnerEnergies(iSpec, iSF, NbrOfParticle)
!===================================================================================================================================
! SideList for SurfaceFlux in BCdata_auxSF is created. Furthermore, the side areas are corrected for Symmetry2D case and finally
! communicated.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
USE MOD_Particle_Vars           ,ONLY: PDM
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_DSMC_Init               ,ONLY: DSMC_SetInternalEnr_LauxVFD
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                        :: iSpec, iSF, NbrOfParticle
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
    IF (SpecDSMC(iSpec)%PolyatomicMol) THEN
      CALL DSMC_SetInternalEnr_Poly(iSpec,iSF,PositionNbr,2)
    ELSE
      CALL DSMC_SetInternalEnr_LauxVFD(iSpec, iSF, PositionNbr,2)
    END IF
  END IF
  iPart = iPart + 1
END DO
END SUBROUTINE SetInnerEnergies

SUBROUTINE ParticleSurfaceflux()
!===================================================================================================================================
! Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
!===================================================================================================================================
! Modules
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_DSMC_Symmetry2D         ,ONLY: CalcRadWeightMPF
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, RadialWeighting
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
USE MOD_MacroBody_Tools         ,ONLY: INSIDEMACROBODY
USE MOD_MacroBody_Vars          ,ONLY: UseMacroBody
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Part_Tools              ,ONLY: GetParticleWeight
USE MOD_Part_Emission_Tools     ,ONLY: SetParticleChargeAndMass, SetParticleMPF
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance, CalcMassflowRate, nPartIn, PartEkinIn
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfFluxSideSize, TriaSurfaceFlux, BCdata_auxSF
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Timedisc_Vars           ,ONLY: RKdtFrac, dt
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping
#endif /*IMPA*/
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
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, currentBC, SideID, NbrOfParticle, ExtraParts, ParticleIndexNbr
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSubSide, iPart, iPartTotal
INTEGER                     :: nReject, allowedRejections, PartsEmitted, Node1, Node2, globElemId
INTEGER                     :: PartInsSideSubSub(1:RadialWeighting%nSubSides)
REAL                        :: Particle_pos(3), RandVal1,  xyzNod(3), RVec(2), minPos(2), xi(2), Vector1(3), Vector2(3)
REAL                        :: ndist(3), midpoint(3)
LOGICAL                     :: AcceptPos
REAL,ALLOCATABLE            :: particle_positions(:), particle_xis(:)
INTEGER,ALLOCATABLE         :: PartInsSubSides(:,:,:)
TYPE(tSurfFluxLink),POINTER :: currentSurfFluxPart => NULL()
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    PartsEmitted = 0
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    NbrOfParticle = 0 ! calculated within (sub)side-Loops!
    iPartTotal=0
    ! Reset the mass flow rate counter for the next time step
    IF(CalcMassflowRate) Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = 0.
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) CALL AdaptiveBoundary_ConstMassflow_Weight(iSpec,iSF)
    END IF
    !Calc Particles for insertion in standard case
    IF ((.NOT.DoPoissonRounding).AND.(.NOT. DoTimeDepInflow).AND.(.NOT.RadialWeighting%DoRadialWeighting) &
        .AND.(.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive)) CALL CalcPartInsSubSidesStandardCase(iSpec,iSF, PartInsSubSides)

!----- 0.: go through (sub)sides if present in proc
    IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) THEN
      CYCLE
    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
#if USE_LOADBALANCE
    CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
    DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
      IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
        IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
      END IF
      BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
      ElemID = SideToElem(1,BCSideID)
      IF (ElemID.LT.1) THEN !not sure if necessary
        ElemID = SideToElem(2,BCSideID)
        iLocSide = SideToElem(4,BCSideID)
      ELSE
        iLocSide = SideToElem(3,BCSideID)
      END IF
      globElemId = ElemID + offSetElem
      SideID=GetGlobalNonUniqueSideID(globElemId,iLocSide)
      IF (TriaSurfaceFlux) xyzNod(1:3) = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1:3)
      DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
        ExtraParts = 0 !set here number of additional to-be-inserted particles in current BCSideID/subsides (e.g. desorption)
        IF(Symmetry2DAxisymmetric.AND.(jSample.EQ.2)) CYCLE
        IF (TriaSurfaceFlux) THEN
          !-- compute parallelogram of triangle
          Node1 = jSample+1     ! normal = cross product of 1-2 and 1-3 for first triangle
          Node2 = jSample+2     !          and 1-3 and 1-4 for second triangle
          Vector1 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
          Vector2 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
          midpoint(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
          ndist(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)
        END IF

        IF (PartBound%Reactive(currentBC)) CALL CalcExtraPartsReactiveBC(iSpec,iSF)
        ! REQUIRED LATER FOR THE POSITION START
        IF(Symmetry2DAxisymmetric) CALL DefineSideDirectVec2D(SideID, xyzNod, minPos, RVec)

        !-- compute number of to be inserted particles
        IF (.NOT.RadialWeighting%DoRadialWeighting) THEN
          IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
            IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
              PartInsSubSide=PartInsSubSides(iSample,jSample,iSide)
            ELSE IF(DoPoissonRounding .AND. .NOT.DoTimeDepInflow)THEN
              CALL CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
            ELSE !DoTimeDepInflow
              CALL RANDOM_NUMBER(RandVal1)
              PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
                             * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR+RandVal1)
            END IF !DoPoissonRounding
          ELSE !Species(iSpec)%Surfaceflux(iSF)%Adaptive
            CALL CalcPartInsAdaptive(iSpec, iSF, BCSideID, iSide, iSample, jSample, PartInsSubSide)
          END IF ! Adaptive SurfaceFlux
        ELSE
          CALL CalcPartInsRadWeight(iSpec, iSF, iSample, jSample, iSide, minPos, RVec, PartInsSubSide, PartInsSideSubSub)
        END IF ! noAdaptive.AND.(.NOT.Symmetry2DAxisymmetric)
        !-- proceed with calculated to be inserted particles
        IF (PartInsSubSide.LT.0) THEN
          CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: PartInsSubSide.LT.0!')
        ELSE IF (PartInsSubSide + ExtraParts.LE.0) THEN
          CYCLE
        END IF
        PartInsSubSide = PartInsSubSide + ExtraParts
        NbrOfParticle = NbrOfParticle + PartInsSubSide
        ALLOCATE(particle_positions(1:PartInsSubSide*3))
        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
          ALLOCATE( particle_xis(1:PartInsSubSide*2))
        END IF !VeloIsNormal
        !-- put particles in subside (rejections are used if contraint reduces actual inserted number)
        iPart=1
        nReject=0
        allowedRejections=0

        !-- Set Positions
        IF(Symmetry2DAxisymmetric) THEN
          CALL CalcPartPosRadWeight(minPos, RVec, PartInsSubSide, PartInsSideSubSub, particle_positions)
        ELSE
          DO WHILE (iPart+allowedRejections .LE. PartInsSubSide)
            IF (TriaSurfaceFlux) THEN
              Particle_pos(1:3) = CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
            ELSE !.NOT.TriaSurfaceFlux
              Particle_pos(1:3) = CalcPartPosBezier(iSpec, iSF, iSample, jSample, iSide, SideID)
            END IF !TriaSurfaceFlux

            AcceptPos=.TRUE.
            IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN !check rmax-rejection
              IF (.NOT.InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)) AcceptPos=.FALSE.
            END IF ! CircularInflow
            IF (UseMacroBody) THEN
              IF (INSIDEMACROBODY(Particle_pos)) AcceptPos=.FALSE.
            END IF
            !-- save position if accepted:
            IF (AcceptPos) THEN
              particle_positions(iPart*3-2) = Particle_pos(1)
              particle_positions(iPart*3-1) = Particle_pos(2)
              particle_positions(iPart*3  ) = Particle_pos(3)
              IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
                particle_xis(iPart*2-1) = xi(1)
                particle_xis(iPart*2  ) = xi(2)
              END IF !VeloIsNormal
              iPart=iPart+1
            ELSE
              nReject=nReject+1
              IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow .OR. UseMacroBody) THEN !check rmax-rejection
                allowedRejections=allowedRejections+1
              END IF
            END IF
          END DO !put particles in subside: WHILE(iPart+allowedRejections .LE. PartInsSubSide)
        END IF
        PartInsSubSide = PartInsSubSide - allowedRejections
        NbrOfParticle = NbrOfParticle - allowedRejections

        !-- Fill Particle Informations (PartState, Partelem, etc.)
        ParticleIndexNbr = 1
        DO iPart=1,PartInsSubSide
          IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) &
              ParticleIndexNbr = PDM%nextFreePosition(iPartTotal + 1 + PDM%CurrentNextFreePosition)
          IF (ParticleIndexNbr .ne. 0) THEN
            PartState(1:3,ParticleIndexNbr) = particle_positions(3*(iPart-1)+1:3*(iPart-1)+3)
            IF (PartBound%Reactive(CurrentBC)) &
               CALL CreateLinkedListReactiveBC(iSpec, iSF, SideID, iSide, iSample, jSample, ParticleIndexNbr, currentSurfFluxPart)
            IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal.AND.(.NOT.TriaSurfaceFlux)) THEN
              PartState(4:5,ParticleIndexNbr) = particle_xis(2*(iPart-1)+1:2*(iPart-1)+2) !use velo as dummy-storage for xi!
            END IF
            LastPartPos(1:3,ParticleIndexNbr)=PartState(1:3,ParticleIndexNbr)
#if defined(IMPA) || defined(ROS)
            IF(DoRefMapping) CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),PartPosRef(1:3,ParticleIndexNbr),globElemId)
#endif /*IMPA*/
            PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
            PDM%dtFracPush(ParticleIndexNbr) = .TRUE.
            PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
            PEM%GlobalElemID(ParticleIndexNbr) = globElemId
            PEM%LastGlobalElemID(ParticleIndexNbr) = globElemId !needed when ParticlePush is not executed, e.g. "delay"
            iPartTotal = iPartTotal + 1
            IF (VarTimeStep%UseVariableTimeStep) THEN
              VarTimeStep%ParticleTimeStep(ParticleIndexNbr) &
                = CalcVarTimeStep(PartState(1,ParticleIndexNbr),PartState(2,ParticleIndexNbr),PEM%LocalElemID(ParticleIndexNbr))
            END IF
            IF (RadialWeighting%DoRadialWeighting) THEN
              PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(2,ParticleIndexNbr), 1,ParticleIndexNbr)
            END IF
            IF(CalcMassflowRate) THEN
              Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow &
                                                                + GetParticleWeight(ParticleIndexNbr)
            END IF
#ifdef CODE_ANALYZE
            CALL AnalyzePartPos(ParticleIndexNbr)
#endif /*CODE_ANALYZE*/
          ELSE
            CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
          END IF
        END DO
        DEALLOCATE(particle_positions)
        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) DEALLOCATE(particle_xis)
!----- 2a.: set velocities if special for each subside
        CALL SetSurfacefluxVelocities(iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartInsSubSide)

        PartsEmitted = PartsEmitted + PartInsSubSide
#if USE_LOADBALANCE
        !used for calculating LoadBalance of tCurrent(LB_SURFFLUX) ==> "2b.: set remaining properties"
        nSurfacefluxPerElem(ElemID)=nSurfacefluxPerElem(ElemID)+PartInsSubSide
#endif /*USE_LOADBALANCE*/

      END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
#if USE_LOADBALANCE
      CALL LBElemSplitTime(ElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
    END DO ! iSide

    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut = 0
    END IF
    IF (NbrOfParticle.NE.iPartTotal) CALL abort(&
__STAMP__&
, 'Error 2 in ParticleSurfaceflux!')
!----- 2b.: set remaining properties
    CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)
    IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(iSpec,NbrOfParticle)
    ! define molecule stuff
    IF (useDSMC.AND.(CollisMode.GT.1)) CALL SetInnerEnergies(iSpec, iSF, NbrOfParticle)
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
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/
    ! Sample Energies on Surfaces when particles are emitted from them
    IF (NbrOfParticle.NE.PartsEmitted) THEN
      ! should be equal for including the following lines in tSurfaceFlux
      CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: NbrOfParticle.NE.PartsEmitted')
    END IF
    IF (PartBound%Reactive(CurrentBC)) CALL SamplingForReactiveBC(iSpec, iSF, PartsEmitted, currentSurfFluxPart)
  END DO !iSF
END DO !iSpec

END SUBROUTINE ParticleSurfaceflux

FUNCTION InSideCircularInflow(iSpec, iSF, iSide, Particle_pos)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER, INTENT(IN)             :: iSpec, iSF, iSide
REAL, INTENT(IN)                :: Particle_pos(3)
LOGICAL                         :: InSideCircularInflow
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: point(2), radius, origin(2)
!===================================================================================================================================
  origin=Species(iSpec)%Surfaceflux(iSF)%origin
  SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide))
  CASE(0) !- RejectType=0 : complete side is inside valid bounds
    InSideCircularInflow=.TRUE.
  CASE(1) !- RejectType=1 : complete side is outside of valid bounds
  CALL abort(&
  __STAMP__&
  ,'side outside of valid bounds was considered although nVFR=0...?!')
                !AcceptPos=.FALSE.
  CASE(2) !- RejectType=2 : side is partly inside valid bounds
    point(1)=Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(2))-origin(1)
    point(2)=Particle_pos(Species(iSpec)%Surfaceflux(iSF)%dir(3))-origin(2)
    radius=SQRT( (point(1))**2+(point(2))**2 )
    IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
      InSideCircularInflow=.TRUE.
    ELSE
      InSideCircularInflow=.FALSE.
    END IF
  CASE DEFAULT
    CALL abort(&
  __STAMP__&
  ,'wrong SurfFluxSideRejectType!')
  END SELECT !SurfFluxSideRejectType

END FUNCTION InSideCircularInflow


FUNCTION CalcPartPosBezier(iSpec, iSF, iSample, jSample, iSide, SideID)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D,BezierSampleXi
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_Mesh_Vars               ,ONLY: NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES)
INTEGER, INTENT(IN)         :: iSpec, iSF, iSide, SideID, iSample, jSample
REAL                        :: CalcPartPosBezier(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal2(2), xiab(1:2,1:2), xi(2), E, F, G, D, gradXiEta2D(1:2,1:2),gradXiEta3D(1:2,1:3), RandVal1
INTEGER                     :: iLoop
!===================================================================================================================================
  iLoop=0
  DO !ARM for xi considering the dA of the Subside in RefSpace
    iLoop = iLoop+1
    CALL RANDOM_NUMBER(RandVal2)
    xiab(1,1:2)=(/BezierSampleXi(iSample-1),BezierSampleXi(iSample)/) !correct order?!?
    xiab(2,1:2)=(/BezierSampleXi(JSample-1),BezierSampleXi(JSample)/) !correct order?!?
    xi=(xiab(:,2)-xiab(:,1))*RandVal2+xiab(:,1)
    IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
      IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
        CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,2 &
          ,Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample &
          ,iSide)%BezierControlPoints2D(1:2,0:NGeo,0:NGeo) &
          ,Gradient=gradXiEta2D)
        E=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(1,1:2))
        F=DOT_PRODUCT(gradXiEta2D(1,1:2),gradXiEta2D(2,1:2))
        G=DOT_PRODUCT(gradXiEta2D(2,1:2),gradXiEta2D(2,1:2))
      ELSE
        CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID) &
          ,Gradient=gradXiEta3D)
        E=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(1,1:3))
        F=DOT_PRODUCT(gradXiEta3D(1,1:3),gradXiEta3D(2,1:3))
        G=DOT_PRODUCT(gradXiEta3D(2,1:3),gradXiEta3D(2,1:3))
      END IF !.NOT.VeloIsNormal
      D=SQRT(E*G-F*F)
      D=D/Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Dmax !scaled Jacobian of xi
      IF (D .GT. 1.01) THEN !arbitrary warning threshold
        IPWRITE(*,'(I4,x,A28,I0,A9,I0,A22,I0)') &
          'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has inaccurate Dmax! ',D
      END IF
      CALL RANDOM_NUMBER(RandVal1)
      IF (RandVal1.LE.D) THEN
        EXIT !accept xi
      ELSE
        IF (MOD(iLoop,100).EQ.0) THEN !arbitrary warning threshold
          IPWRITE(*,'(I4,x,A28,I0,A9,I0,A18,I0)') &
            'WARNING: ARM of SurfaceFlux ',iSF,' of Spec ',iSpec,' has reached loop ',iLoop
          IPWRITE(*,'(I4,x,A19,2(x,E16.8))') &
            '         R, D/Dmax:',RandVal1,D
        END IF
      END IF
    ELSE !no ARM -> accept xi
      EXIT
    END IF
  END DO !Jacobian-based ARM-loop
  IF(MINVAL(XI).LT.-1.)THEN
    IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi<-1',XI
  END IF
  IF(MAXVAL(XI).GT.1.)THEN
    IPWRITE(UNIT_StdOut,'(I0,A,E16.8)') ' Xi>1',XI
  END IF
  CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=CalcPartPosBezier)

END FUNCTION CalcPartPosBezier

FUNCTION CalcPartPosTriaSurface(xyzNod, Vector1, Vector2, ndist, midpoint)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
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
  IF (.NOT.TriaTracking) THEN !prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
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

SUBROUTINE CalcPartPosRadWeight(minPos, RVec, PartInsSubSide, PartInsSideSubSub, particle_positions)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: PartInsSubSide, PartInsSideSubSub(:)
REAL, INTENT(IN)            :: minPos(2), RVec(2)
REAL, INTENT(OUT)           :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal1, PminTemp, PmaxTemp, Particle_pos(2)
INTEGER                     :: iSub, iPart, iPartSub
!===================================================================================================================================
iPart=1
IF (RadialWeighting%DoRadialWeighting.AND.(.NOT.(ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))))) THEN
  IF(RadialWeighting%CellLocalWeighting) THEN
    DO WHILE (iPart .LE. PartInsSubSide)
      CALL RANDOM_NUMBER(RandVal1)
      Particle_pos(2) = minPos(2) + RandVal1 * RVec(2)
      ! x-position depending on the y-location
      Particle_pos(1) = minPos(1) + (Particle_pos(2)-minPos(2)) * RVec(1) / RVec(2)
      particle_positions(iPart*3-2) = Particle_pos(1)
      particle_positions(iPart*3-1) = Particle_pos(2)
      particle_positions(iPart*3  ) = 0.
      iPart = iPart + 1
    END DO
  ELSE
    DO iSub = 1, RadialWeighting%nSubSides
      iPartSub = 1
      DO WHILE (iPartSub.LE.PartInsSideSubSub(iSub))
        CALL RANDOM_NUMBER(RandVal1)
        PminTemp = minPos(2) + RVec(2)/RadialWeighting%nSubSides*(iSub-1.)
        PmaxTemp = minPos(2) + RVec(2)/RadialWeighting%nSubSides*iSub
        Particle_pos(2) = PminTemp + RandVal1 * (PmaxTemp - PminTemp)
        ! x-position depending on the y-location
        Particle_pos(1) = minPos(1) + (Particle_pos(2)-minPos(2)) * RVec(1) / RVec(2)
        particle_positions(iPart*3-2) = Particle_pos(1)
        particle_positions(iPart*3-1) = Particle_pos(2)
        particle_positions(iPart*3  ) = 0.
        iPart = iPart + 1
        iPartSub = iPartSub + 1
      END DO
    END DO
  END IF
ELSE
  DO WHILE (iPart .LE. PartInsSubSide)
    CALL RANDOM_NUMBER(RandVal1)
    IF (ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
      ! y_min = y_max, faces parallel to x-direction, constant distribution
      Particle_pos(1:2) = minPos(1:2) + RVec(1:2) * RandVal1
    ELSE
    ! No RadialWeighting, regular linear distribution of particle positions
      Particle_pos(1:2) = minPos(1:2) + RVec(1:2) &
          * ( SQRT(RandVal1*((minPos(2) + RVec(2))**2-minPos(2)**2)+minPos(2)**2) - minPos(2) ) / (RVec(2))
    END IF
    particle_positions(iPart*3-2) = Particle_pos(1)
    particle_positions(iPart*3-1) = Particle_pos(2)
    particle_positions(iPart*3  ) = 0.
    iPart = iPart + 1
  END DO
END IF

END SUBROUTINE CalcPartPosRadWeight

SUBROUTINE CalcPartInsRadWeight(iSpec, iSF, iSample, jSample, iSide, minPos, RVec, PartInsSubSide, PartInsSideSubSub)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_DSMC_Vars               ,ONLY: RadialWeighting
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSample, jSample, iSide
REAL, INTENT(IN)            :: minPos(2), RVec(2)
INTEGER, INTENT(OUT)        :: PartInsSubSide, PartInsSideSubSub(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: RandVal1
INTEGER                     :: iSub
!===================================================================================================================================
CALL RANDOM_NUMBER(RandVal1)
PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
  * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR + RandVal1)
IF(.NOT.RadialWeighting%CellLocalWeighting) THEN
  IF(.NOT.ALMOSTEQUAL(minPos(2),minPos(2)+RVec(2))) THEN
    PartInsSubSide = 0
    DO iSub = 1, RadialWeighting%nSubSides
      CALL RANDOM_NUMBER(RandVal1)
      PartInsSideSubSub(iSub) = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
              * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%nVFRSub(iSide,iSub)+ RandVal1)
      PartInsSubSide = PartInsSubSide + PartInsSideSubSub(iSub)
    END DO
  END IF
END IF

END SUBROUTINE CalcPartInsRadWeight

SUBROUTINE CalcPartInsPoissonDistr(iSpec, iSF, iSample, jSample, iSide, PartInsSubSide)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Part_Emission_Tools     ,ONLY: SamplePoissonDistri
USE MOD_Particle_Vars           ,ONLY: Species
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, iSample, jSample, iSide
INTEGER, INTENT(OUT)        :: PartInsSubSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: PartIns
!===================================================================================================================================
PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
                      * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR
IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in ParticleSurfaceflux: flux is too large for poisson sampling!')
ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
  CALL SamplePoissonDistri( PartIns , PartInsSubSide )
END IF

END SUBROUTINE CalcPartInsPoissonDistr

SUBROUTINE CalcPartInsAdaptive(iSpec, iSF, BCSideID, iSide, iSample, jSample, PartInsSubSide)
!===================================================================================================================================
! Calculate random normalized vector in 3D (unit space)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi
USE MOD_TimeDisc_Vars           ,ONLY: dt,RKdtFrac
USE MOD_Particle_Vars           ,ONLY: Species, Adaptive_MacroVal
USE MOD_Particle_Surfaces_Vars  ,ONLY: SurfMeshSubSideData
USE MOD_Mesh_Vars               ,ONLY: SideToElem
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: iSpec, iSF, BCSideID, iSide, jSample, iSample
INTEGER, INTENT(OUT)        :: PartInsSubSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                        :: ElemPartDensity, T, pressure, VeloVec(3), vec_nIn(3), veloNormal, VeloIC, VeloVecIC(3)
REAL                        :: projFak, a, v_thermal, vSF, nVFR, RandVal1
INTEGER                     :: ElemID
!===================================================================================================================================
  ElemID = SideToElem(1,BCSideID)
  IF (ElemID.LT.1) ElemID = SideToElem(2,BCSideID)
  SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType)
  CASE(1) ! Pressure inlet (pressure, temperature const)
    ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%PartDensity
    T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(2) ! adaptive Outlet/freestream
    ElemPartDensity = Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec)
    pressure = Species(iSpec)%Surfaceflux(iSF)%AdaptivePressure
    T = pressure / (BoltzmannConst * Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec))
    !T = SQRT(Adaptive_MacroVal(4,ElemID,iSpec)**2+Adaptive_MacroVal(5,ElemID,iSpec)**2 &
    !  + Adaptive_MacroVal(6,ElemID,iSpec)**2)
  CASE(3) ! Mass flow, temperature constant
    VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec)
    VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec)
    VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec)
    vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
    veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
    IF(veloNormal.GT.0.0) THEN
      ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
                        / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
      Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1:3,ElemID) = VeloVec(1:3)
    ELSE
      ! Using the old velocity vector, overwriting the sampled value with the old one
      Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1,ElemID)
      Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(2,ElemID)
      Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(3,ElemID)
      VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec)
      VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec)
      VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec)
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
      IF(veloNormal.GT.0.0) THEN
        ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
          / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
      ELSE
        SWRITE(*,*) 'WARNING: No particles inserted!'
        SWRITE(*,*) 'WARNING: Possibly different adaptive BCs of Type3/4 have been defined next to each other.'
        SWRITE(*,*) 'WARNING: Adaptive BCs sharing a mesh element is currently not supported -> wrong velocity vector!'
        ElemPartDensity = 0
      END IF
    END IF
    T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(4) !Const. massflow inlet after Lei 2017
    T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
  CASE DEFAULT
    SWRITE(*,*) 'Selected adaptive boundary condition type: ', Species(iSpec)%Surfaceflux(iSF)%AdaptiveType
    CALL abort(&
  __STAMP__&
  ,'ERROR Adaptive Inlet: Wrong adaptive type for Surfaceflux!')
  END SELECT
  VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec)
  VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec)
  VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec)
  VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
  IF (ABS(VeloIC).GT.0.) THEN
    VeloVecIC = VeloVec / VeloIC
  ELSE
    VeloVecIC = (/1.,0.,0./)
  END IF
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
  v_thermal = SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC) !thermal speed
  a = 0 !dummy for projected speed ratio in constant v-distri
  !-- compute total volume flow rate through surface
  SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
  CASE('constant')
    vSF = VeloIC * projFak !Velo proj. to inwards normal
    nVFR = MAX(SurfMeshSubSideData(iSample,jSample,BCSideID)%area * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
  CASE('liquid')
    vSF = v_thermal / (2.0*SQRT(PI)) !mean flux velocity through normal sub-face
    nVFR = SurfMeshSubSideData(iSample,jSample,BCSideID)%area * vSF !VFR projected to inwards normal of sub-side
  CASE('maxwell','maxwell_lpn')
    IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
      v_thermal = 1.
    END IF
    a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
    vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
    nVFR = SurfMeshSubSideData(iSample,jSample,BCSideID)%area * vSF !VFR projected to inwards normal of sub-side
  CASE DEFAULT
    CALL abort(&
  __STAMP__&
  ,'wrong velo-distri for adaptive Surfaceflux!')
  END SELECT
  IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) THEN
    CALL RANDOM_NUMBER(RandVal1)
    PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide)     &
                            * (Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow * dt*RKdtFrac    &
                                / (Species(iSpec)%MassIC * Species(iSpec)%MacroParticleFactor)  &
                                + REAL(Species(iSpec)%Surfaceflux(iSF)%AdaptivePartNumOut)) +RandVal1)
  ELSE
    CALL RANDOM_NUMBER(RandVal1)
    PartInsSubSide = INT(ElemPartDensity / Species(iSpec)%MacroParticleFactor * dt*RKdtFrac * nVFR+RandVal1)
  END IF

END SUBROUTINE CalcPartInsAdaptive


SUBROUTINE AdaptiveBoundary_ConstMassflow_Weight(iSpec,iSF)
!===================================================================================================================================
!> Routine calculates the weights of the triangles for AdaptiveType=4 to scale up the number of particles to be inserted
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!                                                                                              ! ----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst, Pi
USE MOD_Particle_Vars          ,ONLY: Species, Adaptive_MacroVal
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfMeshSubSideData, BCdata_auxSF, SurfFluxSideSize
USE MOD_TimeDisc_Vars          ,ONLY: dt, RKdtFrac
USE MOD_Mesh_Vars              ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools    ,ONLY: GetGlobalNonUniqueSideID
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iSF
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSide, BCSideID, ElemID, iLocSide, SideID, currentBC, PartInsSubSum, iSample, jSample
INTEGER, ALLOCATABLE            :: PartInsSubSidesAdapt(:,:,:)
REAL                            :: VeloVec(1:3), vec_nIn(1:3), veloNormal, T, ElemPartDensity, VeloIC, VeloVecIC(1:3), projFak
REAL                            :: v_thermal, a, vSF, nVFR, RandVal1, area
!===================================================================================================================================

currentBC = Species(iSpec)%Surfaceflux(iSF)%BC

SDEALLOCATE(PartInsSubSidesAdapt)
ALLOCATE(PartInsSubSidesAdapt(1:SurfFluxSideSize(1),1:SurfFluxSideSize(2),1:BCdata_auxSF(currentBC)%SideNumber))
PartInsSubSidesAdapt=0

PartInsSubSum = 0

DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
  IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
    IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.1) CYCLE
  END IF
  BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
  ElemID = SideToElem(1,BCSideID)
  IF (ElemID.LT.1) THEN !not sure if necessary
    ElemID = SideToElem(2,BCSideID)
    iLocSide = SideToElem(4,BCSideID)
  ELSE
    iLocSide = SideToElem(3,BCSideID)
  END IF
  SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
  DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
    VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec)
    VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec)
    VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec)
    vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
    veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
    IF(veloNormal.GT.0.0) THEN
      ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
                        / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
      Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1:3,ElemID) = VeloVec(1:3)
    ELSE
      ! Using the old velocity vector, overwriting the sampled value with the old one
      Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(1,ElemID)
      Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(2,ElemID)
      Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec) = Species(iSpec)%Surfaceflux(iSF)%AdaptivePreviousVelocity(3,ElemID)
      VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec)
      VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec)
      VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec)
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      veloNormal = VeloVec(1)*vec_nIn(1) + VeloVec(2)*vec_nIn(2) + VeloVec(3)*vec_nIn(3)
      IF(veloNormal.GT.0.0) THEN
        ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%AdaptiveMassflow &
                        / (veloNormal * Species(iSpec)%Surfaceflux(iSF)%totalAreaSF * Species(iSpec)%MassIC)
      ELSE
        SWRITE(*,*) 'WARNING: Negative/zero velocity at the adaptive boundary, Type 4, no particles inserted! iSF: ', iSF
        ElemPartDensity = 0
      END IF
    END IF
    T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
    VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
    IF (ABS(VeloIC).GT.0.) THEN
      VeloVecIC = VeloVec / VeloIC
    ELSE
      VeloVecIC = (/1.,0.,0./)
    END IF
    projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
    v_thermal = SQRT(2.*BoltzmannConst*T/Species(iSpec)%MassIC) !thermal speed
    a = 0 !dummy for projected speed ratio in constant v-distri
    IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      area = Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(iSample,jSample,iSide)
    ELSE
      area = SurfMeshSubSideData(iSample,jSample,BCSideID)%area
    END IF
    !-- compute total volume flow rate through surface
    SELECT CASE(TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution))
    CASE('constant')
      vSF = VeloIC * projFak !Velo proj. to inwards normal
      nVFR = MAX(area * vSF,0.) !VFR proj. to inwards normal (only positive parts!)
    CASE('maxwell','maxwell_lpn')
      IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
        v_thermal = 1.
      END IF
      a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
      vSF = v_thermal / (2.0*SQRT(PI)) * ( EXP(-(a*a)) + a*SQRT(PI)*(1+ERF(a)) ) !mean flux velocity through normal sub-face
      nVFR = area * vSF !VFR projected to inwards normal of sub-side
    CASE DEFAULT
      CALL abort(&
        __STAMP__&
        ,'wrong velo-distri for adaptive Surfaceflux!')
    END SELECT
    CALL RANDOM_NUMBER(RandVal1)
    PartInsSubSidesAdapt(iSample,jSample,iSide) = INT(ElemPartDensity/Species(iSpec)%MacroParticleFactor*dt*RKdtFrac*nVFR+RandVal1)
    PartInsSubSum = PartInsSubSum + PartInsSubSidesAdapt(iSample,jSample,iSide)
  END DO; END DO
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,PartInsSubSum,1,MPI_INTEGER,MPI_SUM,PartMPI%COMM,IERROR)
#endif

IF(PartInsSubSum.GT.0) THEN
  Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(:,:,:) = REAL(PartInsSubSidesAdapt(:,:,:)) / REAL(PartInsSubSum)
ELSE
  Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(:,:,:) = 0.
END IF

IF(Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
  ! Scaling up the number of particles to be inserted on the triaside
  DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
    BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
    DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
      IF(Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(iSample,jSample,iSide).GT.0.0) THEN
        Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide) = &
          Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide) &
            * SurfMeshSubSideData(iSample,jSample,BCSideID)%area &
            / Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(iSample,jSample,iSide)
      ELSE
        Species(iSpec)%Surfaceflux(iSF)%ConstMassflowWeight(iSample,jSample,iSide) = 0.0
      END IF
    END DO; END DO
  END DO
END IF

SDEALLOCATE(PartInsSubSidesAdapt)

END SUBROUTINE AdaptiveBoundary_ConstMassflow_Weight


SUBROUTINE SetSurfacefluxVelocities(FractNbr,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_Part_Tools,             ONLY : VeloFromDistribution
USE MOD_Particle_Surfaces_Vars, ONLY : SurfMeshSubSideData, TriaSurfaceFlux
USE MOD_Particle_Surfaces,      ONLY : CalcNormAndTangBezier
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,PositionNbr,envelope,currentBC
REAL                             :: Vec3D(3), vec_nIn(1:3), vec_t1(1:3), vec_t2(1:3)
REAL                             :: a,zstar,RandVal1,RandVal2(2),RandVal3(3),u,RandN,RandN_save,Velo1,Velo2,Velosq,T,beta,z
LOGICAL                          :: RandN_in_Mem
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
REAL                             :: projFak                          ! VeloVecIC projected to inwards normal of tria
REAL                             :: Velo_t1                          ! Velo comp. of first orth. vector in tria
REAL                             :: Velo_t2                          ! Velo comp. of second orth. vector in tria
REAL                             :: VeloIC
REAL                             :: VeloVec(1:3)
REAL                             :: VeloVecIC(1:3),v_thermal, pressure
!===================================================================================================================================

IF(PartIns.lt.1) RETURN

IF (TRIM(Species(FractNbr)%Surfaceflux(iSF)%velocityDistribution).EQ.'maxwell' .OR. &
  TRIM(Species(FractNbr)%Surfaceflux(iSF)%velocityDistribution).EQ.'maxwell_lpn') THEN
  velocityDistribution='maxwell_surfaceflux'
ELSE IF (TRIM(Species(FractNbr)%Surfaceflux(iSF)%velocityDistribution).EQ.'liquid' ) THEN
  velocityDistribution='liquid'
ELSE IF (TRIM(Species(FractNbr)%Surfaceflux(iSF)%velocityDistribution).EQ.'constant' ) THEN
  velocityDistribution='constant'
  VeloVecIC(1:3) = Species(FractNbr)%Surfaceflux(iSF)%VeloVecIC(1:3)
  VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
ELSE
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')
END IF
RandN_in_Mem=.FALSE.
envelope=-1
currentBC = Species(FractNbr)%Surfaceflux(iSF)%BC

IF (.NOT.Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
  vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
END IF !.NOT.VeloIsNormal

IF(.NOT.Species(FractNbr)%Surfaceflux(iSF)%Adaptive) THEN
  VeloIC = Species(FractNbr)%Surfaceflux(iSF)%VeloIC
  T = Species(FractNbr)%Surfaceflux(iSF)%MWTemperatureIC
  a = Species(FractNbr)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn
  projFak = Species(FractNbr)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak
  Velo_t1 = Species(FractNbr)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1
  Velo_t2 = Species(FractNbr)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2
ELSE !Species(iSpec)%Surfaceflux(iSF)%Adaptive
  SELECT CASE(Species(FractNbr)%Surfaceflux(iSF)%AdaptiveType)
  CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
    T =  Species(FractNbr)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(2) ! adaptive Outlet/freestream
    pressure = Species(FractNbr)%Surfaceflux(iSF)%AdaptivePressure
    T = pressure / (BoltzmannConst * Adaptive_MacroVal(DSMC_NUMDENS,ElemID,FractNbr))
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'wrong adaptive type for Surfaceflux velocities!')
  END SELECT
  VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,FractNbr)
  VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,FractNbr)
  VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,FractNbr)
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  VeloVec(1:3) = DOT_PRODUCT(VeloVec,vec_nIn)*vec_nIn(1:3)
  VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
  IF (ABS(VeloIC).GT.0.) THEN
    VeloVecIC = VeloVec / VeloIC
  ELSE
    VeloVecIC = (/1.,0.,0./)
  END IF
  projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
  v_thermal = SQRT(2.*BoltzmannConst*T/Species(FractNbr)%MassIC) !thermal speed
  IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
    v_thermal = 1.
  END IF
  a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
  Velo_t1 = VeloIC * DOT_PRODUCT(vec_t1,VeloVecIC) !v in t1-dir
  Velo_t2 = VeloIC * DOT_PRODUCT(vec_t2,VeloVecIC) !v in t2-dir
END IF !Adaptive SurfaceFlux

!-- set velocities
SELECT CASE(TRIM(velocityDistribution))
CASE('constant')

  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
!-- In case of side-normal velocities: calc n-vector at particle position, xi was saved in PartState(4:5)
      IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
        vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
        vec_t1(1:3) = 0. !dummy
        vec_t2(1:3) = 0. !dummy
      ELSE IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
        CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
        vec_nIn(1:3) = -vec_nIn(1:3)
        vec_t1(1:3) = 0. !dummy
        vec_t2(1:3) = 0. !dummy
      ELSE
        vec_nIn(1:3) = VeloVecIC(1:3)
      END IF !VeloIsNormal

!-- build complete velo-vector
      Vec3D(1:3) = vec_nIn(1:3) * Species(FractNbr)%Surfaceflux(iSF)%VeloIC
      PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF !PositionNbr .NE. 0
  END DO !i = ...NbrOfParticle
CASE('liquid')
  !-- 0a.: In case of side-normal velocities: calc n-/t-vectors at particle position, xi was saved in PartState(4:5)
  IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
    vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
    vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
    vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
  ELSE IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
    CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),tang1=vec_t1(1:3),tang2=vec_t2(1:3) &
      ,xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
    vec_nIn(1:3) = -vec_nIn(1:3)
  END IF
  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      Vec3D(1:3) = VeloFromDistribution('liquid_evap',FractNbr,T)
      PartState(4:6,PositionNbr) = vec_nIn(1:3)*(a*SQRT(2*BoltzmannConst*T/Species(FractNbr)%MassIC)+Vec3D(3)) &
                                 + vec_t1(1:3)*(Velo_t1+Vec3D(1)) &
                                 + vec_t2(1:3)*(Velo_t2+Vec3D(2))
    ELSE !PositionNbr .EQ. 0
      CALL abort(&
__STAMP__&
,'!PositionNbr .EQ. 0!')
    END IF !PositionNbr .NE. 0
  END DO
CASE('maxwell_surfaceflux')
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
      IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
        vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
        vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
        vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
      ELSE IF (Species(FractNbr)%Surfaceflux(iSF)%VeloIsNormal) THEN
        CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),tang1=vec_t1(1:3),tang2=vec_t2(1:3) &
          ,xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
        vec_nIn(1:3) = -vec_nIn(1:3)
      END IF !VeloIsNormal
!-- 1.: determine zstar (initial generation of potentially too many RVu is for needed indentities of RVu used multiple times!
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
!            IF (.NOT.DoZigguratSampling) THEN !polar method
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
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
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
!            IF (.NOT.DoZigguratSampling) THEN !polar method
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
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
            zstar = 1./SQRT(2.)*RandN
          END IF
          IF ( (a-zstar)/a .GT. RandVal3(3)) THEN
            EXIT
          END IF
        END DO
      CASE DEFAULT
        CALL abort(&
__STAMP__&
,'wrong enevelope in SetSurfacefluxVelocities!')
      END SELECT

!-- 2.: sample normal directions and build complete velo-vector
      Vec3D(1:3) = vec_nIn(1:3) * SQRT(2.*BoltzmannConst*T/Species(FractNbr)%MassIC)*(a-zstar)
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
!      ELSE !ziggurat method
!        Velo1=rnor()
!        Velo2=rnor()
!      END IF
      Vec3D(1:3) = Vec3D(1:3) + vec_t1(1:3) &
        * ( Velo_t1+Velo1*SQRT(BoltzmannConst*T/Species(FractNbr)%MassIC) )     !t1-Komponente (Gauss)
      Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3) &
        * ( Velo_t2+Velo2*SQRT(BoltzmannConst*T/Species(FractNbr)%MassIC) )     !t2-Komponente (Gauss)

      PartState(4:6,PositionNbr) = Vec3D(1:3)
    ELSE !PositionNbr .EQ. 0
      CALL abort(&
__STAMP__&
,'PositionNbr .EQ. 0!')
    END IF !PositionNbr .NE. 0
  END DO !i = ...NbrOfParticle
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')
END SELECT

END SUBROUTINE SetSurfacefluxVelocities


SUBROUTINE CircularInflow_Area(iSpec,iSF,iSide,BCSideID)
!===================================================================================================================================
!> Routine calculates the partial area of the circular infow per triangle for AdaptiveType=4 (Monte Carlo integration)
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY:Species
USE MOD_Particle_Surfaces_Vars      ,ONLY:SurfMeshSubSideData, BCdata_auxSF, SurfFluxSideSize, TriaSurfaceFlux
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpec, iSF, iSide, BCSideID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iMC, iSample, jSample, dir(3), currentBC, MCVar, Node1, Node2, counter
REAL                            :: CircleArea, point(2), origin(2), radius, Vector1(3), Vector2(3), RandVal2(2), Particle_pos(3)
REAL                            :: PartDistance, TriaNode(1:3), midpoint(1:3), ndist(1:3)
!===================================================================================================================================

MCVar = 1000000

IF (.NOT.TriaSurfaceFlux) THEN
  CALL abort(&
__STAMP__&
,'ERROR: CircularInflow only with TriaSurfaceFlux!')
END IF

currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
dir       = Species(iSpec)%Surfaceflux(iSF)%dir
origin    = Species(iSpec)%Surfaceflux(iSF)%origin

TriaNode(1:3) = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1:3)

DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
  !-- compute parallelogram of triangle
  Node1 = jSample+1     ! normal = cross product of 1-2 and 1-3 for first triangle
  Node2 = jSample+2     !          and 1-3 and 1-4 for second triangle
  Vector1 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node1-1)
  Vector2 = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,Node2-1)
  midpoint(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%midpoint(1:3)
  ndist(1:3) = BCdata_auxSF(currentBC)%TriaSwapGeo(iSample,jSample,iSide)%ndist(1:3)
  IF(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).EQ.0) THEN
    CircleArea = SurfMeshSubSideData(iSample,jSample,BCSideID)%area
  ELSE
    CircleArea = 0.
    counter = 0
    DO iMC = 1,MCVar
      CALL RANDOM_NUMBER(RandVal2)
      Particle_pos(1:3) = TriaNode(1:3) + Vector1(1:3) * RandVal2(1)
      Particle_pos(1:3) = Particle_pos(1:3) + Vector2(1:3) * RandVal2(2)
      PartDistance = ndist(1)*(Particle_pos(1)-midpoint(1)) & !Distance from v1-v2
                    + ndist(2)*(Particle_pos(2)-midpoint(2)) &
                    + ndist(3)*(Particle_pos(3)-midpoint(3))
      IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
        Particle_pos(1:3) = 2*midpoint(1:3)-Particle_pos(1:3)
      END IF
      point(1)=Particle_pos(dir(2))-origin(1)
      point(2)=Particle_pos(dir(3))-origin(2)
      radius=SQRT( (point(1))**2+(point(2))**2 )
      IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
        CircleArea = CircleArea + 1./REAL(MCVar)
        counter = counter + 1
      END IF
    END DO
    CircleArea = CircleArea * SurfMeshSubSideData(iSample,jSample,BCSideID)%area
  END IF
  Species(iSpec)%Surfaceflux(iSF)%CircleAreaPerTriaSide(iSample,jSample,iSide) = CircleArea
END DO; END DO

END SUBROUTINE CircularInflow_Area


END MODULE MOD_surface_flux
