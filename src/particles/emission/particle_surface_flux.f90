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


SUBROUTINE InitializeParticleSurfaceflux()
!===================================================================================================================================
! Init Particle Inserting via Surface Flux
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: PI, BoltzmannConst
USE MOD_ReadInTools
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound,nPartBound, nAdaptiveBC, nPorousBC
USE MOD_Particle_Vars          ,ONLY: Species, nSpecies, DoSurfaceFlux, DoPoissonRounding, nDataBC_CollectCharges, DoTimeDepInflow, &
                                     Adaptive_MacroVal, MacroRestartData_tmp, AdaptiveWeightFac, VarTimeStep
USE MOD_PARTICLE_Vars          ,ONLY: nMacroRestartFiles, UseAdaptive, UseCircularInflow
USE MOD_Particle_Vars          ,ONLY: DoForceFreeSurfaceFlux
USE MOD_DSMC_Vars              ,ONLY: useDSMC, BGGas
USE MOD_Mesh_Vars              ,ONLY: nBCSides, BC, SideToElem, NGeo, nElems, offsetElem
USE MOD_Particle_Surfaces_Vars ,ONLY: BCdata_auxSF, BezierSampleN, SurfMeshSubSideData, SurfMeshSideAreas
USE MOD_Particle_Surfaces_Vars ,ONLY: SurfFluxSideSize, TriaSurfaceFlux, WriteTriaSurfaceFluxDebugMesh, SideType
USE MOD_Particle_Surfaces      ,ONLY: GetBezierSampledAreas, GetSideBoundingBox, CalcNormAndTangTriangle
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Particle_Mesh          ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking, DoRefMapping
USE MOD_IO_HDF5
USE MOD_HDF5_INPUT             ,ONLY: DatasetExists,ReadAttribute,ReadArray,GetDataSize
USE MOD_Restart_Vars           ,ONLY: DoRestart,RestartFile
USE MOD_Particle_Vars          ,ONLY: Symmetry2D, Symmetry2DAxisymmetric
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_DSMC_Symmetry2D        ,ONLY: DSMC_2D_CalcSymmetryArea, DSMC_2D_CalcSymmetryAreaSubSides
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
INTEGER               :: iPartBound,iSpec,iSF,SideID,BCSideID,iSide,ElemID,iLocSide,iSample,jSample,iBC,currentBC,iCount
INTEGER               :: iCopy1, iCopy2, iCopy3, nSides
CHARACTER(32)         :: hilf, hilf2, hilf3
REAL                  :: a, vSF, projFak, v_thermal
REAL                  :: vec_nIn(3), nVFR, vec_t1(3), vec_t2(3), point(2)
LOGICAL               :: noAdaptive
INTEGER               :: MaxSurfacefluxBCs
INTEGER               :: nDataBC                             ! number of different PartBounds used for SFs
INTEGER,ALLOCATABLE   :: TmpMapToBC(:)                       ! PartBC
INTEGER,ALLOCATABLE   :: TmpSideStart(:)                     ! Start of Linked List for Sides in SurfacefluxBC
INTEGER,ALLOCATABLE   :: TmpSideNumber(:)                    ! Number of Particles in Sides in SurfacefluxBC
INTEGER,ALLOCATABLE   :: TmpSideEnd(:)                       ! End of Linked List for Sides in SurfacefluxBC
INTEGER,ALLOCATABLE   :: TmpSideNext(:)                      ! Next Side in same SurfacefluxBC (Linked List)
INTEGER,ALLOCATABLE   :: nType0(:,:), nType1(:,:), nType2(:,:)
REAL                  :: totalArea
REAL,ALLOCATABLE      :: tmp_SubSideAreas(:,:), tmp_SubSideDmax(:,:)
REAL,ALLOCATABLE      :: tmp_Vec_nOut(:,:,:), tmp_Vec_t1(:,:,:), tmp_Vec_t2(:,:,:)
REAL,ALLOCATABLE      :: tmp_BezierControlPoints2D(:,:,:,:,:)
REAL,DIMENSION(1:3,1:8):: BoundingBox
INTEGER,ALLOCATABLE   :: Adaptive_BC_Map(:), tmp_Surfaceflux_BCs(:)
LOGICAL,ALLOCATABLE   :: Adaptive_Found_Flag(:)
INTEGER               :: nAdaptive_Found, iSS, nSurffluxBCs_old, nSurffluxBCs_new, iSFx
REAL,ALLOCATABLE      :: sum_pressurefraction(:)
REAL                  :: Vector1(3),Vector2(3),Vector3(3)
INTEGER               :: dir(3)
REAL                  :: origin(2),xyzNod(3)
REAL                  :: corner(3)
REAL                  :: VecBoundingBox(3)
INTEGER               :: iNode
REAL                  :: vec(2)
REAL                  :: radiusCorner(2,4)
LOGICAL               :: r0inside, intersecExists(2,2)
REAL                  :: corners(2,4),rmin,rmax!,atan2Shift
INTEGER               :: FileID
LOGICAL               :: OutputSurfaceFluxLinked
REAL,ALLOCATABLE      :: ElemData_HDF5(:,:,:)
LOGICAL               :: AdaptiveDataExists, AdaptiveInitDone
INTEGER               :: iElem
#if USE_MPI
INTEGER               :: iProc
REAL, ALLOCATABLE     :: areasLoc(:),areasGlob(:)
REAL                  :: totalAreaSF_global
#endif
REAL                  :: ymin, ymax, VFR_total, yMinTemp, yMaxTemp
TYPE tBCdata_auxSFTemp
  REAL, ALLOCATABLE   :: SubSideWeight(:,:)
  REAL, ALLOCATABLE   :: WeightingFactor(:)
  REAL, ALLOCATABLE   :: SubSideArea(:,:)
END TYPE
TYPE(tBCdata_auxSFTemp),ALLOCATABLE          :: BCdata_auxSFTemp(:)
INTEGER               :: iSub
!===================================================================================================================================

#if USE_MPI
CALL MPI_BARRIER(PartMPI%COMM,iError)
#endif /*USE_MPI*/
OutputSurfaceFluxLinked=GETLOGICAL('OutputSurfaceFluxLinked','.FALSE.')

! global calculations for sampling the faces for area and vector calculations (checks the integration with CODE_ANALYZE)
ALLOCATE (tmp_SubSideAreas(SurfFluxSideSize(1),SurfFluxSideSize(2)), &
  tmp_Vec_nOut(3,SurfFluxSideSize(1),SurfFluxSideSize(2)), &
  tmp_Vec_t1(3,SurfFluxSideSize(1),SurfFluxSideSize(2)), &
  tmp_Vec_t2(3,SurfFluxSideSize(1),SurfFluxSideSize(2)), &
  SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:nBCSides)    )
IF (.NOT.TriaSurfaceFlux) THEN
  ALLOCATE (tmp_SubSideDmax(SurfFluxSideSize(1),SurfFluxSideSize(2)), &
    tmp_BezierControlPoints2D(2,0:NGeo,0:NGeo,SurfFluxSideSize(1),SurfFluxSideSize(2)) )
END IF
ALLOCATE(SurfMeshSideAreas(1:nBCSides))
SurfMeshSideAreas=0.
totalArea=0.
DO BCSideID=1,nBCSides
  ElemID = SideToElem(1,BCSideID)
  IF (ElemID.LT.1) THEN !not sure if necessary
    ElemID = SideToElem(2,BCSideID)
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
        ,TriNum=jSample,ElemID_opt=ElemID,LocSideID_opt=ilocSide)
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

UseCircularInflow=.FALSE.
UseAdaptive=.FALSE.
MaxSurfacefluxBCs=0
nDataBC=nDataBC_CollectCharges !sides may be also used for collectcharges of floating potential!!!
DoSurfaceFlux=.FALSE.
!-- 0.: allocate and initialize aux. data of BCs (SideLists for Surfacefluxes):
!-----moved to end of InitializeVariables in particle_init!!!
!ALLOCATE(BCdata_auxSF(1:nPartBound))
!DO iPartBound=1,nPartBound
!  BCdata_auxSF(iPartBound)%SideNumber=-1 !init value when not used
!END DO

AdaptiveWeightFac = GETREAL('Part-AdaptiveWeightingFactor','0.001')
! auxiliary arrays for defining all Adaptive_BCs
IF (nAdaptiveBC.GT.0) THEN
  ALLOCATE(Adaptive_BC_Map(1:nAdaptiveBC))
  Adaptive_BC_Map(:)=0
  ALLOCATE(Adaptive_Found_Flag(1:nAdaptiveBC))
  iSS = 0
  DO iPartBound=1,nPartBound
    IF(PartBound%Adaptive(iPartBound))THEN
      iSS = iSS + 1
      Adaptive_BC_Map(iSS) = iPartBound
    END IF
  END DO
  Adaptive_Found_Flag(:) = .FALSE.
  ALLOCATE(sum_pressurefraction(1:nAdaptiveBC))
  sum_pressurefraction(:) = 0.
END IF

!-- 1.: read/prepare parameters and determine nec. BCs
DO iSpec=1,nSpecies
  IF (nAdaptiveBC.GT.0) THEN
    Adaptive_Found_Flag(:) = .FALSE.
  END IF
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  Species(iSpec)%nSurfacefluxBCs = GETINT('Part-Species'//TRIM(hilf)//'-nSurfacefluxBCs','0')
  IF (useDSMC) THEN
    IF (BGGas%BGGasSpecies.EQ.iSpec) THEN
      IF (Species(iSpec)%nSurfacefluxBCs.GT.0 .OR. nAdaptiveBC.GT.0) CALL abort(&
__STAMP__&
, 'SurfaceFlux or AdaptiveBCs are not implemented for the BGG-species!')
    END IF
  END IF
  ! if no surfacefluxes defined and only adaptive boundaries then first allocation with adaptive
  IF ((Species(iSpec)%nSurfacefluxBCs.EQ.0) .AND. (nAdaptiveBC.GT.0)) THEN
    Species(iSpec)%nSurfacefluxBCs = nAdaptiveBC
    ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      Species(iSpec)%Surfaceflux(iSF)%BC = Adaptive_BC_Map(iSF)
    END DO
    nAdaptive_Found = nAdaptiveBC
    Adaptive_Found_Flag(:) = .TRUE.
  ! if no surfaceflux needed
  ELSE IF ((Species(iSpec)%nSurfacefluxBCs.EQ.0) .AND. (nAdaptiveBC.EQ.0)) THEN
    CYCLE
  ELSE
    ALLOCATE(Species(iSpec)%Surfaceflux(1:Species(iSpec)%nSurfacefluxBCs))
    ! Initialize Surfaceflux to BC mapping and check if defined Surfacefluxes from init overlap with Adaptive BCs
    Species(iSpec)%Surfaceflux(:)%BC=-1
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      WRITE(UNIT=hilf2,FMT='(I0)') iSF
      hilf2=TRIM(hilf)//'-Surfaceflux'//TRIM(hilf2)
      Species(iSpec)%Surfaceflux(iSF)%BC = GETINT('Part-Species'//TRIM(hilf2)//'-BC','0')
      IF (nAdaptiveBC.GT.0) THEN
        DO iSS=1,nAdaptiveBC
          IF (Adaptive_BC_Map(iSS).EQ.Species(iSpec)%Surfaceflux(iSF)%BC) THEN
            Adaptive_Found_Flag(iSS) = .TRUE.
          END IF
        END DO
      END IF
    END DO
    nAdaptive_Found = 0
    DO iSS=1,nAdaptiveBC
      IF(Adaptive_Found_Flag(iSS)) nAdaptive_Found = nAdaptive_Found + 1
    END DO
  END IF
  ! add missing Adaptive BCs at end of Surfaceflux array and reduce number of constant surfacefluxBC
  ! additionally rearrange surfaceflux array for Adaptive BC being the last entries
  IF (nAdaptiveBC.GT.0) THEN
    ALLOCATE(tmp_Surfaceflux_BCs(1:Species(iSpec)%nSurfacefluxBCs))
    tmp_Surfaceflux_BCs(:) = Species(iSpec)%Surfaceflux(:)%BC
    nSurffluxBCs_old = Species(iSpec)%nSurfacefluxBCs
    Species(iSpec)%nSurfacefluxBCs = Species(iSpec)%nSurfacefluxBCs - nAdaptive_Found
    nSurffluxBCs_new = Species(iSpec)%nSurfacefluxBCs + nAdaptiveBC
    DEALLOCATE(Species(iSpec)%Surfaceflux)
    ALLOCATE(Species(iSpec)%Surfaceflux(1:nSurffluxBCs_new))
    iSFx = 1
    DO iSF=1,nSurffluxBCs_old
      IF (PartBound%Adaptive(tmp_Surfaceflux_BCs(iSF))) CYCLE
      Species(iSpec)%Surfaceflux(iSFx)%BC = tmp_Surfaceflux_BCs(iSF)
      iSFx = iSFx +1
    END DO
    DO iSFx=1,nAdaptiveBC
      Species(iSpec)%Surfaceflux(Species(iSpec)%nSurfacefluxBCs+iSFx)%BC = Adaptive_BC_Map(iSFx)
    END DO
    DEALLOCATE(tmp_Surfaceflux_BCs)
  END IF

  MaxSurfacefluxBCs=MAX(MaxSurfacefluxBCs,Species(iSpec)%nSurfacefluxBCs)
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
    IF (iSF .LE. Species(iSpec)%nSurfacefluxBCs) THEN
      noAdaptive=.TRUE.
    ELSE
      noAdaptive=.FALSE.
    END IF
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
    IF (noAdaptive) THEN
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
      Species(iSpec)%Surfaceflux(iSF)%VeloIC                = GETREAL('Part-Species'//TRIM(hilf2)//'-VeloIC','0.')
      Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal          = GETLOGICAL('Part-Species'//TRIM(hilf2)//'-VeloIsNormal','.FALSE.')
      IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
        Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit=.FALSE.
        Species(iSpec)%Surfaceflux(iSF)%CircularInflow=.FALSE.
      ELSE
        Species(iSpec)%Surfaceflux(iSF)%VeloVecIC          =GETREALARRAY('Part-Species'//TRIM(hilf2)//'-VeloVecIC',3,'1. , 0. , 0.')
        Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit=GETLOGICAL('Part-Species'//TRIM(hilf2)//'-SimpleRadialVeloFit','.FALSE.')
        Species(iSpec)%Surfaceflux(iSF)%CircularInflow=GETLOGICAL('Part-Species'//TRIM(hilf2)//'-CircularInflow','.FALSE.')
        IF (Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit) THEN
          Species(iSpec)%Surfaceflux(iSF)%CircularInflow =.TRUE.
          Species(iSpec)%Surfaceflux(iSF)%preFac       = GETREAL('Part-Species'//TRIM(hilf2)//'-preFac','0.')
          Species(iSpec)%Surfaceflux(iSF)%powerFac     = GETREAL('Part-Species'//TRIM(hilf2)//'-powerFac','0.')
          Species(iSpec)%Surfaceflux(iSF)%shiftFac     = GETREAL('Part-Species'//TRIM(hilf2)//'-shiftFac','0.')
        END IF !Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit
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
    ELSE !Adaptive
      Species(iSpec)%Surfaceflux(iSF)%velocityDistribution  = Species(iSpec)%Init(0)%velocityDistribution
      IF (PartBound%AdaptiveMacroRestartFileID(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.0 .OR. nMacroRestartFiles.EQ.0) THEN
        Species(iSpec)%Surfaceflux(iSF)%VeloIC                = Species(iSpec)%Init(0)%VeloIC
        Species(iSpec)%Surfaceflux(iSF)%VeloVecIC             = Species(iSpec)%Init(0)%VeloVecIC
      END IF
      Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal          = .FALSE.
      Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit   = .FALSE.
      Species(iSpec)%Surfaceflux(iSF)%CircularInflow=.FALSE.
    END IF
    IF (.NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .OR. .NOT.noAdaptive) THEN
      !--- normalize VeloVecIC
      IF (PartBound%AdaptiveMacroRestartFileID(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.0 .OR. nMacroRestartFiles.EQ.0) THEN
        IF (.NOT. ALL(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(:).eq.0.)) THEN
          Species(iSpec)%Surfaceflux(iSF)%VeloVecIC = Species(iSpec)%Surfaceflux(iSF)%VeloVecIC &
            /SQRT(DOT_PRODUCT(Species(iSpec)%Surfaceflux(iSF)%VeloVecIC,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC))
        END IF
      END IF
    END IF
    IF (noAdaptive) THEN
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
    ELSE !Adaptive
      WRITE(UNIT=hilf3,FMT='(I0)') Adaptive_BC_Map(iSF-Species(iSpec)%nSurfacefluxBCs)
      Species(iSpec)%Surfaceflux(iSF)%PressureFraction      = &
        GETREAL('Part-Boundary'//TRIM(hilf3)//'-Species'//TRIM(hilf)//'-Pressurefraction','0.')
      sum_pressurefraction(iSF-Species(iSpec)%nSurfacefluxBCs) = sum_pressurefraction(iSF-Species(iSpec)%nSurfacefluxBCs) &
        + Species(iSpec)%Surfaceflux(iSF)%PressureFraction
      IF (PartBound%AdaptiveMacroRestartFileID(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.0 &
          .OR. PartBound%AdaptiveType(Species(iSpec)%Surfaceflux(iSF)%BC).EQ.1 .OR. nMacroRestartFiles.EQ.0) THEN
        Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC       = PartBound%AdaptiveTemp(Species(iSpec)%Surfaceflux(iSF)%BC)
        Species(iSpec)%Surfaceflux(iSF)%PartDensity           = Species(iSpec)%Surfaceflux(iSF)%PressureFraction &
          * PartBound%AdaptivePressure(Species(iSpec)%Surfaceflux(iSF)%BC) &
          / (BoltzmannConst * Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC)
      END IF
      Species(iSpec)%Surfaceflux(iSF)%ReduceNoise           = .FALSE.
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
      IF (noAdaptive) THEN
        WRITE( hilf3, '(I0.2)') NGeo*NGeo*NGeo !1 for linear elements, this is an arbitray estimation for higher N!
        Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = GETINT('Part-Species'//TRIM(hilf2)//'-ARM_DmaxSampleN',hilf3)
      ELSE !Adaptive
        Species(iSpec)%Surfaceflux(iSF)%ARM_DmaxSampleN = NGeo*NGeo*NGeo
      END IF
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
IF (nAdaptiveBC.GT.0) THEN
  IF( (MINVAL(sum_pressurefraction(:)).LT.0.99).OR.(MAXVAL(sum_pressurefraction(:)).GT.1.01) ) CALL abort( &
__STAMP__&
, 'Sum of all pressurefractions .NE. 1')
END IF

SDEALLOCATE(Adaptive_BC_Map)
SDEALLOCATE(Adaptive_Found_Flag)

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoPoissonRounding,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DoTimeDepInflow,1,MPI_LOGICAL,MPI_LAND,PartMPI%COMM,iError) !set T if this is for all procs
#endif /*USE_MPI*/

!-- 2.: create Side lists for applicable BCs
!--- 2a: temporary (linked) lists
ALLOCATE(TmpMapToBC(1:nDataBC) &
        ,TmpSideStart(1:nDataBC) &
        ,TmpSideNumber(1:nDataBC) &
        ,TmpSideEnd(1:nDataBC) &
        ,TmpSideNext(1:nBCSides)) !Next: Sides of diff. BCs ar not overlapping!
TmpMapToBC = 0
TmpSideStart = 0
TmpSideNumber = 0
TmpSideEnd = 0
TmpSideNext = 0
nDataBC=0
DO iBC=1,nPartBound
  IF (BCdata_auxSF(iBC)%SideNumber.EQ. -1) CYCLE !not set for SFs or CollectCharges
  nDataBC=nDataBC+1
  TmpMapToBC(nDataBC)=iBC
END DO
DO BCSideID=1,nBCSides
  currentBC=0
  DO iBC=1,nDataBC
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
IF (UseCircularInflow) THEN
  ALLOCATE(nType0(1:MaxSurfacefluxBCs,1:nSpecies), &
    nType1(1:MaxSurfacefluxBCs,1:nSpecies), &
    nType2(1:MaxSurfacefluxBCs,1:nSpecies) )
  nType0=0
  nType1=0
  nType2=0
END IF

IF(RadialWeighting%DoRadialWeighting) THEN
  ALLOCATE(BCdata_auxSFTemp(1:nPartBound))
END IF

!--- 2b: save sequential lists in BCdata_auxSF
DO iBC=1,nDataBC
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
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
      IF (TmpMapToBC(iBC).EQ.Species(iSpec)%Surfaceflux(iSF)%BC) THEN !only surfacefluxes with iBC
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),1:TmpSideNumber(iBC)) )
        IF (UseCircularInflow .AND. (iSF .LE. Species(iSpec)%nSurfacefluxBCs)) THEN
          ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(1:TmpSideNumber(iBC)) )
        END IF
      END IF
      IF(RadialWeighting%DoRadialWeighting) THEN
        ALLOCATE(Species(iSpec)%Surfaceflux(iSF)%nVFRSub(1:TmpSideNumber(iBC),1:RadialWeighting%nSubSides))
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
                BCdata_auxSFTemp(TmpMapToBC(iBC))%WeightingFactor(iCount) = (1. + GEO%ElemMidPoint(2,ElemID)&
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
          ,TriNum=jSample,ElemID_opt=ElemID,LocSideID_opt=ilocSide)
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

DEALLOCATE(TmpMapToBC &
          ,TmpSideStart &
          ,TmpSideNumber &
          ,TmpSideEnd &
          ,TmpSideNext)

!-- 3.: initialize Surfaceflux-specific data
! Allocate sampling of near adaptive boundary element values
IF((nAdaptiveBC.GT.0).OR.UseAdaptive.OR.(nPorousBC.GT.0))THEN
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
END IF

DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
    IF (iSF .LE. Species(iSpec)%nSurfacefluxBCs) THEN
      noAdaptive=.TRUE.
    ELSE
      noAdaptive=.FALSE.
    END IF
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
        !-- check where the sides are located relative to rmax (based on corner nodes of bounding box)
        !- RejectType=0 : complete side is inside valid bounds
        !- RejectType=1 : complete side is outside of valid bounds
        !- RejectType=2 : side is partly inside valid bounds
        IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
          CALL GetSideBoundingBox(BCSideID,BoundingBox)
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
            nType1(iSF,iSpec)=nType1(iSF,iSpec)+1
          ELSE IF ( (rmax .LE. Species(iSpec)%Surfaceflux(iSF)%rmax) .AND. (rmin .GE. Species(iSpec)%Surfaceflux(iSF)%rmin) ) THEN
            Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=0
            nType0(iSF,iSpec)=nType0(iSF,iSpec)+1
          ELSE
            Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide)=2
            nType2(iSF,iSpec)=nType2(iSF,iSpec)+1
          END IF !  (rmin > Surfaceflux-rmax) .OR. (rmax < Surfaceflux-rmin)
          IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
            IF((Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide).NE.1)   &
                .AND.(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4)) CALL CircularInflow_Area(iSpec,iSF,iSide,BCSideID)
          END IF
        END IF ! CircularInflow: check r-bounds
        IF (noAdaptive.AND.(.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive)) THEN
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
        END IF
        DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
          IF (Species(iSpec)%Surfaceflux(iSF)%AcceptReject) THEN
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
          END IF
        END DO; END DO !jSample=1,SurfFluxSideSize(2); iSample=1,SurfFluxSideSize(1)
        IF ((.NOT.noAdaptive).OR.(Species(iSpec)%Surfaceflux(iSF)%Adaptive)) THEN
          IF (.NOT.AdaptiveInitDone) THEN
            ! initialize velocity, trans_temperature and density of macrovalues
            IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
              ! ofc, a file for every macrovalue...dirty fix
              FileID = Species(iSpec)%Init(0)%ElemPartDensityFileID
            ELSE
              FileID = PartBound%AdaptiveMacroRestartFileID(Species(iSpec)%Surfaceflux(iSF)%BC)
            END IF
            IF (FileID.GT.0 .AND. FileID.LE.nMacroRestartFiles) THEN
              Adaptive_MacroVal(DSMC_VELOX,ElemID,iSpec) = MacroRestartData_tmp(DSMC_VELOX,ElemID,iSpec,FileID)
              Adaptive_MacroVal(DSMC_VELOY,ElemID,iSpec) = MacroRestartData_tmp(DSMC_VELOY,ElemID,iSpec,FileID)
              Adaptive_MacroVal(DSMC_VELOZ,ElemID,iSpec) = MacroRestartData_tmp(DSMC_VELOZ,ElemID,iSpec,FileID)
              Adaptive_MacroVal(DSMC_TEMPX,ElemID,iSpec) = MAX(0.,MacroRestartData_tmp(DSMC_TEMPX,ElemID,iSpec,FileID))
              Adaptive_MacroVal(DSMC_TEMPY,ElemID,iSpec) = MAX(0.,MacroRestartData_tmp(DSMC_TEMPY,ElemID,iSpec,FileID))
              Adaptive_MacroVal(DSMC_TEMPZ,ElemID,iSpec) = MAX(0.,MacroRestartData_tmp(DSMC_TEMPZ,ElemID,iSpec,FileID))
              Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec) = MacroRestartData_tmp(DSMC_NUMDENS,ElemID,iSpec,FileID)
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
        END IF
      END DO ! iSide

    ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
      CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
    END IF
#ifdef CODE_ANALYZE
    IF (BCdata_auxSF(currentBC)%SideNumber.GT.0 .AND. Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      IPWRITE(*,'(I4,A,2(x,I0),A,3(x,I0))') ' For Surfaceflux/Spec',iSF,iSpec,' are nType0,1,2: ' &
                                            , nType0(iSF,iSpec),nType1(iSF,iSpec),nType2(iSF,iSpec)
    END IF
#endif /*CODE_ANALYZE*/

    !--- 3b: ReduceNoise initialization
    IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
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
    END IF !ReduceNoise
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
  END DO !iSF
END DO !iSpec

! Deallocate auxiliary variable container (no pointers used inside container)
IF(RadialWeighting%DoRadialWeighting.AND.ALLOCATED(BCdata_auxSFTemp)) DEALLOCATE(BCdata_auxSFTemp)

!-- write debug-mesh for tria-surfflux
IF (WriteTriaSurfaceFluxDebugMesh) THEN
  !count sides
  nSides=0
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC !go through sides if present in proc...
      IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
        nSides=nSides+BCdata_auxSF(currentBC)%SideNumber
      ELSE IF (BCdata_auxSF(currentBC)%SideNumber.EQ.-1) THEN
        CALL abort(&
  __STAMP__&
  ,'ERROR in ParticleSurfaceflux: Someting is wrong with SideNumber of BC ',currentBC)
      END IF
    END DO !iSF
  END DO !iSpec
  WRITE(UNIT=hilf,FMT='(I4.4)') myRank
  OPEN(UNIT   = 103, &
         FILE   = 'Tria-Surfflux-debugmesh_'//TRIM(hilf)//'.tec' ,&
         STATUS = 'UNKNOWN')
  WRITE(103,*) 'TITLE="Tria-Surfflux-debugmesh" '
  WRITE(103,'(102a)') 'VARIABLES ="x","y","z","BC","iSF","iSpec"'
  WRITE(103,*) 'ZONE NODES=',4*nSides,', ELEMENTS=',2*nSides,'DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
  ! Write nodes
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC !go through sides if present in proc...
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
          !WRITE(103,'(3(F0.10,1X),3(I0,1X))')GEO%NodeCoords(1:3,1,iLocSide,ElemID),currentBC,iSF,iSpec
          !WRITE(103,'(3(F0.10,1X),3(I0,1X))')GEO%NodeCoords(1:3,2,iLocSide,ElemID),currentBC,iSF,iSpec
          !WRITE(103,'(3(F0.10,1X),3(I0,1X))')GEO%NodeCoords(1:3,3,iLocSide,ElemID),currentBC,iSF,iSpec
          !WRITE(103,'(3(F0.10,1X),3(I0,1X))')GEO%NodeCoords(1:3,4,iLocSide,ElemID),currentBC,iSF,iSpec
          WRITE(103,'(3(F0.10,1X),3(I0,1X))')BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod,currentBC,iSF,iSpec
          WRITE(103,'(3(F0.10,1X),3(I0,1X))')BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod &
            +BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,1),currentBC,iSF,iSpec
          WRITE(103,'(3(F0.10,1X),3(I0,1X))')BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod &
            +BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,2),currentBC,iSF,iSpec
          WRITE(103,'(3(F0.10,1X),3(I0,1X))')BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod &
            +BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%Vectors(:,3),currentBC,iSF,iSpec
        END DO ! iSide
      END IF
    END DO !iSF
  END DO !iSpec
  ! Write sides
  nSides=0
  DO iSpec=1,nSpecies
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC !go through sides if present in proc...
      IF (BCdata_auxSF(currentBC)%SideNumber.GT.0) THEN
        DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
          WRITE(103,'(4(I0,1X))')nSides*4+1,nSides*4+2,nSides*4+3,nSides*4+3 !1. tria
          WRITE(103,'(4(I0,1X))')nSides*4+1,nSides*4+3,nSides*4+4,nSides*4+4 !2. tria
          nSides=nSides+1
        END DO ! iSide
      END IF
    END DO !iSF
  END DO !iSpec
  CLOSE(103)
END IF !TriaSurfaceFlux

! Setting variables required after a restart
IF(DoRestart) THEN
  DO iSpec=1,nSpecies
    DO iSF = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
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

#if USE_MPI
CALL MPI_BARRIER(PartMPI%COMM,iError)
#endif /*USE_MPI*/
END SUBROUTINE InitializeParticleSurfaceflux


SUBROUTINE ParticleSurfaceflux()
!===================================================================================================================================
! Particle Inserting via Surface Flux and (if present) adaptiveBC (Surface Flux adapting part density, velocity or temperature)
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars       ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools              ,ONLY: UpdateNextFreePosition,GetParticleWeight
USE MOD_MacroBody_vars          ,ONLY: UseMacroBody
USE MOD_MacroBody_tools         ,ONLY: INSIDEMACROBODY
USE MOD_DSMC_Vars               ,ONLY: useDSMC, CollisMode, SpecDSMC, DSMC, PartStateIntEn, RadialWeighting
USE MOD_SurfaceModel_Vars       ,ONLY: SurfModel
USE MOD_Particle_Boundary_Tools ,ONLY: CalcWallSample
USE MOD_DSMC_Init               ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_DSMC_Symmetry2D         ,ONLY: CalcRadWeightMPF
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfMesh, PartBound, nAdaptiveBC, nSurfSample
USE MOD_TimeDisc_Vars           ,ONLY: TEnd, time
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance,CalcMassflowRate
USE MOD_Particle_Analyze_Vars   ,ONLY: nPartIn,PartEkinIn
USE MOD_Timedisc_Vars           ,ONLY: RKdtFrac,RKdtFracTotal,Time
USE MOD_Particle_Analyze_Tools  ,ONLY: CalcEkinPart
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem!, ElemBaryNGeo
USE MOD_Mesh_Vars               ,ONLY: NGeo!,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo
USE MOD_Particle_Mesh_Vars      ,ONLY: GEO
USE MOD_Particle_Mesh           ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfMeshSubSideData
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Particle_Tracking_Vars  ,ONLY: TriaTracking
#if defined(IMPA) || defined(ROS)
USE MOD_Particle_Tracking_Vars  ,ONLY: DoRefMapping
#endif /*IMPA*/
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D,BezierSampleXi,SurfFluxSideSize,TriaSurfaceFlux
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_Eval_xyz                ,ONLY: GetPositionInRefElem
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars  ,ONLY: PartOut, MPIRankOut
#if  defined(IMPA) || defined(ROS)
USE MOD_Timedisc_Vars           ,ONLY: iStage,nRKStages
#endif
#endif /*CODE_ANALYZE*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: nSurfacefluxPerElem
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_part_emission_tools     ,ONLY: IntegerDivide,SetParticleChargeAndMass,SetParticleMPF,SamplePoissonDistri
USE MOD_part_pos_and_velo       ,ONLY: SetParticleVelocity
#if CODE_ANALYZE
USE MOD_part_emission_tools     ,ONLY: CalcVectorAdditionCoeffs
#endif /*CODE_ANALYZE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                     :: iSpec , PositionNbr, iSF, iSide, currentBC, SideID, iLoop
INTEGER                     :: NbrOfParticle, ExtraParts
INTEGER                     :: BCSideID, ElemID, iLocSide, iSample, jSample, PartInsSF, PartInsSubSide, iPart, iPartTotal, IntSample
INTEGER                     :: ParticleIndexNbr, allocStat
REAL                        :: PartIns,VFR_total
REAL                        :: Particle_pos(3), RandVal1, RandVal2(2), xNod,yNod,zNod, Vector2D(2), RVec(2), minPos(2)
INTEGER                     :: minVec
REAL,ALLOCATABLE            :: particle_positions(:), particle_xis(:)
INTEGER(KIND=8)             :: inserted_Particle_iter,inserted_Particle_time,inserted_Particle_diff
INTEGER,ALLOCATABLE         :: PartInsProc(:),PartInsSubSides(:,:,:)
REAL                        :: xiab(1:2,1:2),xi(2),E,F,G,D,gradXiEta2D(1:2,1:2),gradXiEta3D(1:2,1:3)
REAL                        :: point(2),origin(2),veloR,vTot,phi,radius,preFac,powerFac,shiftFac
INTEGER                     :: dir(3), nReject, allowedRejections
LOGICAL                     :: AcceptPos, noAdaptive
!variables used for sampling of of energies and impulse of emitted particles from surfaces
INTEGER                     :: PartsEmitted
REAL                        :: TransArray(1:6),IntArray(1:6)
REAL                        :: VelXold, VelYold, VelZold, VeloReal
REAL                        :: EtraOld, EtraWall, EtraNew
REAL                        :: ErotOld, ErotWall, ErotNew
REAL                        :: EvibOld, EvibWall, EVibNew
REAL                        :: Vector1(3),Vector2(3),PartDistance,ndist(3),midpoint(3),AreasTria(2)
INTEGER                     :: p,q,SurfSideID,PartID,Node1,Node2,ExtraPartsTria(2)
REAL                        :: ElemPartDensity, VeloVec(1:3), VeloIC
REAL                        :: VeloVecIC(1:3), ProjFak, v_thermal, a, T, vSF, nVFR,vec_nIn(1:3), pressure, veloNormal
REAL, PARAMETER             :: eps_nontria=1.0E-6 !prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
#if USE_LOADBALANCE
! load balance
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
#ifdef CODE_ANALYZE
REAL                        :: tmpVec(3)
#endif /*CODE_ANALYZE*/
TYPE(tSurfFluxLink),POINTER :: currentSurfFluxPart => NULL()
REAL                        :: PminTemp, PmaxTemp
INTEGER                     :: PartInsSideSubSub(1:RadialWeighting%nSubSides), iSub, iPartSub
!===================================================================================================================================

DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs+nAdaptiveBC
    PartsEmitted = 0
    IF (iSF .LE. Species(iSpec)%nSurfacefluxBCs) THEN
      noAdaptive=.TRUE.
    ELSE
      noAdaptive=.FALSE.
    END IF
    currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
    NbrOfParticle = 0 ! calculated within (sub)side-Loops!
    iPartTotal=0

    IF(CalcMassflowRate) THEN
    ! Reset the mass flow rate counter for the next time step
      Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = 0.
    END IF

    IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN
      dir   =Species(iSpec)%Surfaceflux(iSF)%dir
      origin=Species(iSpec)%Surfaceflux(iSF)%origin
      IF (Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit) THEN
        preFac=Species(iSpec)%Surfaceflux(iSF)%preFac
        powerFac=Species(iSpec)%Surfaceflux(iSF)%powerFac
        shiftFac=Species(iSpec)%Surfaceflux(iSF)%shiftFac
      END IF
    END IF
    IF(Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
      IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) CALL AdaptiveBoundary_ConstMassflow_Weight(iSpec,iSF)
    END IF
    !--- Noise reduction (both ReduceNoise=T (with comm.) and F (proc local), but not for DoPoissonRounding)
    IF (.NOT.DoPoissonRounding .AND. .NOT. DoTimeDepInflow .AND. noAdaptive) THEN
      IF (Species(iSpec)%Surfaceflux(iSF)%ReduceNoise) THEN
        !-- calc global to-be-inserted number of parts and distribute to procs (root)
        SDEALLOCATE(PartInsProc)
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
    END IF !.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow .AND. noAdaptive

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
      SideID=GetGlobalNonUniqueSideID(offsetElem+ElemID,iLocSide)
      IF (TriaSurfaceFlux) THEN
        xNod = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(1)
        yNod = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(2)
        zNod = BCdata_auxSF(currentBC)%TriaSideGeo(iSide)%xyzNod(3)
      END IF
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
        IF (noAdaptive) THEN
          IF ( PartBound%Reactive(currentBC) )  THEN
            IF (SurfMesh%SideIDToSurfID(SideID).GT.0) THEN
              ! sumEvapPart for triatracking is only allocated over (1,1,nsurfsides,nspecies)
              IF (.NOT.TriaSurfaceFlux .OR. (iSample.EQ.1 .AND. jSample.EQ.1)) THEN
                ExtraParts = SurfModel%SumEvapPart(iSample,jSample,SurfMesh%SideIDToSurfID(SideID),iSpec)
                SurfModel%SumEvapPart(iSample,jSample,SurfMesh%SideIDToSurfID(SideID),iSpec) = 0
              END IF
              IF (TriaSurfaceFlux) THEN
                IF (iSample.EQ.1 .AND. jSample.EQ.1) THEN !first tria
                  AreasTria(1)=SurfMeshSubSideData(1,1,BCSideID)%area
                  AreasTria(2)=SurfMeshSubSideData(SurfFluxSideSize(1),SurfFluxSideSize(2),BCSideID)%area
                  ExtraPartsTria(:) = 0
                  CALL IntegerDivide(ExtraParts, 2, AreasTria, ExtraPartsTria)
                  ExtraParts = ExtraPartsTria(1)
                ELSE !second tria
                  ExtraParts = ExtraPartsTria(2)
                END IF
              END IF !TriaSurfaceFlux
              SurfModel%Info(iSpec)%NumOfDes=SurfModel%Info(iSpec)%NumOfDes+ExtraParts
            END IF !SurfMesh%SideIDToSurfID(SideID).GT.0
          END IF !reactive surface
        END IF !noAdaptive

        ! REQUIRED LATER FOR THE POSITION START
        IF(Symmetry2DAxisymmetric) THEN
          !-- compute parallelogram of triangle (only simple 2 value adds/subs, other from init)
          Node1 = 2     ! normal = cross product of 1-2 and 1-3 for first triangle
          Node2 = 4     !          and 1-3 and 1-4 for second triangle
          Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,ElemID)) - xNod
          Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,ElemID)) - yNod
          Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,ElemID)) - zNod
          Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,ElemID)) - xNod
          Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,ElemID)) - yNod
          Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,ElemID)) - zNod
          IF (ABS(Vector1(3)).GT.ABS(Vector2(3))) THEN
            Vector2D(1:2) = Vector2(1:2)
          ELSE
            Vector2D(1:2) = Vector1(1:2)
          END IF
          minVec = MINLOC((/ynod, ynod+Vector2D(2)/),1)
          SELECT CASE(minVec)
          CASE(1)
            minPos(1:2) = (/xnod, ynod/)
            RVec(1:2) =  Vector2D(1:2)
          CASE(2)
            minPos(1:2) = (/xnod,ynod/) + Vector2D(1:2)
            RVec(1:2) = - Vector2D(1:2)
          END SELECT
        END IF
        ! REQUIRED LATER FOR THE POSITION END

!----- 1.: set positions
        !-- compute number of to be inserted particles
        IF (noAdaptive.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
          IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%Adaptive) THEN
            IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
              PartInsSubSide=PartInsSubSides(iSample,jSample,iSide)
            ELSE IF(DoPoissonRounding .AND. .NOT.DoTimeDepInflow)THEN
              PartIns = Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
                      * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR
              IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
                CALL abort(&
  __STAMP__&
  ,'ERROR in ParticleSurfaceflux: flux is too large for poisson sampling!')
              ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
                CALL SamplePoissonDistri( PartIns , PartInsSubSide )
              END IF
            ELSE !DoTimeDepInflow
              CALL RANDOM_NUMBER(RandVal1)
              PartInsSubSide = INT(Species(iSpec)%Surfaceflux(iSF)%PartDensity / Species(iSpec)%MacroParticleFactor &
                             * dt*RKdtFrac * Species(iSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%nVFR+RandVal1)
            END IF !DoPoissonRounding
          ELSE !Species(iSpec)%Surfaceflux(iSF)%Adaptive
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
            CASE(4)
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
          END IF ! Adaptive SurfaceFlux
        ELSE IF(.NOT.noAdaptive) THEN !Adaptive
          SELECT CASE(PartBound%AdaptiveType(currentBC))
          CASE(1) ! Pressure inlet (pressure, temperature const)
            ElemPartDensity = Species(iSpec)%Surfaceflux(iSF)%PartDensity
            T =  Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC
          CASE(2) ! adaptive Outlet/freestream
            ElemPartDensity = Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpec)
            pressure = PartBound%AdaptivePressure(currentBC)
            T = pressure / (BoltzmannConst * SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,:)))
          CASE(3) ! pressure outlet (pressure defined)
          CASE DEFAULT
            CALL abort(&
__STAMP__&
,'wrong adaptive type for Surfaceflux!')
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

          CALL RANDOM_NUMBER(RandVal1)
          PartInsSubSide = INT(ElemPartDensity / Species(iSpec)%MacroParticleFactor * dt*RKdtFrac * nVFR+RandVal1)
        ELSE IF(RadialWeighting%DoRadialWeighting) THEN
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
        ALLOCATE( particle_positions(1:PartInsSubSide*3), STAT=allocStat )
        IF (allocStat .NE. 0) THEN
          CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: cannot allocate particle_positions!')
        END IF
        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
          ALLOCATE( particle_xis(1:PartInsSubSide*2), STAT=allocStat )
          IF (allocStat .NE. 0) THEN
            CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: cannot allocate particle_xis!')
          END IF
        END IF !VeloIsNormal
        !-- put particles in subside (rejections are used if contraint reduces actual inserted number)
        iPart=1
        nReject=0
        allowedRejections=0

        IF(Symmetry2DAxisymmetric) THEN
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
        ELSE
          DO WHILE (iPart+allowedRejections .LE. PartInsSubSide)
            IF (TriaSurfaceFlux) THEN
              CALL RANDOM_NUMBER(RandVal2)
              IF (.NOT.TriaTracking) THEN !prevent inconsistency with non-triatracking by bilinear-routine (tol. might be increased)
                RandVal2 = RandVal2 + eps_nontria*(1 - 2.*RandVal2) !shift randVal off from 0 and 1
                DO WHILE (ABS(RandVal2(1)+RandVal2(2)-1.0).LT.eps_nontria) !sum must not be 1, since this corresponds to third egde
                  CALL RANDOM_NUMBER(RandVal2)
                  RandVal2 = RandVal2 + eps_nontria*(1 - 2.*RandVal2)
                END DO
              END IF
              Particle_pos = (/xNod,yNod,zNod/) + Vector1 * RandVal2(1)
              Particle_pos =       Particle_pos + Vector2 * RandVal2(2)
              PartDistance = ndist(1)*(Particle_pos(1)-midpoint(1)) & !Distance from v1-v2
                          + ndist(2)*(Particle_pos(2)-midpoint(2)) &
                          + ndist(3)*(Particle_pos(3)-midpoint(3))
              IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
                Particle_pos(1:3) = 2*midpoint(1:3)-Particle_pos(1:3)
              END IF
            ELSE !.NOT.TriaSurfaceFlux
              iLoop=0
              DO !ARM for xi considering the dA of the Subside in RefSpace
                iLoop = iLoop+1
                CALL RANDOM_NUMBER(RandVal2)
                xiab(1,1:2)=(/BezierSampleXi(ISample-1),BezierSampleXi(ISample)/) !correct order?!?
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
              CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=Particle_pos)
            END IF !TriaSurfaceFlux

            IF (Species(iSpec)%Surfaceflux(iSF)%CircularInflow) THEN !check rmax-rejection
              SELECT CASE(Species(iSpec)%Surfaceflux(iSF)%SurfFluxSideRejectType(iSide))
              CASE(0) !- RejectType=0 : complete side is inside valid bounds
                AcceptPos=.TRUE.
              CASE(1) !- RejectType=1 : complete side is outside of valid bounds
                CALL abort(&
  __STAMP__&
  ,'side outside of valid bounds was considered although nVFR=0...?!')
                !AcceptPos=.FALSE.
              CASE(2) !- RejectType=2 : side is partly inside valid bounds
                point(1)=Particle_pos(dir(2))-origin(1)
                point(2)=Particle_pos(dir(3))-origin(2)
                radius=SQRT( (point(1))**2+(point(2))**2 )
                IF ((radius.LE.Species(iSpec)%Surfaceflux(iSF)%rmax).AND.(radius.GE.Species(iSpec)%Surfaceflux(iSF)%rmin)) THEN
                  AcceptPos=.TRUE.
                ELSE
                  AcceptPos=.FALSE.
                END IF
              CASE DEFAULT
                CALL abort(&
  __STAMP__&
  ,'wrong SurfFluxSideRejectType!')
              END SELECT !SurfFluxSideRejectType
            ELSE !no check for rmax-rejection
              AcceptPos=.TRUE.
            END IF ! CircularInflow
            IF (UseMacroBody) THEN
              IF (INSIDEMACROBODY(Particle_pos)) THEN
                AcceptPos=.FALSE.
              END IF
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
        
        ParticleIndexNbr = 1
        DO iPart=1,PartInsSubSide
          IF ((iPart.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
            ParticleIndexNbr = PDM%nextFreePosition(iPartTotal + 1 &
              + PDM%CurrentNextFreePosition)
          END IF
          IF (ParticleIndexNbr .ne. 0) THEN
            PartState(1:3,ParticleIndexNbr) = particle_positions(3*(iPart-1)+1:3*(iPart-1)+3)
#ifdef CODE_ANALYZE
            IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
              IF(ParticleIndexNbr.EQ.PARTOUT)THEN
                WRITE(UNIT_stdout,'(A32)') ' ---------------------------------------------------------------'
                IPWRITE(UNIT_stdOut,'(I0,A,3(X,E15.8))') ' SurfFlux Pos:      ', PartState(1:3,ParticleIndexNbr)
                IF (TriaSurfaceFlux) THEN
                  !- the following lines are inverted to recalc. the RandVal
                  !Particle_pos = (/xNod,yNod,zNod/) + Vector1 * RandVal2(1) + Vector2 * RandVal2(2)
                  !PartDistance = ndist(1)*(Particle_pos(1)-midpoint(1)) & !Distance from v1-v2
                  !             + ndist(2)*(Particle_pos(2)-midpoint(2)) &
                  !             + ndist(3)*(Particle_pos(3)-midpoint(3))
                  !IF (PartDistance.GT.0.) THEN !flip into right triangle if outside
                  !  Particle_pos(1:3) = 2*midpoint(1:3)-Particle_pos(1:3)
                  !END IF
                  !- recalc. the RandVal assuming no flip:
                  tmpVec=PartState(1:3,ParticleIndexNbr)
                  tmpVec = tmpVec - (/xNod,yNod,zNod/) != Vector1 * RandVal2(1) + Vector2 * RandVal2(2)
                  IPWRITE(UNIT_stdOut,'(I0,A,2(X,E15.8))') ' SurfFlux RandVals1:', CalcVectorAdditionCoeffs(tmpVec,Vector1,Vector2)
                  !- recalc. the RandVal assuming flip:
                  tmpVec=2*midpoint(1:3)-PartState(1:3,ParticleIndexNbr)
                  tmpVec = tmpVec - (/xNod,yNod,zNod/) != Vector1 * RandVal2(1) + Vector2 * RandVal2(2)
                  IPWRITE(UNIT_stdOut,'(I0,A,2(X,E15.8))') ' SurfFlux RandVals2:', CalcVectorAdditionCoeffs(tmpVec,Vector1,Vector2)
                END IF
                WRITE(UNIT_stdout,'(A32)') ' ---------------------------------------------------------------'
              END IF
            END IF
#endif /*CODE_ANALYZE*/
            IF (noAdaptive) THEN
              ! check if surfaceflux is used for surface sampling (neccessary for desorption and evaporation)
              ! create linked list of surfaceflux-particle-info for sampling case
              IF ((DSMC%CalcSurfaceVal.AND.(Time.GE.(1.-DSMC%TimeFracSamp)*TEnd)) &
                  .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) THEN
                IF (PartBound%TargetBoundCond(CurrentBC).EQ.PartBound%ReflectiveBC) THEN
                  IF ( PartBound%Reactive(CurrentBC) )  THEN
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
                  END IF ! reflective bc
                END IF ! sampling is on (CalcSurfaceVal)
              END IF ! wallmodel or liquidsim
              IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) THEN
                PartState(4:5,ParticleIndexNbr) = particle_xis(2*(iPart-1)+1:2*(iPart-1)+2) !use velo as dummy-storage for xi!
              ELSE IF (Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit) THEN !PartState is used as drift for case of MB-distri!
                point(1)=PartState(dir(2),ParticleIndexNbr)-origin(1)
                point(2)=PartState(dir(3),ParticleIndexNbr)-origin(2)
                radius=SQRT( (point(1))**2+(point(2))**2 )
                phi=ATAN2(point(2),point(1))
                !-- evaluate radial fit
                vTot = Species(iSpec)%Surfaceflux(iSF)%VeloIC
                veloR=-radius*(preFac*exp(powerFac*radius)+shiftFac)
                IF (ABS(veloR).GT.1.) THEN
                  IPWRITE(*,*) 'radius=',radius
                  IPWRITE(*,*) 'veloR-ratio=',veloR
                  CALL abort(__STAMP__,&
                    'ERROR in VeloFit!')
                END IF
                PartState(3+dir(1),ParticleIndexNbr) = SIGN(vTot * SQRT(1.-veloR**2) &
                  ,Species(iSpec)%Surfaceflux(iSF)%VeloVecIC(dir(1)))
                veloR = veloR * vToT
                PartState(3+dir(2),ParticleIndexNbr) = veloR*cos(phi)
                PartState(3+dir(3),ParticleIndexNbr) = veloR*sin(phi)
              END IF !VeloIsNormal or SimpleRadialVeloFit
            END IF

            ! shift lastpartpos minimal into cell for fail-safe tracking
            LastPartPos(1:3,ParticleIndexNbr)=PartState(1:3,ParticleIndexNbr)
            !SELECT CASE(SideType(SideID))
            !CASE(PLANAR_RECT,PLANAR_NONRECT)
            !  LastPartPos(1:3,ParticleIndexNbr)=ElemBaryNGeo(1:3,ElemID) &
            !  + (PartState(1:3,ParticleIndexNbr)-ElemBaryNGeo(1:3,ElemID)) * (0.9999)
            !CASE(BILINEAR,CURVED,PLANAR_CURVED) !to be changed into more efficient method using known xi
            !  CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),Particle_pos(1:3),ElemID) !RefMap PartState
            !  DO iLoop=1,3 !shift border-RefCoords into elem
            !    IF( ABS(Particle_pos(iLoop)) .GT. 0.9999 ) THEN
            !      Particle_pos(iLoop)=SIGN(0.999999,Particle_pos(iLoop))
            !    END IF
            !  END DO
            !  CALL TensorProductInterpolation(Particle_pos(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,ElemID) &
            !    ,LastPartPos(1:3,ParticleIndexNbr)) !Map back into phys. space
            !CASE DEFAULT
            !  CALL abort(&
!__STAMP__&
!,'unknown SideType!')
            !END SELECT

!#ifdef CODE_ANALYZE
!          CALL GetPositionInRefElem(LastPartPos(1:3,ParticleIndexNbr),Particle_pos(1:3),ElemID)
!          IF (ANY(ABS(Particle_pos).GT.1.0)) THEN !maybe 1+epsInCell would be enough...
!            IPWRITE(*,*) 'Particle_pos: ',Particle_pos
!            CALL abort(&
!__STAMP__&
!,'CODE_ANALYZE: RefPos of LastPartPos is outside for ElemID. BC-cells are too deformed for surfaceflux!')
!          END IF
!#endif /*CODE_ANALYZE*/
#if defined(IMPA) || defined(ROS)
            IF(DoRefMapping)THEN
              CALL GetPositionInRefElem(PartState(1:3,ParticleIndexNbr),PartPosRef(1:3,ParticleIndexNbr),ElemID) !RefMap PartState
            END IF
            ! important for implicit, correct norm, etc.
            PartState(1:3,ParticleIndexNbr)=LastPartPos(1:3,ParticleIndexNbr)
#endif /*IMPA*/
#ifdef CODE_ANALYZE
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
#endif /*CODE_ANALYZE*/
            PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
            PDM%dtFracPush(ParticleIndexNbr) = .TRUE.
            PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
            PEM%Element(ParticleIndexNbr) = ElemID
            PEM%lastElement(ParticleIndexNbr) = ElemID !needed when ParticlePush is not executed, e.g. "delay"
            iPartTotal = iPartTotal + 1
            IF (VarTimeStep%UseVariableTimeStep) THEN
              VarTimeStep%ParticleTimeStep(ParticleIndexNbr) &
                = CalcVarTimeStep(PartState(1,ParticleIndexNbr),PartState(2,ParticleIndexNbr),PEM%Element(ParticleIndexNbr))
            END IF
            IF (RadialWeighting%DoRadialWeighting) THEN
              PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(2,ParticleIndexNbr), 1,ParticleIndexNbr)
            END IF
            IF(CalcMassflowRate) THEN
              Species(iSpec)%Surfaceflux(iSF)%SampledMassflow = Species(iSpec)%Surfaceflux(iSF)%SampledMassflow &
                                                                + GetParticleWeight(ParticleIndexNbr)
            END IF
          ELSE
            CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
          END IF
        END DO
        DEALLOCATE(particle_positions)
        IF (Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. .NOT.TriaSurfaceFlux) DEALLOCATE(particle_xis)
!----- 2a.: set velocities if special for each subside
        IF (TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).NE.'constant' &
          .OR. Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
          CALL SetSurfacefluxVelocities(iSpec,iSF,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartInsSubSide)
        END IF

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
    IF (TRIM(Species(iSpec)%Surfaceflux(iSF)%velocityDistribution).EQ.'constant' &
      .AND. .NOT.Species(iSpec)%Surfaceflux(iSF)%SimpleRadialVeloFit &
      .AND. .NOT.Species(iSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      CALL SetParticleVelocity(iSpec,iSF,NbrOfParticle,2)
    END IF
    CALL SetParticleChargeAndMass(iSpec,NbrOfParticle)
    IF (usevMPF.AND.(.NOT.RadialWeighting%DoRadialWeighting)) CALL SetParticleMPF(iSpec,NbrOfParticle)
    ! define molecule stuff
    IF (useDSMC.AND.(CollisMode.GT.1)) THEN
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

    ! instead of an UpdateNextfreePosition we update the particleVecLength only - enough ?!?
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
    PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle
#if USE_LOADBALANCE
    CALL LBPauseTime(LB_SURFFLUX,tLBStart)
#endif /*USE_LOADBALANCE*/
    IF (noAdaptive) THEN
      ! Sample Energies on Surfaces when particles are emitted from them
      IF (NbrOfParticle.NE.PartsEmitted) THEN
        ! should be equal for including the following lines in tSurfaceFlux
        CALL abort(&
__STAMP__&
,'ERROR in ParticleSurfaceflux: NbrOfParticle.NE.PartsEmitted')
      END IF
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
                VeloReal = VECNORM(PartState(4:6,iPart))
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
                CALL LBElemSplitTime(PEM%Element(PartID),tLBStart)
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
    END IF
  END DO !iSF
END DO !iSpec

END SUBROUTINE ParticleSurfaceflux


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
USE MOD_Particle_Mesh          ,ONLY: GetGlobalNonUniqueSideID
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
USE MOD_Particle_Boundary_Vars, ONLY : PartBound
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

IF(iSF.GT.Species(FractNbr)%nSurfacefluxBCs)THEN
  SELECT CASE(PartBound%AdaptiveType(currentBC))
  CASE(1) ! Pressure inlet (pressure, temperature const)
    T =  Species(FractNbr)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(2) ! adaptive Outlet/freestream
    pressure = PartBound%AdaptivePressure(Species(FractNbr)%Surfaceflux(iSF)%BC)
    T = pressure / (BoltzmannConst * SUM(Adaptive_MacroVal(DSMC_NUMDENS,ElemID,:)))
    !T = SQRT(Adaptive_MacroVal(4,ElemID,FractNbr)**2+Adaptive_MacroVal(5,ElemID,FractNbr)**2 &
    !  + Adaptive_MacroVal(6,ElemID,FractNbr)**2)
  CASE(3) ! pressure outlet (pressure defined)
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'wrong adaptive type for Surfaceflux velocities!')
  END SELECT
  VeloVec(1) = Adaptive_MacroVal(DSMC_VELOX,ElemID,FractNbr)
  VeloVec(2) = Adaptive_MacroVal(DSMC_VELOY,ElemID,FractNbr)
  VeloVec(3) = Adaptive_MacroVal(DSMC_VELOZ,ElemID,FractNbr)
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
ELSE
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
END IF

!-- set velocities
SELECT CASE(TRIM(velocityDistribution))
CASE('constant') !constant with normal velocities (for VeloVecIC see SetParticleVelocity!)

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
        CALL abort(&
__STAMP__&
,'this should not happen!')
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
  IF (.NOT.Species(FractNbr)%Surfaceflux(iSF)%SimpleRadialVeloFit) THEN
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
  END IF !.NOT.SimpleRadialVeloFit

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
!-- 0b.: initialize DataTriaSF if particle-dependent (as in case of SimpleRadialVeloFit), drift vector is already in PartState!!!
      ELSE IF (Species(FractNbr)%Surfaceflux(iSF)%SimpleRadialVeloFit) THEN
        VeloIC = VECNORM(PartState(4:6,PositionNbr))
        IF (ALMOSTZERO(VeloIC)) THEN
          projFak = 1. !dummy
          a = 0.
        ELSE
          projFak = DOT_PRODUCT(vec_nIn,PartState(4:6,PositionNbr)) / VeloIC
          a = VeloIC * projFak / SQRT(2.*BoltzmannConst*T/Species(FractNbr)%MassIC) !speed ratio proj. to inwards n (can be negative!)
        END IF
        Velo_t1 = DOT_PRODUCT(vec_t1,PartState(4:6,PositionNbr)) !v in t1-dir
        Velo_t2 = DOT_PRODUCT(vec_t2,PartState(4:6,PositionNbr)) !v in t2-dir
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
      END IF !VeloIsNormal, else if SimpleRadialVeloFit
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
