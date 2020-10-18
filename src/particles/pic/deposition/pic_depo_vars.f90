#include "piclas.h"
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
MODULE MOD_PICDepo_Vars
!===================================================================================================================================
! Contains the variables for the particle deposition
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                         :: DoDeposition              ! flag to switch deposition on/off
LOGICAL                         :: RelaxDeposition           ! relaxation of current PartSource with RelaxFac into PartSourceOld
REAL                            :: RelaxFac

REAL,ALLOCPOINT                 :: PartSource(:,:,:,:,:)     ! PartSource(1:4,PP_N,PP_N,PP_N,nComputeNodeTotalElems) containing 
!                                                            ! current and charge density source terms for Maxwell/Poisson systems
!                                                            ! Access array with CNElemID = GetCNElemID(GlobalElemID)
!                                                            !                            = GetCNElemID(iElem+offSetElem)
#if USE_MPI
REAL, ALLOCATABLE               :: PartSourceProc(:,:,:,:,:)
INTEGER                         :: PartSource_Shared_Win
REAL,ALLOCPOINT                 :: PartSource_Shared(:)
#endif

LOGICAL                         :: PartSourceConstExists
REAL,ALLOCATABLE                :: PartSourceConst(:,:,:,:,:)! PartSource(1:4,PP_N,PP_N,PP_N,nElems) const. part of Source
REAL,ALLOCATABLE                :: PartSourceOld(:,:,:,:,:,:)! PartSource(:,2,PP_N,PP_N,PP_N,nElems) prev. and sec. prev. Source
REAL,ALLOCATABLE                :: GaussBorder(:)            ! 1D coords of gauss points in -1|1 space
INTEGER,ALLOCATABLE             :: GaussBGMIndex(:,:,:,:,:)  ! Background mesh index of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                :: GaussBGMFactor(:,:,:,:,:) ! BGM factor of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                :: GPWeight(:,:,:,:,:,:,:)   ! Weights for splines deposition (check pic_depo for details)
CHARACTER(LEN=256)              :: DepositionType            ! Type of Deposition-Method
INTEGER,ALLOCATABLE             :: PartToFIBGM(:,:)          ! Mapping form Particle to FIBGM
REAL,ALLOCATABLE                :: ElemRadius2_sf(:)         ! elem radius plus radius_sf
REAL, ALLOCATABLE               :: BGMSource(:,:,:,:)
REAL, ALLOCATABLE               :: NodeSourceExt(:)          ! Additional source for cell_volweight_mean (external or surface charge)
!                                                            ! accumulates over time
REAL, ALLOCATABLE               :: NodeSourceExtTmp(:)       ! Additional source for cell_volweight_mean (external or surface charge)
!                                                            ! temporary variable, which is nullified repeatedly
REAL                            :: r_sf                      ! cutoff radius of shape function
REAL                            :: r2_sf                     ! cutoff radius of shape function * cutoff radius of shape function
REAL                            :: r2_sf_inv                 ! 1/cutoff radius of shape function * cutoff radius of shape function
REAL                            :: w_sf                      ! shapefuntion weight
REAL                            :: r_sf0                     ! minimal shape function radius
REAL                            :: r_sf_scale                ! scaling of shape function radius
REAL                            :: BetaFac                   ! betafactor of shape-function || integral =1
INTEGER                         :: sf1d_dir                  ! direction of 1D shape function
LOGICAL                         :: sfDepo3D                  ! when using 1D or 2D deposition, the charge can be deposited over the
LOGICAL                         :: DoSFChargeCons
!                                                            ! volume (3D) or line (1D) / area (2D)
INTEGER                         :: NDepo                     ! polynomial degree of delta distri
REAL,ALLOCATABLE                :: tempcharge(:)             ! temp-charge for epo. kernal
REAL,ALLOCATABLE                :: NDepoChooseK(:,:)         ! array n over n
REAL,ALLOCATABLE                :: wBaryNDepo(:)             ! barycentric weights for deposition
REAL,ALLOCATABLE                :: swGPNDepo(:)              ! integration weights for deposition
REAL,ALLOCATABLE                :: XiNDepo(:)                ! gauss position of barycenters
REAL,ALLOCATABLE                :: Vdm_NDepo_GaussN(:,:)     ! VdM between different polynomial degrees
REAL,ALLOCATABLE                :: DDMassinv(:,:,:,:)        ! inverse mass-matrix for deposition
!LOGICAL                         :: DeltaDistriChangeBasis    ! Change polynomial degree
LOGICAL                         :: DoSFEqui                  ! use equidistant points for SF
LOGICAL                         :: DoSFLocalDepoAtBounds     ! Do not use shape function deposition in elements where a boundary
!                                                            ! would truncate the shape function. Use a local deposition in these
!                                                            ! elements instead of the shape function
INTEGER                         :: SfRadiusInt               ! radius integer for cylindrical and spherical shape function
REAL,ALLOCATABLE                :: ElemDepo_xGP(:,:,:,:,:)   ! element xGPs for deposition
REAL,ALLOCATABLE                :: Vdm_EquiN_GaussN(:,:)     ! Vdm from equidistant points to Gauss Points
INTEGER                         :: alpha_sf                  ! shapefuntion exponent
REAL                            :: BGMdeltas(3)              ! Backgroundmesh size in x,y,z
REAL                            :: FactorBGM(3)              ! Divider for BGM (to allow real numbers)
REAL                            :: BGMVolume                 ! Volume of a BGM Cell
INTEGER                         :: BGMminX                   ! Local minimum BGM Index in x
INTEGER                         :: BGMminY                   ! Local minimum BGM Index in y
INTEGER                         :: BGMminZ                   ! Local minimum BGM Index in z
INTEGER                         :: BGMmaxX                   ! Local maximum BGM Index in x
INTEGER                         :: BGMmaxY                   ! Local maximum BGM Index in y
INTEGER                         :: BGMmaxZ                   ! Local maximum BGM Index in z
LOGICAL                         :: Periodic_Depo             ! Flag for periodic treatment for deposition
!INTEGER                         :: DeltaType                 ! Flag
INTEGER                         :: NKnots
REAL,ALLOCATABLE                :: Knots(:)
LOGICAL                         :: OutputSource              ! write the source to hdf5
REAL,ALLOCATABLE                :: CellVolWeightFac(:)
INTEGER                         :: VolIntOrder
REAL,ALLOCATABLE                :: VolInt_X(:)
REAL,ALLOCATABLE                :: VolInt_W(:)
REAL,ALLOCATABLE                :: CellVolWeight_Volumes(:,:,:,:)
REAL,ALLOCPOINT                 :: NodeVolume(:)
#if USE_MPI
INTEGER                         :: NodeVolume_Shared_Win
REAL,ALLOCPOINT                 :: NodeVolume_Shared(:)
#endif

REAL,ALLOCPOINT                 :: NodeSource(:,:)
#if USE_MPI
REAL, ALLOCATABLE               :: NodeSourceLoc(:,:)
INTEGER                         :: NodeSource_Shared_Win
REAL,ALLOCPOINT                 :: NodeSource_Shared(:,:)

TYPE tNodeMapping
  INTEGER,ALLOCATABLE                   :: RecvNodeUniqueGlobalID(:)
  INTEGER,ALLOCATABLE                   :: SendNodeUniqueGlobalID(:)
  REAL,ALLOCATABLE                      :: RecvNodeSource(:,:)
  REAL,ALLOCATABLE                      :: SendNodeSource(:,:)
  INTEGER                               :: nSendUniqueNodes
  INTEGER                               :: nRecvUniqueNodes
END TYPE
TYPE (tNodeMapping),ALLOCATABLE      :: NodeMapping(:)
#endif


INTEGER                         :: NbrOfSFdepoLayers             ! Number of const. source layer for sf-depo at planar BCs
LOGICAL                         :: PrintSFDepoWarnings           ! flag to print the warnings
LOGICAL                         :: ConstantSFdepoLayers          ! depo just once
LOGICAL                         :: SFdepoLayersAlreadyDone       ! flag for skipping the depo (i.e., when layers are const.)
REAL    , ALLOCATABLE           :: SFdepoLayersGeo(:,:,:)        ! 1:nFixes;1:2(base,normal);1:3(x,y,z) normal outwards!!!
REAL    , ALLOCATABLE           :: SFdepoLayersBounds(:,:,:)     ! 1:nFixes;1:2(min,max);1:3(x,y,z)
LOGICAL , ALLOCATABLE           :: SFdepoLayersUseFixBounds(:)   ! use alls planes of SFdepoFixes as additional bounds?
CHARACTER(LEN=256),ALLOCATABLE  :: SFdepoLayersSpace(:)          ! name of space (cuboid or cylinder)
REAL    , ALLOCATABLE           :: SFdepoLayersBaseVector(:,:,:) ! 1:nFixes;1:2;1:3(x,y,z)
INTEGER , ALLOCATABLE           :: SFdepoLayersSpec(:)           ! species of particles for respective layer
REAL    , ALLOCATABLE           :: SFdepoLayersMPF(:)            ! MPF for layerParts
REAL    , ALLOCATABLE           :: SFdepoLayersPartNum(:)        ! number of particles in volume
REAL    , ALLOCATABLE           :: SFdepoLayersRadius(:)         ! radius for cylinder-space
LOGICAL                         :: SFResampleAnalyzeSurfCollis
TYPE tLastAnalyzeSurfCollis
  INTEGER                       :: PartNumberSamp                ! number of parts from last sampling
  INTEGER                       :: PartNumberReduced             ! max. allowed number of parts to be saved
  LOGICAL                       :: ReducePartNumber              ! reduce PartNumberSamp to PartNumberReduced
  INTEGER                       :: PartNumberDepo                ! number of parts to be inserted in depo
  REAL, ALLOCATABLE             :: WallState(:,:)                ! Pos at wall and velocities from last sampling
  INTEGER, ALLOCATABLE          :: Species(:)                    ! Spec of parts
  REAL                          :: pushTimeStep                  ! timestep for (fractional) euler push from wall into ghost domain
  INTEGER                       :: PartNumThreshold              ! Threshold for checking inserted parts per depo (otherwise abort)
  REAL                          :: NormVecOfWall(3)              ! normVec for pushTimeStep
  REAL                          :: Bounds(1:2,1:3)               ! bounds after push for preventing parts outside of...
                                                                 ! ...extruded domain 1:2(min,max);1:3(x,y,z)
  LOGICAL                       :: UseFixBounds                  ! use alls planes of SFdepoFixes as additional bounds?
  LOGICAL                       :: Restart                       ! read-in old DSMCSurfCollis-file for restart
  CHARACTER(LEN=256)            :: DSMCSurfCollisRestartFile
  INTEGER                       :: NumberOfBCs                   ! Nbr of BC to be analyzed (def.: 1)
  INTEGER, ALLOCATABLE          :: BCs(:)                        ! BCs to be analyzed (def.: 0 = all)
  INTEGER                       :: NbrOfSpeciesForDtCalc         ! Number of species used for SFResample-dt (def.: 1)
  INTEGER, ALLOCATABLE          :: SpeciesForDtCalc(:)           ! Species used for SFResample-dt (def.: 0 = all)
END TYPE
TYPE(tLastAnalyzeSurfCollis)    :: LastAnalyzeSurfCollis

#if USE_MPI
TYPE tShapeMapping
  INTEGER,ALLOCATABLE           :: RecvShapeElemID(:)
  INTEGER                       :: nRecvShapeElems
  REAL,ALLOCATABLE              :: RecvBuffer(:,:,:,:,:)
END TYPE
TYPE(tShapeMapping),ALLOCATABLE :: ShapeMapping(:)

TYPE tCNShapeMapping
  INTEGER,ALLOCATABLE           :: RecvShapeElemID(:)
  INTEGER,ALLOCATABLE           :: SendShapeElemID(:)
  INTEGER                       :: nRecvShapeElems
  INTEGER                       :: nSendShapeElems
  REAL,ALLOCATABLE              :: RecvBuffer(:,:,:,:,:)
  REAL,ALLOCATABLE              :: SendBuffer(:,:,:,:,:)
END TYPE
TYPE(tCNShapeMapping),ALLOCATABLE ::CNShapeMapping(:)


INTEGER                         :: nSendShapeElems            ! number of halo elements on proc to communicate with shape function
INTEGER,ALLOCATABLE             :: SendShapeElemID(:)         ! mapping from CNElemID to ShapeElemID
INTEGER,ALLOCATABLE             :: SendElemShapeID(:)         ! mapping from ShapeElemID to CNElemID
#endif

!===================================================================================================================================
END MODULE MOD_PICDepo_Vars
