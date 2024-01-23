#include "piclas.h"
!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
#if USE_HDG
LOGICAL                         :: DoDirichletDeposition     ! flag to switch deposition on Dirichlet sides on/off
#endif /*USE_HDG*/
LOGICAL                         :: RelaxDeposition           ! relaxation of current PartSource with RelaxFac into PartSourceOld
LOGICAL                         :: DoHaloDepo                ! Flag for enlarging the deposition region (implicit and dielectric)
REAL                            :: RelaxFac

REAL,ALLOCATABLE                 :: PartSource(:,:,:,:,:)    ! PartSource(1:4,PP_N,PP_N,PP_N,nElems) containing
!                                                            ! current and charge density source terms for Maxwell/Poisson systems
REAL,ALLOCATABLE                :: PartSourceGlob(:,:,:,:,:)
!REAL,ALLOCATABLE                :: PartSourceProc(:,:,:,:,:)
!REAL,ALLOCATABLE                :: PartSourceLoc(:,:,:,:,:)
REAL,ALLOCATABLE                :: PartSourceTmp (:,:,:,:)
INTEGER, ALLOCATABLE            :: nDepoDOFPerProc(:)
INTEGER, ALLOCATABLE            :: nDepoOffsetProc(:)

LOGICAL                         :: PartSourceConstExists
REAL,ALLOCATABLE                :: PartSourceConst(:,:,:,:,:)! PartSource(1:4,PP_N,PP_N,PP_N,nElems) const. part of Source
REAL,ALLOCATABLE                :: PartSourceOld(:,:,:,:,:,:)! PartSource(:,2,PP_N,PP_N,PP_N,nElems) prev. and sec. prev. Source
REAL,ALLOCATABLE                :: GaussBorder(:)            ! 1D coords of gauss points in -1|1 space
INTEGER,ALLOCATABLE             :: GaussBGMIndex(:,:,:,:,:)  ! Background mesh index of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                :: GaussBGMFactor(:,:,:,:,:) ! BGM factor of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
REAL,ALLOCATABLE                :: GPWeight(:,:,:,:,:,:,:)   ! Weights for splines deposition (check pic_depo for details)
CHARACTER(LEN=256)              :: DepositionType            ! Type of Deposition-Method
CHARACTER(LEN=50)               :: VerifyChargeStr           ! Auxiliary string for output to std.out for VerifyCharge=T
INTEGER,ALLOCATABLE             :: PartToFIBGM(:,:)          ! Mapping form Particle to FIBGM
REAL,ALLOCATABLE                :: ElemRadius2_sf(:)         ! elem radius plus radius_sf
REAL, ALLOCATABLE               :: BGMSource(:,:,:,:)
REAL                            :: SFAdaptiveDOF             ! Average number of DOF in shape function radius (assuming a Cartesian
!                                                            ! grid with equal elements). Only implemented for PIC-Deposition-Type =
!                                                            ! shape_function_adaptive (2). The maximum number of DOF is limited by
!                                                            ! the polynomial degree and is (4/3)*Pi*(N+1)^3. Default is 33.
LOGICAL                         :: SFAdaptiveSmoothing       ! Enable smooth transition of element-dependent radius when
                                                             ! using shape_function_adaptive, default=FALSE
REAL                            :: r_sf                      ! cutoff radius of shape function
REAL                            :: r2_sf                     ! cutoff radius of shape function * cutoff radius of shape function
REAL                            :: r2_sf_inv                 ! 1/cutoff radius of shape function * cutoff radius of shape function
REAL                            :: w_sf                      ! shapefuntion weight
REAL                            :: r_sf0                     ! minimal shape function radius
REAL                            :: r_sf_scale                ! scaling of shape function radius
REAL                            :: BetaFac                   ! betafactor of shape-function || integral =1
INTEGER                         :: sf1d_dir                  ! direction of 1D shape function
LOGICAL                         :: sfDepo3D                  ! when using 1D or 2D deposition, the charge can be deposited over the
!                                                            ! a line (1D) or area (2D)
REAL                            :: dimFactorSF               ! Scaling factor when using sfDepo3D=F
LOGICAL,ALLOCATABLE             :: ChargeSFDone(:)           ! Element flag for cycling already completed elements
LOGICAL                         :: DoSFChargeCons
!                                                            ! volume (3D) or line (1D) / area (2D)
INTEGER                         :: NDepo                     ! polynomial degree of delta distribution
REAL,ALLOCATABLE                :: tempcharge(:)             ! temp-charge for epo. Kernel
REAL,ALLOCATABLE                :: NDepoChooseK(:,:)         ! array n over n
REAL,ALLOCATABLE                :: wBaryNDepo(:)             ! barycentric weights for deposition
REAL,ALLOCATABLE                :: swGPNDepo(:)              ! integration weights for deposition
REAL,ALLOCATABLE                :: XiNDepo(:)                ! gauss position of bary centers
REAL,ALLOCATABLE                :: Vdm_NDepo_GaussN(:,:)     ! Vandermonde between different polynomial degrees
REAL,ALLOCATABLE                :: DDMassinv(:,:,:,:)        ! inverse mass-matrix for deposition
REAL,ALLOCATABLE                :: Vdm_EquiN_GaussN(:,:)     ! Vandermonde from equidistant points to Gauss Points
INTEGER                         :: alpha_sf                  ! shape function exponent
INTEGER                         :: dim_sf                    ! 1D, 2D or 3D shape function
INTEGER                         :: dim_sf_dir                ! Get shape function direction for 1D (the direction in which the charge
!                                                            ! will be distributed) and 2D (the direction in which the charge will be
!                                                            ! constant)
INTEGER                         :: dim_sf_dir1               ! 1st perpendicular direction used in 2D shape function
INTEGER                         :: dim_sf_dir2               ! 2nd perpendicular direction used in 2D shape function
INTEGER                         :: dim_periodic_vec1         ! 1st periodic vector used in 2D shape function
INTEGER                         :: dim_periodic_vec2         ! 2nd periodic vector used in 2D shape function (if available)
REAL                            :: BGMdeltas(3)              ! Background mesh size in x,y,z
REAL                            :: FactorBGM(3)              ! Divider for BGM (to allow real numbers)
REAL                            :: BGMVolume                 ! Volume of a BGM Cell
INTEGER                         :: BGMminX                   ! Local minimum BGM Index in x
INTEGER                         :: BGMminY                   ! Local minimum BGM Index in y
INTEGER                         :: BGMminZ                   ! Local minimum BGM Index in z
INTEGER                         :: BGMmaxX                   ! Local maximum BGM Index in x
INTEGER                         :: BGMmaxY                   ! Local maximum BGM Index in y
INTEGER                         :: BGMmaxZ                   ! Local maximum BGM Index in z
LOGICAL                         :: Periodic_Depo             ! Flag for periodic treatment for deposition
REAL                            :: totalChargePeriodicSF     ! total charge of particle that is mirrored over periodic boundaries
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

REAL,ALLOCPOINT                 :: SFElemr2_Shared(:,:) ! index 1: radius, index 2: radius squared

REAL,ALLOCATABLE                :: NodeSource(:,:)
INTEGER,ALLOCATABLE             :: NodeSendDepoRankToGlobalRank(:)
INTEGER,ALLOCATABLE             :: NodeRecvDepoRankToGlobalRank(:)
INTEGER,ALLOCATABLE             :: DepoNodetoGlobalNode(:)
INTEGER                         :: nDepoNodes
INTEGER                         :: nDepoNodesTotal
INTEGER                         :: nNodeSendExchangeProcs
INTEGER                         :: nNodeRecvExchangeProcs
! Additional source for cell_volweight_mean (external or surface charge) that accumulates over time in elements adjacent to
! dielectric interfaces.
REAL,ALLOCATABLE                :: NodeSourceExt(:) ! It contains the global, synchronized surface charge contribution that is
!                                                   ! read and written to .h5

INTEGER            :: nTotalPeriodicNodes
INTEGER,ALLOCPOINT :: Periodic_nNodes(:)
INTEGER,ALLOCPOINT :: Periodic_Nodes(:)
INTEGER,ALLOCPOINT :: Periodic_offsetNode(:)

#if USE_MPI
INTEGER            :: Periodic_nNodes_Shared_Win
INTEGER,ALLOCPOINT :: Periodic_nNodes_Shared(:)
INTEGER            :: Periodic_Nodes_Shared_Win
INTEGER,ALLOCPOINT :: Periodic_Nodes_Shared(:)
INTEGER            :: Periodic_offsetNode_Shared_Win
INTEGER,ALLOCPOINT :: Periodic_offsetNode_Shared(:)

REAL,ALLOCATABLE                :: NodeSourceExtTmp(:) ! It contains the local non-synchronized surface charge contribution (does
!                                                      ! not consider the charge contribution from restart files). This
!                                                      ! contribution accumulates over time, but remains local to each processor
!                                                      ! as it is communicated via the container NodeSourceExt.

INTEGER                         :: SFElemr2_Shared_Win

! Send direction of nodes (can be different from number of receive nodes for each processor)
TYPE tNodeMappingSend
  INTEGER,ALLOCATABLE           :: SendNodeUniqueGlobalID(:)
  REAL,ALLOCATABLE              :: SendNodeSourceCharge(:)
  REAL,ALLOCATABLE              :: SendNodeSourceCurrent(:,:)
  REAL,ALLOCATABLE              :: SendNodeSourceExt(:)
  INTEGER                       :: nSendUniqueNodes
END TYPE
TYPE (tNodeMappingSend),ALLOCATABLE      :: NodeMappingSend(:)

! Receive direction of nodes (can be different from number of send nodes for each processor)
TYPE tNodeMappingRecv
  INTEGER,ALLOCATABLE           :: RecvNodeUniqueGlobalID(:)
  REAL,ALLOCATABLE              :: RecvNodeSourceCharge(:)
  REAL,ALLOCATABLE              :: RecvNodeSourceCurrent(:,:)
  REAL,ALLOCATABLE              :: RecvNodeSourceExt(:)
  INTEGER                       :: nRecvUniqueNodes
END TYPE
TYPE (tNodeMappingRecv),ALLOCATABLE      :: NodeMappingRecv(:)

TYPE tShapeMapping
  INTEGER,ALLOCATABLE           :: RecvShapeElemID(:)
  INTEGER                       :: nRecvShapeElems
  INTEGER,ALLOCATABLE           :: SendShapeElemID(:)
  INTEGER                       :: nSendShapeElems
  REAL,ALLOCATABLE              :: RecvBuffer(:,:,:,:,:)
  REAL,ALLOCATABLE              :: SendBuffer(:,:,:,:,:)
  LOGICAL,ALLOCATABLE           :: DoSendElem(:)
  INTEGER                       :: nNonZeroSendElems
  INTEGER                       :: Rank
END TYPE
TYPE(tShapeMapping),ALLOCATABLE :: ShapeMapping(:)

TYPE tCNShapeMapping
  INTEGER,ALLOCATABLE           :: RecvShapeElemID(:)
  INTEGER,ALLOCATABLE           :: RecvShapeProcElemID(:,:)
  INTEGER,ALLOCATABLE           :: SendShapeElemID(:)
  INTEGER,ALLOCATABLE           :: SendShapeProcElemID(:,:)
  INTEGER                       :: nRecvShapeElems(2)
  INTEGER                       :: nSendShapeElems(2)
  REAL,ALLOCATABLE              :: RecvBuffer(:,:,:,:,:)
  REAL,ALLOCATABLE              :: SendBuffer(:,:,:,:,:)
END TYPE
TYPE(tCNShapeMapping),ALLOCATABLE ::CNShapeMapping(:)

INTEGER                         :: ShapeElemProcSend_Shared_Win
LOGICAL,ALLOCPOINT              :: ShapeElemProcSend_Shared(:,:)
LOGICAL,ALLOCATABLE             :: FlagShapeElem(:)
INTEGER,ALLOCATABLE             :: SendElemShapeID(:)         ! mapping from ShapeElemID to CNElemID
INTEGER                         :: nShapeExchangeProcs

!INTEGER             :: SendRequest
INTEGER,ALLOCATABLE :: RecvRequest(:), SendRequest(:), CNRankToSendRank(:)
#endif

!===================================================================================================================================
END MODULE MOD_PICDepo_Vars
