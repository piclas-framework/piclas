!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_Photon_TrackingOutput
!===================================================================================================================================
! Module for the main radiation transport routines
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
PUBLIC :: WritePhotonSurfSampleToHDF5,WritePhotonVolSampleToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WritePhotonVolSampleToHDF5()
!===================================================================================================================================
! Writes Radiation values to HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars            ,ONLY: nElems,MeshFile,offSetElem,sJ
USE MOD_RayTracing_Vars      ,ONLY: Ray,nVarRay,U_N_Ray,U_N_Ray_loc,N_DG_Ray,PREF_VDM_Ray,N_DG_Ray_loc,N_Inter_Ray
USE MOD_RayTracing_Vars      ,ONLY: RayElemPassedEnergyLoc1st,RayElemPassedEnergyLoc2nd
USE MOD_RayTracing_Vars      ,ONLY: RaySecondaryVectorX,RaySecondaryVectorY,RaySecondaryVectorZ,ElemVolume,RayElemEmission
USE MOD_HDF5_output          ,ONLY: GatheredWriteArray
#if USE_MPI
USE MOD_RayTracing_Vars      ,ONLY: RayElemPassedEnergy_Shared,RayElemOffset,RayElemPassedEnergyHO_Shared
#else
USE MOD_RayTracing_Vars      ,ONLY: RayElemPassedEnergy
#endif /*USE_MPI*/
USE MOD_io_HDF5
USE MOD_HDF5_output          ,ONLY: GenerateFileSkeleton
USE MOD_HDF5_Output_ElemData ,ONLY: WriteAdditionalElemData
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nGlobalElems
USE MOD_ChangeBasis          ,ONLY: ChangeBasis3D
USE MOD_Interpolation_Vars   ,ONLY: NodeType
USE MOD_Interpolation        ,ONLY: GetVandermonde
USE MOD_Mesh_Tools           ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Vars   ,ONLY: ElemVolume_Shared
USE MOD_Photon_TrackingVars  ,ONLY: RadiationVolState
USE MOD_HDF5_Output_State    ,ONLY: WriteElemDataToSeparateContainer
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: iElem,iGlobalElem,iCNElem,Nloc,k,l,m
CHARACTER(LEN=255), ALLOCATABLE     :: StrVarNames(:)
REAL                                :: UNMax(nVarRay,0:Ray%NMax,0:Ray%NMax,0:Ray%NMax,PP_nElems)
#if USE_MPI
INTEGER                             :: NlocOffset
#endif /*USE_MPI*/
REAL                                :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL                                :: J_Nmax(1:1,0:Ray%NMax,0:Ray%NMax,0:Ray%NMax)
REAL                                :: J_Nloc(1:1,0:Ray%NMax,0:Ray%NMax,0:Ray%NMax)
REAL                                :: IntegrationWeight
REAL                                :: Vdm_GaussN_NMax(0:PP_N,0:Ray%NMax)    !< for interpolation to Analyze points (from NodeType nodes to Gauss-Lobatto nodes)
REAL, ALLOCATABLE                   :: Vdm_GaussN_Nloc(:,:)    !< for interpolation to Analyze points (from NodeType nodes to Gauss-Lobatto nodes)
REAL, PARAMETER                     :: tolerance=1e-2
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') ' WRITE Radiation TO HDF5 FILE...'

ALLOCATE(RayElemPassedEnergyLoc1st(1:nElems))
ALLOCATE(RayElemPassedEnergyLoc2nd(1:nElems))
ALLOCATE(ElemVolume(1:nElems))
RayElemPassedEnergyLoc1st=-1.0
RayElemPassedEnergyLoc2nd=-1.0
ElemVolume=-1.0
CALL AddToElemData(ElementOut,'RayElemPassedEnergy1st',RealArray=RayElemPassedEnergyLoc1st)
CALL AddToElemData(ElementOut,'RayElemPassedEnergy2nd',RealArray=RayElemPassedEnergyLoc2nd)
CALL AddToElemData(ElementOut,'RaySecondaryVectorX'   ,RealArray=RaySecondaryVectorX)
CALL AddToElemData(ElementOut,'RaySecondaryVectorY'   ,RealArray=RaySecondaryVectorY)
CALL AddToElemData(ElementOut,'RaySecondaryVectorZ'   ,RealArray=RaySecondaryVectorZ)
CALL AddToElemData(ElementOut,'ElemVolume'            ,RealArray=ElemVolume)
CALL AddToElemData(ElementOut,'Nloc',IntArray=N_DG_Ray_loc)

! Copy data from shared array
DO iElem = 1, nElems
  N_DG_Ray_loc(iElem) = N_DG_Ray(iElem + offsetElem)
END DO ! iElem = 1, nElems

! Sanity check
IF(ANY(N_DG_Ray_loc.LE.0)) CALL abort(__STAMP__,'N_DG_Ray_loc cannot contain zeros!')

ALLOCATE(StrVarNames(1:nVarRay))
StrVarNames(1)='RayElemPassedEnergy1st'
StrVarNames(2)='RayElemPassedEnergy2nd'
StrVarNames(3)='ElemVolume'
StrVarNames(4)='PhotonEnergyDensity1st'
StrVarNames(5)='PhotonEnergyDensity2nd'

#if USE_MPI
CALL ExchangeRayVolInfo()
#endif /*USE_MPI*/

! p-refinement: Interpolate lower degree to higher degree (other way around would require model=T)
CALL GetVandermonde(PP_N, NodeType, Ray%NMax, Ray%NodeType, Vdm_GaussN_NMax, modal=.FALSE.)

#if USE_MPI
ASSOCIATE( RayElemPassedEnergy => RayElemPassedEnergy_Shared )
#endif /*USE_MPI*/
  DO iElem=1,PP_nElems
    iGlobalElem = iElem+offSetElem
    iCNElem = GetCNElemID(iGlobalElem)

    ! 1. Elem-constant data
    ! Primary energy
    RayElemPassedEnergyLoc1st(iElem) = RayElemPassedEnergy(1,iGlobalElem) !/ ElemVolume_Shared(iCNElem)
    ! Secondary energy
    RayElemPassedEnergyLoc2nd(iElem) = RayElemPassedEnergy(2,iGlobalElem) !/ ElemVolume_Shared(iCNElem)
    ! Activate volume emission
    IF((RayElemPassedEnergyLoc1st(iElem).GT.0.0)) RayElemEmission(1,iElem) = .TRUE.
    IF((RayElemPassedEnergyLoc2nd(iElem).GT.0.0)) RayElemEmission(2,iElem) = .TRUE.
    ! Check if secondary energy is greater than zero
    IF(RayElemPassedEnergyLoc2nd(iElem).GT.0.0)THEN
      IF(RayElemPassedEnergy(6,iGlobalElem).LE.0.0) CALL abort(__STAMP__,'Secondary ray counter is zero but energy is not!')
      ! x-, y- and z-direction of secondary energy
      RaySecondaryVectorX(iElem) = RayElemPassedEnergy(3,iGlobalElem) / RayElemPassedEnergy(6,iGlobalElem)
      RaySecondaryVectorY(iElem) = RayElemPassedEnergy(4,iGlobalElem) / RayElemPassedEnergy(6,iGlobalElem)
      RaySecondaryVectorZ(iElem) = RayElemPassedEnergy(5,iGlobalElem) / RayElemPassedEnergy(6,iGlobalElem)
    ELSE
      RaySecondaryVectorX(iElem) = 0.
      RaySecondaryVectorY(iElem) = 0.
      RaySecondaryVectorZ(iElem) = 0.
    END IF ! RayElemPassedEnergyLoc2nd(iElem).GT.0

    ! Store element volume for later calculation of photon energy density
    ElemVolume(iElem) = ElemVolume_Shared(iCNElem)

    ! 2. Variable polynomial degree data
    Nloc = N_DG_Ray_loc(iElem)
    !U_N_Ray(iElem)%U(1:1,:,:,:) = RayElemPassedEnergy(1,iGlobalElem)
    !U_N_Ray(iElem)%U(2:2,:,:,:) = RayElemPassedEnergy(2,iGlobalElem)
#if USE_MPI
    IF(nProcessors.GT.1)THEN
      ! Get data from shared array
      NlocOffset = (Nloc+1)**3
      ASSOCIATE( i => RayElemOffset(iGlobalElem))
        U_N_Ray(iGlobalElem)%U(1,:,:,:) = RESHAPE(RayElemPassedEnergyHO_Shared(1, i+1:i+NlocOffset), (/Nloc+1, Nloc+1, Nloc+1/))
        U_N_Ray(iGlobalElem)%U(2,:,:,:) = RESHAPE(RayElemPassedEnergyHO_Shared(2, i+1:i+NlocOffset), (/Nloc+1, Nloc+1, Nloc+1/))
      END ASSOCIATE
    END IF ! nProcessors.GT.1
#endif /*USE_MPI*/

    ! Get the element volumes on Nloc
    ! Apply integration weights and the Jacobian
    ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the Jacobian ;-)
    J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
    ! p-refinement: Interpolate lower degree to higher degree (other way around would require model=T)
    J_Nloc = 0.
    IF(PP_N.EQ.Nloc)THEN
      J_Nloc(1,0:PP_N,0:PP_N,0:PP_N) = J_N(1,0:PP_N,0:PP_N,0:PP_N)
    ELSE
      ALLOCATE(Vdm_GaussN_Nloc(0:PP_N,0:Nloc))
      IF(Nloc.LT.PP_N)THEN
        CALL GetVandermonde(PP_N, NodeType, Nloc, Ray%NodeType, Vdm_GaussN_Nloc, modal=.TRUE.)
      ELSE
        CALL GetVandermonde(PP_N, NodeType, Nloc, Ray%NodeType, Vdm_GaussN_Nloc, modal=.FALSE.)
      END IF ! Nloc.LT.PP_N
        CALL ChangeBasis3D(1,PP_N,Nloc,Vdm_GaussN_Nloc,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_Nloc(1:1,0:Nloc,0:Nloc,0:Nloc))
    END IF ! PP_N.EQ.Nloc

    ! Calculate the sub-volumes
    DO m=0,Nloc
      DO l=0,Nloc
        DO k=0,Nloc
          IntegrationWeight = N_Inter_Ray(Nloc)%wGP(k)*&
                              N_Inter_Ray(Nloc)%wGP(l)*&
                              N_Inter_Ray(Nloc)%wGP(m)*J_Nloc(1,k,l,m)
          !UNMax(1:2,k,l,m,iElem) = UNMax(1:2,k,l,m,iElem) / IntegrationWeight
          U_N_Ray(iGlobalElem)%U(3,k,l,m) = IntegrationWeight
        END DO ! k
      END DO ! l
    END DO ! m
    IF(PP_N.NE.Nloc) DEALLOCATE(Vdm_GaussN_Nloc)

    ! Sanity checks: Low order
    ! 1.) compare sum of sub-volumes with cell-const value and abort
    ! 2.) compare sum of sub-cell energies with cell-const value and abort (1st and 2nd energies)
    ! 1st energy
    IF(RayElemPassedEnergyLoc1st(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(RayElemPassedEnergyLoc1st(iElem), SUM(U_N_Ray(iGlobalElem)%U(1,:,:,:)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem                    = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "RayElemPassedEnergyLoc1st(iElem)     = ", RayElemPassedEnergyLoc1st(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(U_N_Ray(iGlobalElem)%U(1,:,:,:)) = ", SUM(U_N_Ray(iGlobalElem)%U(1,:,:,:))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(U_N_Ray(iGlobalElem)%U(1,:,:,:))/RayElemPassedEnergyLoc1st(iElem)
        CALL abort(__STAMP__,'Before: RayElemPassedEnergyLoc1st does not match U_N_Ray%U(1) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF
    ! 2nd energy
    IF(RayElemPassedEnergyLoc2nd(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(RayElemPassedEnergyLoc2nd(iElem), SUM(U_N_Ray(iGlobalElem)%U(2,:,:,:)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem                    = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "RayElemPassedEnergyLoc2nd(iElem)     = ", RayElemPassedEnergyLoc2nd(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(U_N_Ray(iGlobalElem)%U(2,:,:,:)) = ", SUM(U_N_Ray(iGlobalElem)%U(2,:,:,:))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(U_N_Ray(iGlobalElem)%U(2,:,:,:))/RayElemPassedEnergyLoc2nd(iElem)
        CALL abort(__STAMP__,'Before: RayElemPassedEnergyLoc1st does not match U_N_Ray%U(2) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF
    ! volume
    IF(ElemVolume(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(ElemVolume(iElem), SUM(U_N_Ray(iGlobalElem)%U(3,:,:,:)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem                    = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "ElemVolume(iElem)                    = ", ElemVolume(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(U_N_Ray(iGlobalElem)%U(3,:,:,:)) = ", SUM(U_N_Ray(iGlobalElem)%U(3,:,:,:))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(U_N_Ray(iGlobalElem)%U(3,:,:,:))/ElemVolume(iElem)
        CALL abort(__STAMP__,'Before: ElemVolume(iElem) does not match U_N_Ray%U(3) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF

    ! Calculate the photon energy density on Nloc
    U_N_Ray(iGlobalElem)%U(4,:,:,:) = U_N_Ray(iGlobalElem)%U(1,:,:,:)/U_N_Ray(iGlobalElem)%U(3,:,:,:)
    U_N_Ray(iGlobalElem)%U(5,:,:,:) = U_N_Ray(iGlobalElem)%U(2,:,:,:)/U_N_Ray(iGlobalElem)%U(3,:,:,:)

    ! Map from Nloc to NMax for output to .h5 on the highest polynomial degree NMax
    ! The higher order element volume UNMax(3,:,:,:,:) is over-written later on
    IF(Nloc.EQ.Ray%Nmax)THEN
      UNMax(:,:,:,:,iElem) = U_N_Ray(iGlobalElem)%U(:,:,:,:)
    ELSE
      CALL ChangeBasis3D(nVarRay, Nloc, Ray%NMax, PREF_VDM_Ray(Nloc,Ray%NMax)%Vdm, U_N_Ray(iGlobalElem)%U(:,:,:,:), UNMax(:,:,:,:,iElem))
    END IF ! Nloc.Eq.Nmax

    ! Copy data from global array (later used for emission)
    U_N_Ray_loc(iElem)%U(:,:,:,:) = U_N_Ray(iGlobalElem)%U(:,:,:,:)

    ! Apply integration weights and the Jacobian
    ! Interpolate the Jacobian to the analyze grid: be careful we interpolate the inverse of the inverse of the Jacobian ;-)
    J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
    CALL ChangeBasis3D(1,PP_N,Ray%NMax,Vdm_GaussN_NMax,J_N(1:1,0:PP_N,0:PP_N,0:PP_N),J_Nmax(1:1,:,:,:))
    DO m=0,Ray%NMax
      DO l=0,Ray%NMax
        DO k=0,Ray%NMax
          IntegrationWeight = N_Inter_Ray(Ray%NMax)%wGP(k)*&
                              N_Inter_Ray(Ray%NMax)%wGP(l)*&
                              N_Inter_Ray(Ray%NMax)%wGP(m)*J_Nmax(1,k,l,m)
          !UNMax(1:2,k,l,m,iElem) = UNMax(1:2,k,l,m,iElem) / IntegrationWeight
          UNMax(3,k,l,m,iElem) = IntegrationWeight
        END DO ! k
      END DO ! l
    END DO ! m

    ! Sanity checks: High order
    ! 1.) compare sum of sub-volumes with cell-const value and abort
    ! 2.) compare sum of sub-cell energies with cell-const value and abort (1st and 2nd energies)
    ! 1st energy
    IF(RayElemPassedEnergyLoc1st(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(RayElemPassedEnergyLoc1st(iElem), SUM(UNMax(1,:,:,:,iElem)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem                = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "RayElemPassedEnergyLoc1st(iElem) = ", RayElemPassedEnergyLoc1st(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(UNMax(1,:,:,:,iElem))        = ", SUM(UNMax(1,:,:,:,iElem))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(UNMax(1,:,:,:,iElem))/RayElemPassedEnergyLoc1st(iElem)
        CALL abort(__STAMP__,'After: RayElemPassedEnergyLoc1st does not match UNMax(1) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF
    ! 2nd energy
    IF(RayElemPassedEnergyLoc2nd(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(RayElemPassedEnergyLoc2nd(iElem), SUM(UNMax(2,:,:,:,iElem)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem                = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "RayElemPassedEnergyLoc2nd(iElem) = ", RayElemPassedEnergyLoc2nd(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(UNMax(2,:,:,:,iElem))        = ", SUM(UNMax(2,:,:,:,iElem))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(UNMax(1,:,:,:,iElem))/RayElemPassedEnergyLoc2nd(iElem)
        CALL abort(__STAMP__,'After: RayElemPassedEnergyLoc1st does not match UNMax(2) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF
    ! volume
    IF(ElemVolume(iElem).GT.0.0)THEN
      IF(.NOT.ALMOSTEQUALRELATIVE(ElemVolume(iElem), SUM(UNMax(3,:,:,:,iElem)), tolerance))THEN
        IPWRITE(UNIT_StdOut,*) "iElem,iGlobalElem         = ", iElem,iGlobalElem
        IPWRITE(UNIT_StdOut,*) "ElemVolume(iElem)         = ", ElemVolume(iElem)
        IPWRITE(UNIT_StdOut,*) "SUM(UNMax(3,:,:,:,iElem)) = ", SUM(UNMax(3,:,:,:,iElem))
        IPWRITE(UNIT_StdOut,*) "ratio =", SUM(UNMax(3,:,:,:,iElem))/ElemVolume(iElem)
        CALL abort(__STAMP__,'After: ElemVolume(iElem) does not match UNMax(3) for tolerance = ',RealInfoOpt=tolerance)
      END IF
    END IF

  END DO ! iElem=1,PP_nElems

  ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  ! Write file after last abort to prevent a corrupt output file (which might be used when restarting the simulation)
  IF(MPIRoot) CALL GenerateFileSkeleton('RadiationVolState',nVarRay,StrVarNames,TRIM(MeshFile),0.,FileNameIn=RadiationVolState,NIn=Ray%NMax,NodeType_in=Ray%NodeType)
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

  ! Sanity check
  IF(ANY(ISNAN(UNMax))) CALL abort(__STAMP__,'Found one or more NaN in the array UNMax!')

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nVarRay           => INT(nVarRay,IK)            ,&
        NMax              => INT(Ray%NMax,IK)           ,&
        nGlobalElems      => INT(nGlobalElems,IK)       ,&
        PP_nElems         => INT(PP_nElems,IK)          ,&
        offsetElem        => INT(offsetElem,IK)         )
    CALL GatheredWriteArray(RadiationVolState,create=.FALSE.,&
         DataSetName='DG_Solution', rank=5,&
         nValGlobal=(/nVarRay     , NMax+1_IK , NMax+1_IK , NMax+1_IK , nGlobalElems/) , &
         nVal=      (/nVarRay     , NMax+1_IK , NMax+1_IK , NMax+1_IK , PP_nElems/)    , &
         offset=    (/0_IK        , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
         collective=.TRUE.,RealArray=UNMax)
  END ASSOCIATE
#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

! Write all 'ElemData' arrays to a single container in the state.h5 file
CALL WriteAdditionalElemData(RadiationVolState,ElementOut)

! Write Nloc and the reflected vector to separate containers for element- and process-wise read-in during restart or loadbalance
CALL WriteElemDataToSeparateContainer(RadiationVolState,ElementOut,'Nloc')
CALL WriteElemDataToSeparateContainer(RadiationVolState,ElementOut,'RaySecondaryVectorX')
CALL WriteElemDataToSeparateContainer(RadiationVolState,ElementOut,'RaySecondaryVectorY')
CALL WriteElemDataToSeparateContainer(RadiationVolState,ElementOut,'RaySecondaryVectorZ')

SWRITE(*,*) 'DONE'
END SUBROUTINE WritePhotonVolSampleToHDF5


SUBROUTINE WritePhotonSurfSampleToHDF5()
!===================================================================================================================================
!> write the final values of the surface sampling to a HDF5 state file
!> additional performs all the final required computations
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Particle_Boundary_Vars ,ONLY: nComputeNodeSurfOutputSides,nGlobalOutputSides, nSurfBC
USE MOD_Particle_Boundary_Vars ,ONLY: offsetComputeNodeSurfOutputSide, SurfBCName, nComputeNodeSurfSides
USE MOD_Particle_Boundary_Vars ,ONLY: SurfSide2GlobalSide, GlobalSide2SurfSide
USE MOD_HDF5_Output            ,ONLY: WriteAttributeToHDF5,WriteArrayToHDF5,WriteHDF5Header
USE MOD_Mesh_Vars              ,ONLY: MeshFile
USE MOD_Particle_Mesh_Vars     ,ONLY: SideInfo_Shared
USE MOD_MPI_Shared_Vars        ,ONLY: mySurfRank
#if USE_MPI
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_LEADERS_SURF
USE MOD_Particle_Boundary_Vars ,ONLY: nGlobalSurfSides
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSurfSideArea_Shared
#else
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSurfSideArea
#endif /*USE_MPI*/
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSampWall
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Photon_TrackingVars    ,ONLY: RadiationSurfState
USE MOD_RayTracing_Vars        ,ONLY: Ray
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                  :: Statedummy
CHARACTER(LEN=255)                  :: H5_Name, H5_Name2
CHARACTER(LEN=255),ALLOCATABLE      :: Str2DVarNames(:)
INTEGER                             :: GlobalSideID, iSurfSide, OutputCounter, SurfSideNb, p, q
INTEGER,PARAMETER                   :: nVar2D=3
REAL                                :: tstart,tend
REAL, ALLOCATABLE                   :: helpArray(:,:,:,:)
INTEGER, ALLOCATABLE                :: helpArray2(:)
!===================================================================================================================================
#if USE_MPI
CALL ExchangeRadiationSurfData()
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.EQ.MPI_COMM_NULL) RETURN
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)

! Return if no sampling sides
IF (nGlobalSurfSides.EQ.0) RETURN
#endif /*USE_MPI*/
IF (mySurfRank.EQ.0) THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE Radiation SurfSTATE TO HDF5 FILE...'
  tstart=LOCALTIME()
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF (mySurfRank.EQ.0) THEN
#endif /*USE_MPI*/
  CALL OpenDataFile(RadiationSurfState,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
  Statedummy = 'RadiationSurfState'
  ! Write file header
  CALL WriteHDF5Header(Statedummy,File_ID)
  CALL WriteAttributeToHDF5(File_ID , 'DSMC_nSurfSample' , 1       , IntegerScalar = Ray%nSurfSample        )
  CALL WriteAttributeToHDF5(File_ID , 'MeshFile'         , 1       , StrScalar     = (/TRIM(MeshFile)/) )
  CALL WriteAttributeToHDF5(File_ID , 'BC_Surf'          , nSurfBC , StrArray      = SurfBCName         )
  CALL WriteAttributeToHDF5(File_ID , 'N'                , 1       , IntegerScalar = Ray%nSurfSample        )
  CALL WriteAttributeToHDF5(File_ID , 'NodeType'         , 1       , StrScalar     = (/Ray%NodeType/)   )
  CALL WriteAttributeToHDF5(File_ID , 'Time'             , 1       , RealScalar    = 0.                 )

  ALLOCATE(Str2DVarNames(1:nVar2D))
  ! fill varnames for total values
  Str2DVarNames(1) ='PhotonCount'
  Str2DVarNames(2) ='HeatFlux'
  Str2DVarNames(3) ='iBC'

  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurface',nVar2D,StrArray=Str2DVarNames)

   CALL CloseDataFile()
  DEALLOCATE(Str2DVarNames)
#if USE_MPI
END IF
CALL MPI_BARRIER(MPI_COMM_LEADERS_SURF,iERROR)
CALL OpenDataFile(RadiationSurfState,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
#else
CALL OpenDataFile(RadiationSurfState,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
#endif /*USE_MPI*/


WRITE(H5_Name,'(A)') 'SurfaceData'
WRITE(H5_Name2,'(A)') 'SurfaceDataGlobalSideIndex'
#if USE_MPI
ASSOCIATE(PhotonSurfSideArea   => PhotonSurfSideArea_Shared)
#endif

ASSOCIATE (&
      nSurfSample    => INT(Ray%nSurfSample                 , IK)  , &
      nGlobalSides   => INT(nGlobalOutputSides              , IK)  , &
      LocalnBCSides  => INT(nComputeNodeSurfOutputSides     , IK)  , &
      offsetSurfSide => INT(offsetComputeNodeSurfOutputSide , IK)  , &
      nVar2D         => INT(nVar2D                          , IK))

  ALLOCATE(helpArray(nVar2D,1:nSurfSample,1:nSurfSample,LocalnBCSides))
  ALLOCATE(helpArray2(LocalnBCSides))
  OutputCounter = 0
  DO iSurfSide = 1,nComputeNodeSurfSides
    GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    IF(SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID).GT.0) THEN
      IF(GlobalSideID.LT.SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID)) THEN
        SurfSideNb = GlobalSide2SurfSide(SURF_SIDEID,SideInfo_Shared(SIDE_NBSIDEID,GlobalSideID))
        ! Add your contribution to my inner BC
        PhotonSampWall(:,:,:,iSurfSide) = PhotonSampWall(:,:,:,iSurfSide) + PhotonSampWall(:,:,:,SurfSideNb)
      ELSE
        CYCLE
      END IF
    END IF
    OutputCounter = OutputCounter + 1
    helpArray(1,1:nSurfSample,1:nSurfSample,OutputCounter) = PhotonSampWall(1,1:nSurfSample,1:nSurfSample,iSurfSide)
    helpArray2(OutputCounter) = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    !  SurfaceArea should be changed to 1:SurfMesh%nSides if inner sampling sides exist...
    DO p = 1, INT(nSurfSample)
      DO q = 1, INT(nSurfSample)
        helpArray(2,p,q,OutputCounter) = PhotonSampWall(2,p,q,iSurfSide)/PhotonSurfSideArea(p,q,iSurfSide)
        helpArray(3,p,q,OutputCounter) = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))
      END DO ! q = 1, nSurfSample
    END DO ! p = 1, nSurfSample
  END DO
  ! WARNING: Only the sampling leaders write the data to .h5
  CALL WriteArrayToHDF5(DataSetName=H5_Name  , rank=4      ,                                  &
                        nValGlobal =(/nVar2D , nSurfSample , nSurfSample , nGlobalSides/)   , &
                        nVal       =(/nVar2D , nSurfSample , nSurfSample , LocalnBCSides/)  , &
                        offset     =(/0_IK   , 0_IK        , 0_IK        , offsetSurfSide/) , &
                        collective =.FALSE.  ,                                                &
                        RealArray=helpArray(1:nVar2D,1:nSurfSample,1:nSurfSample,1:LocalnBCSides))
  CALL WriteArrayToHDF5(DataSetName = H5_Name2             , rank = 1 , &
                        nValGlobal  = (/nGlobalSides/)     , &
                        nVal        = (/LocalnBCSides/)      , &
                        offset      = (/offsetSurfSide/) , &
                        collective  = .FALSE.  , IntegerArray_i4 = helpArray2(1:INT(LocalnBCSides,4)))
  DEALLOCATE(helpArray)
  DEALLOCATE(helpArray2)
END ASSOCIATE

#if USE_MPI
END ASSOCIATE
#endif /*USE_MPI*/

CALL CloseDataFile()

IF (mySurfRank.EQ.0) THEN
  tend=LOCALTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',tend-tstart,'s]'
END IF

END SUBROUTINE WritePhotonSurfSampleToHDF5

#if USE_MPI
SUBROUTINE ExchangeRadiationSurfData()
!===================================================================================================================================
! exchange the surface data
! only processes with samling sides in their halo region and the original process participate on the communication
! structure is similar to particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars ,ONLY: SurfTotalSideOnNode, SurfMapping, nComputeNodeSurfTotalSides, GlobalSide2SurfSide
USE MOD_Particle_MPI_Vars      ,ONLY: SurfSendBuf,SurfRecvBuf
USE MOD_Photon_TrackingVars    ,ONLY: PhotonSampWallProc, PhotonSampWall_Shared, PhotonSampWall_Shared_Win
USE MOD_RayTracing_Vars        ,ONLY: Ray
USE MOD_MPI_Shared_Vars        ,ONLY: MPI_COMM_LEADERS_SURF, MPI_COMM_SHARED, nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_MPI_Shared
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID, SideID
INTEGER                         :: iPos,iProc,p,q
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

MessageSize = 2*nComputeNodeSurfTotalSides*(Ray%nSurfSample**2)
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(PhotonSampWallProc, PhotonSampWall_Shared, MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
ELSE
  CALL MPI_REDUCE(PhotonSampWallProc, 0                    , MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
ENDIF
! Nullify process-local array
PhotonSampWallProc = 0.0

! Update
CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win,MPI_COMM_SHARED)

! prepare buffers for surf leader communication
IF (myComputeNodeRank.EQ.0) THEN
  nValues = 2

  ! open receive buffer
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvSurfSides * nValues
    CALL MPI_IRECV( SurfRecvBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE
    ! Only assemble message if we are expecting sides to send to this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    SurfSendBuf(iProc)%content = 0.
    DO iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      SideID     = SurfMapping(iProc)%SendSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      ! Assemble message
      DO q = 1,Ray%nSurfSample
        DO p = 1,Ray%nSurfSample
        SurfSendBuf(iProc)%content(iPos+1:iPos+2) = PhotonSampWall_Shared(:,p,q,SurfSideID)
        iPos = iPos + 2
        END DO ! p=0,Ray%nSurfSample
      END DO ! q=0,Ray%nSurfSample
      PhotonSampWall_Shared(:,:,:,SurfSideID)=0.
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE
    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendSurfSides * nValues
    CALL MPI_ISEND( SurfSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE
    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
      SideID     = SurfMapping(iProc)%RecvSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      DO q = 1,Ray%nSurfSample
        DO p = 1,Ray%nSurfSample
          PhotonSampWall_Shared(:,p,q,SurfSideID) = PhotonSampWall_Shared(:,p,q,SurfSideID) &
                                                  + SurfRecvBuf(iProc)%content(iPos+1:iPos+2)
          iPos = iPos + 2
        END DO ! p=0,Ray%nSurfSample
      END DO ! q=0,Ray%nSurfSample
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
     ! Nullify buffer
    SurfRecvBuf(iProc)%content = 0.
  END DO ! iProc
END IF

CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win,MPI_COMM_SHARED)

END SUBROUTINE ExchangeRadiationSurfData


!===================================================================================================================================
! Exchanges and add up the volume ray tracing data between all processes to have the global data available on each process
! 1. Exchange the low-order data
! 2. Exchange the high-order data
!===================================================================================================================================
SUBROUTINE ExchangeRayVolInfo()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_RayTracing_Vars ,ONLY: RayElemPassedEnergy,RayElemPassedEnergy_Shared,RayElemPassedEnergy_Shared_Win,RayElemSize
USE MOD_RayTracing_Vars ,ONLY: nVarRay,RayElemPassedEnergyHO_Shared,RayElemPassedEnergyHO_Shared_Win,N_DG_Ray,U_N_Ray,RayElemOffset
USE MOD_Mesh_Vars       ,ONLY: nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER           :: MessageSize,offset,NlocOffset
INTEGER(KIND=8)   :: nGlobalEntries
REAL, ALLOCATABLE :: RayElemPassedEnergyHO(:,:) ! <
INTEGER           :: iElem,Nloc
CHARACTER(LEN=255):: hilf
!===================================================================================================================================
! Collect the information from the process-local shadow arrays in the compute-node shared array
MessageSize = RayElemSize*nGlobalElems

! 1. Exchange the low-order data
! Reduce data to each node leader
IF (myComputeNodeRank.EQ.0) THEN
  CALL MPI_REDUCE(RayElemPassedEnergy, RayElemPassedEnergy_Shared, MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
ELSE
  CALL MPI_REDUCE(RayElemPassedEnergy, 0                         , MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
ENDIF
CALL BARRIER_AND_SYNC(RayElemPassedEnergy_Shared_Win, MPI_COMM_SHARED)

! Synchronize data between node leaders with all-reduce
IF(nLeaderGroupProcs.GT.1)THEN
  IF(myComputeNodeRank.EQ.0)THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RayElemPassedEnergy_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,iError)
  END IF

  ! Synchronize data from node leaders to node workers
  CALL BARRIER_AND_SYNC(RayElemPassedEnergy_Shared_Win, MPI_COMM_SHARED)
END IF

! 2. Exchange the high-order data
! Only duplicate and reduce the data when more than one process are used
IF(nProcessors.GT.1)THEN

  nGlobalEntries = 0
  DO iElem = 1, nGlobalElems
    Nloc = N_DG_Ray(iElem)
    nGlobalEntries = nGlobalEntries + INT((Nloc+1)**3,8)
  END DO ! iElem = 1, nGlobalElems

  ! Sanity check
  IF(nGlobalEntries * INT(nVarRay,8).GT.INT(HUGE(1_4),8))THEN
    IF(MPIRoot)THEN
      WRITE(UNIT=hilf,FMT='(A,I0,A,I0)') "Number of entries in RayElemPassedEnergyHO(1:nVarRay,1:nGlobalEntries) "&
          ,nGlobalEntries * INT(nVarRay,8)," is larger than ",HUGE(1_4)
      CALL abort(__STAMP__,TRIM(hilf))
    END IF
  END IF

  ALLOCATE(RayElemPassedEnergyHO(nVarRay,nGlobalEntries))
  RayElemPassedEnergyHO=0.
  ALLOCATE(RayElemOffset(nGlobalElems))
  !> Shared arrays for high-order volume sampling
  CALL Allocate_Shared((/nVarRay,INT(nGlobalEntries,4)/),RayElemPassedEnergyHO_Shared_Win,RayElemPassedEnergyHO_Shared)
  CALL MPI_WIN_LOCK_ALL(0,RayElemPassedEnergyHO_Shared_Win,IERROR)
  CALL BARRIER_AND_SYNC(RayElemPassedEnergyHO_Shared_Win,MPI_COMM_SHARED)

  ! Store data in local array
  offset = 0
  DO iElem = 1, nGlobalElems
    Nloc = N_DG_Ray(iElem)
    NlocOffset = (Nloc+1)**3
    RayElemOffset(iElem) = offset
    RayElemPassedEnergyHO(1,offset+1:offset+NlocOffset) = RESHAPE(U_N_Ray(iElem)%U(1,:,:,:),(/(Nloc+1)**3/))
    RayElemPassedEnergyHO(2,offset+1:offset+NlocOffset) = RESHAPE(U_N_Ray(iElem)%U(2,:,:,:),(/(Nloc+1)**3/))
    offset = offset + NlocOffset
  END DO ! iElem = 1, nGlobalElems

  MessageSize = nVarRay * INT(nGlobalEntries,4)

  IF (myComputeNodeRank.EQ.0) THEN
    CALL MPI_REDUCE(RayElemPassedEnergyHO, RayElemPassedEnergyHO_Shared, MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
  ELSE
    CALL MPI_REDUCE(RayElemPassedEnergyHO, 0                           , MessageSize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_SHARED, IERROR)
  ENDIF
  CALL BARRIER_AND_SYNC(RayElemPassedEnergyHO_Shared_Win, MPI_COMM_SHARED)

  ! Synchronize data between node leaders with all-reduce
  IF(nLeaderGroupProcs.GT.1)THEN
    IF(myComputeNodeRank.EQ.0)THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RayElemPassedEnergyHO_Shared,MessageSize,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_LEADERS_SHARED,iError)
    END IF

    ! Synchronize data from node leaders to node workers
    CALL BARRIER_AND_SYNC(RayElemPassedEnergyHO_Shared_Win, MPI_COMM_SHARED)
  END IF

  DEALLOCATE(RayElemPassedEnergyHO)

END IF ! nProcessors.GT.1

END SUBROUTINE ExchangeRayVolInfo
#endif /*USE_MPI*/

END MODULE MOD_Photon_TrackingOutput
