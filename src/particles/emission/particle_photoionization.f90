!==================================================================================================================================
! Copyright (c) 2023-2025 boltzplatz - numerical plasma dynamics GmbH
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

MODULE MOD_Particle_Photoionization
!===================================================================================================================================
!> Module for particle insertion through photo-ionization
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: PhotoIonization_RayTracing_SEE, PhotoIonization_RayTracing_Volume
PUBLIC :: CalcPhotoIonizationNumber, PhotoIonization_InsertProducts
!===================================================================================================================================
CONTAINS

!===================================================================================================================================
!> Calculate the time-dependent scaling factor for the maximum intensity, depending on the pulse type
!===================================================================================================================================
SUBROUTINE GetPulseIntensity(TimeScalingFactor)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Timedisc_Vars           ,ONLY: dt,time
USE MOD_RayTracing_Vars         ,ONLY: Ray
#ifdef LSERK
USE MOD_Timedisc_Vars           ,ONLY: iStage, RK_c, nRKStages
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)     :: TimeScalingFactor        !< Scaling factor of the maximum intensity: I(t) = I_max * factor
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: t_1, t_2
INTEGER               :: NbrOfRepetitions
!===================================================================================================================================
! Initialize
TimeScalingFactor = 0.

! Determine the time interval
#ifdef LSERK
IF (iStage.EQ.1) THEN
  t_1 = Time
  t_2 = Time + RK_c(2) * dt
ELSE
  IF (iStage.NE.nRKStages) THEN
    t_1 = Time + RK_c(iStage) * dt
    t_2 = Time + RK_c(iStage+1) * dt
  ELSE
    t_1 = Time + RK_c(iStage) * dt
    t_2 = Time + dt
  END IF
END IF
#else
t_1 = Time
t_2 = Time + dt
#endif

! Calculate the current pulse, starting at zero
NbrOfRepetitions = INT(Time/Ray%Period)

! Leave the subroutine if the maximum number of pulses has been reached
IF((NbrOfRepetitions+1).GT.Ray%NbrOfPulses) THEN
  TimeScalingFactor = 0.
  RETURN
END IF

! Gaussian pulse only: Add time shift of -SQRT(8) * Ray%PulseDuration so that I_max is after half of the pulse duration
! For the square pulse, this reduces to dt, Ray%tShift = 0.
! In case of multiple pulses, shift the time interval accordingly
t_1 = t_1 - Ray%tShift - NbrOfRepetitions * Ray%Period
t_2 = t_2 - Ray%tShift - NbrOfRepetitions * Ray%Period

! Calculate the time-dependent ray intensity
SELECT CASE(Ray%PulseType)
CASE('Gaussian')
  ! check if t_2 is outside of the pulse
  IF(t_2.GT.2.0*Ray%tShift) t_2 = 2.0*Ray%tShift
  ! Determine the time scaling factor
  TimeScalingFactor = 0.5 * SQRT(PI) * Ray%PulseDuration * (ERF(t_2/Ray%PulseDuration)-ERF(t_1/Ray%PulseDuration))
CASE('square')
  IF(t_1.LT.Ray%PulseDuration) THEN
    ! check if t_2 is outside of the pulse
    IF(t_2.GT.Ray%PulseDuration) t_2 = Ray%PulseDuration
    TimeScalingFactor = t_2 - t_1
  ELSE
    TimeScalingFactor = 0.
  END IF
CASE DEFAULT
  CALL Abort(__STAMP__,'Unknown pulse type for ray tracing: '//TRIM(Ray%PulseType)//'. Select square or Gaussian!')
END SELECT

! Additional scaling in case the power density has been changed
TimeScalingFactor = TimeScalingFactor * Ray%IntensityAmplitudeFactor

END SUBROUTINE GetPulseIntensity


SUBROUTINE PhotoIonization_RayTracing_SEE()
!===================================================================================================================================
!> Routine calculates the number of secondary electrons to be emitted and inserts them on the surface, utilizing the cell-local
!> photon energy from the raytracing
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Particle_Boundary_Vars  ,ONLY: Partbound,DoBoundaryParticleOutputRay
USE MOD_Particle_Vars           ,ONLY: Species, PartState, usevMPF,SpeciesOffsetVDL
USE MOD_RayTracing_Vars         ,ONLY: Ray,UseRayTracing,RayElemEmission
USE MOD_part_emission_tools     ,ONLY: CalcPhotonEnergy
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared,UseBezierControlPoints
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D, BezierSampleXi
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient, CalcNormAndTangBezier
USE MOD_Mesh_Vars               ,ONLY: NGeo,offsetElem,nElems
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_FromWorkFuncSEE
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_part_operations         ,ONLY: CreateParticle
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfSides, SurfSide2GlobalSide
#ifdef LSERK
USE MOD_Timedisc_Vars           ,ONLY: RK_c, nRKStages
#endif
USE MOD_Photon_TrackingVars     ,ONLY: PhotonSampWall_loc,PhotonSurfSideArea
#if USE_HDG
USE MOD_HDG_Vars                ,ONLY: UseFPC,FPC,UseEPC,EPC
USE MOD_Mesh_Vars               ,ONLY: BoundaryType
USE MOD_Particle_Boundary_Vars  ,ONLY: DoVirtualDielectricLayer,Do2DSurfaceCharge
USE MOD_Particle_Vars           ,ONLY: LastPartPos,PartSpecies
USE MOD_PICDepo_Tools           ,ONLY: DepositParticleOnSurface
#if USE_PETSC
USE MOD_HDG_Vars                ,ONLY: UseCircuitModel
#endif /*USE_PETSC*/
#endif /*USE_HDG*/
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SEE,CalcPhotonSEE
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemBaryNGeo
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_Particle_Vars           ,ONLY: nSpecies
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: E_Intensity,vec(3)
INTEGER               :: SideID, GlobElemID, PartID, locElemID, iSurfSide, CNElemID
INTEGER               :: p, q, iPartBound, SpecID, iPart, NbrOfSEE, iSEEBC
REAL                  :: RealNbrOfSEE, TimeScalingFactor, MPF,PhotonEnergy
REAL                  :: PartPos(1:3), PartPosSurf(1:3), xi(2)
REAL                  :: RandVal, RandVal2(2), xiab(1:2,1:2), nVec(3), tang1(3), tang2(3), Velo3D(3)
#if USE_HDG
REAL                  :: ChargeHole
INTEGER               :: iBC,iUniqueFPCBC,iUniqueEPCBC,BCState
#endif /*USE_HDG*/
#if USE_LOADBALANCE
REAL                  :: tLBStart
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================
! Check if ray tracing based SEE is active
! 1) Boundary from which rays are emitted
IF(.NOT.UseRayTracing) RETURN
! 2) SEE yield for any BC greater than zero
IF(.NOT.ANY(PartBound%PhotonSEEYield(:).GT.0.)) RETURN

! Get the time-dependent pulse intensity factor based on the pulse type (square, Gaussian, etc.)
CALL GetPulseIntensity(TimeScalingFactor)

! Get the photon energy (currently: single, fixed wave length -> calculate energy once)
PhotonEnergy = CalcPhotonEnergy(Ray%WaveLength)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! Loop over all surface sides (to include inner BCs, which are not part of nBCSides)
DO iSurfSide = 1, nComputeNodeSurfSides
  SideID    = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  ! Determine which element the particles are going to be inserted
  GlobElemID = SideInfo_Shared(SIDE_ELEMID ,SideID)
  locElemID = GlobElemID - offsetElem
  ! Cycle non-local elements
  IF((locElemID.LE.0).OR.(locElemID.GT.nElems)) CYCLE
  ! Skip elements without ionization
  IF(.NOT.RayElemEmission(1,locElemID)) CYCLE
  ! Skip non-reflective BC sides
  iPartBound = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  IF(PartBound%TargetBoundCond(iPartBound).NE.PartBound%ReflectiveBC) CYCLE
  ! Skip BC sides with zero yield
  IF(PartBound%PhotonSEEYield(iPartBound).LE.0.) CYCLE
  ! Determine which species is to be inserted
  SpecID = PartBound%PhotonSEEElectronSpecies(iPartBound)
  ! Sanity check
  IF(SpecID.EQ.0)THEN
    IPWRITE(UNIT_StdOut,*) "iPartBound =", iPartBound
    IPWRITE(UNIT_StdOut,*) "PartBound%PhotonSEEElectronSpecies(iPartBound) =", PartBound%PhotonSEEElectronSpecies(iPartBound)
    CALL abort(__STAMP__,'Electron species index cannot be zero!')
  END IF ! SpecID.eq.0
  IF(SpecID.GT.nSpecies)THEN
    IPWRITE(UNIT_StdOut,*) "iPartBound =", iPartBound
    IPWRITE(UNIT_StdOut,*) "PartBound%PhotonSEEElectronSpecies(iPartBound) =", PartBound%PhotonSEEElectronSpecies(iPartBound)
    IPWRITE(UNIT_StdOut,*) "nSpecies   =", nSpecies
    CALL abort(__STAMP__,'Electron species index cannot greater than nSpecies!')
  END IF ! SpecID.eq.0
  ! Determine the weighting factor of the electron species
  IF(usevMPF)THEN
    MPF = PartBound%PhotonSEEMacroParticleFactor(iPartBound) ! Use SEE-specific MPF
  ELSE
    MPF = Species(SpecID)%MacroParticleFactor ! Use species MPF
  END IF ! usevMPF
  ! Loop over the subsides
  DO p = 1, Ray%nSurfSample
    DO q = 1, Ray%nSurfSample
      IF(PhotonSampWall_loc(p,q,iSurfSide).LT.0.0) CALL abort(__STAMP__,'ERROR in PhotoIonization_RayTracing_SEE: PhotonSampWall_loc not defined!')
      ! Calculate the number of SEEs per subside
      E_Intensity = PhotonSampWall_loc(p,q,iSurfSide) * PhotonSurfSideArea(p,q,iSurfSide) * TimeScalingFactor
      RealNbrOfSEE = E_Intensity / PhotonEnergy * PartBound%PhotonSEEYield(iPartBound) / MPF
      ! Add random number to calculated real/float value of SEE particles and user INT()for lower-bound cut-off
      CALL RANDOM_NUMBER(RandVal)
      NbrOfSEE = INT(RealNbrOfSEE+RandVal)
      ! NbrOfSEE: Number of particles to be inserted
      IF(NbrOfSEE.GT.0)THEN
        ! Check if photon SEE electric current is to be measured
        IF(CalcPhotonSEE)THEN
          ! Note that the negative value of the charge -q is used below
          iSEEBC = SEE%BCIDToSEEBCID(iPartBound)
          SEE%RealElectronOutPhoton(iSEEBC) = SEE%RealElectronOutPhoton(iSEEBC) - MPF*NbrOfSEE*Species(SpecID)%ChargeIC
        END IF ! (NbrOfSEE.GT.0).AND.(CalcPhotonSEE)
        ! Calculate the normal & tangential vectors
        IF(UseBezierControlPoints)THEN
          ! Use Bezier polynomial
          xi(1)=(BezierSampleXi(p-1)+BezierSampleXi(p))/2. ! (a+b)/2
          xi(2)=(BezierSampleXi(q-1)+BezierSampleXi(q))/2. ! (a+b)/2
          xiab(1,1:2)=(/BezierSampleXi(p-1),BezierSampleXi(p)/)
          xiab(2,1:2)=(/BezierSampleXi(q-1),BezierSampleXi(q)/)
          CALL CalcNormAndTangBezier(nVec,tang1,tang2,xi(1),xi(2),SideID)
        ELSE
          ! Sanity check
          CALL abort(__STAMP__,'Photoionization with ray tracing requires BezierControlPoints3D')
        END IF ! UseBezierControlPoints
        ! Normal vector provided by the routine points outside of the domain
        nVec = -nVec
        ! Loop over number of particles to be inserted
        DO iPart = 1, NbrOfSEE
          ! Determine particle position within the sub-side
          CALL RANDOM_NUMBER(RandVal2)
          IF(UseBezierControlPoints)THEN
            ! Use Bezier polynomial
            xi=(xiab(:,2)-xiab(:,1))*RandVal2+xiab(:,1)
            CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=PartPosSurf(1:3))
          ELSE
            ! Sanity check
            CALL abort(__STAMP__,'Photoionization with ray tracing requires BezierControlPoints3D')
          END IF ! UseBezierControlPoints
          ! Determine particle velocity
          CALL CalcVelocity_FromWorkFuncSEE(PartBound%PhotonSEEWorkFunction(iPartBound), Species(SpecID)%MassIC, tang1, nVec, Velo3D)
          ! Move particle slightly into the domain away from the surface because TriaTracking loses the particle during restart
          ! as the InsideQuad3D test returns "not inside" for these particles and they are deleted
          CNElemID = GetCNElemID(GlobElemID)
          vec(1:3) = ElemBaryNGeo(1:3,CNElemID) - PartPosSurf(1:3)
          PartPos(1:3) = PartPosSurf(1:3) + 1e-7 * vec(1:3)
          ! Create new particle
          ! Create with PEM%LastGlobalElemID = 0 to prevent tracking directly after creation
          CALL CreateParticle(SpecID,PartPos(1:3),GlobElemID,0,Velo3D(1:3),0.,0.,0.,NewPartID=PartID,NewMPF=MPF)
          ! 1. Store the particle information in PartStateBoundary.h5
          IF(DoBoundaryParticleOutputRay) CALL StoreBoundaryParticleProperties(PartID,SpecID,PartState(1:3,PartID),&
                                                   UNITVECTOR(PartState(4:6,PartID)),nVec,iPartBound=iPartBound,mode=2,MPF_optIN=MPF)
#if USE_HDG
          ! 2. Check if floating boundary conditions (FPC) are used and consider electron holes
          IF(UseFPC)THEN
            iBC = PartBound%MapToFieldBC(iPartBound)
            IF(iBC.LE.0) CALL abort(__STAMP__,'UseFPC=T: iBC = PartBound%MapToFieldBC(PartBCIndex) must be >0',IntInfoOpt=iBC)
            IF(BoundaryType(iBC,BC_TYPE).EQ.20)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
              BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
              iUniqueFPCBC = FPC%Group(BCState,2)
              FPC%ChargeProc(iUniqueFPCBC) = FPC%ChargeProc(iUniqueFPCBC) - Species(SpecID)%ChargeIC * MPF ! Use negative charge!
            END IF ! BCType.EQ.20
          END IF ! UseFPC

          ! 3. Check if electric potential condition (EPC) are used and consider electron holes
          IF(UseEPC)THEN
            iBC = PartBound%MapToFieldBC(iPartBound)
            IF(iBC.LE.0) CALL abort(__STAMP__,'UseEPC=T: iBC = PartBound%MapToFieldBC(PartBCIndex) must be >0',IntInfoOpt=iBC)
            IF(BoundaryType(iBC,BC_TYPE).EQ.8)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
              BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
              iUniqueEPCBC = EPC%Group(BCState,2)
              EPC%ChargeProc(iUniqueEPCBC) = EPC%ChargeProc(iUniqueEPCBC) - Species(SpecID)%ChargeIC * MPF ! Use negative charge!
            END IF ! BCType.EQ.8
          END IF ! UseEPC

#if USE_PETSC
          ! 4. Check if circuit model boundary condition (CMBC) are used and consider electron holes
          IF (UseCircuitModel) THEN
            CALL abort(__STAMP__,'Circuit model (CMBC) not implemented for electron holes due to photon emission')
          END IF ! UseCircuitModel
#endif /*USE_PETSC*/

          ! 3. Check if SEE holes are to be deposited
          ! 3a. VDL
          IF(DoVirtualDielectricLayer)THEN
            IF(ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0)THEN
              ! Set velocity to zero as these virtual particles are deleted after the tracking/MPI communication step
              Velo3D(1:3) = 0.

              ! Create particle (with opposite charge by setting a nagative species index)
              ! Create with PEM%LastGlobalElemID = GlobElemID to trigger tracking directly after creation to find a possible new host
              ! element
              CALL CreateParticle(SpecID,PartPosSurf(1:3),GlobElemID,GlobElemID,Velo3D(1:3),0.,0.,0.,NewPartID=PartID,NewMPF=MPF)

              ! Set new particle position: Shift away from face by ratio of real VDL thickness and permittivity.
              ! Note the negative sign that is required as the normal vector points outwards has already been applied above
              PartState(1:3,PartID) = LastPartPos(1:3,PartID) &
                                    + (PartBound%ThicknessVDL(iPartBound)/PartBound%PermittivityVDL(iPartBound))*nVec

              ! Encode species index: Set a temporarily invalid number, which holds the information that the particle has interacted with a VDL.
              ! The Particle is removed after MPI communication because the new position might be on a different process due to the displacement
              PartSpecies(PartID) = PartSpecies(PartID) + SpeciesOffsetVDL
              ! Invert species index
              PartSpecies(PartID) = -PartSpecies(PartID)
            END IF ! ABS(PartBound%PermittivityVDL(iPartBound)).GT.0.0
          END IF ! DoVirtualDielectricLayer

          ! 3b. 2D surface charge deposition
          IF(Do2DSurfaceCharge) THEN
            iBC = PartBound%MapToFieldBC(iPartBound)
            IF(iBC.LE.0) CALL abort(__STAMP__,'Do2DSurfaceCharge=T: iBC = PartBound%MapToFieldBC(PartBCIndex) must be >0',IntInfoOpt=iBC)
            ! Method 3: 2D Surface Charging
            IF (PartBound%UseSurfaceCharge(iBC)) THEN
              ! Calculate the opposite charge
              ChargeHole = -Species(SpecID)%ChargeIC*MPF
              ! Deposit the charge(s). Use the previously calcualted reference coordinates xi(1:2) directly
              CALL DepositParticleOnSurface(ChargeHole, PartPosSurf(1:3), GlobElemID, SideID, PartID=0, xi=xi(1:2))
            END IF ! PartBound%UseSurfaceCharge(locBCID)
          END IF ! Do2DSurfaceCharge
#endif /*USE_HDG*/
        END DO ! iPart = 1, NbrOfSEE
      END IF ! NbrOfSEE.GT.0
    END DO ! q = 1, Ray%nSurfSample
  END DO ! p = 1, Ray%nSurfSample
#if USE_LOADBALANCE
CALL LBElemSplitTime(locElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iSurfSide = 1, nComputeNodeSurfSides

END SUBROUTINE PhotoIonization_RayTracing_SEE


SUBROUTINE PhotoIonization_RayTracing_Volume()
!===================================================================================================================================
!> Routine calculates the number of photo-ionization reactions, utilizing the cell-local photon energy from the raytracing
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
! Variables
USE MOD_Globals_Vars            ,ONLY: PI, c
USE MOD_Timedisc_Vars           ,ONLY: dt
USE MOD_Mesh_Vars               ,ONLY: nElems, offsetElem
USE MOD_Mesh_Vars               ,ONLY: NGeo,wBaryCL_NGeo,XiCL_NGeo,XCL_NGeo
USE MOD_RayTracing_Vars         ,ONLY: UseRayTracing, Ray,RayElemEmission
USE MOD_RayTracing_Vars         ,ONLY: U_N_Ray_loc,N_DG_Ray_loc,N_Inter_Ray
USE MOD_RayTracing_Vars         ,ONLY: RaySecondaryVectorX,RaySecondaryVectorY,RaySecondaryVectorZ
USE MOD_Particle_Vars           ,ONLY: Species, PartState, usevMPF, PartMPF, PDM, PEM, PartSpecies
USE MOD_DSMC_Vars               ,ONLY: ChemReac, DSMC, BGGas, Coll_pData, CollisMode
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
! Functions/Subroutines
USE MOD_Eval_xyz                ,ONLY: TensorProductInterpolation
USE MOD_part_emission_tools     ,ONLY: CalcPhotonEnergy
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_DSMC_PolyAtomicModel    ,ONLY: DSMC_SetInternalEnr
USE MOD_part_tools              ,ONLY: CalcVelocity_maxwell_particle
#if defined(LSERK)
USE MOD_TimeDisc_Vars           ,ONLY: nRKStages,RK_c
#endif /*defined(LSERK)*/
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Boundary_Vars  ,ONLY: DoBoundaryParticleOutputRay
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers      ,ONLY: LBStartTime, LBElemSplitTime
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Tools     ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iElem,k,l,m,iReac,iPair,iGlobalElem,iVar
INTEGER               :: SpecID,nPair,NRayLoc,BGGSpecID
INTEGER               :: PartID,newPartID
REAL                  :: E_Intensity, TimeScalingFactor, PhotonEnergy
REAL                  :: density, NbrOfPhotons, NbrOfReactions
REAL                  :: RandNum,RandVal(3),Xi(3)
REAL                  :: RandomPos(1:3),MPF
REAL                  :: RayDirection(1:3),RayBaseVector1IC(1:3),RayBaseVector2IC(1:3)
#if USE_LOADBALANCE
REAL                  :: tLBStart
#endif /*USE_LOADBALANCE*/
LOGICAL               :: positionIsInside
!===================================================================================================================================

IF(.NOT.UseRayTracing) RETURN

! Photo-ionization is currently only supported with CollisMode = 3
IF(CollisMode.NE.3) RETURN

IF(.NOT.ChemReac%AnyPhIonReaction) RETURN

! Get the time-dependent pulse intensity factor based on the pulse type (square, Gaussian, etc.)
CALL GetPulseIntensity(TimeScalingFactor)

! Get the photon energy (currently: single, fixed wave length -> calculate energy once)
PhotonEnergy = CalcPhotonEnergy(Ray%WaveLength)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/

! Loop over the first and secondary rays
DO iVar = 1, 2
  ! Set the direction for the first ray direction before the element loop since it does not change
  IF(iVar.EQ.1) THEN
    RayDirection(1:3) = Ray%Direction(1:3)
    RayBaseVector1IC(1:3) = Ray%BaseVector1IC(1:3)
    RayBaseVector2IC(1:3) = Ray%BaseVector2IC(1:3)
  END IF
  ! Loop over all elements
  DO iElem=1, nElems
    ! Skip elements without ionization
    IF(.NOT.RayElemEmission(iVar,iElem)) CYCLE
    iGlobalElem = iElem+offSetElem
    ! Set the ray-direction in case of secondary rays
    IF(iVar.EQ.2) THEN
      ! Generate two base vectors perpendicular to the ray direction
      RayDirection(1:3) = (/RaySecondaryVectorX(iElem),RaySecondaryVectorY(iElem),RaySecondaryVectorZ(iElem)/)
      ! Check whether the direction is zero
      IF(.NOT.ALL(RayDirection(1:3).EQ.0.)) THEN
        CALL OrthoNormVec(RayDirection(1:3),RayBaseVector1IC(1:3),RayBaseVector2IC(1:3))
      ELSE
        CYCLE
      END IF
    END IF
    ! Set the element local order and loop over all sub-elements
    NRayLoc = N_DG_Ray_loc(iElem)
    DO m=0,NRayLoc
      DO l=0,NRayLoc
        DO k=0,NRayLoc
          E_Intensity = U_N_Ray_loc(iElem)%U(iVar,k,l,m) * TimeScalingFactor
          ! Number of photons (spectrum not implemented)
          NbrOfPhotons = E_Intensity / (PhotonEnergy * c * dt)
          DO iReac = 1, ChemReac%NumOfReact
            SpecID = ChemReac%Reactants(iReac,1)
            ! FEATURE: Background gas density distribution
            BGGSpecID = BGGas%MapSpecToBGSpec(SpecID)
            density = BGGas%NumberDensity(BGGSpecID)
            ! Determine the number of particles to insert
            ! Collision number: Z = n_gas * n_ph * sigma_reac * v (in the case of photons its speed of light)
            ! Number of reactions: N = Z * dt * V (number of photons cancels out the volume)
            ! Number of reactions: N = n_gas * N_ph * sigma_reac * v * dt
            NbrOfReactions = density * NbrOfPhotons * ChemReac%CrossSection(iReac) * c * dt / Species(SpecID)%MacroParticleFactor
            CALL RANDOM_NUMBER(RandNum)
            nPair = INT(NbrOfReactions+RandNum)
            ALLOCATE(Coll_pData(nPair))
            Coll_pData%Ec = 0.
            ! Loop over all newly created particles
            DO iPair = 1, nPair
              ! When using Triatracking, the reference position may not be inside the element.
              ! Create random numbers until the position is inside the element
              positionIsInside = .FALSE.
              DO WHILE(.NOT.positionIsInside)
                ! Get a random position in the subelement
                CALL RANDOM_NUMBER(RandVal)
                Xi(1) = -1.0 + SUM(N_Inter_Ray(NRayLoc)%wGP(0:k-1)) + N_Inter_Ray(NRayLoc)%wGP(k) * RandVal(1)
                Xi(2) = -1.0 + SUM(N_Inter_Ray(NRayLoc)%wGP(0:l-1)) + N_Inter_Ray(NRayLoc)%wGP(l) * RandVal(2)
                Xi(3) = -1.0 + SUM(N_Inter_Ray(NRayLoc)%wGP(0:m-1)) + N_Inter_Ray(NRayLoc)%wGP(m) * RandVal(3)
                IF(ANY(Xi.GT.1.0).OR.ANY(Xi.LT.-1.0))THEN
                  IPWRITE(UNIT_StdOut,*) "Xi =", Xi
                  CALL abort(__STAMP__,'xi out of range')
                END IF ! ANY(Xi.GT.1.0).OR.ANY(Xi.LT.-1.0)
                ! Get the physical coordinates that correspond to the reference coordinates
                CALL TensorProductInterpolation(Xi(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),RandomPos(1:3))
                IF(TrackingMethod.EQ.TRIATRACKING) THEN
                  CALL ParticleInsideQuad3D(RandomPos,iGlobalElem,positionIsInside)
                ELSE
                  positionIsInside = .TRUE.
                END IF
              END DO

              ! Create new particle from the background gas
              PartID = GetNextFreePosition()
              ! Set the position
              PartState(1:3,PartID) = RandomPos(1:3)
              ! Set the species
              PartSpecies(PartID) = SpecID
              ! Set the velocity (required for the collision energy, although relatively small compared to the photon energy)
              IF(BGGas%UseDistribution) THEN
                PartState(4:6,PartID) = CalcVelocity_maxwell_particle(SpecID,BGGas%Distribution(BGGSpecID,4:6,iElem)) &
                                              + BGGas%Distribution(BGGSpecID,1:3,iElem)
              ELSE
                CALL CalcVelocity_maxwell_lpn(FractNbr=SpecID, Vec3D=PartState(4:6,PartID), iInit=1)
              END IF
              ! Ambipolar diffusion
              IF (DSMC%DoAmbipolarDiff) THEN
                newAmbiParts = newAmbiParts + 1
                iPartIndx_NodeNewAmbi(newAmbiParts) = PartID
              END IF
              ! Set the internal energies
              IF(CollisMode.GT.1) THEN
                CALL DSMC_SetInternalEnr(SpecID,1,PartID,1)
              END IF
              ! Particle flags
              PDM%ParticleInside(PartID)  = .TRUE.
              PDM%IsNewPart(PartID)       = .TRUE.
              PDM%dtFracPush(PartID)      = .FALSE.
              PEM%GlobalElemID(PartID)     = iGlobalElem
              PEM%LastGlobalElemID(PartID) = 0 ! Initialize with invalid value
              ! Create second particle (only the index and the flags/elements needs to be set)
              newPartID = GetNextFreePosition()
              IF (DSMC%DoAmbipolarDiff) THEN
                newAmbiParts = newAmbiParts + 1
                iPartIndx_NodeNewAmbi(newAmbiParts) = newPartID
              END IF
              ! Set the position
              PartState(1:3,newPartID) = RandomPos(1:3)
              ! Particle flags
              PDM%ParticleInside(newPartID)  = .TRUE.
              PDM%IsNewPart(newPartID)       = .TRUE.
              PDM%dtFracPush(newPartID)      = .FALSE.
              PEM%GlobalElemID(newPartID)     = iGlobalElem
              PEM%LastGlobalElemID(newPartID) = 0 ! Initialize with invalid value
              ! Pairing (first particle is the background gas species)
              Coll_pData(iPair)%iPart_p1 = PartID
              Coll_pData(iPair)%iPart_p2 = newPartID
              ! Relative velocity is not required as the relative translational energy will not be considered
              Coll_pData(iPair)%CRela2 = 0.
              ! Weighting factor: use the weighting factor of the emission init
              MPF = Species(SpecID)%MacroParticleFactor
              IF(usevMPF) THEN
                PartMPF(PartID)    = MPF
                PartMPF(newPartID) = MPF
              END IF
              ! Store the ion particle information in PartStateBoundary.h5
              IF(DoBoundaryParticleOutputRay) THEN
                CALL StoreBoundaryParticleProperties(PartID,SpecID,PartState(1:3,PartID),&
                    UNITVECTOR(PartState(4:6,PartID)),(/0.,0.,1./),iPartBound=0,mode=2,MPF_optIN=MPF)
                CALL StoreBoundaryParticleProperties(NewPartID,SpecID,PartState(1:3,NewPartID),&
                    UNITVECTOR(PartState(4:6,NewPartID)),(/0.,0.,1./),iPartBound=0,mode=2,MPF_optIN=MPF)
              END IF ! DoBoundaryParticleOutputRay
              ! Velocity (set it to zero, as it will be subtracted in the chemistry module)
              PartState(4:6,newPartID) = 0.
              ! Insert the products and distribute the reaction energy (Requires: Pair indices, Coll_pData(iPair)%iPart_p1/2)
              IF(DoBoundaryParticleOutputRay) THEN
                ! Store the electron particle information in PartStateBoundary.h5
                CALL PhotoIonization_InsertProducts(iPair, iReac, RayBaseVector1IC, RayBaseVector2IC, RayDirection, PartBCIndex=0)
              ELSE
                ! DO NOT store the electron particle information in PartStateBoundary.h5
                CALL PhotoIonization_InsertProducts(iPair, iReac, RayBaseVector1IC, RayBaseVector2IC, RayDirection, PartBCIndex=-1)
              END IF
            END DO  ! iPart = 1, nPair
            DEALLOCATE(Coll_pData)
          END DO    ! iReac = 1, ChemReac%NumOfReact
        END DO      ! k
      END DO        ! l
    END DO          ! m
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO            ! iElem = 1, nElems
END DO              ! iVar = 1, 2

END SUBROUTINE PhotoIonization_RayTracing_Volume


!===================================================================================================================================
!> Calculate the number of photo-ionization reactions to perform based on a given number of photons and cross-section
!===================================================================================================================================
SUBROUTINE CalcPhotoIonizationNumber(iSpec,NbrOfPhotons,NbrOfReactions)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: c
USE MOD_Particle_Vars ,ONLY: Species
USE MOD_DSMC_Vars     ,ONLY: BGGas,ChemReac
USE MOD_MCC_Vars      ,ONLY: NbrOfPhotonXsecReactions,SpecPhotonXSecInterpolated
USE MOD_MCC_Vars      ,ONLY: PhotoIonFirstLine,PhotoIonLastLine,PhotonDistribution,PhotoReacToReac
USE MOD_TimeDisc_Vars ,ONLY: dt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec              !< Species index
REAL, INTENT(IN)              :: NbrOfPhotons       !< Simulated number of photons
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: NbrOfReactions     !< Number of reactions to be performed
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iReac,iPhotoReac,iLine
REAL                          :: density
!===================================================================================================================================

NbrOfReactions = 0.

IF(NbrOfPhotonXsecReactions.GT.0)THEN
  ! Distribute the photons according to the distribution function
  PhotonDistribution = SpecPhotonXSecInterpolated(:,2) * NbrOfPhotons

  DO iPhotoReac = 1, NbrOfPhotonXsecReactions
    iReac          = PhotoReacToReac(iPhotoReac)
    density        = BGGas%NumberDensity(BGGas%MapSpecToBGSpec(ChemReac%Reactants(iReac,1)))
    ! Only consider lines with cross-section larger than zero
    DO iLine = PhotoIonFirstLine, PhotoIonLastLine
      ASSOCIATE( CrossSection => SpecPhotonXSecInterpolated(iLine,2+iPhotoReac) )
        ! Consider the ratio of the cross-section to the sum of al cross-sections
        NbrOfReactions = NbrOfReactions + PhotonDistribution(iLine) * CrossSection
      END ASSOCIATE
    END DO ! PhotoIonFirstLine, PhotoIonLastLine
  END DO ! iPhotoReac = 1, NbrOfPhotonXsecReactions
  NbrOfReactions = NbrOfReactions * density * c * dt / Species(iSpec)%MacroParticleFactor
END IF ! NbrOfPhotonXsecReactions.GT.0

! Photoionization reactions with constant cross sections
DO iReac = 1, ChemReac%NumOfReact
  ! Only treat photoionization reactions
  IF(TRIM(ChemReac%ReactModel(iReac)).NE.'phIon') CYCLE
  ! First reactant of the reaction is the actual heavy particle species
  ASSOCIATE( density      => BGGas%NumberDensity(BGGas%MapSpecToBGSpec(ChemReac%Reactants(iReac,1))) ,&
             CrossSection => ChemReac%CrossSection(iReac))
    ! Collision number: Z = n_gas * n_ph * sigma_reac * v (in the case of photons its speed of light)
    ! Number of reactions: N = Z * dt * V (number of photons cancels out the volume)
    NbrOfReactions = NbrOfReactions + density * NbrOfPhotons * CrossSection * c * dt / Species(iSpec)%MacroParticleFactor
  END ASSOCIATE
END DO

END SUBROUTINE CalcPhotoIonizationNumber


!===================================================================================================================================
!> Routine performing the photo-ionization reaction: initializing the heavy species at the background gas temperature (first
!> reactant) and distributing the remaining collision energy onto the electrons
!===================================================================================================================================
SUBROUTINE PhotoIonization_InsertProducts(iPair, iReac, b1, b2, normal, iLineOpt, PartBCIndex)
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: eV2Joule
USE MOD_DSMC_Vars               ,ONLY: Coll_pData, DSMC, SpecDSMC, CollInf
USE MOD_DSMC_Vars               ,ONLY: ChemReac,PartIntEn
USE MOD_DSMC_Vars               ,ONLY: newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_MCC_Vars                ,ONLY: ReacToPhotoReac,NbrOfPhotonXsecReactions,SpecPhotonXSecInterpolated
USE MOD_Particle_Vars           ,ONLY: PartSpecies, PartState, PDM, PEM, PartPosRef, Species, PartMPF, usevMPF
USE MOD_Particle_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_Tracking_Vars  ,ONLY: TrackingMethod
USE MOD_part_tools              ,ONLY: GetParticleWeight, DiceUnitVector,CalcEVib_particle, RotInitPolyRoutineFuncPTR
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition, CalcEElec_particle
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_Particle_Analyze_Vars   ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Analyze_Pure   ,ONLY: CalcEkinPart
USE MOD_Particle_Boundary_Vars  ,ONLY: DoBoundaryParticleOutputHDF5
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_DSMC_CollisVec          ,ONLY: PostCollVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iPair, iReac
REAL, INTENT(IN), OPTIONAL    :: b1(3),b2(3),normal(3)
INTEGER, INTENT(IN), OPTIONAL :: iLineOpt
INTEGER, INTENT(IN), OPTIONAL :: PartBCIndex
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iPart, iSpec, iProd, NumProd
INTEGER                       :: ReactInx(1:4), EductReac(1:3), ProductReac(1:4)
REAL                          :: Weight(1:4), SumWeightProd, Mass_Electron, CRela2_Electron, RandVal, NumElec
REAL                          :: VeloCOM(1:3), Temp_Trans, Temp_Rot, Temp_Vib, Temp_Elec,EForm
REAL                          :: FracMassCent1, FracMassCent2, MassRed, cRelaNew(1:3),MPF
LOGICAL                       :: IonizationReaction
!===================================================================================================================================

EductReac(1:3) = ChemReac%Reactants(iReac,1:3)
ProductReac(1:4) = ChemReac%Products(iReac,1:4)

Weight = 0.
NumProd = 2; SumWeightProd = 0.

! Reaction is only defined with a single reactant (that of the background gas species)
IF (PartSpecies(Coll_pData(iPair)%iPart_p1).EQ.ChemReac%Reactants(iReac,1)) THEN
  ReactInx(1) = Coll_pData(iPair)%iPart_p1
  ReactInx(2) = Coll_pData(iPair)%iPart_p2
ELSE IF (PartSpecies(Coll_pData(iPair)%iPart_p2).EQ.ChemReac%Reactants(iReac,1)) THEN
  ReactInx(2) = Coll_pData(iPair)%iPart_p1
  ReactInx(1) = Coll_pData(iPair)%iPart_p2
ELSE
  CALL abort(__STAMP__,'ERROR in PhotoIonization_InsertProducts: Pair does not correspond to the reactants!')
END IF

! Do not perform the reaction in case the reaction is to be calculated at a constant gas composition (DSMC%ReservoirSimuRate = T)
IF (DSMC%ReservoirSimu.AND.DSMC%ReservoirSimuRate) THEN
  ! Count the number of reactions to determine the actual reaction rate
  IF (DSMC%ReservoirRateStatistic) ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + GetParticleWeight(ReactInx(1))
  ! Leave the routine again
  RETURN
END IF

! Set the particle weights to the same as the background species
Weight(1:2) = GetParticleWeight(ReactInx(1))
! Set the particle species for the first two products (the other products require an index first)
DO iProd = 1, NumProd
  PartSpecies(ReactInx(iProd)) = ProductReac(iProd)
END DO

! Calculate the sum of the weights of the products
SumWeightProd = Weight(1) + Weight(2)

IF(EductReac(3).EQ.0) THEN
  IF(ProductReac(3).NE.0) THEN
    ! === Get free particle index for the 3rd product
    ReactInx(3) = GetNextFreePosition()
    PDM%ParticleInside(ReactInx(3)) = .true.
    PDM%IsNewPart(ReactInx(3)) = .true.
    PDM%dtFracPush(ReactInx(3)) = .FALSE.
    ! Set species index of new particle
    PartSpecies(ReactInx(3)) = ProductReac(3)
    PartState(1:3,ReactInx(3)) = PartState(1:3,ReactInx(1))
    IF(TrackingMethod.EQ.REFMAPPING) THEN
      PartPosRef(1:3,ReactInx(3))=PartPosRef(1:3,ReactInx(1))
    END IF
    IF((Species(ProductReac(3))%InterID.EQ.2).OR.(Species(ProductReac(3))%InterID.EQ.20)) THEN
      ALLOCATE(PartIntEn(ReactInx(3))%ERot(1), PartIntEn(ReactInx(3))%EVib(1))
      PartIntEn(ReactInx(3))%EVib(1)= 0.
      PartIntEn(ReactInx(3))%ERot(1) = 0.
    END IF
    IF(DSMC%ElectronicModel.GT.0) THEN
      IF((Species(ProductReac(3))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(3))%FullyIonized)) THEN
        ALLOCATE(PartIntEn(ReactInx(3))%EElec(1))
        PartIntEn(ReactInx(3))%EElec(1) = 0.
      END IF
    END IF
    PEM%GlobalElemID(ReactInx(3)) = PEM%GlobalElemID(ReactInx(1))
    PEM%LastGlobalElemID(ReactInx(3)) = PEM%GlobalElemID(ReactInx(3))
    IF(usevMPF) PartMPF(ReactInx(3)) = PartMPF(ReactInx(1))
    IF(UseVarTimeStep) PartTimeStep(ReactInx(3)) = PartTimeStep(ReactInx(1))
    Weight(3) = Weight(1)
    NumProd = 3
    SumWeightProd = SumWeightProd + Weight(3)
    IF (DSMC%DoAmbipolarDiff) THEN
      newAmbiParts = newAmbiParts + 1
      iPartIndx_NodeNewAmbi(newAmbiParts) = ReactInx(3)
    END IF
  END IF
END IF

IF(ProductReac(4).NE.0) THEN
  ! === Get free particle index for the 4th product
  ReactInx(4) = GetNextFreePosition()
  PDM%ParticleInside(ReactInx(4)) = .true.
  PDM%IsNewPart(ReactInx(4)) = .true.
  PDM%dtFracPush(ReactInx(4)) = .FALSE.
  ! Set species index of new particle
  PartSpecies(ReactInx(4)) = ProductReac(4)
  PartState(1:3,ReactInx(4)) = PartState(1:3,ReactInx(1))
  IF(TrackingMethod.EQ.REFMAPPING) THEN ! here Nearst-GP is missing
    PartPosRef(1:3,ReactInx(4))=PartPosRef(1:3,ReactInx(1))
  END IF
  IF((Species(ProductReac(4))%InterID.EQ.2).OR.(Species(ProductReac(4))%InterID.EQ.20)) THEN
    ALLOCATE(PartIntEn(ReactInx(4))%ERot(1), PartIntEn(ReactInx(4))%EVib(1))
    PartIntEn(ReactInx(4))%EVib(1)= 0.
    PartIntEn(ReactInx(4))%ERot(1) = 0.
  END IF
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(ProductReac(4))%InterID.NE.4).AND.(.NOT.SpecDSMC(ProductReac(4))%FullyIonized)) THEN
      ALLOCATE(PartIntEn(ReactInx(4))%EElec(1))
      PartIntEn(ReactInx(4))%EElec(1) = 0.
    END IF
  END IF
  PEM%GlobalElemID(ReactInx(4)) = PEM%GlobalElemID(ReactInx(1))
  PEM%LastGlobalElemID(ReactInx(4)) = PEM%GlobalElemID(ReactInx(4))
  IF(usevMPF) PartMPF(ReactInx(4)) = PartMPF(ReactInx(1))
  IF(UseVarTimeStep) PartTimeStep(ReactInx(4)) = PartTimeStep(ReactInx(1))
  Weight(4) = Weight(1)
  NumProd = 4
  SumWeightProd = SumWeightProd + Weight(4)
  IF (DSMC%DoAmbipolarDiff) THEN
    newAmbiParts = newAmbiParts + 1
    iPartIndx_NodeNewAmbi(newAmbiParts) = ReactInx(4)
  END IF
END IF

! Set the formation energy
IF(NbrOfPhotonXsecReactions.GT.0)THEN
  IF(SpecPhotonXSecInterpolated(iLineOpt,2+ReacToPhotoReac(iReac)).LE.0)THEN
    IPWRITE(UNIT_StdOut,'(I6,3(A,I3))') " (sigma=0) iLine =",iLineOpt," iPhotoReac =",ReacToPhotoReac(iReac)," iReac =",iReac
    CALL abort(__STAMP__,'Cross-section for this reaction is zero')
  END IF
  ! Add photon energy to the formation energy of the reaction
  EForm = ChemReac%EForm(iReac) + SpecPhotonXSecInterpolated(iLineOpt,1)*eV2Joule
  IF(EForm.LE.0)THEN
    IPWRITE(UNIT_StdOut,'(I6,3(A,I3))') " (EForm=0) iLine =",iLineOpt," iPhotoReac =",ReacToPhotoReac(iReac)," iReac =",iReac
    IPWRITE(UNIT_StdOut,*) "Photon energy [J] =", SpecPhotonXSecInterpolated(iLineOpt,1)*eV2Joule ,&
        "and in [eV]",SpecPhotonXSecInterpolated(iLineOpt,1)
    CALL abort(__STAMP__,'Energy of formation for photoionization is zero or negative: ',RealInfoOpt=EForm)
  END IF
ELSE
  EForm = ChemReac%EForm(iReac)
END IF ! NbrOfPhotonXsecReactions.GT.0

! Consider the energy of the background gas particle and the remaining energy from the photo-reaction
Coll_pData(iPair)%Ec = 0.5 * Species(PartSpecies(ReactInx(1)))%MassIC * DOTPRODUCT(PartState(4:6,ReactInx(1))) * Weight(1) &
                      + EForm*SumWeightProd/NumProd
! Adding the vibrational and rotational energy to the collision energy
IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) &
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + (PartIntEn(ReactInx(1))%EVib(1) + PartIntEn(ReactInx(1))%ERot(1))*Weight(1)
! Addition of the electronic energy to the collision energy
IF (DSMC%ElectronicModel.GT.0) THEN
  IF((Species(EductReac(1))%InterID.NE.4).AND.(.NOT.SpecDSMC(EductReac(1))%FullyIonized)) &
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartIntEn(ReactInx(1))%EElec(1)*Weight(1)
END IF

! Saving the velocity of the background particle as the centre of mass velocity
VeloCOM(1:3) = PartState(4:6,ReactInx(1))
! Get the properties of the background species used for the photo-ionization reaction
Temp_Trans = Species(EductReac(1))%Init(1)%MWTemperatureIC
IF((Species(EductReac(1))%InterID.EQ.2).OR.(Species(EductReac(1))%InterID.EQ.20)) THEN
  Temp_Vib   = SpecDSMC(EductReac(1))%Init(1)%TVib
  Temp_Rot   = SpecDSMC(EductReac(1))%Init(1)%TRot
ELSE
  Temp_Vib   = Temp_Trans
  Temp_Rot   = Temp_Trans
END IF
IF(DSMC%ElectronicModel.GT.0) Temp_Elec = SpecDSMC(EductReac(1))%Init(1)%TElec
!-------------------------------------------------------------------------------------------------------------------------------
! Insert the heavy species at the properties of the background gas
!-------------------------------------------------------------------------------------------------------------------------------
NumElec = 0.
IonizationReaction = .FALSE.
DO iProd = 1, NumProd
  iPart = ReactInx(iProd)
  iSpec = ProductReac(iProd)
  IF(Species(iSpec)%InterID.EQ.4) THEN
    NumElec = NumElec + Weight(iProd)
    Mass_Electron = Species(iSpec)%MassIC
    IonizationReaction = .TRUE.
    CYCLE
  END IF
  ! Set the internal energies
  IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) THEN
    IF (.NOT.ALLOCATED(PartIntEn(iPart)%EVib)) ALLOCATE(PartIntEn(iPart)%EVib(1))
    IF (.NOT.ALLOCATED(PartIntEn(iPart)%ERot)) ALLOCATE(PartIntEn(iPart)%ERot(1))
    PartIntEn(iPart)%EVib = CalcEVib_particle(iSpec,Temp_Vib,iPart)
    PartIntEn(iPart)%ERot = RotInitPolyRoutineFuncPTR(iSpec,Temp_Rot,iPart)
  ELSE
    SDEALLOCATE(PartIntEn(iPart)%EVib)
    SDEALLOCATE(PartIntEn(iPart)%ERot)
  END IF
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) THEN
      IF (.NOT.ALLOCATED(PartIntEn(iPart)%EElec)) ALLOCATE(PartIntEn(iPart)%EElec(1))
      PartIntEn(iPart)%EElec = CalcEElec_particle(iSpec,Temp_Elec,iPart)
    ELSE
      SDEALLOCATE(PartIntEn(iPart)%EElec)
    END IF
  END IF
  ! Determine the particle velocity (is going to be added to the PartState)
  CALL CalcVelocity_maxwell_lpn(FractNbr=iSpec, Vec3D=PartState(4:6,iPart), Temperature=Temp_Trans)
  ! Remove the distributed energy from the available collision energy
  Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - 0.5 * Species(iSpec)%MassIC * DOTPRODUCT(PartState(4:6,iPart)) * Weight(iProd)
  IF((Species(iSpec)%InterID.EQ.2).OR.(Species(iSpec)%InterID.EQ.20)) &
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - (PartIntEn(iPart)%EVib(1) + PartIntEn(iPart)%ERot(1)) * Weight(iProd)
  IF(DSMC%ElectronicModel.GT.0) THEN
    IF((Species(iSpec)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpec)%FullyIonized)) &
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartIntEn(iPart)%EElec(1) * Weight(iProd)
  END IF
  IF(Coll_pData(iPair)%Ec.LE.0)THEN
    IF(NbrOfPhotonXsecReactions.GT.0)THEN
      IPWRITE(UNIT_StdOut,'(I6,3(A,I3))') " (%Ec=0)   iLine =",iLineOpt," iPhotoReac =",ReacToPhotoReac(iReac)," iReac =",iReac
      IPWRITE(UNIT_StdOut,'(I6,A,E24.12,A,E15.2)') " Photon energy [J] =", SpecPhotonXSecInterpolated(iLineOpt,1)*eV2Joule ,&
          "and in [eV]",SpecPhotonXSecInterpolated(iLineOpt,1)
    ELSE
      IPWRITE(UNIT_StdOut,*) "iReac =", iReac
    END IF ! NbrOfPhotonXsecReactions.GT.0
    CALL abort(__STAMP__,'Coll_pData(iPair)%Ec is zero or negative: ',RealInfoOpt=Coll_pData(iPair)%Ec)
  END IF
END DO
!--------------------------------------------------------------------------------------------------!
! Calculation of new electron velocities OR distribute remaining energy onto the heavy species (currently only for 2 products)
!--------------------------------------------------------------------------------------------------!
IF(IonizationReaction) THEN
  CRela2_Electron = NumElec * Mass_Electron
  CRela2_Electron = 2. * Coll_pData(iPair)%Ec / CRela2_Electron
  DO iProd = 1, NumProd
    iPart = ReactInx(iProd)
    iSpec = ProductReac(iProd)
    ! Check if particle is an electron
    IF(Species(iSpec)%InterID.EQ.4) THEN
      PartState(4:6,iPart) = VeloCOM(1:3) + SQRT(CRela2_Electron) * DiceUnitVector()
      ! Change the direction of its velocity vector (randomly) to be perpendicular to the photon's path
      ! Get random vector b3 in b1-b2-plane
      CALL RANDOM_NUMBER(RandVal)
      PartState(4:6,iPart) = GetRandomVectorInPlane(b1,b2,PartState(4:6,iPart),RandVal)
      ! Rotate the resulting vector in the b3-NormalIC-plane
      PartState(4:6,iPart) = GetRotatedVector(PartState(4:6,iPart),normal)
      ! Store the particle information in PartStateBoundary.h5
      IF(DoBoundaryParticleOutputHDF5) THEN
        IF(usevMPF)THEN
          MPF = PartMPF(iPart) ! Use emission-specific MPF
        ELSE
          MPF = Species(iSpec)%MacroParticleFactor ! Use species MPF
        END IF ! usevMPF
        ! Only store volume-emitted particle data in PartStateBoundary.h5 if the PartBCIndex is greater/equal zero
        IF(PartBCIndex.GE.0) CALL StoreBoundaryParticleProperties(iPart,iSpec,PartState(1:3,iPart),&
                                    UNITVECTOR(PartState(4:6,iPart)),normal,iPartBound=PartBCIndex,mode=2,MPF_optIN=MPF)
      END IF ! DoBoundaryParticleOutputHDF5
    END IF
  END DO
ELSE
  IF (usevMPF) THEN
    FracMassCent1 = Species(ProductReac(1))%MassIC *Weight(1) &
        /(Species(ProductReac(1))%MassIC* Weight(1) + Species(ProductReac(2))%MassIC * Weight(2))
    FracMassCent2 = Species(ProductReac(2))%MassIC *Weight(2) &
        /(Species(ProductReac(1))%MassIC* Weight(1) + Species(ProductReac(2))%MassIC * Weight(2))
    MassRed = Species(ProductReac(1))%MassIC *Weight(1)* Species(ProductReac(2))%MassIC *Weight(2) &
        / (Species(ProductReac(1))%MassIC*Weight(1) + Species(ProductReac(2))%MassIC *Weight(2))
  ELSE
    ! Scattering of (AB)
    FracMassCent1 = CollInf%FracMassCent(ProductReac(1),CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
    FracMassCent2 = CollInf%FracMassCent(ProductReac(2),CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
    MassRed = CollInf%MassRed(CollInf%Coll_Case(ProductReac(1),ProductReac(2)))
  END IF

  Coll_pData(iPair)%cRela2 = 2 * Coll_pData(iPair)%Ec / MassRed
  cRelaNew(1:3) = PostCollVec(iPair)

  !deltaV particle 1
  PartState(4:6,ReactInx(1)) = VeloCOM(1:3) + FracMassCent2*cRelaNew(1:3)
  !deltaV particle 2
  PartState(4:6,ReactInx(2)) = VeloCOM(1:3) - FracMassCent1*cRelaNew(1:3)
END IF

IF(CalcPartBalance) THEN
  DO iProd = 1, NumProd
    iSpec = ProductReac(iProd)
    nPartIn(iSpec) = nPartIn(iSpec) + 1
    PartEkinIn(iSpec) = PartEkinIn(iSpec) + CalcEkinPart(ReactInx(iProd))
  END DO
END IF

END SUBROUTINE PhotoIonization_InsertProducts


!===================================================================================================================================
!> Pick random vector in a plane set up by the basis vectors b1 and b2
!===================================================================================================================================
PPURE FUNCTION GetRandomVectorInPlane(b1,b2,VeloVec,RandVal)
! MODULES
USE MOD_Globals      ,ONLY: VECNORM3D
USE MOD_Globals_Vars ,ONLY: PI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: b1(1:3),b2(1:3) !< Basis vectors (normalized)
REAL,INTENT(IN)    :: VeloVec(1:3)    !< Velocity vector before the random direction selection within the plane defined by b1 and b2
REAL,INTENT(IN)    :: RandVal         !< Random number (given from outside to render this function PPURE)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL               :: GetRandomVectorInPlane(1:3) ! Output velocity vector
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vabs ! Absolute velocity
REAL               :: phi ! random angle between 0 and 2*PI
!===================================================================================================================================
Vabs = VECNORM3D(VeloVec)
phi = RandVal * 2.0 * PI
GetRandomVectorInPlane = Vabs*(b1*COS(phi) + b2*SIN(phi))
END FUNCTION GetRandomVectorInPlane


!===================================================================================================================================
!> Rotate the vector in the plane set up by VeloVec and NormVec by choosing an angle from a 4.0 / PI * COS(Theta)**2
!> distribution via the ARM
!===================================================================================================================================
FUNCTION GetRotatedVector(VeloVec,NormVec)
! MODULES
USE MOD_Globals      ,ONLY: VECNORM3D, UNITVECTOR
USE MOD_Globals_Vars ,ONLY: PI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: NormVec(1:3) !< Basis vector (normalized)
REAL,INTENT(IN)    :: VeloVec(1:3) !< Velocity vector before the random direction selection within the plane defined by b1 and b2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL               :: GetRotatedVector(1:3) !< Output velocity vector
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vabs ! Absolute velocity
REAL               :: RandVal, v(1:3)
REAL               :: Theta, Theta_temp
REAL               :: PDF_temp
REAL, PARAMETER    :: PDF_max=4./PI
LOGICAL            :: ARM_SEE_PDF
!===================================================================================================================================
v = UNITVECTOR(VeloVec)
Vabs = VECNORM3D(VeloVec)

! ARM for angular distribution
ARM_SEE_PDF=.TRUE.
DO WHILE(ARM_SEE_PDF)
  CALL RANDOM_NUMBER(RandVal)
  Theta_temp = PI*(RandVal-0.5)
  PDF_temp = 4.0 / PI * COS(Theta_temp)**2
  CALL RANDOM_NUMBER(RandVal)
  IF ((PDF_temp/PDF_max).GT.RandVal) ARM_SEE_PDF = .FALSE.
END DO
Theta = Theta_temp

! Rotate original vector Vabs*v
GetRotatedVector = Vabs*(v*COS(Theta) + NormVec*SIN(Theta))
END FUNCTION GetRotatedVector


END MODULE MOD_Particle_Photoionization