!==================================================================================================================================
! Copyright (c) 2023 boltzplatz - numerical plasma dynamics GmbH
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
PUBLIC :: PhotoIonization_RayTracing_SEE
!===================================================================================================================================
CONTAINS

SUBROUTINE PhotoIonization_RayTracing_SEE()
!===================================================================================================================================
!> Routine calculates the number of secondary electrons to be emitted and inserts them on the surface, utilizing the cell-local
!> photon energy from the raytracing
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI
USE MOD_Timedisc_Vars           ,ONLY: dt,time
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample, Partbound, SurfSide2GlobalSide, DoBoundaryParticleOutputHDF5
USE MOD_Particle_Vars           ,ONLY: Species, PartState, usevMPF
USE MOD_RayTracing_Vars         ,ONLY: Ray,RayPartBound
USE MOD_part_emission_tools     ,ONLY: CalcPhotonEnergy
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Surfaces_Vars  ,ONLY: BezierControlPoints3D, BezierSampleXi
USE MOD_Particle_Surfaces       ,ONLY: EvaluateBezierPolynomialAndGradient, CalcNormAndTangBezier
USE MOD_Mesh_Vars               ,ONLY: NGeo
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_FromWorkFuncSEE
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_part_operations         ,ONLY: CreateParticle
#ifdef LSERK
USE MOD_Timedisc_Vars           ,ONLY: iStage, RK_c, nRKStages
#endif
#if USE_MPI
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfSideArea_Shared,nComputeNodeSurfTotalSides
USE MOD_Photon_TrackingVars     ,ONLY: PhotonSampWall_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: nComputeNodeProcessors,myComputeNodeRank
#else
USE MOD_Photon_TrackingVars     ,ONLY: PhotonSampWall
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfTotalSides
#endif /*USE_MPI*/
#if USE_HDG
USE MOD_HDG_Vars                ,ONLY: UseFPC,FPC,UseEPC,EPC
USE MOD_Mesh_Vars               ,ONLY: BoundaryType
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: t_1, t_2, E_Intensity
INTEGER               :: NbrOfRepetitions, firstSide, lastSide, SideID, iSample, GlobElemID, PartID
INTEGER               :: iSurfSide, p, q, BCID, SpecID, iPart, NbrOfSEE
REAL                  :: RealNbrOfSEE, TimeScalingFactor, MPF
REAL                  :: Particle_pos(1:3), xi(2)
REAL                  :: RandVal, RandVal2(2), xiab(1:2,1:2), nVec(3), tang1(3), tang2(3), Velo3D(3)
#if USE_HDG
INTEGER               :: iBC,iUniqueFPCBC,iUniqueEPCBC,BCState
#endif /*USE_HDG*/
!===================================================================================================================================
! Check if ray tracing based SEE is active
! 1) Boundary from which rays are emitted
IF(RayPartBound.LE.0) RETURN
! 2) SEE yield for any BC greater than zero
IF(.NOT.ANY(PartBound%PhotonSEEYield(:).GT.0.)) RETURN

! TODO: Copied here from InitParticleMesh, which is only build if not TriaSurfaceFlux
IF(.NOT.ALLOCATED(BezierSampleXi)) ALLOCATE(BezierSampleXi(0:nSurfSample))
DO iSample=0,nSurfSample
  BezierSampleXi(iSample)=-1.+2.0/nSurfSample*iSample
END DO

! Surf sides are shared, array calculation can be distributed
#if USE_MPI
firstSide = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nSurfTotalSides
#endif /*USE_MPI*/

ASSOCIATE( tau         => Ray%PulseDuration      ,&
           tShift      => Ray%tShift             ,&
           lambda      => Ray%WaveLength         ,&
           Period      => Ray%Period)

  ! Temporal bound of integration
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

! Calculate the current pulse
NbrOfRepetitions = INT(Time/Period)

! Add arbitrary time shift (-4 sigma_t) so that I_max is not at t=0s
! Note that sigma_t = tau / sqrt(2)
t_1 = t_1 - tShift - NbrOfRepetitions * Period
t_2 = t_2 - tShift - NbrOfRepetitions * Period

! check if t_2 is outside of the pulse
IF(t_2.GT.2.0*tShift) t_2 = 2.0*tShift

TimeScalingFactor = 0.5 * SQRT(PI) * tau * (ERF(t_2/tau)-ERF(t_1/tau))

DO iSurfSide = firstSide, lastSide
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
  ! TODO: Skip sides which are not mine in the MPI case
  BCID = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  ! Skip non-reflective BC sides
  IF(PartBound%TargetBoundCond(BCID).NE.PartBound%ReflectiveBC) CYCLE
  ! Skip BC sides with zero yield
  IF(PartBound%PhotonSEEYield(BCID).LE.0.) CYCLE
  ! Determine which species is to be inserted
  SpecID = PartBound%PhotonSEEElectronSpecies(BCID)
  ! Determine which element the particles are going to be inserted
  GlobElemID = SideInfo_Shared(SIDE_ELEMID ,SideID)
  ! Determine the weighting factor of the electron species
  IF(usevMPF)THEN
    MPF = PartBound%PhotonSEEMacroParticleFactor(BCID) ! Use SEE-specific MPF
  ELSE
    MPF = Species(SpecID)%MacroParticleFactor ! Use species MPF
  END IF ! usevMPF
  ! Loop over the subsides
  DO p = 1, nSurfSample
    DO q = 1, nSurfSample
      ! Calculate the number of SEEs per subside
#if USE_MPI
      E_Intensity = PhotonSampWall_Shared(2,p,q,iSurfSide) * TimeScalingFactor
#else
      E_Intensity = PhotonSampWall(2,p,q,iSurfSide) * TimeScalingFactor
#endif /*USE_MPI*/
      RealNbrOfSEE = E_Intensity / CalcPhotonEnergy(lambda) * PartBound%PhotonSEEYield(BCID) / MPF
      CALL RANDOM_NUMBER(RandVal)
      NbrOfSEE = INT(RealNbrOfSEE+RandVal)
      ! Calculate the normal & tangential vectors
      xi(1)=(BezierSampleXi(p-1)+BezierSampleXi(p))/2. ! (a+b)/2
      xi(2)=(BezierSampleXi(q-1)+BezierSampleXi(q))/2. ! (a+b)/2
      xiab(1,1:2)=(/BezierSampleXi(p-1),BezierSampleXi(p)/)
      xiab(2,1:2)=(/BezierSampleXi(q-1),BezierSampleXi(q)/)
      CALL CalcNormAndTangBezier(nVec,tang1,tang2,xi(1),xi(2),SideID)
      ! Normal vector provided by the routine points outside of the domain
      nVec = -nVec
      ! Loop over number of particles to be inserted
      DO iPart = 1, NbrOfSEE
        ! Determine particle position within the sub-side
        CALL RANDOM_NUMBER(RandVal2)
        xi=(xiab(:,2)-xiab(:,1))*RandVal2+xiab(:,1)
        CALL EvaluateBezierPolynomialAndGradient(xi,NGeo,3,BezierControlPoints3D(1:3,0:NGeo,0:NGeo,SideID),Point=Particle_pos(1:3))
        ! Determine particle velocity
        CALL CalcVelocity_FromWorkFuncSEE(PartBound%PhotonSEEWorkFunction(BCID), Species(SpecID)%MassIC, tang1, nVec, Velo3D)
        ! Create new particle
        CALL CreateParticle(SpecID,Particle_pos(1:3),GlobElemID,Velo3D(1:3),0.,0.,0.,NewPartID=PartID,NewMPF=MPF)
        ! 1. Store the particle information in PartStateBoundary.h5
        IF(DoBoundaryParticleOutputHDF5) THEN
          CALL StoreBoundaryParticleProperties(PartID,SpecID,PartState(1:3,PartID),&
                UNITVECTOR(PartState(4:6,PartID)),nVec,iPartBound=BCID,mode=2,MPF_optIN=MPF)
        END IF ! DoBoundaryParticleOutputHDF5
#if USE_HDG
        ! 2. Check if floating boundary conditions (FPC) are used and consider electron holes
        IF(UseFPC)THEN
          iBC = PartBound%MapToFieldBC(BCID)
          IF(iBC.LE.0) CALL abort(__STAMP__,'iBC = PartBound%MapToFieldBC(PartBCIndex) must be >0',IntInfoOpt=iBC)
          IF(BoundaryType(iBC,BC_TYPE).EQ.20)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
            BCState = BoundaryType(iBC,BC_STATE) ! State is iFPC
            iUniqueFPCBC = FPC%Group(BCState,2)
            FPC%ChargeProc(iUniqueFPCBC) = FPC%ChargeProc(iUniqueFPCBC) - Species(SpecID)%ChargeIC * MPF ! Use negative charge!
          END IF ! BCType.EQ.20
        END IF ! UseFPC
        ! 3. Check if electric potential condition (EPC) are used and consider electron holes
        IF(UseEPC)THEN
          iBC = PartBound%MapToFieldBC(BCID)
          IF(iBC.LE.0) CALL abort(__STAMP__,'iBC = PartBound%MapToFieldBC(PartBCIndex) must be >0',IntInfoOpt=iBC)
          IF(BoundaryType(iBC,BC_TYPE).EQ.8)THEN ! BCType = BoundaryType(iBC,BC_TYPE)
            BCState = BoundaryType(iBC,BC_STATE) ! State is iEPC
            iUniqueEPCBC = EPC%Group(BCState,2)
            EPC%ChargeProc(iUniqueEPCBC) = EPC%ChargeProc(iUniqueEPCBC) - Species(SpecID)%ChargeIC * MPF ! Use negative charge!
          END IF ! BCType.EQ.8
        END IF ! UseEPC
#endif /*USE_HDG*/
      END DO
    END DO ! q = 1, nSurfSample
  END DO ! p = 1, nSurfSample
END DO

END ASSOCIATE

END SUBROUTINE PhotoIonization_RayTracing_SEE

END MODULE MOD_Particle_Photoionization
