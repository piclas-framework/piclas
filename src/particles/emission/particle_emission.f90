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

MODULE MOD_part_emission
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
INTERFACE ParticleInserting
  MODULE PROCEDURE ParticleInserting
END INTERFACE

INTERFACE AdaptiveBCAnalyze
  MODULE PROCEDURE AdaptiveBCAnalyze
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: ParticleInserting, AdaptiveBCAnalyze
!===================================================================================================================================
CONTAINS

SUBROUTINE ParticleInserting()
!===================================================================================================================================
!> Particle emission at every iteration (ParticleEmissionType.GT.0)
!===================================================================================================================================
! Modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: c,PI
USE MOD_Timedisc_Vars          ,ONLY: dt,time
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac,RKdtFracTotal
USE MOD_Particle_Vars
USE MOD_PIC_Vars
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: useDSMC, CollisMode, SpecDSMC
USE MOD_part_emission_tools    ,ONLY: DSMC_SetInternalEnr_LauxVFD
USE MOD_DSMC_PolyAtomicModel   ,ONLY: DSMC_SetInternalEnr_Poly
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcPartBalance,nPartIn,PartEkinIn
USE MOD_Particle_Analyze_Tools ,ONLY: CalcEkinPart
USE MOD_part_emission_tools    ,ONLY: SetParticleChargeAndMass,SetParticleMPF,SamplePoissonDistri,SetParticleTimeStep,CalcNbrOfPhotons
USE MOD_part_pos_and_velo      ,ONLY: SetParticlePosition,SetParticleVelocity
USE MOD_DSMC_BGGas             ,ONLY: BGGas_PhotoIonization
USE MOD_DSMC_ChemReact         ,ONLY: CalcPhotoIonizationNumber
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_ReadInTools            ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Local variable declaration
INTEGER                          :: i , iPart, PositionNbr, iInit, IntSample
INTEGER                          :: NbrOfParticle
INTEGER(KIND=8)                  :: inserted_Particle_iter,inserted_Particle_time
INTEGER(KIND=8)                  :: inserted_Particle_diff
REAL                             :: PartIns, RandVal1
REAL                             :: RiseFactor, RiseTime,NbrOfPhotons
#if USE_MPI
INTEGER                          :: InitGroup
#endif
REAL                             :: NbrOfReactions,NbrOfParticlesReal
!===================================================================================================================================

!---  Emission at time step
DO i=1,nSpecies
  DO iInit = 1, Species(i)%NumberOfInits
    ! Reset the number of particles per species AND init region
    NbrOfParticle = 0
    ! Only use inits defined for emission (every time step)
    IF (Species(i)%Init(iInit)%ParticleEmissionType.GT.0) THEN
      SELECT CASE(Species(i)%Init(iInit)%ParticleEmissionType)
      CASE(1) ! Emission Type: Particles per !!!!!SECOND!!!!!!!! (not per ns)
        IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dt*RKdtFrac  ! emitted particles during time-slab
          inserted_Particle_iter = INT(PartIns,8)                                     ! number of particles to be inserted
          PartIns=Species(i)%Init(iInit)%ParticleNumber * (Time + dt*RKdtFracTotal) ! total number of emitted particle over
                                                                                      ! simulation
          !-- random-round the inserted_Particle_time for preventing periodicity
          ! PO & SC: why, sometimes we do not want this add, TB is bad!
          IF (inserted_Particle_iter.GE.1) THEN
            CALL RANDOM_NUMBER(RandVal1)
            inserted_Particle_time = INT(PartIns + RandVal1,8) ! adds up to ONE
          ELSE IF((inserted_Particle_iter.GE.0).AND.(inserted_Particle_iter.LT.1)) THEN
                                                      !needed, since InsertedParticleSurplus can increase
                                                      !and _iter>1 needs to be possible for preventing periodicity
            IF (ALMOSTEQUAL(PartIns,0.)) THEN !dummy
              inserted_Particle_time = INT(PartIns,8)
            ELSE !poisson-distri of PartIns-INT(PartIns)
              CALL SamplePoissonDistri( PartIns-INT(PartIns) , IntSample )
              inserted_Particle_time = INT(INT(PartIns)+IntSample,8) !INT(PartIns) + POISDISTRI( PartIns-INT(PartIns) )
            END IF
          ELSE !dummy
            inserted_Particle_time = INT(PartIns,8)
          END IF
          !-- evaluate inserted_Particle_time and inserted_Particle_iter
          inserted_Particle_diff = inserted_Particle_time - Species(i)%Init(iInit)%InsertedParticle &
            - inserted_Particle_iter - Species(i)%Init(iInit)%InsertedParticleSurplus &
            + Species(i)%Init(iInit)%InsertedParticleMisMatch
          Species(i)%Init(iInit)%InsertedParticleSurplus = ABS(MIN(inserted_Particle_iter + inserted_Particle_diff,0))
          NbrOfParticle = MAX(INT(inserted_Particle_iter + inserted_Particle_diff,4),0)
          !-- if maxwell velo dist and less than 5 parts: skip (to ensure maxwell dist)
          IF (TRIM(Species(i)%Init(iInit)%velocityDistribution).EQ.'maxwell') THEN
            IF (NbrOfParticle.LT.5) NbrOfParticle=0
          END IF
        ELSE IF (DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
          ! linear rise of inflow
          RiseTime=Species(i)%Init(iInit)%InflowRiseTime
          IF(RiseTime.GT.0.)THEN
            IF(Time-DelayTime.LT.RiseTime)THEN
              RiseFactor=(time-DelayTime)/RiseTime
            ELSE
              RiseFactor=1.
            END IF
          ELSE
            RiseFactor=1.
          EnD IF
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dt*RKdtFrac * RiseFactor  ! emitted particles during time-slab
          CALL RANDOM_NUMBER(RandVal1)
          IF (EXP(-PartIns).LE.TINY(PartIns)) THEN
            IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
            NbrOfParticle = INT(PartIns + RandVal1)
            DoPoissonRounding = .FALSE.
          ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
            CALL SamplePoissonDistri( PartIns , NbrOfParticle , DoPoissonRounding)
          END IF
        ELSE ! DoTimeDepInflow
          ! linear rise of inflow
          RiseTime=Species(i)%Init(iInit)%InflowRiseTime
          IF(RiseTime.GT.0.)THEN
            IF(Time-DelayTime.LT.RiseTime)THEN
              RiseFactor=(time-DelayTime)/RiseTime
            ELSE
              RiseFactor=1.
            END IF
          ELSE
            RiseFactor=1.
          EnD IF
          ! emitted particles during time-slab
          PartIns=Species(i)%Init(iInit)%ParticleNumber * dt*RKdtFrac * RiseFactor &
                  + Species(i)%Init(iInit)%InsertedParticleMisMatch
          CALL RANDOM_NUMBER(RandVal1)
          NbrOfParticle = INT(PartIns + RandVal1)
        END IF
#if USE_MPI
        ! communicate number of particles with all procs in the same init group
        InitGroup=Species(i)%Init(iInit)%InitCOMM
        IF(PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
          ! only procs which are part of group take part in the communication
            !NbrOfParticle based on RandVals!
          CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
        ELSE
          NbrOfParticle=0
        END IF
        !CALL MPI_BCAST(NbrOfParticle, 1, MPI_INTEGER,0,PartMPI%COMM,IERROR) !NbrOfParticle based on RandVals!
#endif
        Species(i)%Init(iInit)%InsertedParticle = Species(i)%Init(iInit)%InsertedParticle + INT(NbrOfParticle,8)
      CASE(2)    ! Emission Type: Particles per Iteration
        IF (RKdtFracTotal .EQ. 1.) THEN !insert in last stage only, so that no reconstruction is nec. and number/iter matches
          NbrOfParticle = INT(Species(i)%Init(iInit)%ParticleNumber)
        ELSE
          NbrOfParticle = 0
        END IF
      CASE(7) ! SEE based on photon impact and photo-ionization in the volume
        ASSOCIATE( tShift => Species(i)%Init(iInit)%tShift )
          ! Check if all pulses have terminated
          IF(Time.LE.Species(i)%Init(iInit)%tActive)THEN
            ! Check if pulse is currently active of in between two pulses (in the latter case, do nothing)
            IF(MOD(MERGE(Time-tShift, Time, Time.GE.tShift), Species(i)%Init(iInit)%Period).LE.2.0*tShift)THEN
              ! Calculate the number of currently active photons (both surface SEE and volumetric emission)
              CALL CalcNbrOfPhotons(i, iInit, NbrOfPhotons)

              ! Check if only particles in the first quadrant are to be inserted that is spanned by the vectors
              ! x=BaseVector1IC and y=BaseVector2IC in the interval x,y in [0,R] and reduce the number of photon accordingly
              IF(Species(i)%Init(iInit)%FirstQuadrantOnly) NbrOfPhotons = NbrOfPhotons / 4.0

              ! Select surface SEE or volumetric emission
              IF(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'photon_SEE_disc')THEN
                ! SEE based on photon impact
                NbrOfPhotons = Species(i)%Init(iInit)%YieldSEE * NbrOfPhotons / Species(i)%MacroParticleFactor &
                              + Species(i)%Init(iInit)%NINT_Correction
                NbrOfParticle = NINT(NbrOfPhotons)
                Species(i)%Init(iInit)%NINT_Correction = NbrOfPhotons - REAL(NbrOfParticle)
              ELSE
                ! Photo-ionization in the volume
                ! Calculation of the number of photons (using actual number and applying the weighting factor on the number of reactions)
                NbrOfPhotons = Species(i)%Init(iInit)%EffectiveIntensityFactor * NbrOfPhotons
                ! Calculation of the number of photons depending on the cylinder height (ratio of actual to virtual cylinder height, which
                ! is spanned by the disk and the length given by c*dt)
                NbrOfPhotons = NbrOfPhotons * Species(i)%Init(iInit)%CylinderHeightIC / (c*dt)
                ! Calculation of the number of electron resulting from the chemical reactions in the photoionization region
                CALL CalcPhotoIonizationNumber(NbrOfPhotons,NbrOfReactions)
                NbrOfReactions = NbrOfReactions + Species(i)%Init(iInit)%NINT_Correction
                NbrOfParticle = NINT(NbrOfReactions)
                Species(i)%Init(iInit)%NINT_Correction = NbrOfReactions - REAL(NbrOfParticle)
              END IF
            ELSE
              NbrOfParticle = 0
            END IF ! MOD(MERGE(Time-T0/2., Time, Time.GE.T0/2.), Period).GT.T0
          ELSE
            NbrOfParticle = 0
          END IF ! Time.LE.Species(i)%Init(iInit)%tActive
        END ASSOCIATE
      CASE(8) ! SpaceIC='2D_landmark','2D_landmark_copy'
              ! Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
              ! for low-temperature partially magnetized plasmas (2019)
       ASSOCIATE( x2 => 1.0e-2  ,& ! m
                  x1 => 0.25e-2 ,& ! m
                  Ly => 1.28e-2 ,& ! m
                  S0 => 5.23e23 )  ! m^-3 s^-1
         NbrOfParticlesReal = Ly*dt*S0*2.0*(x2-x1)/PI ! yields 1.60E+08 m^-1 (i.e. per metre in z-direction)
         NbrOfParticlesReal = NbrOfParticlesReal*(GEO%zmaxglob-GEO%zminglob)/Species(i)%MacroParticleFactor &
                              + Species(i)%Init(iInit)%NINT_Correction
         NbrOfParticle = NINT(NbrOfParticlesReal)
         Species(i)%Init(iInit)%NINT_Correction = NbrOfParticlesReal - REAL(NbrOfParticle)
         IF(.NOT.ALLOCATED(PartPosLandmark))THEN
           NbrOfParticleLandmarkMax = MAX(10,NINT(NbrOfParticlesReal*1.2)) ! Add 20% safety
           CALL PrintOption('Landmark volume emission allocation size NbrOfParticleLandmarkMax' , 'CALCUL.' , IntOpt=NbrOfParticleLandmarkMax)
           ALLOCATE(PartPosLandmark(1:3,1:NbrOfParticleLandmarkMax))
           PartPosLandmark=HUGE(1.)
           IF(NbrOfParticleLandmarkMax.LE.0)THEN
             IPWRITE(UNIT_StdOut,*) "NbrOfParticleLandmarkMax =", NbrOfParticleLandmarkMax
             CALL abort(&
                 __STAMP__&
                 ,'NbrOfParticleLandmarkMax.LE.0')
           END IF
         ELSE
           IF(NbrOfParticleLandmarkMax.LT.NbrOfParticle)THEN
             IPWRITE(UNIT_StdOut,*) "NbrOfParticleLandmarkMax,NbrOfParticle =", NbrOfParticleLandmarkMax,NbrOfParticle
             CALL abort(&
             __STAMP__&
             ,'NbrOfParticleLandmarkMax.LT.NbrOfParticle is not allowed! Allocate PartPosLandmark to the appropriate size.')
           END IF ! NbrOfParticleLandmarkMax.LE.NbrOfParticle
         END IF ! .NOT.ALLOCATED()
       END ASSOCIATE
     CASE(9) ! '2D_landmark_neutralization'
#if USE_MPI
       ! Communicate number of particles with all procs in the same init group
       InitGroup=Species(i)%Init(iInit)%InitCOMM
       IF(PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL) THEN
         ! Only processors which are part of group take part in the communication
         CALL MPI_ALLREDUCE(NeutralizationBalance,NeutralizationBalanceGlobal,1,MPI_INTEGER,MPI_SUM,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
       ELSE
         NeutralizationBalanceGlobal=0
       END IF
#else
       NeutralizationBalanceGlobal = NeutralizationBalance
#endif
       ! Insert electrons only when the number is greater than zero
       IF(NeutralizationBalanceGlobal.GT.0)THEN
         ! Insert only when positive
         NbrOfParticle = NeutralizationBalanceGlobal
         ! Reset the counter
         NeutralizationBalance = 0
       ELSE
         NbrOfParticle = 0
       END IF ! NeutralizationBalance.GT.0

      CASE DEFAULT
        NbrOfParticle = 0
      END SELECT

      CALL SetParticlePosition(i,iInit,NbrOfParticle)
      ! Pairing of "electrons" with the background species and performing the reaction
      IF(TRIM(Species(i)%Init(iInit)%SpaceIC).EQ.'photon_cylinder') THEN
        CALL BGGas_PhotoIonization(i,iInit,NbrOfParticle)
        CYCLE
      END IF

      CALL SetParticleVelocity(i,iInit,NbrOfParticle)
      CALL SetParticleChargeAndMass(i,NbrOfParticle)
      IF (usevMPF) CALL SetParticleMPF(i,NbrOfParticle)
      IF (VarTimeStep%UseVariableTimeStep) CALL SetParticleTimeStep(NbrOfParticle)
      ! define molecule stuff
      IF (useDSMC.AND.(CollisMode.GT.1)) THEN
        iPart = 1
        DO WHILE (iPart.LE.NbrOfParticle)
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr.NE.0) THEN
            IF (SpecDSMC(i)%PolyatomicMol) THEN
              CALL DSMC_SetInternalEnr_Poly(i,iInit,PositionNbr,1)
            ELSE
              CALL DSMC_SetInternalEnr_LauxVFD(i,iInit,PositionNbr,1)
            END IF
          END IF
          iPart = iPart + 1
        END DO
      END IF
      ! Compute number of input particles and energy
      IF(CalcPartBalance.AND.(NbrOfParticle.GT.0)) THEN
        nPartIn(i)=nPartIn(i) + NbrOfparticle
        DO iPart=1,NbrOfParticle
          PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
          IF (PositionNbr .NE. 0) PartEkinIn(i) = PartEkinIn(i) + CalcEkinPart(PositionNbr)
        END DO ! iPart
      END IF ! CalcPartBalance
      ! Update the current next free position and increase the particle vector length
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + NbrOfParticle
      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfParticle

      ! Complete check if all particles were emitted successfully
#if USE_MPI
      InitGroup=Species(i)%Init(iInit)%InitCOMM
      IF (PartMPI%InitGroup(InitGroup)%COMM.NE.MPI_COMM_NULL .AND. Species(i)%Init(iInit)%sumOfRequestedParticles.GT.0) THEN
        CALL MPI_WAIT(PartMPI%InitGroup(InitGroup)%Request, MPI_STATUS_IGNORE, iError)

        IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
        ! add number of matching error to particle emission to fit
        ! number of added particles
        Species(i)%Init(iInit)%InsertedParticleMisMatch = Species(i)%Init(iInit)%sumOfRequestedParticles - Species(i)%Init(iInit)%sumOfMatchedParticles
        IF (Species(i)%Init(iInit)%sumOfRequestedParticles .GT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
          WRITE(UNIT_StdOut,'(A)')      'WARNING in ParticleEmission_parallel:'
          WRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
          WRITE(UNIT_StdOut,'(A,I0,A)') 'matched only '   , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
          WRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
        ELSE IF (Species(i)%Init(iInit)%sumOfRequestedParticles .LT. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
          WRITE(UNIT_StdOut,'(A)')      'ERROR in ParticleEmission_parallel:'
          WRITE(UNIT_StdOut,'(A,I0)')   'Fraction Nbr: '  , i
          WRITE(UNIT_StdOut,'(A,I0,A)') 'matched '        , Species(i)%Init(iInit)%sumOfMatchedParticles  , ' particles'
          WRITE(UNIT_StdOut,'(A,I0,A)') 'when '           , Species(i)%Init(iInit)%sumOfRequestedParticles, ' particles were required!'
        ! ELSE IF (nbrOfParticle .EQ. Species(i)%Init(iInit)%sumOfMatchedParticles) THEN
        !  WRITE(UNIT_stdOut,'(A,I0)')   'Fraction Nbr: '  , FractNbr
        !  WRITE(UNIT_stdOut,'(A,I0,A)') 'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
        END IF
#if USE_MPI
        END IF ! PartMPI%iProc.EQ.0
      END IF
#endif
    END IF ! Species(iSpec)%Init(iInit)%ParticleEmissionType.GT.0
  END DO  ! iInit
END DO  ! i=1,nSpecies

END SUBROUTINE ParticleInserting


SUBROUTINE AdaptiveBCAnalyze(initSampling_opt)
!===================================================================================================================================
! Sampling of variables (part-density, velocity and energy) for Adaptive BC elements
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_DSMC_Analyze           ,ONLY: CalcTVibPoly,CalcTelec
USE MOD_DSMC_Vars              ,ONLY: PartStateIntEn, DSMC, CollisMode, SpecDSMC, useDSMC, RadialWeighting
USE MOD_Mesh_Vars              ,ONLY: nElems, offsetElem
USE MOD_Part_Tools             ,ONLY: GetParticleWeight
USE MOD_SurfaceModel_Vars      ,ONLY: PorousBCSampIter, PorousBCMacroVal
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared, SideInfo_Shared
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared
USE MOD_Particle_Vars          ,ONLY: PartState, PDM, PartSpecies, Species, nSpecies, PEM, Adaptive_MacroVal, AdaptiveWeightFac
USE MOD_Particle_Vars          ,ONLY: usevMPF
USE MOD_Timedisc_Vars          ,ONLY: iter
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime, LBElemSplitTime, LBPauseTime
USE MOD_LoadBalance_vars       ,ONLY: nPartsPerBCElem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL, INTENT(IN), OPTIONAL   :: initSampling_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: ElemID, AdaptiveElemID, i, iSpec, SamplingIteration
REAL                            :: TVib_TempFac, TTrans_TempFac, RelaxationFactor
REAL, ALLOCATABLE               :: Source(:,:,:)
LOGICAL                         :: initSampling, isBCElem
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER                         :: nlocSides, GlobalSideID, GlobalElemID, iLocSide
!===================================================================================================================================
ALLOCATE(Source(1:11,1:nElems,1:nSpecies))
Source=0.0

! Optional flag for the utilization of the routine for an initial sampling of the density and pressure distribution before simstart
IF(PRESENT(initSampling_opt)) THEN
  initSampling = initSampling_opt
ELSE
 initSampling = .FALSE.
END IF

! If no particles are present during the initial sampling for a porous BC, leave the routine, otherwise initial variables for the
! adaptive inlet surface flux will be overwritten by zero's.
IF (PDM%ParticleVecLength.LT.1) RETURN

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
DO i=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(i)) THEN
    GlobalElemID = PEM%GlobalElemID(i)
    ElemID = GlobalElemID - offsetElem
    IF ((ElemID.LT.1).OR.(ElemID.GT.nElems)) CYCLE
    ! not a BC element
    ! IF (ElemToBCSides(ELEM_NBR_BCSIDES,ElemID).EQ.-1) CYCLE
    ! ======================================
    ! ElemToBCSides is only built for RefMapping, eigenes Array beschreiben statt untere schleife jedes mal abzulaufen
    isBCElem = .FALSE.
    nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
    DO iLocSide = 1,nlocSides
      GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
      IF (SideInfo_Shared(SIDE_BCID,GlobalSideID).GT.0) isBCElem = .TRUE.
    END DO
    IF(.NOT.isBCElem) CYCLE
    ! ======================================
#if USE_LOADBALANCE
    nPartsPerBCElem(ElemID) = nPartsPerBCElem(ElemID) + 1
#endif /*USE_LOADBALANCE*/
    iSpec = PartSpecies(i)
    Source(1:3,ElemID, iSpec) = Source(1:3,ElemID,iSpec) + PartState(4:6,i) * GetParticleWeight(i)
    Source(4:6,ElemID, iSpec) = Source(4:6,ElemID,iSpec) + PartState(4:6,i)**2 * GetParticleWeight(i)
    Source(7,ElemID, iSpec) = Source(7,ElemID, iSpec) + 1.0  ! simulation particle number
    IF(useDSMC)THEN
      IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
        IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
          Source(8:9,ElemID, iSpec) = Source(8:9,ElemID, iSpec) + PartStateIntEn(1:2,i) * GetParticleWeight(i)
        END IF
      END IF
      IF (DSMC%ElectronicModel.GT.0) THEN
        Source(10,ElemID, iSpec) = Source(10,ElemID, iSpec) + PartStateIntEn(3,i) * GetParticleWeight(i)
      END IF
    END IF
    Source(11,ElemID, iSpec) = Source(11,ElemID, iSpec) + GetParticleWeight(i)
  END IF
END DO
#if USE_LOADBALANCE
CALL LBPauseTime(LB_ADAPTIVE,tLBStart)
#endif /*USE_LOADBALANCE*/

IF(initSampling) THEN
  RelaxationFactor = 1
  SamplingIteration = 1
ELSE
  RelaxationFactor = AdaptiveWeightFac
  SamplingIteration = PorousBCSampIter
END IF

DO AdaptiveElemID = 1,nElems
  GlobalElemID = AdaptiveElemID + offsetElem
  ! not a BC element
  ! IF (ElemToBCSides(ELEM_NBR_BCSIDES,AdaptiveElemID).EQ.-1) CYCLE
  ! ======================================
  ! ElemToBCSides is only built for RefMapping, eigenes Array beschreiben statt untere schleife jedes mal abzulaufen?
  isBCElem = .FALSE.
  nlocSides = ElemInfo_Shared(ELEM_LASTSIDEIND,GlobalElemID) -  ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID)
  DO iLocSide = 1,nlocSides
    GlobalSideID = ElemInfo_Shared(ELEM_FIRSTSIDEIND,GlobalElemID) + iLocSide
    IF (SideInfo_Shared(SIDE_BCID,GlobalSideID).GT.0) isBCElem = .TRUE.
  END DO
  IF(.NOT.isBCElem) CYCLE
  ! ======================================
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart)
#endif /*USE_LOADBALANCE*/
  DO iSpec = 1,nSpecies
    ! write timesample particle values of bc elements in global macrovalues of bc elements
    IF (Source(7,AdaptiveElemID,iSpec).GT.0.0) THEN
      IF(.NOT.initSampling) THEN
        ! compute flow velocity (during computation, not for the initial distribution, where the velocity from the ini is used)
        Adaptive_MacroVal(1:3,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(1:3,AdaptiveElemID,iSpec) &
            + RelaxationFactor*Source(1:3,AdaptiveElemID, iSpec) / Source(11,AdaptiveElemID,iSpec)
      END IF
      ! compute flow Temperature
      IF (Source(7,AdaptiveElemID,iSpec).GT.1.0) THEN
        Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) &
          + RelaxationFactor &
            * (Source(7,AdaptiveElemID,iSpec)/(Source(7,AdaptiveElemID,iSpec)-1.0)) &
            * Species(iSpec)%MassIC/ BoltzmannConst &
            * ( Source(4:6,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec) &
            - (Source(1:3,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2)
      ELSE
        Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(4:6,AdaptiveElemID,iSpec)
      END IF
      ! ================================================================
      IF(PorousBCSampIter.GT.0) THEN
        ! Sampling the number density and pressure every given number of iterations and RESETTING it after calculation
        PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) = PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) + Source(1:6,AdaptiveElemID, iSpec)
        ! Sampling the particle weights
        PorousBCMacroVal(7,AdaptiveElemID,iSpec) = PorousBCMacroVal(7,AdaptiveElemID,iSpec) + Source(11,AdaptiveElemID,iSpec)
        ! Sampling the number of simulation particles
        PorousBCMacroVal(8,AdaptiveElemID,iSpec) = PorousBCMacroVal(8,AdaptiveElemID,iSpec) + Source(7,AdaptiveElemID,iSpec)
        IF(MOD(iter,SamplingIteration).EQ.0) THEN
          IF(PorousBCMacroVal(8,AdaptiveElemID,iSpec).GT.1) THEN
            ! number density
            IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
              Adaptive_MacroVal(7,AdaptiveElemID,iSpec)=PorousBCMacroVal(7,AdaptiveElemID,iSpec)/SamplingIteration   &
                                                        /ElemVolume_Shared(GlobalElemID)
            ELSE
              Adaptive_MacroVal(7,AdaptiveElemID,iSpec)=PorousBCMacroVal(7,AdaptiveElemID,iSpec)/SamplingIteration   &
                                                        /ElemVolume_Shared(GlobalElemID) * Species(iSpec)%MacroParticleFactor
            END IF
            ! instantaneous temperature WITHOUT 1/BoltzmannConst
            PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) = PorousBCMacroVal(1:6,AdaptiveElemID,iSpec) &
                                                          / PorousBCMacroVal(7,AdaptiveElemID,iSpec)
            TTrans_TempFac = (PorousBCMacroVal(8,AdaptiveElemID,iSpec)/(PorousBCMacroVal(8,AdaptiveElemID,iSpec)-1.0)) &
                *Species(iSpec)%MassIC*(PorousBCMacroVal(4,AdaptiveElemID,iSpec) - PorousBCMacroVal(1,AdaptiveElemID,iSpec)**2   &
                                      + PorousBCMacroVal(5,AdaptiveElemID,iSpec) - PorousBCMacroVal(2,AdaptiveElemID,iSpec)**2   &
                                      + PorousBCMacroVal(6,AdaptiveElemID,iSpec) - PorousBCMacroVal(3,AdaptiveElemID,iSpec)**2) / 3.
            ! pressure (BoltzmannConstant canceled out in temperature calculation)
            Adaptive_MacroVal(12,AdaptiveElemID,iSpec)=Adaptive_MacroVal(7,AdaptiveElemID,iSpec)*TTrans_TempFac
            ! Resetting the sampling values
            PorousBCMacroVal(1:8,AdaptiveElemID,iSpec) = 0.0
          END IF
        END IF
      ELSE
        ! Calculation of the number density and pressure with the relaxation factor
        ! compute instantaneous temperature WITHOUT 1/BoltzmannConst
        IF (Source(7,AdaptiveElemID,iSpec).GT.1.0) THEN
          TTrans_TempFac = (Source(7,AdaptiveElemID,iSpec)/(Source(7,AdaptiveElemID,iSpec)-1.0)) &
                          * Species(iSpec)%MassIC * (Source(4,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)       &
                          - (Source(1,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2   &
                          + Source(5,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)      &
                          - (Source(2,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2   &
                          + Source(6,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec)      &
                          - (Source(3,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec))**2) / 3.
        ELSE
          TTrans_TempFac = 0.0
        END IF
        IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) THEN
          ! compute density
          Adaptive_MacroVal(7,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(7,AdaptiveElemID,iSpec) &
            + RelaxationFactor*Source(11,AdaptiveElemID,iSpec) /ElemVolume_Shared(GlobalElemID)
          ! Pressure with relaxation factor
          Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
          +RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/ElemVolume_Shared(GlobalElemID)*TTrans_TempFac
        ELSE
          ! compute density
          Adaptive_MacroVal(7,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(7,AdaptiveElemID,iSpec) &
            + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)/ElemVolume_Shared(GlobalElemID)*Species(iSpec)%MacroParticleFactor
          ! Pressure with relaxation factor
          Adaptive_MacroVal(12,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(12,AdaptiveElemID,iSpec) &
                                                      + RelaxationFactor*Source(11,AdaptiveElemID,iSpec)  &
                                                      / ElemVolume_Shared(GlobalElemID)*Species(iSpec)%MacroParticleFactor*TTrans_TempFac
        END IF
      END IF
      ! !==================================================================================================
      IF(useDSMC)THEN
        IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
        IF ((SpecDSMC(iSpec)%InterID.EQ.2).OR.(SpecDSMC(iSpec)%InterID.EQ.20)) THEN
            IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
              IF( (Source(8,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec)) .GT. SpecDSMC(iSpec)%EZeroPoint) THEN
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec) &
                  + RelaxationFactor*CalcTVibPoly(Source(8,AdaptiveElemID,iSpec) / Source(11,AdaptiveElemID,iSpec),iSpec)
              ELSE
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec)
              END IF
            ELSE
              TVib_TempFac=Source(8,AdaptiveElemID,iSpec)/ (Source(11,AdaptiveElemID,iSpec) &
                *BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
              IF (TVib_TempFac.LE.DSMC%GammaQuant) THEN
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec)
              ELSE
                Adaptive_MacroVal(8,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(8,AdaptiveElemID,iSpec) &
                  + RelaxationFactor*SpecDSMC(iSpec)%CharaTVib / LOG(1 + 1/(TVib_TempFac-DSMC%GammaQuant))
              END IF
            END IF
            Adaptive_MacroVal(9,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(9,AdaptiveElemID,iSpec) &
                + RelaxationFactor*Source(9,AdaptiveElemID,iSpec)/(Source(11,AdaptiveElemID,iSpec)*BoltzmannConst)
            IF (DSMC%ElectronicModel.GT.0) THEN
              Adaptive_MacroVal(10,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(10,AdaptiveElemID,iSpec) &
                + RelaxationFactor*CalcTelec( Source(10,AdaptiveElemID,iSpec)/Source(11,AdaptiveElemID,iSpec),iSpec)
            END IF
          END IF
        END IF
      END IF
    ELSE
      Adaptive_MacroVal(1:10,AdaptiveElemID,iSpec) = (1-RelaxationFactor)*Adaptive_MacroVal(1:10,AdaptiveElemID,iSpec)
    END IF
  END DO
#if USE_LOADBALANCE
  CALL LBElemSplitTime(AdaptiveElemID,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO

END SUBROUTINE AdaptiveBCAnalyze

END MODULE MOD_part_emission
