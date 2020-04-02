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

MODULE MOD_part_pos_and_velo
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

INTERFACE SetParticlePosition
  MODULE PROCEDURE SetParticlePosition
END INTERFACE

INTERFACE SetParticleVelocity
  MODULE PROCEDURE SetParticleVelocity
END INTERFACE

!===================================================================================================================================
PUBLIC         :: SetParticleVelocity,SetParticlePosition
!===================================================================================================================================
CONTAINS


SUBROUTINE SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species,Symmetry2DAxisymmetric
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_part_Emission_Tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: chunkSize
LOGICAL                                  :: DoExactPartNumInsert
#if USE_MPI
INTEGER                                  :: InitGroup
REAL,ALLOCATABLE                         :: ProcMeshVol(:)
INTEGER,ALLOCATABLE                      :: ProcNbrOfParticle(:)
#endif
!===================================================================================================================================/

  DoExactPartNumInsert =  .FALSE.
  ! check if particle inserting during simulation or initial inserting and also if via partdensity or exact particle number
  ! nbrOfParticles is set for initial inserting if initialPartNum or partdensity is set in ini
  ! ParticleEmission and Partdensity not working together
  IF (NbrofParticle.EQ.0.AND.(Species(FractNbr)%Init(iInit)%ParticleEmission.EQ.0)) RETURN
  IF ((NbrofParticle.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.LE.0.)) THEN
    DoExactPartNumInsert =  .TRUE.
    IF(Symmetry2DAxisymmetric) CALL abort(&
__STAMP__&
,'Axisymmetric: Particle insertion only possible with PartDensity!')
  END IF
  chunksize = 0
#if USE_MPI
! emission group communicator
  InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
  IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
    NbrofParticle=0
    RETURN
  END IF
  IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1 .AND. Species(FractNbr)%Init(iInit)%ElemPartDensityFileID.EQ.0) THEN
    IF (DoExactPartNumInsert) THEN !###$ ToDo
      IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
        ALLOCATE(ProcMeshVol(0:PartMPI%InitGroup(InitGroup)%nProcs-1))
        ALLOCATE(ProcNbrOfParticle(0:PartMPI%InitGroup(InitGroup)%nProcs-1))
        ProcMeshVol=0.
        ProcNbrOfParticle=0
      ELSE ! to reduce global memory allocation if a lot of procs are used
        ALLOCATE(ProcMeshVol(1))
        ALLOCATE(ProcNbrOfParticle(1))
        ProcMeshVol=0.
        ProcNbrOfParticle=0
      END IF !InitGroup%MPIroot
      CALL MPI_GATHER(LocalVolume,1,MPI_DOUBLE_PRECISION &
                     ,ProcMeshVol,1,MPI_DOUBLE_PRECISION,0,PartMPI%InitGroup(InitGroup)%COMM,iError)
      IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
        CALL IntegerDivide(NbrOfParticle,PartMPI%InitGroup(InitGroup)%nProcs,ProcMeshVol,ProcNbrOfParticle)
      END IF
      CALL MPI_SCATTER(ProcNbrOfParticle, 1, MPI_INTEGER, chunksize, 1, MPI_INTEGER, 0, PartMPI%InitGroup(InitGroup)%COMM, IERROR)
      SDEALLOCATE(ProcMeshVol)
      SDEALLOCATE(ProcNbrOfParticle)
    END IF
  ELSE
    chunksize = NbrOfParticle
  END IF
#else
  IF (DoExactPartNumInsert) chunksize = NbrOfParticle
#endif /*USE_MPI*/
  IF ((chunksize.GT.0).OR.(Species(FractNbr)%Init(iInit)%PartDensity.GT.0.)) THEN
    CALL SetCellLocalParticlePosition(chunkSize,FractNbr,iInit,DoExactPartNumInsert)
  END IF
  NbrOfParticle = chunksize

END SUBROUTINE SetParticlePositionCellLocal


SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PartState, Symmetry2DAxisymmetric
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_part_emission_tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition,SetParticlePositionPoint
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionEquidistLine, SetParticlePositionLine, SetParticlePositionDisk
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionCuboidCylinder, SetParticlePositionGyrotronCircle,SetParticlePositionCircle
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionSphere, SetParticlePositionSinDeviation
#if USE_MPI
USE MOD_Particle_MPI_Emission  ,ONLY: SendEmissionParticlesToProcs
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                         :: particle_positions(:)
INTEGER                                  :: i,ParticleIndexNbr,allocStat,nChunks, chunkSize
INTEGER                                  :: mySumOfMatchedParticles, sumOfMatchedParticles, DimSend
LOGICAL                                  :: DoExactPartNumInsert
#if USE_MPI
INTEGER                                  :: InitGroup
REAL,ALLOCATABLE                         :: ProcMeshVol(:)
INTEGER,ALLOCATABLE                      :: ProcNbrOfParticle(:)
#endif
!===================================================================================================================================
IF (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
  CALL SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
  RETURN
END IF
IF ( (NbrOfParticle .LE. 0).AND. (ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.) ) RETURN

DimSend=3 !save (and send) only positions
nChunks = 1                   ! Standard: Nicht-MPI
sumOfMatchedParticles = 0
mySumOfMatchedParticles = 0
chunkSize = nbrOfParticle
! emission group communicator
#if USE_MPI
InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
IF ( (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle_equidistant'                 ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'sin_deviation'                      ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle'                             ) .OR.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'line_with_equidistant_distribution' )) THEN
  nChunks = 1
ELSE IF (nbrOfParticle.GT.(PartMPI%InitGroup(InitGroup)%nProcs*10)) THEN
  nChunks = PartMPI%InitGroup(InitGroup)%nProcs
ELSE
  nChunks = 1
END IF
chunkSize = INT(nbrOfParticle/nChunks)
IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
  chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
END IF
! all proc taking part in particle inserting
IF (PartMPI%InitGroup(InitGroup)%MPIROOT.OR.nChunks.GT.1) THEN
#endif
  ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
  IF (allocStat .NE. 0) &
    CALL abort(__STAMP__,'ERROR in SetParticlePosition: cannot allocate particle_positions!')

  !------------------SpaceIC-cases: start-----------------------------------------------------------!
  SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
  CASE ('point')
    CALL SetParticlePositionPoint(FractNbr,iInit,chunkSize,particle_positions)
  CASE ('line_with_equidistant_distribution')
    CALL SetParticlePositionEquidistLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE ('line')
    CALL SetParticlePositionLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE('disc')
    CALL SetParticlePositionDisk(FractNbr,iInit,chunkSize,particle_positions)
  CASE('circle', 'circle_equidistant')
    CALL SetParticlePositionCircle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('gyrotron_circle')
    CALL SetParticlePositionGyrotronCircle(FractNbr,iInit,chunkSize,particle_positions, NbrOfParticle)
  CASE('cuboid','cylinder')
    CALL SetParticlePositionCuboidCylinder(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sphere')
    CALL SetParticlePositionSphere(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sin_deviation')
    CALL SetParticlePositionSinDeviation(FractNbr,iInit,chunkSize,particle_positions)
  END SELECT
  !------------------SpaceIC-cases: end-------------------------------------------------------------------------------------------
#if USE_MPI
ELSE !no mpi root, nchunks=1
  chunkSize=0
END IF
! Need to open MPI communication regardless of the chunk number. Make it only dependent on the number of procs
IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1) THEN
  CALL SendEmissionParticlesToProcs(chunkSize,DimSend,FractNbr,iInit,mySumOfMatchedParticles,particle_positions)
! Finish emission on local proc
ELSE
  mySumOfMatchedParticles = 0
  ParticleIndexNbr        = 1
  DO i=1,chunkSize
    ! Find a free position in the PDM array
    IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
      ParticleIndexNbr = PDM%nextFreePosition(mySumOfMatchedParticles + 1 + PDM%CurrentNextFreePosition)
    END IF
    IF (ParticleIndexNbr.NE.0) THEN
      PartState(1:DimSend,ParticleIndexNbr) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
      PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
      CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.FALSE.)
      IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
        mySumOfMatchedParticles = mySumOfMatchedParticles + 1
        PDM%IsNewPart(ParticleIndexNbr) = .TRUE.
        PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
      END IF
    ELSE
          CALL abort(&
    __STAMP__&
    ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO
END IF
! we want always warnings to know if the emission has failed. if a timedisc does not require this, this
! timedisc has to be handled separately
! check the sum of the matched particles: did each particle find its "home"-CPU?
CALL MPI_ALLREDUCE( mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM &
                  , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
#else
! in the seriell case, particles are only emitted on the current proc
sumOfMatchedParticles = mySumOfMatchedParticles
#endif

#if USE_MPI
IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
  ! add number of matching error to particle emission to fit
  ! number of added particles
  Species(FractNbr)%Init(iInit)%InsertedParticleMisMatch = nbrOfParticle  - sumOfMatchedParticles
  IF (nbrOfParticle .GT. sumOfMatchedParticles) THEN
    SWRITE(UNIT_StdOut,'(A)')'WARNING in ParticleEmission_parallel:'
    SWRITE(UNIT_StdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
    SWRITE(UNIT_StdOut,'(A,I0,A)')'matched only ', sumOfMatchedParticles, ' particles'
    SWRITE(UNIT_StdOut,'(A,I0,A)')'when ', NbrOfParticle, ' particles were required!'
  ELSE IF (nbrOfParticle .LT. sumOfMatchedParticles) THEN
        SWRITE(UNIT_StdOut,'(A)')'ERROR in ParticleEmission_parallel:'
        SWRITE(UNIT_StdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
        SWRITE(UNIT_StdOut,'(A,I0,A)')'matched ', sumOfMatchedParticles, ' particles'
        SWRITE(UNIT_StdOut,'(A,I0,A)')'when ', NbrOfParticle, ' particles were required!'
  ELSE IF (nbrOfParticle .EQ. sumOfMatchedParticles) THEN
  !  WRITE(UNIT_stdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
  !  WRITE(UNIT_stdOut,'(A,I0,A)')'ParticleEmission_parallel: matched all (',NbrOfParticle,') particles!'
  END IF
#if USE_MPI
END IF ! PartMPI%iProc.EQ.0
#endif
! Return the *local* NbrOfParticle so that the following Routines only fill in
! the values for the local particles
NbrOfParticle = mySumOfMatchedParticles

IF (chunkSize.GT.0) THEN
  DEALLOCATE(particle_positions, STAT=allocStat)
  IF (allocStat .NE. 0) &
    CALL ABORT(__STAMP__,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
END IF

END SUBROUTINE SetParticlePosition

SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle,init_or_sf)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PIC_Vars
USE MOD_Particle_Vars
USE MOD_Globals_Vars          ,ONLY: BoltzmannConst
USE MOD_Timedisc_Vars         ,ONLY: dt
USE MOD_Equation_Vars         ,ONLY: c,c2
USE MOD_PICInterpolation_vars ,ONLY: externalField
USE MOD_part_emission_tools   ,ONLY: CalcVelocity_maxwell_lpn,BessK,DEVI,SYNGE,QUASIREL,CalcVelocity_taylorgreenvortex
USE MOD_part_emission_tools   ,ONLY: CalcVelocity_emmert
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr,iInit,init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,j,PositionNbr
REAL                             :: Radius(3), n_vec(3), tan_vec(3), Velo1, Angle, Velo2, f
REAL                             :: Vec3D(3), RandVal(3), Vec1D
REAL                             :: II(3,3),JJ(3,3),NN(3,3)
INTEGER                          :: distnum,Rotation
REAL                             :: r1,r2,x_1,x_2,y_1,y_2,a,b,e,g,x_01,x_02,y_01,y_02, RandVal1
REAL                             :: Velosq, v_sum(3), v2_sum, maxwellfac
LOGICAL                          :: Is_ElemMacro
REAL                             :: sigma(3), ftl, PartVelo
REAL                             :: RandN_save
LOGICAL                          :: RandN_in_Mem
CHARACTER(30)                    :: velocityDistribution             ! specifying keyword for velocity distribution
REAL                             :: RadiusIC                         ! Radius for IC circle
REAL                             :: RadiusICGyro                     ! Radius for Gyrotron gyro radius
REAL                             :: NormalIC(3)                      ! Normal / Orientation of circle
REAL                             :: BasePointIC(3)                   ! base point for IC cuboid and IC sphere
REAL                             :: VeloIC                           ! velocity for inital Data
REAL                             :: VeloIC2                          ! square of velocity for inital Data
REAL                             :: VeloVecIC(3)                     ! normalized velocity vector
REAL                             :: WeibelVeloPar                    ! Parrallel velocity component for Weibel
REAL                             :: WeibelVeloPer                    ! Perpendicular velocity component for Weibel
REAL                             :: OneDTwoStreamVelo                ! Stream Velocity for the Two Stream Instability
REAL                             :: OneDTwoStreamTransRatio          ! Ratio between perpendicular and parallel velocity
REAL                             :: Alpha                            ! WaveNumber for sin-deviation initiation.
REAL                             :: MWTemperatureIC                  ! Temperature for Maxwell Distribution
REAL                             :: MJRatio(3)                       ! momentum to temperature ratio
! Maxwell-Juettner
REAL                             :: eps, anta, BesselK2,  gamm_k, max_val, qq, u_max, value, velabs, xixi, f_gamm
REAL                             :: VelocitySpread                         ! widening of init velocity
REAL                             :: vMag2                                  ! magnitude of velocity
!===================================================================================================================================

IF(NbrOfParticle.lt.1) RETURN
   IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PIC%maxParticleNumber!')
   END IF
RandN_in_Mem=.FALSE.
Is_ElemMacro = .FALSE.
SELECT CASE (init_or_sf)
CASE(1) !iInit
  IF (Species(FractNbr)%Init(iInit)%ElemVelocityICFileID.GT.0 .OR. Species(FractNbr)%Init(iInit)%ElemTemperatureFileID.GT.0) THEN
    Is_ElemMacro = .TRUE.
  END IF

  velocityDistribution=Species(FractNbr)%Init(iInit)%velocityDistribution
  VeloVecIC=Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
  VeloIC=Species(FractNbr)%Init(iInit)%VeloIC
  BasePointIC=Species(FractNbr)%Init(iInit)%BasePointIC(1:3)
  NormalIC=Species(FractNbr)%Init(iInit)%NormalIC(1:3)
  RadiusIC=Species(FractNbr)%Init(iInit)%RadiusIC
  Alpha=Species(FractNbr)%Init(iInit)%alpha
  RadiusICGyro=Species(FractNbr)%Init(iInit)%RadiusICGyro
  MWTemperatureIC=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  WeibelVeloPar=Species(FractNbr)%Init(iInit)%WeibelVeloPar
  WeibelVeloPer=Species(FractNbr)%Init(iInit)%WeibelVeloPer
  OneDTwoStreamVelo=Species(FractNbr)%Init(iInit)%OneDTwoStreamVelo
  OneDTwoStreamTransRatio=Species(FractNbr)%Init(iInit)%OneDTwoStreamTransRatio
  MJRatio(1)=Species(FractNbr)%Init(iInit)%MJxRatio
  MJRatio(2)=Species(FractNbr)%Init(iInit)%MJyRatio
  MJRatio(3)=Species(FractNbr)%Init(iInit)%MJzRatio
  SELECT CASE(TRIM(velocityDistribution))
  CASE('tangential_constant')
    Rotation       = Species(FractNbr)%Init(iInit)%Rotation
    VelocitySpread = Species(FractNbr)%Init(iInit)%VelocitySpread
    IF(VelocitySpread.GT.0)THEN
      IF(Species(FractNbr)%Init(iInit)%VelocitySpreadMethod.EQ.0)THEN
        ! sigma of normal Distribution, Kostas proposal
        VelocitySpread = VelocitySpread * VeloIC   !/(2.*SQRT(2.*LOG(10.)))
      ELSE IF(Species(FractNbr)%Init(iInit)%VelocitySpreadMethod.EQ.1)THEN
        ! sigma is defined by changing the width of the distribution function at 10% of its maxima
        ! the input value is the spread in percent, hence, 5% => v = v +- 0.05*v at 10% of maximum value
        ! width of the velocity spread, deltaV:
        VelocitySpread = 2.0*VelocitySpread * VeloIC
        ! computing the corresponding sigma
        VelocitySpread = VelocitySpread / (2.*SQRT(2.*LOG(10.)))
      ELSE
     CALL abort(&
__STAMP__&
,' This method for the velocity spread is not implemented.')
      END IF
      IF(alpha.GT.0) THEN
        vMag2 = (1.0+1./(alpha*alpha)) * VeloIC*VeloIC
      ELSE
        vMag2 = VeloIC*VeloIC
      END IF
    END IF
    VeloIC2        = VeloIC*VeloIC
  END SELECT
CASE(2) !SurfaceFlux
  IF (TRIM(Species(FractNbr)%Surfaceflux(iInit)%velocityDistribution).EQ.'constant') THEN
    velocityDistribution=Species(FractNbr)%Surfaceflux(iInit)%velocityDistribution
  ELSE
    CALL abort(&
__STAMP__&
,'only constant velo-distri implemented in SetParticleVelocity for surfaceflux!') !other distris in SetSurfacefluxVelocities!!!
  END IF
  VeloVecIC=Species(FractNbr)%Surfaceflux(iInit)%VeloVecIC(1:3)
  VeloIC=Species(FractNbr)%Surfaceflux(iInit)%VeloIC
  MWTemperatureIC=Species(FractNbr)%Surfaceflux(iInit)%MWTemperatureIC

CASE DEFAULT
  CALL abort(&
__STAMP__&
,'neither iInit nor Surfaceflux defined as reference!')
END SELECT

SELECT CASE(TRIM(velocityDistribution))
CASE('random')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .ne. 0) THEN
      CALL RANDOM_NUMBER(RandVal)
      RandVal(:) = RandVal(:) - 0.5
      RandVal(:) = RandVal(:)/SQRT(RandVal(1)**2+RandVal(2)**2+RandVal(3)**2)
      PartState(4:6,PositionNbr) = RandVal(1:3) * VeloIC
    END IF
    i = i + 1
  END DO
!  CASE('EOC_Test')
!     ! for leapfrog EOC test velo has to be set half dt before ICPos.
!     Radius = Species(1)%BasePointIC
!     IF (Species(1)%BasePointIC(1) > 0. ) THEN
!       n_vec = (/0.,0.,-1./)
!       Angle = 0.
!       Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!       Radius(1) = cos(Angle)
!       Radius(2) = sin(Angle)
!     ELSEIF (Species(1)%BasePointIC(3) > 0. ) THEN
!       n_vec = (/0.,1.,0./)
!       Angle = PI
!       Angle = Angle - PIC%GyrationFrequency * dt * 0.5
!       Radius(1) = sin(Angle)
!       Radius(3) = cos(Angle)
!     END IF
!     tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
!     tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
!     tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
!     PartState(4:6,1) = tan_vec(1:3) * Species(1)%VeloIC
CASE('constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        IF (Is_ElemMacro) THEN
          IF (Species(FractNbr)%Init(iInit)%ElemVelocityICFileID.GT.0) THEN
            PartState(4:6,PositionNbr) = Species(FractNbr)%Init(iInit)%ElemVelocityIC(1:3,PEM%Element(PositionNbr))
          ELSE
            PartState(4:6,PositionNbr) = VeloVecIC(1:3) * VeloIC
          END IF
        ELSE
          PartState(4:6,PositionNbr) = VeloVecIC(1:3) * VeloIC
        END IF
     END IF
     i = i + 1
  END DO
CASE('radial_constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        Radius(1:3) = PartState(1:3,PositionNbr) - BasePointIC(1:3)
        !  Unity radius
        !Radius(1:3) = Radius(1:3) / RadiusIC
        Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
        PartState(4:6,PositionNbr) = Radius(1:3) * VeloIC
     END IF
     i = i + 1
  END DO
CASE('tangential_constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        Radius(1:3) = PartState(1:3,PositionNbr) - BasePointIC(1:3)
        !  Normal Vector of circle
        n_vec(1:3) = NormalIC(1:3)
        ! If we're doing Leapfrog, then use velocities from half-timestep before
        IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
          Angle = 0.5 * dt * VeloIC / RadiusIC ! 0.5*dt*(v/r)
          JJ(1,1:3) = (/   0.,-n_vec(3), n_vec(2)/)
          JJ(2,1:3) = (/ n_vec(3),   0.,-n_vec(1)/)
          JJ(3,1:3) = (/-n_vec(2), n_vec(1),   0./)
          II(1,1:3) = (/1.,0.,0./)
          II(2,1:3) = (/0.,1.,0./)
          II(3,1:3) = (/0.,0.,1./)
          forall(j=1:3) NN(:,j) = n_vec(:)*n_vec(j)
          Radius = MATMUL( NN+cos(Angle)*(II-NN)+sin(Angle)*JJ , Radius )
        END IF
        !  Unity radius
        Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
        !  Vector Product rxn
        tan_vec(1) = Radius(2)*n_vec(3) - Radius(3)*n_vec(2)
        tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3)
        tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)

        IF(VelocitySpread.GT.0.)THEN
          IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
            Vec1D = RandN_save
            RandN_in_Mem=.FALSE.
          ELSE
            CALL RANDOM_NUMBER(RandVal)
            Velo1 = 2.0*RandVal(1)-1.0
            Velo2 = 2.0*RandVal(2)-1.0
            Velosq= Velo1**2+Velo2**2
            DO WHILE ((Velosq.LE.0).OR.(Velosq.GE.1))
              CALL RANDOM_NUMBER(RandVal)
              Velo1 = 2.0*RandVal(1)-1.0
              Velo2 = 2.0*RandVal(2)-1.0
              Velosq= Velo1**2+Velo2**2
            END DO
            Vec1D = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
            RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
            RandN_in_Mem=.TRUE.
          END IF
          ! velocity spread of tangential velocity
          IF(Rotation.EQ.1)THEN
            Vec3D  = tan_vec(1:3) * (VeloIC+Vec1D*VelocitySpread)
          ELSE
            Vec3D = -tan_vec(1:3) * (VeloIC+Vec1D*VelocitySpread)
          END IF
          ! compute axial velocity
          Vec1D = vMag2  - DOT_PRODUCT(Vec3D,Vec3D)
          IF(Vec1D.LT.0) CALL abort(&
__STAMP__&
,' Error in set velocity!',PositionNbr)
          Vec1D=SQRT(Vec1D)
          PartState(4:6,PositionNbr) = Vec3D+n_vec(1:3) * Vec1D
        ELSE ! no velocity spread
          ! If Gyrotron resonator: Add velocity in normal direction!
          IF (Alpha .gt. 0.) THEN
            n_vec = n_vec * ( 1 / Alpha )
          ELSE
            n_vec = 0
          END IF
          !  And finally the velocities
          IF(Rotation.EQ.1)THEN
            PartState(4:6,PositionNbr) = tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
          ELSE
            PartState(4:6,PositionNbr) = -tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
          END IF
        END IF
     END IF
     i = i + 1
  END DO

CASE('gyrotron_circle')
  i = 1
  IF (externalField(6).NE.0) THEN
    PIC%GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))
  ELSE
    PIC%GyroVecDirSIGN = -1
  END IF
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
     !! Position of particle on gyro circle changed in SetParticlePosition.F90: Problem
     !! We don't have the radius-vector any more. Thus transport the radius vector from there to here.
     ! Or do Alternative way: Hack the radius by intersecting two circles (big IC and small gyro circle)
       r1 = RadiusIC
       r2 = RadiusICGyro
       x_1 = 0.
       y_1 = 0.
       x_2 = PartState(1,PositionNbr)
       y_2 = PartState(2,PositionNbr)
       IF (x_1 .eq. x_2) THEN
         a = (x_1 - x_2)/(y_2-y_1)
         b = ((r1**2-r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
             & /(2.*y_2-2.*y_1)
         e = (a**2+1.)
         f = (2.*a*(b-y_1))-2.*x_1
         g = (b-y_1)**2-r1**2+x_1**2
         ! intersection points
         x_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
         x_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
         y_01 = x_01 * a + b
         y_02 = x_02 * a + b
       ELSE
         a = (y_1 - y_2)/(x_2-x_1)
         b = ((r1**2 - r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
              & /(2.*x_2 - 2. * x_1)
         e = (a**2 + 1.)
         f = 2. * a * (b - x_1) -2 *y_1
         g = (b-x_1)**2 - r1**2 + y_1**2
         y_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
         y_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
         x_01 = y_01 * a + b
         x_02 = y_02 * a + b
       END IF
       CALL RANDOM_NUMBER(RandVal1)
       IF (RandVal1 .ge. 0.5) THEN
         Radius(1) = PartState(1,PositionNbr) - x_01
         Radius(2) = PartState(2,PositionNbr) - y_01
       ELSE
         Radius(1) = PartState(1,PositionNbr) - x_02
         Radius(2) = PartState(2,PositionNbr) - y_02
       END IF

        Radius(3) = 0.
        !Check if Radius has correct length
        IF ((SQRT(Radius(1)**2+Radius(2)**2)-r1).ge.1E-15) THEN
          IPWRITE(UNIT_stdOut,*)"Error in setparticle velocity, gyrotron circle. &
                    & Radius too big after intersection."
        END IF
        !  Normal Vector of circle
        n_vec(1:3) = NormalIC(1:3)

        ! If we're doing Leapfrog, then use velocities from half-timestep before. This only applies in
        ! x- and y-direction. z has allways same velo.
!           IF (ParticlePushMethod.EQ.'boris_leap_frog_scheme') THEN
!             ! get angle of part on gyrocircle
!             Angle = ACOS(Radius(1)/Species(1)%RadiusICGyro)
!             IF (Radius(2).LE.0) THEN
!               Angle = 2*PI-Angle
!             END IF
!             ! shift position angle half dt back in time (as particle moves clockwise,
!             ! we add dalpha in ccw direction)
!             Angle = Angle + PIC%GyrationFrequency * dt * 0.5 * PIC%GyroVecDirSIGN
!             Radius(1) = cos(Angle)
!             Radius(2) = sin(Angle)
!           END IF
           !  Unity radius
           Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
           !  Vector Product rxn
           tan_vec(1) = Radius(2)*n_vec(3)*PIC%GyroVecDirSIGN - Radius(3)*n_vec(2)
           tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3) *PIC%GyroVecDirSIGN
           tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
           ! If Gyrotron resonator: Add velocity in normal direction!
           IF (Alpha .gt. 0.) THEN
             n_vec = n_vec * ( 1. / Alpha )
           ELSE
             n_vec = 0.
           END IF
           !  And finally the velocities
           PartState(4:6,PositionNbr) = (tan_vec(1:3) + n_vec(1:3)) * VeloIC
           IF (ABS(SQRT(PartState(4,PositionNbr)*PartState(4,PositionNbr) &
                      + PartState(5,PositionNbr)*PartState(5,PositionNbr))&
                      - VeloIC) .GT. 10.) THEN
             SWRITE(*,'(A,3(E21.14,X))') 'Velocity=', PartState(4:6,PositionNbr)
             CALL abort(&
__STAMP__&
,'ERROR in gyrotron_circle spaceIC!',PositionNbr)
           END If
           IF (PartState(4,PositionNbr).NE.PartState(4,PositionNbr) .OR. &
               PartState(5,PositionNbr).NE.PartState(5,PositionNbr) .OR. &
               PartState(6,PositionNbr).NE.PartState(6,PositionNbr)     ) THEN
             SWRITE(*,'(A,3(E21.14,X))') 'WARNING:! NaN: Velocity=', PartState(4:6,PositionNbr)
           END If
        END IF
        i = i + 1
     END DO
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
       IF (Is_ElemMacro) THEN
         CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit, Element=PEM%Element(PositionNbr))
       ELSE
         CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit)
       END IF
       PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
CASE('taylorgreenvortex')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
       CALL CalcVelocity_taylorgreenvortex(FractNbr, Vec3D, iInit=iInit, Element=PEM%Element(PositionNbr))
       PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
CASE('emmert')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      CALL CalcVelocity_emmert(FractNbr, iInit, Vec3D)
    END IF
    PartState(4:6,PositionNbr) = Vec3D(1:3)
  END DO
CASE('maxwell')
  v_sum(1:3) = 0.0
  v2_sum = 0.0

  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        DO distnum = 1, 3
!          IF (.NOT.DoZigguratSampling) THEN !polar method
            IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
              Vec3D(distnum) = RandN_save
              RandN_in_Mem=.FALSE.
            ELSE
              CALL RANDOM_NUMBER(RandVal)
              Velo1 = 2.0*RandVal(1)-1.0
              Velo2 = 2.0*RandVal(2)-1.0
              Velosq= Velo1**2+Velo2**2
              DO WHILE ((Velosq.LE.0).OR.(Velosq.GE.1))
                CALL RANDOM_NUMBER(RandVal)
                Velo1 = 2.0*RandVal(1)-1.0
                Velo2 = 2.0*RandVal(2)-1.0
                Velosq= Velo1**2+Velo2**2
              END DO
              Vec3D(distnum) = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_in_Mem=.TRUE.
            END IF
!          ELSE !ziggurat method
!            Vec3D(distnum)=rnor()
!          END IF
        END DO
        PartState(4:6,PositionNbr) = Vec3D(1:3)
        v_sum(1:3) = v_sum(1:3) + Vec3D(1:3)
        v2_sum = v2_sum + Vec3D(1)**2+Vec3D(2)**2+Vec3D(3)**2
     END IF
     i = i + 1
  END DO
  v_sum(1:3) = v_sum(1:3) / NbrOfParticle
  v2_sum = v2_sum / NbrOfParticle
  maxwellfac = SQRT(3. * BoltzmannConst * MWTemperatureIC / &              ! velocity of maximum
                 (Species(FractNbr)%MassIC*v2_sum))

  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
       PartState(4:6,PositionNbr) = (PartState(4:6,PositionNbr) - v_sum(1:3)) * maxwellfac &
                                    + VeloIC *VeloVecIC(1:3)
     END IF
     i = i + 1
  END DO

CASE('maxwell-juettner')
  xixi = Species(FractNbr)%MassIC*c2/ &
         (BoltzmannConst*MWTemperatureIC)
  BesselK2 = BessK(2,xixi)

  ! Find initial value for Newton Algorithm
  IF (xixi .LT. 4.d0) THEN
    gamm_k = 5.d0 * BoltzmannConst*MWTemperatureIC/ &
                    (Species(FractNbr)%MassIC*c2)
  ELSE
    gamm_k = 1.d0 + BoltzmannConst*MWTemperatureIC/ &
                    (Species(FractNbr)%MassIC*c2)
  END IF
  f_gamm = DEVI(Species(FractNbr)%MassIC, MWTemperatureIC, gamm_k)

  ! Newton Algorithm to find maximum value of distribution function
  ! (valid for both the relativistic and quasi relativistic distribution)
  i = 0
  eps=1e-8
  DO WHILE (abs(f_gamm) .GT. eps )
    i = i+1
    gamm_k = gamm_k - f_gamm/(xixi*(3._8*gamm_k**2._8-1._8)-10._8*gamm_k)
    f_gamm = DEVI(Species(FractNbr)%MassIC, MWTemperatureIC, gamm_k)
    IF(i.EQ.101) &
      CALL abort(&
__STAMP__&
,' Newton Algorithm to find maximum value of Maxwell-Juettner distribution has not been successfull!')
  END DO

  u_max = sqrt(1.d0-1.d0/(gamm_k*gamm_k))*c
  IF (xixi .LT. 692.5_8) THEN                  ! due to numerical precision
        max_val = SYNGE(u_max, MWTemperatureIC, &
                              Species(FractNbr)%MassIC, BesselK2)
      ELSE
        max_val = QUASIREL(u_max, MWTemperatureIC, &
                                 Species(FractNbr)%MassIC)
      END IF

  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    anta  = 1._8
    value = 0._8

    ! acception rejection method for velocity's absolute value
    DO WHILE (anta .GT. value)
      CALL RANDOM_NUMBER(velabs)
      CALL RANDOM_NUMBER(anta)
      velabs = velabs*c
      anta = anta*max_val
      IF (xixi .LT. 692.5_8) THEN
        value = SYNGE(velabs, MWTemperatureIC, &
                              Species(FractNbr)%MassIC, BesselK2)
      ELSE
        value = QUASIREL(velabs, MWTemperatureIC, &
                                 Species(FractNbr)%MassIC)
      END IF
    END DO

    ! polar method for velocity's x&y direction
    ! (required to generate elliptical random distribution)
    qq = 2._8
    DO WHILE ((qq .GT. 1._8) .OR. (qq .EQ. 0._8))
      CALL RANDOM_NUMBER(RandVal)
      RandVal = 2._8*RandVal-1._8
      qq = RandVal(1)*RandVal(1) + RandVal(2)*RandVal(2)
    END DO
    qq = sqrt(-2._8*log(qq)/qq)
    Vec3D(1) = RandVal(1)*qq*MJRatio(1)
    Vec3D(2) = RandVal(2)*qq*MJRatio(2)

    ! polar method for velocity's z direction
    qq = 2._8
    DO WHILE ((qq .GT. 1._8) .OR. (qq .EQ. 0._8))
      CALL RANDOM_NUMBER(RandVal)
      RandVal(:) = 2*RandVal(:)-1._8
      qq = RandVal(1)*RandVal(1) + RandVal(2)*RandVal(2)
    END DO
    qq = sqrt(-2._8*log(qq)/qq)
    Vec3D(3) = RandVal(1)*qq*MJRatio(3)

    Velosq  = sqrt(Vec3D(1)*Vec3D(1)+Vec3D(2)*Vec3D(2)+Vec3D(3)*Vec3D(3))
    PartState(4:6,PositionNbr) = velabs/Velosq*Vec3D
  END DO


CASE('weibel')
  v_sum(:)  = 0.0
  sigma(:) = 0.0

  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
!      IF (.NOT.DoZigguratSampling) THEN !polar method
        Velosq = 2.
        DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
          CALL RANDOM_NUMBER(RandVal)
          RandVal(:) = 2*RandVal(:)-1
          Velosq = RandVal(1)**2 + RandVal(2)**2
        END DO
        Velosq = sqrt(-2*log(Velosq)/Velosq)
        Vec3D(1) = RandVal(1)*Velosq
        Vec3D(2) = RandVal(2)*Velosq

        IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
          Vec3D(3) = RandN_save
          RandN_in_Mem=.FALSE.
        ELSE
          Velosq = 2.
          DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
            CALL RANDOM_NUMBER(RandVal)
            RandVal(:) = 2*RandVal(:)-1
            Velosq = RandVal(1)**2 + RandVal(2)**2
          END DO
          Velosq = sqrt(-2*log(Velosq)/Velosq)
          Vec3D(3) = RandVal(1)*Velosq
          RandN_save = RandVal(2)*Velosq
          RandN_in_Mem=.TRUE.
        END IF
!      ELSE !ziggurat method
!        Vec3D(1) = rnor()
!        Vec3D(2) = rnor()
!        Vec3D(3) = rnor()
!      END IF
      v_sum(:) = v_sum(:)  + Vec3D(:)
      sigma(:)   = sigma(:)    + Vec3D(:)**2
      PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
!  WRITE(*,*) PartVeloX(1), PartVeloY(1), PartVeloZ(1)

!    IF (NbrOfParticle .GT. 0) THEN
!      v_sum(:)  = 0.0
!      sigma(:) = 0.0

!      DO i=1,NbrOfParticle
!  !       elemNbr = PartToElem%Element(i)
!  v_sum(1)  = v_sum(1)  + PartVeloX(i)
!  v_sum(2)  = v_sum(2)  + PartVeloY(i)
!  v_sum(3)  = v_sum(3)  + PartVeloZ(i)
!  sigma(1) = sigma(1) + PartVeloX(i)**2
!  sigma(2) = sigma(2) + PartVeloY(i)**2
!  sigma(3) = sigma(3) + PartVeloZ(i)**2
!  !       NPart(elemNbr) = NPart(elemNbr) + 1
!      END DO


  IF (NbrOfParticle .GT. 1) THEN
    v_sum(:)  = v_sum(:)/NbrOfParticle
    sigma(:) = (NbrOfParticle/(NbrOfParticle-1))*(sigma(:)/NbrOfParticle-v_sum(:)**2)
        ! Verschiebungssatz der korrigierten Stichprobenkovarianz:
  ELSE ! s^2(X)=1/(N-1)(N*E(X^2) - N*E(X)^2)
    v_sum(:) = 0.
    sigma(:) = 1.
  END IF

  ftl = 0
  DO i=1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(4,PositionNbr)   = (PartState(4,PositionNbr)  -v_sum(1)) * &
                                    SQRT(WeibelVeloPar**2/sigma(1)) *c
      PartState(5:6,PositionNbr) = (PartState(5:6,PositionNbr)-v_sum(2:3)) * &
                                    SQRT(WeibelVeloPer**2/sigma(2:3)) *c
      PartVelo = SQRT(PartState(4,PositionNbr)**2+PartState(5,PositionNbr)**2+PartState(6,PositionNbr)**2)

      DO WHILE (PartVelo .GE. c)
        ftl = ftl+1
        IPWRITE(UNIT_stdOut,*) 'Number of Particles FTL:', ftl
!        IF (.NOT.DoZigguratSampling) THEN !polar method
          Velosq = 2.
          DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
            CALL RANDOM_NUMBER(RandVal)
            RandVal(:) = 2*RandVal(:)-1
            Velosq = RandVal(1)**2 + RandVal(2)**2
          END DO
          Velosq = sqrt(-2*log(Velosq)/Velosq)
          Vec3D(1) = RandVal(1)*Velosq
          Vec3D(2) = RandVal(2)*Velosq

          IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
            Vec3D(3) = RandN_save
            RandN_in_Mem=.FALSE.
          ELSE
            Velosq = 2.
            DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
              CALL RANDOM_NUMBER(RandVal)
              RandVal(:) = 2*RandVal(:)-1
              Velosq = RandVal(1)**2 + RandVal(2)**2
            END DO
            Velosq = sqrt(-2*log(Velosq)/Velosq)
            Vec3D(3) = RandVal(1)*Velosq
            RandN_save = RandVal(2)*Velosq
            RandN_in_Mem=.TRUE.
          END IF
!        ELSE !ziggurat method
!          Vec3D(1) = rnor()
!          Vec3D(2) = rnor()
!          Vec3D(3) = rnor()
!        END IF

        PartState(4:6,PositionNbr) = Vec3D(1:3)

        PartState(4,PositionNbr)   = (Vec3D(1)  -v_sum(1)) * &
            SQRT(WeibelVeloPar**2/sigma(1)) *c
        PartState(5:6,PositionNbr) = (Vec3D(2:3)-v_sum(2:3)) * &
            SQRT(WeibelVeloPer**2/sigma(2:3)) *c
        PartVelo = SQRT(PartState(4,PositionNbr)**2+PartState(5,PositionNbr)**2+PartState(6,PositionNbr)**2)
      END DO
    END IF
  END DO

CASE('OneD-twostreaminstabilty')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(4,PositionNbr) = OneDTwoStreamVelo
      PartState(5:6,PositionNbr) = OneDTwoStreamTransRatio
!      IF (.NOT.DoZigguratSampling) THEN !polar method
        Velosq = 2.
        DO WHILE ((Velosq .GT. 1.) .OR. (Velosq .EQ. 0.))
          CALL RANDOM_NUMBER(RandVal)
          RandVal(:) = 2*RandVal(:)-1
          Velosq = RandVal(1)**2 + RandVal(2)**2
        END DO
        RandVal(1:2) = RandVal(1:2)*sqrt(-2*log(Velosq)/Velosq)
!      ELSE
!        RandVal(1) = rnor()
!        RandVal(2) = rnor()
!      END IF
      PartState(5,PositionNbr) = RandVal(1)*OneDTwoStreamTransRatio* &
                                                   OneDTwoStreamVelo
      PartState(6,PositionNbr) = RandVal(2)*OneDTwoStreamTransRatio* &
                                                   OneDTwoStreamVelo
    END IF
  END DO

CASE('IMD') ! read IMD particle velocity from *.chkpt file -> velocity space has already been read when particles position was done
  ! do nothing
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')

END SELECT
END SUBROUTINE SetParticleVelocity


END  MODULE MOD_part_pos_and_velo
