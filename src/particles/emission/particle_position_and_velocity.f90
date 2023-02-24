!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

INTERFACE SetPartPosAndVeloEmissionDistribution
  MODULE PROCEDURE SetPartPosAndVeloEmissionDistribution
END INTERFACE

!===================================================================================================================================
PUBLIC :: SetParticleVelocity, SetParticlePosition, SetPartPosAndVeloEmissionDistribution
!===================================================================================================================================
CONTAINS


SUBROUTINE SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position for processor-local particles (only in processor elements)
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_part_Emission_Tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
#endif /*USE_MPI*/
USE MOD_Symmetry_Vars          ,ONLY: Symmetry
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
! ParticleNumber and PartDensity not working together
IF (NbrofParticle.EQ.0.AND.(Species(FractNbr)%Init(iInit)%ParticleNumber.EQ.0)) RETURN
IF ((NbrofParticle.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.LE.0.)) THEN
  DoExactPartNumInsert =  .TRUE.
  IF(Symmetry%Axisymmetric) CALL abort(__STAMP__,'Axisymmetric: Particle insertion only possible with PartDensity!')
END IF
chunksize = 0
#if USE_MPI
! emission group communicator
InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1) THEN
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
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PartState,FractNbrOld,chunkSizeOld,NeutralizationBalance
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_part_emission_tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition,SetParticlePositionPoint
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionEquidistLine, SetParticlePositionLine, SetParticlePositionDisk
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionCuboidCylinder, SetParticlePositionGyrotronCircle,SetParticlePositionCircle
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionSphere, SetParticlePositionSinDeviation
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionPhotonSEEDisc, SetParticlePositionPhotonCylinder
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionPhotonSEERectangle, SetParticlePositionPhotonRectangle
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionPhotonHoneycomb, SetParticlePositionPhotonSEEHoneycomb
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionLandmark,SetParticlePositionLandmarkNeutralization
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionLiu2010Neutralization,SetParticlePositionLiu2010Neutralization3D
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionLiu2010SzaboNeutralization
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
INTEGER                                  :: DimSend
#if USE_MPI
INTEGER                                  :: InitGroup
#endif
!===================================================================================================================================
IF (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
  CALL SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
  RETURN
END IF

Species(FractNbr)%Init(iInit)%sumOfRequestedParticles = NbrOfParticle
IF((NbrOfParticle.LE.0).AND.(ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.)) RETURN

DimSend  = 3 ! save (and send) only positions
nChunks  = 1 ! Standard: Nicht-MPI
Species(FractNbr)%Init(iInit)%sumOfMatchedParticles   = 0
Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles = 0
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

! Set default chunkSize
chunkSize = INT(nbrOfParticle/nChunks)
IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
  chunkSize = chunkSize*(1-nChunks) + nbrOfParticle
END IF

#endif
! Set special chunkSize (also for MPI=OFF)
SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
CASE('2D_Liu2010_neutralization_Szabo','3D_Liu2010_neutralization_Szabo')
  ! Override the chunkSize with the processor-local sum of the required number of emitted particles
  chunkSize = NeutralizationBalance ! Sum over all elements of each processor (not global over all procs)
  nChunks   = 2 ! dummy value that is greater than 1
CASE DEFAULT
END SELECT
#if USE_MPI

! all proc taking part in particle inserting
IF (PartMPI%InitGroup(InitGroup)%MPIROOT.OR.nChunks.GT.1) THEN
#endif
  ! Allocate part pos buffer
  ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
  ! Sanity check
  IF (allocStat .NE. 0) CALL abort(__STAMP__,'ERROR in SetParticlePosition: cannot allocate particle_positions!')

  !------------------SpaceIC-cases: start-----------------------------------------------------------!
  SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
  CASE('point')
    CALL SetParticlePositionPoint(FractNbr,iInit,chunkSize,particle_positions)
  CASE('line_with_equidistant_distribution')
    CALL SetParticlePositionEquidistLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE('line')
    CALL SetParticlePositionLine(FractNbr,iInit,chunkSize,particle_positions)
  CASE('disc')
    CALL SetParticlePositionDisk(FractNbr,iInit,chunkSize,particle_positions)
  CASE('circle', 'circle_equidistant')
    CALL SetParticlePositionCircle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('gyrotron_circle')
    CALL SetParticlePositionGyrotronCircle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('cuboid','cylinder')
    CALL SetParticlePositionCuboidCylinder(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sphere')
    CALL SetParticlePositionSphere(FractNbr,iInit,chunkSize,particle_positions)
  CASE('sin_deviation')
    CALL SetParticlePositionSinDeviation(FractNbr,iInit,particle_positions)
  CASE('photon_SEE_disc') ! disc case for surface distribution
    CALL SetParticlePositionPhotonSEEDisc(FractNbr,iInit,chunkSize,particle_positions)
  CASE('photon_SEE_rectangle') ! rectangle case for surface distribution
    CALL SetParticlePositionPhotonSEERectangle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('photon_SEE_honeycomb') ! Honeycomb disc case for surface distribution
    CALL SetParticlePositionPhotonSEEHoneycomb(FractNbr,iInit,chunkSize,particle_positions)
  CASE('photon_cylinder') ! cylinder case for photonionization
    CALL SetParticlePositionPhotonCylinder(FractNbr,iInit,chunkSize,particle_positions)
  CASE('photon_rectangle') ! rectangle case for photonionization
    CALL SetParticlePositionPhotonRectangle(FractNbr,iInit,chunkSize,particle_positions)
  CASE('photon_honeycomb') ! Honeycomb case for photonionization
    CALL SetParticlePositionPhotonHoneycomb(FractNbr,iInit,chunkSize,particle_positions)
  CASE('2D_landmark','2D_landmark_copy')
    ! Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
    ! for low-temperature partially magnetized plasmas (2019)
    IF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'2D_landmark')THEN
      FractNbrOld  = FractNbr
      chunkSizeOld = chunkSize
      CALL SetParticlePositionLandmark(chunkSize,particle_positions,1)
    ELSEIF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'2D_landmark_copy')THEN
      CALL SetParticlePositionLandmark(chunkSizeOld,particle_positions,2)
    END IF ! TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'2D_landmark')
  CASE('2D_landmark_neutralization')
    ! Neutralization at const. x-position from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark
    ! for low-temperature partially magnetized plasmas (2019)
    CALL SetParticlePositionLandmarkNeutralization(chunkSize,particle_positions)
  CASE('2D_Liu2010_neutralization')
    ! Neutralization at right BC (max. x-position) H. Liu "Particle-in-cell simulation of a Hall thruster" (2010) - 2D case
    CALL SetParticlePositionLiu2010Neutralization(chunkSize,particle_positions)
  CASE('3D_Liu2010_neutralization')
    ! Neutralization at right BC (max. z-position) H. Liu "Particle-in-cell simulation of a Hall thruster" (2010) - 3D case
    CALL SetParticlePositionLiu2010Neutralization3D(FractNbr,iInit,chunkSize,particle_positions)
  CASE('2D_Liu2010_neutralization_Szabo','3D_Liu2010_neutralization_Szabo')
    ! Neutralization at right BC (max. x-position) H. Liu "Particle-in-cell simulation of a Hall thruster" (2010) - 2D and 3D case
    ! Some procs might have nothing to emit (cells are quasi neutral or negatively charged)
    IF(chunkSize.GT.0) CALL SetParticlePositionLiu2010SzaboNeutralization(chunkSize,particle_positions)
  END SELECT
  !------------------SpaceIC-cases: end-------------------------------------------------------------------------------------------
#if USE_MPI
ELSE !no mpi root, nchunks=1
  chunkSize = 0
END IF
! Need to open MPI communication regardless of the chunk number. Make it only dependent on the number of procs
IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1) THEN
  CALL SendEmissionParticlesToProcs(chunkSize,DimSend,FractNbr,iInit,Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles,particle_positions)
! Finish emission on local proc
ELSE
#endif /*USE_MPI*/
  Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles = 0
  ParticleIndexNbr = 1
  DO i = 1,chunkSize
    ! Find a free position in the PDM array
    IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
      ParticleIndexNbr = PDM%nextFreePosition(Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles + 1 + PDM%CurrentNextFreePosition)
    END IF
    IF (ParticleIndexNbr.NE.0) THEN
      PartState(1:DimSend,ParticleIndexNbr) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
      PDM%ParticleInside( ParticleIndexNbr) = .TRUE.
      CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.FALSE.)
      IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
        Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles + 1
        PDM%IsNewPart(ParticleIndexNbr)  = .TRUE.
        PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
      END IF
    ELSE
      CALL ABORT(__STAMP__,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO
#if USE_MPI
END IF

! Start communicating matched particles. This routine is finished in particle_emission.f90
CALL MPI_IREDUCE( Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles &
                , Species(FractNbr)%Init(iInit)%sumOfMatchedParticles   &
                , 1                                    &
                , MPI_INTEGER                          &
                , MPI_SUM                              &
                , 0                                    &
                , PartMPI%InitGroup(InitGroup)%COMM    &
                , PartMPI%InitGroup(InitGroup)%Request &
                , IERROR)
#else
! in the serial case, particles are only emitted on the current processor
Species(FractNbr)%Init(iInit)%sumOfMatchedParticles = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles
#endif /*USE_MPI*/

! Return the *local* NbrOfParticle so that the following Routines only fill in
! the values for the local particles
NbrOfParticle = Species(FractNbr)%Init(iInit)%mySumOfMatchedParticles

IF (chunkSize.GT.0) THEN
  DEALLOCATE(particle_positions, STAT=allocStat)
  IF (allocStat .NE. 0) CALL ABORT(__STAMP__,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
END IF

END SUBROUTINE SetParticlePosition


SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst,eV2Kelvin
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn, CalcVelocity_taylorgreenvortex, CalcVelocity_FromWorkFuncSEE
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_gyrotroncircle
USE MOD_Particle_Boundary_Vars  ,ONLY: DoBoundaryParticleOutputHDF5
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_part_tools              ,ONLY: BuildTransGaussNums, InRotRefFrameCheck
USE MOD_Particle_Vars           ,ONLY: CalcBulkElectronTemp,BulkElectronTemp
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: FractNbr,iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)           :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i, PositionNbr
CHARACTER(30)                   :: velocityDistribution
REAL                            :: VeloIC, VeloVecIC(3), maxwellfac, VeloVecNorm
REAL                            :: iRanPart(3, NbrOfParticle), Vec3D(3),MPF
!===================================================================================================================================

IF(NbrOfParticle.LT.1) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber) CALL abort(__STAMP__,'NbrOfParticle > PDM%maxParticleNumber! '//&
                                                                'Increase Part-maxParticleNumber or use more processors.')

velocityDistribution=Species(FractNbr)%Init(iInit)%velocityDistribution
VeloIC=Species(FractNbr)%Init(iInit)%VeloIC
VeloVecIC=Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
VeloVecNorm = VECNORM(VeloVecIC(1:3))
IF (VeloVecNorm.GT.0.0) THEN
  VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
END IF

SELECT CASE(TRIM(velocityDistribution))
CASE('constant')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      PartState(4:6,PositionNbr) = VeloVecIC(1:3) * VeloIC
    END IF
  END DO
CASE('gyrotron_circle')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_gyrotroncircle(FractNbr, Vec3D, iInit, PositionNbr)
      PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
CASE('maxwell_lpn','2D_landmark','2D_landmark_copy','2D_landmark_neutralization')
  ! maxwell_lpn: Maxwell low particle number
  ! 2D_landmark: Ionization profile from T. Charoy, 2D axial-azimuthal particle-in-cell benchmark for low-temperature partially
  !              magnetized plasmas (2019)
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit)
      PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
CASE('2D_Liu2010_neutralization','3D_Liu2010_neutralization','2D_Liu2010_neutralization_Szabo','3D_Liu2010_neutralization_Szabo')
  IF(.NOT.CalcBulkElectronTemp) CALL abort(__STAMP__,&
      'Velocity distribution 2D_Liu2010_neutralization requires CalcBulkElectronTemp=T')
  ! Use the global electron temperature if available
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, Temperature=BulkElectronTemp*eV2Kelvin)
      PartState(4:6,PositionNbr) = Vec3D(1:3)
      ! Mirror the x-velocity component (force all electrons to travel in -x direction)
      PartState(4,PositionNbr) = -ABS(PartState(4,PositionNbr))
    END IF
  END DO
CASE('taylorgreenvortex')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
       CALL CalcVelocity_taylorgreenvortex(FractNbr, Vec3D, iInit=iInit, Element=PEM%GlobalElemID(PositionNbr))
       PartState(4:6,PositionNbr) = Vec3D(1:3)
    END IF
  END DO
CASE('maxwell')
  CALL BuildTransGaussNums(NbrOfParticle, iRanPart)
  maxwellfac = SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC)
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
       PartState(4:6,PositionNbr) = VeloIC *VeloVecIC(1:3) + iRanPart(1:3,i)*maxwellfac
    END IF
  END DO
CASE('IMD') ! read IMD particle velocity from *.chkpt file -> velocity space has already been read when particles position was done
  ! do nothing
CASE('photon_SEE_energy')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
        CALL CalcVelocity_FromWorkFuncSEE(FractNbr, Vec3D, iInit=iInit)
        PartState(4:6,PositionNbr) = Vec3D(1:3)
        ! Store the particle information in PartStateBoundary.h5
        IF(DoBoundaryParticleOutputHDF5) THEN
          IF(usevMPF)THEN
            MPF = Species(FractNbr)%Init(iInit)%MacroParticleFactor ! Use emission-specific MPF
          ELSE
            MPF = Species(FractNbr)%MacroParticleFactor ! Use species MPF
          END IF ! usevMPF
          CALL StoreBoundaryParticleProperties(PositionNbr,FractNbr,PartState(1:3,PositionNbr),&
               UNITVECTOR(PartState(4:6,PositionNbr)),Species(FractNbr)%Init(iInit)%NormalIC,&
               iBC=Species(FractNbr)%Init(iInit)%PartBCIndex,mode=2,MPF_optIN=MPF)
        END IF ! DoBoundaryParticleOutputHDF5
    END IF
  END DO
CASE DEFAULT
  CALL abort(__STAMP__,'wrong velo-distri! velocityDistribution='//TRIM(velocityDistribution))
END SELECT

IF(UseRotRefFrame) THEN
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      PDM%InRotRefFrame(PositionNbr) = InRotRefFrameCheck(PositionNbr)
    END IF
  END DO
END IF

END SUBROUTINE SetParticleVelocity


!===================================================================================================================================
!> Initialize particle position and velocity from a distribution given by .h5 file
!> the .h5 file contains n, T, vr and vz for each species.
!> Only the electron temperature is currently used as the temperature of all heavy species is initialized with 300K.
!> Each processor creates all species randomly in each element using the sub-volumes defined by the Gaussian quadrature and
!> guarantees that at least one particle is created for each sub volume (depending of course on the MPF).
!===================================================================================================================================
SUBROUTINE SetPartPosAndVeloEmissionDistribution(iSpec,iInit,NbrOfParticle)
! modules
!USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals                ,ONLY: myrank,UNIT_StdOut,MPI_COMM_WORLD,abort
USE MOD_part_tools             ,ONLY: InitializeParticleMaxwell,InterpolateEmissionDistribution2D
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Vars          ,ONLY: Species, PDM, PartState, PEM, LastPartPos
USE MOD_Particle_Tracking      ,ONLY: ParticleInsideCheck
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Emission_Vars ,ONLY: EmissionDistributionDim, EmissionDistributionN
USE MOD_Interpolation          ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_Basis                  ,ONLY: BarycentricWeights
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Equation               ,ONLY: ExactFunc
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP,sJ
USE MOD_Interpolation_Vars     ,ONLY: NodeTypeVISU,NodeType
USE MOD_Eval_xyz               ,ONLY: TensorProductInterpolation
USE MOD_Mesh_Vars              ,ONLY: NGeo,XCL_NGeo,XiCL_NGeo,wBaryCL_NGeo,offsetElem
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemVolume_Shared,BoundsOfElem_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Dielectric_Vars        ,ONLY: DoDielectric,isDielectricElem,DielectricNoParticles
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: iSpec, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,iPart,GlobalElemID,PositionNbr
REAL              :: iRan, RandomPos(3), MPF
REAL              :: PartDens(1:3) ! dummy vector because the routine can only return vector values
LOGICAL           :: InsideFlag
REAL              :: IntegrationWeight
INTEGER           :: k,l,m
REAL              :: RandVal(3),Xi(3)
INTEGER           :: nPart,Nloc,Nred
REAL              :: nPartCell,nPartPerSubVol
REAL              :: BBoxVolume,SimPartDens

TYPE Interpolation
  REAL,ALLOCATABLE  :: densityVISU(:,:,:,:)
  REAL,ALLOCATABLE  :: J_NAnalyze(:,:,:,:)
  REAL,ALLOCATABLE  :: Coords_NAnalyze(:,:,:,:)
  REAL,ALLOCATABLE  :: xIP_VISU(:),wIP_VISU(:)
  REAL,ALLOCATABLE  :: Vdm_N_EQ_emission(:,:) ! < Vandermonde mapping from NodeType to equidistant (visu) node set
END TYPE Interpolation

TYPE(Interpolation),ALLOCATABLE :: NED(:)
!===================================================================================================================================/
! Sanity check
IF(PP_N.GT.EmissionDistributionN) CALL abort(__STAMP__,'PP_N > EmissionDistributionN')

! Allocate type for all polynomial degrees from PP_N to EmissionDistributionN
ALLOCATE(NED(PP_N:EmissionDistributionN))

DO Nloc = PP_N, EmissionDistributionN
  ! Allocate all arrays for Nloc
  ALLOCATE(NED(Nloc)%densityVISU(1,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NED(Nloc)%J_NAnalyze(1,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NED(Nloc)%Coords_NAnalyze(3,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(NED(Nloc)%xIP_VISU(0:Nloc))
  ALLOCATE(NED(Nloc)%wIP_VISU(0:Nloc))
  ALLOCATE(NED(Nloc)%Vdm_N_EQ_emission(0:Nloc,0:Nloc))

  ! Allocate and determine Vandermonde mapping from NodeType to equidistant (visu) node set
  CALL GetVandermonde(PP_N, NodeType, Nloc, NodeTypeVISU, NED(Nloc)%Vdm_N_EQ_emission, modal=.FALSE.)
  CALL GetNodesAndWeights(Nloc, NodeTypeVISU, NED(Nloc)%xIP_VISU, wIP=NED(Nloc)%wIP_VISU)
END DO ! Nloc = 1, EmissionDistributionN

NbrOfParticle = 0
MPF = Species(iSpec)%MacroParticleFactor

DO iElem = 1, nElems
  ! Do not emit particles in dielectrics
  IF(DoDielectric)THEN
    IF(DielectricNoParticles.AND.isDielectricElem(iElem)) CYCLE
  END IF ! DoDielectric

  GlobalElemID = iElem + offsetElem

  SELECT CASE(EmissionDistributionDim)
  CASE(1) ! 1D
    CALL abort(__STAMP__,'EmissionDistributionDim=1 is not implemented')
  CASE(2) ! 2D

    ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
    CALL ChangeBasis3D(3, PP_N, EmissionDistributionN, NED(EmissionDistributionN)%Vdm_N_EQ_emission, Elem_xGP(1:3,:,:,:,iElem), &
                                                       NED(EmissionDistributionN)%Coords_NAnalyze(1:3,:,:,:))

    ! Interpolate the Jacobian to the equidistant grid: be careful we interpolate the inverse of the inverse of the jacobian ;-)
    CALL ChangeBasis3D(1, PP_N, EmissionDistributionN, NED(EmissionDistributionN)%Vdm_N_EQ_emission, 1./sJ(:,:,:,iElem), &
                                                       NED(EmissionDistributionN)%J_NAnalyze(1:1,:,:,:))

    ! Integrate the density with the highest polynomial degree and calculate the average number of particles for the whole cell
    nPartCell = 0.
    NED(EmissionDistributionN)%densityVISU = -1.
    DO m=0,EmissionDistributionN
      DO l=0,EmissionDistributionN
        DO k=0,EmissionDistributionN
          ! Evaluate the equidistant solution from .h5 data file
          PartDens(1:3) = InterpolateEmissionDistribution2D(iSpec, iInit, NED(EmissionDistributionN)%Coords_NAnalyze(1:3,k,l,m),&
              dimLower=3, dimUpper=3, transformation=.FALSE.)

          ! Get density at interpolation point
          NED(EmissionDistributionN)%densityVISU(1,k,l,m) = PartDens(1)

          ! Integrate the Jacobian
          IntegrationWeight = NED(EmissionDistributionN)%wIP_VISU(k)*&
                              NED(EmissionDistributionN)%wIP_VISU(l)*&
                              NED(EmissionDistributionN)%wIP_VISU(m)*&
                              NED(EmissionDistributionN)%J_NAnalyze(1,k,l,m)

          ! Get total number of particles per cell
          nPartCell = nPartCell + NED(EmissionDistributionN)%densityVISU(1,k,l,m)*IntegrationWeight
        END DO ! k
      END DO ! l
    END DO ! m

    ! Apply MPF after summation
    nPartCell = nPartCell/REAL(MPF)

    ! Skip empty elements: Cut-off at arbitrarily small number
    IF(nPartCell.LT.1e-4) CYCLE

    ! Get number of particles per sub volume
    nPartPerSubVol = nPartCell/REAL((EmissionDistributionN+1)**3)

    ! Check if number of particles per sub cell drops below 1 (this would cause artificial noise in the density distribution)
    IF(nPartPerSubVol.LT.1.0)THEN
      ! Approximation of polynomial degree for 1 particle per sub volume, use lower value INT() to increase the number of particles
      Nred = INT(nPartCell**(1./3.) - 1.0)
    ELSE
      ! Use full emission polynomial
      Nred = EmissionDistributionN
    END IF

    ! If less than 1 particle per element is emitted, force ARM emission
    IF(nPartCell.LT.1.0) Nred = 0

    ! Check if the emission polynomial degree has to be reduced and if it falls below PP_N use cell-const. emission (via ARM)
    IF(Nred.GE.PP_N)THEN
      ! Emit all particles in the sub volumes

      ! Check if already calculated
      IF(Nred.NE.EmissionDistributionN)THEN
        ! Interpolate the physical position Elem_xGP to the analyze position, needed for exact function
        CALL ChangeBasis3D(3, PP_N, Nred, NED(Nred)%Vdm_N_EQ_emission, Elem_xGP(1:3,:,:,:,iElem), NED(Nred)%Coords_NAnalyze(1:3,:,:,:))

        ! Interpolate the Jacobian to the equidistant grid: be careful we interpolate the inverse of the inverse of the jacobian ;-)
        CALL ChangeBasis3D(1, PP_N, Nred, NED(Nred)%Vdm_N_EQ_emission, 1./sJ(:,:,:,iElem), NED(Nred)%J_NAnalyze(1:1,:,:,:))
      END IF ! Nred.NE.EmissionDistributionN

      ! Loop over all interpolation points
      NED(Nred)%densityVISU = -1.
      DO m=0,Nred
        DO l=0,Nred
          DO k=0,Nred
            ! Density field from .h5 data that is interpolated to the equidistant interpolation points (bilinear interpolation)
            PartDens(1:3) = InterpolateEmissionDistribution2D(iSpec,iInit,NED(Nred)%Coords_NAnalyze(1:3,k,l,m),dimLower=3,dimUpper=3,&
                transformation=.FALSE.)
            NED(Nred)%densityVISU(1,k,l,m) = PartDens(1)

            ! Integrate the Jacobian
            IntegrationWeight = NED(Nred)%wIP_VISU(k)*NED(Nred)%wIP_VISU(l)*NED(Nred)%wIP_VISU(m)*NED(Nred)%J_NAnalyze(1,k,l,m)

            ! Add noise via random number
            CALL RANDOM_NUMBER(iRan)
            nPart = INT(NED(Nred)%densityVISU(1,k,l,m)*IntegrationWeight/MPF + iRan)

            ! Loop over all newly created particles
            DO iPart = 1, nPart
              CALL RANDOM_NUMBER(RandVal)
              Xi(1) = -1.0 + SUM(NED(Nred)%wIP_VISU(0:k-1)) + NED(Nred)%wIP_VISU(k) * RandVal(1)
              Xi(2) = -1.0 + SUM(NED(Nred)%wIP_VISU(0:l-1)) + NED(Nred)%wIP_VISU(l) * RandVal(2)
              Xi(3) = -1.0 + SUM(NED(Nred)%wIP_VISU(0:m-1)) + NED(Nred)%wIP_VISU(m) * RandVal(3)
              IF(ANY(Xi.GT.1.0).OR.ANY(Xi.LT.-1.0))THEN
                IPWRITE(UNIT_StdOut,*) "Xi =", Xi
                CALL abort(__STAMP__,'xi out of range')
              END IF ! ANY(Xi.GT.1.0).OR.ANY(Xi.LT.-1.0)
              ! Get the physical coordinates that correspond to the reference coordinates
              CALL TensorProductInterpolation(Xi(1:3),3,NGeo,XiCL_NGeo,wBaryCL_NGeo,XCL_NGeo(1:3,0:NGeo,0:NGeo,0:NGeo,iElem),&
                  RandomPos(1:3)) !Map into phys. space

              !InsideFlag = .FALSE.
              !InsideFlag = ParticleInsideCheck(RandomPos,iPart,GlobalElemID)
              InsideFlag = .TRUE.

              ! Exclude particles outside of the element
              IF (InsideFlag) THEN
                NbrOfParticle = NbrOfParticle + 1
                PositionNbr                     = PDM%nextFreePosition(NbrOfParticle+PDM%CurrentNextFreePosition)
                PEM%GlobalElemID(PositionNbr)   = GlobalElemID
                PDM%ParticleInside(PositionNbr) = .TRUE.
                IF((PositionNbr.GE.PDM%maxParticleNumber).OR.&
                    (PositionNbr.EQ.0)) CALL abort(__STAMP__,'Emission: Increase maxParticleNumber!',PositionNbr)
                PartState(1:3,PositionNbr) = RandomPos(1:3)
                CALL InitializeParticleMaxwell(PositionNbr,iSpec,iElem,Mode=2,iInit=iInit)
              END IF
            END DO ! nPart
          END DO ! k
        END DO ! l
      END DO ! m

    ELSE ! Nred < PP_N
      ! Emit all particles in the element (ARM)
      ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,GlobalElemID) ) ! 1-2: Min, Max value; 1-3: x,y,z

        ! nPartCell is the REAL number of particles for the element with the correct volume of the element
        ! get number of particles for BB volume (outside particles are removed via ARM)
        BBoxVolume = (Bounds(2,3) - Bounds(1,3))*(Bounds(2,2) - Bounds(1,2))*(Bounds(2,1) - Bounds(1,1))
        ! Reconstruct density of simulation particles (without considering MPF)
        SimPartDens = nPartCell/ElemVolume_Shared(GetCNElemID(GlobalElemID))
        CALL RANDOM_NUMBER(iRan)
        nPart = INT(SimPartDens * BBoxVolume + iRan)
        !nPart=0

        DO iPart = 1, nPart
          InsideFlag = .FALSE. ! default

          ! Get random position in bounding box
          CALL RANDOM_NUMBER(RandomPos)
          RandomPos(1:3) = Bounds(1,1:3) + RandomPos(1:3)*(Bounds(2,1:3)-Bounds(1,1:3))

          ! Check if particle is inside the element (refmapping and tracing require a virtual particle with correct properties)
          SELECT CASE(TrackingMethod)
          CASE(REFMAPPING,TRACING)
            ! Attention: NbrOfParticle+PDM%CurrentNextFreePosition + 1
            PositionNbr                   = PDM%nextFreePosition(NbrOfParticle+PDM%CurrentNextFreePosition + 1)
            PEM%GlobalElemID(PositionNbr) = GlobalElemID
            LastPartPos(1:3,PositionNbr)  = RandomPos(1:3)
            InsideFlag                    = ParticleInsideCheck(RandomPos,PositionNbr,GlobalElemID)
          CASE(TRIATRACKING)
            InsideFlag  = ParticleInsideCheck(RandomPos,-1,GlobalElemID)
          CASE DEFAULT
            CALL abort(__STAMP__,'TrackingMethod not implemented! TrackingMethod =',IntInfoOpt=TrackingMethod)
          END SELECT

          ! Exclude particles outside of the element
          IF (InsideFlag) THEN
            NbrOfParticle = NbrOfParticle + 1
            PositionNbr                     = PDM%nextFreePosition(NbrOfParticle+PDM%CurrentNextFreePosition)
            PEM%GlobalElemID(PositionNbr)   = GlobalElemID
            PDM%ParticleInside(PositionNbr) = .TRUE.
            IF((PositionNbr.GE.PDM%maxParticleNumber).OR.&
               (PositionNbr.EQ.0)) CALL abort(__STAMP__,'Emission: Increase maxParticleNumber!',PositionNbr)
            PartState(1:3,PositionNbr) = RandomPos(1:3)
            CALL InitializeParticleMaxwell(PositionNbr,iSpec,iElem,Mode=2,iInit=iInit)
          END IF
        END DO ! nPart

      END ASSOCIATE

    END IF ! Nred.GE.PP_N

  CASE(3) ! 3D
    CALL abort(__STAMP__,'EmissionDistributionDim=3 is not implemented')
  END SELECT

END DO ! iElem = 1, nElems
END SUBROUTINE SetPartPosAndVeloEmissionDistribution


END  MODULE MOD_part_pos_and_velo
