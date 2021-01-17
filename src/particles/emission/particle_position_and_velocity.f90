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

INTERFACE AD_SetInitElectronVelo
  MODULE PROCEDURE AD_SetInitElectronVelo
END INTERFACE

!===================================================================================================================================
PUBLIC         :: SetParticleVelocity,SetParticlePosition, AD_SetInitElectronVelo
!===================================================================================================================================
CONTAINS


SUBROUTINE SetParticlePositionCellLocal(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species,Symmetry
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
  ! ParticleNumber and PartDensity not working together
  IF (NbrofParticle.EQ.0.AND.(Species(FractNbr)%Init(iInit)%ParticleNumber.EQ.0)) RETURN
  IF ((NbrofParticle.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.LE.0.)) THEN
    DoExactPartNumInsert =  .TRUE.
    IF(Symmetry%Axisymmetric) CALL abort(&
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
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PartState
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_part_emission_tools    ,ONLY: IntegerDivide,SetCellLocalParticlePosition,SetParticlePositionPoint
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionEquidistLine, SetParticlePositionLine, SetParticlePositionDisk
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionCuboidCylinder, SetParticlePositionGyrotronCircle,SetParticlePositionCircle
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionSphere, SetParticlePositionSinDeviation
USE MOD_part_emission_tools    ,ONLY: SetParticlePositionPhotonSEEDisc, SetParticlePositionPhotonCylinder
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
IF ((NbrOfParticle .LE. 0).AND. (ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.)) RETURN

DimSend  = 3                   !save (and send) only positions
nChunks  = 1                   ! Standard: Nicht-MPI
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
chunkSize = INT(nbrOfParticle/nChunks)
IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
  chunkSize = chunkSize*(1-nChunks) + nbrOfParticle
END IF
! all proc taking part in particle inserting
IF (PartMPI%InitGroup(InitGroup)%MPIROOT.OR.nChunks.GT.1) THEN
#endif
  ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
  IF (allocStat .NE. 0) &
    CALL abort(__STAMP__,'ERROR in SetParticlePosition: cannot allocate particle_positions!')

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
  CASE('photon_cylinder') ! cylinder case for photonionization
    CALL SetParticlePositionPhotonCylinder(FractNbr,iInit,chunkSize,particle_positions)
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
  IF (allocStat .NE. 0) &
    CALL ABORT(__STAMP__,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
END IF

END SUBROUTINE SetParticlePosition

SUBROUTINE SetParticleVelocity(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Determine the particle velocity of each inserted particle
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn, CalcVelocity_taylorgreenvortex, CalcVelocity_FromWorkFuncSEE
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_gyrotroncircle
USE MOD_Particle_Boundary_Vars  ,ONLY: DoBoundaryParticleOutputHDF5
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_FPFlow_Init             ,ONLY: FP_BuildTransGaussNums
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
REAL                            :: iRanPart(3, NbrOfParticle), Vec3D(3)
!===================================================================================================================================

IF(NbrOfParticle.LT.1) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PDM%maxParticleNumber!')
END IF

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
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit)
      PartState(4:6,PositionNbr) = Vec3D(1:3)
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
  CALL FP_BuildTransGaussNums(NbrOfParticle, iRanPart)
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
        IF(DoBoundaryParticleOutputHDF5) CALL StoreBoundaryParticleProperties(PositionNbr,FractNbr,PartState(1:3,PositionNbr),&
                                          UNITVECTOR(PartState(4:6,PositionNbr)),Species(FractNbr)%Init(iInit)%NormalIC,mode=2,&
                                          usevMPF_optIN=.FALSE.)
    END IF
  END DO
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'wrong velo-distri!')
END SELECT
END SUBROUTINE SetParticleVelocity

SUBROUTINE AD_SetInitElectronVelo(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
! Deletes all background gas particles and updates the particle index list
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_FPFlow_Init             ,ONLY: FP_BuildTransGaussNums
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: FractNbr,iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)           :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! LOCAL VARIABLES
INTEGER                         :: i, PositionNbr
CHARACTER(30)                   :: velocityDistribution
REAL                            :: VeloIC, VeloVecIC(3), maxwellfac, VeloVecNorm
REAL                            :: iRanPart(3, NbrOfParticle), Vec3D(3)
!===================================================================================================================================
IF(NbrOfParticle.LT.1) RETURN
IF(Species(FractNbr)%ChargeIC.LE.0.0) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PDM%maxParticleNumber!')
END IF

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
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloVecIC(1:3) * VeloIC 
    END IF
  END DO
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      CALL CalcVelocity_maxwell_lpn(DSMC%AmbiDiffElecSpec, Vec3D, Temperature=Species(FractNbr)%Init(iInit)%MWTemperatureIC)
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + Vec3D(1:3)
    END IF
  END DO
CASE('maxwell')
  CALL FP_BuildTransGaussNums(NbrOfParticle, iRanPart)
  maxwellfac = SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(DSMC%AmbiDiffElecSpec)%MassIC)
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr.GT.0) THEN
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + iRanPart(1:3,i)*maxwellfac
    END IF
  END DO
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'Velo-Distri not implemented for ambipolar diffusion!')
END SELECT

END SUBROUTINE AD_SetInitElectronVelo

END  MODULE MOD_part_pos_and_velo
