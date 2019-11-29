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


#if USE_MPI
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle,mode)
#else
SUBROUTINE SetParticlePosition(FractNbr,iInit,NbrOfParticle)
#endif /*USE_MPI*/
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI,PartMPIInsert
USE MOD_Particle_Vars          ,ONLY: DoPoissonRounding,DoTimeDepInflow
#endif /*USE_MPI*/
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: BoltzmannConst
USE MOD_Particle_Vars          ,ONLY: IMDTimeScale,IMDLengthScale,IMDNumber,IMDCutOff,IMDCutOffxValue,IMDAtomFile
USE MOD_PIC_Vars
USE MOD_Particle_Vars          ,ONLY: Species,PDM,PartState,OutputVpiWarnings, Symmetry2DAxisymmetric
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Globals_Vars           ,ONLY: PI, TwoepsMach
USE MOD_Timedisc_Vars          ,ONLY: dt
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping, TriaTracking
USE MOD_Part_tools             ,ONLY: DICEUNITVECTOR
USE MOD_MacroBody_Vars         ,ONLY: UseMacroBody
USE MOD_MacroBody_tools        ,ONLY: INSIDEMACROBODY
USE MOD_PICInterpolation       ,ONLY: InterpolateVariableExternalField
USE MOD_PICInterpolation_Vars  ,ONLY: VariableExternalField
USE MOD_PICInterpolation_vars  ,ONLY: useVariableExternalField
USE MOD_Equation_vars          ,ONLY: c_inv
USE MOD_ReadInTools            ,ONLY: PrintOption
USE MOD_part_emission_tools    ,ONLY: IntegerDivide,CalcVelocity_maxwell_lpn,SamplePoissonDistri,SetCellLocalParticlePosition
USE MOD_part_emission_tools    ,ONLY: InsideExcludeRegionCheck
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
#if USE_MPI
INTEGER                                  :: mode
INTEGER                                  :: iProc,tProc, CellX, CellY, CellZ
INTEGER                                  :: msg_status(1:MPI_STATUS_SIZE)
INTEGER                                  :: MessageSize
LOGICAL                                  :: InsideMyBGM
#endif
REAL,ALLOCATABLE                         :: particle_positions(:)
INTEGER                                  :: allocStat
INTEGER                                  :: i,j,k,ParticleIndexNbr
INTEGER                                  :: mySumOfMatchedParticles, sumOfMatchedParticles
INTEGER                                  :: nChunks, chunkSize, chunkSize2
REAL                                     :: lineVector(3),VectorGap(3)
REAL                                     :: RandVal(3), Particle_pos(3),lineVector2(3)
REAL                                     :: n(3) , radius_vec(3)
REAL                                     :: II(3,3),JJ(3,3),NN(3,3)
REAL                                     :: RandVal1
REAL                                     :: radius, argumentTheta
REAL                                     :: rgyrate, Bintpol, pilen
REAL                                     :: x_step, y_step, z_step,  x_pos , y_pos
REAL                                     :: xlen, ylen, zlen
REAL                                     :: IMD_array(12),xMin,xMax,yMin,yMax,zMin,zMax
INTEGER                                  :: Nshift,ioUnit,io_error,IndNum
CHARACTER(LEN=255)                       :: StrTmp
INTEGER                                  :: iPart
REAL                                     :: Vec3D(3), l_ins, v_line, delta_l, v_drift_line, A_ins, PartIns
REAL                                     :: v_drift_BV(2), lrel_ins_BV(4), BV_lengths(2), v_BV(2), delta_lBV(2)
REAL                                     :: intersecPoint(3), orifice_delta, lPeri, ParaCheck(3)
INTEGER                                  :: DimSend, orificePeriodic
LOGICAL                                  :: orificePeriodicLog(2), insideExcludeRegion
LOGICAL                                  :: DoExactPartNumInsert
#if USE_MPI
INTEGER                                  :: InitGroup,nChunksTemp,mySumOfRemovedParticles
INTEGER,ALLOCATABLE                      :: PartFoundInProc(:,:) ! 1 proc id, 2 local part id
REAL,ALLOCATABLE                         :: ProcMeshVol(:)
INTEGER,ALLOCATABLE                      :: ProcNbrOfParticle(:)
#endif
!===================================================================================================================================
! emission group communicator
#if USE_MPI
InitGroup=Species(FractNbr)%Init(iInit)%InitCOMM
IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
  NbrofParticle=0
  RETURN
END IF
#endif /*USE_MPI*/

IF (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'cell_local') THEN
  DoExactPartNumInsert =  .FALSE.
  ! check if particle inserting during simulation or initial inserting and also if via partdensity or exact particle number
  ! nbrOfParticles is set for initial inserting if initialPartNum or partdensity is set in ini
  ! ParticleEmission and Partdensity not working together
  IF (NbrofParticle.EQ.0.AND.(Species(FractNbr)%Init(iInit)%ParticleEmission.EQ.0)) RETURN
  IF ((NbrofParticle.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.LE.0.)) THEN
    DoExactPartNumInsert = .TRUE.
    IF(Symmetry2DAxisymmetric) THEN
      CALL abort(&
__STAMP__&
,'Axisymmetric: Particle insertion only possible with PartDensity!')
    END IF
  END IF
  !IF ((Species(FractNbr)%Init(iInit)%ParticleEmission.GT.0).AND.(Species(FractNbr)%Init(iInit)%PartDensity.GT.0.)) CALL abort(&
!__STAMP__&
!,'ParticleEmission>0 and PartDensity>0. Can not be set at the same time for cell_local inserting. Set both for species: ',FractNbr)
  chunksize = 0
#if USE_MPI
  IF (mode.EQ.2) RETURN
  IF (PartMPI%InitGroup(InitGroup)%nProcs.GT.1 .AND. Species(FractNbr)%Init(iInit)%ElemPartDensityFileID.EQ.0) THEN
    IF (DoExactPartNumInsert) THEN
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
      CALL MPI_GATHER(GEO%LocalVolume,1,MPI_DOUBLE_PRECISION &
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
  !------------------SpaceIC-case: cell_local-------------------------------------------------------------------------------------
  IF ((chunksize.GT.0).OR.(Species(FractNbr)%Init(iInit)%PartDensity.GT.0.)) THEN
    CALL SetCellLocalParticlePosition(chunkSize,FractNbr,iInit,DoExactPartNumInsert)
  END IF
  NbrOfParticle = chunksize
  RETURN
END IF


PartIns=0.
lineVector = 0.0
A_ins = 0.0
l_ins = 0.0
lrel_ins_BV = 0.0
BV_lengths = 0.0
Particle_pos = 0.0
orificePeriodic = 0
orificePeriodicLog = .FALSE.
IF(Species(FractNbr)%Init(iInit)%VirtPreInsert) THEN ! (SpaceIC.EQ.'cuboid_vpi').OR.(SpaceIC.EQ.'cylinder_vpi')
  DimSend=6 !save (and send) velocities and positions
  ! the following is here, and not inside the select case, as it is needed to be excecuted just once and, most importantly, it
  ! calculates the virt. insertion height defining the virt. NbrOfParticle (which is chunked next) when PartDensity is used
  ! (-> could also be moved to particle_init?)
  lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
  lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
  lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
  IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
    CALL abort(&
__STAMP__&
,'BaseVectors are parallel!')
  ELSE
    A_ins = SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + lineVector(3) * lineVector(3))
    lineVector = lineVector / A_ins
  END IF
  v_drift_line = Species(FractNbr)%Init(iInit)%VeloIC * &
    ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*lineVector(1) + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*lineVector(2) &
    + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*lineVector(3) ) !lineVector component of drift-velocity
  l_ins=dt*RKdtFrac * ( v_drift_line + Species(FractNbr)%Init(iInit)%NSigma & !virt. insertion height
    * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) )
  IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'freestream') .OR. &
      (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') ) THEN
    BV_lengths(1) = SQRT(Species(FractNbr)%Init(iInit)%BaseVector1IC(1)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector1IC(2)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector1IC(3)**2)
    v_drift_BV(1) = Species(FractNbr)%Init(iInit)%VeloIC / BV_lengths(1) * &
      ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3) ) !BV1 component of drift-velocity
    BV_lengths(2) = SQRT(Species(FractNbr)%Init(iInit)%BaseVector2IC(1)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector2IC(2)**2 + &
                         Species(FractNbr)%Init(iInit)%BaseVector2IC(3)**2)
    v_drift_BV(2) = Species(FractNbr)%Init(iInit)%VeloIC / BV_lengths(2) * &
      ( Species(FractNbr)%Init(iInit)%VeloVecIC(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
      + Species(FractNbr)%Init(iInit)%VeloVecIC(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3) ) !BV2 component of drift-velocity
    IF ( .NOT.ALMOSTEQUAL(A_ins,BV_lengths(1)*BV_lengths(2)) ) THEN
      SWRITE(*,'(A72,2(x,I0),A1)') 'cross product and product of theirs lenghts for BaseVectors of Spec/Init',&
        FractNbr, iInit, ':'
      SWRITE(*,*) A_ins, BV_lengths(1)*BV_lengths(2)
      CALL abort(&
__STAMP__&
,' BaseVectors of the current SpaceIC are not parallel?')
    END IF
    IF ( .NOT.ALMOSTEQUAL(SQRT(v_drift_BV(1)**2+v_drift_BV(2)**2+v_drift_line**2),ABS(Species(FractNbr)%Init(iInit)%VeloIC)) ) THEN
      SWRITE(*,'(A60,2(x,I0),A1)') 'v_drift_BV1, v_drift_BV2, v_drift_line, VeloIC for Spec/Init',&
        FractNbr, iInit, ':'
      SWRITE(*,*) v_drift_BV(1),v_drift_BV(2),v_drift_line,Species(FractNbr)%Init(iInit)%VeloIC
      CALL abort(&
__STAMP__&
,' Something is wrong with the Basis of the current SpaceIC!')
    END IF
    IF ( TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'cuboid_vpi' ) THEN
      lrel_ins_BV=dt*RKdtFrac*( (/v_drift_BV(1),-v_drift_BV(1),v_drift_BV(2),-v_drift_BV(2)/)+Species(FractNbr)%Init(iInit)%NSigma &
        * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) )!rel. virt. insertion height:
      lrel_ins_BV(1:2)=lrel_ins_BV(1:2)/BV_lengths(1)                                                        !... in -BV1/+BV1 dir.
      lrel_ins_BV(3:4)=lrel_ins_BV(3:4)/BV_lengths(2)                                                        !... in -BV2/+BV2 dir.
      DO i=1,4
        IF (.NOT.Species(FractNbr)%Init(iInit)%vpiBVBuffer(i)) lrel_ins_BV(i)=0.
      END DO
    ELSE IF ( TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'orifice' ) THEN !cylinder-orifice
      lrel_ins_BV(1:4)=dt*RKdtFrac * ( SQRT(v_drift_BV(1)**2+v_drift_BV(2)**2) + Species(FractNbr)%Init(iInit)%NSigma &
        * SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(FractNbr)%MassIC) ) &
        / Species(FractNbr)%Init(iInit)%RadiusIC !rel. virt. insertion height is a single, maximum value for whole circumference
    END IF
    IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') .AND. &
        (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC) .EQ. 'cuboid_vpi') ) THEN !needs further devel.!
      SELECT CASE (GEO%nPeriodicVectors)
      CASE (0)
      CASE (1)
        ParaCheck(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * GEO%PeriodicVectors(3,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * GEO%PeriodicVectors(2,1)
        ParaCheck(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * GEO%PeriodicVectors(1,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * GEO%PeriodicVectors(3,1)
        ParaCheck(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * GEO%PeriodicVectors(2,1) - &
          Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * GEO%PeriodicVectors(1,1)
        IF ( .NOT.(SQRT(ParaCheck(1)**2+ParaCheck(2)**2+ParaCheck(3)**2).GT.TwoepsMach) ) orificePeriodic=1 !parallel with BV1
        IF (orificePeriodic .EQ. 0) THEN
          ParaCheck(1) = Species(FractNbr)%Init(iInit)%BaseVector2IC(2) * GEO%PeriodicVectors(3,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(3) * GEO%PeriodicVectors(2,1)
          ParaCheck(2) = Species(FractNbr)%Init(iInit)%BaseVector2IC(3) * GEO%PeriodicVectors(1,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(1) * GEO%PeriodicVectors(3,1)
          ParaCheck(3) = Species(FractNbr)%Init(iInit)%BaseVector2IC(1) * GEO%PeriodicVectors(2,1) - &
            Species(FractNbr)%Init(iInit)%BaseVector2IC(2) * GEO%PeriodicVectors(1,1)
          IF ( .NOT.(SQRT(ParaCheck(1)**2+ParaCheck(2)**2+ParaCheck(3)**2).GT.TwoepsMach) ) THEN
            orificePeriodic=2 !parallel with BV2
          ELSE
            CALL abort(&
__STAMP__&
,' PeriodicVector is not parallel to any orifice BV -> not implemented yet!')
          END IF
        ELSE IF (orificePeriodic .EQ. 1) THEN
          !PerVec cannot be parallel with BV2, as BV1 already is
        ELSE
          CALL abort(&
__STAMP__&
,' Something is wrong with the PeriodicVector and the orifice BV!')
        END IF
        lPeri=SQRT(GEO%PeriodicVectors(1,1)**2+GEO%PeriodicVectors(2,1)**2+GEO%PeriodicVectors(3,1)**2)
        IF ( .NOT.ALMOSTEQUAL(lPeri,BV_lengths(orificePeriodic)) ) THEN
          SWRITE(*,'(A22,I1,x,A1)') 'lPeri and length of BV',orificePeriodic,':'
          SWRITE(*,'(G0,x,G0)') lPeri,BV_lengths(orificePeriodic)
          CALL abort(&
__STAMP__&
,' PeriodicVector and its parallel BV ar not of same length! ')
        END IF
        orificePeriodicLog(1)=(orificePeriodic.EQ.1)
        orificePeriodicLog(2)=(orificePeriodic.EQ.2)
      CASE DEFAULT
        CALL abort(&
__STAMP__&
,' orifice region only implemented for 0 or 1 PeriodicVector!')
      END SELECT
    END IF !cuboid-orifice
  END IF !freestream or orifice
  !--calculation of (virtual) NbrOfParticle from virt. insertion height
  IF(Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) THEN
    SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    CASE ('cylinder_vpi')
      IF( TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType) .EQ. 'orifice' ) THEN
        A_ins = PI * ( Species(FractNbr)%Init(iInit)%RadiusIC*(1.+lrel_ins_BV(1)) )**2
      ELSE
        A_ins = PI * ( Species(FractNbr)%Init(iInit)%RadiusIC**2 - Species(FractNbr)%Init(iInit)%Radius2IC**2 )
      END IF
    CASE ('cuboid_vpi')
      IF( (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'freestream') .OR. &
          (TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType).EQ.'orifice') ) THEN
        A_ins = BV_lengths(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) * BV_lengths(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4))
      END IF
    CASE DEFAULT
      CALL abort(&
__STAMP__&
,'wrong SpaceIC for virtual Pre-Inserting region!')
    END SELECT
    PartIns = Species(FractNbr)%Init(iInit)%PartDensity * l_ins * A_ins / (Species(FractNbr)%MacroParticleFactor)
    IF(PartIns.GT.0.) THEN
      NbrOfParticle = INT(PartIns)
    ELSE
      NbrOfParticle = 0
    END IF
  END IF
ELSE
  DimSend=3 !save (and send) only positions
END IF

IF ( (NbrOfParticle .LE. 0).AND.(PartIns .LE. 0.).AND. (ABS(Species(FractNbr)%Init(iInit)%PartDensity).LE.0.) ) &
  RETURN !0<Partins<1: statistical handling of exact REAL-INT-conv. below!

nChunks = 1                   ! Standard: Nicht-MPI
sumOfMatchedParticles = 0
mySumOfMatchedParticles = 0

chunkSize = nbrOfParticle

! process myRank=0 generates the complete list of random positions for all emitted particles
#if USE_MPI
IF(( (nbrOfParticle.GT.PartMPI%InitGroup(InitGroup)%nProcs*10                             ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'circle_equidistant'                 ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'sin_deviation'                      ) .AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'cuboid_with_equidistant_distribution').AND.  &
     (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).NE.'line_with_equidistant_distribution' )))THEN
   nChunks = PartMPI%InitGroup(InitGroup)%nProcs
ELSE
   nChunks = 1
END IF

! communication
IF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle') nChunks=1

IF (mode.EQ.1) THEN
  chunkSize = INT(nbrOfParticle/nChunks)
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
    IF( Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) ) THEN
      ! statistical handling of exact REAL-INT-conversion -> values in send(1)- and receive(2)-mode might differ for VPI+PartDens
      ! (NbrOf Particle can differ from root to other procs and, thus, need to be communicated of calculated later again)
      CALL RANDOM_NUMBER(RandVal1)
      IF (.NOT.DoPoissonRounding .AND. .NOT.DoTimeDepInflow) THEN
        NbrOfParticle = INT(PartIns + RandVal1)
      ELSE IF (EXP(-PartIns).LE.TINY(PartIns) .AND.DoPoissonRounding) THEN
        IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
        NbrOfParticle = INT(PartIns + RandVal1)
        DoPoissonRounding = .FALSE.
      ELSE IF (DoPoissonRounding) THEN
        !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
        CALL SamplePoissonDistri( PartIns , NbrOfParticle , DoPoissonRounding)
      ELSE ! DoTimeDepInflow
        NbrOfParticle = INT(PartIns + RandVal1)
      END IF
    END IF
    chunkSize = chunkSize + ( nbrOfParticle - (nChunks*chunkSize) )
  END IF
  IF (PartMPI%InitGroup(InitGroup)%MPIROOT .OR. nChunks.GT.1) THEN
#endif
    ALLOCATE( particle_positions(1:chunkSize*DimSend), STAT=allocStat )
    IF (allocStat .NE. 0) THEN
      CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition: cannot allocate particle_positions!')
    END IF

    chunkSize2=chunkSize !will be changed during insertion for:
                         !  1.: vpi with PartDensity (orig. chunksize is for buffer region)
                         !  2.: excludeRegions (orig. chunksize is for SpaceIC without taking excludeRegions into account)
    !------------------SpaceIC-cases: start-----------------------------------------------------------!
    SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    !------------------SpaceIC-case: point------------------------------------------------------------------------------------------
    CASE ('point')
       Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC
       DO i=1,chunkSize
          particle_positions(i*3-2) = Particle_pos(1)
          particle_positions(i*3-1) = Particle_pos(2)
          particle_positions(i*3  ) = Particle_pos(3)
       END DO
    !------------------SpaceIC-case: line_with_equidistant_distribution-------------------------------------------------------------
    CASE ('line_with_equidistant_distribution')
      IF(NbrOfParticle.EQ.1)THEN
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + 0.5 * Species(FractNbr)%Init(iInit)%BaseVector1IC
      ELSE
        VectorGap = Species(FractNbr)%Init(iInit)%BaseVector1IC/(NbrOfParticle-1)
        DO i=1,chunkSize
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (i-1)*VectorGap
          particle_positions(i*3-2) = Particle_pos(1)
          particle_positions(i*3-1) = Particle_pos(2)
          particle_positions(i*3  ) = Particle_pos(3)
        END DO
      END IF
    !------------------SpaceIC-case: line-------------------------------------------------------------------------------------------
    CASE ('line')
      DO i=1,chunkSize
        CALL RANDOM_NUMBER(RandVal1)
        Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC*RandVal1
        particle_positions(i*3-2) = Particle_pos(1)
        particle_positions(i*3-1) = Particle_pos(2)
        particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: disc-------------------------------------------------------------------------------------------
    CASE('disc')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
        lineVector(1) = 1.0
        lineVector(2) = 1.0
        lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                         Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
        IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
          lineVector(1) = 1.0
          lineVector(3) = 1.0
          lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                            Species(FractNbr)%Init(iInit)%NormalIC(2)
        ELSE
          IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
            lineVector(2) = 1.0
            lineVector(3) = 1.0
            lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                 Species(FractNbr)%Init(iInit)%NormalIC(1)
          ELSE
            CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC Vektor darf nicht Nullvektor sein')
          END IF
        END IF
      END IF

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      DO i=1,chunkSize
         radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
         DO WHILE(radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC)
            CALL RANDOM_NUMBER(RandVal)
            RandVal = RandVal * 2. - 1.
            Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
                     (RandVal(1) * lineVector + RandVal(2) *lineVector2)

            radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
                           (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
                           (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
                           (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
                           (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
                           (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
         END DO
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    CASE('circle')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
         lineVector(1) = 1.0
         lineVector(2) = 1.0
         lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                           Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
         IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
            lineVector(1) = 1.0
            lineVector(3) = 1.0
            lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                              Species(FractNbr)%Init(iInit)%NormalIC(2)
         ELSE
            IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
               lineVector(2) = 1.0
               lineVector(3) = 1.0
               lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                 Species(FractNbr)%Init(iInit)%NormalIC(1)
            ELSE
              CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC should not be zero')
            END IF
         END IF
      END IF

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      radius = Species(FractNbr)%Init(iInit)%RadiusIC
      DO i=1,chunkSize
         CALL RANDOM_NUMBER(RandVal1)
         argumentTheta = 2.*pi*RandVal1
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                        linevector * cos(argumentTheta) * radius +  &
                        linevector2 * sin(argumentTheta) * radius
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: gyrotron_circle--------------------------------------------------------------------------------
    CASE('gyrotron_circle')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
         lineVector(1) = 1.0
         lineVector(2) = 1.0
         lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                           Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
         IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
            lineVector(1) = 1.0
            lineVector(3) = 1.0
            lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                              Species(FractNbr)%Init(iInit)%NormalIC(2)
         ELSE
            IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
               lineVector(2) = 1.0
               lineVector(3) = 1.0
               lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                 Species(FractNbr)%Init(iInit)%NormalIC(1)
            ELSE
              CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC should not be zero')
            END IF
         END IF
      END IF

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      radius = Species(FractNbr)%Init(iInit)%RadiusIC
      DO i=1,chunkSize
         CALL RANDOM_NUMBER(RandVal1)
         argumentTheta = 2.*pi*RandVal1
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                        linevector * cos(argumentTheta) * radius +  &
                        linevector2 * sin(argumentTheta) * radius
         ! Change position of particle on the small gyro circle
         ! take normal vecotr of the circle
         n(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
         ! generate radius vector (later it will be multiplied by the length of the
         ! gyro circles. For now we just need the vector)
         radius_vec(1) = Particle_pos(1) - Species(FractNbr)%Init(iInit)%BasePointIC(1)
         radius_vec(2) = Particle_pos(2) - Species(FractNbr)%Init(iInit)%BasePointIC(2)
         radius_vec(3) = Particle_pos(3) - Species(FractNbr)%Init(iInit)%BasePointIC(3)
         !rotate radius vector with random angle
         CALL RANDOM_NUMBER(RandVal1)
         argumentTheta=2.*pi*RandVal1
         JJ(1,1:3) = (/   0.,-n(3), n(2)/)
         JJ(2,1:3) = (/ n(3),   0.,-n(1)/)
         JJ(3,1:3) = (/-n(2), n(1),   0./)
         II(1,1:3) = (/1.,0.,0./)
         II(2,1:3) = (/0.,1.,0./)
         II(3,1:3) = (/0.,0.,1./)
         forall(j=1:3) NN(:,j) = n(:)*n(j)

!        1. determine the z-position in order to get the interpolated curved B-field
         CALL RANDOM_NUMBER(RandVal1)
         IF (NbrOfParticle.EQ.Species(FractNbr)%Init(iInit)%initialParticleNumber) THEN
           particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                           + RandVal1 * Species(FractNbr)%Init(iInit)%CuboidHeightIC
         ELSE
           particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                           + RandVal1 * dt*RKdtFrac &
                                           * Species(FractNbr)%Init(iInit)%VeloIC/Species(FractNbr)%Init(iInit)%alpha
         END IF

!        2. calculate curved B-field at z-position in order to determine size of gyro radius
         IF (useVariableExternalField) THEN
            IF(particle_positions(i*3).LT.VariableExternalField(1,1))THEN ! assume particles travel in positive z-direction
              CALL abort(&
__STAMP__&
,'SetParticlePosition: particle_positions(i*3) cannot be smaller than VariableExternalField(1,1). Fix *.csv data or emission!')
            END IF
            Bintpol = InterpolateVariableExternalField(particle_positions(i*3))
            rgyrate = 1./ SQRT ( 1 - (Species(FractNbr)%Init(iInit)%VeloIC**2 * (1 + 1./Species(FractNbr)%Init(iInit)%alpha**2)) &
                                * c_inv * c_inv ) * Species(FractNbr)%MassIC * Species(FractNbr)%Init(iInit)%VeloIC / &
                      ( Bintpol * abs( Species(FractNbr)%ChargeIC) )
         ELSE
           rgyrate =  Species(FractNbr)%Init(iInit)%RadiusICGyro
         END IF

         radius_vec = MATMUL( NN+cos(argumentTheta)*(II-NN)+sin(argumentTheta)*JJ , radius_vec )
         radius_vec(1:3) = radius_vec(1:3) / SQRT(radius_vec(1)**2+radius_vec(2)**2+radius_vec(3)**2) &
                       * rgyrate !Species(1)%RadiusICGyro
         ! Set new particles position:
         particle_positions(i*3-2) = Particle_pos(1) + radius_vec(1)
         particle_positions(i*3-1) = Particle_pos(2) + radius_vec(2)
         !particle_positions(i*3  )=0.
      END DO
    !------------------SpaceIC-case: circle_equidistant-----------------------------------------------------------------------------
    CASE('circle_equidistant')
      IF (Species(FractNbr)%Init(iInit)%NormalIC(3).NE.0) THEN
         lineVector(1) = 1.0
         lineVector(2) = 1.0
         lineVector(3) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(2))/ &
                           Species(FractNbr)%Init(iInit)%NormalIC(3)
      ELSE
         IF (Species(FractNbr)%Init(iInit)%NormalIC(2).NE.0) THEN
            lineVector(1) = 1.0
            lineVector(3) = 1.0
            lineVector(2) = -(Species(FractNbr)%Init(iInit)%NormalIC(1)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                              Species(FractNbr)%Init(iInit)%NormalIC(2)
         ELSE
            IF (Species(FractNbr)%Init(iInit)%NormalIC(1).NE.0) THEN
               lineVector(2) = 1.0
               lineVector(3) = 1.0
               lineVector(1) = -(Species(FractNbr)%Init(iInit)%NormalIC(2)+Species(FractNbr)%Init(iInit)%NormalIC(3))/ &
                                 Species(FractNbr)%Init(iInit)%NormalIC(1)
            ELSE
              CALL abort(&
__STAMP__&
,'Error in SetParticlePosition, NormalIC should not be zero')
            END IF
         END IF
      END IF

      lineVector2(1) = Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(3) - &
           Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(2)
      lineVector2(2) = Species(FractNbr)%Init(iInit)%NormalIC(3) * lineVector(1) - &
           Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(3)
      lineVector2(3) = Species(FractNbr)%Init(iInit)%NormalIC(1) * lineVector(2) - &
           Species(FractNbr)%Init(iInit)%NormalIC(2) * lineVector(1)

      lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * &
           lineVector(2) + lineVector(3) * lineVector(3))

      lineVector2 = lineVector2 / SQRT(lineVector2(1) * lineVector2(1) + lineVector2(2) * &
           lineVector2(2) + lineVector2(3) * lineVector2(3))

      radius = Species(FractNbr)%Init(iInit)%RadiusIC
      DO i=1,chunkSize
         argumentTheta = 2.*pi*i/chunkSize
         Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                        linevector * cos(argumentTheta) * radius +  &
                        linevector2 * sin(argumentTheta) * radius
         particle_positions(i*3-2) = Particle_pos(1)
         particle_positions(i*3-1) = Particle_pos(2)
         particle_positions(i*3  ) = Particle_pos(3)
      END DO
    !------------------SpaceIC-case: cuboid or cylinder-----------------------------------------------------------------------------
    CASE('cuboid','cylinder')
      lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
      lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
      lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
        Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
      IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
        CALL abort(&
__STAMP__&
,'BaseVectors are parallel!')
      ELSE
        lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
          lineVector(3) * lineVector(3))
      END IF
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
        SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
        CASE ('cuboid')
          CALL RANDOM_NUMBER(RandVal)
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
          IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
            Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
          ELSE
            Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3)
          END IF
        CASE ('cylinder')
          radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
          DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
             CALL RANDOM_NUMBER(RandVal)
             Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2-1) &
                          + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)
             radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                            Particle_pos(2) * Particle_pos(2) + &
                            Particle_pos(3) * Particle_pos(3) )
          END DO
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
          IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
            Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
          ELSE
            Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
          END IF
        END SELECT
        IF (UseMacroBody) THEN
          IF (INSIDEMACROBODY(Particle_pos)) THEN
            i=i+1
            CYCLE !particle is inside MacroParticle
          END IF
        END IF
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
        particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
        particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: sphere-----------------------------------------------------------------------------------------
    CASE('sphere')
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
        CALL RANDOM_NUMBER(RandVal1)
        radius = Species(FractNbr)%Init(iInit)%RadiusIC*RandVal1**(1./3.)
        Particle_pos = DICEUNITVECTOR()*radius + Species(FractNbr)%Init(iInit)%BasePointIC
        IF (UseMacroBody) THEN
          IF (INSIDEMACROBODY(Particle_pos)) THEN
            i=i+1
            CYCLE !particle is inside MacroParticle
          END IF
        END IF
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
        particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
        particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: cuboid_vpi-------------------------------------------------------------------------------------
    CASE('cuboid_vpi')
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
        ! Check if particle would reach comp. domain in one timestep
        CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit)
        CALL RANDOM_NUMBER(RandVal)
        v_line = Vec3D(1)*lineVector(1) + Vec3D(2)*lineVector(2) + Vec3D(3)*lineVector(3) !lineVector component of velocity
        delta_l = dt*RKdtFrac * v_line - l_ins * RandVal(3)
        IF (delta_l .LT. 0.) THEN
          IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
          CYCLE !particle would not reach comp. domain -> try new velo
        END IF
        SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType))
        CASE('perpendicular_extrusion')
          ! set particle positions depending on SpaceIC
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1) &
                                                                   + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2) &
                                                                   + lineVector * delta_l
        CASE('freestream')
          v_BV(1) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3)
          v_BV(1) = v_BV(1) / BV_lengths(1) !BV1 component of velocity
          v_BV(2) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
          v_BV(2) = v_BV(2) / BV_lengths(2) !BV2 component of velocity
          delta_lBV = dt*RKdtFrac * v_BV
          delta_lBV(1) = delta_lBV(1) + BV_lengths(1) * ( RandVal(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) - lrel_ins_BV(1) )
          delta_lBV(2) = delta_lBV(2) + BV_lengths(2) * ( RandVal(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4)) - lrel_ins_BV(3) )
          IF ( (delta_lBV(1).LT.0.) .OR. (delta_lBV(1).GT.BV_lengths(1)) ) THEN
            IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
            CYCLE !particle would not reach comp. domain in direction of BaseVector1 -> try new velo
          END IF
          IF ( (delta_lBV(2).LT.0.) .OR. (delta_lBV(2).GT.BV_lengths(2)) ) THEN
            IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
            CYCLE !particle would not reach comp. domain in direction of BaseVector2 -> try new velo
          END IF
          ! set particle positions depending on SpaceIC
          Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC &
                       + Species(FractNbr)%Init(iInit)%BaseVector1IC/BV_lengths(1) * delta_lBV(1) &
                       + Species(FractNbr)%Init(iInit)%BaseVector2IC/BV_lengths(2) * delta_lBV(2) &
                       + lineVector * delta_l
        CASE('orifice')
          ! set particle position (to be tried) depending on SpaceIC
          v_BV(1) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3)
          v_BV(1) = v_BV(1) / BV_lengths(1) !BV1 component of velocity
          v_BV(2) = Vec3D(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + Vec3D(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + Vec3D(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
          v_BV(2) = v_BV(2) / BV_lengths(2) !BV2 component of velocity
          delta_lBV = dt*RKdtFrac * v_BV
          delta_lBV(1) = delta_lBV(1) + BV_lengths(1) * ( RandVal(1)*(lrel_ins_BV(1)+1.0+lrel_ins_BV(2)) - lrel_ins_BV(1) )
          delta_lBV(2) = delta_lBV(2) + BV_lengths(2) * ( RandVal(2)*(lrel_ins_BV(3)+1.0+lrel_ins_BV(4)) - lrel_ins_BV(3) )
          Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC/BV_lengths(1) * delta_lBV(1) &
                       + Species(FractNbr)%Init(iInit)%BaseVector2IC/BV_lengths(2) * delta_lBV(2) &
                       + lineVector * delta_l
          IntersecPoint = Particle_pos - Vec3D * delta_l/v_line !Vector from BP to Intersec point of virt. path with orifice plane
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3))/BV_lengths(1)
          IF ( orificePeriodicLog(1) ) THEN
            IF ( (delta_lBV(1).LT.0.) .OR. (delta_lBV(1).GT.BV_lengths(1)) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. domain in direction of BaseVector1 -> try new velo
            END IF
          ELSE
            IF ( (orifice_delta.GT.BV_lengths(1)) .OR. (orifice_delta.LT.0.) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          END IF
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3))/BV_lengths(2)
          IF ( orificePeriodicLog(2) ) THEN
            IF ( (delta_lBV(2).LT.0.) .OR. (delta_lBV(2).GT.BV_lengths(2)) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. domain in direction of BaseVector2 -> try new velo
            END IF
          ELSE
            IF ( (orifice_delta.GT.BV_lengths(2)) .OR. (orifice_delta.LT.0.) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          END IF
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
        CASE DEFAULT
          CALL abort(&
__STAMP__&
,'wrong vpiDomainType for virtual Pre-Inserting region!')
        END SELECT
        IF (UseMacroBody) THEN
          IF (INSIDEMACROBODY(Particle_pos)) THEN
            i=i+1
            CYCLE !particle is inside MacroParticle
          END IF
        END IF
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        !store determined values or go to next particle
        particle_positions((chunkSize2+1)*6-5) = Particle_pos(1)
        particle_positions((chunkSize2+1)*6-4) = Particle_pos(2)
        particle_positions((chunkSize2+1)*6-3) = Particle_pos(3)
        particle_positions((chunkSize2+1)*6-2) = Vec3D(1)
        particle_positions((chunkSize2+1)*6-1) = Vec3D(2)
        particle_positions((chunkSize2+1)*6  ) = Vec3D(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: cylinder_vpi-----------------------------------------------------------------------------------
    CASE('cylinder_vpi')
      i=1
      chunkSize2=0
      DO WHILE (i .LE. chunkSize)
        ! Check if particle would reach comp. domain in one timestep
        CALL CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit=iInit)
        CALL RANDOM_NUMBER(RandVal)
        v_line = Vec3D(1)*lineVector(1) + Vec3D(2)*lineVector(2) + Vec3D(3)*lineVector(3) !lineVector component of velocity
        delta_l = dt*RKdtFrac * v_line - l_ins * RandVal(3)
        IF (delta_l .LT. 0.) THEN
          IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
          CYCLE !particle would not reach comp. domain -> try new velo
        END IF
        SELECT CASE(TRIM(Species(FractNbr)%Init(iInit)%vpiDomainType))
        CASE('perpendicular_extrusion')
          ! set particle positions depending on SpaceIC
          radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
          DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2-1)
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
          END DO
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC + lineVector * delta_l
        CASE('orifice')
          ! set particle position (to be tried) depending on SpaceIC
          radius = Species(FractNbr)%Init(iInit)%RadiusIC * (1.+lrel_ins_BV(1)) + 1. !lrel_ins_BV(1)=lrel_ins_BV(2)=...
          DO WHILE ( radius .GT. Species(FractNbr)%Init(iInit)%RadiusIC * (1.+lrel_ins_BV(1)) )
            CALL RANDOM_NUMBER(RandVal)
            Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC*(1.+lrel_ins_BV(1)) * (RandVal(1)*2-1) &
                         + Species(FractNbr)%Init(iInit)%BaseVector2IC*(1.+lrel_ins_BV(2)) * (RandVal(2)*2-1) !BV_lengths(:)=R
            radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                           Particle_pos(2) * Particle_pos(2) + &
                           Particle_pos(3) * Particle_pos(3) )
          END DO
          Particle_pos = Particle_pos + dt*RKdtFrac*Vec3D - lineVector*(dt*RKdtFrac*v_line-delta_l) !get old RandVal(3) for l_ins*R3

          IntersecPoint = Particle_pos - Vec3D * delta_l/v_line !Vector from BP to Intersec point of virt. path with orifice plane
          orifice_delta = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector1IC(1) &
                         + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector1IC(2) &
                         + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector1IC(3))/BV_lengths(1)
          radius = (IntersecPoint(1)*Species(FractNbr)%Init(iInit)%BaseVector2IC(1) &
                  + IntersecPoint(2)*Species(FractNbr)%Init(iInit)%BaseVector2IC(2) &
                  + IntersecPoint(3)*Species(FractNbr)%Init(iInit)%BaseVector2IC(3))/BV_lengths(2)
          radius = SQRT( orifice_delta*orifice_delta + radius*radius )
            IF ( radius.GT.BV_lengths(1) ) THEN
              IF (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) i=i+1
              CYCLE !particle would not reach comp. through orifice -> try new velo
            END IF
          Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
        CASE DEFAULT
          CALL abort(&
__STAMP__&
,'wrong vpiDomainType for virtual Pre-Inserting region!')
        END SELECT
        IF (UseMacroBody) THEN
          IF (INSIDEMACROBODY(Particle_pos)) THEN
            i=i+1
            CYCLE !particle is inside MacroParticle
          END IF
        END IF
        IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
          CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
          IF (insideExcludeRegion) THEN
            i=i+1
            CYCLE !particle is in excluded region
          END IF
        END IF
        !store determined values or go to next particle
        particle_positions((chunkSize2+1)*6-5) = Particle_pos(1)
        particle_positions((chunkSize2+1)*6-4) = Particle_pos(2)
        particle_positions((chunkSize2+1)*6-3) = Particle_pos(3)
        particle_positions((chunkSize2+1)*6-2) = Vec3D(1)
        particle_positions((chunkSize2+1)*6-1) = Vec3D(2)
        particle_positions((chunkSize2+1)*6  ) = Vec3D(3)
        i=i+1
        chunkSize2=chunkSize2+1
      END DO
    !------------------SpaceIC-case: cuboid_equal-----------------------------------------------------------------------------------
    CASE('cuboid_equal')
#if USE_MPI
      IF (PartMPI%InitGroup(InitGroup)%nProcs.GT. 1) THEN
        SWRITE(UNIT_stdOut,*)'WARNING in SetParticlePosition:'
        SWRITE(UNIT_stdOut,*)'cannot fully handle Particle Initial Condition \"cuboid equal\"'
        SWRITE(UNIT_stdOut,*)'in parallel mode (with more than one CPU)!'
        SWRITE(UNIT_stdOut,*)'USE WITH CARE!!!'
      END IF
      j=0
      mySumOfMatchedParticles = 0
      DO i=1,PDM%ParticleVecLength
         j=j+1
         ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
         IF (ParticleIndexNbr .ne. 0) THEN
            PartState(ParticleIndexNbr,1:3) = PartState(j,1:3)
            CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.FALSE.)
            IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
               mySumOfMatchedParticles = mySumOfMatchedParticles + 1
            ELSE
               PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
            END IF
            IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
              PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
              PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
            END IF
         ELSE
           CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition: ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
         END IF
      END DO
      CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER &
          , MPI_SUM, PartMPI%InitGroup(InitGroup)%COMM, IERROR)
      nbrOfParticle = NbrOfParticle - sumOfMatchedParticles
      IF (nbrOfParticle .NE. 0) THEN
        IPWRITE(UNIT_stdOut,*)'ERROR in ParticleEmission_parallel:'
        IPWRITE(UNIT_stdOut,'(I4,A,I8,A)')'matched ', sumOfMatchedParticles, ' particles'
        IPWRITE(UNIT_stdOut,'(I4,A,I8,A)')'when ', NbrOfParticle+sumOfMatchedParticles, ' particles were required!'
        CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel')
      END IF
      NbrOfParticle = mySumOfMatchedParticles
      DEALLOCATE( particle_positions, STAT=allocStat )
      IF (allocStat .NE. 0) THEN
        CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
      END IF
      RETURN
#else
      DO i=1,chunkSize
        ParticleIndexNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
        particle_positions(i*3-2 : i*3) = PartState(ParticleIndexNbr-Species(FractNbr)%Init(iInit)%initialParticleNumber,1:3)
      END DO
#endif
    !------------------SpaceIC-case: cuboid_with_equidistant_distribution-----------------------------------------------------------
    CASE ('cuboid_with_equidistant_distribution')
       IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
            (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
            * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(&
__STAMP__&
,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector1IC(1)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector1IC(2)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector1IC(3)**2 )
       ylen = SQRT(Species(FractNbr)%Init(iInit)%BaseVector2IC(1)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector2IC(2)**2 &
            + Species(FractNbr)%Init(iInit)%BaseVector2IC(3)**2 )
       zlen = ABS(Species(FractNbr)%Init(iInit)%CuboidHeightIC)

       ! make sure the vectors correspond to x,y,z-dir
       IF ((xlen.NE.Species(FractNbr)%Init(iInit)%BaseVector1IC(1)).OR. &
          (ylen.NE.Species(FractNbr)%Init(iInit)%BaseVector2IC(2)).OR. &
          (zlen.NE.Species(FractNbr)%Init(iInit)%CuboidHeightIC)) THEN
         CALL abort(&
__STAMP__&
,'Basevectors1IC,-2IC and CuboidHeightIC have to be in x,y,z-direction, respectively for emission condition')
        END IF
       x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
       y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
       z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
       iPart = 1
       DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
         x_pos = (i-0.5) * x_step + Species(FractNbr)%Init(iInit)%BasePointIC(1)
         DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
           y_pos =  Species(FractNbr)%Init(iInit)%BasePointIC(2) + (j-0.5) * y_step
           DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
             particle_positions(iPart*3-2) = x_pos
             particle_positions(iPart*3-1) = y_pos
             particle_positions(iPart*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                  + (k-0.5) * z_step
             iPart = iPart + 1
           END DO
         END DO
       END DO
    !------------------SpaceIC-case: sin_deviation----------------------------------------------------------------------------------
    CASE('sin_deviation')
       IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
            (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
            * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
         SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
         CALL abort(&
         __STAMP__&
         ,'ERROR: Number of particles in init / emission region',iInit)
       END IF
       xlen = abs(GEO%xmaxglob  - GEO%xminglob)
       ylen = abs(GEO%ymaxglob  - GEO%yminglob)
       zlen = abs(GEO%zmaxglob  - GEO%zminglob)
       pilen=2.0*PI/xlen
       x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
       y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
       z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
       iPart = 1
       DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
          x_pos = (i * x_step - x_step*0.5)
          x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Init(iInit)%Amplitude &
                  * sin(Species(FractNbr)%Init(iInit)%WaveNumber * pilen * x_pos)
          DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
            y_pos =  GEO%yminglob + j * y_step - y_step * 0.5
            DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
              particle_positions(iPart*3-2) = x_pos
              particle_positions(iPart*3-1) = y_pos
              particle_positions(iPart*3  ) = GEO%zminglob &
                                        + k * z_step - z_step * 0.5
              iPart = iPart + 1
            END DO
          END DO
       END DO
    !------------------SpaceIC-case: IMD--------------------------------------------------------------------------------------------
    CASE('IMD') ! read IMD particle position from *.chkpt file
      ! set velocity distribution to read external data
      SWRITE(UNIT_stdOut,'(A,A)') " Reading IMD atom data from file (IMDAtomFile): ",TRIM(IMDAtomFile)
      IF(TRIM(IMDAtomFile).NE.'no file found')THEN
        Species(FractNbr)%Init(iInit)%velocityDistribution='IMD'
#if USE_MPI
        IF(.NOT.PartMPI%InitGroup(InitGroup)%MPIROOT)THEN
          CALL abort(__STAMP__&
          ,'ERROR: Cannot SetParticlePosition in multi-core environment for SpaceIC=IMD!')
        END IF
#endif /*USE_MPI*/
        ! Read particle data from file
        ioUnit=GETFREEUNIT()
        OPEN(UNIT=ioUnit,FILE=TRIM(IMDAtomFile),STATUS='OLD',ACTION='READ',IOSTAT=io_error)
        IF(io_error.NE.0)THEN
          CALL abort(__STAMP__&
          ,'ERROR in particle_emission.f90: Cannot open specified File (particle position) for SpaceIC=IMD!')
        END IF
        ! IMD Data Format (ASCII)
        !   1      2    3         4           5         6         7         8         9         10       11     12
        !#C number type mass      x           y         z         vx        vy        vz        Epot     Z      eam_rho
        !   2294   0    26.981538 3589.254381 46.066405 91.985804 -1.576543 -0.168184 -0.163417 0.000000 2.4332 0.000000
        IndNum=INDEX(IMDAtomFile, '/',BACK = .TRUE.)
        IF(IndNum.GT.0)THEN
          !IndNum=INDEX(IMDAtomFile,'/',BACK = .TRUE.) ! get path without binary
          StrTmp=TRIM(IMDAtomFile(IndNum+1:LEN(IMDAtomFile)))
          IndNum=INDEX(StrTmp,'.',BACK = .TRUE.)
          IF(IndNum.GT.0)THEN
            StrTmp=StrTmp(1:IndNum-1)
            IndNum=INDEX(StrTmp,'.')
            IF(IndNum.GT.0)THEN
              StrTmp=StrTmp(IndNum+1:LEN(StrTmp))
            END IF
          END IF
        END IF
        read(StrTmp,*,iostat=io_error)  IMDNumber
        CALL PrintOption('IMD *.chkpt file','OUTPUT',StrOpt=StrTmp)
        CALL PrintOption('IMDNumber','OUTPUT',IntOpt=IMDNumber)
        Nshift=0
        xMin=HUGE(1.)
        yMin=HUGE(1.)
        zMin=HUGE(1.)
        xMax=-HUGE(1.)
        yMax=-HUGE(1.)
        zMax=-HUGE(1.)
        DO i=1,9
          READ(ioUnit,'(A)',IOSTAT=io_error)StrTmp
          IF(io_error.NE.0)THEN
             SWRITE(UNIT_stdOut,'(A,I5,A3,A)') 'Error in line ',i,' : ',TRIM(StrTmp)
          END IF
        END DO
        DO i=1,chunkSize
          READ(ioUnit,*,IOSTAT=io_error) IMD_array(1:12)
          IF(io_error>0)THEN
            CALL abort(__STAMP__&
            ,'ERROR in particle_emission.f90: Error reading specified File (particle position) for SpaceIC=IMD!')
          ELSE IF(io_error<0)THEN
            EXIT
          ELSE
            IF(1.EQ.2)THEN ! transformation
              ! 0.) multiply by unit system factor (1e-10)
              ! 1.) switch X and Z axis and invert Z
              ! 2.) shift origin in X- and Y-direction by -10nm
              Particle_pos = (/  IMD_array(6)       *1.E-10-10.13E-9,&
                                 IMD_array(5)       *1.E-10-10.13E-9,&
                               -(IMD_array(4)-10500)*1.E-10/)
            ELSE ! no transformation
              Particle_pos = (/  IMD_array(4)*IMDLengthScale,&
                                 IMD_array(5)*IMDLengthScale,&
                                 IMD_array(6)*IMDLengthScale/)
            END IF
            particle_positions((i-Nshift)*3-2) = Particle_pos(1)
            particle_positions((i-Nshift)*3-1) = Particle_pos(2)
            particle_positions((i-Nshift)*3  ) = Particle_pos(3)

            PartState(i-Nshift,4:6) =&
            (/IMD_array(7)*IMDLengthScale/IMDTimeScale,&
              IMD_array(8)*IMDLengthScale/IMDTimeScale,&
              IMD_array(9)*IMDLengthScale/IMDTimeScale/)

            xMin=MIN(Particle_pos(1),xMin)
            yMin=MIN(Particle_pos(2),yMin)
            zMin=MIN(Particle_pos(3),zMin)
            xMax=MAX(Particle_pos(1),xMax)
            yMax=MAX(Particle_pos(2),yMax)
            zMax=MAX(Particle_pos(3),zMax)
            ! check cutoff
            SELECT CASE(TRIM(IMDCutOff))
            CASE('no_cutoff') ! nothing to do
            CASE('Epot') ! kill particles that have Epot (i.e. they are in the solid body)
              IF(ABS(IMD_array(10)).GT.0.0)THEN ! IMD_array(10) is Epot
                Nshift=Nshift+1
              END IF
            CASE('coordinates') ! kill particles that are below a certain threshold in z-direction
              IF(IMD_array(4)*IMDLengthScale.GT.IMDCutOffxValue)THEN
                Nshift=Nshift+1
              END IF
            CASE('velocity') ! kill particles that are below a certain velocity threshold
              CALL abort(__STAMP__&
              ,'ERROR in particle_emission.f90: Error reading specified File (particle position) for SpaceIC=IMD!')
            END SELECT
          END IF
        END DO
        CLOSE(ioUnit)
        SWRITE(UNIT_stdOut,'(A,I15)')  "Particles Read: chunkSize = NbrOfParticle = ",(i-Nshift)-1
        chunkSize     = (i-Nshift)-1 ! don't change here, change at velocity
        NbrOfParticle = (i-Nshift)-1 ! don't change here, change at velocity
        SWRITE(UNIT_stdOut,'(A)') 'Min-Max particle positions from IMD source file:'
        SWRITE(UNIT_stdOut,'(A25,A25)')  "x-Min [nm]","x-Max [nm]"
        SWRITE(UNIT_stdOut,'(ES25.14E3,ES25.14E3)') xMin*1.e9,xMax*1.e9
        SWRITE(UNIT_stdOut,'(A25,A25)')  "y-Min [nm]","y-Max [nm]"
        SWRITE(UNIT_stdOut,'(ES25.14E3,ES25.14E3)') yMin*1.e9,yMax*1.e9
        SWRITE(UNIT_stdOut,'(A25,A25)')  "z-Min [nm]","z-Max [nm]"
        SWRITE(UNIT_stdOut,'(ES25.14E3,ES25.14E3)') zMin*1.e9,zMax*1.e9
        CALL PrintOption('IMD Particles Found','OUTPUT',IntOpt=(i-Nshift)-1)
      ELSE ! TRIM(IMDAtomFile) = 'no file found' -> exit
        Species(FractNbr)%Init(iInit)%velocityDistribution=''
      END IF
    END SELECT
    !------------------SpaceIC-cases: end-------------------------------------------------------------------------------------------
    chunkSize=chunkSize2

#if USE_MPI
 ELSE !no mpi root, nchunks=1
   chunkSize=0
 END IF
 IF(nChunks.GT.1) THEN
   ALLOCATE( PartMPIInsert%nPartsSend  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%nPartsRecv  (0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   ALLOCATE( PartMPIInsert%SendRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%RecvRequest (0:PartMPI%InitGroup(InitGroup)%nProcs-1,1:2), STAT=allocStat )
   ALLOCATE( PartMPIInsert%send_message(0:PartMPI%InitGroup(InitGroup)%nProcs-1), STAT=allocStat )
   PartMPIInsert%nPartsSend(:)=0
   DO i=1,chunkSize
     CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
     CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
     CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
     InsideMyBGM=.TRUE.
     IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
         (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
         (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
       InsideMyBGM=.FALSE.
     END If
     IF (InsideMyBGM) THEN
       IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
     END IF
     IF (InsideMyBGM) THEN
       DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
         iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
         tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
         IF(tProc.EQ.-1)CYCLE
         !IF(PartMPI%InitGroup(InitGroup)%COMM.EQ.MPI_COMM_NULL) THEN
         PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
       END DO
       PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)=&
              PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
     ELSE
       DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
!         IF (iProc.EQ.PartMPI%iProc) CYCLE
         PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
       END DO
     END IF
   END DO
 ELSE
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
      ALLOCATE( PartMPIInsert%send_message(0:0), STAT=allocStat )
      MessageSize=DimSend*chunkSize
      ALLOCATE( PartMPIInsert%send_message(0)%content(1:MessageSize), STAT=allocStat )
      PartMPIInsert%send_message(0)%content(:)=particle_positions(1:DimSend*chunkSize)
      DEALLOCATE(particle_positions, STAT=allocStat)
    END IF
 END IF
 IF (nChunks.GT.1) THEN
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      ! sent particles
      !--- MPI_ISEND lengths of lists of particles leaving local mesh
      CALL MPI_ISEND(PartMPIInsert%nPartsSend(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%SendRequest(iProc,1), IERROR)
      !--- MPI_IRECV lengths of lists of particles entering local mesh
      CALL MPI_IRECV(PartMPIInsert%nPartsRecv(iProc), 1, MPI_INTEGER, iProc, 1011+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, &
                     PartMPIInsert%RecvRequest(iProc,1), IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        ALLOCATE( PartMPIInsert%send_message(iProc)%content(1:DimSend*PartMPIInsert%nPartsSend(iProc)), STAT=allocStat )
      END IF
    END DO
    PartMPIInsert%nPartsSend(:)=0
    DO i=1,chunkSize
      CellX = INT((particle_positions(DimSend*(i-1)+1)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
      CellY = INT((particle_positions(DimSend*(i-1)+2)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
      CellZ = INT((particle_positions(DimSend*(i-1)+3)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
      InsideMyBGM=.TRUE.
      IF ((CellX.GT.GEO%FIBGMimax).OR.(CellX.LT.GEO%FIBGMimin) .OR. &
          (CellY.GT.GEO%FIBGMjmax).OR.(CellY.LT.GEO%FIBGMjmin) .OR. &
          (CellZ.GT.GEO%FIBGMkmax).OR.(CellZ.LT.GEO%FIBGMkmin)) THEN
        InsideMyBGM=.FALSE.
      END If
      IF (InsideMyBGM) THEN
        IF (.NOT.ALLOCATED(GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs)) InsideMyBGM=.FALSE.
      END IF
      IF (InsideMyBGM) THEN
        DO j=2,GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(1)+1
          iProc=GEO%FIBGM(CellX,CellY,CellZ)%ShapeProcs(j)
          tProc=PartMPI%InitGroup(InitGroup)%CommToGroup(iProc)
          IF(tProc.EQ.-1)CYCLE
          PartMPIInsert%nPartsSend(tProc)=PartMPIInsert%nPartsSend(tProc)+1
          k=PartMPIInsert%nPartsSend(tProc)
          PartMPIInsert%send_message(tProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
        PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)= &
            PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)+1
        k=PartMPIInsert%nPartsSend(PartMPI%InitGroup(InitGroup)%MyRank)
        PartMPIInsert%send_message(PartMPI%InitGroup(InitGroup)%MyRank)%content(DimSend*(k-1)+1:DimSend*k)=&
                                                                          particle_positions(DimSend*(i-1)+1:DimSend*i)
      ELSE
        DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
 !         IF (iProc.EQ.PartMPI%iProc) CYCLE
          PartMPIInsert%nPartsSend(iProc)=PartMPIInsert%nPartsSend(iProc)+1
          k=PartMPIInsert%nPartsSend(iProc)
          PartMPIInsert%send_message(iProc)%content(DimSend*(k-1)+1:DimSend*k)=particle_positions(DimSend*(i-1)+1:DimSend*i)
        END DO
      END IF
    END DO
    DEALLOCATE(particle_positions, STAT=allocStat)
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      !--- (non-blocking:) send messages to all procs receiving particles from myself
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_ISEND(PartMPIInsert%send_message(iProc)%content, DimSend*PartMPIInsert%nPartsSend(iProc),&
         MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr, PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%SendRequest(iProc,2), IERROR)
      END IF
    END DO
  END IF
ELSE ! mode.NE.1:
!--- RECEIVE:
  nChunksTemp=0
  IF(nChunks.EQ.1) THEN
    IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN !chunkSize can be 1 higher than NbrOfParticle for VPI+PartDens
       chunkSize=INT( REAL(SIZE(PartMPIInsert%send_message(0)%content)) / REAL(DimSend) )
       ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
       particle_positions(:)=PartMPIInsert%send_message(0)%content(:)
       DEALLOCATE( PartMPIInsert%send_message(0)%content )
       DEALLOCATE( PartMPIInsert%send_message )
    END IF
    !IF( (Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity.GT.0.)) .OR. &
    !    (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) ) THEN
      CALL MPI_BCAST(chunkSize, 1, MPI_INTEGER,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
    !ELSE
    !  chunkSize=NbrOfParticle
    !END IF
    IF(.NOT.PartMPI%InitGroup(InitGroup)%MPIROOT) THEN
      ALLOCATE(particle_positions(1:chunkSize*DimSend), STAT=allocStat)
    END IF
    CALL MPI_BCAST(particle_positions, chunkSize*DimSend, MPI_DOUBLE_PRECISION,0,PartMPI%InitGroup(InitGroup)%COMM,IERROR)
    nChunksTemp=1
  ELSE
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,1),msg_status(:),IERROR)
    END DO
    k=SUM(PartMPIInsert%nPartsRecv)
    ALLOCATE(particle_positions(1:k*DimSend), STAT=allocStat)
    k=0
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      IF (PartMPIInsert%nPartsRecv(iProc).GT.0) THEN
      !--- MPI_IRECV lengths of lists of particles entering local mesh
        CALL MPI_IRECV(particle_positions(k*DimSend+1), DimSend*PartMPIInsert%nPartsRecv(iProc),&
                                                  MPI_DOUBLE_PRECISION, iProc, 1022+FractNbr,   &
                                                  PartMPI%InitGroup(InitGroup)%COMM, PartMPIInsert%RecvRequest(iProc,2), IERROR)
        CALL MPI_WAIT(PartMPIInsert%RecvRequest(iProc,2),msg_status(:),IERROR)
        k=k+PartMPIInsert%nPartsRecv(iProc)
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsRecv )
    DEALLOCATE( PartMPIInsert%RecvRequest )
    DO iProc=0,PartMPI%InitGroup(InitGroup)%nProcs-1
      CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,1),msg_status(:),IERROR)
      IF (PartMPIInsert%nPartsSend(iProc).GT.0) THEN
        CALL MPI_WAIT(PartMPIInsert%SendRequest(iProc,2),msg_status(:),IERROR)
        DEALLOCATE( PartMPIInsert%send_message(iProc)%content )
      END IF
    END DO
    DEALLOCATE( PartMPIInsert%nPartsSend )
    DEALLOCATE( PartMPIInsert%send_message )
    DEALLOCATE( PartMPIInsert%SendRequest )
    chunkSize=k
    nChunks=1
  END IF
#endif
   ! each process checks which particle can be matched to its elements, counting the elements inside (local particles)
!   WRITE(*,*)'locating',chunkSize,'*',nChunks,' particles...'
!   WRITE(UNIT=debugFileName,FMT='(A,I2.2)')'prtcls_',PartMPI%iProc
!   OPEN(UNIT=130+PartMPI%iProc,FILE=debugFileName)
!   DO i=1,chunkSize*nChunks
!      WRITE(130+PartMPI%iProc,'(3(ES15.8))')particle_positions(i*3-2:i*3)
!   END DO
!   CLOSE(130+PartMPI%iProc)

#if USE_MPI
  ! in order to remove duplicated particles
  IF(nChunksTemp.EQ.1) THEN
    ALLOCATE(PartFoundInProc(1:2,1:ChunkSize),STAT=ALLOCSTAT)
      IF (ALLOCSTAT.NE.0) THEN
        CALL abort(&
__STAMP__,&
"abort: Error during emission in PartFoundInProc allocation")
      END IF
    PartFoundInProc=-1
  END IF
#endif /*USE_MPI*/

  mySumOfMatchedParticles=0
  ParticleIndexNbr = 1
  DO i=1,chunkSize*nChunks
    IF ((i.EQ.1).OR.PDM%ParticleInside(ParticleIndexNbr)) THEN
       ParticleIndexNbr = PDM%nextFreePosition(mySumOfMatchedParticles + 1 &
                                             + PDM%CurrentNextFreePosition)
    END IF
    IF (ParticleIndexNbr .ne. 0) THEN
       PartState(ParticleIndexNbr,1:DimSend) = particle_positions(DimSend*(i-1)+1:DimSend*(i-1)+DimSend)
       PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
       CALL LocateParticleInElement(ParticleIndexNbr,doHALO=.FALSE.)
       IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
          mySumOfMatchedParticles = mySumOfMatchedParticles + 1
          ! IF (VarTimeStep%UseVariableTimeStep) THEN
          !   VarTimeStep%ParticleTimeStep(ParticleIndexNbr) = &
          !     CalcVarTimeStep(PartState(ParticleIndexNbr,1), PartState(ParticleIndexNbr,2),PEM%Element(ParticleIndexNbr))
          ! END IF
          ! IF(RadialWeighting%DoRadialWeighting) THEN
          !    PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(ParticleIndexNbr,2),FractNbr,ParticleIndexNbr)
          ! END IF
#if USE_MPI
          IF(nChunksTemp.EQ.1) THEN
            ! mark elements with Rank and local found particle index
            PartFoundInProc(1,i)=MyRank
            PartFoundInProc(2,i)=mySumOfMatchedParticles
          END IF ! nChunks.EQ.1
#endif /*USE_MPI*/
       ELSE
          PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
       END IF
       IF (PDM%ParticleInside(ParticleIndexNbr)) THEN
         PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
         PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
       END IF
    ELSE
      CALL abort(&
__STAMP__&
,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
    END IF
  END DO

! we want always warnings to know if the emission has failed. if a timedisc does not require this, this
! timedisc has to be handled separately
#if USE_MPI
  mySumOfRemovedParticles=0
  IF(nChunksTemp.EQ.1) THEN
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,PartfoundInProc(1,:), ChunkSize, MPI_INTEGER, MPI_MAX &
                                                        , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
    ! loop over all particles and check, if particle is found in my proc
    ! proc with LARGES id gets the particle, all other procs remove the duplicated
    ! particle from their list
    DO i=1,chunkSize
      IF(PartFoundInProc(2,i).GT.-1)THEN ! particle has been previously found by MyRank
        IF(PartFoundInProc(1,i).NE.MyRank)THEN ! particle should not be found by MyRank
          !ParticleIndexNbr = PartFoundInProc(2,i)
          ParticleIndexNbr = PDM%nextFreePosition(PartFoundInProc(2,i) + PDM%CurrentNextFreePosition)
          IF(.NOT.PDM%ParticleInside(ParticleIndexNbr)) WRITE(UNIT_stdOut,*) ' Error in emission in parallel!!'
          PDM%ParticleInside(ParticleIndexNbr) = .FALSE.
          PDM%IsNewPart(ParticleIndexNbr)=.FALSE.
          ! correct number of found particles
          mySumOfRemovedParticles = mySumOfRemovedParticles +1
          ! set update next free position to zero for removed particle
          PDM%nextFreePosition(PartFoundInProc(2,i) + PDM%CurrentNextFreePosition) = 0
          !mySumOfMatchedParticles = mySumOfMatchedParticles -1
        END IF
      END IF
    END DO ! i=1,chunkSize
    DEALLOCATE(PartFoundInProc)
    mySumOfMatchedParticles = mySumOfMatchedParticles - mySumOfRemovedParticles
  END IF

  ! check the sum of the matched particles: did each particle find its "home"-CPU?
  CALL MPI_ALLREDUCE(mySumOfMatchedParticles, sumOfMatchedParticles, 1, MPI_INTEGER, MPI_SUM &
                                           , PartMPI%InitGroup(InitGroup)%COMM, IERROR)
#else
  ! im seriellen Fall kommen alle Partikel auf einen CPU,
  ! daher ist PIC%maxParticleNumber die harte Grenze
  sumOfMatchedParticles = mySumOfMatchedParticles
#endif

#if USE_MPI
  IF(PartMPI%InitGroup(InitGroup)%MPIRoot) THEN
#endif
    IF( Species(FractNbr)%Init(iInit)%VirtPreInsert .AND. (Species(FractNbr)%Init(iInit)%PartDensity .GT. 0.) ) THEN
      IF ((nbrOfParticle .NE. sumOfMatchedParticles).AND.OutputVpiWarnings) THEN
        SWRITE(UNIT_StdOut,'(A)')'WARNING in ParticleEmission_parallel:'
        SWRITE(UNIT_StdOut,'(A,I0)')'Fraction Nbr: ', FractNbr
        SWRITE(UNIT_StdOut,'(I0,A)') sumOfMatchedParticles, ' particles reached the domain when'
        SWRITE(UNIT_StdOut,'(I0,A)') NbrOfParticle, '(+1) velocities were calculated with vpi+PartDens'
      END IF
    ELSE
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
    END IF
#if USE_MPI
  END IF ! PartMPI%iProc.EQ.0
#endif

  ! Return the *local* NbrOfParticle so that the following Routines only fill in
  ! the values for the local particles
#if USE_MPI
  NbrOfParticle = mySumOfMatchedParticles + mySumOfRemovedParticles
#else
  NbrOfParticle = mySumOfMatchedParticles
#endif

  DEALLOCATE( particle_positions, STAT=allocStat )
  IF (allocStat .NE. 0) THEN
    CALL abort(&
__STAMP__&
,'ERROR in ParticleEmission_parallel: cannot deallocate particle_positions!')
  END IF
#if USE_MPI
END IF ! mode 1/2
#endif

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
  IF(Species(FractNbr)%Init(iInit)%VirtPreInsert) RETURN !velocities already set in SetParticlePosition!

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
      PartState(PositionNbr,4:6) = RandVal(1:3) * VeloIC
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
!     PartState(1,4:6) = tan_vec(1:3) * Species(1)%VeloIC
CASE('constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        IF (Is_ElemMacro) THEN
          IF (Species(FractNbr)%Init(iInit)%ElemVelocityICFileID.GT.0) THEN
            PartState(PositionNbr,4:6) = Species(FractNbr)%Init(iInit)%ElemVelocityIC(1:3,PEM%Element(PositionNbr))
          ELSE
            PartState(PositionNbr,4:6) = VeloVecIC(1:3) * VeloIC
          END IF
        ELSE
          PartState(PositionNbr,4:6) = VeloVecIC(1:3) * VeloIC
        END IF
     END IF
     i = i + 1
  END DO
CASE('radial_constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        Radius(1:3) = PartState(PositionNbr,1:3) - BasePointIC(1:3)
        !  Unity radius
        !Radius(1:3) = Radius(1:3) / RadiusIC
        Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
        PartState(PositionNbr,4:6) = Radius(1:3) * VeloIC
     END IF
     i = i + 1
  END DO
CASE('tangential_constant')
  i = 1
  DO WHILE (i .le. NbrOfParticle)
     PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
     IF (PositionNbr .ne. 0) THEN
        Radius(1:3) = PartState(PositionNbr,1:3) - BasePointIC(1:3)
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
          PartState(PositionNbr,4:6) = Vec3D+n_vec(1:3) * Vec1D
        ELSE ! no velocity spread
          ! If Gyrotron resonator: Add velocity in normal direction!
          IF (Alpha .gt. 0.) THEN
            n_vec = n_vec * ( 1 / Alpha )
          ELSE
            n_vec = 0
          END IF
          !  And finally the velocities
          IF(Rotation.EQ.1)THEN
            PartState(PositionNbr,4:6) = tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
          ELSE
            PartState(PositionNbr,4:6) = -tan_vec(1:3) * VeloIC + n_vec(1:3) * VeloIC
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
       x_2 = PartState(PositionNbr,1)
       y_2 = PartState(PositionNbr,2)
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
         Radius(1) = PartState(PositionNbr,1) - x_01
         Radius(2) = PartState(PositionNbr,2) - y_01
       ELSE
         Radius(1) = PartState(PositionNbr,1) - x_02
         Radius(2) = PartState(PositionNbr,2) - y_02
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
           PartState(PositionNbr,4:6) = (tan_vec(1:3) + n_vec(1:3)) * VeloIC
           IF (ABS(SQRT(PartState(PositionNbr,4)*PartState(PositionNbr,4) &
                      + PartState(PositionNbr,5)*PartState(PositionNbr,5))&
                      - VeloIC) .GT. 10.) THEN
             SWRITE(*,'(A,3(E21.14,X))') 'Velocity=', PartState(PositionNbr,4:6)
             CALL abort(&
__STAMP__&
,'ERROR in gyrotron_circle spaceIC!',PositionNbr)
           END If
           IF (PartState(PositionNbr,4).NE.PartState(PositionNbr,4) .OR. &
               PartState(PositionNbr,5).NE.PartState(PositionNbr,5) .OR. &
               PartState(PositionNbr,6).NE.PartState(PositionNbr,6)     ) THEN
             SWRITE(*,'(A,3(E21.14,X))') 'WARNING:! NaN: Velocity=', PartState(PositionNbr,4:6)
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
       PartState(PositionNbr,4:6) = Vec3D(1:3)
    END IF
  END DO
CASE('taylorgreenvortex')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
       CALL CalcVelocity_taylorgreenvortex(FractNbr, Vec3D, iInit=iInit, Element=PEM%Element(PositionNbr))
       PartState(PositionNbr,4:6) = Vec3D(1:3)
    END IF
  END DO
CASE('emmert')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      CALL CalcVelocity_emmert(FractNbr, iInit, Vec3D)
    END IF
    PartState(PositionNbr,4:6) = Vec3D(1:3)
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
        PartState(PositionNbr,4:6) = Vec3D(1:3)
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
       PartState(PositionNbr,4:6) = (PartState(PositionNbr,4:6) - v_sum(1:3)) * maxwellfac &
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
    PartState(PositionNbr,4:6) = velabs/Velosq*Vec3D
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
      PartState(PositionNbr,4:6) = Vec3D(1:3)
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
  ELSE                                                                            ! s^2(X)=1/(N-1)(N*E(X^2) - N*E(X)^2)
    v_sum(:) = 0.
    sigma(:) = 1.
  END IF

  ftl = 0
  DO i=1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(PositionNbr,4)   = (PartState(PositionNbr,4)  -v_sum(1)) * &
                                    SQRT(WeibelVeloPar**2/sigma(1)) *c
      PartState(PositionNbr,5:6) = (PartState(PositionNbr,5:6)-v_sum(2:3)) * &
                                    SQRT(WeibelVeloPer**2/sigma(2:3)) *c
      PartVelo = SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)

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

        PartState(PositionNbr,4:6) = Vec3D(1:3)

        PartState(PositionNbr,4)   = (Vec3D(1)  -v_sum(1)) * &
            SQRT(WeibelVeloPar**2/sigma(1)) *c
        PartState(PositionNbr,5:6) = (Vec3D(2:3)-v_sum(2:3)) * &
            SQRT(WeibelVeloPer**2/sigma(2:3)) *c
        PartVelo = SQRT(PartState(PositionNbr,4)**2+PartState(PositionNbr,5)**2+PartState(PositionNbr,6)**2)
      END DO
    END IF
  END DO

CASE('OneD-twostreaminstabilty')
  DO i = 1,NbrOfParticle
    PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
    IF (PositionNbr .NE. 0) THEN
      PartState(PositionNbr,4) = OneDTwoStreamVelo
      PartState(PositionNbr,5:6) = OneDTwoStreamTransRatio
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
      PartState(PositionNbr,5) = RandVal(1)*OneDTwoStreamTransRatio* &
                                                   OneDTwoStreamVelo
      PartState(PositionNbr,6) = RandVal(2)*OneDTwoStreamTransRatio* &
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
