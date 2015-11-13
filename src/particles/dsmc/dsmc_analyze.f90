#include "boltzplatz.h"

MODULE MOD_DSMC_Analyze
!===================================================================================================================================
! Module for DSMC Sampling and Output
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_data_sampling
  MODULE PROCEDURE DSMC_data_sampling
END INTERFACE

INTERFACE WriteDSMCToHDF5
  MODULE PROCEDURE WriteDSMCToHDF5
END INTERFACE

INTERFACE DSMC_output_calc
  MODULE PROCEDURE DSMC_output_calc
END INTERFACE

INTERFACE CalcTVib
  MODULE PROCEDURE CalcTVib
END INTERFACE

INTERFACE CalcSurfaceValues
  MODULE PROCEDURE CalcSurfaceValues
END INTERFACE

INTERFACE CalcTelec
  MODULE PROCEDURE CalcTelec
END INTERFACE

INTERFACE CalcTVibPoly
  MODULE PROCEDURE CalcTVibPoly
END INTERFACE

INTERFACE CalcMeanFreePath
  MODULE PROCEDURE CalcMeanFreePath
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_data_sampling,  DSMCHO_data_sampling, CalcMeanFreePath,WriteDSMCToHDF5
PUBLIC :: DSMC_output_calc, CalcTVib, CalcSurfaceValues, CalcTelec, CalcTVibPoly, InitHODSMC, WriteDSMCHOToHDF5

!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_data_sampling()
!===================================================================================================================================
! Sampling of variables velocity and energy for DSMC
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, PartStateIntEn, CollisMode, SpecDSMC, useDSMC
  USE MOD_Particle_Vars,          ONLY : PEM, PDM, PartSpecies, PartState, PartMPF, usevMPF
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iElem
!===================================================================================================================================
  DSMC%SampNum = DSMC%SampNum + 1
  DO iPart = 1, PDM%ParticleVecLength
    IF (PDM%ParticleInside(ipart)) THEN
      iElem = PEM%Element(iPart)
      IF (usevMPF) THEN
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartState(iPart,4:6) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) &
                                                       + PartState(iPart,4:6)**2 * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%SimPartNum  = SampDSMC(iElem,PartSpecies(iPart))%SimPartNum + 1
        ! real number of particles in SampDSMC(iElem,PartSpecies(iPart))%PartNum for vMPF
        IF (useDSMC) THEN
          IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.(SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
            SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                         + PartStateIntEn(iPart,1) * PartMPF(iPart)
            SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                         + PartStateIntEn(iPart,2) * PartMPF(iPart)
          END IF
          IF (DSMC%ElectronicState) THEN
            SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
              + PartStateIntEn(iPart,3) * PartMPF(iPart)
          END IF
        END IF
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartState(iPart,4:6)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) &
                                                       + PartState(iPart,4:6)**2
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
        IF (useDSMC) THEN
          IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
            IF (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2) THEN
              SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                           + PartStateIntEn(iPart,1)
              SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                           + PartStateIntEn(iPart,2)
            END IF
          END IF
          IF (DSMC%ElectronicState) THEN
            SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
                                                         + PartStateIntEn(iPart,3)
          END IF
        END IF
      END IF
    END IF
  END DO
END SUBROUTINE DSMC_data_sampling


SUBROUTINE DSMC_output_calc
!===================================================================================================================================
! Compute the macroscopic values from samples
!===================================================================================================================================
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, MacroDSMC, CollisMode, SpecDSMC, realtime, useDSMC
  USE MOD_Mesh_Vars,              ONLY : nElems,MeshFile
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst, Species, usevMPF
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_Particle_Vars,          ONLY : Time
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, iter, dt
  USE MOD_Restart_Vars,           ONLY : RestartTime
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES     
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem, iSpec
REAL                           :: TVib_TempFac, MolecPartNum, HeavyPartNum     
!===================================================================================================================================

  ALLOCATE(MacroDSMC(nElems,nSpecies + 1))
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV(1)  = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV(2)  = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV(3)  = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV(4)  = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(1) = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(2) = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartV2(3) = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%Temp(1)   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%Temp(2)   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%Temp(3)   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%Temp(4)   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%PartNum   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%NumDens   = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%TVib      = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%TRot      = 0
  MacroDSMC(1:nElems,1:nSpecies+1)%TElec     = 0
  IF (useDSMC) DSMC%CollProbOut(1:nElems,2) = 0.0       ! resetting the time-averaged mean collision probability

  DO iSpec = 1, nSpecies
    DO iElem = 1, nElems ! element/cell main loop    
      IF(SampDSMC(iElem,iSpec)%PartNum.gt. 0) THEN
! compute flow velocity
        MacroDSMC(iElem,iSpec)%PartV(1) = SampDSMC(iElem,iSpec)%PartV(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(2) = SampDSMC(iElem,iSpec)%PartV(2) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(3) = SampDSMC(iElem,iSpec)%PartV(3) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(4) = SQRT( MacroDSMC(iElem,iSpec)%PartV(1)**2 &
                                              + MacroDSMC(iElem,iSpec)%PartV(2)**2 &
                                              + MacroDSMC(iElem,iSpec)%PartV(3)**2 )
! compute flow Temperature
        MacroDSMC(iElem,iSpec)%PartV2(1) = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV2(2) = SampDSMC(iElem,iSpec)%PartV2(2) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV2(3) = SampDSMC(iElem,iSpec)%PartV2(3) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%Temp(1:3) = Species(iSpec)%MassIC / BoltzmannConst &
                                             * (MacroDSMC(iElem,iSpec)%PartV2(1:3) - MacroDSMC(iElem,iSpec)%PartV(1:3)**2)
        MacroDSMC(iElem,iSpec)%Temp(4) = (MacroDSMC(iElem,iSpec)%Temp(1) &
                                        + MacroDSMC(iElem,iSpec)%Temp(2) &
                                        + MacroDSMC(iElem,iSpec)%Temp(3)) / 3
! compute density
        MacroDSMC(iElem,iSpec)%PartNum = SampDSMC(iElem,iSpec)%PartNum / REAL(DSMC%SampNum)
        ! if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
        IF (usevMPF) THEN
          MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
        ELSE 
          MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
        END IF
! compute internal energies / has to be changed for vfd 
        IF (useDSMC) THEN
          IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
            IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
              IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
                TVib_TempFac=SampDSMC(iElem,iSpec)%EVib / (SampDSMC(iElem,iSpec)%PartNum*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
                IF (TVib_TempFac.LE.DSMC%GammaQuant) THEN
                  MacroDSMC(iElem,iSpec)%TVib = 0.0
                ELSE
                  MacroDSMC(iElem,iSpec)%TVib = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(TVib_TempFac-DSMC%GammaQuant))
                END IF
              ELSE                                            ! TSHO-model
                MacroDSMC(iElem,iSpec)%TVib = CalcTVib(SpecDSMC(iSpec)%CharaTVib &
                  , SampDSMC(iElem,iSpec)%EVib/SampDSMC(iElem,iSpec)%PartNum, SpecDSMC(iSpec)%MaxVibQuant)
              END IF
              MacroDSMC(iElem,iSpec)%TRot = SampDSMC(iElem, iSpec)%ERot/(BoltzmannConst*SampDSMC(iElem,iSpec)%PartNum)
              IF (DSMC%ElectronicState) THEN
                MacroDSMC(iElem,iSpec)%TElec= CalcTelec( SampDSMC(iElem,iSpec)%EElec/SampDSMC(iElem,iSpec)%PartNum, iSpec)
              END IF
            END IF
          END IF
        END IF
      END IF
    END DO
  END DO
! compute total values
  DO iElem = 1, nElems ! element/cell main loop 
    MolecPartNum = 0
    HeavyPartNum = 0
    DO iSpec = 1, nSpecies
      IF(SampDSMC(iElem,iSpec)%PartNum.GT. 0) THEN
        MacroDSMC(iElem,nSpecies + 1)%PartV(1) = MacroDSMC(iElem,nSpecies + 1)%PartV(1) &
                                               + MacroDSMC(iElem,iSpec)%PartV(1) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%PartV(2) = MacroDSMC(iElem,nSpecies + 1)%PartV(2) &
                                               + MacroDSMC(iElem,iSpec)%PartV(2) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%PartV(3) = MacroDSMC(iElem,nSpecies + 1)%PartV(3) &
                                               + MacroDSMC(iElem,iSpec)%PartV(3) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = MacroDSMC(iElem,nSpecies + 1)%Temp(1) &
                                               + MacroDSMC(iElem,iSpec)%Temp(1) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = MacroDSMC(iElem,nSpecies + 1)%Temp(2) &
                                               + MacroDSMC(iElem,iSpec)%Temp(2) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = MacroDSMC(iElem,nSpecies + 1)%Temp(3) &
                                               + MacroDSMC(iElem,iSpec)%Temp(3) * MacroDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,nSpecies + 1)%PartNum  = MacroDSMC(iElem,nSpecies + 1)%PartNum + MacroDSMC(iElem,iSpec)%PartNum
        IF (useDSMC) THEN
          IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
            IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
              MacroDSMC(iElem,nSpecies + 1)%TVib = MacroDSMC(iElem,nSpecies + 1)%TVib &
                                                 + MacroDSMC(iElem,iSpec)%TVib * MacroDSMC(iElem,iSpec)%PartNum
              MacroDSMC(iElem,nSpecies + 1)%TRot = MacroDSMC(iElem,nSpecies + 1)%TRot &
                                                 + MacroDSMC(iElem,iSpec)%TRot * MacroDSMC(iElem,iSpec)%PartNum
              MolecPartNum                       = MolecPartNum + MacroDSMC(iElem,iSpec)%PartNum
            END IF
            IF ( DSMC%ElectronicState .AND. (SpecDSMC(iSpec)%InterID.NE.4) ) THEN
              MacroDSMC(iElem,nSpecies + 1)%TElec = MacroDSMC(iElem, nSpecies+1)%TElec &
                                                  + MacroDSMC(iElem,iSpec)%TElec * MacroDSMC(iElem,iSpec)%PartNum
              HeavyPartNum                        = HeavyPartNum + MacroDSMC(iElem,iSpec)%PartNum
            END IF
          END IF
        END IF
      END IF
    END DO
    IF(MacroDSMC(iElem,nSpecies + 1)%PartNum.GT. 0) THEN
      MacroDSMC(iElem,nSpecies + 1)%PartV(1) = MacroDSMC(iElem,nSpecies + 1)%PartV(1) / MacroDSMC(iElem,nSpecies + 1)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%PartV(2) = MacroDSMC(iElem,nSpecies + 1)%PartV(2) / MacroDSMC(iElem,nSpecies + 1)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%PartV(3) = MacroDSMC(iElem,nSpecies + 1)%PartV(3) / MacroDSMC(iElem,nSpecies + 1)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = MacroDSMC(iElem,nSpecies + 1)%Temp(1)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = MacroDSMC(iElem,nSpecies + 1)%Temp(2)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
      MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = MacroDSMC(iElem,nSpecies + 1)%Temp(3)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (MolecPartNum.GT. 0)) THEN
                MacroDSMC(iElem,nSpecies + 1)%TVib = MacroDSMC(iElem,nSpecies + 1)%TVib / MolecPartNum
                MacroDSMC(iElem,nSpecies + 1)%TRot = MacroDSMC(iElem,nSpecies + 1)%TRot / MolecPartNum
      END IF
      IF ( DSMC%ElectronicState) THEN
        MacroDSMC(iElem,nSpecies + 1)%TElec = MacroDSMC(iElem, nSpecies+1)%TElec / HeavyPartNum
      END IF
    END IF
    MacroDSMC(iElem,nSpecies + 1)%PartV(4) = SQRT(MacroDSMC(iElem,nSpecies + 1)%PartV(1)**2 &
                                            + MacroDSMC(iElem,nSpecies + 1)%PartV(2)**2 &
                                            + MacroDSMC(iElem,nSpecies + 1)%PartV(3)**2)
    MacroDSMC(iElem,nSpecies + 1)%Temp(4)  = (MacroDSMC(iElem,nSpecies + 1)%Temp(1) &
                                            + MacroDSMC(iElem,nSpecies + 1)%Temp(2) &
                                            + MacroDSMC(iElem,nSpecies + 1)%Temp(3)) /3
    MacroDSMC(iElem,nSpecies + 1)%NumDens = MacroDSMC(iElem,nSpecies + 1)%PartNum * Species(1)%MacroParticleFactor & 
                                   / GEO%Volume(iElem) ! the calculation is limitied for MPF = const.
    IF (usevMPF) THEN
      MacroDSMC(iElem,nSpecies + 1)%PartNum = 0
      DO iSpec = 1, nSpecies
        MacroDSMC(iElem,iSpec)%PartNum = SampDSMC(iElem,iSpec)%SimPartNum / REAL(DSMC%SampNum)
        MacroDSMC(iElem,nSpecies + 1)%PartNum = MacroDSMC(iElem,nSpecies + 1)%PartNum + MacroDSMC(iElem,iSpec)%PartNum
      END DO
    END IF
    IF (useDSMC) THEN
      IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
        DSMC%CollProbOut(iElem,2) = DSMC%CollProbSamp(iElem) / iter
      ELSE
        DSMC%CollProbOut(iElem,2) = DSMC%CollProbSamp(iElem)*dt / (Time-(1-DSMC%TimeFracSamp)*TEnd)
      END IF
    END IF
  END DO
  CALL WriteDSMCToHDF5(TRIM(MeshFile),realtime)
  DEALLOCATE(MacroDSMC)
END SUBROUTINE DSMC_output_calc


SUBROUTINE WriteDSMCSurfToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
! Writes DSMC surface values to HDF5
!===================================================================================================================================
! MODULES
   USE MOD_Globals
   USE MOD_PreProc
   USE MOD_io_HDF5
   USE MOD_HDF5_output,   ONLY:WriteArrayToHDF5,WriteAttributeToHDF5,WriteHDF5Header
   USE MOD_PARTICLE_Vars, ONLY:nSpecies
   USE MOD_Mesh_Vars,     ONLY:nGlobalElems,offsetSurfElem
   USE MOD_DSMC_Vars,     ONLY:SurfMesh, MacroSurfaceVal , CollisMode
   USE MOD_Output_Vars,   ONLY:ProjectName
#ifdef MPI
   USE MOD_MPI_Vars,      ONLY:offsetSurfElemMPI
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*),INTENT(IN)          :: MeshFileName
  REAL,INTENT(IN)                       :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  CHARACTER(LEN=255)                 :: FileName,FileString,Statedummy
  INTEGER                             :: nVal
!===================================================================================================================================
  SWRITE(*,*) ' WRITE DSMCSurfSTATE TO HDF5 FILE...'
  FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
  FileString=TRIM(FileName)//'.h5'
#ifdef MPI
  CALL OpenDataFile(FileString,create=.TRUE.,single=.FALSE.)
#else
  CALL OpenDataFile(FileString,create=.TRUE.)
#endif

  Statedummy = 'DSMCSurfState'
  CALL WriteHDF5Header(Statedummy,File_ID)

  nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
#ifdef MPI
  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceX', rank=1,&
                        nValGlobal=(/offsetSurfElemMPI(nProcessors)/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem  /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(1))

  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceY', rank=1,&
                        nValGlobal=(/offsetSurfElemMPI(nProcessors)/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem  /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(2))

  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceZ', rank=1,&
                        nValGlobal=(/offsetSurfElemMPI(nProcessors)/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem  /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(3))

  CALL WriteArrayToHDF5(DataSetName='DSMC_Heatflux', rank=1,&
                        nValGlobal=(/offsetSurfElemMPI(nProcessors)/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem  /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Heatflux)

  CALL WriteArrayToHDF5(DataSetName='DSMC_Counter', rank=1,&
                        nValGlobal=(/offsetSurfElemMPI(nProcessors)/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem  /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%CounterOut)
#else
  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceX', rank=1,&
                        nValGlobal=(/SurfMesh%nSurfaceBCSides/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(1))

  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceY', rank=1,&
                        nValGlobal=(/SurfMesh%nSurfaceBCSides/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(2))

  CALL WriteArrayToHDF5(DataSetName='DSMC_ForceZ', rank=1,&
                        nValGlobal=(/SurfMesh%nSurfaceBCSides/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Force(3))

  CALL WriteArrayToHDF5(DataSetName='DSMC_Heatflux', rank=1,&
                        nValGlobal=(/SurfMesh%nSurfaceBCSides/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%Heatflux)

  CALL WriteArrayToHDF5(DataSetName='DSMC_Counter', rank=1,&
                        nValGlobal=(/SurfMesh%nSurfaceBCSides/),&
                        nVal=      (/SurfMesh%nSurfaceBCSides/),&
                        offset=    (/ offsetSurfElem /),&
                        collective=.FALSE., existing=.FALSE., RealArray=MacroSurfaceVal(:)%CounterOut)
#endif
  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)

  CALL CloseDataFile()

END SUBROUTINE WriteDSMCSurfToHDF5


SUBROUTINE CalcSurfaceValues
!===================================================================================================================================
! Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_DSMC_Vars,       ONLY:SurfMesh,  MacroSurfaceVal, DSMC, realtime!,SampWall
  USE MOD_Particle_Vars,   ONLY:Time, WriteMacroValues, nSpecies, MacroValSampTime
  USE MOD_TimeDisc_Vars,   ONLY:TEnd
  USE MOD_Mesh_Vars,       ONLY:MeshFile
  USE MOD_Restart_Vars,    ONLY:RestartTime  
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                            :: iSpec
  REAL                               :: TimeSample
  INTEGER, ALLOCATABLE               :: CounterTotal(:), SumCounterTotal(:)              ! Total Wall-Collision counter
!===================================================================================================================================
 
STOP
 
#ifdef MPI
  CALL MPISurfaceValuesSend()  
#endif

  ALLOCATE(MacroSurfaceVal(SurfMesh%nSurfaceBCSides))
  IF (DSMC%CalcSurfCollis_Output) THEN
    ALLOCATE(CounterTotal(1:nSpecies+1))
    ALLOCATE(SumCounterTotal(1:nSpecies+1))
    CounterTotal(1:nSpecies)=0
    SumCounterTotal(1:nSpecies+1)=0
  END IF

  IF (WriteMacroValues) THEN
    TimeSample = Time - MacroValSampTime !elapsed time since last sampling (variable dt's possible!)
    MacroValSampTime = Time
  ELSE IF (RestartTime.GT.(1-DSMC%TimeFracSamp)*TEnd) THEN
    TimeSample = Time - RestartTime
  ELSE
    TimeSample = (Time-(1-DSMC%TimeFracSamp)*TEnd)
  END IF
!
!  DO iElem=1,SurfMesh%nSurfaceBCSides
!    ALLOCATE(MacroSurfaceVal(iElem)%Counter(1:nSpecies))
!    MacroSurfaceVal(iElem)%Counter(1:nSpecies)=0.0
!    MacroSurfaceVal(iElem)%CounterOut=0.0
!    MacroSurfaceVal(iElem)%Heatflux = (SampWall(iElem)%Energy(1)+SampWall(iElem)%Energy(4)+SampWall(iElem)%Energy(7) &
!                                      -SampWall(iElem)%Energy(3)-SampWall(iElem)%Energy(6)-SampWall(iElem)%Energy(9))&
!                                     /(SurfMesh%SurfaceArea(iElem) * TimeSample)
!    MacroSurfaceVal(iElem)%Force(1) = SampWall(iElem)%Force(1) /(SurfMesh%SurfaceArea(iElem) * TimeSample)
!    MacroSurfaceVal(iElem)%Force(2) = SampWall(iElem)%Force(2) /(SurfMesh%SurfaceArea(iElem) * TimeSample)
!    MacroSurfaceVal(iElem)%Force(3) = SampWall(iElem)%Force(3) / (SurfMesh%SurfaceArea(iElem) * TimeSample)
!    DO iSpec=1,nSpecies
!      MacroSurfaceVal(iElem)%Counter(iSpec) = SampWall(iElem)%Counter(iSpec) / TimeSample
!      IF (DSMC%CalcSurfCollis_Output) CounterTotal(iSpec) = CounterTotal(iSpec) + INT(SampWall(iElem)%Counter(iSpec))
!      IF (DSMC%CalcSurfCollis_SpeciesFlags(iSpec)) THEN !Sum up all Collisions with SpeciesFlags for output
!        MacroSurfaceVal(iElem)%CounterOut = MacroSurfaceVal(iElem)%CounterOut &
!                                          + MacroSurfaceVal(iElem)%Counter(iSpec)
!      END IF
!    END DO
!  END DO

  IF (DSMC%CalcSurfCollis_Output) THEN
#ifdef MPI
    CALL MPI_REDUCE(CounterTotal,SumCounterTotal(1:nSpecies),nSpecies,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
#endif
    DO iSpec=1,nSpecies
      IF (DSMC%CalcSurfCollis_SpeciesFlags(iSpec)) THEN !Sum up all Collisions with SpeciesFlags for output
        SumCounterTotal(nSpecies+1) = SumCounterTotal(nSpecies+1) + SumCounterTotal(iSpec)
      END IF
    END DO
    SWRITE(*,*) ' The following species swaps at walls have been sampled:'
    DO iSpec=1,nSpecies
      SWRITE(*,'(A9,I2,A2,E16.9,A6)') ' Species ',iSpec,': ',REAL(SumCounterTotal(iSpec)) / TimeSample,' MP/s;'
    END DO
    SWRITE(*,'(A23,E16.9,A6)') ' All with SpeciesFlag: ',REAL(SumCounterTotal(nSpecies+1)) / TimeSample,' MP/s.'
    DEALLOCATE(CounterTotal)
    DEALLOCATE(SumCounterTotal)
  END IF

  CALL WriteDSMCSurfToHDF5(TRIM(MeshFile),realtime)
  DEALLOCATE(MacroSurfaceVal)

END SUBROUTINE CalcSurfaceValues


#ifdef MPI
SUBROUTINE MPISurfaceValuesSend()
!===================================================================================================================================
! Sends surface values of the halo cells to the respective processor
!===================================================================================================================================
! MODULES
!  USE MOD_part_MPI_Vars, ONLY : PartMPI, MPIGEO
  USE mpi
!  USE MOD_Globals,           ONLY:IERROR, MPISTATUS
!  USE MOD_DSMC_Vars,         ONLY:SurfMesh, SampWall, SampWallHaloCell
!  USE MOD_Mesh_Vars,         ONLY:ElemToSide
!  USE MOD_Particle_Vars,     ONLY:nSpecies
!  USE MOD_Particle_MPI_Vars, ONLY:PartMPI
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!  TYPE tMPISurfContent
!    REAL, POINTER                    :: content(:) =>NULL()
!  END TYPE
!  INTEGER                            :: iSide, iProc, Element, iLocSide, iCount, iSampWallAlloc
!  INTEGER                            :: SurfSideID    
!  INTEGER, ALLOCATABLE               :: RecvMsgSurfs(:), SendMsgSurfs(:), PosCount(:)
!  TYPE(tMPISurfContent), POINTER     :: SendContent(:) => NULL()
!  TYPE(tMPISurfContent), POINTER     :: RecvContent(:) => NULL()
!===================================================================================================================================

!  ALLOCATE(SendMsgSurfs(0:PartMPI%nProcs-1))
!  ALLOCATE(RecvMsgSurfs(0:PartMPI%nProcs-1))
!  ALLOCATE(SendContent (0:PartMPI%nProcs-1))
!  ALLOCATE(RecvContent (0:PartMPI%nProcs-1))
!  ALLOCATE(PosCount    (0:PartMPI%nProcs-1))
!  PosCount    (0:PartMPI%nProcs-1) = 0
!  SendMsgSurfs(0:PartMPI%nProcs-1) = 0
!  RecvMsgSurfs(0:PartMPI%nProcs-1) = 0

!---- calculation of send-massages number for each proc
STOP
!  DO iSide = 1, SIZE(MPIGEO%BC,2)  !---- loop over all halosides that are defined as BC
!    IF(SurfMesh%HaloSideIDToSurfSideMap(iSide).NE.0) THEN  !---- only surfaces (=wall-sides)   
!      ! get halo cells ElemID. Not quite sure if ELEM_ID or NB_ELEM_ID. 
!      ! maybe check this later and delete IF-query here to save cpu time
!      IF (MPIGEO%SideToElem(1,iSide).NE.-1) THEN
!       Element = MPIGEO%SideToElem(1,iSide)
!      ELSE
!       Element = MPIGEO%SideToElem(2,iSide)
!      END IF
!      iProc = MPIGEO%ElemMPIID(Element)
!      SendMsgSurfs(iProc) = SendMsgSurfs(iProc) + 1
!    END IF
!  END DO
!---- comunicate send-massages number (=number of surfaces) for each proc
!  DO iProc=0, PartMPI%nProcs-1 
!    IF (PartMPI%iProc.LT.iProc) THEN
!      CALL MPI_SEND(SendMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)    
!      CALL MPI_RECV(RecvMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)      
!    ELSE IF (PartMPI%iProc.GT.iProc) THEN
!      CALL MPI_RECV(RecvMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)
!      CALL MPI_SEND(SendMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PartMPI%COMM,IERROR)
!    END IF
!  END DO
!!---- allocate (send and recv) massages size (15 * send-massages number)
!!---- 1=iLocSide; 2=target element; 3-11=Energy(1-9); 12-14=Force(1-3); 14+nSpecies=Counter(1:nSpecies)  per surface
!  DO iProc=0, PartMPI%nProcs-1
!    IF (SendMsgSurfs(iProc).NE.0) THEN
!      ALLOCATE(SendContent(iProc)%content((14+nSpecies)*SendMsgSurfs(iProc)))
!    END IF
!    IF (RecvMsgSurfs(iProc).NE.0) THEN
!      ALLOCATE(RecvContent(iProc)%content((14+nSpecies)*RecvMsgSurfs(iProc)))
!    END IF
!  END DO
!---- build massage
!  DO iSide = 1, SIZE(MPIGEO%BC,2)
!    IF(SurfMesh%HaloSideIDToSurfSideMap(iSide).NE.0) THEN     
!      ! get halo cells ElemID. Not quite sure if ELEM_ID or NB_ELEM_ID. 
!      ! maybe check this later and delete IF-query here to save cpu time
!      IF (MPIGEO%SideToElem(1,iSide).NE.-1) THEN
!        Element = MPIGEO%SideToElem(1,iSide)
!      ELSE
!        Element = MPIGEO%SideToElem(2,iSide)
!      END IF
!      iProc = MPIGEO%ElemMPIID(Element)
!      PosCount(iProc) = PosCount(iProc) + 1
!      DO iLocSide=1,6  !---- search for iLocSide
!        IF (MPIGEO%ElemToSide(1,iLocSide,Element).EQ.iSide) THEN
!          SendContent(iProc)%content(PosCount(iProc))= REAL(iLocSide)  !---- 1=iLocSide
!          EXIT
!        END IF
!      END DO
!      PosCount(iProc) = PosCount(iProc) + 1
!      SendContent(iProc)%content(PosCount(iProc))= REAL(MPIGEO%NativeElemID(Element))  !---- 2=target element
!      PosCount(iProc) = PosCount(iProc) + 1
!      SendContent(iProc)%content(PosCount(iProc):PosCount(iProc)+8)= &  !---- 3-11=Energy(1-9)
!                      SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Energy(1:9)
!      SendContent(iProc)%content(PosCount(iProc)+9:PosCount(iProc)+11)= &  !---- 12-14=Force(1-3)
!                      SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Force(1:3)
!      SendContent(iProc)%content(PosCount(iProc)+12:PosCount(iProc)+11+nSpecies)= &  !---- 14+nSpecies=Counter(1:nSpecies)
!                      SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Counter(1:nSpecies)
!      PosCount(iProc) = PosCount(iProc)+11+nSpecies
!    END IF
!  END DO
!!---- comunication
!  DO iProc=0, PartMPI%nProcs-1          
!    IF (PartMPI%iProc.LT.iProc) THEN
!      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
!          (14+nSpecies)*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PartMPI%COMM,IERROR)    
!      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, &
!          (14+nSpecies)*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR) 
!    ELSE IF (PartMPI%iProc.GT.iProc) THEN
!      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, & 
!          (14+nSpecies)*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PartMPI%COMM,MPISTATUS,IERROR)  
!      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
!          (14+nSpecies)*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PartMPI%COMM,IERROR)    
!    END IF
!  END DO
!
!!---- sum up surface values  
!  DO iProc=0, PartMPI%nProcs-1
!!    IF (RecvMsgSurfs(iProc).ne.0) print*,RecvMsgSurfs(iProc)
!    IF ((iProc.NE.PartMPI%iProc).AND.(RecvMsgSurfs(iProc).NE.0)) THEN
!      DO iCount = 1, RecvMsgSurfs(iProc)
!        SurfSideID = SurfMesh%GlobSideToSurfSideMap( &  !---- get local surfaceID (iLocSide(target) = iLocSide(local))
!        ElemToSide(1,INT(RecvContent(iProc)%content(iCount*(14+nSpecies)-13-nSpecies)), &
!        INT(RecvContent(iProc)%content(iCount*(14+nSpecies)-12-nSpecies)))) 
!        SampWall(SurfSideID)%Energy(1:9) = SampWall(SurfSideID)%Energy(1:9) &
!        + RecvContent(iProc)%content(iCount*(14+nSpecies)-11-nSpecies:iCount*(14+nSpecies)-3-nSpecies)
!        SampWall(SurfSideID)%Force(1:3) = SampWall(SurfSideID)%Force(1:3) &
!        + RecvContent(iProc)%content(iCount*(14+nSpecies)-2-nSpecies:iCount*(14+nSpecies)-nSpecies)
!        SampWall(SurfSideID)%Counter(1:nSpecies) = SampWall(SurfSideID)%Counter(1:nSpecies) &
!        + RecvContent(iProc)%content(iCount*(14+nSpecies)+1-nSpecies:iCount*(14+nSpecies))
!      END DO
!    END IF
!  END DO  
!
!  DO iProc=0, PartMPI%nProcs-1
!    IF (ASSOCIATED(SendContent(iProc)%content)) DEALLOCATE(SendContent(iProc)%content)
!    IF (ASSOCIATED(RecvContent(iProc)%content)) DEALLOCATE(RecvContent(iProc)%content)
!  END DO
!  DEALLOCATE(SendMsgSurfs)
!  DEALLOCATE(RecvMsgSurfs)
!  DEALLOCATE(SendContent)
!  DEALLOCATE(RecvContent)
!  DEALLOCATE(PosCount)
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(1) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(2) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(3) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(4) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(5) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(6) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(7) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(8) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Energy(9) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(1) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(2) = 0.0
!  SampWallHaloCell(1:SurfMesh%nHaloSurfaceBCSides)%Force(3) = 0.0
!  DO iSampWallAlloc=1,SurfMesh%nHaloSurfaceBCSides
!    SampWallHaloCell(iSampWallAlloc)%Counter(1:nSpecies) = 0.0
!  END DO


END SUBROUTINE MPISurfaceValuesSend

#endif


REAL FUNCTION CalcTVib(ChaTVib,MeanEVib,nMax)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
  USE MOD_Globals,                ONLY : abort
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : DSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL, INTENT(IN)                :: ChaTVib,MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
  INTEGER, INTENT(IN)             :: nMax              ! INT(CharaTDisss/CharaTVib) + 1 
  REAL(KIND=8)                    :: LowerVal, UpperVal, MiddleVal, MaxPosiVal  ! upper and lower value of zero point search 
  REAl(KIND=8)                    :: eps_prec=1.0e-5   ! precision of zero point search
  REAL(KIND=8)                    :: ZeroVal1, ZeroVal2 ! both fuction values to compare
!===================================================================================================================================

  IF (MeanEVib.GT.0) THEN
    !.... Initial limits for a: lower limit = very small value
    !                           upper limit = max. value allowed by system
    !     zero point = CharaTVib / TVib
    LowerVal  = 1.0/(2.0*nMax)                                    ! Tvib is max for nMax => lower limit = 1.0/nMax
    UpperVal  = LOG(HUGE(MiddleVal*nMax))/nMax-1.0/(2.0 * nMax)   ! upper limit = for max possible EXP(nMax*MiddleVal)-value
    MaxPosiVal = LOG(HUGE(MaxPosiVal))  ! maximum value possible in system
    DO WHILE (ABS(LowerVal-UpperVal).GT.eps_prec)                      !  Let's search the zero point by bisection
      MiddleVal = 0.5*(LowerVal+UpperVal)

      IF ((LowerVal.GT.MaxPosiVal).OR.(MiddleVal.GT.MaxPosiVal)) THEN
         CALL Abort(&
           __STAMP__,&
          'Cannot find zero point in TVib Calculation Function! CharTVib:',RealInfoOpt=ChaTVib)
      END IF
      
      ! Calc of actual function values
      ZeroVal1 = DSMC%GammaQuant + 1/(EXP(LowerVal)-1) - nMax/(EXP(nMax*LowerVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
      ZeroVal2 = DSMC%GammaQuant + 1/(EXP(MiddleVal)-1) - nMax/(EXP(nMax*MiddleVal)-1) - MeanEVib/(ChaTVib*BoltzmannConst)
      ! decision of direction of bisection
      IF (ZeroVal1*ZeroVal2.LT.0) THEN
        UpperVal = MiddleVal
      ELSE
        LowerVal = MiddleVal
      END IF
    END DO
    CalcTVib = ChaTVib/LowerVal ! LowerVal = CharaTVib / TVib
  ELSE
    CalcTVib = 0
  END IF  

  RETURN

END FUNCTION CalcTVib

!-----------------------------------------------------------------------------------------------------------------------------------

REAL FUNCTION CalcTelec(MeanEelec, iSpec)
!===================================================================================================================================
! Calculation of the electronic temperature (zero-point search)
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: MeanEelec  ! Charak TVib, mean vibrational Energy of all molecules
  INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: ii
  REAL(KIND=8)                    :: LowerTemp, UpperTemp, MiddleTemp ! upper and lower value of modified zero point search
  REAL(KIND=8)                    :: eps_prec=1.0e-5   ! precision of zero point search
  REAL(KIND=8)                    :: SumOne, SumTwo    ! both summs
!===================================================================================================================================

  ! lower limit: very small value or lowest temperature if ionized
  ! upper limit: highest possible temperature
  IF ( MeanEelec .GT. 0 ) THEN
    IF ( SpecDSMC(iSpec)%ElectronicState(2,0) .EQ. 0 ) THEN
      LowerTemp = 1.0
    ELSE
      LowerTemp = SpecDSMC(iSpec)%ElectronicState(2,0)
    END IF
    UpperTemp = SpecDSMC(iSpec)%ElectronicState(2,SpecDSMC(iSpec)%MaxElecQuant-1)
    DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
      MiddleTemp = 0.5*( LowerTemp + UpperTemp)
      SumOne = 0.0
      SumTwo = 0.0
      DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant-1
        SumOne = SumOne + SpecDSMC(iSpec)%ElectronicState(1,ii) * &
                  exp( - SpecDSMC(iSpec)%ElectronicState(2,ii) / MiddleTemp )
        SumTwo = SumTwo + SpecDSMC(iSpec)%ElectronicState(1,ii) * SpecDSMC(iSpec)%ElectronicState(2,ii) * &
                  exp( - SpecDSMC(iSpec)%ElectronicState(2,ii) / MiddleTemp )
      END DO
      IF ( SumTwo / SumOne .GT. MeanEelec / BoltzmannConst ) THEN
        UpperTemp = MiddleTemp
      ELSE
        LowerTemp = MiddleTemp
      END IF
    END DO
    CalcTelec = UpperTemp ! or 0.5*( Tmax + Tmin)
  ELSE
    CalcTelec = 0. ! sup
  END IF

  RETURN

END FUNCTION CalcTelec


REAL FUNCTION CalcTVibPoly(MeanEVib, iSpec)
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for polyatomic molecules
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PolyatomMolDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
  INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: iDOF,iPolyatMole
  REAL(KIND=8)                    :: LowerTemp, UpperTemp, MiddleTemp ! upper and lower value of modified zero point search
  REAl(KIND=8)                    :: eps_prec=1.0E-5,JToEv   ! precision of zero point search
  REAL(KIND=8)                    :: SumOne    ! both summs
!===================================================================================================================================

  ! lower limit: very small value or lowest temperature if ionized
  ! upper limit: highest possible temperature
  JToEv = 1.602176565E-19  
  iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
  IF ( MeanEVib .GT. 0 ) THEN
    LowerTemp = 1.0
    UpperTemp = 5.0*SpecDSMC(iSpec)%Ediss_eV*JToEv/BoltzmannConst
    DO WHILE ( ABS( UpperTemp - LowerTemp ) .GT. eps_prec )
      MiddleTemp = 0.5*( LowerTemp + UpperTemp)
      SumOne = 0.0
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        SumOne = SumOne + 0.5*BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
              + BoltzmannConst * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) &
              / (EXP(PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)/MiddleTemp) -1.0)
      END DO
      IF ( SumOne .GT. MeanEVib) THEN
        UpperTemp = MiddleTemp
      ELSE
        LowerTemp = MiddleTemp
      END IF
    END DO
    CalcTVibPoly = UpperTemp ! or 0.5*( Tmax + Tmin)
  ELSE
    CalcTVibPoly = 0. ! sup
  END IF
  RETURN

END FUNCTION CalcTVibPoly


REAL FUNCTION CalcMeanFreePath(SpecPartNum, nPart, Volume, opt_omega, opt_temp)
!===================================================================================================================================
! Calculation of the mean free path for the hard sphere and variable hard sphere (if omega and temperature are given)
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Globals_Vars,           ONLY : Pi
  USE MOD_Particle_Vars,          ONLY : Species, nSpecies
  USE MOD_DSMC_Vars,              ONLY : SpecDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                :: Volume
  INTEGER, INTENT(IN)             :: SpecPartNum(:), nPart
  REAL, OPTIONAL, INTENT(IN)      :: opt_omega, opt_temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                         :: iSpec, jSpec
  REAL                            :: DrefMixture, omega, Temp, MFP_Tmp
!===================================================================================================================================
  DrefMixture = 0.0
  CalcMeanFreePath = 0.0

  ! Calculation of mixture reference diameter

  DO iSpec = 1, nSpecies
    DrefMixture = DrefMixture + REAL(SpecPartNum(iSpec))*SpecDSMC(iSpec)%DrefVHS / REAL(nPart)
  END DO
  ! Calculation of mean free path for a gas mixture (Bird 1986, p. 96, Eq. 4.77)
  ! (only defined for a single weighting factor, if omega is present calculation of the mean free path with the VHS model)
  IF(PRESENT(opt_omega).AND.PRESENT(opt_temp)) THEN
    omega=opt_omega
    Temp = opt_temp
    DO iSpec = 1, nSpecies
      MFP_Tmp = 0.0
      DO jSpec = 1, nSpecies
        MFP_Tmp = MFP_Tmp + (Pi*DrefMixture**2.*REAL(SpecPartNum(jSpec))*Species(jSpec)%MacroParticleFactor / Volume &
                                * (SpecDSMC(iSpec)%TrefVHS/Temp)**(omega) &
                                * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
      END DO
      CalcMeanFreePath = CalcMeanFreePath + REAL(SpecPartNum(iSpec)) / REAL(nPart) / MFP_Tmp
    END DO
  ELSE
    DO iSpec = 1, nSpecies
      MFP_Tmp = 0.0
      DO jSpec = 1, nSpecies
        MFP_Tmp = MFP_Tmp + (Pi*DrefMixture**2.*REAL(SpecPartNum(jSpec))*Species(jSpec)%MacroParticleFactor / Volume &
                                * SQRT(1+Species(iSpec)%MassIC/Species(jSpec)%MassIC))
      END DO
      CalcMeanFreePath = CalcMeanFreePath + REAL(SpecPartNum(iSpec)) / REAL(nPart) / MFP_Tmp
    END DO
  END IF
  RETURN

END FUNCTION CalcMeanFreePath



SUBROUTINE InitHODSMC()
!===================================================================================================================================
! Calculates macroscopic surface values from samples
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,          ONLY:nElems, Elem_xGP, sJ, nBCSides, SideToElem,XCL_NGeo
USE MOD_DSMC_Vars,          ONLY:DSMCSampVolWe, HODSMC,DSMCSampNearInt, DSMCSampCellVolW
USE MOD_Globals
USE MOD_ReadInTools
#ifndef MPI
USE MOD_Particle_Mesh_Vars, ONLY:GEO
#endif /*NOT MPI*/
USE MOD_PreProc,            ONLY:PP_N
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_Basis,              ONLY:LegendreGaussNodesAndWeights, LegGaussLobNodesAndWeights
USE MOD_Basis,              ONLY:BarycentricWeights,InitializeVandermonde
USE MOD_Interpolation_Vars, ONLY:xGP, wBary
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL          :: xmin, ymin, zmin, xmax, ymax, zmax
  INTEGER       :: iElem, i, ALLOCSTAT, j, k, m, l, iSide, jj,kk,mm
  REAL,ALLOCATABLE                        :: Vdm_ElemxgpN_DSMCNOut(:,:)
  REAL,ALLOCATABLE                        :: xGP_tmp(:)
  REAL, ALLOCATABLE                       :: DetJacGauss_N(:,:,:,:), DetLocal(:,:,:,:)!, Volumes(:,:,:)
  LOGICAL, ALLOCATABLE                    :: VolumeDone(:,:,:)
#ifndef MPI
  INTEGER       :: k2,m2,l2
#endif /*NOT MPI*/
!===================================================================================================================================
 
SWRITE(UNIT_stdOut,'(A)') ' INIT High Order DSMC Sampling...'

ALLOCATE( Vdm_ElemxgpN_DSMCNOut(0:HODSMC%nOutputDSMC,0:PP_N) &
          , xGP_tmp(0:HODSMC%nOutputDSMC))
ALLOCATE(HODSMC%DSMC_wGP(0:HODSMC%nOutputDSMC))

HODSMC%NodeType = GETSTR('DSMC-HOSampling-NodeType','visu')
SELECT CASE(TRIM(HODSMC%NodeType))
CASE('visu')
  DO i=0,HODSMC%nOutputDSMC
    xGP_tmp(i) = 2./REAL(HODSMC%nOutputDSMC) * REAL(i) - 1. 
    HODSMC%DSMC_wGP(i) = 2./REAL(HODSMC%nOutputDSMC)
  END DO
  HODSMC%DSMC_wGP(0) = HODSMC%DSMC_wGP(0) * 0.5
  HODSMC%DSMC_wGP(HODSMC%nOutputDSMC) = HODSMC%DSMC_wGP(HODSMC%nOutputDSMC) * 0.5
CASE('gauss')
  CALL LegendreGaussNodesAndWeights(HODSMC%nOutputDSMC,xGP_tmp,HODSMC%DSMC_wGP)
CASE('gauss-lobatto')
  CALL LegGaussLobNodesAndWeights(HODSMC%nOutputDSMC,xGP_tmp,HODSMC%DSMC_wGP)
CASE DEFAULT
  CALL abort(__STAMP__, &
      'Unknown HODSMCNodeType in dsmc_analyze.f90')
END SELECT

CALL InitializeVandermonde(PP_N,HODSMC%nOutputDSMC,wBary,xGP,xGP_tmp,Vdm_ElemxgpN_DSMCNOut)

SELECT CASE(TRIM(HODSMC%SampleType))
CASE('cartmesh_volumeweighting')
  ! ChangeBasis3D to lower or higher polynomial degree
  ALLOCATE(HODSMC%DSMC_xGP(1:3,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems))
  DO iElem=1,nElems
    CALL ChangeBasis3D(3,PP_N, HODSMC%nOutputDSMC,Vdm_ElemxgpN_DSMCNOut, Elem_xGP(:,:,:,:,iElem),HODSMC%DSMC_xGP(:,:,:,:,iElem))
  END DO ! iElem
  ! read in background mesh size
  DSMCSampVolWe%BGMdeltas(1:3) = GETREALARRAY('DSMCSampVolWe-BGMdeltas',3,'0. , 0. , 0.')
  DSMCSampVolWe%FactorBGM(1:3) = GETREALARRAY('DSMCSampVolWe-FactorBGM',3,'1. , 1. , 1.')
  DSMCSampVolWe%BGMdeltas(1:3) = 1./DSMCSampVolWe%FactorBGM(1:3)*DSMCSampVolWe%BGMdeltas(1:3)
  IF (ANY(DSMCSampVolWe%BGMdeltas.EQ.0.0)) THEN
    CALL abort(__STAMP__, &
      'ERROR: DSMCSampVolWe-BGMdeltas: No size for the cartesian background mesh definded.')
  END IF

  DSMCSampVolWe%OrderVolInt = GETINT('DSMCSampVolWe-VolIntOrd','50')
  ALLOCATE(DSMCSampVolWe%x_VolInt(0:DSMCSampVolWe%OrderVolInt),DSMCSampVolWe%w_VolInt(0:DSMCSampVolWe%OrderVolInt))
  CALL LegendreGaussNodesAndWeights(DSMCSampVolWe%OrderVolInt,DSMCSampVolWe%x_VolInt,DSMCSampVolWe%w_VolInt)  

  ! calc min and max coordinates for local mesh
  xmin = 1.0E200
  ymin = 1.0E200
  zmin = 1.0E200
  xmax = -1.0E200
  ymax = -1.0E200
  zmax = -1.0E200
  DO iElem = 1, nElems
    xmin=MIN(xmin,MINVAL(XCL_NGeo(1,:,:,:,iElem)))
    xmax=MAX(xmax,MAXVAL(XCL_NGeo(1,:,:,:,iElem)))
    ymin=MIN(ymin,MINVAL(XCL_NGeo(2,:,:,:,iElem)))
    ymax=MAX(ymax,MAXVAL(XCL_NGeo(2,:,:,:,iElem)))
    zmin=MIN(zmin,MINVAL(XCL_NGeo(3,:,:,:,iElem)))
    zmax=MAX(zmax,MAXVAL(XCL_NGeo(3,:,:,:,iElem)))
  END DO
  ! define minimum and maximum backgroundmesh index, compute volume
  DSMCSampVolWe%BGMVolume = DSMCSampVolWe%BGMdeltas(1)*DSMCSampVolWe%BGMdeltas(2)*DSMCSampVolWe%BGMdeltas(3)
  DSMCSampVolWe%BGMminX = FLOOR(xmin/DSMCSampVolWe%BGMdeltas(1)-0.0001)
  DSMCSampVolWe%BGMminY = FLOOR(ymin/DSMCSampVolWe%BGMdeltas(2)-0.0001)
  DSMCSampVolWe%BGMminZ = FLOOR(zmin/DSMCSampVolWe%BGMdeltas(3)-0.0001)
  DSMCSampVolWe%BGMmaxX = CEILING(xmax/DSMCSampVolWe%BGMdeltas(1)+0.0001)
  DSMCSampVolWe%BGMmaxY = CEILING(ymax/DSMCSampVolWe%BGMdeltas(2)+0.0001)
  DSMCSampVolWe%BGMmaxZ = CEILING(zmax/DSMCSampVolWe%BGMdeltas(3)+0.0001)

  ! mapping from gausspoints to BGM
  ALLOCATE(DSMCSampVolWe%GaussBGMIndex(1:3,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMIndex!') 
  END IF
  ALLOCATE(DSMCSampVolWe%GaussBGMFactor(1:3,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in pic_depo.f90: Cannot allocate GaussBGMFactor!')
  END IF
  DO iElem = 1, nElems
    DO j = 0, HODSMC%nOutputDSMC
      DO k = 0, HODSMC%nOutputDSMC
        DO m = 0, HODSMC%nOutputDSMC
          DSMCSampVolWe%GaussBGMIndex(1,j,k,m,iElem) = FLOOR(HODSMC%DSMC_xGP(1,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(1))
          DSMCSampVolWe%GaussBGMIndex(2,j,k,m,iElem) = FLOOR(HODSMC%DSMC_xGP(2,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(2))
          DSMCSampVolWe%GaussBGMIndex(3,j,k,m,iElem) = FLOOR(HODSMC%DSMC_xGP(3,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(3))
          DSMCSampVolWe%GaussBGMFactor(1,j,k,m,iElem) = (HODSMC%DSMC_xGP(1,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(1)) &
              -REAL(DSMCSampVolWe%GaussBGMIndex(1,j,k,m,iElem))
          DSMCSampVolWe%GaussBGMFactor(2,j,k,m,iElem) = (HODSMC%DSMC_xGP(2,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(2)) &
              -REAL(DSMCSampVolWe%GaussBGMIndex(2,j,k,m,iElem))
          DSMCSampVolWe%GaussBGMFactor(3,j,k,m,iElem) = (HODSMC%DSMC_xGP(3,j,k,m,iElem)/DSMCSampVolWe%BGMdeltas(3)) &
              -REAL(DSMCSampVolWe%GaussBGMIndex(3,j,k,m,iElem))
        END DO
      END DO
    END DO
  END DO

  ALLOCATE(DSMCSampVolWe%isBoundBGCell(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, &
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ))
  DSMCSampVolWe%isBoundBGCell = .false.
  DO iSide = 1, nBCSides
    iElem = SideToElem(S2E_ELEM_ID, iSide)
    !IF (iElem.EQ.-1) iElem = SideToElem(S2E_NB_ELEM_ID, iSide)
    DO j = 0, HODSMC%nOutputDSMC
      DO k = 0, HODSMC%nOutputDSMC
        DO m = 0, HODSMC%nOutputDSMC
         DO jj=-1,1; DO kk=-1,1; DO mm=-1,1
           IF (.NOT.(DSMCSampVolWe%GaussBGMIndex(1,j,k,m,iElem)+jj.LT.DSMCSampVolWe%BGMminX &
                .OR.DSMCSampVolWe%GaussBGMIndex(2,j,k,m,iElem)+kk.LT.DSMCSampVolWe%BGMminY &
                .OR.DSMCSampVolWe%GaussBGMIndex(3,j,k,m,iElem)+mm.LT.DSMCSampVolWe%BGMminZ)) THEN
             DSMCSampVolWe%isBoundBGCell(DSMCSampVolWe%GaussBGMIndex(1,j,k,m,iElem)+jj, &
                DSMCSampVolWe%GaussBGMIndex(2,j,k,m,iElem)+kk, DSMCSampVolWe%GaussBGMIndex(3,j,k,m,iElem)+mm) = .true.
           END IF
          END DO; END DO; END DO
        END DO
      END DO
    END DO
  END DO

  ALLOCATE(DSMCSampVolWe%BGMVolumes(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, &
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ))  
  ALLOCATE(VolumeDone(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, &
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ))
  ALLOCATE(DSMCSampVolWe%BGMVolumes2(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, &
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ))
  VolumeDone = .false.
  DSMCSampVolWe%BGMVolumes = DSMCSampVolWe%BGMVolume !0.0
  DSMCSampVolWe%BGMVolumes2  = 0.0

  DO j = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
    DO k = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
      DO m = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
        IF (DSMCSampVolWe%isBoundBGCell(j,k,m)) THEN
          CALL VolumeBoundBGMCInt(j, k, m, DSMCSampVolWe%BGMVolumes(j,k,m))
        END IF
      END DO      
    END DO
  END DO  

#ifdef MPI
  CALL MPIBackgroundMeshInitDSMCHO()
  CALL MPIVolumeExchangeBGMDSMCHO()
#else

  IF(GEO%nPeriodicVectors.GT.0)THEN
    ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
    ALLOCATE(DSMCSampVolWe%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(__STAMP__, &
      'ERROR in MPIBackgroundMeshInitDSMCHO: cannot allocate DSMCSampVolWe%PeriodicBGMVectors!')
    END IF
    DO i = 1, GEO%nPeriodicVectors
      DSMCSampVolWe%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/DSMCSampVolWe%BGMdeltas(1))
      IF(ABS(GEO%PeriodicVectors(1,i)/DSMCSampVolWe%BGMdeltas(1)-REAL(DSMCSampVolWe%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      DSMCSampVolWe%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/DSMCSampVolWe%BGMdeltas(2))
      IF(ABS(GEO%PeriodicVectors(2,i)/DSMCSampVolWe%BGMdeltas(2)-REAL(DSMCSampVolWe%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
      DSMCSampVolWe%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/DSMCSampVolWe%BGMdeltas(3))
      IF(ABS(GEO%PeriodicVectors(3,i)/DSMCSampVolWe%BGMdeltas(3)-REAL(DSMCSampVolWe%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
        CALL abort(__STAMP__, &
      'ERROR: Periodic Vector ist not multiple of background mesh delta')
      END IF
    END DO
  END IF
  DO i = 1,GEO%nPeriodicVectors
    DO k = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
      k2 = k + DSMCSampVolWe%PeriodicBGMVectors(1,i)
      DO l = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
        l2 = l + DSMCSampVolWe%PeriodicBGMVectors(2,i)
        DO m = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
          m2 = m + DSMCSampVolWe%PeriodicBGMVectors(3,i)
          IF ((k2.GE.DSMCSampVolWe%BGMminX).AND.(k2.LE.DSMCSampVolWe%BGMmaxX)) THEN
            IF ((l2.GE.DSMCSampVolWe%BGMminY).AND.(l2.LE.DSMCSampVolWe%BGMmaxY)) THEN
              IF ((m2.GE.DSMCSampVolWe%BGMminZ).AND.(m2.LE.DSMCSampVolWe%BGMmaxZ)) THEN
                DSMCSampVolWe%BGMVolumes(k,l,m) = DSMCSampVolWe%BGMVolumes(k,l,m) + DSMCSampVolWe%BGMVolumes(k2,l2,m2)
                DSMCSampVolWe%BGMVolumes(k2,l2,m2) = DSMCSampVolWe%BGMVolumes(k,l,m)
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
#endif
  DEALLOCATE(HODSMC%DSMC_xGP,HODSMC%DSMC_wGP)
CASE('nearest_gausspoint')
  
  ALLOCATE(HODSMC%sJ(0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,1:nElems))
  ALLOCATE( DetJacGauss_N(1,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC) &
          , DetLocal(1,0:PP_N,0:PP_N,0:PP_N))
  DO iElem=1,nElems
    DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
      DetLocal(1,j,k,l)=1./sJ(j,k,l,iElem)
    END DO; END DO; END DO
    CALL ChangeBasis3D(1,PP_N, HODSMC%nOutputDSMC,Vdm_ElemxgpN_DSMCNOut, DetLocal(:,:,:,:),DetJacGauss_N(:,:,:,:))
    DO j=0, HODSMC%nOutputDSMC; DO k=0, HODSMC%nOutputDSMC; DO l=0, HODSMC%nOutputDSMC
      HODSMC%sJ(j,k,l,iElem)=1./DetJacGauss_N(1,j,k,l)
    END DO; END DO; END DO
  END DO ! iElem    
  ! compute the borders of the virtual volumes around the gauss points in -1|1 space
  ALLOCATE(DSMCSampNearInt%GaussBorder(1:HODSMC%nOutputDSMC),STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) THEN
    CALL abort(__STAMP__, &
      'ERROR in dsmc_analyze.f90: Cannot allocate Mapped Gauss Border Coords!')
  END IF
  DO i = 1,HODSMC%nOutputDSMC
    DSMCSampNearInt%GaussBorder(i) = (xGP_tmp(i) + xGP_tmp(i-1))/2
  END DO
CASE('cell_mean')    
  DEALLOCATE(HODSMC%DSMC_wGP)
CASE('cell_volweight')    
  ALLOCATE(DSMCSampCellVolW%xGP(0:HODSMC%nOutputDSMC))
  DSMCSampCellVolW%xGP(0:HODSMC%nOutputDSMC) = xGP_tmp(0:HODSMC%nOutputDSMC)
  DSMCSampCellVolW%xGP(0:HODSMC%nOutputDSMC) = (DSMCSampCellVolW%xGP(0:HODSMC%nOutputDSMC)+1.0)/2.0
  DEALLOCATE(HODSMC%DSMC_wGP)
CASE DEFAULT
  CALL abort(__STAMP__, &
      'Unknown DSMCHOSampleType in dsmc_analyze.f90')
END SELECT

END SUBROUTINE InitHODSMC



SUBROUTINE DSMCHO_data_sampling()
!===================================================================================================================================
! Sampling of variables velocity and energy for DSMC
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY:PartStateIntEn, DSMCSampVolWe, DSMC, CollisMode, SpecDSMC, HODSMC, DSMC_HOSolution
  USE MOD_DSMC_Vars,              ONLY:DSMCSampNearInt, DSMCSampCellVolW
  USE MOD_Particle_Vars,          ONLY:PartState, PDM, PartSpecies, Species, nSpecies, PEM,PartPosRef
  USE MOD_Mesh_Vars,              ONLY:nElems
  USE MOD_Particle_Mesh_Vars,     ONLY:Geo
  USE MOD_Particle_Tracking_vars, ONLY:DoRefMapping
  USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
  !USE MOD_part_MPFtools,          ONLY:GeoCoordToMap
  USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iElem, iLoopx, iLoopy, iLoopz, k, l, m, i, kk, ll, mm, iSpec, a, b, ii
  REAL, ALLOCATABLE             :: BGMSource(:,:,:,:,:), alphaSum(:,:,:,:),BGMSourceCellVol(:,:,:,:,:,:)
  REAL, ALLOCATABLE             :: alphaSumCellVol(:,:,:,:,:), Source(:,:,:,:,:,:)
  REAL                          :: alpha1, alpha2, alpha3, TSource(1:11)
!===================================================================================================================================
DSMC%SampNum = DSMC%SampNum + 1
SELECT CASE(TRIM(HODSMC%SampleType))
 CASE('cartmesh_volumeweighting')
  ! Step 1: Deposition of all particles onto background mesh -> densities
  ALLOCATE(BGMSource(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, &
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ,1:11, 1:nSpecies), &
          alphaSum(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY:DSMCSampVolWe%BGMmaxY, & 
          DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ, 1:nSpecies))

  BGMSource(:,:,:,:,:) = 0.0
  alphaSum(:,:,:,:) = 0.0
  DO i = 1, PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      iSpec = PartSpecies(i)
      k = FLOOR(PartState(i,1)/DSMCSampVolWe%BGMdeltas(1))
      l = FLOOR(PartState(i,2)/DSMCSampVolWe%BGMdeltas(2))
      m = FLOOR(PartState(i,3)/DSMCSampVolWe%BGMdeltas(3))
      alpha1 = (PartState(i,1) / DSMCSampVolWe%BGMdeltas(1)) - k
      alpha2 = (PartState(i,2) / DSMCSampVolWe%BGMdeltas(2)) - l
      alpha3 = (PartState(i,3) / DSMCSampVolWe%BGMdeltas(3)) - m
      TSource(:) = 0.0
      TSource(1:3) = PartState(i,4:6)
      TSource(4:6) = PartState(i,4:6)**2
      TSource(7) = 1.0  !density
      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(PartSpecies(i))%InterID.EQ.2)) THEN
        TSource(8:9)      =  PartStateIntEn(i,1:2)
      ELSE
        TSource(8:9) = 0.0
      END IF
      IF (DSMC%ElectronicState) THEN
        TSource(10)     =  PartStateIntEn(i,3)
      ELSE
        TSource(10) = 0.0
      END IF
      TSource(11) = 1.0 

      BGMSource(k,l,m,1:11,iSpec)       = BGMSource(k,l,m,1:11,iSpec) + (TSource(1:11) * (1-alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k,l,m+1,1:11,iSpec)     = BGMSource(k,l,m+1,1:11,iSpec) + (TSource(1:11) * (1-alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k,l+1,m,1:11,iSpec)     = BGMSource(k,l+1,m,1:11,iSpec) + (TSource(1:11) * (1-alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k,l+1,m+1,1:11,iSpec)   = BGMSource(k,l+1,m+1,1:11,iSpec) + (TSource(1:11) * (1-alpha1)*(alpha2)*(alpha3))
      BGMSource(k+1,l,m,1:11,iSpec)     = BGMSource(k+1,l,m,1:11,iSpec) + (TSource(1:11) * (alpha1)*(1-alpha2)*(1-alpha3))
      BGMSource(k+1,l,m+1,1:11,iSpec)   = BGMSource(k+1,l,m+1,1:11,iSpec) + (TSource(1:11) * (alpha1)*(1-alpha2)*(alpha3))
      BGMSource(k+1,l+1,m,1:11,iSpec)   = BGMSource(k+1,l+1,m,1:11,iSpec) + (TSource(1:11) * (alpha1)*(alpha2)*(1-alpha3))
      BGMSource(k+1,l+1,m+1,1:11,iSpec) = BGMSource(k+1,l+1,m+1,1:11,iSpec) + (TSource(1:11) * (alpha1)*(alpha2)*(alpha3))

      alphaSum(k,l,m,iSpec)       = alphaSum(k,l,m,iSpec) + (1-alpha1)*(1-alpha2)*(1-alpha3)
      alphaSum(k,l,m+1,iSpec)     = alphaSum(k,l,m+1,iSpec) + (1-alpha1)*(1-alpha2)*(alpha3)
      alphaSum(k,l+1,m,iSpec)     = alphaSum(k,l+1,m,iSpec) + (1-alpha1)*(alpha2)*(1-alpha3)
      alphaSum(k,l+1,m+1,iSpec)   = alphaSum(k,l+1,m+1,iSpec) + (1-alpha1)*(alpha2)*(alpha3)
      alphaSum(k+1,l,m,iSpec)     = alphaSum(k+1,l,m,iSpec) + (alpha1)*(1-alpha2)*(1-alpha3)
      alphaSum(k+1,l,m+1,iSpec)   = alphaSum(k+1,l,m+1,iSpec) + (alpha1)*(1-alpha2)*(alpha3)
      alphaSum(k+1,l+1,m,iSpec)   = alphaSum(k+1,l+1,m,iSpec) + (alpha1)*(alpha2)*(1-alpha3)
      alphaSum(k+1,l+1,m+1,iSpec) = alphaSum(k+1,l+1,m+1,iSpec) + (alpha1)*(alpha2)*(alpha3)
     END IF
  END DO

#ifdef MPI
  CALL MPISourceExchangeBGMDSMCHO(BGMSource, alphaSum)
#else
  IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchangeDSMCHO(BGMSource, alphaSum)
#endif

  DO iSpec = 1, nSpecies
    DO iLoopx = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
      DO iLoopy = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
        DO iLoopz = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
          IF (alphaSum(iLoopx, iLoopy, iLoopz,iSpec).GT.0.0) THEN        
            BGMSource(iLoopx,iLoopy,iLoopz,1:6,iSpec) = BGMSource(iLoopx,iLoopy,iLoopz,1:6,iSpec) &
                    / alphaSum(iLoopx,iLoopy,iLoopz,iSpec)
            BGMSource(iLoopx,iLoopy,iLoopz,8:10,iSpec) = BGMSource(iLoopx,iLoopy,iLoopz,8:10,iSpec) &
                   / alphaSum(iLoopx,iLoopy,iLoopz,iSpec)
          ELSE
            BGMSource(iLoopx,iLoopy,iLoopz,1:6,iSpec) = 0.0
            BGMSource(iLoopx,iLoopy,iLoopz,8:10,iSpec) = 0.0
          END IF
          IF (DSMCSampVolWe%BGMVolumes(iLoopx,iLoopy,iLoopz).GT.0.0) THEN
              BGMSource(iLoopx,iLoopy,iLoopz,7,iSpec) = BGMSource(iLoopx,iLoopy,iLoopz,7,iSpec) &
                / DSMCSampVolWe%BGMVolumes(iLoopx,iLoopy,iLoopz)* Species(iSpec)%MacroParticleFactor   
          ELSE
              BGMSource(iLoopx,iLoopy,iLoopz,7,iSpec) = 0.0
          END IF
        END DO
      END DO
    END DO
  END DO

  ! Step 2: Interpolation of densities onto grid
  DO iSpec = 1, nSpecies
    DO iElem = 1, nElems
      DO kk = 0, HODSMC%nOutputDSMC
        DO ll = 0, HODSMC%nOutputDSMC
          DO mm = 0, HODSMC%nOutputDSMC
           k = DSMCSampVolWe%GaussBGMIndex(1,kk,ll,mm,iElem)
           l = DSMCSampVolWe%GaussBGMIndex(2,kk,ll,mm,iElem)
           m = DSMCSampVolWe%GaussBGMIndex(3,kk,ll,mm,iElem)
           alpha1 = DSMCSampVolWe%GaussBGMFactor(1,kk,ll,mm,iElem)
           alpha2 = DSMCSampVolWe%GaussBGMFactor(2,kk,ll,mm,iElem)
           alpha3 = DSMCSampVolWe%GaussBGMFactor(3,kk,ll,mm,iElem)     
           DSMC_HOSolution(:,kk,ll,mm,iElem,iSpec) = (DSMC_HOSolution(:,kk,ll,mm,iElem, iSpec) * (REAL(DSMC%SampNum) - 1.0) & 
             +  (BGMSource(k,l,m,:,iSpec) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k,l,m+1,:,iSpec) * (1-alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k,l+1,m,:,iSpec) * (1-alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k,l+1,m+1,:,iSpec) * (1-alpha1) * (alpha2) * (alpha3) + &
                BGMSource(k+1,l,m,:,iSpec) * (alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSource(k+1,l,m+1,:,iSpec) * (alpha1) * (1-alpha2) * (alpha3) + &
                BGMSource(k+1,l+1,m,:,iSpec) * (alpha1) * (alpha2) * (1-alpha3) + &
                BGMSource(k+1,l+1,m+1,:,iSpec) * (alpha1) * (alpha2) * (alpha3))) / REAL(DSMC%SampNum)
         END DO !mm
       END DO !ll
     END DO !kk
   END DO !iElem
 END DO
 DEALLOCATE(BGMSource, alphasum)
CASE('nearest_gausspoint')
  ALLOCATE(Source(1:11,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC, &
          0:HODSMC%nOutputDSMC,1:nElems, 1:nSpecies))
  IF(MOD(HODSMC%nOutputDSMC,2).EQ.0) THEN
    a = HODSMC%nOutputDSMC/2
    b = a
  ELSE
    a = (HODSMC%nOutputDSMC+1)/2
    b = a-1
  END IF
  Source=0.0
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      iSpec = PartSpecies(i)
      iElem = PEM%Element(i)
      ! Map Particle to -1|1 space (re-used in interpolation)
      ! check with depositions and PartPosRef already mapped
      IF(.NOT.DoRefMapping)THEN
        CALL Eval_xyz_ElemCheck(PartState(i,1:3),PartPosRef(1:3,i),iElem,i)
      END IF
      !CALL GeoCoordToMap(PartState(i,1:3),PartPosRef(1:3),iElem)
      ! Find out which gausspoint is closest and add up charges and currents
      !! x-direction
      k = a
      DO ii = 0,b-1
        IF(ABS(PartPosRef(1,i)).GE.DSMCSampNearInt%GaussBorder(HODSMC%nOutputDSMC-ii))THEN
          k = HODSMC%nOutputDSMC-ii
          EXIT
        END IF
      END DO
      k = NINT((HODSMC%nOutputDSMC+SIGN(2.0*k-HODSMC%nOutputDSMC,PartPosRef(1,i)))/2)
      !! y-direction
      l = a
      DO ii = 0,b-1
        IF(ABS(PartPosRef(2,i)).GE.DSMCSampNearInt%GaussBorder(HODSMC%nOutputDSMC-ii))THEN
          l = HODSMC%nOutputDSMC-ii
          EXIT
        END IF
      END DO
      l = NINT((HODSMC%nOutputDSMC+SIGN(2.0*l-HODSMC%nOutputDSMC,PartPosRef(2,i)))/2)
      !! z-direction
      m = a
      DO ii = 0,b-1
        IF(ABS(PartPosRef(3,i)).GE.DSMCSampNearInt%GaussBorder(HODSMC%nOutputDSMC-ii))THEN
          m = HODSMC%nOutputDSMC-ii
          EXIT
        END IF
      END DO
      m = NINT((HODSMC%nOutputDSMC+SIGN(2.0*m-HODSMC%nOutputDSMC,PartPosRef(3,i)))/2)
      Source(1:3,k,l,m,iElem, iSpec) = Source(1:3,k,l,m,iElem, iSpec) + PartState(i,4:6)
      Source(4:6,k,l,m,iElem, iSpec) = Source(4:6,k,l,m,iElem, iSpec) + PartState(i,4:6)**2
      Source(7,k,l,m,iElem, iSpec) = Source(7,k,l,m,iElem, iSpec) + 1.0  !density
      IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(PartSpecies(i))%InterID.EQ.2)) THEN
        Source(8:9,k,l,m,iElem, iSpec) = Source(8:9,k,l,m,iElem, iSpec) + PartStateIntEn(i,1:2)
      END IF
      IF (DSMC%ElectronicState) THEN
        Source(10,k,l,m,iElem, iSpec) = Source(10,k,l,m,iElem, iSpec) + PartStateIntEn(i,3)
      END IF
      Source(11,k,l,m,iElem, iSpec) = Source(11,k,l,m,iElem, iSpec) + 1.0 
    END IF
  END DO
  DSMC_HOSolution(:,:,:,:,:,:) = (DSMC_HOSolution(:,:,:,:,:,:) * (REAL(DSMC%SampNum) - 1.0) &
        + Source(:,:,:,:,:,:))/REAL(DSMC%SampNum)
CASE('cell_mean')
  ALLOCATE(Source(1:11,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC, &
          0:HODSMC%nOutputDSMC,1:nElems, 1:nSpecies))
  Source=0.0
  DO i=1,PDM%ParticleVecLength
    IF (PDM%ParticleInside(i)) THEN
      iSpec = PartSpecies(i)
      iElem = PEM%Element(i)
      DO kk = 0, HODSMC%nOutputDSMC
      DO ll = 0, HODSMC%nOutputDSMC
      DO mm = 0, HODSMC%nOutputDSMC
        Source(1:3,kk,ll,mm,iElem, iSpec) = Source(1:3,kk,ll,mm,iElem, iSpec) + PartState(i,4:6)
        Source(4:6,kk,ll,mm,iElem, iSpec) = Source(4:6,kk,ll,mm,iElem, iSpec) + PartState(i,4:6)**2
        Source(7,kk,ll,mm,iElem, iSpec) = Source(7,kk,ll,mm,iElem, iSpec) + 1.0  !density
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(i))%InterID.EQ.2)) THEN
          Source(8:9,kk,ll,mm,iElem, iSpec) = Source(8:9,kk,ll,mm,iElem, iSpec) + PartStateIntEn(i,1:2)
        END IF
        IF (DSMC%ElectronicState) THEN
          Source(10,kk,ll,mm,iElem, iSpec) = Source(10,kk,ll,mm,iElem, iSpec) + PartStateIntEn(i,3)
        END IF
        Source(11,kk,ll,mm,iElem, iSpec) = Source(11,kk,ll,mm,iElem, iSpec) + 1.0 
      END DO; END DO; END DO
    END IF
  END DO
  DSMC_HOSolution(:,:,:,:,:,:) = (DSMC_HOSolution(:,:,:,:,:,:) * (REAL(DSMC%SampNum) - 1.0) &
        + Source(:,:,:,:,:,:))/REAL(DSMC%SampNum)
CASE('cell_volweight') 
  ALLOCATE(BGMSourceCellVol(0:1,0:1,0:1,1:nElems,1:11, 1:nSpecies), &
          alphaSumCellVol(0:1,0:1,0:1,1:nElems, 1:nSpecies))
  BGMSourceCellVol(:,:,:,:,:,:) = 0.0
  alphaSumCellVol(:,:,:,:,:) = 0.0

  DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN  
    iElem = PEM%Element(iPart)
    iSpec = PartSpecies(iPart)
    IF(.NOT.DoRefMapping)THEN
      CALL Eval_xyz_ElemCheck(PartState(iPart,1:3),PartPosRef(1:3,iPart),iElem,iPart)
    END IF
    !CALL GeoCoordToMap(PartState(iPart,1:3), TempPartPos(1:3), iElem)
    TSource(:) = 0.0
    TSource(1:3) = PartState(iPart,4:6)
    TSource(4:6) = PartState(iPart,4:6)**2
    TSource(7) = 1.0  !density
    IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
        (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
      TSource(8:9)      =  PartStateIntEn(iPart,1:2)
    ELSE
      TSource(8:9) = 0.0
    END IF
    IF (DSMC%ElectronicState) THEN
      TSource(10)     =  PartStateIntEn(iPart,3)
    ELSE
      TSource(10) = 0.0
    END IF
    TSource(11) = 1.0 
    alpha1=(PartPosRef(1,iPart)+1.0)/2.0
    alpha2=(PartPosRef(2,iPart)+1.0)/2.0
    alpha3=(PartPosRef(3,iPart)+1.0)/2.0
    BGMSourceCellVol(0,0,0,iElem,1:11,iSpec) = BGMSourceCellVol(0,0,0,iElem,1:11,iSpec) &
            + (TSource(1:11) * (1-alpha1)*(1-alpha2)*(1-alpha3))
    BGMSourceCellVol(0,0,1,iElem,1:11,iSpec) = BGMSourceCellVol(0,0,1,iElem,1:11,iSpec) &
            + (TSource(1:11) * (1-alpha1)*(1-alpha2)*(alpha3))
    BGMSourceCellVol(0,1,0,iElem,1:11,iSpec) = BGMSourceCellVol(0,1,0,iElem,1:11,iSpec) &
            + (TSource(1:11) * (1-alpha1)*(alpha2)*(1-alpha3))
    BGMSourceCellVol(0,1,1,iElem,1:11,iSpec) = BGMSourceCellVol(0,1,1,iElem,1:11,iSpec) &
            + (TSource(1:11) * (1-alpha1)*(alpha2)*(alpha3))
    BGMSourceCellVol(1,0,0,iElem,1:11,iSpec) = BGMSourceCellVol(1,0,0,iElem,1:11,iSpec) &
            + (TSource(1:11) * (alpha1)*(1-alpha2)*(1-alpha3))
    BGMSourceCellVol(1,0,1,iElem,1:11,iSpec) = BGMSourceCellVol(1,0,1,iElem,1:11,iSpec) &
            + (TSource(1:11) * (alpha1)*(1-alpha2)*(alpha3))
    BGMSourceCellVol(1,1,0,iElem,1:11,iSpec) = BGMSourceCellVol(1,1,0,iElem,1:11,iSpec) &
            + (TSource(1:11) * (alpha1)*(alpha2)*(1-alpha3))
    BGMSourceCellVol(1,1,1,iElem,1:11,iSpec) = BGMSourceCellVol(1,1,1,iElem,1:11,iSpec) &
            + (TSource(1:11) * (alpha1)*(alpha2)*(alpha3))

    alphaSumCellVol(0,0,0,iElem,iSpec) = alphaSumCellVol(0,0,0,iElem,iSpec) + (1-alpha1)*(1-alpha2)*(1-alpha3)
    alphaSumCellVol(0,0,1,iElem,iSpec) = alphaSumCellVol(0,0,1,iElem,iSpec) + (1-alpha1)*(1-alpha2)*(alpha3)
    alphaSumCellVol(0,1,0,iElem,iSpec) = alphaSumCellVol(0,1,0,iElem,iSpec) + (1-alpha1)*(alpha2)*(1-alpha3)
    alphaSumCellVol(0,1,1,iElem,iSpec) = alphaSumCellVol(0,1,1,iElem,iSpec) + (1-alpha1)*(alpha2)*(alpha3)
    alphaSumCellVol(1,0,0,iElem,iSpec) = alphaSumCellVol(1,0,0,iElem,iSpec) + (alpha1)*(1-alpha2)*(1-alpha3)
    alphaSumCellVol(1,0,1,iElem,iSpec) = alphaSumCellVol(1,0,1,iElem,iSpec) + (alpha1)*(1-alpha2)*(alpha3)
    alphaSumCellVol(1,1,0,iElem,iSpec) = alphaSumCellVol(1,1,0,iElem,iSpec) + (alpha1)*(alpha2)*(1-alpha3)
    alphaSumCellVol(1,1,1,iElem,iSpec) = alphaSumCellVol(1,1,1,iElem,iSpec) + (alpha1)*(alpha2)*(alpha3)   
  END IF
  END DO
  
  DO iSpec = 1, nSpecies
    DO iElem=1, nElems
      DO iLoopx = 0,1
        DO iLoopy = 0,1
          DO iLoopz = 0,1
            IF (alphaSumCellVol(iLoopx, iLoopy, iLoopz,iElem,iSpec).GT.0.0) THEN        
              BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,1:6,iSpec) = BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,1:6,iSpec) &
                      / alphaSumCellVol(iLoopx,iLoopy,iLoopz,iElem,iSpec)
              BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,8:10,iSpec) = BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,8:10,iSpec) &
                     / alphaSumCellVol(iLoopx,iLoopy,iLoopz,iElem,iSpec)
            ELSE
              BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,1:6,iSpec) = 0.0
              BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,8:10,iSpec) = 0.0
            END IF
            BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,7,iSpec) = BGMSourceCellVol(iLoopx,iLoopy,iLoopz,iElem,7,iSpec) &
                / GEO%Volume(iElem) * Species(iSpec)%MacroParticleFactor   
          END DO
        END DO
      END DO
    END DO
  END DO

  DO iSpec = 1, nSpecies
    DO iElem = 1, nElems
      DO kk = 0, HODSMC%nOutputDSMC
        DO ll = 0, HODSMC%nOutputDSMC
          DO mm = 0, HODSMC%nOutputDSMC
           alpha1 = DSMCSampCellVolW%xGP(kk)
           alpha2 = DSMCSampCellVolW%xGP(ll)
           alpha3 = DSMCSampCellVolW%xGP(mm)
           DSMC_HOSolution(:,kk,ll,mm,iElem,iSpec) = (DSMC_HOSolution(:,kk,ll,mm,iElem, iSpec) * (REAL(DSMC%SampNum) - 1.0) & 
             +  (BGMSourceCellVol(0,0,0,iElem,:,iSpec) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSourceCellVol(0,0,1,iElem,:,iSpec) * (1-alpha1) * (1-alpha2) * (alpha3) + &
                BGMSourceCellVol(0,1,0,iElem,:,iSpec) * (1-alpha1) * (alpha2) * (1-alpha3) + &
                BGMSourceCellVol(0,1,1,iElem,:,iSpec) * (1-alpha1) * (alpha2) * (alpha3) + &
                BGMSourceCellVol(1,0,0,iElem,:,iSpec) * (alpha1) * (1-alpha2) * (1-alpha3) + &
                BGMSourceCellVol(1,0,1,iElem,:,iSpec) * (alpha1) * (1-alpha2) * (alpha3) + &
                BGMSourceCellVol(1,1,0,iElem,:,iSpec) * (alpha1) * (alpha2) * (1-alpha3) + &
                BGMSourceCellVol(1,1,1,iElem,:,iSpec) * (alpha1) * (alpha2) * (alpha3))) / REAL(DSMC%SampNum)
         END DO !mm
       END DO !ll
     END DO !kk
   END DO !iElem
 END DO
 DEALLOCATE(BGMSourceCellVol, alphaSumCellVol)
CASE DEFAULT
 CALL abort(__STAMP__, &
      'Unknown DepositionType in pic_depo.f90')
END SELECT
END SUBROUTINE DSMCHO_data_sampling


SUBROUTINE WriteDSMCHOToHDF5(MeshFileName,OutputTime, FutureTime)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,            ONLY:HODSMC, DSMC_HOSolution, CollisMode, SpecDSMC, DSMC
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,          ONLY:ProjectName
USE MOD_Mesh_Vars,            ONLY:offsetElem,nGlobalElems, nElems
USE MOD_io_HDF5
USE MOD_HDF5_output,          ONLY: WriteArrayToHDF5
USE MOD_Particle_Vars,        ONLY: Species, BoltzmannConst, nSpecies
USE MOD_Particle_Mesh_Vars,   ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVal, iElem, kk , ll, mm, iSpec
REAL,ALLOCATABLE               :: DSMC_MacroVal(:,:,:,:,:), TVib_TempFac(:,:,:)
#ifdef MPI
REAL                           :: StartT,EndT
#endif
!===================================================================================================================================
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE DSMC-HO TO HDF5 FILE...'
#ifdef MPI
  StartT=MPI_WTIME()
#endif

ALLOCATE(DSMC_MacroVal(1:11,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,nElems))
ALLOCATE(TVib_TempFac(0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC,0:HODSMC%nOutputDSMC))
! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:11))
StrVarNames(1)='VeloX'
StrVarNames(2)='VeloY'
StrVarNames(3)='VeloZ'
StrVarNames(4)='TempX'
StrVarNames(5)='TempY'
StrVarNames(6)='TempZ'
StrVarNames(7)='Density'       
StrVarNames(8)='TVib'
StrVarNames(9)='TRot'
StrVarNames(10)='TElec'
StrVarNames(11)='PointWeight'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_DSMCHOState',OutputTime))//'.h5'
! PO:
! excahnge PP_N through Nout
IF(MPIRoot) CALL GenerateDSMCHOFileSkeleton('DSMCHOState',11,StrVarNames,MeshFileName,OutputTime,FutureTime)
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

#ifdef MPI
CALL OpenDataFile(FileName,create=.false.,single=.FALSE.)
#else
CALL OpenDataFile(FileName,create=.false.)
#endif

! Write DG solution ----------------------------------------------------------------------------------------------------------------
nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
DO iSpec = 1, nSpecies
  DSMC_MacroVal = 0.0
  DO iElem = 1, nElems ! element/cell main loop    
    DO kk = 0, HODSMC%nOutputDSMC; DO ll = 0, HODSMC%nOutputDSMC; DO mm = 0, HODSMC%nOutputDSMC
      SELECT CASE(TRIM(HODSMC%SampleType))
      CASE('cartmesh_volumeweighting','cell_volweight')
        ! compute flow velocity
        DSMC_MacroVal(1:3,kk,ll,mm, iElem) = DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec)
        ! compute flow Temperature
        DSMC_MacroVal(4:6,kk,ll,mm, iElem) = Species(iSpec)%MassIC/ BoltzmannConst &
                            * (DSMC_HOSolution(4:6,kk,ll,mm, iElem, iSpec)- DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec)**2)
        ! compute density
        DSMC_MacroVal(7,kk,ll,mm, iElem) = DSMC_HOSolution(7,kk,ll,mm, iElem, iSpec)
          !       if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
        !      IF (usevMPF) THEN
        !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
        !      ELSE 
        !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * &
        !         Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
        !      END IF
        ! compute internal energies / has to be changed for vfd 
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
            TVib_TempFac(kk,ll,mm)=DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec)/ (BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
            IF (TVib_TempFac(kk,ll,mm).LE.DSMC%GammaQuant) THEN
              DSMC_MacroVal(8,kk,ll,mm, iElem) = 0.0           
            ELSE
              DSMC_MacroVal(8,kk,ll,mm, iElem) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(TVib_TempFac(kk,ll,mm)-DSMC%GammaQuant))
            END IF
          ELSE                                            ! TSHO-model
            DSMC_MacroVal(8,kk,ll,mm, iElem)  = CalcTVib(SpecDSMC(iSpec)%CharaTVib & 
                , DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec), SpecDSMC(iSpec)%MaxVibQuant) 
          END IF       
          DSMC_MacroVal(9,kk,ll,mm, iElem) = DSMC_HOSolution(9,kk,ll,mm, iElem, iSpec)/(BoltzmannConst)
          IF (DSMC%ElectronicState) THEN
            DSMC_MacroVal(10,kk,ll,mm, iElem)= CalcTelec( DSMC_HOSolution(10,kk,ll,mm, iElem, iSpec), iSpec)
          END IF
        END IF
        DSMC_MacroVal(11,kk,ll,mm, iElem) = DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec) 
      CASE('nearest_gausspoint') 
        IF (DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec).GT.0.0) THEN
          ! compute flow velocity
          DSMC_MacroVal(1:3,kk,ll,mm, iElem) = DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec) &
                /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)
          ! compute flow Temperature
          DSMC_MacroVal(4:6,kk,ll,mm, iElem) = Species(iSpec)%MassIC/ BoltzmannConst &
                              * (DSMC_HOSolution(4:6,kk,ll,mm, iElem, iSpec) /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec) &
                            - (DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec)/DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec))**2)
          ! compute density
          DSMC_MacroVal(7,kk,ll,mm, iElem) = DSMC_HOSolution(7,kk,ll,mm, iElem, iSpec)*HODSMC%sJ(kk,ll,mm,iElem) &
            /(HODSMC%DSMC_wGP(kk)*HODSMC%DSMC_wGP(ll)*HODSMC%DSMC_wGP(mm))*Species(iSpec)%MacroParticleFactor
            !       if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
          !      IF (usevMPF) THEN
          !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
          !      ELSE 
          !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * &
          !         Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
          !      END IF
          ! compute internal energies / has to be changed for vfd 
          IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
            IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
              TVib_TempFac(kk,ll,mm)=DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec)/ (DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec) &
                    *BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
              IF (TVib_TempFac(kk,ll,mm).LE.DSMC%GammaQuant) THEN
                DSMC_MacroVal(8,kk,ll,mm, iElem) = 0.0           
              ELSE
                DSMC_MacroVal(8,kk,ll,mm, iElem) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(TVib_TempFac(kk,ll,mm)-DSMC%GammaQuant))
              END IF
            ELSE                                            ! TSHO-model
              DSMC_MacroVal(8,kk,ll,mm, iElem)  = CalcTVib(SpecDSMC(iSpec)%CharaTVib & 
                  , DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec) &
                  /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec),SpecDSMC(iSpec)%MaxVibQuant)
            END IF       
            DSMC_MacroVal(9,kk,ll,mm, iElem) = DSMC_HOSolution(9,kk,ll,mm, iElem, iSpec) &
               /(DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)*BoltzmannConst)
            IF (DSMC%ElectronicState) THEN
              DSMC_MacroVal(10,kk,ll,mm, iElem)= CalcTelec( DSMC_HOSolution(10,kk,ll,mm, iElem, iSpec)&
                  /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec), iSpec)
            END IF
          END IF
        ELSE
          DSMC_MacroVal(1:10,kk,ll,mm, iElem) = 0.0
        END IF
        DSMC_MacroVal(11,kk,ll,mm, iElem) = DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)    
      CASE('cell_mean')         
        IF (DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec).GT.0.0) THEN
          ! compute flow velocity
          DSMC_MacroVal(1:3,kk,ll,mm, iElem) = DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec) &
                /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)
          ! compute flow Temperature
          DSMC_MacroVal(4:6,kk,ll,mm, iElem) = Species(iSpec)%MassIC/ BoltzmannConst &
                              * (DSMC_HOSolution(4:6,kk,ll,mm, iElem, iSpec) /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec) &
                            - (DSMC_HOSolution(1:3,kk,ll,mm, iElem, iSpec)/DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec))**2)
          ! compute density
          DSMC_MacroVal(7,kk,ll,mm, iElem) = DSMC_HOSolution(7,kk,ll,mm, iElem, iSpec) &
            /GEO%Volume(iElem)*Species(iSpec)%MacroParticleFactor
          !       if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
          !      IF (usevMPF) THEN
          !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
          !      ELSE 
          !        MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * &
          !         Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
          !      END IF
          ! compute internal energies / has to be changed for vfd 
          IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
            IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
              TVib_TempFac(kk,ll,mm)=DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec)/ (DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec) &
                    *BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
              IF (TVib_TempFac(kk,ll,mm).LE.DSMC%GammaQuant) THEN
                DSMC_MacroVal(8,kk,ll,mm, iElem) = 0.0           
              ELSE
                DSMC_MacroVal(8,kk,ll,mm, iElem) = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(TVib_TempFac(kk,ll,mm)-DSMC%GammaQuant))
              END IF
            ELSE                                            ! TSHO-model
              DSMC_MacroVal(8,kk,ll,mm, iElem)  = CalcTVib(SpecDSMC(iSpec)%CharaTVib & 
                  , DSMC_HOSolution(8,kk,ll,mm, iElem, iSpec)&
                  /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec),SpecDSMC(iSpec)%MaxVibQuant)
            END IF       
            DSMC_MacroVal(9,kk,ll,mm, iElem) = DSMC_HOSolution(9,kk,ll,mm, iElem, iSpec) &
               /(DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)*BoltzmannConst)
            IF (DSMC%ElectronicState) THEN
              DSMC_MacroVal(10,kk,ll,mm, iElem)= CalcTelec( DSMC_HOSolution(10,kk,ll,mm, iElem, iSpec)&
                  /DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec), iSpec)
            END IF
          END IF
        ELSE
          DSMC_MacroVal(1:10,kk,ll,mm, iElem) = 0.0
        END IF
        DSMC_MacroVal(11,kk,ll,mm, iElem) = DSMC_HOSolution(11,kk,ll,mm, iElem, iSpec)              
      END SELECT
    END DO; END DO; END DO
  END DO
      
  WRITE(H5_Name,'(A,I3.3)') 'DSMCHO_Spec',iSpec
  CALL WriteArrayToHDF5(DataSetName=H5_Name, rank=5,&
                      nValGlobal=(/11,HODSMC%nOutputDSMC+1,HODSMC%nOutputDSMC+1,HODSMC%nOutputDSMC+1,nGlobalElems/),&
                      nVal=      (/11,HODSMC%nOutputDSMC+1,HODSMC%nOutputDSMC+1,HODSMC%nOutputDSMC+1,PP_nElems/),&
                      offset=    (/0,      0,     0,     0,     offsetElem/),&
                      collective=.false., existing=.FALSE., RealArray=DSMC_MacroVal)
END DO

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
#ifdef MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
END SUBROUTINE WriteDSMCHOToHDF5


SUBROUTINE GenerateDSMCHOFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY: ProjectName
USE MOD_ReadInTools,ONLY: GetParameters
USE MOD_PreProcFlags
USE MOD_io_HDF5
USE MOD_DSMC_Vars,            ONLY:HODSMC
USE MOD_HDF5_Output,       ONLY: WriteAttributeToHDF5, WriteHDF5Header
USE MOD_Particle_Vars,      ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
INTEGER,INTENT(IN)             :: nVar
CHARACTER(LEN=255)             :: StrVarNames(nVar)
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,MeshFile255
CHARACTER(LEN=255),ALLOCATABLE :: params(:)
CHARACTER(LEN=255)             :: NodeTypeTemp
!===================================================================================================================================
! Create file
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
#ifdef MPI
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.)
#else
CALL OpenDataFile(TRIM(FileName),create=.TRUE.)
#endif

SELECT CASE(TRIM(HODSMC%NodeType))
CASE('visu')
  NodeTypeTemp = 'VISU'
CASE('gauss')
  NodeTypeTemp = 'GAUSS'
CASE('gauss-lobatto')
  NodeTypeTemp = 'GAUSS-LOBATTO'
CASE DEFAULT
  CALL abort(__STAMP__, &
      'Unknown HODSMCNodeType in dsmc_analyze.f90')
END SELECT

! Write file header
CALL WriteHDF5Header(TRIM(TypeString),File_ID)

! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=N)
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
IF(PRESENT(FutureTime))THEN
  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
END IF
!CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeTypeTemp/))
CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)

CALL WriteAttributeToHDF5(File_ID,'NComputation',1,IntegerScalar=HODSMC%nOutputDSMC)
CALL WriteAttributeToHDF5(File_ID,'NSpecies',1,IntegerScalar=nSpecies)

! Write ini file parameters and compile flags
CALL GetParameters(params)
CALL WriteAttributeToHDF5(File_ID,'Parameters',SIZE(params),StrArray=params)
CALL WriteAttributeToHDF5(File_ID,'Compile',1,StrScalar=(/PREPROC_FLAGS/))
DEALLOCATE(params)

CALL CloseDataFile()
END SUBROUTINE GenerateDSMCHOFileSkeleton


#ifndef MPI
SUBROUTINE PeriodicSourceExchangeDSMCHO(BGMSource, alphasum)    
!============================================================================================================================
! Exchange sources in periodic case
!============================================================================================================================
! use MODULES                                                    
  USE MOD_Particle_MPI_Vars, ONLY:PartMPI
  !USE MOD_PICDepo_Vars
  USE MOD_Particle_Mesh_Vars, ONLY:Geo
  USE MOD_Particle_Vars
  USE MOD_DSMC_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL,INTENT(INOUT)         :: BGMSource(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY &
                            :DSMCSampVolWe%BGMmaxY,DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ,1:11, 1:nSpecies)
    REAL,INTENT(INOUT)         :: alphasum(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY &
                            :DSMCSampVolWe%BGMmaxY,DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ, 1:nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                           
    INTEGER                     :: i,k,l,m,k2,l2,m2
!-----------------------------------------------------------------------------------------------------------------------------------

DO i = 1,GEO%nPeriodicVectors
  DO k = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
    k2 = k + DSMCSampVolWe%PeriodicBGMVectors(1,i)
    DO l = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
      l2 = l + DSMCSampVolWe%PeriodicBGMVectors(2,i)
      DO m = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
        m2 = m + DSMCSampVolWe%PeriodicBGMVectors(3,i)
        IF ((k2.GE.DSMCSampVolWe%BGMminX).AND.(k2.LE.DSMCSampVolWe%BGMmaxX)) THEN
          IF ((l2.GE.DSMCSampVolWe%BGMminY).AND.(l2.LE.DSMCSampVolWe%BGMmaxY)) THEN
            IF ((m2.GE.DSMCSampVolWe%BGMminZ).AND.(m2.LE.DSMCSampVolWe%BGMmaxZ)) THEN
              BGMSource(k,l,m,:,:) = BGMSource(k,l,m,:,:) + BGMSource(k2,l2,m2,:,:)
              BGMSource(k2,l2,m2,:,:) = BGMSource(k,l,m,:,:)
              alphasum(k,l,m,:) = alphasum(k,l,m,:) + alphasum(k2,l2,m2,:)
              alphasum(k2,l2,m2,:) = alphasum(k,l,m,:)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO
RETURN
END SUBROUTINE PeriodicSourceExchangeDSMCHO
#else /*MPI*/
SUBROUTINE MPISourceExchangeBGMDSMCHO(BGMSource, alphasum)            
!============================================================================================================================
! Exchange sources in periodic case for MPI
!============================================================================================================================
! use MODULES                                                         
  USE MOD_Globals
  USE MOD_Particle_Vars
  USE MOD_Particle_Mesh_Vars, ONLY: GEO
  USE MOD_Particle_MPI_Vars,  ONLY: PartMPI,tMPIMessage
  USE MOD_DSMC_Vars,          ONLY: DSMCSampVolWe
!-----------------------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL,INTENT(INOUT)        :: BGMSource(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY &
                          :DSMCSampVolWe%BGMmaxY,DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ,1:11, 1:nSpecies)
    REAL,INTENT(INOUT)         :: alphasum(DSMCSampVolWe%BGMminX:DSMCSampVolWe%BGMmaxX,DSMCSampVolWe%BGMminY &
                            :DSMCSampVolWe%BGMmaxY,DSMCSampVolWe%BGMminZ:DSMCSampVolWe%BGMmaxZ, 1:nSpecies)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
    TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)                              
    TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)                              
    INTEGER                     :: send_request(0:PartMPI%nProcs-1)                              
    INTEGER                     :: recv_request(0:PartMPI%nProcs-1)                              
    INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
    INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
    INTEGER                     :: i,k,l,m,n, ppp, Counter        
    INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)                              
    INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)                                
    INTEGER                     :: allocStat, Counter2                                     
    INTEGER                     :: messageCounterS, messageCounterR                                
    INTEGER                     :: myRealKind, k2,l2,m2, iSpec
    REAL                        :: myRealTestValue                                                
!-----------------------------------------------------------------------------------------------------------------------------------

     myRealKind = KIND(myRealTestValue)
     IF (myRealKind.EQ.4) THEN
       myRealKind = MPI_REAL
     ELSE IF (myRealKind.EQ.8) THEN
       myRealKind = MPI_DOUBLE_PRECISION
     ELSE
       myRealKind = MPI_REAL
     END IF
         
   
    !--- Assemble actual sources to send
     DO i = 0,PartMPI%nProcs-1
      ! sourcelength muss noch für periodisch angepasst werden
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.&
            (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          SourceLength(i)=(DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1) + 1) &
            * (DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2) + 1) &
            * (DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3) + 1)
          IF(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
            DO k = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
              SourceLength(i) = SourceLength(i) + (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,1) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
                  (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,2) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
                  (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,3) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
            END DO
          END IF
          ALLOCATE(send_message(i)%content(1:SourceLength(i)*12*nSpecies), STAT=allocStat)
          IF (allocStat .NE. 0) THEN
             CALL abort(__STAMP__, &
                'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
          END IF
       END IF
       Counter = 0
       Counter2 = 0
       IF (DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor) THEN
          DO k = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1)
          DO l = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2)
          DO m = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3)             
            Counter2 = Counter2 + 1
            DO iSpec = 1, nSpecies
              DO n = 1,11
                 send_message(i)%content((Counter2-1)*11*nSpecies +n + (iSpec-1)*11) = BGMSource(k,l,m,n,iSpec)
              END DO
              send_message(i)%content(SourceLength(i)*11*nSpecies + (Counter2-1)*nSpecies+iSpec) = alphasum(k,l,m,iSpec)
            END DO
          END DO
          END DO
          END DO
       END IF  
       IF (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
          DO n = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
             DO k = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,1),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter2 = Counter2 + 1
               DO iSpec = 1, nSpecies
                 DO ppp = 1,11
                    send_message(i)%content((Counter2-1)*11*nSpecies +ppp +(iSpec-1)*11) = BGMSource(k,l,m,ppp,iSpec)
                 END DO
                 send_message(i)%content(SourceLength(i)*11*nSpecies + (Counter2-1)*nSpecies+iSpec) = alphasum(k,l,m,iSpec)
               END DO
             END DO
             END DO
             END DO
          END DO
       END IF
    END DO

    !--- allocate actual source receive buffer
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          Counter = SourceLength(i)
          RecvLength(i) = Counter
          ALLOCATE(recv_message(i)%content(1:Counter*12*nSpecies), STAT=allocStat)
          IF (allocStat .NE. 0) THEN
             CALL abort(__STAMP__, &
               'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
          END IF
       END IF
    END DO
    !--- communicate
    messageCounterS = 0
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.&
            (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          ! MPI_ISEND true/false list for all border BGM points
          messageCounterS = messageCounterS + 1
          CALL MPI_ISEND(send_message(i)%content,SourceLength(i)*12*nSpecies,myRealKind,i,1,PartMPI%COMM, &
                         send_request(messageCounterS), IERROR)
       END IF
    END DO
    messageCounterR = 0
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.&
            (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
          messageCounterR = messageCounterR + 1
          CALL MPI_IRECV(recv_message(i)%content,RecvLength(i)*12*nSpecies,myRealKind,i,1,PartMPI%COMM, &
                         recv_request(messageCounterR), IERROR)
       END IF
    END DO
    ! MPI_WAITALL for the non-blocking MPI-communication to be finished
    IF (messageCounterS .GE. 1) THEN
       CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
    END IF
    IF (messageCounterR .GE. 1) THEN
       CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
    END IF
    !--- Deallocate Send Message Buffers
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
           DEALLOCATE(send_message(i)%content, STAT=allocStat)
           IF (allocStat .NE. 0) THEN
              CALL abort(__STAMP__, &
                'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
           END IF
       END IF
    END DO

    !--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
    !---           received sources have been added!
    IF ((GEO%nPeriodicVectors.GT.0).AND.(DSMCSampVolWe%SelfPeriodic)) THEN
       DO i = 1, GEO%nPeriodicVectors
          DO k = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
             k2 = k + DSMCSampVolWe%PeriodicBGMVectors(1,i)
          DO l = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
             l2 = l + DSMCSampVolWe%PeriodicBGMVectors(2,i)
          DO m = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
             m2 = m + DSMCSampVolWe%PeriodicBGMVectors(3,i)
             IF ((k2.GE.DSMCSampVolWe%BGMminX).AND.(k2.LE.DSMCSampVolWe%BGMmaxX)) THEN
             IF ((l2.GE.DSMCSampVolWe%BGMminY).AND.(l2.LE.DSMCSampVolWe%BGMmaxY)) THEN
             IF ((m2.GE.DSMCSampVolWe%BGMminZ).AND.(m2.LE.DSMCSampVolWe%BGMmaxZ)) THEN
                BGMSource(k,l,m,:,:) = BGMSource(k,l,m,:,:) + BGMSource(k2,l2,m2,:,:)
                BGMSource(k2,l2,m2,:,:) = BGMSource(k,l,m,:,:)
                alphasum(k,l,m,:) = alphasum(k,l,m,:) + alphasum(k2,l2,m2,:)
                alphasum(k2,l2,m2,:) = alphasum(k,l,m,:)
             END IF
             END IF
             END IF
          END DO
          END DO
          END DO
       END DO
    END IF

    !--- Add Sources and Deallocate Receive Message Buffers
    DO i = 0,PartMPI%nProcs-1
       IF (RecvLength(i).GT.0) THEN
          Counter = 0
          Counter2 = 0
          IF (DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor) THEN
             DO k = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1)
             DO l = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2)
             DO m = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3)
               Counter2 = Counter2 + 1
               DO iSpec =1, nSpecies
                 DO n = 1,11
                   BGMSource(k,l,m,n, iSpec) = BGMSource(k,l,m,n,iSpec)  &
                      + recv_message(i)%content((Counter2-1)*11*nSpecies+n +(iSpec-1)*11)
                 END DO
                 alphasum(k,l,m,iSpec) = alphasum(k,l,m,iSpec) & 
                      + recv_message(i)%content(SourceLength(i)*11*nSpecies + (Counter2-1)*nSpecies+iSpec)
               END DO
             END DO
             END DO
             END DO
          END IF
          IF (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
             DO n = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
                DO k = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,1),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,1)
                DO l = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,2),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,2)
                DO m = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,3),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter2 = Counter2 + 1
                  DO iSpec = 1, nSpecies
                    DO ppp = 1,11
                       BGMSource(k,l,m,ppp,iSpec) = BGMSource(k,l,m,ppp, iSpec) &
                          + recv_message(i)%content((Counter2-1)*11*nSpecies+ppp+(iSpec-1)*11)
                    END DO
                    alphasum(k,l,m,iSpec) = alphasum(k,l,m,iSpec) & 
                      + recv_message(i)%content(SourceLength(i)*11*nSpecies + (Counter2-1)*nSpecies+iSpec)
                  END DO
                END DO
                END DO
                END DO
             END DO
          END IF
          IF ((DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor).OR.(DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor)) THEN
             DEALLOCATE(recv_message(i)%content, STAT=allocStat)
             IF (allocStat .NE. 0) THEN
                CALL abort(__STAMP__, &
                  'ERROR in MPISourceExchangeBGMDSMCHO: cannot deallocate recv_message')
             END IF
          END IF
       END IF       
    END DO
END SUBROUTINE MPISourceExchangeBGMDSMCHO



SUBROUTINE MPIVolumeExchangeBGMDSMCHO()            
!============================================================================================================================
! Exchange sources in periodic case for MPI
!============================================================================================================================
! use MODULES                                                         
  USE MOD_Particle_MPI_Vars, ONLY: PartMPI,tMPIMessage
  USE MOD_Particle_Mesh_Vars,ONLY: GEO
  USE MOD_Globals
  USE MOD_Particle_Vars
  USE MOD_DSMC_Vars,        ONLY: DSMCSampVolWe
!-----------------------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES

!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
    TYPE(tMPIMessage)           :: send_message(0:PartMPI%nProcs-1)                              
    TYPE(tMPIMessage)           :: recv_message(0:PartMPI%nProcs-1)                              
    INTEGER                     :: send_request(0:PartMPI%nProcs-1)                              
    INTEGER                     :: recv_request(0:PartMPI%nProcs-1)                              
    INTEGER                     :: send_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
    INTEGER                     :: recv_status_list(1:MPI_STATUS_SIZE,0:PartMPI%nProcs-1)        
    INTEGER                     :: i,k,l,m,n, Counter
    INTEGER                     :: SourceLength(0:PartMPI%nProcs-1)                              
    INTEGER                     :: RecvLength(0:PartMPI%nProcs-1)                                
    INTEGER                     :: allocStat, Counter2                                     
    INTEGER                     :: messageCounterS, messageCounterR                                
    INTEGER                     :: myRealKind, k2,l2,m2
    REAL                        :: myRealTestValue                                                
!-----------------------------------------------------------------------------------------------------------------------------------

     myRealKind = KIND(myRealTestValue)
     IF (myRealKind.EQ.4) THEN
       myRealKind = MPI_REAL
     ELSE IF (myRealKind.EQ.8) THEN
       myRealKind = MPI_DOUBLE_PRECISION
     ELSE
       myRealKind = MPI_REAL
     END IF
         
        !--- Assemble actual sources to send
     DO i = 0,PartMPI%nProcs-1
       IF (DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor.OR.DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
          DO k = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1)
          DO l = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2)
          DO m = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3) 
            IF (.NOT.(DSMCSampVolWe%isBoundBGCell(k,l,m))) THEN
              CALL VolumeBoundBGMCInt(k, l, m, DSMCSampVolWe%BGMVolumes(k,l,m))      
              DSMCSampVolWe%isBoundBGCell(k,l,m) = .true.    
            END IF
          END DO
          END DO
          END DO
       END IF  
    END DO
   
    !--- Assemble actual sources to send
     DO i = 0,PartMPI%nProcs-1
      ! sourcelength muss noch für periodisch angepasst werden
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
          SourceLength(i)=(DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1) + 1) &
            * (DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2) + 1) &
            * (DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3) - DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3) + 1)         
          IF(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
            DO k = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
              SourceLength(i) = SourceLength(i) + (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,1) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,1) + 1) * &
                  (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,2) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,2) + 1) * &
                  (DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(2,3) -&
                   DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1,3) + 1)
            END DO
          END IF
          ALLOCATE(send_message(i)%content(1:SourceLength(i)), STAT=allocStat)
          IF (allocStat .NE. 0) THEN
             CALL abort(__STAMP__, &
                'ERROR in MPISourceExchangeBGM: cannot allocate send_message')
          END IF
       END IF
       Counter = 0
       Counter2 = 0
       IF (DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor) THEN
          DO k = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1)
          DO l = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2)
          DO m = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3)             
            Counter2 = Counter2 + 1
            send_message(i)%content(Counter2) = DSMCSampVolWe%BGMVolumes(k,l,m)
          END DO
          END DO
          END DO
       END IF  
       IF (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
          DO n = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
             DO k = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,1),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,1)
             DO l = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,2),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,2)
             DO m = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,3),&
                    DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,3)
               Counter2 = Counter2 + 1
               send_message(i)%content(Counter2) = DSMCSampVolWe%BGMVolumes(k,l,m)
             END DO
             END DO
             END DO
          END DO
       END IF
    END DO

    !--- allocate actual source receive buffer
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          Counter = SourceLength(i)
          RecvLength(i) = Counter
          ALLOCATE(recv_message(i)%content(1:Counter), STAT=allocStat)
          IF (allocStat .NE. 0) THEN
             CALL abort(__STAMP__, &
               'ERROR in MPISourceExchangeBGM: cannot allocate recv_message')
          END IF
       END IF
    END DO
    !--- communicate
    messageCounterS = 0
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.&
            (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          ! MPI_ISEND true/false list for all border BGM points
          messageCounterS = messageCounterS + 1
          CALL MPI_ISEND(send_message(i)%content,SourceLength(i),myRealKind,i,1,PartMPI%COMM, &
                         send_request(messageCounterS), IERROR)
       END IF
    END DO
    messageCounterR = 0
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.&
            (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
          ! MPI_IRECV true/false list for all border BGM points from neighbor CPUs
          messageCounterR = messageCounterR + 1
          CALL MPI_IRECV(recv_message(i)%content,RecvLength(i),myRealKind,i,1,PartMPI%COMM, &
                         recv_request(messageCounterR), IERROR)
       END IF
    END DO
    ! MPI_WAITALL for the non-blocking MPI-communication to be finished
    IF (messageCounterS .GE. 1) THEN
       CALL MPI_WAITALL(messageCounterS,send_request(1:messageCounterS),send_status_list(:,1:messageCounterS),IERROR)
    END IF
    IF (messageCounterR .GE. 1) THEN
       CALL MPI_WAITALL(messageCounterR,recv_request(1:messageCounterR),recv_status_list(:,1:messageCounterR),IERROR)
    END IF
    !--- Deallocate Send Message Buffers
    DO i = 0,PartMPI%nProcs-1
       IF ((DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor).OR.(DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
           DEALLOCATE(send_message(i)%content, STAT=allocStat)
           IF (allocStat .NE. 0) THEN
              CALL abort(__STAMP__, &
                'ERROR in MPISourceExchangeBGM: cannot deallocate send_message')
           END IF
       END IF
    END DO

    !--- add selfperiodic sources, if any (needs to be done after send message is compiled and before
    !---           received sources have been added!
    IF ((GEO%nPeriodicVectors.GT.0).AND.(DSMCSampVolWe%SelfPeriodic)) THEN
       DO i = 1, GEO%nPeriodicVectors
          DO k = DSMCSampVolWe%BGMminX, DSMCSampVolWe%BGMmaxX
             k2 = k + DSMCSampVolWe%PeriodicBGMVectors(1,i)
          DO l = DSMCSampVolWe%BGMminY, DSMCSampVolWe%BGMmaxY
             l2 = l + DSMCSampVolWe%PeriodicBGMVectors(2,i)
          DO m = DSMCSampVolWe%BGMminZ, DSMCSampVolWe%BGMmaxZ
             m2 = m + DSMCSampVolWe%PeriodicBGMVectors(3,i)
             IF ((k2.GE.DSMCSampVolWe%BGMminX).AND.(k2.LE.DSMCSampVolWe%BGMmaxX)) THEN
             IF ((l2.GE.DSMCSampVolWe%BGMminY).AND.(l2.LE.DSMCSampVolWe%BGMmaxY)) THEN
             IF ((m2.GE.DSMCSampVolWe%BGMminZ).AND.(m2.LE.DSMCSampVolWe%BGMmaxZ)) THEN
                DSMCSampVolWe%BGMVolumes(k,l,m) = DSMCSampVolWe%BGMVolumes(k,l,m) + DSMCSampVolWe%BGMVolumes(k2,l2,m2)
                DSMCSampVolWe%BGMVolumes(k2,l2,m2) = DSMCSampVolWe%BGMVolumes(k,l,m)             
             END IF
             END IF
             END IF
          END DO
          END DO
          END DO
       END DO
    END IF

    !--- Add Sources and Deallocate Receive Message Buffers
    DO i = 0,PartMPI%nProcs-1
       IF (RecvLength(i).GT.0) THEN
          Counter = 0
          Counter2 = 0
          IF (DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor) THEN
             DO k = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,1), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,1)
             DO l = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,2), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,2)
             DO m = DSMCSampVolWe%MPIConnect(i)%BGMBorder(1,3), DSMCSampVolWe%MPIConnect(i)%BGMBorder(2,3)
               Counter2 = Counter2 + 1
               DSMCSampVolWe%BGMVolumes(k,l,m) = DSMCSampVolWe%BGMVolumes(k,l,m)  &
                   + recv_message(i)%content(Counter2)
             END DO
             END DO
             END DO
          END IF       
          IF (DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor) THEN
             DO n = 1, DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount
                DO k = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,1),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,1)
                DO l = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,2),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,2)
                DO m = DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(1,3),&
                       DSMCSampVolWe%MPIConnect(i)%Periodic(n)%BGMPeriodicBorder(2,3)
                  Counter2 = Counter2 + 1
                  DSMCSampVolWe%BGMVolumes(k,l,m) = DSMCSampVolWe%BGMVolumes(k,l,m)  &
                   + recv_message(i)%content(Counter2)
                END DO
                END DO
                END DO
             END DO
          END IF   
          IF ((DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor)) THEN
             DEALLOCATE(recv_message(i)%content, STAT=allocStat)
             IF (allocStat .NE. 0) THEN
                CALL abort(__STAMP__, &
                  'ERROR in MPISourceExchangeBGMDSMCHO: cannot deallocate recv_message')
             END IF
          END IF
       END IF       
    END DO
END SUBROUTINE MPIVolumeExchangeBGMDSMCHO

SUBROUTINE MPIBackgroundMeshInitDSMCHO()  
!============================================================================================================================
! initialize MPI background mesh
!============================================================================================================================
! use MODULES          
  USE MOD_Particle_Vars    
  USE MOD_Globals
  USE MOD_Particle_MPI_Vars,  ONLY: PartMPI
  USE MOD_DSMC_Vars,          ONLY: DSMCSampVolWe
  USE MOD_Particle_Mesh_Vars, ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE                                                                                  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
    INTEGER                     :: i,k,m,n                                    
    INTEGER                     :: localminmax(6), maxofmin, minofmax                              
    INTEGER                     :: completeminmax(6*PartMPI%nProcs)                              
    INTEGER                     :: allocStat, NeighCount                                                                 
    INTEGER                     :: TempBorder(1:2,1:3)            
    INTEGER                     :: Periodicminmax(6), coord, PeriodicVec(1:3)                   
    INTEGER                     :: TempPeriBord(1:26,1:2,1:3)                                      
    LOGICAL                     :: CHECKNEIGHBOR     
!-----------------------------------------------------------------------------------------------------------------------------------
 ! Periodic Init stuff
 IF(GEO%nPeriodicVectors.GT.0)THEN
   ! Compute PeriodicBGMVectors (from PeriodicVectors and BGMdeltas)
   ALLOCATE(DSMCSampVolWe%PeriodicBGMVectors(1:3,1:GEO%nPeriodicVectors),STAT=allocStat)
   IF (allocStat .NE. 0) THEN
     CALL abort(__STAMP__, &
        'ERROR in MPIBackgroundMeshInitDSMCHO: cannot allocate DSMCSampVolWe%PeriodicBGMVectors!')
   END IF
   DO i = 1, GEO%nPeriodicVectors
     DSMCSampVolWe%PeriodicBGMVectors(1,i) = NINT(GEO%PeriodicVectors(1,i)/DSMCSampVolWe%BGMdeltas(1))
     IF(ABS(GEO%PeriodicVectors(1,i)/DSMCSampVolWe%BGMdeltas(1)-REAL(DSMCSampVolWe%PeriodicBGMVectors(1,i))).GT.1E-10)THEN
       CALL abort(__STAMP__, &
        'ERROR: Periodic Vector ist not multiple of background mesh delta')
     END IF
     DSMCSampVolWe%PeriodicBGMVectors(2,i) = NINT(GEO%PeriodicVectors(2,i)/DSMCSampVolWe%BGMdeltas(2))
     IF(ABS(GEO%PeriodicVectors(2,i)/DSMCSampVolWe%BGMdeltas(2)-REAL(DSMCSampVolWe%PeriodicBGMVectors(2,i))).GT.1E-10)THEN
       CALL abort(__STAMP__, &
        'ERROR: Periodic Vector ist not multiple of background mesh delta')
     END IF
     DSMCSampVolWe%PeriodicBGMVectors(3,i) = NINT(GEO%PeriodicVectors(3,i)/DSMCSampVolWe%BGMdeltas(3))
     IF(ABS(GEO%PeriodicVectors(3,i)/DSMCSampVolWe%BGMdeltas(3)-REAL(DSMCSampVolWe%PeriodicBGMVectors(3,i))).GT.1E-10)THEN
       CALL abort(__STAMP__, &
        'ERROR: Periodic Vector ist not multiple of background mesh delta')
     END IF
   END DO
   ! Check whether process is periodic with itself
   DSMCSampVolWe%SelfPeriodic = .FALSE.
   !--- virtually move myself according to periodic vectors in order to find overlapping areas
   !--- 26 possibilities,
   localminmax(1) = DSMCSampVolWe%BGMminX
   localminmax(2) = DSMCSampVolWe%BGMminY
   localminmax(3) = DSMCSampVolWe%BGMminZ
   localminmax(4) = DSMCSampVolWe%BGMmaxX
   localminmax(5) = DSMCSampVolWe%BGMmaxY
   localminmax(6) = DSMCSampVolWe%BGMmaxZ
   DO k = -1,1
     DO m = -1,1
       DO n = -1,1
         PeriodicVec = k*DSMCSampVolWe%PeriodicBGMVectors(:,1) + m*DSMCSampVolWe%PeriodicBGMVectors(:,1) & 
            + n*DSMCSampVolWe%PeriodicBGMVectors(:,1)
         IF (ALL(PeriodicVec(:).EQ.0)) CYCLE
         periodicminmax(1) = localminmax(1) + PeriodicVec(1)
         periodicminmax(2) = localminmax(2) + PeriodicVec(2)
         periodicminmax(3) = localminmax(3) + PeriodicVec(3)
         periodicminmax(4) = localminmax(4) + PeriodicVec(1)
         periodicminmax(5) = localminmax(5) + PeriodicVec(2)
         periodicminmax(6) = localminmax(6) + PeriodicVec(3)
         !--- find overlap
         DO coord = 1,3           ! x y z direction
           maxofmin = MAX(periodicminmax(coord),localminmax(coord))
           minofmax = MIN(periodicminmax(3+coord),localminmax(3+coord))
           IF (maxofmin.LE.minofmax) DSMCSampVolWe%SelfPeriodic = .TRUE.      ! overlapping
         END DO
       END DO
     END DO
   END DO
 END IF

 !--- send and receive min max indices to and from all processes

    !--- enter local min max vector (xmin, ymin, zmin, xmax, ymax, zmax)
    localminmax(1) = DSMCSampVolWe%BGMminX
    localminmax(2) = DSMCSampVolWe%BGMminY
    localminmax(3) = DSMCSampVolWe%BGMminZ
    localminmax(4) = DSMCSampVolWe%BGMmaxX
    localminmax(5) = DSMCSampVolWe%BGMmaxY
    localminmax(6) = DSMCSampVolWe%BGMmaxZ
    !--- do allgather into complete min max vector
    CALL MPI_ALLGATHER(localminmax,6,MPI_INTEGER,completeminmax,6,MPI_INTEGER,PartMPI%COMM,IERROR)
    ! Allocate MPIConnect
    SDEALLOCATE(DSMCSampVolWe%MPIConnect)
    ALLOCATE(DSMCSampVolWe%MPIConnect(0:PartMPI%nProcs-1),STAT=allocStat)
    IF (allocStat .NE. 0) THEN
      CALL abort(__STAMP__, &
        'ERROR in MPIBackgroundMeshInit: cannot allocate DSMCSampVolWe%MPIConnect')
    END IF

    !--- determine borders indices (=overlapping BGM mesh points) with each process
    DO i = 0,PartMPI%nProcs-1
      DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor = .FALSE.
      DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount = 0
       IF (i.EQ.PartMPI%MyRank) THEN
          DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor = .FALSE.
       ELSE
          DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor = .TRUE.
          DO k = 1,3           ! x y z direction
             maxofmin = MAX(localminmax(k),completeminmax((i*6)+k))
             minofmax = MIN(localminmax(3+k),completeminmax((i*6)+3+k))
             IF (maxofmin.LE.minofmax) THEN           ! overlapping
                TempBorder(1,k) = maxofmin
                TempBorder(2,k) = minofmax
             ELSE
                DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor = .FALSE.
             END IF
          END DO          
       END IF
       IF(DSMCSampVolWe%MPIConnect(i)%isBGMNeighbor)THEN
          SDEALLOCATE(DSMCSampVolWe%MPIConnect(i)%BGMBorder)
          ALLOCATE(DSMCSampVolWe%MPIConnect(i)%BGMBorder(1:2,1:3),STAT=allocStat)
          IF (allocStat .NE. 0) THEN
             CALL abort(__STAMP__, &
              'ERROR in MPIBackgroundMeshInit: cannot allocate DSMCSampVolWe%MPIConnect')
          END IF
          DSMCSampVolWe%MPIConnect(i)%BGMBorder(1:2,1:3) = TempBorder(1:2,1:3)
       END IF
    END DO

    
    !--- determine border indices for periodic meshes  
    IF (GEO%nPeriodicVectors.GT.0) THEN   
      DO i = 0,PartMPI%nProcs-1
        IF (i.EQ.PartMPI%MyRank) THEN
          DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor = .FALSE.
          DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount = 0
        ELSE
          !--- virtually move myself according to periodic vectors in order to find overlapping areas
          !--- 26 possibilities, processes need to work through them in opposite direction in order
          !--- to get matching areas.
          !--- Example for 2D:  I am process #3, I compare myself with #7
          !--- Periodic Vectors are p1 and p2.
          !--- I check p1, p2, p1+p2, p1-p2, -p1+p2, -p1-p2, -p2, -p1
          !--- #7 has to check -p1, -p2, -p1-p2, -p1+p2, p1-p2, p1+p1, p2, p1
          !--- This is done by doing 3 loops from -1 to 1 (for the higher process number)
          !--- or 1 to -1 (for the lower process number) and multiplying
          !--- these numbers to the periodic vectors
          NeighCount = 0   !-- counter: how often is the process my periodic neighbor?
          DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor = .FALSE.
          DO k = -SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i)
            DO m = -SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i)
              DO n = -SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i),SIGN(1,PartMPI%MyRank-i)
                IF ((k.EQ.0).AND.(m.EQ.0).AND.(n.EQ.0)) CYCLE !this is not periodic and already done above
                CHECKNEIGHBOR = .TRUE.
                PeriodicVec = k*DSMCSampVolWe%PeriodicBGMVectors(:,1)
                IF (GEO%nPeriodicVectors.GT.1) THEN
                  PeriodicVec = PeriodicVec + m*DSMCSampVolWe%PeriodicBGMVectors(:,2)
                END IF
                IF (GEO%nPeriodicVectors.GT.2) THEN
                  PeriodicVec = PeriodicVec + n*DSMCSampVolWe%PeriodicBGMVectors(:,3)
                END IF
                periodicminmax(1) = localminmax(1) + PeriodicVec(1)
                periodicminmax(2) = localminmax(2) + PeriodicVec(2)
                periodicminmax(3) = localminmax(3) + PeriodicVec(3)
                periodicminmax(4) = localminmax(4) + PeriodicVec(1)
                periodicminmax(5) = localminmax(5) + PeriodicVec(2)
                periodicminmax(6) = localminmax(6) + PeriodicVec(3)
                !--- find overlap
                DO coord = 1,3           ! x y z direction
                  maxofmin = MAX(periodicminmax(coord),completeminmax((i*6)+coord))
                  minofmax = MIN(periodicminmax(3+coord),completeminmax((i*6)+3+coord))
                  IF (maxofmin.LE.minofmax) THEN           ! overlapping
                    TempBorder(1,coord) = maxofmin
                    TempBorder(2,coord) = minofmax
                  ELSE
                    CHECKNEIGHBOR = .FALSE.
                  END IF
                END DO
                IF(CHECKNEIGHBOR)THEN
                  NeighCount = NeighCount + 1
                  TempBorder(:,1) = TempBorder(:,1) - PeriodicVec(1)
                  TempBorder(:,2) = TempBorder(:,2) - PeriodicVec(2)
                  TempBorder(:,3) = TempBorder(:,3) - PeriodicVec(3)
                  TempPeriBord(NeighCount,1:2,1:3) = TempBorder(1:2,1:3)
                  DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor = .TRUE.
                END IF
              END DO
            END DO
          END DO
          DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount = NeighCount
          ALLOCATE(DSMCSampVolWe%MPIConnect(i)%Periodic(1:DSMCSampVolWe%MPIConnect(i)%BGMPeriodicBorderCount),STAT=allocStat)
          IF (allocStat .NE. 0) THEN
            CALL abort(__STAMP__,  &
              'ERROR in MPIBackgroundMeshInit: cannot allocate DSMCSampVolWe%MPIConnect')
          END IF
          DO k = 1,NeighCount
            ALLOCATE(DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1:2,1:3),STAT=allocStat)
            IF (allocStat .NE. 0) THEN
              CALL abort(__STAMP__, &
                'ERROR in MPIBackgroundMeshInit: cannot allocate DSMCSampVolWe%MPIConnect')
            END IF
            DSMCSampVolWe%MPIConnect(i)%Periodic(k)%BGMPeriodicBorder(1:2,1:3) = TempPeriBord(k,1:2,1:3)
          END DO
        END IF
      END DO
    ELSE
      !--- initialize to FALSE for completely non-periodic cases
      DO i = 0,PartMPI%nProcs-1
        DSMCSampVolWe%MPIConnect(i)%isBGMPeriodicNeighbor = .FALSE.
      END DO
    END IF
  RETURN
END SUBROUTINE MPIBackgroundMeshInitDSMCHO
#endif /*MPI*/

SUBROUTINE VolumeBoundBGMCInt(i, j, k, Volume)                                                         !
  USE MOD_Particle_Vars
  USE MOD_Particle_Mesh_Vars,     ONLY:GEO,epsOneCell
  USE MOD_DSMC_Vars,              ONLY:DSMCSampVolWe
  USE MOD_Eval_xyz,               ONLY:eval_xyz_elemcheck
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
   INTEGER                          :: i,j, k, Element, CellX,CellY,CellZ, iElem
   INTEGER                          :: stepx, stepy, stepz
   REAL                             :: xi(3)
   REAL                             :: GuessPos(3), Volume, Found
   REAL       :: alpha1, alpha2, alpha3
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
Found = 0

DO stepx=0, DSMCSampVolWe%OrderVolInt
DO stepy=0, DSMCSampVolWe%OrderVolInt
DO stepz=0, DSMCSampVolWe%OrderVolInt
  GuessPos(1) = DSMCSampVolWe%BGMdeltas(1)*(i) + DSMCSampVolWe%x_VolInt(stepx)*DSMCSampVolWe%BGMdeltas(1)
  GuessPos(2) = DSMCSampVolWe%BGMdeltas(2)*(j) + DSMCSampVolWe%x_VolInt(stepy)*DSMCSampVolWe%BGMdeltas(2)
  GuessPos(3) = DSMCSampVolWe%BGMdeltas(3)*(k) + DSMCSampVolWe%x_VolInt(stepz)*DSMCSampVolWe%BGMdeltas(3)
  IF ( (GuessPos(1).LT.GEO%xmin).OR.(GuessPos(1).GT.GEO%xmax).OR. &
      (GuessPos(2).LT.GEO%ymin).OR.(GuessPos(2).GT.GEO%ymax).OR. &
      (GuessPos(3).LT.GEO%zmin).OR.(GuessPos(3).GT.GEO%zmax)) THEN
    CYCLE
  END IF
   !--- get background mesh cell of particle
   CellX = CEILING((GuessPos(1)-GEO%xminglob)/GEO%FIBGMdeltas(1)) 
   CellX = MIN(GEO%FIBGMimax,CellX)                             
   CellY = CEILING((GuessPos(2)-GEO%yminglob)/GEO%FIBGMdeltas(2))
   CellY = MIN(GEO%FIBGMjmax,CellY) 
   CellZ = CEILING((GuessPos(3)-GEO%zminglob)/GEO%FIBGMdeltas(3))
   CellZ = MIN(GEO%FIBGMkmax,CellZ)
   !--- check all cells associated with this beckground mesh cell
   DO iElem = 1, GEO%FIBGM(CellX,CellY,CellZ)%nElem
      Element = GEO%FIBGM(CellX,CellY,CellZ)%Element(iElem)
      CALL Eval_xyz_ElemCheck(GuessPos,Xi,Element)
      IF(MAXVAL(ABS(Xi)).GT.epsOneCell)THEN ! particle outside
        alpha1 = (GuessPos(1) / DSMCSampVolWe%BGMdeltas(1)) - i
        alpha2 = (GuessPos(2) / DSMCSampVolWe%BGMdeltas(2)) - j
        alpha3 = (GuessPos(3) / DSMCSampVolWe%BGMdeltas(3)) - k
         Found = Found + (1.-ABS(alpha1))*(1.-ABS(alpha2))*(1.-ABS(alpha3)) &
             *DSMCSampVolWe%w_VolInt(stepx)*DSMCSampVolWe%w_VolInt(stepy)*DSMCSampVolWe%w_VolInt(stepz)
        EXIT
      END IF
   END DO
END DO; END DO; END DO;
Volume = REAL(Found)*DSMCSampVolWe%BGMVolume

END SUBROUTINE VolumeBoundBGMCInt


SUBROUTINE WriteDSMCToHDF5(MeshFileName,OutputTime)
!===================================================================================================================================
! Writes DSMC state values to HDF5
!===================================================================================================================================
! MODULES
   USE MOD_Globals
   USE MOD_PreProc
   USE MOD_io_HDF5
   USE MOD_HDF5_output,   ONLY:WriteArrayToHDF5,WriteAttributeToHDF5,WriteHDF5Header
   USE MOD_PARTICLE_Vars, ONLY:nSpecies
   USE MOD_Mesh_Vars,     ONLY:offsetElem,nGlobalElems
   USE MOD_DSMC_Vars,     ONLY:MacroDSMC, CollisMode, DSMC, useDSMC
   USE MOD_Output_Vars,   ONLY:ProjectName
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
  REAL,INTENT(IN)                 :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  CHARACTER(LEN=255)                 :: FileName,FileString,Statedummy
  INTEGER                             :: nVal
!===================================================================================================================================
  SWRITE(*,*) ' WRITE DSMCSTATE TO HDF5 FILE...'
  FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCState',OutputTime)
  FileString=TRIM(FileName)//'.h5'
#ifdef MPI
  CALL OpenDataFile(FileString,create=.TRUE.,single=.FALSE.)
#else
  CALL OpenDataFile(FileString,create=.TRUE.)
#endif
  Statedummy = 'DSMCState'
  CALL WriteHDF5Header(Statedummy,File_ID)

  nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)

  CALL WriteArrayToHDF5(DataSetName='DSMC_velx', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV(1))

  CALL WriteArrayToHDF5(DataSetName='DSMC_vely', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV(2))

  CALL WriteArrayToHDF5(DataSetName='DSMC_velz', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV(3))

  CALL WriteArrayToHDF5(DataSetName='DSMC_vel', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV(4))

  CALL WriteArrayToHDF5(DataSetName='DSMC_velx2', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV2(1))

  CALL WriteArrayToHDF5(DataSetName='DSMC_vely2', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV2(2))

  CALL WriteArrayToHDF5(DataSetName='DSMC_velz2', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartV2(3))

  CALL WriteArrayToHDF5(DataSetName='DSMC_tempx', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Temp(1))

  CALL WriteArrayToHDF5(DataSetName='DSMC_tempy', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Temp(2))

  CALL WriteArrayToHDF5(DataSetName='DSMC_tempz', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Temp(3))

  CALL WriteArrayToHDF5(DataSetName='DSMC_temp', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Temp(4))

  CALL WriteArrayToHDF5(DataSetName='DSMC_dens', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%NumDens)

  CALL WriteArrayToHDF5(DataSetName='DSMC_partnum', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%PartNum)

  IF (useDSMC) CALL WriteArrayToHDF5(DataSetName='DSMC_collprob', rank=2,&
                        nValGlobal=(/nGlobalElems, nSpecies+1/),&
                        nVal=      (/PP_nElems,    nSpecies+1/),&
                        offset=    (/offsetElem, 0  /),&
                        collective=.TRUE., existing=.FALSE., RealArray=DSMC%CollProbOut(:,:))

  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
    CALL WriteArrayToHDF5(DataSetName='DSMC_tvib', rank=2,&
                          nValGlobal=(/nGlobalElems, nSpecies+1/),&
                          nVal=      (/PP_nElems,    nSpecies+1/),&
                          offset=    (/offsetElem, 0  /),&
                          collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Tvib)

    CALL WriteArrayToHDF5(DataSetName='DSMC_trot', rank=2,&
                          nValGlobal=(/nGlobalElems, nSpecies+1/),&
                          nVal=      (/PP_nElems,    nSpecies+1/),&
                          offset=    (/offsetElem, 0  /),&
                          collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Trot)
  END IF

  IF (DSMC%ElectronicState) THEN
    CALL WriteArrayToHDF5(DataSetName='DSMC_telec', rank=2,&
                          nValGlobal=(/nGlobalElems, nSpecies+1/),&
                          nVal=      (/PP_nElems,    nSpecies+1/),&
                          offset=    (/offsetElem, 0  /),&
                          collective=.TRUE., existing=.FALSE., RealArray=MacroDSMC(:,:)%Telec)
  END IF

  CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
  CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)

  CALL CloseDataFile()

END SUBROUTINE WriteDSMCToHDF5

END MODULE MOD_DSMC_Analyze
