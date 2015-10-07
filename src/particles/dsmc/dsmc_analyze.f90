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

INTERFACE WriteOutputMesh
  MODULE PROCEDURE WriteOutputMesh
END INTERFACE

INTERFACE WriteDSMCToHDF5
  MODULE PROCEDURE WriteDSMCToHDF5
END INTERFACE

INTERFACE WriteOutputMeshSamp
  MODULE PROCEDURE WriteOutputMeshSamp
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

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_data_sampling, WriteOutputMesh, WriteDSMCToHDF5, WriteOutputMeshSamp
PUBLIC :: DSMC_output_calc, CalcTVib, CalcSurfaceValues, CalcTelec, CalcTVibPoly
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

SUBROUTINE WriteOutputMesh()
!===================================================================================================================================
! Subroutine writes the mesh for debugging  
!===================================================================================================================================
! MODULES
   USE MOD_Particle_Vars
   USE MOD_DSMC_Vars,      ONLY : CollisMode, useDSMC  
   USE MOD_Mesh_Vars,     ONLY : nElems!, nNodes
#ifdef MPI
   USE MOD_Particle_MPI_Vars, ONLY: PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                                
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  CHARACTER(LEN=26)                  :: myFileName
  INTEGER                            :: iElem, iNode, iSpec, iInit
  INTEGER                            :: withMolecules
!===================================================================================================================================

  IF (useDSMC) THEN
    IF (CollisMode.EQ.1) THEN
      withMolecules = 0  
    ELSE
      withMolecules = 1
    END IF
  ELSE
    withMolecules = 0
  END IF
  STOP

!#ifdef MPI
!  WRITE(myFileName,'(A8,I5.5,A4)')'DSMCMesh',PartMPI%MyRank,'.vtk'
!#else
!  WRITE(myFileName,'(A12)')'DSMCMesh.vtk'
!#endif  
!  OPEN(1112,FILE=myFileName,STATUS='replace')
!  WRITE(1112,'(A36,I0,A11,I0)')'# vtk DataFile Version 2.0 Species: ',nSpecies+1,' Molcules: ',withMolecules 
!  WRITE(1112,'(A)')'Debug Mesh '
!  WRITE(1112,'(A)')'ASCII'
!  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
!  WRITE(1112,'(A)')''
!  WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
!  DO iNode=1, nNodes
!    WRITE(1112,*) GEO%NodeCoords(1:3, iNode)
!  END DO
!  WRITE(1112,*)''
!  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nElems,9*nElems
!  DO iElem=1, nElems
!    WRITE(1112,'(I0)',ADVANCE="NO")8
!    DO iNode=1, 8
!    WRITE(1112,'(1X,I0)',ADVANCE="NO") GEO%ElemToNodeID(iNode,iElem) -1
!    END DO
!    WRITE(1112,*)''
!  END DO
!  WRITE(1112,*)''
!  WRITE(1112,'(A,I0)')'CELL_TYPES ',nElems
!  DO iElem=1,nElems
!    WRITE(1112,'(1X,I0)',ADVANCE="NO")12
!  END DO  
!  WRITE(1112,*)''
!  WRITE(1112,*)''
!  WRITE(1112,'(A,I0)')'CELL_DATA ',nElems
!  DO iSpec=1, nSpecies
!    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
!      IF ((Species(iSpec)%Init(iInit)%ParticleEmissionType.GE.3).AND.(Species(iSpec)%Init(iInit)%ParticleEmissionType.LE.6)) THEN
!        WRITE(1112,'(A32,I3.3,A5,I3.3,A)')'SCALARS PressureElemType_Species', iSpec, '_Init', iInit, ' FLOAT'
!        WRITE(1112,'(A)')'LOOKUP_TABLE default'
!        DO iElem = 1, nElems
!          WRITE(1112,*) Species(iSpec)%Init(iInit)%ConstPress%ElemStat(iElem)
!        END DO
!      END IF
!    END DO
!  END DO
!  CLOSE(1112)
END SUBROUTINE WriteOutputMesh


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
  SWRITE(*,*) 'aetsch'
  RETURN


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
  USE MOD_DSMC_Vars,       ONLY:SurfMesh, SampWall, MacroSurfaceVal, DSMC, realtime
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
  INTEGER                            :: iElem, iSpec
  REAL                               :: TimeSample
  INTEGER, ALLOCATABLE               :: CounterTotal(:), SumCounterTotal(:)              ! Total Wall-Collision counter
!===================================================================================================================================
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
!  USE MOD_part_MPI_Vars, ONLY : PMPIVAR, MPIGEO
  USE mpi
  USE MOD_Globals,           ONLY:IERROR, MPISTATUS
  USE MOD_DSMC_Vars,         ONLY:SurfMesh, SampWall, SampWallHaloCell
  USE MOD_Mesh_Vars,         ONLY:ElemToSide
  USE MOD_Particle_Vars,     ONLY:nSpecies
  USE MOD_Particle_MPI_Vars, ONLY:PartMPI
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  TYPE tMPISurfContent
    REAL, POINTER                    :: content(:) =>NULL()
  END TYPE
  INTEGER                            :: iSide, iProc, Element, iLocSide, iCount, iSampWallAlloc
  INTEGER                            :: SurfSideID    
  INTEGER, ALLOCATABLE               :: RecvMsgSurfs(:), SendMsgSurfs(:), PosCount(:)
  TYPE(tMPISurfContent), POINTER     :: SendContent(:) => NULL()
  TYPE(tMPISurfContent), POINTER     :: RecvContent(:) => NULL()
!===================================================================================================================================

  ALLOCATE(SendMsgSurfs(0:PartMPI%nProcs-1))
  ALLOCATE(RecvMsgSurfs(0:PartMPI%nProcs-1))
  ALLOCATE(SendContent (0:PartMPI%nProcs-1))
  ALLOCATE(RecvContent (0:PartMPI%nProcs-1))
  ALLOCATE(PosCount    (0:PartMPI%nProcs-1))
  PosCount    (0:PartMPI%nProcs-1) = 0
  SendMsgSurfs(0:PartMPI%nProcs-1) = 0
  RecvMsgSurfs(0:PartMPI%nProcs-1) = 0

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
!    ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
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
!  DO iProc=0, PMPIVAR%nProcs-1          
!    IF (PMPIVAR%iProc.LT.iProc) THEN
!      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
!          (14+nSpecies)*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,IERROR)    
!      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, &
!          (14+nSpecies)*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR) 
!    ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
!      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, & 
!          (14+nSpecies)*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)  
!      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
!          (14+nSpecies)*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,IERROR)    
!    END IF
!  END DO
!
!!---- sum up surface values  
!  DO iProc=0, PMPIVAR%nProcs-1
!!    IF (RecvMsgSurfs(iProc).ne.0) print*,RecvMsgSurfs(iProc)
!    IF ((iProc.NE.PMPIVAR%iProc).AND.(RecvMsgSurfs(iProc).NE.0)) THEN
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
!  DO iProc=0, PMPIVAR%nProcs-1
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


SUBROUTINE WriteOutputMeshSamp()
!===================================================================================================================================
! Subroutine to write the mesh with sampling values of emistype6-BC as data
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars
  USE MOD_Mesh_Vars,     ONLY : nElems!, nNodes
  USE MOD_Globals
#ifdef MPI
  !USE MOD_part_MPI_Vars, ONLY : PMPIVAR
  USE MOD_Particle_MPI_Vars, ONLY: PartMPI
#endif
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                          :: iInit                                                          
  CHARACTER(LEN=255)              :: myFileName
  INTEGER                          :: iElem, iNode, iSpec, Elem
  REAL, ALLOCATABLE               :: ElemSampOutput(:,:)
!===================================================================================================================================

  STOP 
  ALLOCATE (ElemSampOutput(nElems,0:6))
#ifdef MPI
  WRITE(myFileName,'(A13,I5.5)')'DSMCMesh_Samp',PartMPI%MyRank
  myFileName=TRIM(TIMESTAMP(myFileName,Time))//'.vtk'
#else
  myFileName=TRIM(TIMESTAMP('DSMCMesh_Samp',Time))//'.vtk'
#endif
  OPEN(1503,FILE=myFileName,STATUS='replace')
  WRITE(1503,'(A26)')'# vtk DataFile Version 2.0'
  WRITE(1503,'(A)')'Debug Mesh2 '
  WRITE(1503,'(A)')'ASCII'
  WRITE(1503,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1503,'(A)')''
! WRITE(1503,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
!  DO iNode=1, nNodes
!    WRITE(1503,*) GEO%NodeCoords(1:3, iNode)
!  END DO
  WRITE(1503,*)''
  WRITE(1503,'(A,I0,1X,I0)')'CELLS ',nElems,9*nElems
  DO iElem=1, nElems
    WRITE(1503,'(I0)',ADVANCE="NO")8
    !DO iNode=1, 8
    !  !WRITE(1503,'(1X,I0)',ADVANCE="NO") GEO%ElemToNodeID(iNode,iElem) -1
    !END DO
    WRITE(1503,*)''
  END DO
  WRITE(1503,*)''
  WRITE(1503,'(A,I0)')'CELL_TYPES ',nElems
  DO iElem=1,nElems
    WRITE(1503,'(1X,I0)',ADVANCE="NO")12
  END DO
  WRITE(1503,*)''
  WRITE(1503,*)''
  WRITE(1503,'(A,I0)')'CELL_DATA ',nElems
  DO iSpec=1, nSpecies
    !write Samplingvalues in array of all elements
    ElemSampOutput(:,:)=0.
    DO iInit = Species(iSpec)%StartnumberOfInits, Species(iSpec)%NumberOfInits
      IF (Species(iSpec)%Init(iInit)%ParticleEmissionType.EQ.6) THEN
        DO iElem = 1,Species(iSpec)%Init(iInit)%ConstPress%nElemTotalInside
          Elem = Species(iSpec)%Init(iInit)%ConstPress%ElemTotalInside(iElem)
          IF (ElemSampOutput(Elem,0).EQ.0.) THEN
            ElemSampOutput(Elem,0)=1.
            ElemSampOutput(Elem,1:6)=Species(iSpec)%Init(iInit)%ConstPress%ConstPressureSamp(iElem,1:6)
          END IF
        END DO
      END IF
    END DO
    !write data
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_status_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,0)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_vx_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,1)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_vy_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,2)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_vz_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,3)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_n_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,4)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_p_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,5)
    END DO
    WRITE(1503,'(A,I3.3,A)')'SCALARS PressureElem_a2_', iSpec, ' FLOAT'
    WRITE(1503,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      WRITE(1503,*) ElemSampOutput(iElem,6)
    END DO
  END DO
  CLOSE(1503)
  DEALLOCATE(ElemSampOutput)
END SUBROUTINE WriteOutputMeshSamp

END MODULE MOD_DSMC_Analyze
