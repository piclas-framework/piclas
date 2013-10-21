MODULE MOD_DSMC_Analyze
!==================================================================================================
! module for DSMC Sampling and Output
!==================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!--------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!--------------------------------------------------------------------------------------------------
! Private Part ------------------------------------------------------------------------------------
! Public Part -------------------------------------------------------------------------------------
PUBLIC :: DSMC_data_sampling, WriteOutputMesh, WriteDSMCToHDF5
PUBLIC :: DSMC_output_calc, CalcTVib, CalcSurfaceValues, OutputMaxCollProb, CalcTelec
!===================================================================================================

CONTAINS

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_data_sampling()

  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, PartStateIntEn, CollisMode, SpecDSMC
  USE MOD_Particle_Vars,          ONLY : PEM, PDM, PartSpecies, PartState, PartMPF, usevMPF
!--------------------------------------------------------------------------------------------------!
! statistical pairing method                                                                       !
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: iPart, iElem
!--------------------------------------------------------------------------------------------------!
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
        ! if usevMPF SampDSMC(iElem,PartSpecies(iPart))%PartNum == real number of particles
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1) * PartMPF(iPart)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2) * PartMPF(iPart)
        END IF
        IF (DSMC%ElectronicState) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
                                                       + PartStateIntEn(iPart,3) * PartMPF(iPart)
        END IF
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartState(iPart,4:6)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1:3) &
                                                       + PartState(iPart,4:6)**2
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2)
        END IF
        IF (DSMC%ElectronicState) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EElec     = SampDSMC(iElem,PartSpecies(iPart))%EElec &
                                                       + PartStateIntEn(iPart,3)
        END IF
      END IF
    END IF
  END DO

END SUBROUTINE DSMC_data_sampling

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_output_calc(nOutput)

  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, MacroDSMC, CollisMode, SpecDSMC, realtime
  USE MOD_Mesh_Vars,              ONLY : nElems,MeshFile
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst, Species, GEO, PartSpecies, usevMPF
!--------------------------------------------------------------------------------------------------!
! statistical pairing method                                                                       !
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: iPart, iElem, iSpec, MolecPartNum, HeavyPartNum                   !
REAL                          :: TVib_TempFac
! input variable declaration                                                                       !
INTEGER, INTENT(IN)           :: nOutput                                                           !
!--------------------------------------------------------------------------------------------------!

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
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
          IF (DSMC%VibEnergyModel.EQ.0) THEN              ! SHO-model
              TVib_TempFac=SampDSMC(iElem,iSpec)%EVib &
                      /(SampDSMC(iElem,iSpec)%PartNum*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
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
!              MacroDSMC(iElem,iSpec)%TVib = SpecDSMC(iSpec)%CharaTVib / LOG(1.0 + &
!                                SpecDSMC(iSpec)%CharaTVib * BoltzmannConst * SampDSMC(iElem, iSpec)%PartNum &
!                                 /SampDSMC(iElem,iSpec)%EVib)
!              print*,'TVIB:',MacroDSMC(iElem,iSpec)%TVib
        IF (DSMC%ElectronicState) THEN
          MacroDSMC(iElem,iSpec)%TElec= CalcTelec( SampDSMC(iElem,iSpec)%EElec/SampDSMC(iElem,iSpec)%PartNum, iSpec)
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
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
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
      IF ( DSMC%ElectronicState .AND. (SpecDSMC(iSpec)%InterID.NE.4) ) THEN
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
  END DO
  CALL WriteDSMCToHDF5(TRIM(MeshFile),realtime)
  CALL WriteOutputDSMC(nOutput)
  DEALLOCATE(MacroDSMC)
  
END SUBROUTINE DSMC_output_calc

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE WriteOutputMesh()
   USE MOD_Particle_Vars,  ONLY : GEO,nSpecies
   USE MOD_DSMC_Vars,      ONLY : CollisMode, useDSMC  
   USE MOD_Mesh_Vars,     ONLY : nElems, nNodes
#ifdef MPI
   USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                      !                                       
! Local variable declaration                                                                       !
  CHARACTER(LEN=26)                  :: myFileName
  INTEGER                            :: iElem, iNode
  INTEGER                            :: withMolecules                   

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  IF (useDSMC) THEN
    IF (CollisMode.EQ.1) THEN
      withMolecules = 0  
    ELSE
      withMolecules = 1
    END IF
  ELSE
    withMolecules = 0
  END IF

#ifdef MPI
  WRITE(myFileName,'(A8,I5.5,A4)')'DSMCMesh',PMPIVAR%iProc,'.vtk'
#else
  WRITE(myFileName,'(A12)')'DSMCMesh.vtk'
#endif  
  OPEN(1112,FILE=myFileName,STATUS='replace')
  WRITE(1112,'(A36,I0,A11,I0)')'# vtk DataFile Version 2.0 Species: ',nSpecies+1,' Molcules: ',withMolecules 
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
  DO iNode=1, nNodes
    WRITE(1112,*) GEO%NodeCoords(1:3, iNode)
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nElems,9*nElems
  DO iElem=1, nElems
    WRITE(1112,'(I0)',ADVANCE="NO")8
    DO iNode=1, 8
    WRITE(1112,'(1X,I0)',ADVANCE="NO") GEO%ElemToNodeID(iNode,iElem) -1
    END DO
    WRITE(1112,*)''
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',nElems
  DO iElem=1,nElems
    WRITE(1112,'(1X,I0)',ADVANCE="NO")12
  END DO  
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_DATA ',nElems
  CLOSE(1112)

END SUBROUTINE WriteOutputMesh

!--------------------------------------------------------------------------------------------------!

SUBROUTINE WriteOutputDSMC(nOutput)
   USE MOD_Globals
   USE MOD_Particle_Vars
   USE MOD_Mesh_Vars,     ONLY : nElems, nNodes
   USE MOD_DSMC_Vars,     ONLY : MacroDSMC, CollisMode, SpecDSMC, useDSMC
#ifdef MPI
   USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  CHARACTER(LEN=26)                  :: myFileName
  CHARACTER(LEN=43)                  :: copyCommand 
  INTEGER                            :: iSpec, iElem, iNode
  INTEGER, INTENT(IN)                :: nOutput                                         
  INTEGER                            :: withMolecules   
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  IF(MPIroot)THEN
    WRITE(UNIT_StdOut,'(132("-"))')      
    WRITE(UNIT_stdOut,'(A)') 'WRITE DSMC OUTPUT'
  END IF

  IF (useDSMC) THEN
    IF (CollisMode.EQ.1) THEN
      withMolecules = 0  
    ELSE
      withMolecules = 1
    END IF
  ELSE
    withMolecules = 0
  END IF

#ifdef MPI !MPI 
  WRITE(copyCommand,'(A8,I5.5,A1,I4.4,A4)')'DSMCOut_',PMPIVAR%iProc,'_',nOutput,'.vtk'
#else
! copy of the Outputmesh to add the macroscopic values
  write(copyCommand,'(A8,I4.4,A4)') 'DSMCOut_',nOutput,'.vtk'
#endif
  OPEN (1112, file=copyCommand, STATUS='replace')
  WRITE(1112,'(A36,I0,A11,I0)')'# vtk DataFile Version 2.0 Species: ',nSpecies+1,' Molcules: ',withMolecules 
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',nNodes,' FLOAT'
  DO iNode=1, nNodes
    WRITE(1112,*) GEO%NodeCoords(1:3, iNode)
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',nElems,9*nElems
  DO iElem=1, nElems
    WRITE(1112,'(I0)',ADVANCE="NO")8
    DO iNode=1, 8
    WRITE(1112,'(1X,I0)',ADVANCE="NO") GEO%ElemToNodeID(iNode,iElem) -1
    END DO
    WRITE(1112,*)''
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',nElems
  DO iElem=1,nElems
    WRITE(1112,'(1X,I0)',ADVANCE="NO")12
  END DO  
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_DATA ',nElems
  DO iSpec=1, nSpecies + 1
    WRITE(1112,'(A,I3.3,A)')'SCALARS VelocityX_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%PartV(1)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS VelocityY_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%PartV(2)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS VelocityZ_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%PartV(3)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS VelocityTotal_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%PartV(4)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS TempX_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%Temp(1)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS TempY_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%Temp(2)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS TempZ_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%Temp(3)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS TempTotal_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%Temp(4)
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS PartNum_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%PartNum
      END IF
    END DO
    WRITE(1112,'(A,I3.3,A)')'SCALARS Dens_Species_',iSpec,' FLOAT'
    WRITE(1112,'(A)')'LOOKUP_TABLE default'
    DO iElem = 1, nElems
      IF (MacroDSMC(iElem,iSpec)%PartNum .eq.0) THEN
        WRITE(1112,*) '0'
      ELSE
        WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%NumDens
      END IF
    END DO
    IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
      WRITE(1112,'(A,I3.3,A)')'SCALARS TVib_Species_',iSpec,' FLOAT'
      WRITE(1112,'(A)')'LOOKUP_TABLE default'
      DO iElem = 1, nElems
        IF (MacroDSMC(iElem,iSpec)%PartNum.eq.0) THEN
          WRITE(1112,*) '0'
        ELSE
          WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%TVib
        END IF
      END DO
      WRITE(1112,'(A,I3.3,A)')'SCALARS TRot_Species_',iSpec,' FLOAT'
      WRITE(1112,'(A)')'LOOKUP_TABLE default'
      DO iElem = 1, nElems
        IF (MacroDSMC(iElem,iSpec)%PartNum.eq.0) THEN
          WRITE(1112,*) '0'
        ELSE
          WRITE(1112,'(E22.15)') MacroDSMC(iElem,iSpec)%TRot
        END IF
      END DO
    END IF
   END DO  
      
  IF(MPIroot)THEN
    WRITE(UNIT_stdOut,'(A)') 'DSMC OUTPUT......DONE'
    WRITE(UNIT_StdOut,'(132("-"))')
  END IF

END SUBROUTINE WriteOutputDSMC

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE WriteDSMCToHDF5(MeshFileName,OutputTime)
   USE MOD_Globals
   USE MOD_PreProc
   USE MOD_io_HDF5
   USE MOD_HDF5_output
   USE MOD_PARTICLE_Vars, ONLY:nSpecies
   USE MOD_Mesh_Vars,ONLY:offsetElem,nGlobalElems
   USE MOD_DSMC_Vars,     ONLY :MacroDSMC, CollisMode, DSMC
   USE MOD_Output_Vars,   ONLY:ProjectName
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  CHARACTER(LEN=255)                 :: FileName,FileString,MeshFile255,Statedummy
  INTEGER                        :: nVal
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
WRITE(*,*) ' WRITE DSMCSTATE TO HDF5 FILE...'
FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCState',OutputTime)
FileString=TRIM(FileName)//'.h5'
CALL OpenDataFile(Filestring,.TRUE.)

Statedummy = 'DSMCState'
CALL WriteHDF5Header(Statedummy,File_ID)
! Write DG solution ----------------------------------------------------------------------------------------------------------------
nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
CALL WriteArrayToHDF5('DSMC_velx',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(1))
CALL WriteArrayToHDF5('DSMC_vely',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(2))
CALL WriteArrayToHDF5('DSMC_velz',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(3))
CALL WriteArrayToHDF5('DSMC_vel',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(4))
CALL WriteArrayToHDF5('DSMC_velx2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(1))
CALL WriteArrayToHDF5('DSMC_vely2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(2))
CALL WriteArrayToHDF5('DSMC_velz2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(3))
CALL WriteArrayToHDF5('DSMC_tempx',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(1))
CALL WriteArrayToHDF5('DSMC_tempy',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(2))
CALL WriteArrayToHDF5('DSMC_tempz',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(3))
CALL WriteArrayToHDF5('DSMC_temp',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(4))
CALL WriteArrayToHDF5('DSMC_dens',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%NumDens)
CALL WriteArrayToHDF5('DSMC_partnum',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartNum)
IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
  CALL WriteArrayToHDF5('DSMC_tvib',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%TVib)  
  CALL WriteArrayToHDF5('DSMC_trot',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%TRot)  
END IF
IF (DSMC%ElectronicState) THEN
  CALL WriteArrayToHDF5('DSMC_telec',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Telec)  
END IF
CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
MeshFile255=TRIM(MeshFileName)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=MeshFile255)
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)

CALL CloseDataFile()

END SUBROUTINE WriteDSMCToHDF5

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE WriteDSMCSurfToHDF5(MeshFileName,OutputTime)
   USE MOD_Globals
   USE MOD_PreProc
   USE MOD_io_HDF5
   USE MOD_HDF5_output
   USE MOD_PARTICLE_Vars, ONLY:nSpecies
   USE MOD_Mesh_Vars,ONLY:offsetElem,nGlobalElems
   USE MOD_DSMC_Vars,     ONLY :SurfMesh, MacroSurfaceVal , CollisMode, DSMC
   USE MOD_Output_Vars,   ONLY:ProjectName
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  CHARACTER(LEN=255)                 :: FileName,FileString,MeshFile255,Statedummy
  INTEGER                        :: nVal
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
WRITE(*,*) ' WRITE DSMCSurfSTATE TO HDF5 FILE...'
FileName=TIMESTAMP(TRIM(ProjectName)//'_DSMCSurfState',OutputTime)
FileString=TRIM(FileName)//'.h5'
CALL OpenDataFile(Filestring,.TRUE.)

Statedummy = 'DSMCSurfState'
CALL WriteHDF5Header(Statedummy,File_ID)
! Write DG solution ----------------------------------------------------------------------------------------------------------------
nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
!CALL WriteArrayToHDF5('DSMC_velx',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(1))
!CALL WriteArrayToHDF5('DSMC_vely',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(2))
!CALL WriteArrayToHDF5('DSMC_velz',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(3))
!CALL WriteArrayToHDF5('DSMC_vel',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV(4))
!CALL WriteArrayToHDF5('DSMC_velx2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(1))
!CALL WriteArrayToHDF5('DSMC_vely2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(2))
!CALL WriteArrayToHDF5('DSMC_velz2',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartV2(3))
!CALL WriteArrayToHDF5('DSMC_tempx',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(1))
!CALL WriteArrayToHDF5('DSMC_tempy',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(2))
!CALL WriteArrayToHDF5('DSMC_tempz',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(3))
!CALL WriteArrayToHDF5('DSMC_temp',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%Temp(4))
!CALL WriteArrayToHDF5('DSMC_dens',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%NumDens)
!CALL WriteArrayToHDF5('DSMC_partnum',nGlobalElems,2,(/PP_nElems, nSpecies+1/),offsetElem,1,RealArray=MacroDSMC(:,:)%PartNum)

CALL WriteAttributeToHDF5(File_ID,'DSMC_nSpecies',1,IntegerScalar=nSpecies)
CALL WriteAttributeToHDF5(File_ID,'DSMC_CollisMode',1,IntegerScalar=CollisMode)
MeshFile255=TRIM(MeshFileName)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=MeshFile255)
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)

!CALL WriteAttributeToHDF5(File_ID,'DSMCSurf_nSurfNode',1,IntegerScalar=SurfMesh%nSurfaceNode)
!CALL WriteArrayToHDF5('DSMC_partnum',nGlobalElems,2,(/PP_nElems, SurfMesh%nSurfaceNode/), &
!      offsetElem,1,IntegerArray=SurfMesh%BCSurfNodes(:)

CALL CloseDataFile()

END SUBROUTINE WriteDSMCSurfToHDF5

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!



SUBROUTINE CalcSurfaceValues(nOutput)
   USE MOD_DSMC_Vars,      ONLY : SurfMesh, SampWall, MacroSurfaceVal, DSMC
   USE MOD_Particle_Vars,  ONLY : Time
   USE MOD_TimeDisc_Vars,  ONLY : TEnd  
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                      !      
  INTEGER, INTENT(INOUT)             :: nOutput                                 
! Local variable declaration                                                                       !
  INTEGER                            :: iElem             
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
#ifdef MPI
  CALL MPISurfaceValuesSend()  
#endif

  ALLOCATE(MacroSurfaceVal(SurfMesh%nSurfaceBCSides))
  DO iElem=1,SurfMesh%nSurfaceBCSides
    MacroSurfaceVal(iElem)%Heatflux = (SampWall(iElem)%Energy(1)+SampWall(iElem)%Energy(4)+SampWall(iElem)%Energy(7) &
                                      -SampWall(iElem)%Energy(3)-SampWall(iElem)%Energy(6)-SampWall(iElem)%Energy(9))&
                                     /(SurfMesh%SurfaceArea(iElem) * (Time-(1-DSMC%TimeFracSamp)*TEnd))
    MacroSurfaceVal(iElem)%Force(1) = SampWall(iElem)%Force(1) /(SurfMesh%SurfaceArea(iElem) * (Time-(1-DSMC%TimeFracSamp)*TEnd))
    MacroSurfaceVal(iElem)%Force(2) = SampWall(iElem)%Force(2) /(SurfMesh%SurfaceArea(iElem) * (Time-(1-DSMC%TimeFracSamp)*TEnd))
    MacroSurfaceVal(iElem)%Force(3) = SampWall(iElem)%Force(3) / (SurfMesh%SurfaceArea(iElem) * (Time-(1-DSMC%TimeFracSamp)*TEnd))
    MacroSurfaceVal(iElem)%Counter(1) = SampWall(iElem)%Counter(1) / (Time-(1-DSMC%TimeFracSamp)*TEnd)
  END DO

  CALL WriteOutputSurfaceMesh(nOutput)
  DEALLOCATE(MacroSurfaceVal)

END SUBROUTINE CalcSurfaceValues

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE WriteOutputSurfaceMesh(nOutput)
   USE MOD_Globals
   USE MOD_Particle_Vars
   USE MOD_DSMC_Vars,      ONLY : SurfMesh, MacroSurfaceVal  
#ifdef MPI
   USE MOD_part_MPI_Vars, ONLY : PMPIVAR
#endif
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                      !                                       
  INTEGER, INTENT(IN)                :: nOutput
! Local variable declaration                                                                       !
  CHARACTER(LEN=255)                  :: myFileName
  INTEGER                            :: iElem, iNode             
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

  IF(MPIroot)THEN
    WRITE(UNIT_StdOut,'(132("-"))')      
    WRITE(UNIT_stdOut,'(A)') 'WRITE SURFACE OUTPUT'
  END IF


#ifdef MPI
  WRITE(myFileName,'(A17,I5.5,A1,I4.4,A4)')'DSMCSurfaceValues',PMPIVAR%iProc,'_',nOutput,'.vtk'
#else
  WRITE(myFileName,'(A18,I4.4,A4)')'DSMCSurfaceValues_',nOutput,'.vtk'
#endif  
  OPEN(1112,FILE=myFileName,STATUS='replace')
  WRITE(1112,'(A)')'# vtk DataFile Version 2.0'
  WRITE(1112,'(A)')'Debug Mesh '
  WRITE(1112,'(A)')'ASCII'
  WRITE(1112,'(A)')'DATASET UNSTRUCTURED_GRID'
  WRITE(1112,'(A)')''
  WRITE(1112,'(A,I0,A)')'POINTS ',SurfMesh%nSurfaceNode,' FLOAT'
  DO iNode=1, SurfMesh%nSurfaceNode
    WRITE(1112,*) GEO%NodeCoords(1:3, SurfMesh%BCSurfNodes(iNode))
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0,1X,I0)')'CELLS ',SurfMesh%nSurfaceBCSides,5*SurfMesh%nSurfaceBCSides
  DO iElem=1, SurfMesh%nSurfaceBCSides
    WRITE(1112,'(I0)',ADVANCE="NO")4
    DO iNode=1, 4
    WRITE(1112,'(1X,I0)',ADVANCE="NO") SurfMesh%SideSurfNodeMap(iNode,iElem) -1
    END DO
    WRITE(1112,*)''
  END DO
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_TYPES ',SurfMesh%nSurfaceBCSides
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(1X,I0)',ADVANCE="NO")9
  END DO  
  WRITE(1112,*)''
  WRITE(1112,*)''
  WRITE(1112,'(A,I0)')'CELL_DATA ',SurfMesh%nSurfaceBCSides
  WRITE(1112,'(A)')'SCALARS Heatflux FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(E22.15)') MacroSurfaceVal(iElem)%Heatflux
  END DO
  WRITE(1112,'(A)')'SCALARS ForceX FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(E22.15)') MacroSurfaceVal(iElem)%Force(1)
  END DO
  WRITE(1112,'(A)')'SCALARS ForceY FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(E22.15)') MacroSurfaceVal(iElem)%Force(2)
  END DO
  WRITE(1112,'(A)')'SCALARS ForceZ FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(E22.15)') MacroSurfaceVal(iElem)%Force(3)
  END DO
  WRITE(1112,'(A)')'SCALARS ImpactCount FLOAT'
  WRITE(1112,'(A)')'LOOKUP_TABLE default'
  DO iElem=1,SurfMesh%nSurfaceBCSides
    WRITE(1112,'(E22.15)') MacroSurfaceVal(iElem)%Counter(1)
  END DO
  CLOSE(1112)

  IF(MPIroot)THEN    
    WRITE(UNIT_stdOut,'(A)') 'SURFACE OUTPUT DONE'
    WRITE(UNIT_StdOut,'(132("-"))')  
  END IF

END SUBROUTINE WriteOutputSurfaceMesh

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE OutputMaxCollProb(Time)
   USE MOD_Particle_Vars,  ONLY : GEO
   USE MOD_DSMC_Vars,      ONLY : DSMC  
#ifdef MPI
   USE MOD_part_MPI_Vars,  ONLY : PMPIVAR
#endif
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
  REAL,INTENT(IN)                     :: Time
! Local variable declaration                                                                       !
  LOGICAL                             :: isOpen, FileExists                                        !
  CHARACTER(LEN=255)                  :: outfile
  INTEGER                             :: unit_index
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  unit_index = 5351
  DSMC%CollMean = DSMC%CollMean / DSMC%CollMeanCount
  INQUIRE(UNIT   = unit_index , OPENED = isOpen)
  IF (.NOT.isOpen) THEN
#ifdef MPI	
    WRITE(outfile,'(A11,I5.5,A4)')'CollProbMax',PMPIVAR%iProc,'.csv'
#else
    WRITE(outfile,'(A)')'CollProbMax.csv'
#endif
    INQUIRE(file=TRIM(outfile),EXIST=FileExists)
    OPEN(unit_index,file=TRIM(outfile))
         !--- insert header
    WRITE(unit_index,'(A14)',ADVANCE='NO') 'Time'
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(A14)',ADVANCE='NO') 'CollProbMax'    
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,'(A14)',ADVANCE='NO') 'CollProbMean'
    WRITE(unit_index,'(A1)') ' '
  END IF
  WRITE(unit_index,104,ADVANCE='NO') Time
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,104,ADVANCE='NO') DSMC%CollProbMax
  WRITE(unit_index,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index,104,ADVANCE='NO') DSMC%CollMean
  WRITE(unit_index,'(A1)') ' '

104    FORMAT (e25.14)

END SUBROUTINE OutputMaxCollProb

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
#ifdef MPI

SUBROUTINE MPISurfaceValuesSend()
  USE MOD_Globals,       ONLY : myRank, IERROR, MPISTATUS
  USE MOD_part_MPI_Vars, ONLY : PMPIVAR, MPIGEO
  USE MOD_DSMC_Vars,     ONLY : SurfMesh, SampWall, SampWallHaloCell
  USE MOD_Mesh_Vars,     ONLY : ElemToSide
  !----------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                  !
  !----------------------------------------------------------------------------------------------!
  INCLUDE 'mpif.h'                                                                               !
  !----------------------------------------------------------------------------------------------!
! argument list declaration                                      !      
                           
! Local variable declaration                                                                       !
  TYPE tMPISurfContent
    REAL, POINTER                    :: content(:) =>NULL()
  END TYPE
  INTEGER                            :: iSide, iProc, Element, iLocSide, iCount  
  INTEGER                            :: SurfSideID    
  INTEGER, ALLOCATABLE               :: RecvMsgSurfs(:), SendMsgSurfs(:), PosCount(:)
  TYPE(tMPISurfContent), POINTER     :: SendContent(:) => NULL()
  TYPE(tMPISurfContent), POINTER     :: RecvContent(:) => NULL()
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  ALLOCATE(SendMsgSurfs(0:PMPIVAR%nProcs-1))
  ALLOCATE(RecvMsgSurfs(0:PMPIVAR%nProcs-1))
  ALLOCATE(SendContent(0:PMPIVAR%nProcs-1))
  ALLOCATE(RecvContent(0:PMPIVAR%nProcs-1))
  ALLOCATE(PosCount(0:PMPIVAR%nProcs-1))
  PosCount(0:PMPIVAR%nProcs-1) = 0
  SendMsgSurfs(0:PMPIVAR%nProcs-1) = 0
  RecvMsgSurfs(0:PMPIVAR%nProcs-1) = 0
!  print*,'SendMsgSurfs(0)',SendMsgSurfs(0)
!  print*,'SendMsgSurfs(1)',SendMsgSurfs(1)
!  print*,'SendMsgSurfs(2)',SendMsgSurfs(2)
!---- calculation of send-massages number for each proc
  DO iSide = 1, SIZE(MPIGEO%BC,2)  !---- loop over all halosides that are defined as BC
    IF(SurfMesh%HaloSideIDToSurfSideMap(iSide).NE.0) THEN  !---- only surfaces (=wall-sides)   
       ! get halo cells ElemID. Not quite sure if ELEM_ID or NB_ELEM_ID. 
       ! maybe check this later and delete IF-query here to save cpu time
       IF (MPIGEO%SideToElem(1,iSide).NE.-1) THEN
         Element = MPIGEO%SideToElem(1,iSide)
       ELSE
         Element = MPIGEO%SideToElem(2,iSide)
       END IF
       iProc = MPIGEO%ElemMPIID(Element)
       SendMsgSurfs(iProc) = SendMsgSurfs(iProc) + 1
    END IF
  END DO
!---- comunicate send-massages number (=number of surfaces) for each proc
  DO iProc=0, PMPIVAR%nProcs-1 
    IF (PMPIVAR%iProc.LT.iProc) THEN
      CALL MPI_SEND(SendMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,IERROR)    
      CALL MPI_RECV(RecvMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)      
    ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
      CALL MPI_RECV(RecvMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)
      CALL MPI_SEND(SendMsgSurfs(iProc),1,MPI_INTEGER,iProc,1101,PMPIVAR%COMM,IERROR)
    END IF
  END DO
!---- allocate (send and recv) massages size (15 * send-massages number)
!---- 1=iLocSide; 2=target element; 3-11=Energy(1-9); 12-14=Force(1-3); 15=Counter(1)  per surface  
  DO iProc=0, PMPIVAR%nProcs-1
    IF (SendMsgSurfs(iProc).NE.0) THEN
      ALLOCATE(SendContent(iProc)%content(15*SendMsgSurfs(iProc)))
    END IF
    IF (RecvMsgSurfs(iProc).NE.0) THEN
      ALLOCATE(RecvContent(iProc)%content(15*RecvMsgSurfs(iProc)))
    END IF
  END DO
!---- build massage
  DO iSide = 1, SIZE(MPIGEO%BC,2)
    IF(SurfMesh%HaloSideIDToSurfSideMap(iSide).NE.0) THEN     
       ! get halo cells ElemID. Not quite sure if ELEM_ID or NB_ELEM_ID. 
       ! maybe check this later and delete IF-query here to save cpu time
       IF (MPIGEO%SideToElem(1,iSide).NE.-1) THEN
         Element = MPIGEO%SideToElem(1,iSide)
       ELSE
         Element = MPIGEO%SideToElem(2,iSide)
       END IF
       iProc = MPIGEO%ElemMPIID(Element)
       PosCount(iProc) = PosCount(iProc) + 1
       DO iLocSide=1,6  !---- search for iLocSide
         IF (MPIGEO%ElemToSide(1,iLocSide,Element).EQ.iSide) THEN
          SendContent(iProc)%content(PosCount(iProc))= REAL(iLocSide)  !---- 1=iLocSide
          EXIT
         END IF
       END DO
       PosCount(iProc) = PosCount(iProc) + 1
       SendContent(iProc)%content(PosCount(iProc))= REAL(MPIGEO%NativeElemID(Element))  !---- 2=target element
       PosCount(iProc) = PosCount(iProc) + 1
       SendContent(iProc)%content(PosCount(iProc):PosCount(iProc)+8)= &  !---- 3-11=Energy(1-9)
                        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Energy(1:9)
       SendContent(iProc)%content(PosCount(iProc)+9:PosCount(iProc)+11)= &  !---- 12-14=Force(1-3)
                        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Force(1:3)
       SendContent(iProc)%content(PosCount(iProc)+12)= &  !---- 15=Counter(1)
                        SampWallHaloCell(SurfMesh%HaloSideIDToSurfSideMap(iSide))%Counter(1)
       PosCount(iProc) = PosCount(iProc)+12
    END IF
  END DO
!---- comunication
  DO iProc=0, PMPIVAR%nProcs-1          
    IF (PMPIVAR%iProc.LT.iProc) THEN
      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
          15*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,IERROR)    
      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, &
          15*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR) 
    ELSE IF (PMPIVAR%iProc.GT.iProc) THEN
      IF (RecvMsgSurfs(iProc).NE.0) CALL MPI_RECV(RecvContent(iProc)%content, & 
          15*RecvMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,MPISTATUS,IERROR)  
      IF (SendMsgSurfs(iProc).NE.0) CALL MPI_SEND(SendContent(iProc)%content, &
          15*SendMsgSurfs(iProc),MPI_DOUBLE_PRECISION,iProc,1101,PMPIVAR%COMM,IERROR)    
    END IF
  END DO

!---- sum up surface values  
  DO iProc=0, PMPIVAR%nProcs-1
!    IF (RecvMsgSurfs(iProc).ne.0) print*,RecvMsgSurfs(iProc)
    IF ((iProc.NE.PMPIVAR%iProc).AND.(RecvMsgSurfs(iProc).NE.0)) THEN
      DO iCount = 1, RecvMsgSurfs(iProc)
      SurfSideID = SurfMesh%GlobSideToSurfSideMap( &  !---- get local surfaceID (iLocSide(target) = iLocSide(local))
        ElemToSide(1,INT(RecvContent(iProc)%content(iCount*15-14)), &
        INT(RecvContent(iProc)%content(iCount*15-13)))) 
      SampWall(SurfSideID)%Energy(1:9) = SampWall(SurfSideID)%Energy(1:9) &
        + RecvContent(iProc)%content(iCount*15-12:iCount*15-4)
      SampWall(SurfSideID)%Force(1:3) = SampWall(SurfSideID)%Force(1:3) &
        + RecvContent(iProc)%content(iCount*15-3:iCount*15-1)
!      IF (PMPIVAR%iProc.eq.0) then
!        print*, SurfSideID
!        print*, SampWall(SurfSideID)%Counter(1),RecvContent(iProc)%content(iCount*15)
!      end if
      SampWall(SurfSideID)%Counter(1) = SampWall(SurfSideID)%Counter(1) &
        + RecvContent(iProc)%content(iCount*15)      
      END DO
    END IF
  END DO  

  DO iProc=0, PMPIVAR%nProcs-1
    IF (ASSOCIATED(SendContent(iProc)%content)) DEALLOCATE(SendContent(iProc)%content)
    IF (ASSOCIATED(RecvContent(iProc)%content)) DEALLOCATE(RecvContent(iProc)%content)
  END DO
  DEALLOCATE(SendMsgSurfs)
  DEALLOCATE(RecvMsgSurfs)
  DEALLOCATE(SendContent)
  DEALLOCATE(RecvContent)
  DEALLOCATE(PosCount)

END SUBROUTINE MPISurfaceValuesSend

#endif
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!


!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcTVib(ChaTVib,MeanEVib,nMax)
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : DSMC
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  REAL, INTENT(IN)                :: ChaTVib,MeanEVib  ! Charak TVib, mean vibrational Energy of all molecules
  INTEGER, INTENT(IN)             :: nMax              ! INT(CharaTDisss/CharaTVib) + 1 
  REAL(KIND=8)                    :: LowerVal, UpperVal, MiddleVal, MaxPosiVal  ! upper and lower value of zero point search 
  REAl(KIND=8)                    :: eps_prec=1.0e-5   ! precision of zero point search
  REAL(KIND=8)                    :: ZeroVal1, ZeroVal2 ! both fuction values to compare

!--------------------------------------------------------------------------------------------------!
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
        WRITE(*,*) 'Cannot find zero point in TVib Calculation Function!'
        WRITE(*,*) 'System stopped!'
        STOP
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

!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcTelec(MeanEelec, iSpec)
  USE MOD_Particle_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : DSMC, SpecDSMC
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                              !
! Local variable declaration                                                                       !
  REAL, INTENT(IN)                :: MeanEelec  ! Charak TVib, mean vibrational Energy of all molecules
  INTEGER, INTENT(IN)             :: iSpec      ! Number of Species
  INTEGER                         :: ii
  REAL(KIND=8)                    :: LowerTemp, UpperTemp, MiddleTemp ! upper and lower value of modified
                                                                      ! zero point search
  REAl(KIND=8)                    :: eps_prec=1.0e-5   ! precision of zero point search
  REAL(KIND=8)                    :: SumOne, SumTwo    ! both summs
  REAL(KIND=8)                    :: ZeroVal1, ZeroVal2 ! both fuction values to compare

!--------------------------------------------------------------------------------------------------!
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
      DO ii = 0, SpecDSMC(iSpec)%MaxElecQuant
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

!--------------------------------------------------------------------------------------------------!

END MODULE MOD_DSMC_Analyze
