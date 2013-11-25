MODULE MOD_LD_Analyze
!==================================================================================================
! module for LD Sampling and Output
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
PUBLIC :: LD_data_sampling, LD_output_calc
!===================================================================================================

CONTAINS

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_data_sampling()

  USE MOD_LD_Vars
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, PartStateIntEn, CollisMode, SpecDSMC
  USE MOD_Particle_Vars,          ONLY : PEM, PDM, PartSpecies, PartState, PartMPF, usevMPF
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
                                                       + PartStateBulkValues(iPart,1:3) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1)   = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1) &  ! V2 stands for BulkTemp
                                                       + PartStateBulkValues(iPart,4) * PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%SimPartNum  = SampDSMC(iElem,PartSpecies(iPart))%SimPartNum + 1
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartStateBulkValues(iPart,1:3)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1)   = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1) & ! V2 stands for BulkTemp
                                                       + PartStateBulkValues(iPart,4)
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
      END IF
    END IF
  END DO

END SUBROUTINE LD_data_sampling

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_output_calc(nOutput)

  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, MacroDSMC, CollisMode, SpecDSMC, realtime
  USE MOD_Mesh_Vars,              ONLY : nElems,MeshFile
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst, Species, GEO, PartSpecies, usevMPF
  USE MOD_DSMC_Analyze,           ONLY : WriteDSMCToHDF5
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
      IF(SampDSMC(iElem,iSpec)%PartNum.GT. 0) THEN
! compute flow velocity
        MacroDSMC(iElem,iSpec)%PartV(1) = SampDSMC(iElem,iSpec)%PartV(1) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(2) = SampDSMC(iElem,iSpec)%PartV(2) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(3) = SampDSMC(iElem,iSpec)%PartV(3) / SampDSMC(iElem,iSpec)%PartNum
        MacroDSMC(iElem,iSpec)%PartV(4) = SQRT( MacroDSMC(iElem,iSpec)%PartV(1)**2 &
                                              + MacroDSMC(iElem,iSpec)%PartV(2)**2 &
                                              + MacroDSMC(iElem,iSpec)%PartV(3)**2 )
! compute flow Temperature
        MacroDSMC(iElem,iSpec)%Temp(1:3) = 0.0 
        MacroDSMC(iElem,iSpec)%Temp(4)   = SampDSMC(iElem,iSpec)%PartV2(1) / SampDSMC(iElem,iSpec)%PartNum
! compute density
        MacroDSMC(iElem,iSpec)%PartNum = SampDSMC(iElem,iSpec)%PartNum / REAL(DSMC%SampNum)
        ! Comment: if usevMPF MacroDSMC(iElem,iSpec)%PartNum == real number of particles
        IF (usevMPF) THEN
          MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum / GEO%Volume(iElem)
        ELSE 
          MacroDSMC(iElem,iSpec)%NumDens = MacroDSMC(iElem,iSpec)%PartNum * Species(iSpec)%MacroParticleFactor / GEO%Volume(iElem)
        END IF
! compute internal energies / has to be changed for vfd 
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
          (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
          IF (DSMC%VibEnergyModel.EQ.0) THEN
            IF (TVib_TempFac.LE.DSMC%GammaQuant) THEN
              MacroDSMC(iElem,iSpec)%TVib =  MacroDSMC(iElem,iSpec)%Temp(4)          
            ELSE
              MacroDSMC(iElem,iSpec)%TVib = MacroDSMC(iElem,iSpec)%Temp(4)
            END IF
          ELSE                                            ! TSHO-model
              MacroDSMC(iElem,iSpec)%TVib = MacroDSMC(iElem,iSpec)%Temp(4)
          END IF       
          MacroDSMC(iElem,iSpec)%TRot = MacroDSMC(iElem,iSpec)%Temp(4)
        IF (DSMC%ElectronicState) THEN
          MacroDSMC(iElem,iSpec)%TElec= MacroDSMC(iElem,iSpec)%Temp(4)
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
        MacroDSMC(iElem,nSpecies + 1)%Temp(1)  = 0.0
        MacroDSMC(iElem,nSpecies + 1)%Temp(2)  = 0.0
        MacroDSMC(iElem,nSpecies + 1)%Temp(3)  = 0.0
        MacroDSMC(iElem,nSpecies + 1)%Temp(4)  = MacroDSMC(iElem,nSpecies + 1)%Temp(4) &
                                               + MacroDSMC(iElem,iSpec)%Temp(4) * MacroDSMC(iElem,iSpec)%PartNum
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
      MacroDSMC(iElem,nSpecies + 1)%Temp(4)  = MacroDSMC(iElem,nSpecies + 1)%Temp(4)  / MacroDSMC(iElem,nSpecies + 1)%PartNum
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
!  CALL WriteOutputDSMC(nOutput)
  DEALLOCATE(MacroDSMC)
  
END SUBROUTINE LD_output_calc

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_Analyze
