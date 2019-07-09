!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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


MODULE MOD_LD_Analyze
!===================================================================================================================================
! module for LD Sampling and Output Calculation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE LD_data_sampling
  MODULE PROCEDURE LD_data_sampling
END INTERFACE
INTERFACE LD_output_calc
  MODULE PROCEDURE LD_output_calc
END INTERFACE

PUBLIC :: LD_data_sampling, LD_output_calc
!===================================================================================================================================
CONTAINS

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_data_sampling()
!===================================================================================================================================
! Sample celldata for output
!===================================================================================================================================
! MODULES
  USE MOD_LD_Vars
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, PartStateIntEn, CollisMode, SpecDSMC
  USE MOD_Particle_Vars,          ONLY : PEM, PDM, PartSpecies, PartMPF, usevMPF
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
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
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1) * PartMPF(iPart)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2) * PartMPF(iPart)
        END IF
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + PartMPF(iPart)
        SampDSMC(iElem,PartSpecies(iPart))%SimPartNum  = SampDSMC(iElem,PartSpecies(iPart))%SimPartNum + 1
      ELSE ! normal sampling without weighting
        SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3)  = SampDSMC(iElem,PartSpecies(iPart))%PartV(1:3) &
                                                       + PartStateBulkValues(iPart,1:3)
        SampDSMC(iElem,PartSpecies(iPart))%PartV2(1)   = SampDSMC(iElem,PartSpecies(iPart))%PartV2(1) & ! V2 stands for BulkTemp
                                                       + PartStateBulkValues(iPart,4)
        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
            (SpecDSMC(PartSpecies(iPart))%InterID.EQ.2)) THEN
          SampDSMC(iElem,PartSpecies(iPart))%EVib      = SampDSMC(iElem,PartSpecies(iPart))%EVib &
                                                       + PartStateIntEn(iPart,1)
          SampDSMC(iElem,PartSpecies(iPart))%ERot      = SampDSMC(iElem,PartSpecies(iPart))%ERot &
                                                       + PartStateIntEn(iPart,2)
        END IF
        SampDSMC(iElem,PartSpecies(iPart))%PartNum     = SampDSMC(iElem,PartSpecies(iPart))%PartNum + 1
      END IF
    END IF
  END DO

END SUBROUTINE LD_data_sampling

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_output_calc()
!===================================================================================================================================
! Calculation of outputdata on the basis of sampled values
!===================================================================================================================================
! MODULES
  USE MOD_Timedisc_Vars,          ONLY : time
  USE MOD_Globals_Vars,           ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,              ONLY : DSMC, SampDSMC, MacroDSMC, CollisMode, SpecDSMC
  USE MOD_Mesh_Vars,              ONLY : nElems,MeshFile
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species,  usevMPF
  USE MOD_Particle_Mesh_Vars,     ONLY : GEO
  USE MOD_DSMC_Analyze,           ONLY : WriteDSMCToHDF5, CalcTVib
  USE MOD_LD_Vars,                ONLY : LD_Residual, LD_CalcResidual
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: iElem, iSpec                  !
REAL                          :: TVib_TempFac, HeavyPartNum, MolecPartNum
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
        END IF

!!!!!!!!!!!!! compute internal energies / has to be changed for vfd
!!!!!!!!!!!!        IF (((CollisMode.EQ.2).OR.(CollisMode.EQ.3)).AND.&
!!!!!!!!!!!!          (SpecDSMC(iSpec)%InterID.EQ.2)) THEN
!!!!!!!!!!!!          MacroDSMC(iElem,iSpec)%TVib = MacroDSMC(iElem,iSpec)%Temp(4)
!!!!!!!!!!!!          MacroDSMC(iElem,iSpec)%TRot = MacroDSMC(iElem,iSpec)%Temp(4)
!!!!!!!!!!!!          IF (DSMC%ElectronicModel) THEN
!!!!!!!!!!!!            MacroDSMC(iElem,iSpec)%TElec= MacroDSMC(iElem,iSpec)%Temp(4)
!!!!!!!!!!!!          END IF
!!!!!!!!!!!!        END IF

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
        IF ( DSMC%ElectronicModel .AND. (SpecDSMC(iSpec)%InterID.NE.4) ) THEN
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
      IF ( DSMC%ElectronicModel ) THEN
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

    IF (LD_CalcResidual) THEN! residual calculation
      IF(LD_Residual(iElem,1).NE.0.0) THEN
        IF(MacroDSMC(iElem,nSpecies + 1)%PartV(1).GT. 1E-3) THEN
          LD_Residual(iElem,1) = ABS( (LD_Residual(iElem,1) - MacroDSMC(iElem,nSpecies + 1)%PartV(1)) &
                                      /MacroDSMC(iElem,nSpecies + 1)%PartV(1) )
          ! LD_Residual(iElem,1) = LD_Residual(iElem,1).....keep residual due to very small values
        END IF
        IF(MacroDSMC(iElem,nSpecies + 1)%PartV(2).GT. 1E-3) THEN
          LD_Residual(iElem,2) = ABS( (LD_Residual(iElem,2) - MacroDSMC(iElem,nSpecies + 1)%PartV(2)) &
                                      /MacroDSMC(iElem,nSpecies + 1)%PartV(2) )
        END IF
        IF(MacroDSMC(iElem,nSpecies + 1)%PartV(3).GT. 1E-3) THEN
          LD_Residual(iElem,3) = ABS( (LD_Residual(iElem,2) - MacroDSMC(iElem,nSpecies + 1)%PartV(3)) &
                                      /MacroDSMC(iElem,nSpecies + 1)%PartV(3) )
        END IF
        IF(MacroDSMC(iElem,nSpecies + 1)%PartV(4).GT. 1E-3) THEN
          LD_Residual(iElem,4) = ABS( (LD_Residual(iElem,3) - MacroDSMC(iElem,nSpecies + 1)%PartV(4)) &
                                      /MacroDSMC(iElem,nSpecies + 1)%PartV(4) )
        END IF
        LD_Residual(iElem,5) = ABS( (LD_Residual(iElem,4) - MacroDSMC(iElem,nSpecies + 1)%PartNum) &
                                    /MacroDSMC(iElem,nSpecies + 1)%PartNum )
        LD_Residual(iElem,6) = ABS( (LD_Residual(iElem,5) - MacroDSMC(iElem,nSpecies + 1)%Temp(4)) &
                                    /MacroDSMC(iElem,nSpecies + 1)%Temp(4) )
      ELSE
        LD_Residual(iElem,1) = 1.0
        LD_Residual(iElem,2) = 1.0
        LD_Residual(iElem,3) = 1.0
        LD_Residual(iElem,4) = 1.0
        LD_Residual(iElem,5) = 1.0
        LD_Residual(iElem,6) = 1.0
      END IF
    END IF ! end residual calculation
  END DO

  IF (LD_CalcResidual) THEN! residual output
    CALL LD_ResidualOutout
    DO iElem = 1, nElems ! element/cell main loop
      LD_Residual(iElem,1) = MacroDSMC(iElem,nSpecies + 1)%PartV(1)
      LD_Residual(iElem,2) = MacroDSMC(iElem,nSpecies + 1)%PartV(2)
      LD_Residual(iElem,3) = MacroDSMC(iElem,nSpecies + 1)%PartV(3)
      LD_Residual(iElem,4) = MacroDSMC(iElem,nSpecies + 1)%PartV(4)
      LD_Residual(iElem,5) = MacroDSMC(iElem,nSpecies + 1)%PartNum
      LD_Residual(iElem,6) = MacroDSMC(iElem,nSpecies + 1)%Temp(4)
    END DO
  END IF


  CALL WriteDSMCToHDF5(TRIM(MeshFile),time)
  DEALLOCATE(MacroDSMC)

END SUBROUTINE LD_output_calc
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE LD_ResidualOutout
!===================================================================================================================================
! Calculation of LD residuum values
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_TimeDisc_Vars,         ONLY : Time
  USE MOD_LD_Vars,               ONLY : LD_Residual
  USE MOD_Mesh_Vars,             ONLY : nElems
!--------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------!
! Local variable declaration                                                                       !
!--------------------------------------------------------------------------------------------------!
  LOGICAL                     :: isOpen
  CHARACTER(LEN=350)          :: outfile
  INTEGER                     :: unit_index
  REAL                        :: LD_ResidualTemp_Linf(6)
  REAL                        :: LD_ResidualTemp_L1(6)
!--------------------------------------------------------------------------------------------------!

    unit_index = 12345
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
      outfile = 'LD_Residual.csv'
      OPEN(unit_index,file=TRIM(outfile))
      !--- insert header
      WRITE(unit_index,'(A8)',ADVANCE='NO') '001-TIME'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_MAX_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloX_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloY_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloZ_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Velo_Abs_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Dens_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Temp_Linf'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_MAX_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloX_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloY_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_VeloZ_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Velo_Abs_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Dens_L1'
      WRITE(unit_index,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index,'(A14)',ADVANCE='NO') 'L_Temp_L1'
      WRITE(unit_index,'(A14)') ' '
    END IF

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') 'Write LD_Residual....'
    LD_ResidualTemp_Linf = MAXVAL(LD_Residual,1)
    LD_ResidualTemp_L1  = SUM(LD_Residual,1) / nElems

    WRITE(unit_index,104,ADVANCE='NO') Time
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') MAXVAL(LD_ResidualTemp_Linf)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(1)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(2)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(3)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(4)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(5)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_Linf(6)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') MAXVAL(LD_ResidualTemp_L1)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(1)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(2)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(3)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(4)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(5)
    WRITE(unit_index,'(A1)',ADVANCE='NO') ','
    WRITE(unit_index,104,ADVANCE='NO') LD_ResidualTemp_L1(6)
    WRITE(unit_index,'(A14)') ' '

104    FORMAT (e25.14)

END SUBROUTINE LD_ResidualOutout

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_Analyze
