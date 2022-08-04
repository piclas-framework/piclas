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

MODULE MOD_Particle_Restart
!===================================================================================================================================
! Module to handle PICLas's restart
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ParticleRestart
  MODULE PROCEDURE ParticleRestart
END INTERFACE

INTERFACE FinalizeParticleRestart
  MODULE PROCEDURE FinalizeParticleRestart
END INTERFACE

PUBLIC :: ParticleRestart
PUBLIC :: FinalizeParticleRestart
!===================================================================================================================================

CONTAINS

SUBROUTINE ParticleRestart()
!===================================================================================================================================
! Restart particle tracking
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Readin
USE MOD_Particle_Restart_Vars
! DSMC
USE MOD_DSMC_Vars              ,ONLY: UseDSMC,CollisMode,PartStateIntEn,DSMC,VibQuantsPar,PolyatomMolDSMC,SpecDSMC,RadialWeighting
USE MOD_DSMC_Vars              ,ONLY: ElectronicDistriPart, AmbipolElecVelo
! Localization
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemEpsOneCell
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod,NbrOfLostParticles,CountNbrOfLostParts
USE MOD_Particle_Tracking_Vars ,ONLY: NbrOfLostParticlesTotal,TotalNbrOfMissingParticlesSum,NbrOfLostParticlesTotal_old
! Interpolation
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_Eval_XYZ               ,ONLY: GetPositionInRefElem
! Mesh
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Mesh_Vars              ,ONLY: offsetElem
! Particles
USE MOD_HDF5_Input_Particles   ,ONLY: ReadEmissionVariablesFromHDF5
USE MOD_Part_Operations        ,ONLY: RemoveAllElectrons
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition,StoreLostParticleProperties
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Vars          ,ONLY: PartInt,PartData,PartState,PartSpecies,PEM,PDM,usevMPF,PartMPF,PartPosRef,SpecReset,Species
! Restart
USE MOD_Restart_Vars           ,ONLY: DoMacroscopicRestart
! HDG
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: UseBRElectronFluid,BRConvertElectronsToFluid,BRConvertFluidToElectrons
USE MOD_Part_BR_Elecron_Fluid  ,ONLY: CreateElectronsFromBRFluid
#endif /*USE_HDG*/
! MPI
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
USE MOD_Particle_Vars          ,ONLY: PartDataSize
#endif /*USE_MPI*/
USE MOD_Particle_Vars          ,ONLY: VibQuantData,ElecDistriData,AD_Data
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! Parameters
INTEGER,PARAMETER                  :: ELEM_FirstPartInd = 1
INTEGER,PARAMETER                  :: ELEM_LastPartInd  = 2
! Counters
INTEGER                            :: iElem
INTEGER                            :: FirstElemInd,LastelemInd
INTEGER(KIND=IK)                   :: locnPart,offsetnPart,iLoop
INTEGER                            :: CounterPoly
INTEGER                            :: iPolyatMole,iPart,CounterElec,CounterAmbi,SpecID
! Localization
LOGICAL                            :: InElementCheck
REAL                               :: xi(3)
REAL                               :: det(6,2)
INTEGER                            :: NbrOfMissingParticles
! MPI
#if USE_MPI
INTEGER,ALLOCATABLE                :: IndexOfFoundParticles(:),CompleteIndexOfFoundParticles(:)
INTEGER                            :: CompleteNbrOfLost,CompleteNbrOfFound,CompleteNbrOfDuplicate
REAL, ALLOCATABLE                  :: RecBuff(:,:)
INTEGER                            :: TotalNbrOfMissingParticles(0:PartMPI%nProcs-1), Displace(0:PartMPI%nProcs-1),CurrentPartNum
INTEGER                            :: OffsetTotalNbrOfMissingParticles(0:PartMPI%nProcs-1)
INTEGER                            :: NbrOfFoundParts, RecCount(0:PartMPI%nProcs-1)
INTEGER, ALLOCATABLE               :: SendBuffPoly(:), RecBuffPoly(:)
REAL, ALLOCATABLE                  :: SendBuffAmbi(:), RecBuffAmbi(:), SendBuffElec(:), RecBuffElec(:)
INTEGER                            :: LostPartsPoly(0:PartMPI%nProcs-1), DisplacePoly(0:PartMPI%nProcs-1)
INTEGER                            :: LostPartsElec(0:PartMPI%nProcs-1), DisplaceElec(0:PartMPI%nProcs-1)
INTEGER                            :: LostPartsAmbi(0:PartMPI%nProcs-1), DisplaceAmbi(0:PartMPI%nProcs-1)
INTEGER                            :: iProc
#endif /*USE_MPI*/
!===================================================================================================================================

! ===========================================================================
! Distribute or read the particle solution
! ===========================================================================
CALL ParticleReadin()

IF(.NOT.DoMacroscopicRestart) THEN
  IF(PartIntExists)THEN
    IF(PartDataExists)THEN
      iPart        = 0
      FirstElemInd = offsetElem+1
      LastElemInd  = offsetElem+PP_nElems
      locnPart     = PartInt(ELEM_LastPartInd,LastElemInd)-PartInt(ELEM_FirstPartInd,FirstElemInd)
      offsetnPart  = PartInt(ELEM_FirstPartInd,FirstElemInd)

      DO iLoop = 1_IK,locnPart
        ! Sanity check: SpecID > 0
        SpecID=INT(PartData(7,offsetnPart+iLoop),4)
        IF(SpecID.LE.0)THEN
          IPWRITE(UNIT_StdOut,'(I0,A,I0,A,I0,A,3(ES25.14E3),A,I0)') "Warning: Found particle in restart file with SpecID =", &
              SpecID,"for iLoop=",iLoop,"with pos: ",PartData(1:3,offsetnPart+iLoop)," and offsetnPart:",offsetnPart
          CALL abort(__STAMP__,'Found particle in restart file with species ID zero, which indicates a corrupted restart file.')
          !IPWRITE(UNIT_StdOut,'(I0,A,I0,A,I0)') "Warning: Found particle in restart file with SpecID =", SpecID,", which will be deleted. iLoop=",iLoop
          CYCLE
        END IF ! SpecID

        ! Check if species is to be removed during restart
        IF(SpecReset(SpecID)) CYCLE

        iPart = iPart + 1
        PartState(1,iPart) = PartData(1,offsetnPart+iLoop)
        PartState(2,iPart) = PartData(2,offsetnPart+iLoop)
        PartState(3,iPart) = PartData(3,offsetnPart+iLoop)
        PartState(4,iPart) = PartData(4,offsetnPart+iLoop)
        PartState(5,iPart) = PartData(5,offsetnPart+iLoop)
        PartState(6,iPart) = PartData(6,offsetnPart+iLoop)
        PartSpecies(iPart) = SpecID

        IF (useDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)) THEN
            PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
            PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
            PartStateIntEn(3,iPart)=PartData(10,offsetnPart+iLoop)
            ! Check if MPF was read from .h5 (or restarting with vMPF from non-vMPF restart file)
            IF(readVarFromState(11))THEN
              PartMPF(iPart)=PartData(11,offsetnPart+iLoop)
            ELSE
              PartMPF(iPart)=Species(SpecID)%MacroParticleFactor
            END IF ! readVarFromState(11)
          ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
            PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
            PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
            ! Check if MPF was read from .h5 (or restarting with vMPF from non-vMPF restart file)
            IF(readVarFromState(10))THEN
              PartMPF(iPart)=PartData(10,offsetnPart+iLoop)
            ELSE
              PartMPF(iPart)=Species(SpecID)%MacroParticleFactor
            END IF ! readVarFromState(10)
          ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel.GT.0)) THEN
            PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
            PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
            PartStateIntEn(3,iPart)=PartData(10,offsetnPart+iLoop)
          ELSE IF (CollisMode.GT.1) THEN
            IF (readVarFromState(8).AND.readVarFromState(9)) THEN
              PartStateIntEn(1,iPart)=PartData(8,offsetnPart+iLoop)
              PartStateIntEn(2,iPart)=PartData(9,offsetnPart+iLoop)
            ELSE IF ((SpecDSMC(PartSpecies(iPart))%InterID.EQ.1).OR.&
                     (SpecDSMC(PartSpecies(iPart))%InterID.EQ.10).OR.&
                     (SpecDSMC(PartSpecies(iPart))%InterID.EQ.15)) THEN
              !- setting inner DOF to 0 for atoms
              PartStateIntEn(1,iPart)=0.
              PartStateIntEn(2,iPart)=0.
            ELSE
              IPWRITE(UNIT_StdOut,*) "SpecDSMC(PartSpecies(iPart))%InterID =", SpecDSMC(PartSpecies(iPart))%InterID
              IPWRITE(UNIT_StdOut,*) "SpecID =", SpecID
              IPWRITE(UNIT_StdOut,*) "iPart =", iPart
              CALL Abort(__STAMP__,"resetting inner DOF for molecules is not implemented yet!")
            END IF ! readVarFromState(8).AND.readVarFromState(9)
          ELSE IF (usevMPF) THEN
            PartMPF(iPart)=PartData(8,offsetnPart+iLoop)
          END IF ! (CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)

          ! Polyatomic
          IF (DSMC%NumPolyatomMolecs.GT.0) THEN
            IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
              iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
              SDEALLOCATE(VibQuantsPar(iPart)%Quants)
              ALLOCATE(   VibQuantsPar(iPart)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
              VibQuantsPar(iPart)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)= &
                  VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,offsetnPart+iLoop)
            END IF ! SpecDSMC(PartSpecies(iPart))%PolyatomicMol
          END IF ! DSMC%NumPolyatomMolecs.GT.0

          ! Electronic
          IF (DSMC%ElectronicModel.EQ.2) THEN
            IF (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) THEN
              SDEALLOCATE(ElectronicDistriPart(iPart)%DistriFunc)
              ALLOCATE(   ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(PartSpecies(iPart))%MaxElecQuant))
              ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(PartSpecies(iPart))%MaxElecQuant)= &
              ElecDistriData(1:SpecDSMC(PartSpecies(iPart))%MaxElecQuant,offsetnPart+iLoop)
            END IF
          END IF

          ! Ambipolar Diffusion
          IF (DSMC%DoAmbipolarDiff) THEN
            IF (Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
              SDEALLOCATE(AmbipolElecVelo(iPart)%ElecVelo)
              ALLOCATE(   AmbipolElecVelo(iPart)%ElecVelo(1:3))
              AmbipolElecVelo(iPart)%ElecVelo(1:3)= AD_Data(1:3,offsetnPart+iLoop)
            END IF
          END IF
        ELSE IF (usevMPF) THEN
          PartMPF(iPart)=PartData(8,offsetnPart+iLoop)
        END IF ! useDSMC

        PDM%ParticleInside(iPart) = .TRUE.
      END DO ! iLoop = 1_IK,locnPart

      iPart = 0

      DO iElem=FirstElemInd,LastElemInd
        IF (PartInt(ELEM_LastPartInd,iElem).GT.PartInt(ELEM_FirstPartInd,iElem)) THEN
          DO iLoop = PartInt(ELEM_FirstPartInd,iElem)-offsetnPart+1_IK , PartInt(ELEM_LastPartInd,iElem)- offsetnPart

            ! Sanity check: SpecID > 0
            SpecID=INT(PartData(7,offsetnPart+iLoop),4)

            IF(SpecID.LE.0) CYCLE

            ! Check if species is to be removed during restart
            IF(SpecReset(SpecID)) CYCLE

            iPart                       = iPart +1
            PEM%GlobalElemID(iPart)     = iElem
            PEM%LastGlobalElemID(iPart) = iElem
          END DO ! iLoop
        END IF ! PartInt(ELEM_LastPartInd,iElem).GT.PartInt(ELEM_FirstPartInd,iElem)
      END DO ! iElem=FirstElemInd,LastElemInd

      DEALLOCATE(PartData)
      IF (useDSMC) THEN
        ! Polyatomic
        IF (DSMC%NumPolyatomMolecs.GT.0) THEN
          SDEALLOCATE(VibQuantData)
        END IF
        ! Electronic
        IF (DSMC%ElectronicModel.EQ.2) THEN
          SDEALLOCATE(ElecDistriData)
        END IF
        ! Ambipolar Diffusion
        IF (DSMC%DoAmbipolarDiff) THEN
          SDEALLOCATE(AD_Data)
        END IF
      END IF ! useDSMC
    ELSE ! not PartDataExists
      SWRITE(UNIT_stdOut,*)'PartData does not exists in restart file'
    END IF ! PartDataExists
    DEALLOCATE(PartInt)
    SDEALLOCATE(readVarFromState)

    PDM%ParticleVecLength = PDM%ParticleVecLength + iPart
    CALL UpdateNextFreePosition()
    LBWRITE(UNIT_stdOut,*)' DONE!'

    ! if ParticleVecLength GT maxParticleNumber: Stop
    IF (PDM%ParticleVecLength.GT.PDM%maxParticleNumber) THEN
      SWRITE (UNIT_stdOut,*) "PDM%ParticleVecLength =", PDM%ParticleVecLength
      SWRITE (UNIT_stdOut,*) "PDM%maxParticleNumber =", PDM%maxParticleNumber
      CALL abort(__STAMP__&
          ,' Number of Particles in Restart file is higher than MaxParticleNumber! Increase MaxParticleNumber!')
    END IF ! PDM%ParticleVecLength.GT.PDM%maxParticleNumber

    ! Since the elementside-local node number are NOT persistant and dependent on the location
    ! of the MPI borders, all particle-element mappings need to be checked after a restart
    ! Step 1: Identify particles that are not in the element in which they were before the restart
    NbrOfMissingParticles = 0
    NbrOfLostParticles    = 0
    CounterPoly           = 0
    CounterElec           = 0
    CounterAmbi           = 0

    SELECT CASE(TrackingMethod)
      CASE(TRIATRACKING)
        DO iPart = 1,PDM%ParticleVecLength
          ! Check if particle is inside the correct element
          CALL ParticleInsideQuad3D(PartState(1:3,iPart),PEM%GlobalElemID(iPart),InElementCheck,det)

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
#if !(USE_MPI)
              IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%GlobalElemID(iPart), UsePartState_opt=.TRUE.)
#endif /*!(USE_MPI)*/

              IF (useDSMC) THEN
                ! Polyatomic
                IF (DSMC%NumPolyatomMolecs.GT.0) THEN
                  IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
                    iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
                    CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                  END IF
                END IF
                ! Electronic
                IF (DSMC%ElectronicModel.EQ.2) THEN
                  IF (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) &
                    CounterElec = CounterElec + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
                END IF
                ! Ambipolar Diffusion
                IF (DSMC%DoAmbipolarDiff) THEN
                  IF (Species(PartSpecies(iPart))%ChargeIC.GT.0.0) &
                    CounterAmbi = CounterAmbi + 3
                END IF
              END IF ! useDSMC
            ELSE
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
            END IF
          END IF
        END DO ! iPart = 1,PDM%ParticleVecLength

      CASE(TRACING)
        DO iPart = 1,PDM%ParticleVecLength
          ! Check if particle is inside the correct element
          CALL GetPositionInRefElem(PartState(1:3,iPart),Xi,PEM%GlobalElemID(iPart))
          IF (ALL(ABS(Xi).LE.1.0)) THEN ! particle inside
            InElementCheck = .TRUE.
            IF(ALLOCATED(PartPosRef)) PartPosRef(1:3,iPart)=Xi
          ELSE
            InElementCheck = .FALSE.
          END IF

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
#if !(USE_MPI)
              IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%GlobalElemID(iPart), UsePartState_opt=.TRUE.)
#endif /*!(USE_MPI)*/

              IF (useDSMC) THEN
                ! Polyatomic
                IF (DSMC%NumPolyatomMolecs.GT.0) THEN
                  IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
                    iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
                    CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                  END IF
                END IF
                ! Electronic
                IF (DSMC%ElectronicModel.EQ.2) THEN
                  IF (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) &
                    CounterElec = CounterElec + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
                END IF
                ! Ambipolar Diffusion
                IF (DSMC%DoAmbipolarDiff) THEN
                  IF (Species(PartSpecies(iPart))%ChargeIC.GT.0.0) &
                    CounterAmbi = CounterAmbi + 3
                END IF
              END IF ! useDSMC
            ELSE
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
            END IF ! .NOT.PDM%ParticleInside(iPart)
          END IF ! .NOT.InElementCheck
        END DO ! iPart = 1,PDM%ParticleVecLength

      CASE(REFMAPPING)
        DO iPart = 1,PDM%ParticleVecLength
          ! Check if particle is inside the correct element
          CALL GetPositionInRefElem(PartState(1:3,iPart),Xi,PEM%GlobalElemID(iPart))
          IF (ALL(ABS(Xi).LE.ElemEpsOneCell(GetCNElemID(PEM%GlobalElemID(iPart))))) THEN ! particle inside
            InElementCheck    = .TRUE.
            PartPosRef(1:3,iPart) = Xi
          ELSE
            InElementCheck    = .FALSE.
          END IF

          ! Particle not in correct element, try to find them within MyProc
          IF (.NOT.InElementCheck) THEN
            NbrOfMissingParticles = NbrOfMissingParticles + 1
            CALL LocateParticleInElement(iPart,doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (.NOT.PDM%ParticleInside(iPart)) THEN
              NbrOfLostParticles = NbrOfLostParticles + 1
#if !(USE_MPI)
              IF (CountNbrOfLostParts) CALL StoreLostParticleProperties(iPart, PEM%GlobalElemID(iPart), UsePartState_opt=.TRUE.)
#endif /*!(USE_MPI)*/

              IF (useDSMC) THEN
                ! Polyatomic
                IF (DSMC%NumPolyatomMolecs.GT.0) THEN
                  IF (SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
                    iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
                    CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
                  END IF
                END IF
                ! Electronic
                IF (DSMC%ElectronicModel.EQ.2) THEN
                  IF (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) &
                    CounterElec = CounterElec + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
                END IF
                ! Ambipolar Diffusion
                IF (DSMC%DoAmbipolarDiff) THEN
                  IF (Species(PartSpecies(iPart))%ChargeIC.GT.0.0) &
                    CounterAmbi = CounterAmbi + 3
                END IF
              END IF ! useDSMC
              PartPosRef(1:3,iPart) = -888.
            ELSE
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
            END IF
          END IF
        END DO ! iPart = 1,PDM%ParticleVecLength
    END SELECT

#if USE_MPI
    ! Step 2: All particles that are not found within MyProc need to be communicated to the others and located there
    ! Combine number of lost particles of all processes and allocate variables
    ! Note: Particles that are lost on MyProc are also searched for here again
    CALL MPI_ALLGATHER(NbrOfLostParticles, 1, MPI_INTEGER, TotalNbrOfMissingParticles, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
    NbrOfLostParticles=0
    IF (useDSMC) THEN
      IF (DSMC%NumPolyatomMolecs.GT.0) CALL MPI_ALLGATHER(CounterPoly, 1, MPI_INTEGER, LostPartsPoly, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
      IF (DSMC%ElectronicModel.EQ.2)   CALL MPI_ALLGATHER(CounterElec, 1, MPI_INTEGER, LostPartsElec, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
      IF (DSMC%DoAmbipolarDiff)        CALL MPI_ALLGATHER(CounterAmbi, 1, MPI_INTEGER, LostPartsAmbi, 1, MPI_INTEGER, PartMPI%COMM, IERROR)
    END IF ! useDSMC

    !TotalNbrOfMissingParticlesSum = SUM(INT(TotalNbrOfMissingParticles,8))
    TotalNbrOfMissingParticlesSum = SUM(TotalNbrOfMissingParticles)

    ! Check total number of missing particles and start re-locating them on other procs
    IF (TotalNbrOfMissingParticlesSum.GT.0) THEN

      ! Set offsets
      OffsetTotalNbrOfMissingParticles(0) = 0
      DO iProc = 1, PartMPI%nProcs-1
        OffsetTotalNbrOfMissingParticles(iProc) = OffsetTotalNbrOfMissingParticles(iProc-1) + TotalNbrOfMissingParticles(iProc-1)
      END DO ! iProc = 0, PartMPI%nProcs-1

      ALLOCATE(RecBuff(PartDataSize,1:TotalNbrOfMissingParticlesSum))
      IF (useDSMC) THEN
        ! Polyatomic
        IF (DSMC%NumPolyatomMolecs.GT.0) THEN
          ALLOCATE(SendBuffPoly(1:CounterPoly))
          ALLOCATE(RecBuffPoly( 1:SUM(LostPartsPoly)))
        END IF
        ! Electronic
        IF (DSMC%ElectronicModel.EQ.2) THEN
          ALLOCATE(SendBuffElec(1:CounterElec))
          ALLOCATE(RecBuffElec( 1:SUM(LostPartsElec)))
        END IF
        ! Ambipolar Diffusion
        IF (DSMC%DoAmbipolarDiff) THEN
          ALLOCATE(SendBuffAmbi(1:CounterAmbi))
          ALLOCATE(RecBuffAmbi( 1:SUM(LostPartsAmbi)))
        END IF
      END IF ! useDSMC

      ! Fill SendBuffer
      NbrOfMissingParticles = OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + 1
      CounterPoly = 0
      CounterAmbi = 0
      CounterElec = 0
      DO iPart = 1, PDM%ParticleVecLength
        IF (.NOT.PDM%ParticleInside(iPart)) THEN
          RecBuff(1:6,NbrOfMissingParticles) = PartState(1:6,iPart)
          RecBuff(7,NbrOfMissingParticles)   = REAL(PartSpecies(iPart))
          IF (useDSMC) THEN
            IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)) THEN
              RecBuff(8,NbrOfMissingParticles)  = PartStateIntEn(1,iPart)
              RecBuff(9,NbrOfMissingParticles)  = PartStateIntEn(2,iPart)
              RecBuff(10,NbrOfMissingParticles) = PartMPF(iPart)
              RecBuff(11,NbrOfMissingParticles) = PartStateIntEn(3,iPart)
            ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
              RecBuff(8,NbrOfMissingParticles)  = PartStateIntEn(1,iPart)
              RecBuff(9,NbrOfMissingParticles)  = PartStateIntEn(2,iPart)
              RecBuff(10,NbrOfMissingParticles) = PartMPF(iPart)
            ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel.GT.0)) THEN
              RecBuff(8,NbrOfMissingParticles)  = PartStateIntEn(1,iPart)
              RecBuff(9,NbrOfMissingParticles)  = PartStateIntEn(2,iPart)
              RecBuff(10,NbrOfMissingParticles) = PartStateIntEn(3,iPart)
            ELSE IF (CollisMode.GT.1) THEN
              RecBuff(8,NbrOfMissingParticles)  = PartStateIntEn(1,iPart)
              RecBuff(9,NbrOfMissingParticles)  = PartStateIntEn(2,iPart)
            ELSE IF (usevMPF) THEN
              RecBuff(8,NbrOfMissingParticles)  = PartMPF(iPart)
            END IF
          ELSE IF (usevMPF) THEN
            RecBuff(8,NbrOfMissingParticles) = PartMPF(iPart)
          END IF
          NbrOfMissingParticles = NbrOfMissingParticles + 1

          !--- receive the polyatomic vibquants per particle at the end of the message
          ! Polyatomic
          IF(useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
            IF(SpecDSMC(PartSpecies(iPart))%PolyatomicMol) THEN
              iPolyatMole = SpecDSMC(PartSpecies(iPart))%SpecToPolyArray
              SendBuffPoly(CounterPoly+1:CounterPoly+PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                  = VibQuantsPar(iPart)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
              CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
            END IF ! SpecDSMC(PartSpecies(iPart))%PolyatomicMol
          END IF ! useDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)
          !--- receive the polyatomic vibquants per particle at the end of the message
          ! Electronic
          IF(useDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
            IF (.NOT.((SpecDSMC(PartSpecies(iPart))%InterID.EQ.4).OR.SpecDSMC(PartSpecies(iPart))%FullyIonized)) THEN
              SendBuffElec(CounterElec+1:CounterElec+SpecDSMC(PartSpecies(iPart))%MaxElecQuant) &
                  = ElectronicDistriPart(iPart)%DistriFunc(1:SpecDSMC(PartSpecies(iPart))%MaxElecQuant)
              CounterElec = CounterElec + SpecDSMC(PartSpecies(iPart))%MaxElecQuant
            END IF !
          END IF !
          ! Ambipolar Diffusion
          IF(useDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
            IF (Species(PartSpecies(iPart))%ChargeIC.GT.0.0)  THEN
              SendBuffAmbi(CounterAmbi+1:CounterAmbi+3) = AmbipolElecVelo(iPart)%ElecVelo(1:3)
              CounterAmbi = CounterAmbi + 3
            END IF !
          END IF !

        END IF ! .NOT.PDM%ParticleInside(iPart)
      END DO ! iPart = 1, PDM%ParticleVecLength

      ! Distribute lost particles to all procs
      NbrOfMissingParticles = 0
      CounterPoly = 0
      CounterElec = 0
      CounterAmbi = 0

      DO iProc = 0, PartMPI%nProcs-1
        RecCount(iProc) = TotalNbrOfMissingParticles(iProc)
        Displace(iProc) = NbrOfMissingParticles
        NbrOfMissingParticles = NbrOfMissingParticles + TotalNbrOfMissingParticles(iProc)

        IF (useDSMC) THEN
          ! Polyatomic
          IF (DSMC%NumPolyatomMolecs.GT.0) THEN
            DisplacePoly(iProc) = CounterPoly
            CounterPoly = CounterPoly + LostPartsPoly(iProc)
          END IF
          ! Electronic
          IF (DSMC%ElectronicModel.EQ.2) THEN
            DisplaceElec(iProc) = CounterElec
            CounterElec = CounterElec + LostPartsElec(iProc)
          END IF
          ! Ambipolar Diffusion
          IF (DSMC%DoAmbipolarDiff) THEN
            DisplaceAmbi(iProc) = CounterAmbi
            CounterAmbi = CounterAmbi + LostPartsAmbi(iProc)
          END IF
        END IF ! useDSMC
      END DO ! iProc = 0, PartMPI%nProcs-1

      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                         , 0                                                &
                         , MPI_DATATYPE_NULL                                &
                         , RecBuff                                          &
                         , PartDataSize*TotalNbrOfMissingParticles(:)       &
                         , PartDataSize*OffsetTotalNbrOfMissingParticles(:) &
                         , MPI_DOUBLE_PRECISION                             &
                         , PartMPI%COMM                                     &
                         , IERROR)

      IF (useDSMC) THEN
        ! Polyatomic
        IF (DSMC%NumPolyatomMolecs.GT.0) CALL MPI_ALLGATHERV(SendBuffPoly, LostPartsPoly(PartMPI%MyRank), MPI_INTEGER, &
            RecBuffPoly, LostPartsPoly, DisplacePoly, MPI_INTEGER, PartMPI%COMM, IERROR)
        ! Electronic
        IF (DSMC%ElectronicModel.EQ.2)   CALL MPI_ALLGATHERV(SendBuffElec, LostPartsElec(PartMPI%MyRank), MPI_INTEGER, &
            RecBuffElec, LostPartsElec, DisplaceElec, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)
        ! Ambipolar Diffusion
        IF (DSMC%DoAmbipolarDiff)        CALL MPI_ALLGATHERV(SendBuffAmbi, LostPartsAmbi(PartMPI%MyRank), MPI_INTEGER, &
            RecBuffAmbi, LostPartsAmbi, DisplaceAmbi, MPI_DOUBLE_PRECISION, PartMPI%COMM, IERROR)
      END IF

      ! Keep track which particles are found on the current proc
      ALLOCATE(IndexOfFoundParticles          (TotalNbrOfMissingParticlesSum))
      IF (MPIRoot) &
        ALLOCATE(CompleteIndexOfFoundParticles(TotalNbrOfMissingParticlesSum))
      IndexOfFoundParticles = -1

      ! Free lost particle positions in local array to make room for missing particles that are tested
      CALL UpdateNextFreePosition()

      ! Add them to particle list and check if they are in MyProcs domain
      NbrOfFoundParts = 0
      CurrentPartNum  = PDM%ParticleVecLength+1
      CounterPoly = 0
      CounterElec = 0
      CounterAmbi = 0

      DO iPart = 1, TotalNbrOfMissingParticlesSum
        ! Sanity check
        IF(CurrentPartNum.GT.PDM%maxParticleNumber)THEn
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " CurrentPartNum        = ",  CurrentPartNum
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " PDM%maxParticleNumber = ",  PDM%maxParticleNumber
          CALL abort(__STAMP__,'Missing particle ID > PDM%maxParticleNumber. Increase Part-MaxParticleNumber!')
        END IF !CurrentPartNum.GT.PDM%maxParticleNumber

        ! Do not search particles twice: Skip my own particles, because these have already been searched for before they are
        ! sent to all other procs
        ASSOCIATE( myFirst => OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + 1 ,&
                   myLast  => OffsetTotalNbrOfMissingParticles(PartMPI%MyRank) + TotalNbrOfMissingParticles(PartMPI%MyRank))
          IF((iPart.GE.myFirst).AND.(iPart.LE.myLast))THEN
            IndexOfFoundParticles(iPart) = 0
            CYCLE
          END IF
        END ASSOCIATE

        PartState(     1:6,CurrentPartNum) = RecBuff(1:6,iPart)
        PDM%ParticleInside(CurrentPartNum) = .true.

        CALL LocateParticleInElement(CurrentPartNum,doHALO=.FALSE.)
        IF (PDM%ParticleInside(CurrentPartNum)) THEN
          IndexOfFoundParticles(iPart) = 1
          PEM%LastGlobalElemID(CurrentPartNum) = PEM%GlobalElemID(CurrentPartNum)

          ! Set particle properties (if the particle is lost, it's properties are written to a .h5 file)
          PartSpecies(CurrentPartNum) = INT(RecBuff(7,iPart))
          IF (useDSMC) THEN
            IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0)) THEN
              PartStateIntEn(1,CurrentPartNum) = RecBuff(8,iPart)
              PartStateIntEn(2,CurrentPartNum) = RecBuff(9,iPart)
              PartStateIntEn(3,CurrentPartNum) = RecBuff(11,iPart)
              PartMPF(CurrentPartNum)          = RecBuff(10,iPart)
            ELSE IF ((CollisMode.GT.1).AND. (usevMPF)) THEN
              PartStateIntEn(1,CurrentPartNum) = RecBuff(8,iPart)
              PartStateIntEn(2,CurrentPartNum) = RecBuff(9,iPart)
              PartMPF(CurrentPartNum)          = RecBuff(10,iPart)
            ELSE IF ((CollisMode.GT.1).AND. (DSMC%ElectronicModel.GT.0)) THEN
              PartStateIntEn(1,CurrentPartNum) = RecBuff(8,iPart)
              PartStateIntEn(2,CurrentPartNum) = RecBuff(9,iPart)
              PartStateIntEn(3,CurrentPartNum) = RecBuff(10,iPart)
            ELSE IF (CollisMode.GT.1) THEN
              PartStateIntEn(1,CurrentPartNum) = RecBuff(8,iPart)
              PartStateIntEn(2,CurrentPartNum) = RecBuff(9,iPart)
            ELSE IF (usevMPF) THEN
              PartMPF(CurrentPartNum)          = RecBuff(8,iPart)
            END IF
          ELSE IF (usevMPF) THEN
            PartMPF(CurrentPartNum)          = RecBuff(8,iPart)
          END IF ! useDSMC
          NbrOfFoundParts = NbrOfFoundParts + 1

          ! Check if particle was found inside of an element
          IF (useDSMC) THEN
            ! Polyatomic
            IF (DSMC%NumPolyatomMolecs.GT.0) THEN
              IF(SpecDSMC(PartSpecies(CurrentPartNum))%PolyatomicMol) THEN
                iPolyatMole = SpecDSMC(PartSpecies(CurrentPartNum))%SpecToPolyArray
                SDEALLOCATE(VibQuantsPar(CurrentPartNum)%Quants)
                ALLOCATE(VibQuantsPar(CurrentPartNum)%Quants(PolyatomMolDSMC(iPolyatMole)%VibDOF))
                VibQuantsPar(CurrentPartNum)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
                    = RecBuffPoly(CounterPoly+1:CounterPoly+PolyatomMolDSMC(iPolyatMole)%VibDOF)
                CounterPoly = CounterPoly + PolyatomMolDSMC(iPolyatMole)%VibDOF
              END IF
            END IF
            ! Electronic
            IF (DSMC%ElectronicModel.EQ.2) THEN
              IF (.NOT.((SpecDSMC(PartSpecies(CurrentPartNum))%InterID.EQ.4) &
                  .OR.SpecDSMC(PartSpecies(CurrentPartNum))%FullyIonized)) THEN
                SDEALLOCATE(ElectronicDistriPart(CurrentPartNum)%DistriFunc)
                ALLOCATE(ElectronicDistriPart(CurrentPartNum)%DistriFunc(1:SpecDSMC(PartSpecies(CurrentPartNum))%MaxElecQuant))
                ElectronicDistriPart(CurrentPartNum)%DistriFunc(1:SpecDSMC(PartSpecies(CurrentPartNum))%MaxElecQuant)= &
                  RecBuffElec(CounterElec+1:CounterElec+SpecDSMC(PartSpecies(CurrentPartNum))%MaxElecQuant)
                CounterElec = CounterElec +SpecDSMC(PartSpecies(CurrentPartNum))%MaxElecQuant
              END IF
            END IF
            ! Ambipolar Diffusion
            IF (DSMC%DoAmbipolarDiff) THEN
              IF (Species(PartSpecies(CurrentPartNum))%ChargeIC.GT.0.0) THEN
                SDEALLOCATE(AmbipolElecVelo(CurrentPartNum)%ElecVelo)
                ALLOCATE(AmbipolElecVelo(CurrentPartNum)%ElecVelo(1:3))
                AmbipolElecVelo(CurrentPartNum)%ElecVelo(1:3)= RecBuffAmbi(CounterAmbi+1:CounterAmbi+3)
                CounterAmbi = CounterAmbi + 3
              END IF
            END IF
          END IF ! useDSMC

          CurrentPartNum = CurrentPartNum + 1
        ELSE ! Lost
          IndexOfFoundParticles(iPart) = 0
        END IF

        ! Sanity Check
        IF(IndexOfFoundParticles(iPart).EQ.-1)THEN
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " iPart                        : ",  iPart
          IPWRITE(UNIT_StdOut,'(I0,A,I0)') " IndexOfFoundParticles(iPart) : ",  IndexOfFoundParticles(iPart)
          CALL abort(__STAMP__,'IndexOfFoundParticles(iPart) was not set correctly)')
        END IF ! IndexOfFoundParticles(iPart)
      END DO ! iPart = 1, TotalNbrOfMissingParticlesSum

      PDM%ParticleVecLength = PDM%ParticleVecLength + NbrOfFoundParts

      ! Combine number of found particles to make sure none are lost completely or found twice
      IF(MPIroot)THEN
        CALL MPI_REDUCE(IndexOfFoundParticles,CompleteIndexOfFoundParticles,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
      ELSE
        CALL MPI_REDUCE(IndexOfFoundParticles,0                            ,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,PartMPI%COMM,IERROR)
      END IF

      CompleteNbrOfFound      = 0
      CompleteNbrOfLost       = 0
      CompleteNbrOfDuplicate  = 0
      ! Only mpi root: Check if particle are not found or found twice
      IF (MPIRoot) THEN
        DO iPart = 1,TotalNbrOfMissingParticlesSum
          IF    (CompleteIndexOfFoundParticles(iPart).EQ.0) THEN ! permanently lost
            CompleteNbrOfLost       = CompleteNbrOfLost      + 1
          ELSEIF(CompleteIndexOfFoundParticles(iPart).EQ.1) THEN ! lost and found on one other processor
            CompleteNbrOfFound      = CompleteNbrOfFound     + 1
          ELSE ! lost and found multiple times on other processors
            CompleteNbrOfDuplicate  = CompleteNbrOfDuplicate + 1
          END IF

          ! Store the particle info
          IF(CountNbrOfLostParts)THEN
            CurrentPartNum = PDM%ParticleVecLength+1

            ! Set properties of the "virtual" particle (only for using the routine StoreLostParticleProperties to store this info
            ! in the .h5 container)
            PartState(1:6,CurrentPartNum)        = RecBuff(1:6,iPart)
            PartSpecies(CurrentPartNum)          = INT(RecBuff(7,iPart))
            PEM%LastGlobalElemID(CurrentPartNum) = -1
            PDM%ParticleInside(CurrentPartNum)   = .FALSE.
            IF(usevMPF) PartMPF(CurrentPartNum)  = RecBuff(8,iPart) ! only required when using vMPF

            CALL StoreLostParticleProperties(CurrentPartNum, PEM%GlobalElemID(CurrentPartNum), &
                                             UsePartState_opt=.TRUE., PartMissingType_opt=CompleteIndexOfFoundParticles(iPart))
          END IF ! CountNbrOfLostParts
        END DO

        WRITE(UNIT_stdOut,'(A,I0)') ' Particles initially lost during restart    : ',TotalNbrOfMissingParticlesSum
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found on (other) procs : ',CompleteNbrOfFound
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles permanently lost       : ',CompleteNbrOfLost
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found multiple times   : ',CompleteNbrOfDuplicate
        NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + CompleteNbrOfLost

        DEALLOCATE(CompleteIndexOfFoundParticles)
      END IF ! MPIRoot

      CALL MPI_BCAST(NbrOfLostParticlesTotal,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
      NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal
    END IF ! TotalNbrOfMissingParticlesSum.GT.0

#else /*not USE_MPI*/
    IF(NbrOfMissingParticles.GT.0)THEN
      NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + NbrOfLostParticles
      TotalNbrOfMissingParticlesSum = NbrOfMissingParticles
      NbrOfLostParticlesTotal_old = NbrOfLostParticlesTotal
      WRITE(UNIT_stdOut,'(A,I0)') ' Particles initially lost during restart : ',NbrOfMissingParticles
      WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles permanently lost    : ',NbrOfLostParticles
      WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles relocated           : ',NbrOfMissingParticles-NbrOfLostParticles
      NbrOfLostParticles = 0
    END IF ! NbrOfMissingParticles.GT.0
#endif /*USE_MPI*/

    CALL UpdateNextFreePosition()

    ! Read-in the stored cloned particles
    IF (RadialWeighting%PerformCloning) CALL RestartClones()
  ELSE ! not PartIntExists
    SWRITE(UNIT_stdOut,*)'PartInt does not exists in restart file'
  END IF ! PartIntExists
ELSE ! DoMacroscopicRestart
  CALL MacroscopicRestart()
  CALL UpdateNextFreePosition()
END IF ! .NOT.DoMacroscopicRestart

! Read-in the cell-local wall temperature
IF (ANY(PartBound%UseAdaptedWallTemp)) CALL RestartAdaptiveWallTemp()

#if USE_HDG
! Remove electron species when using BR electron fluid model
IF(UseBRElectronFluid.AND.BRConvertElectronsToFluid) CALL RemoveAllElectrons()
#endif /*USE_HDG*/

! This routines opens the data file in single mode (i.e. only the root opens and closes the data file)
! ------------------------------------------------
! Particle Emission Parameters
! ------------------------------------------------
CALL ReadEmissionVariablesFromHDF5()

#if USE_HDG
  ! Create electrons from BR fluid properties
  IF(BRConvertFluidToElectrons) CALL CreateElectronsFromBRFluid(.TRUE.)
#endif /*USE_HDG*/

END SUBROUTINE ParticleRestart


SUBROUTINE RestartClones()
!===================================================================================================================================
! Axisymmetric 2D simulation with particle weighting: Read-in of clone particles saved during output of particle data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_input
USE MOD_io_hdf5
USE MOD_Mesh_Vars         ,ONLY: offsetElem, nElems
USE MOD_DSMC_Vars         ,ONLY: useDSMC, CollisMode, DSMC, PolyatomMolDSMC, SpecDSMC
USE MOD_DSMC_Vars         ,ONLY: RadialWeighting, ClonedParticles
USE MOD_Particle_Vars     ,ONLY: nSpecies, usevMPF, Species
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nDimsClone, CloneDataSize, ClonePartNum, iPart, iDelay, maxDelay, iElem, tempDelay
INTEGER(HSIZE_T), POINTER :: SizeClone(:)
REAL,ALLOCATABLE          :: CloneData(:,:)
INTEGER                   :: iPolyatmole, MaxQuantNum, iSpec, compareDelay, MaxElecQuant
INTEGER,ALLOCATABLE       :: pcount(:), VibQuantData(:,:)
REAL, ALLOCATABLE         :: ElecDistriData(:,:), AD_Data(:,:)
LOGICAL                   :: CloneExists
!===================================================================================================================================

CALL DatasetExists(File_ID,'CloneData',CloneExists)
IF(.NOT.CloneExists) THEN
  LBWRITE(*,*) 'No clone data found! Restart without cloning.'
  IF(RadialWeighting%CloneMode.EQ.1) THEN
    RadialWeighting%CloneDelayDiff = 1
  ELSEIF (RadialWeighting%CloneMode.EQ.2) THEN
    RadialWeighting%CloneDelayDiff = 0
  END IF ! RadialWeighting%CloneMode.EQ.1
  RETURN
END IF ! CloneExists

CALL GetDataSize(File_ID,'CloneData',nDimsClone,SizeClone)

CloneDataSize = INT(SizeClone(1),4)
ClonePartNum = INT(SizeClone(2),4)
DEALLOCATE(SizeClone)

IF(ClonePartNum.GT.0) THEN
  ALLOCATE(CloneData(1:CloneDataSize,1:ClonePartNum))
  ASSOCIATE(ClonePartNum  => INT(ClonePartNum,IK)  ,&
            CloneDataSize => INT(CloneDataSize,IK) )
    CALL ReadArray('CloneData',2,(/CloneDataSize,ClonePartNum/),0_IK,2,RealArray=CloneData)
  END ASSOCIATE
  LBWRITE(*,*) 'Read-in of cloned particles complete. Total clone number: ', ClonePartNum
  ! Determing the old clone delay
  maxDelay = INT(MAXVAL(CloneData(9,:)))
  IF(RadialWeighting%CloneMode.EQ.1) THEN
    ! Array is allocated from 0 to maxDelay
    compareDelay = maxDelay + 1
  ELSE
    compareDelay = maxDelay
  END IF
  IF(compareDelay.GT.RadialWeighting%CloneInputDelay) THEN
    LBWRITE(*,*) 'Old clone delay is greater than the new delay. Old delay:', compareDelay
    RadialWeighting%CloneDelayDiff = RadialWeighting%CloneInputDelay + 1
  ELSEIF(compareDelay.EQ.RadialWeighting%CloneInputDelay) THEN
    LBWRITE(*,*) 'The clone delay has not been changed.'
    RadialWeighting%CloneDelayDiff = RadialWeighting%CloneInputDelay + 1
  ELSE
    LBWRITE(*,*) 'New clone delay is greater than the old delay. Old delay:', compareDelay
    RadialWeighting%CloneDelayDiff = compareDelay + 1
  END IF
  IF(RadialWeighting%CloneMode.EQ.1) THEN
    tempDelay = RadialWeighting%CloneInputDelay - 1
  ELSE
    tempDelay = RadialWeighting%CloneInputDelay
  END IF
  ALLOCATE(pcount(0:tempDelay))
  pcount(0:tempDelay) = 0
  ! Polyatomic clones: determining the size of the VibQuant array
  IF (UseDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
    MaxQuantNum = 0
    DO iSpec = 1, nSpecies
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
      END IF
    END DO
    ALLOCATE(VibQuantData(1:MaxQuantNum,1:ClonePartNum))
    ASSOCIATE(ClonePartNum => INT(ClonePartNum,IK),MaxQuantNum => INT(MaxQuantNum,IK))
      CALL ReadArray('CloneVibQuantData',2,(/MaxQuantNum,ClonePartNum/),0_IK,2,IntegerArray_i4=VibQuantData)
    END ASSOCIATE
  END IF
  IF (UseDSMC.AND.(DSMC%ElectronicModel.EQ.2)) THEN
    MaxElecQuant = 0
    DO iSpec = 1, nSpecies
      IF (.NOT.((SpecDSMC(iSpec)%InterID.EQ.4).OR.SpecDSMC(iSpec)%FullyIonized)) THEN
        IF (SpecDSMC(iSpec)%MaxElecQuant.GT.MaxElecQuant) MaxElecQuant = SpecDSMC(iSpec)%MaxElecQuant
      END IF
    END DO
    ALLOCATE(ElecDistriData(1:MaxElecQuant,1:ClonePartNum))
    ASSOCIATE(ClonePartNum => INT(ClonePartNum,IK),MaxElecQuant => INT(MaxElecQuant,IK))
      CALL ReadArray('CloneElecDistriData',2,(/MaxElecQuant,ClonePartNum/),0_IK,2,RealArray=ElecDistriData)
    END ASSOCIATE
  END IF
  IF (UseDSMC.AND.DSMC%DoAmbipolarDiff) THEN
    ALLOCATE(AD_Data(1:3,1:ClonePartNum))
    ASSOCIATE(ClonePartNum => INT(ClonePartNum,IK))
      CALL ReadArray('CloneADVeloData',2,(/INT(3,IK),ClonePartNum/),0_IK,2,RealArray=AD_Data)
    END ASSOCIATE
  END IF
  ! Copying particles into ClonedParticles array
  DO iPart = 1, ClonePartNum
    iDelay = INT(CloneData(9,iPart))
    iElem = INT(CloneData(8,iPart)) - offsetElem
    IF((iElem.LE.nElems).AND.(iElem.GT.0)) THEN
      IF(iDelay.LE.tempDelay) THEN
        pcount(iDelay) = pcount(iDelay) + 1
        RadialWeighting%ClonePartNum(iDelay) = pcount(iDelay)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(1) = CloneData(1,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(2) = CloneData(2,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(3) = CloneData(3,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(4) = CloneData(4,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(5) = CloneData(5,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%PartState(6) = CloneData(6,iPart)
        ClonedParticles(pcount(iDelay),iDelay)%Species = INT(CloneData(7,iPart))
        ClonedParticles(pcount(iDelay),iDelay)%Element = INT(CloneData(8,iPart))
        ClonedParticles(pcount(iDelay),iDelay)%lastPartPos(1:3) = CloneData(1:3,iPart)
        IF (UseDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(3) = CloneData(12,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor   = CloneData(13,iPart)
          ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor   = CloneData(12,iPart)
          ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel.GT.0) ) THEN
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(3) = CloneData(12,iPart)
          ELSE IF (CollisMode.GT.1) THEN
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(1) = CloneData(10,iPart)
            ClonedParticles(pcount(iDelay),iDelay)%PartStateIntEn(2) = CloneData(11,iPart)
          ELSE IF (usevMPF) THEN
            ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor = CloneData(10,iPart)
          END IF
        ELSE IF (usevMPF) THEN
            ClonedParticles(pcount(iDelay),iDelay)%WeightingFactor = CloneData(10,iPart)
        END IF
        IF (UseDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
          IF (SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%PolyatomicMol) THEN
            iPolyatMole = SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%SpecToPolyArray
            ALLOCATE(ClonedParticles(pcount(iDelay),iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF))
            ClonedParticles(pcount(iDelay),iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF) &
              = VibQuantData(1:PolyatomMolDSMC(iPolyatMole)%VibDOF,iPart)
          END IF
        END IF
        IF (UseDSMC.AND.(DSMC%ElectronicModel.EQ.2))  THEN
          IF (.NOT.((SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%InterID.EQ.4) &
              .OR.SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%FullyIonized)) THEN
            ALLOCATE(ClonedParticles(pcount(iDelay),iDelay)%DistriFunc( &
                    1:SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%MaxElecQuant))
            ClonedParticles(pcount(iDelay),iDelay)%DistriFunc(1:SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%MaxElecQuant) &
              = ElecDistriData(1:SpecDSMC(ClonedParticles(pcount(iDelay),iDelay)%Species)%MaxElecQuant,iPart)
          END IF
        END IF
        IF (UseDSMC.AND.DSMC%DoAmbipolarDiff)  THEN
          IF (Species(ClonedParticles(pcount(iDelay),iDelay)%Species)%ChargeIC.GT.0.0) THEN
            ALLOCATE(ClonedParticles(pcount(iDelay),iDelay)%AmbiPolVelo(1:3))
            ClonedParticles(pcount(iDelay),iDelay)%AmbiPolVelo(1:3) = AD_Data(1:3,iPart)
          END IF
        END IF
      END IF
    END IF
  END DO
ELSE
  LBWRITE(*,*) 'Read-in of cloned particles complete. No clones detected.'
END IF

END SUBROUTINE RestartClones


SUBROUTINE RestartAdaptiveWallTemp()
!===================================================================================================================================
!> Read-in of the adaptive side-local wall temperature and the corresponding global side index
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_input
USE MOD_io_hdf5
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, nSurfTotalSides
USE MOD_Particle_Boundary_Vars    ,ONLY: BoundaryWallTemp, GlobalSide2SurfSide
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars           ,ONLY: MPI_COMM_LEADERS_SURF, MPI_COMM_SHARED
USE MOD_Particle_Boundary_Vars    ,ONLY: BoundaryWallTemp_Shared_Win
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE      :: tmpGlobalSideInx(:)
REAL, ALLOCATABLE         :: tmpWallTemp(:,:,:)
INTEGER                   :: iSide, tmpSide, iSurfSide
LOGICAL                   :: AdaptiveWallTempExists
!===================================================================================================================================

CALL DatasetExists(File_ID,'AdaptiveBoundaryWallTemp',AdaptiveWallTempExists)
IF (.NOT.AdaptiveWallTempExists) THEN
  SWRITE(*,*) 'No side-local temperature found. The wall temperature will be adapted during the next macroscopic output.'
  RETURN
END IF

CALL DatasetExists(File_ID,'AdaptiveBoundaryGlobalSideIndx',AdaptiveWallTempExists)
IF (.NOT.AdaptiveWallTempExists) THEN
  CALL Abort(__STAMP__,&
    'ERROR during Restart: AdaptiveBoundaryWallTemp was found in the restart file but not the GlobalSideIndx array!')
END IF

IF (nSurfTotalSides.EQ.0) RETURN

ALLOCATE(tmpGlobalSideInx(nSurfTotalSides), &
      tmpWallTemp(nSurfSample,nSurfSample,nSurfTotalSides))
! Associate construct for integer KIND=8 possibility
#if USE_MPI
! Return if not a sampling leader
IF (MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
#endif
  ASSOCIATE (&
        nSurfSample          => INT(nSurfSample,IK)                     , &
        nGlobalSides         => INT(nSurfTotalSides,IK))
    CALL ReadArray('AdaptiveBoundaryGlobalSideIndx',1,(/nGlobalSides/),0_IK,1,IntegerArray_i4=tmpGlobalSideInx)
    CALL ReadArray('AdaptiveBoundaryWallTemp',3,(/nSurfSample, nSurfSample, nGlobalSides/),0_IK,1,RealArray=tmpWallTemp)
  END ASSOCIATE

  DO iSide = 1, nSurfTotalSides
    tmpSide = tmpGlobalSideInx(iSide)
    IF (GlobalSide2SurfSide(SURF_SIDEID,tmpSide).EQ.-1) CYCLE
    iSurfSide = GlobalSide2SurfSide(SURF_SIDEID,tmpSide)
    BoundaryWallTemp(:,:,iSurfSide) = tmpWallTemp(:,:,iSide)
  END DO
#if USE_MPI
ELSE
  ASSOCIATE (&
        nSurfSample          => INT(0,IK)                     , &
        nGlobalSides         => INT(0,IK))
    CALL ReadArray('AdaptiveBoundaryGlobalSideIndx',1,(/nGlobalSides/),0_IK,1,IntegerArray_i4=tmpGlobalSideInx)
    CALL ReadArray('AdaptiveBoundaryWallTemp',3,(/nSurfSample, nSurfSample, nGlobalSides/),0_IK,1,RealArray=tmpWallTemp)
  END ASSOCIATE
END IF

CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
#endif

END SUBROUTINE RestartAdaptiveWallTemp


SUBROUTINE MacroscopicRestart()
!===================================================================================================================================
!> Read-in of the element data from a DSMC state and insertion of particles based on the macroscopic values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_HDF5_Input    ,ONLY: OpenDataFile,ReadArray,GetDataSize,ReadAttribute
USE MOD_HDF5_Input    ,ONLY: nDims,HSize,File_ID
USE MOD_Restart_Vars  ,ONLY: MacroRestartFileName, MacroRestartValues
USE MOD_Mesh_Vars     ,ONLY: offsetElem, nElems
USE MOD_Particle_Vars ,ONLY: nSpecies
USE MOD_Macro_Restart ,ONLY: MacroRestart_InsertParticles
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: nVar_HDF5, iVar, iSpec, iElem
REAL, ALLOCATABLE                 :: ElemData_HDF5(:,:)
CHARACTER(LEN=255)                :: File_Type
!===================================================================================================================================

SWRITE(UNIT_stdOut,*) 'Using macroscopic values from file: ',TRIM(MacroRestartFileName)

CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

! Check if the provided file is a DSMC state file.
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
IF(TRIM(File_Type).NE.'DSMCState') THEN
  SWRITE(*,*) 'ERROR: The given file type is: ', TRIM(File_Type)
  CALL abort(__STAMP__,&
      'ERROR: Given file for the macroscopic restart is not of the type "DSMCState", please check the input file!')
END IF

CALL GetDataSize(File_ID,'ElemData',nDims,HSize,attrib=.FALSE.)
nVar_HDF5=INT(HSize(1),4)
DEALLOCATE(HSize)

ALLOCATE(MacroRestartValues(1:nElems,1:nSpecies+1,1:DSMC_NVARS))
MacroRestartValues = 0.

ALLOCATE(ElemData_HDF5(1:nVar_HDF5,1:nElems))
! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
  nVar_HDF5  => INT(nVar_HDF5,IK) ,&
  offsetElem => INT(offsetElem,IK),&
  nElems     => INT(nElems,IK)    )
  CALL ReadArray('ElemData',2,(/nVar_HDF5,nElems/),offsetElem,2,RealArray=ElemData_HDF5(:,:))
END ASSOCIATE

iVar = 1
DO iSpec = 1, nSpecies
  DO iElem = 1, nElems
    MacroRestartValues(iElem,iSpec,:) = ElemData_HDF5(iVar:iVar-1+DSMC_NVARS,iElem)
  END DO
  iVar = iVar + DSMC_NVARS
END DO

CALL MacroRestart_InsertParticles()

DEALLOCATE(MacroRestartValues)
DEALLOCATE(ElemData_HDF5)

END SUBROUTINE MacroscopicRestart


SUBROUTINE FinalizeParticleRestart()
!===================================================================================================================================
! Finalizes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
#if USE_HDG
USE MOD_HDG_Vars     ,ONLY: BRConvertFluidToElectrons,BRConvertElectronsToFluid
#endif /*USE_HDG*/
! IMPLICIT VARIABLE HANDLINGDGInitIsDone
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#if USE_HDG
! Avoid converting BR electron fluid to actual particles or vice versa during an automatic load balance restart
BRConvertFluidToElectrons = .FALSE.
BRConvertElectronsToFluid = .FALSE.
#endif /*USE_HDG*/

END SUBROUTINE FinalizeParticleRestart

END MODULE MOD_Particle_Restart
