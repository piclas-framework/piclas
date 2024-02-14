!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
USE MOD_DSMC_Vars              ,ONLY: UseDSMC,CollisMode,PartStateIntEn,DSMC,VibQuantsPar,PolyatomMolDSMC,SpecDSMC
USE MOD_DSMC_Vars              ,ONLY: ElectronicDistriPart, AmbipolElecVelo
! Localization
USE MOD_Particle_Localization  ,ONLY: LocateParticleInElement,SinglePointToElement
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
USE MOD_Part_Operations        ,ONLY: RemoveAllElectrons, RemoveParticle
USE MOD_Part_Tools             ,ONLY: UpdateNextFreePosition,StoreLostParticleProperties, MergeCells
USE MOD_Particle_Boundary_Vars ,ONLY: PartBound
USE MOD_Particle_Vars          ,ONLY: PartInt,PartData,PartState,PartSpecies,PEM,PDM,usevMPF,PartMPF,PartPosRef,SpecReset,Species
USE MOD_Particle_Vars          ,ONLY: DoVirtualCellMerge
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
! Restart
USE MOD_Restart_Vars           ,ONLY: DoMacroscopicRestart, MacroRestartValues
! HDG
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: UseBRElectronFluid,BRConvertElectronsToFluid,BRConvertFluidToElectrons
USE MOD_Part_BR_Elecron_Fluid  ,ONLY: CreateElectronsFromBRFluid
#endif /*USE_HDG*/
! MPI
#if USE_MPI
USE MOD_Particle_Vars          ,ONLY: PartDataSize
#endif /*USE_MPI*/
USE MOD_Particle_Vars          ,ONLY: VibQuantData,ElecDistriData,AD_Data
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! Rotational frame of reference
USE MOD_Particle_Vars          ,ONLY: UseRotRefFrame, PartVeloRotRef, RotRefFrameOmega
USE MOD_Part_Tools             ,ONLY: InRotRefFrameCheck, IncreaseMaxParticleNumber
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
INTEGER                            :: iElem,iPos
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
INTEGER                            :: TotalNbrOfMissingParticles(0:nProcessors-1), Displace(0:nProcessors-1),CurrentPartNum
INTEGER                            :: OffsetTotalNbrOfMissingParticles(0:nProcessors-1)
INTEGER                            :: NbrOfFoundParts, RecCount(0:nProcessors-1)
INTEGER, ALLOCATABLE               :: SendBuffPoly(:), RecBuffPoly(:)
REAL, ALLOCATABLE                  :: SendBuffAmbi(:), RecBuffAmbi(:), SendBuffElec(:), RecBuffElec(:)
INTEGER                            :: LostPartsPoly(0:nProcessors-1), DisplacePoly(0:nProcessors-1)
INTEGER                            :: LostPartsElec(0:nProcessors-1), DisplaceElec(0:nProcessors-1)
INTEGER                            :: LostPartsAmbi(0:nProcessors-1), DisplaceAmbi(0:nProcessors-1)
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
      CALL IncreaseMaxParticleNumber(INT(locnPart))

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
        PartState(1:6,iPart) = PartData(1:6,offsetnPart+iLoop)
        PartSpecies(iPart) = SpecID

        iPos = 7
        ! Rotational frame of reference: initialize logical and velocity
        IF(UseRotRefFrame) THEN
          PDM%InRotRefFrame(iPart) = InRotRefFrameCheck(iPart)
          IF(PDM%InRotRefFrame(iPart)) THEN
            IF(readVarFromState(1+iPos).AND.readVarFromState(2+iPos).AND.readVarFromState(3+iPos)) THEN
              PartVeloRotRef(1:3,iPart) = PartData(MapPartDataToReadin(1+iPos):MapPartDataToReadin(3+iPos),offsetnPart+iLoop)
            ELSE
              PartVeloRotRef(1:3,iPart) = PartState(4:6,iPart) - CROSS(RotRefFrameOmega(1:3),PartState(1:3,iPart))
            END IF
          ELSE
            PartVeloRotRef(1:3,iPart) = 0.
          END IF
          iPos=iPos+3
        END IF

        IF(useDSMC) THEN
          IF(CollisMode.GT.1) THEN
            IF(readVarFromState(1+iPos).AND.readVarFromState(2+iPos)) THEN
              PartStateIntEn(1:2,iPart)=PartData(MapPartDataToReadin(1+iPos):MapPartDataToReadin(2+iPos),offsetnPart+iLoop)
              iPos=iPos+2
            ELSE IF((SpecDSMC(SpecID)%InterID.EQ.1).OR.(SpecDSMC(SpecID)%InterID.EQ.10).OR.(SpecDSMC(SpecID)%InterID.EQ.15)) THEN
              !- setting inner DOF to 0 for atoms
              PartStateIntEn(1:2,iPart) = 0.
            ELSE
              IPWRITE(UNIT_StdOut,*) "SpecDSMC(PartSpecies(iPart))%InterID =", SpecDSMC(PartSpecies(iPart))%InterID
              IPWRITE(UNIT_StdOut,*) "SpecID =", SpecID
              IPWRITE(UNIT_StdOut,*) "iPart =", iPart
              CALL Abort(__STAMP__,"resetting inner DOF for molecules is not implemented yet!")
            END IF ! readVarFromState
            IF(DSMC%ElectronicModel.GT.0) THEN
              PartStateIntEn(3,iPart)=PartData(MapPartDataToReadin(1+iPos),offsetnPart+iLoop)
              iPos=iPos+1
            END IF
          END IF
        END IF
        IF(usevMPF) THEN
          ! Check if MPF was read from .h5 (or restarting with vMPF from non-vMPF restart file)
          IF(readVarFromState(1+iPos))THEN
            PartMPF(iPart)=PartData(MapPartDataToReadin(1+iPos),offsetnPart+iLoop)
            iPos=iPos+1
          ELSE
            PartMPF(iPart)=Species(SpecID)%MacroParticleFactor
          END IF ! readVarFromState
        END IF
        IF(useDSMC) THEN
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
        END IF
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
    SDEALLOCATE(MapPartDataToReadin)

    PDM%ParticleVecLength = PDM%ParticleVecLength + iPart
#ifdef CODE_ANALYZE
    IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
    DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
      IF (PDM%ParticleInside(iPart)) THEN
        IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
        CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
      END IF
    END DO
#endif
    CALL UpdateNextFreePosition()
    LBWRITE(UNIT_stdOut,*)' DONE!'

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
            PEM%GlobalElemID(iPart) = SinglePointToElement(PartState(1:3,iPart),doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (PEM%GlobalElemID(iPart).GT.0) THEN
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
              PDM%ParticleInside(iPart)=.TRUE.
              IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),PEM%GlobalElemID(iPart))
            ELSE
              PDM%ParticleInside(iPart)=.FALSE.
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
            PEM%GlobalElemID(iPart) = SinglePointToElement(PartState(1:3,iPart),doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (PEM%GlobalElemID(iPart).GT.0) THEN
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
              PDM%ParticleInside(iPart)=.TRUE.
              IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),PEM%GlobalElemID(iPart))
            ELSE
              PDM%ParticleInside(iPart)=.FALSE.
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
            PEM%GlobalElemID(iPart) = SinglePointToElement(PartState(1:3,iPart),doHALO=.FALSE.)

            ! Particle not found within MyProc
            IF (PEM%GlobalElemID(iPart).GT.0) THEN
              PEM%LastGlobalElemID(iPart) = PEM%GlobalElemID(iPart)
              PDM%ParticleInside(iPart)=.TRUE.
              IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),PEM%GlobalElemID(iPart))
            ELSE
              PDM%ParticleInside(iPart)=.FALSE.
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
            END IF
          END IF
        END DO ! iPart = 1,PDM%ParticleVecLength
    END SELECT

#if USE_MPI
    ! Step 2: All particles that are not found within MyProc need to be communicated to the others and located there
    ! Combine number of lost particles of all processes and allocate variables
    ! Note: Particles that are lost on MyProc are also searched for here again
    CALL MPI_ALLGATHER(NbrOfLostParticles, 1, MPI_INTEGER, TotalNbrOfMissingParticles, 1, MPI_INTEGER, MPI_COMM_PICLAS, IERROR)
    NbrOfLostParticles=0
    IF (useDSMC) THEN
      IF (DSMC%NumPolyatomMolecs.GT.0) CALL MPI_ALLGATHER(CounterPoly, 1, MPI_INTEGER, LostPartsPoly, 1, MPI_INTEGER, MPI_COMM_PICLAS, IERROR)
      IF (DSMC%ElectronicModel.EQ.2)   CALL MPI_ALLGATHER(CounterElec, 1, MPI_INTEGER, LostPartsElec, 1, MPI_INTEGER, MPI_COMM_PICLAS, IERROR)
      IF (DSMC%DoAmbipolarDiff)        CALL MPI_ALLGATHER(CounterAmbi, 1, MPI_INTEGER, LostPartsAmbi, 1, MPI_INTEGER, MPI_COMM_PICLAS, IERROR)
    END IF ! useDSMC

    !TotalNbrOfMissingParticlesSum = SUM(INT(TotalNbrOfMissingParticles,8))
    TotalNbrOfMissingParticlesSum = SUM(TotalNbrOfMissingParticles)

    ! Check total number of missing particles and start re-locating them on other procs
    IF (TotalNbrOfMissingParticlesSum.GT.0) THEN

      ! Set offsets
      OffsetTotalNbrOfMissingParticles(0) = 0
      DO iProc = 1, nProcessors-1
        OffsetTotalNbrOfMissingParticles(iProc) = OffsetTotalNbrOfMissingParticles(iProc-1) + TotalNbrOfMissingParticles(iProc-1)
      END DO ! iProc = 0, nProcessors-1

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
      NbrOfMissingParticles = OffsetTotalNbrOfMissingParticles(myRank) + 1
      CounterPoly = 0
      CounterAmbi = 0
      CounterElec = 0
      DO iPart = 1, PDM%ParticleVecLength
        IF (.NOT.PDM%ParticleInside(iPart)) THEN
          RecBuff(1:6,NbrOfMissingParticles) = PartState(1:6,iPart)
          RecBuff(7,NbrOfMissingParticles)   = REAL(PartSpecies(iPart))
          iPos=7
          ! Rotational frame of reference
          IF(UseRotRefFrame) THEN
            RecBuff(1+iPos:3+iPos,NbrOfMissingParticles) = PartVeloRotRef(1:3,iPart)
            iPos = iPos + 3
          END IF
          IF (useDSMC) THEN
            IF (CollisMode.GT.1) THEN
              RecBuff(1+iPos:2+iPos,NbrOfMissingParticles)  = PartStateIntEn(1:2,iPart)
              iPos=iPos+2
              IF(DSMC%ElectronicModel.GT.0) THEN
                RecBuff(1+iPos,NbrOfMissingParticles) = PartStateIntEn(3,iPart)
                iPos=iPos+1
              END IF
            END IF
          END IF
          IF (usevMPF) THEN
            RecBuff(1+iPos,NbrOfMissingParticles) = PartMPF(iPart)
            iPos=iPos+1
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

      DO iProc = 0, nProcessors-1
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
      END DO ! iProc = 0, nProcessors-1

      CALL MPI_ALLGATHERV( MPI_IN_PLACE                                     &
                         , 0                                                &
                         , MPI_DATATYPE_NULL                                &
                         , RecBuff                                          &
                         , PartDataSize*TotalNbrOfMissingParticles(:)       &
                         , PartDataSize*OffsetTotalNbrOfMissingParticles(:) &
                         , MPI_DOUBLE_PRECISION                             &
                         , MPI_COMM_PICLAS                                     &
                         , IERROR)

      IF (useDSMC) THEN
        ! Polyatomic
        IF (DSMC%NumPolyatomMolecs.GT.0) CALL MPI_ALLGATHERV(SendBuffPoly, LostPartsPoly(myRank), MPI_INTEGER, &
            RecBuffPoly, LostPartsPoly, DisplacePoly, MPI_INTEGER, MPI_COMM_PICLAS, IERROR)
        ! Electronic
        IF (DSMC%ElectronicModel.EQ.2)   CALL MPI_ALLGATHERV(SendBuffElec, LostPartsElec(myRank), MPI_INTEGER, &
            RecBuffElec, LostPartsElec, DisplaceElec, MPI_DOUBLE_PRECISION, MPI_COMM_PICLAS, IERROR)
        ! Ambipolar Diffusion
        IF (DSMC%DoAmbipolarDiff)        CALL MPI_ALLGATHERV(SendBuffAmbi, LostPartsAmbi(myRank), MPI_INTEGER, &
            RecBuffAmbi, LostPartsAmbi, DisplaceAmbi, MPI_DOUBLE_PRECISION, MPI_COMM_PICLAS, IERROR)
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
        ASSOCIATE( myFirst => OffsetTotalNbrOfMissingParticles(myRank) + 1 ,&
                   myLast  => OffsetTotalNbrOfMissingParticles(myRank) + TotalNbrOfMissingParticles(myRank))
          IF((iPart.GE.myFirst).AND.(iPart.LE.myLast))THEN
            IndexOfFoundParticles(iPart) = 0
            CYCLE
          END IF
        END ASSOCIATE

        PartState(     1:6,CurrentPartNum) = RecBuff(1:6,iPart)

        PEM%GlobalElemID(CurrentPartNum) = SinglePointToElement(PartState(1:3,CurrentPartNum),doHALO=.FALSE.)

        IF (PEM%GlobalElemID(CurrentPartNum).GT.0) THEN
          PEM%LastGlobalElemID(CurrentPartNum) = PEM%GlobalElemID(CurrentPartNum)
          PDM%ParticleInside(CurrentPartNum)=.TRUE.
          IF(TrackingMethod.EQ.REFMAPPING) CALL GetPositionInRefElem(PartState(1:3,CurrentPartNum),PartPosRef(1:3,CurrentPartNum),PEM%GlobalElemID(iPart))
          IndexOfFoundParticles(iPart) = 1
          PEM%LastGlobalElemID(CurrentPartNum) = PEM%GlobalElemID(CurrentPartNum)

          ! Set particle properties (if the particle is lost, it's properties are written to a .h5 file)
          PartSpecies(CurrentPartNum) = INT(RecBuff(7,iPart))
          iPos = 7
          ! Rotational frame of reference
          IF(UseRotRefFrame) THEN
            PDM%InRotRefFrame(CurrentPartNum) = InRotRefFrameCheck(CurrentPartNum)
            IF(PDM%InRotRefFrame(CurrentPartNum)) THEN
              PartVeloRotRef(1:3,CurrentPartNum) = RecBuff(1+iPos:3+iPos,iPart)
            ELSE
              PartVeloRotRef(1:3,CurrentPartNum) = 0.
            END IF
            iPos = iPos + 3
          END IF
          ! DSMC-specific variables
          IF (useDSMC) THEN
            IF (CollisMode.GT.1) THEN
              PartStateIntEn(1:2,CurrentPartNum) = RecBuff(1+iPos:2+iPos,iPart)
              iPos = iPos + 2
              IF(DSMC%ElectronicModel.GT.0) THEN
                PartStateIntEn(3,CurrentPartNum) = RecBuff(1+iPos,iPart)
                iPos = iPos + 1
              END IF
            END IF
          END IF
          ! Variable particle weighting
          IF (usevMPF) THEN
            PartMPF(CurrentPartNum) = RecBuff(1+iPos,iPart)
            iPos = iPos + 1
          END IF
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
          PDM%ParticleInside(iPart)=.FALSE.
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
#ifdef CODE_ANALYZE
      IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL Abort(__STAMP__,'PDM%ParticleVeclength exceeds PDM%maxParticleNumber, Difference:',IntInfoOpt=PDM%ParticleVeclength-PDM%maxParticleNumber)
      DO iPart=PDM%ParticleVecLength+1,PDM%maxParticleNumber
        IF (PDM%ParticleInside(iPart)) THEN
          IPWRITE(*,*) iPart,PDM%ParticleVecLength,PDM%maxParticleNumber
          CALL Abort(__STAMP__,'Particle outside PDM%ParticleVeclength',IntInfoOpt=iPart)
        END IF
      END DO
#endif
      ! IF(PDM%ParticleVecLength.GT.PDM%maxParticleNumber) CALL IncreaseMaxParticleNumber(PDM%ParticleVecLength*CEILING(1+0.5*PDM%MaxPartNumIncrease)-PDM%maxParticleNumber)

      ! Combine number of found particles to make sure none are lost completely or found twice
      IF(MPIroot)THEN
        CALL MPI_REDUCE(IndexOfFoundParticles,CompleteIndexOfFoundParticles,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
      ELSE
        CALL MPI_REDUCE(IndexOfFoundParticles,0                            ,TotalNbrOfMissingParticlesSum,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
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
            CALL RemoveParticle(CurrentPartNum)
          END IF ! CountNbrOfLostParts
        END DO

        WRITE(UNIT_stdOut,'(A,I0)') ' Particles initially lost during restart    : ',TotalNbrOfMissingParticlesSum
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found on (other) procs : ',CompleteNbrOfFound
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles permanently lost       : ',CompleteNbrOfLost
        WRITE(UNIT_stdOut,'(A,I0)') ' Number of particles found multiple times   : ',CompleteNbrOfDuplicate
        NbrOfLostParticlesTotal = NbrOfLostParticlesTotal + CompleteNbrOfLost

        DEALLOCATE(CompleteIndexOfFoundParticles)
      END IF ! MPIRoot

      CALL MPI_BCAST(NbrOfLostParticlesTotal,1,MPI_INTEGER,0,MPI_COMM_PICLAS,iError)
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
  ELSE ! not PartIntExists
    SWRITE(UNIT_stdOut,*)'PartInt does not exists in restart file'
  END IF ! PartIntExists
ELSE ! DoMacroscopicRestart
  CALL MacroscopicRestart()
  CALL UpdateNextFreePosition()
END IF ! .NOT.DoMacroscopicRestart

! Read-in the cell-local wall temperature
IF (ANY(PartBound%UseAdaptedWallTemp)) CALL RestartAdaptiveWallTemp()

! Read-in of adaptive BC sampling values
IF(UseAdaptiveBC.OR.nPorousBC.GT.0) CALL RestartAdaptiveBCSampling()

#if USE_HDG
! Remove electron species when using BR electron fluid model
IF(UseBRElectronFluid.AND.BRConvertElectronsToFluid) CALL RemoveAllElectrons()
#endif /*USE_HDG*/

! This routines opens the data file in single mode (i.e. only the root opens and closes the data file)
! ------------------------------------------------
! Particle Emission Parameters
! ------------------------------------------------
CALL ReadEmissionVariablesFromHDF5()

IF (DoVirtualCellMerge) CALL MergeCells()

#if USE_HDG
  ! Create electrons from BR fluid properties
  IF(BRConvertFluidToElectrons) CALL CreateElectronsFromBRFluid(.TRUE.)
#endif /*USE_HDG*/

! Deallocate the read-in macroscopic values (might have been utilized in RestartAdaptiveBCSampling)
SDEALLOCATE(MacroRestartValues)

END SUBROUTINE ParticleRestart


SUBROUTINE RestartAdaptiveWallTemp()
!===================================================================================================================================
!> Read-in of the adaptive side-local wall temperature and the corresponding global side index
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_HDF5_input
USE MOD_io_hdf5
USE MOD_Restart_Vars              ,ONLY: RestartFile
USE MOD_Particle_Boundary_Vars    ,ONLY: nSurfSample, nGlobalSurfSides
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

! Leave routine if no surface sides have been defined in the domain
IF (nGlobalSurfSides.EQ.0) RETURN

#if USE_MPI
! Only the surface leaders open the file
IF (MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
  CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_LEADERS_SURF)
END IF
#else
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
#endif

! Associate construct for integer KIND=8 possibility
#if USE_MPI
! Only the surface leaders read the array
IF (MPI_COMM_LEADERS_SURF.NE.MPI_COMM_NULL) THEN
#endif
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

  ALLOCATE(tmpGlobalSideInx(nGlobalSurfSides),tmpWallTemp(nSurfSample,nSurfSample,nGlobalSurfSides))

  ASSOCIATE (nSurfSample          => INT(nSurfSample,IK), &
             nGlobalSides         => INT(nGlobalSurfSides,IK))
    CALL ReadArray('AdaptiveBoundaryGlobalSideIndx',1,(/nGlobalSides/),0_IK,1,IntegerArray_i4=tmpGlobalSideInx)
    CALL ReadArray('AdaptiveBoundaryWallTemp',3,(/nSurfSample, nSurfSample, nGlobalSides/),0_IK,1,RealArray=tmpWallTemp)
  END ASSOCIATE
  ! Mapping of the temperature on the global side to the node-local surf side
  DO iSide = 1, nGlobalSurfSides
    tmpSide = tmpGlobalSideInx(iSide)
    IF (GlobalSide2SurfSide(SURF_SIDEID,tmpSide).EQ.-1) CYCLE
    iSurfSide = GlobalSide2SurfSide(SURF_SIDEID,tmpSide)
    BoundaryWallTemp(:,:,iSurfSide) = tmpWallTemp(:,:,iSide)
  END DO
#if USE_MPI
END IF
! Distribute the temperature distribution onto the shared array
CALL BARRIER_AND_SYNC(BoundaryWallTemp_Shared_Win,MPI_COMM_SHARED)
#endif

CALL CloseDataFile()

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

CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)

! Check if the provided file is a DSMC state file.
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=File_Type)
IF(TRIM(File_Type).NE.'DSMCState') THEN
  SWRITE(*,*) 'ERROR: The given file type is: ', TRIM(File_Type)
  CALL abort(__STAMP__,'ERROR: Given file for the macroscopic restart is not of the type "DSMCState", please check the input file!')
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

CALL CloseDataFile()

iVar = 1
DO iSpec = 1, nSpecies
  DO iElem = 1, nElems
    MacroRestartValues(iElem,iSpec,:) = ElemData_HDF5(iVar:iVar-1+DSMC_NVARS,iElem)
  END DO
  iVar = iVar + DSMC_NVARS
END DO

! Insert simulation particles based on the macroscopic values
CALL MacroRestart_InsertParticles()

DEALLOCATE(ElemData_HDF5)

END SUBROUTINE MacroscopicRestart


SUBROUTINE RestartAdaptiveBCSampling()
!===================================================================================================================================
!> 1) If a restart is performed,
!>  1a) Check if AdaptiveInfo exists in state, read it in and write to AdaptBCMacroValues
!>  1b) If TruncateRunningAverage: read-in of AdaptiveRunningAverage to continue the sample
!>  1c) SurfaceFlux, Type=4: read-in of AdaptBCPartNumOut to avoid mass flux jumps
!> 2) Adaptive Type = 4: Read-in of the number of particles leaving the domain through the BC (required for the calculation of the massflow)
!> 3) Fall-back: Initialize the macroscopic values from the macroscopic values or surface flux parameter input values (if no values have been read-in)
!> 4) Approximation of particles leaving the domain, assuming zero bulk velocity, using the macrorestart values or init sampling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_INPUT              ,ONLY: ReadArray, ReadAttribute, DatasetExists, GetDataSize
USE MOD_Particle_Sampling_Adapt ,ONLY: AdaptiveBCSampling
USE MOD_Particle_Sampling_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst, Pi
USE MOD_TimeDisc_Vars           ,ONLY: ManualTimeStep, RKdtFrac
USE MOD_Mesh_Vars               ,ONLY: offsetElem, nElems, SideToElem
USE MOD_Particle_Vars           ,ONLY: Species, nSpecies, VarTimeStep
USE MOD_Particle_Surfaces_Vars  ,ONLY: BCdata_auxSF, SurfFluxSideSize, SurfMeshSubSideData
USE MOD_Restart_Vars            ,ONLY: RestartFile, DoMacroscopicRestart, MacroRestartValues, MacroRestartFileName
USE MOD_SurfaceModel_Vars       ,ONLY: nPorousBC
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                           :: AdaptiveDataExists, RunningAverageExists, UseAdaptiveType4, AdaptBCPartNumOutExists
REAL                              :: TimeStepOverWeight, v_thermal, dtVar, WeightingFactor(1:nSpecies)
REAL,ALLOCATABLE                  :: ElemData_HDF5(:,:), ElemData2_HDF5(:,:,:,:)
INTEGER                           :: iElem, iSpec, iSF, iSide, ElemID, SampleElemID, nVar, GlobalElemID, currentBC
INTEGER                           :: jSample, iSample, BCSideID, nElemReadin, nVarTotal, iVar, nVarArrayStart, nVarArrayEnd
INTEGER                           :: SampIterEnd, nSurfacefluxBCs
INTEGER,ALLOCATABLE               :: GlobalElemIndex(:)
!===================================================================================================================================

AdaptiveDataExists = .FALSE.
RunningAverageExists = .FALSE.
UseAdaptiveType4 = .FALSE.
AdaptBCPartNumOutExists = .FALSE.

DO iSpec=1,nSpecies
  DO iSF=1,Species(iSpec)%nSurfacefluxBCs
    IF(Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) UseAdaptiveType4 = .TRUE.
  END DO
END DO

! 1) Check if AdaptiveInfo exists in state, read it in and write to AdaptBCMacroValues
CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
! read local ParticleInfo from HDF5
CALL DatasetExists(File_ID,'AdaptiveInfo',AdaptiveDataExists)
IF(AdaptiveDataExists)THEN
  CALL GetDataSize(File_ID,'AdaptiveInfo',nDims,HSize)
  nVarTotal=INT(HSize(1),4)
  DEALLOCATE(HSize)
  ALLOCATE(ElemData_HDF5(1:nVarTotal,1:nElems))
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        offsetElem => INT(offsetElem,IK),&
        nElems     => INT(nElems,IK)    ,&
        nVarTotal  => INT(nVarTotal,IK)    )
    CALL ReadArray('AdaptiveInfo',2,(/nVarTotal, nElems/),offsetElem,2,RealArray=ElemData_HDF5(:,:))
  END ASSOCIATE
  nVar = 7
  iVar = 1
  DO iSpec = 1,nSpecies
    DO SampleElemID = 1,AdaptBCSampleElemNum
      ElemID = AdaptBCMapSampleToElem(SampleElemID)
      AdaptBCMacroVal(1:7,SampleElemID,iSpec)   = ElemData_HDF5(iVar:iVar-1+nVar,ElemID)
    END DO
    iVar = iVar + nVar
  END DO
  SDEALLOCATE(ElemData_HDF5)
  LBWRITE(*,*) '| Macroscopic values successfully read-in from restart file.'
END IF
! Read-in the running average values from the state file
IF(AdaptBCTruncAverage) THEN
  ! Avoid deleting the sampling iteration after a restart during a later load balacing step
  IF(.NOT.PerformLoadBalance) AdaptBCSampIterReadIn = 0
  CALL DatasetExists(File_ID,'AdaptiveRunningAverage',RunningAverageExists)
  IF(RunningAverageExists)THEN
    ! Read-in the number of sampling iterations from the restart file (might differ from the current number)
    IF(.NOT.PerformLoadBalance) CALL ReadAttribute(File_ID,'AdaptBCSampIter',1,IntScalar=AdaptBCSampIterReadIn)
    ! Get the data size of the read-in array
    CALL GetDataSize(File_ID,'AdaptiveRunningAverage',nDims,HSize)
    nVar=INT(HSize(2),4)
    nElemReadin = INT(HSize(3),4)
    DEALLOCATE(HSize)
    ! Skip the read-in if the array size does not correspond to the current adaptive BC configuration (e.g. a new adaptive BC was added)
    IF(AdaptBCSampleElemNumGlobal.NE.nElemReadin) THEN
      CALL abort(__STAMP__,&
        'TruncateRunningAverage: Number of read-in elements does not correspond to current number of sample elements!')
    END IF
    ! Treatment of different number of iterations between read-in and parameter input
    IF(AdaptBCSampIter.EQ.nVar) THEN
      nVarArrayStart = 1
      nVarArrayEnd = nVar
      SampIterEnd = AdaptBCSampIter
      IF(AdaptBCSampIterReadIn.LT.nVar.AND..NOT.PerformLoadBalance) THEN
        LBWRITE(*,*) '| TruncateRunningAverage: Array not filled in previous simulation run. Continuing at: ', AdaptBCSampIterReadIn + 1
      END IF
    ELSE IF(AdaptBCSampIter.GT.nVar) THEN
      nVarArrayStart = 1
      nVarArrayEnd = nVar
      SampIterEnd = nVar
      LBWRITE(*,*) '| TruncateRunningAverage: Smaller number of sampling iterations in restart file. Continuing at: ', AdaptBCSampIterReadIn + 1
    ELSE
      nVarArrayStart = nVar - AdaptBCSampIter + 1
      nVarArrayEnd = nVar
      SampIterEnd = AdaptBCSampIter
      AdaptBCSampIterReadIn = AdaptBCSampIter
      LBWRITE(*,*) '| TruncateRunningAverage: Greater number of sampling iterations in restart file. Using the last ', AdaptBCSampIterReadIn, ' sample iterations.'
    END IF
    ALLOCATE(ElemData2_HDF5(1:8,1:nVar,1:nElemReadin,1:nSpecies))
    ALLOCATE(GlobalElemIndex(1:nElemReadin))
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nSpecies    => INT(nSpecies,IK) ,&
          nElemReadin => INT(nElemReadin,IK)    ,&
          nVar        => INT(nVar,IK)    )
      CALL ReadArray('AdaptiveRunningAverage',4,(/8_IK, nVar, nElemReadin, nSpecies/),0_IK,3,RealArray=ElemData2_HDF5(:,:,:,:))
      CALL ReadArray('AdaptiveRunningAverageIndex',1,(/nElemReadin/),0_IK,1,IntegerArray_i4=GlobalElemIndex(:))
    END ASSOCIATE
    ! Map the read-in values to the sampling array (GlobalElemID -> LocalElemID -> SampleElemID)
    IF(AdaptBCSampleElemNum.GT.0) THEN
      DO iElem = 1,nElemReadin
        GlobalElemID = GlobalElemIndex(iElem)
        ! Skip elements outside my local region
        IF((GlobalElemID.LT.1+offsetElem).OR.(GlobalElemID.GT.nElems+offsetElem)) CYCLE
        ! Get the sample element ID
        SampleElemID = AdaptBCMapElemToSample(GlobalElemID-offsetElem)
        IF(SampleElemID.GT.0) AdaptBCAverage(1:8,1:SampIterEnd,SampleElemID,1:nSpecies) = ElemData2_HDF5(1:8,nVarArrayStart:nVarArrayEnd,iElem,1:nSpecies)
      END DO
      ! Scaling of the weighted particle number in case of a macroscopic restart with a particle weighting change
      IF(DoMacroscopicRestart.AND..NOT.PerformLoadBalance) THEN
        CALL ReadAttribute(File_ID,'AdaptBCWeightingFactor',nSpecies,RealArray=WeightingFactor)
        DO iSpec = 1, nSpecies
          IF(WeightingFactor(iSpec).NE.Species(iSpec)%MacroParticleFactor) THEN
            AdaptBCAverage(7:8,1:SampIterEnd,:,iSpec) = AdaptBCAverage(7:8,1:SampIterEnd,:,iSpec) * WeightingFactor(iSpec) &
                                                                                              / Species(iSpec)%MacroParticleFactor
          END IF
        END DO
        LBWRITE(*,*) '| TruncateRunningAverage: Sample successfully initiliazed from restart file and scaled due to MacroscopicRestart.'
      ELSE
        LBWRITE(*,*) '| TruncateRunningAverage: Sample successfully initiliazed from restart file.'
      END IF
    END IF
    IF(.NOT.AdaptiveDataExists) THEN
      ! Calculate the macro values initially from the sample for the first iteration
      CALL AdaptiveBCSampling(initTruncAverage_opt=.TRUE.)
      ! Avoid overwriting it with the intial sampling from the current particle state
      AdaptiveDataExists = .TRUE.
      LBWRITE(*,*) '| TruncateRunningAverage: AdaptiveInfo not found in state file. Macroscopic values calculated from sample.'
    END IF
    SDEALLOCATE(ElemData2_HDF5)
    SDEALLOCATE(GlobalElemIndex)
  ELSE
    LBWRITE(*,*) '| TruncateRunningAverage: No running average values found. Values initiliazed with zeros.'
  END IF
END IF
CALL CloseDataFile()

! 2) Adaptive Type = 4: Read-in of the number of particles leaving the domain through the BC (required for the calculation of the massflow)
IF(UseAdaptiveType4) THEN
  IF(PerformLoadBalance) THEN
    ! Array is not deallocated during loadbalance
    AdaptBCPartNumOutExists = .TRUE.
  ELSE IF(DoMacroscopicRestart) THEN
    ! Reset of the number due to a potentially new weighting factor
    AdaptBCPartNumOutExists = .FALSE.
  ELSE
    ! Read-in array during restart only with the root as it is distributed onto all procs later
    IF(MPIRoot)THEN
      CALL OpenDataFile(RestartFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
      CALL DatasetExists(File_ID,'AdaptBCPartNumOut',AdaptBCPartNumOutExists)
      IF(AdaptBCPartNumOutExists) THEN
        ! Associate construct for integer KIND=8 possibility
        ASSOCIATE (&
              nSpecies     => INT(nSpecies,IK) ,&
              nSurfFluxBCs => INT(nSurfacefluxBCs,IK)  )
          CALL ReadArray('AdaptBCPartNumOut',2,(/nSpecies,nSurfFluxBCs/),0_IK,1,IntegerArray_i4=AdaptBCPartNumOut(:,:))
        END ASSOCIATE
        LBWRITE(*,*) '| Surface Flux, Type=4: Number of particles leaving the domain successfully read-in from restart file.'
      END IF
      CALL CloseDataFile()
    END IF
#if USE_MPI
    CALL MPI_BCAST(AdaptBCPartNumOutExists,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,IERROR)
#endif /*USE_MPI*/
  END IF
END IF

! 3) Fall-back: Initialize the macroscopic values from the macroscopic values or surface flux parameter input values (if no values have been read-in)
IF(.NOT.AdaptiveDataExists) THEN
  IF (DoMacroscopicRestart) THEN
    DO SampleElemID = 1,AdaptBCSampleElemNum
      ElemID = AdaptBCMapSampleToElem(SampleElemID)
      AdaptBCMacroVal(DSMC_VELOX,SampleElemID,1:nSpecies) = MacroRestartValues(ElemID,1:nSpecies,DSMC_VELOX)
      AdaptBCMacroVal(DSMC_VELOY,SampleElemID,1:nSpecies) = MacroRestartValues(ElemID,1:nSpecies,DSMC_VELOY)
      AdaptBCMacroVal(DSMC_VELOZ,SampleElemID,1:nSpecies) = MacroRestartValues(ElemID,1:nSpecies,DSMC_VELOZ)
      AdaptBCMacroVal(4,SampleElemID,1:nSpecies)          = MacroRestartValues(ElemID,1:nSpecies,DSMC_NUMDENS)
    END DO
    LBWRITE(*,*) '| Marcroscopic values have been initialized from: ', TRIM(MacroRestartFileName)
    IF(nPorousBC.GT.0) THEN
      CALL abort(__STAMP__,&
        'Macroscopic restart with porous BC and without state file including adaptive BC info not implemented!')
    END IF
  ELSE
    IF(.NOT.PerformLoadBalance) THEN
      CALL AdaptiveBCSampling(initSampling_opt=.TRUE.)
      LBWRITE(*,*) '| Sampling of inserted particles has been performed for an initial distribution.'
    END IF
  END IF
END IF

! 4) Adaptive Type = 4: Approximation of particles leaving the domain, using the values from AdaptBCMacroVal for velocity and number density
IF(UseAdaptiveType4.AND..NOT.AdaptBCPartNumOutExists) THEN
  DO iSpec=1,nSpecies
    ! Species-specific time step
    IF(VarTimeStep%UseSpeciesSpecific) THEN
      dtVar = ManualTimeStep * RKdtFrac * Species(iSpec)%TimeStepFactor
    ELSE
      dtVar = ManualTimeStep * RKdtFrac
    END IF
    DO iSF=1,Species(iSpec)%nSurfacefluxBCs
      currentBC = Species(iSpec)%Surfaceflux(iSF)%BC
      ! Skip processors without a surface flux
      IF (BCdata_auxSF(currentBC)%SideNumber.EQ.0) CYCLE
      ! Skip other regular surface flux and other types
      IF(.NOT.Species(iSpec)%Surfaceflux(iSF)%AdaptiveType.EQ.4) CYCLE
      ! Calculate the velocity for the particles leaving the domain with the thermal velocity assuming a zero bulk velocity
      v_thermal = SQRT(2.*BoltzmannConst*Species(iSpec)%Surfaceflux(iSF)%MWTemperatureIC/Species(iSpec)%MassIC) / (2.0*SQRT(PI))
      TimeStepOverWeight = dtVar / Species(iSpec)%MacroParticleFactor
      ! Loop over sides on the surface flux
      DO iSide=1,BCdata_auxSF(currentBC)%SideNumber
        BCSideID=BCdata_auxSF(currentBC)%SideList(iSide)
        ElemID = SideToElem(S2E_ELEM_ID,BCSideID)
        SampleElemID = AdaptBCMapElemToSample(ElemID)
        IF(SampleElemID.GT.0) THEN
          DO jSample=1,SurfFluxSideSize(2); DO iSample=1,SurfFluxSideSize(1)
            AdaptBCPartNumOut(iSpec,iSF) = AdaptBCPartNumOut(iSpec,iSF) + INT(AdaptBCMacroVal(4,SampleElemID,iSpec) &
              * TimeStepOverWeight * SurfMeshSubSideData(iSample,jSample,BCSideID)%area * v_thermal)
          END DO; END DO
        END IF  ! SampleElemID.GT.0
      END DO    ! iSide=1,BCdata_auxSF(currentBC)%SideNumber
    END DO      ! iSF=1,Species(iSpec)%nSurfacefluxBCs
  END DO        ! iSpec=1,nSpecies
END IF

END SUBROUTINE RestartAdaptiveBCSampling


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
