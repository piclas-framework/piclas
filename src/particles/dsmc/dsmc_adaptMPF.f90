!==================================================================================================================================
! Copyright (c) 2024 Simone Lauterbach
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

MODULE MOD_CellLocalWeighting
!===================================================================================================================================
!> Routines for the adaption of the cell-local particle weights and mapping of cell-local values to the nodes
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: PerformCellLocalWeighting
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadinForCellLocalWeighting()
!===================================================================================================================================
!> Read-in of DSMCState file through the macroscopic restart file name
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_io_hdf5
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_DSMC_Vars               ,ONLY: AdaptMPFInfo_Shared
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems
USE MOD_Particle_Mesh_Vars      ,ONLY: nComputeNodeElems
USE MOD_HDF5_Input              ,ONLY: OpenDataFile, CloseDataFile, ReadArray, ReadAttribute, GetDataProps
USE MOD_HDF5_Input              ,ONLY: GetDataSize, nDims, HSize, File_ID
USE MOD_Restart_Vars            ,ONLY: MacroRestartFileName
USE MOD_StringTools             ,ONLY: STRICMP
#if USE_MPI
USE MOD_DSMC_Vars               ,ONLY: AdaptMPFInfo_Shared_Win
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: nVar_HDF5, N_HDF5, iVar
INTEGER                             :: nVar_TotalPartNum, nVar_TotalDens, nVar_DSMC, nVar_BGK, nVar_FP, nVar_AdaptMPF
INTEGER                             :: offSetLocal, nVar_Ratio_FP, nVar_Ratio_BGK
INTEGER                             :: iElem, ReadInElems, iCNElem, firstElem, lastElem
REAL, ALLOCATABLE                   :: ElemData_HDF5(:,:)
CHARACTER(LEN=255),ALLOCATABLE      :: VarNames_tmp(:)
!===================================================================================================================================
nVar_TotalPartNum = 0; nVar_TotalDens = 0; nVar_Ratio_FP = 0; nVar_Ratio_BGK = 0; nVar_DSMC = 0; nVar_BGK = 0; nVar_AdaptMPF = 0

! Open DSMC state file
CALL OpenDataFile(MacroRestartFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
CALL GetDataProps('ElemData',nVar_HDF5,N_HDF5,nGlobalElems)

IF(nVar_HDF5.LE.0) THEN
  SWRITE(*,*) 'ERROR: Something is wrong with the MacroscopicRestart file:', TRIM(MacroRestartFileName)
  CALL abort(__STAMP__, 'ERROR: Number of variables in the ElemData array appears to be zero!')
END IF

! Get the variable names from the DSMC state and find the position of required quality factors
ALLOCATE(VarNames_tmp(1:nVar_HDF5))
CALL ReadAttribute(File_ID,'VarNamesAdd',nVar_HDF5,StrArray=VarNames_tmp(1:nVar_HDF5))

DO iVar=1,nVar_HDF5
  IF (STRICMP(VarNames_tmp(iVar),"Total_SimPartNum")) nVar_TotalPartNum = iVar
  IF (STRICMP(VarNames_tmp(iVar),"Total_NumberDensity")) nVar_TotalDens = iVar
  IF (STRICMP(VarNames_tmp(iVar),"DSMC_MCS_over_MFP")) nVar_DSMC = iVar
  IF (STRICMP(VarNames_tmp(iVar),"BGK_DSMC_Ratio")) nVar_Ratio_BGK = iVar
  IF (STRICMP(VarNames_tmp(iVar),"FP_DSMC_Ratio")) nVar_Ratio_FP = iVar
  IF (STRICMP(VarNames_tmp(iVar),"BGK_MaxRelaxationFactor")) nVar_BGK = iVar
  IF (STRICMP(VarNames_tmp(iVar),"FP_MaxRelaxationFactor")) nVar_FP = iVar
  IF (STRICMP(VarNames_tmp(iVar),"WeightingFactorCell")) nVar_AdaptMPF = iVar
END DO

IF(Symmetry%Axisymmetric) THEN
  IF(nVar_AdaptMPF.EQ.0) CALL abort(__STAMP__, 'ERROR: Restart of an axisymmetric simulation with the weighting type cell_local'//&
    'from a DSMCState without WeightingFactorCell is not supported!')
END IF

#if USE_MPI
firstElem = INT(REAL(myComputeNodeRank)*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
lastElem = INT(REAL(myComputeNodeRank+1)*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))
offsetLocal = GetGlobalElemID(firstElem)-1
ReadInElems = lastElem - firstElem +1
#else
firstElem = 1
lastElem = nGlobalElems
offSetLocal = 0
ReadInElems = nGlobalElems
#endif

ALLOCATE(ElemData_HDF5(1:nVar_HDF5,1:ReadInElems))
! Associate construct for integer KIND=8 possibility
ASSOCIATE (nVar_HDF5     => INT(nVar_HDF5,IK) ,&
            offSetLocal   => INT(offSetLocal,IK) ,&
            ReadInElems   => INT(ReadInElems,IK))
  CALL ReadArray('ElemData',2,(/nVar_HDF5,ReadInElems/),offSetLocal,2,RealArray=ElemData_HDF5(:,:))
END ASSOCIATE

#if USE_MPI
CALL Allocate_Shared((/7,nComputeNodeElems/),AdaptMPFInfo_Shared_Win,AdaptMPFInfo_Shared)
CALL MPI_WIN_LOCK_ALL(0,AdaptMPFInfo_Shared_Win,iError)
#else
ALLOCATE(AdaptMPFInfo_Shared(7,nComputeNodeElems))
#endif

! Read in of the quality factors from a previous simulation
DO iCNElem=firstElem, lastElem
  iElem = iCNElem - firstElem +1
  AdaptMPFInfo_Shared(1,iCNElem) = ElemData_HDF5(nVar_TotalPartNum,iElem)
  AdaptMPFInfo_Shared(2,iCNElem) = ElemData_HDF5(nVar_TotalDens,iElem)
  IF (nVar_DSMC.NE.0) THEN
    AdaptMPFInfo_Shared(3,iCNElem) = ElemData_HDF5(nVar_DSMC,iElem)
    IF (ISNAN(AdaptMPFInfo_Shared(3,iCNElem))) THEN
      AdaptMPFInfo_Shared(3,iCNElem) = 0.0
    END IF
  ELSE
    AdaptMPFInfo_Shared(3,iCNElem) = 0.
  END IF
  ! BGK-DSMC Ratio
  IF (nVar_BGK.NE.0) THEN
    AdaptMPFInfo_Shared(4,iCNElem) = ElemData_HDF5(nVar_BGK,iElem)
  ELSE
    AdaptMPFInfo_Shared(4,iCNElem) = 0.
  END IF
  ! FP-DSMC Ratio
  IF (nVar_FP.NE.0) THEN
    AdaptMPFInfo_Shared(4,iCNElem) = ElemData_HDF5(nVar_FP,iElem)
  ELSE
    AdaptMPFInfo_Shared(4,iCNElem) = 0.
  END IF
  ! BGK quality factors
  IF (nVar_Ratio_BGK.NE.0) THEN
    AdaptMPFInfo_Shared(5,iCNElem) = REAL(NINT(ElemData_HDF5(nVar_Ratio_BGK,iElem)))
  ELSE IF (nVar_BGK.NE.0) THEN
    AdaptMPFInfo_Shared(5,iCNElem) = 1.
  ! FP quality factors
  ELSE IF (nVar_Ratio_FP.NE.0) THEN
    AdaptMPFInfo_Shared(5,iCNElem) = REAL(2*NINT(ElemData_HDF5(nVar_Ratio_FP,iElem)))
  ELSE IF (nVar_FP.NE.0) THEN
    AdaptMPFInfo_Shared(5,iCNElem) = 2.
  ELSE
    AdaptMPFInfo_Shared(5,iCNElem) = 0.
  END IF
  AdaptMPFInfo_Shared(6,iCNElem) = 0.
  ! Adapted MPF from a previous simulation
  IF (nVar_AdaptMPF.NE.0) THEN
    AdaptMPFInfo_Shared(7,iCNElem) = ElemData_HDF5(nVar_AdaptMPF,iElem)
  ELSE
    AdaptMPFInfo_Shared(7,iCNElem) = 0.
  END IF
END DO

#if USE_MPI
CALL BARRIER_AND_SYNC(AdaptMPFInfo_Shared_Win,MPI_COMM_SHARED)
#endif
CALL CloseDataFile()

SDEALLOCATE(ElemData_HDF5)

END SUBROUTINE ReadinForCellLocalWeighting


SUBROUTINE PerformCellLocalWeighting()
!===================================================================================================================================
!> Routine for the automatic adaption of the particles weights in each simulation cell based on the read-in of the particle
!> number density and simulation particle number from a previous simulation
!> 1) Read-in of DSMCState file through the macroscopic restart file name
!> 2) Initialize the node mapping (requires NodeToElemMapping and thus FindNeighbourElems = T)
!> 3) Determine cell-local weighting factors
!> 3A) SkipAdaption = F
!>    a) Determine the scaling parameters (such as minimum required particle number) based on the method (DSMC/BGK/FP)
!>    b) Calculate the new weighting factor in OptimalMPF_Shared
!> 3B) SkipAdaption = T
!>    a) Read-in weights from the State file (ElemLocalWeight) and store in OptimalMPF_Shared
!> 4) Map the weights from the element to the nodes and average at each node
!> 5) (Optional) Average the MPF distribution by the neighbour values
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_PreProc
USE MOD_io_hdf5
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars
USE MOD_DSMC_Symmetry
USE MOD_DSMC_Vars
USE MOD_Restart_Vars            ,ONLY: DoCatalyticRestart, RestartFile, DoMacroscopicRestart
USE MOD_HDF5_Input              ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,ReadArray
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems, offSetElem, SideToElem, nBCSides, BC
USE MOD_Globals_Vars            ,ONLY: Pi
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Symmetry_Vars           ,ONLY: Symmetry
USE MOD_Particle_Vars           ,ONLY: Species
USE MOD_Particle_Mesh_Vars      ,ONLY: ElemVolume_Shared, ElemMidPoint_Shared, nComputeNodeElems, GEO
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_part_tools              ,ONLY: CalcVarWeightMPF
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER                           :: iCNElem, firstElem, lastElem, offSetLocal, ReadInElems, iRefine, iGlobalElem
INTEGER                           :: iBound, iSide, ElemID, CNElemId
REAL                              :: MinPartNum, MaxPartNum
LOGICAL, ALLOCATABLE              :: RefineCatElem(:)
LOGICAL                           :: MPFExists
REAL, ALLOCATABLE                 :: MPFData_HDF5(:)
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' PERFORM CELL-LOCAL WEIGHTING...'

! 1) If adaption is enabled, read-in cell-local values from the given DSMCState
IF (.NOT.CellLocalWeight%SkipAdaption) CALL ReadinForCellLocalWeighting()

! 2) Initialize node mapping
CALL InitNodeMapping()

#if USE_MPI
firstElem = INT(REAL(myComputeNodeRank)*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))+1
lastElem = INT(REAL(myComputeNodeRank+1)*REAL(nComputeNodeElems)/REAL(nComputeNodeProcessors))
offsetLocal = GetGlobalElemID(firstElem)-1
ReadInElems = lastElem - firstElem +1
CALL Allocate_Shared((/nComputeNodeElems/),OptimalMPF_Shared_Win,OptimalMPF_Shared)
CALL MPI_WIN_LOCK_ALL(0,OptimalMPF_Shared_Win,iError)
#else
firstElem = 1
lastElem = nGlobalElems
offSetLocal = 0
ReadInElems = nGlobalElems
ALLOCATE(OptimalMPF_Shared(nGlobalElems))
#endif

! 3) Determine cell-local weighting factors
IF (.NOT.CellLocalWeight%SkipAdaption) THEN
  ! Further refinement of elements close to a catalytic surface
  IF (DoCatalyticRestart) THEN
    ALLOCATE(RefineCatElem(1:nComputeNodeElems))
    RefineCatElem = .FALSE.
    DO iSide=1,nBCSides
      iBound = PartBound%MapToPartBC(BC(iSide))
      IF (PartBound%TargetBoundCond(iBound).EQ.PartBound%ReflectiveBC) THEN
        ElemID = SideToElem(S2E_ELEM_ID,iSide)
        CNElemId = GetCNElemID(ElemID + offSetElem)
        RefineCatElem(CNElemID) = .TRUE.
      END IF
    END DO
  END IF

  ! Determine the MPF based on the particle number from the reference simulation
  DO iCNElem = firstElem, lastElem
    ! Determine the reference MPF
    IF (AdaptMPFInfo_Shared(7,iCNElem).NE.0.) THEN
      AdaptMPFInfo_Shared(6,iCNElem) = AdaptMPFInfo_Shared(7,iCNElem)
    ELSE
      AdaptMPFInfo_Shared(6,iCNElem) = Species(1)%MacroParticleFactor
    END IF

    IF (AdaptMPFInfo_Shared(5,iCNElem).EQ.1.) THEN
      ! Adaption based on the BGK quality factor
      IF (AdaptMPFInfo_Shared(4,iCNElem).GT.CellLocalWeight%QualityFactor) THEN
        OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)*(CellLocalWeight%QualityFactor/AdaptMPFInfo_Shared(4,iCNElem))
      ! Adaption based on the particle number per simulation cell
      ELSE ! BGKQualityFactors
        ! Further refinement for the elements close to the symmetry axis or catalytic surfaces
        IF (DoCatalyticRestart) THEN
          IF (RefineCatElem(iCNElem)) THEN
            MinPartNum = CellLocalWeight%Cat_MinPartNum
          ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
            MinPartNum = CellLocalWeight%SymAxis_MinPartNum
          ELSE
            MinPartNum = CellLocalWeight%MinPartNum
          END IF
        ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
          MinPartNum = CellLocalWeight%SymAxis_MinPartNum
        ELSE
          MinPartNum = CellLocalWeight%MinPartNum
        END IF
        IF (MinPartNum.GT.CellLocalWeight%MaxPartNum) THEN
          MaxPartNum = MinPartNum*10.
        ELSE
          MaxPartNum = CellLocalWeight%MaxPartNum
        END IF
        IF(AdaptMPFInfo_Shared(1,iCNElem).LT.MinPartNum) THEN
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/MinPartNum
        ELSE IF(CellLocalWeight%IncludeMaxPartNum) THEN ! Check if the particle number should be decreased
          IF(AdaptMPFInfo_Shared(1,iCNElem).GT.MaxPartNum) THEN
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/CellLocalWeight%MaxPartNum
          ELSE ! Further refinement BGK
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem) * CellLocalWeight%BGKFactor
          END IF
        ELSE ! Further refinement BGK
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem) * CellLocalWeight%BGKFactor
        END IF
      END IF ! BGKQualityFactors

    ELSE IF (AdaptMPFInfo_Shared(5,iCNElem).EQ.2.) THEN
      ! Adaption based on the FP quality factor
      IF (AdaptMPFInfo_Shared(4,iCNElem).GT.CellLocalWeight%QualityFactor) THEN
        OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)*(CellLocalWeight%QualityFactor/AdaptMPFInfo_Shared(4,iCNElem))
      ! Adaption based on the particle number per simulation cell
      ELSE ! FPQualityFactors
        ! Further refinement for the elements close to the symmetry axis or catalytic surfaces
        IF (DoCatalyticRestart) THEN
          IF (RefineCatElem(iCNElem)) THEN
            MinPartNum = CellLocalWeight%Cat_MinPartNum
          ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
            MinPartNum = CellLocalWeight%SymAxis_MinPartNum
          ELSE
            MinPartNum = CellLocalWeight%MinPartNum
          END IF
        ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
          MinPartNum = CellLocalWeight%SymAxis_MinPartNum
        ELSE
          MinPartNum = CellLocalWeight%MinPartNum
        END IF
        IF (MinPartNum.GT.CellLocalWeight%MaxPartNum) THEN
          MaxPartNum = MinPartNum*10.
        ELSE
          MaxPartNum = CellLocalWeight%MaxPartNum
        END IF
        IF(AdaptMPFInfo_Shared(1,iCNElem).LT.MinPartNum) THEN
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/MinPartNum
        ELSE IF(CellLocalWeight%IncludeMaxPartNum) THEN ! Check if the particle number should be decreased
          IF(AdaptMPFInfo_Shared(1,iCNElem).GT.MaxPartNum) THEN
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/CellLocalWeight%MaxPartNum
          ELSE ! Further refinement FP
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem) * CellLocalWeight%FPFactor
          END IF
        ELSE ! Further refinement FP
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem) * CellLocalWeight%FPFactor
        END IF
      END IF ! FPQualityFactors

    ELSE

      ! Adaption based on the DSMC quality factor
      IF (AdaptMPFInfo_Shared(3,iCNElem).GT.CellLocalWeight%QualityFactor) THEN
        IF (Symmetry%Order.EQ.2) THEN
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)*(CellLocalWeight%QualityFactor/AdaptMPFInfo_Shared(3,iCNElem))**2
        ELSE
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)*(CellLocalWeight%QualityFactor/AdaptMPFInfo_Shared(3,iCNElem))**3
        END IF
      ! Adaption based on the particle number per simulation cell
      ELSE ! DSMCQualityFactors
        ! Further refinement for the elements close to the symmetry axis or catalytic surfaces
        IF (DoCatalyticRestart) THEN
          IF (RefineCatElem(iCNElem)) THEN
            MinPartNum = CellLocalWeight%Cat_MinPartNum
          ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
            MinPartNum = CellLocalWeight%SymAxis_MinPartNum
          ELSE
            MinPartNum = CellLocalWeight%MinPartNum
          END IF
        ELSE IF ((Symmetry%Axisymmetric).AND.(ElemMidPoint_Shared(2,iCNElem).LE.(GEO%ymaxglob*0.05))) THEN
          MinPartNum = CellLocalWeight%SymAxis_MinPartNum
        ELSE
          MinPartNum = CellLocalWeight%MinPartNum
        END IF
        IF (MinPartNum.GT.CellLocalWeight%MaxPartNum) THEN
          MaxPartNum = MinPartNum*10.
        ELSE
          MaxPartNum = CellLocalWeight%MaxPartNum
        END IF
        IF(AdaptMPFInfo_Shared(1,iCNElem).LT.MinPartNum) THEN
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/MinPartNum
        ELSE IF(CellLocalWeight%IncludeMaxPartNum) THEN ! Check if the particle number should be decreased
          IF(AdaptMPFInfo_Shared(1,iCNElem).GT.MaxPartNum) THEN
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(2,iCNElem)*ElemVolume_Shared(iCNElem)/CellLocalWeight%MaxPartNum
          ELSE
            OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)
          END IF
        ELSE
          OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)
        END IF
      END IF ! DSMCQualityFactors
    END IF !BGK_DSMC_Ratio

    ! If not defined, determine the optimal MPF from the previous simulation
    IF (OptimalMPF_Shared(iCNElem).LE.0.) THEN
      OptimalMPF_Shared(iCNElem) = AdaptMPFInfo_Shared(6,iCNElem)
    END IF
  END DO ! iGlobalElem

ELSE ! Skip Adaption
  IF(PerformLoadBalance.AND..NOT.UseH5IOLoadBalance) THEN
    CALL abort(__STAMP__, 'ERROR: Cell-local weighting requires UseH5IOLoadBalance!')
  END IF
  ! Read-in of the MPF per element from a previous adaption process
  ALLOCATE(MPFData_HDF5(1:nGlobalElems))
  MPFData_HDF5 = 0.
  CALL OpenDataFile(TRIM(RestartFile),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL DatasetExists(File_ID,'ElemLocalWeight',MPFExists)
  IF(MPFExists) THEN
    ASSOCIATE(nGlobalElems    => INT(nGlobalElems,IK))
      CALL ReadArray('ElemLocalWeight',2,(/nGlobalElems, 1_IK/),0_IK,1,RealArray=MPFData_HDF5(1:nGlobalElems))
    END ASSOCIATE
    DO iCNElem=firstElem, lastElem
      iGlobalElem = GetGlobalElemID(iCNElem)
      OptimalMPF_Shared(iCNElem) = MPFData_HDF5(iGlobalElem)
    END DO

    LBWRITE(UNIT_stdOut,*)' | Cell-local weighting: Read-in of particle weight distribution from state file.'

  ELSE IF(.NOT.DoMacroscopicRestart) THEN
    CALL abort(__STAMP__, 'ERROR: Cell-local weight requires a given particle weight distribution or -DoMacroscopicRestart=T!')
  END IF
  CALL CloseDataFile()
END IF

#if USE_MPI
CALL BARRIER_AND_SYNC(OptimalMPF_Shared_Win,MPI_COMM_SHARED)
#endif

! 4) Mapping of the particle MPF to the nodes of an element (and average)
CALL MappingCellLocalWeightToNode()

! 5) (Optional) Average the MPF distribution by the neighbour values
IF (CellLocalWeight%UseMedianFilter) THEN
  IF (CellLocalWeight%SkipAdaption) THEN
    LBWRITE(UNIT_stdOut,*) ' | ApplyMedianFilter is not possible with SkipAdaption, filtering process is skipped.'
  ELSE
    DO iRefine=1, CellLocalWeight%nRefine
      CALL NodeMappingFilterMPF()
    END DO
  END IF
END IF ! UseMedianFilter

CALL FinalizeNodeMapping()

! Disabling the adaption for load-balance steps
CellLocalWeight%SkipAdaption = .TRUE.

LBWRITE(UNIT_stdOut,'(A)') ' PERFORM CELL-LOCAL WEIGHTING DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE PerformCellLocalWeighting


SUBROUTINE InitNodeMapping()
!===================================================================================================================================
!> Mapping of the adapted particle weights from the elements to the nodes and interpolation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Particle_Vars
USE MOD_Particle_Mesh_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_Mesh_Vars              ,ONLY: offsetElem
USE MOD_MPI_Shared             ,ONLY: BARRIER_AND_SYNC
USE MOD_Mesh_Tools             ,ONLY: GetGlobalElemID, GetCNElemID
USE MOD_MPI_Shared_Vars        ,ONLY: nComputeNodeTotalElems, nLeaderGroupProcs, nProcessors_Global
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: iNode
#if USE_MPI
INTEGER                   :: iElem, ElemID
INTEGER                   :: UniqueNodeID
INTEGER                   :: jElem, NonUniqueNodeID
INTEGER                   :: SendNodeCount, GlobalElemRank, iProc
INTEGER                   :: TestElemID, GlobalElemRankOrig, iRank
LOGICAL,ALLOCATABLE       :: NodeMapping(:,:), DoNodeMapping(:), SendNode(:), IsMappedNode(:)
LOGICAL                   :: bordersMyrank
TYPE(MPI_Request)         :: SendRequestNonSym(0:nProcessors_Global-1)      , RecvRequestNonSym(0:nProcessors_Global-1)
INTEGER                   :: nSendUniqueNodesNonSym(0:nProcessors_Global-1) , nRecvUniqueNodesNonSym(0:nProcessors_Global-1)
INTEGER                   :: GlobalRankToNodeSendRank(0:nProcessors_Global-1)
LOGICAL,ALLOCATABLE       :: FlagShapeElemAdapt(:)
#endif
!===================================================================================================================================
IF (.NOT.ALLOCATED(PartWeightAtNode)) THEN
  ALLOCATE(PartWeightAtNode(1:2,1:nUniqueGlobalNodes))
ELSE
  SDEALLOCATE(PartWeightAtNode)
  ALLOCATE(PartWeightAtNode(1:2,1:nUniqueGlobalNodes))
END IF
PartWeightAtNode=0.0

#if USE_MPI
! Initialization
ALLOCATE(FlagShapeElemAdapt(nComputeNodeTotalElems))
FlagShapeElemAdapt = .FALSE.

DO iElem = 1,nComputeNodeTotalElems
  ElemID    = GetGlobalElemID(iElem)
  IF (ElemInfo_Shared(ELEM_HALOFLAG,ElemID).NE.4) FlagShapeElemAdapt(iElem) = .TRUE.
END DO

ALLOCATE(RecvRequestCN(0:nLeaderGroupProcs-1), SendRequestCN(0:nLeaderGroupProcs-1))
ALLOCATE(DoNodeMapping(0:nProcessors_Global-1),SendNode(1:nUniqueGlobalNodes))
DoNodeMapping = .FALSE.
SendNode = .FALSE.
DO iElem = 1,nComputeNodeTotalElems
  IF (FlagShapeElemAdapt(iElem)) THEN
    bordersMyrank = .FALSE.
    ! Loop all local nodes
    TestElemID = GetGlobalElemID(iElem)
    GlobalElemRankOrig = ElemInfo_Shared(ELEM_RANK,TestElemID)

    DO iNode = 1, 8
      NonUniqueNodeID = ElemNodeID_Shared(iNode,iElem)
      UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
      ! Loop 1D array [offset + 1 : offset + NbrOfElems]
      ! (all CN elements that are connected to the local nodes)
      DO jElem = NodeToElemMapping(1,UniqueNodeID) + 1, NodeToElemMapping(1,UniqueNodeID) + NodeToElemMapping(2,UniqueNodeID)
        TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
        GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
        ! check if element for this side is on the current compute-node. Alternative version to the check above
        IF (GlobalElemRank.EQ.myRank) THEN
          bordersMyrank = .TRUE.
          SendNode(UniqueNodeID) = .TRUE.
        END IF
      END DO
      IF (bordersMyrank) DoNodeMapping(GlobalElemRankOrig) = .TRUE.
    END DO
  END IF
END DO

nMapNodes = 0
ALLOCATE(IsMappedNode(1:nUniqueGlobalNodes))
IsMappedNode = .FALSE.
DO iElem =1, nElems
  TestElemID = GetCNElemID(iElem + offsetElem)
  DO iNode = 1, 8
    NonUniqueNodeID = ElemNodeID_Shared(iNode,TestElemID)
    UniqueNodeID = NodeInfo_Shared(NonUniqueNodeID)
    IsMappedNode(UniqueNodeID) = .TRUE.
  END DO
END DO
nMapNodes = COUNT(IsMappedNode)
nMapNodesTotal = nMapNodes
DO iNode=1, nUniqueGlobalNodes
  IF (.NOT.IsMappedNode(iNode).AND.SendNode(iNode)) THEN
    nMapNodesTotal = nMapNodesTotal + 1
  END IF
END DO

ALLOCATE(NodetoGlobalNode(1:nMapNodesTotal))
nMapNodesTotal = 0
DO iNode=1, nUniqueGlobalNodes
  IF (IsMappedNode(iNode)) THEN
    nMapNodesTotal = nMapNodesTotal + 1
    NodetoGlobalNode(nMapNodesTotal) = iNode
  END IF
END DO
DO iNode=1, nUniqueGlobalNodes
  IF (.NOT.IsMappedNode(iNode).AND.SendNode(iNode)) THEN
    nMapNodesTotal = nMapNodesTotal + 1
    NodetoGlobalNode(nMapNodesTotal) = iNode
  END IF
END DO

GlobalRankToNodeSendRank = -1
nNodeSendExchangeProcs = COUNT(DoNodeMapping)
ALLOCATE(NodeSendRankToGlobalRank(1:nNodeSendExchangeProcs))
NodeSendRankToGlobalRank = 0
nNodeSendExchangeProcs = 0
DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (DoNodeMapping(iRank)) THEN
    nNodeSendExchangeProcs = nNodeSendExchangeProcs + 1
    GlobalRankToNodeSendRank(iRank) = nNodeSendExchangeProcs
    NodeSendRankToGlobalRank(nNodeSendExchangeProcs) = iRank
  END IF
END DO
ALLOCATE(NodeMapping(1:nNodeSendExchangeProcs, 1:nUniqueGlobalNodes))
NodeMapping = .FALSE.

DO iNode = 1, nUniqueGlobalNodes
  IF (SendNode(iNode)) THEN
    DO jElem = NodeToElemMapping(1,iNode) + 1, NodeToElemMapping(1,iNode) + NodeToElemMapping(2,iNode)
        TestElemID = GetGlobalElemID(NodeToElemInfo(jElem))
        GlobalElemRank = ElemInfo_Shared(ELEM_RANK,TestElemID)
        ! check if element for this side is on the current compute-node. Alternative version to the check above
        IF (GlobalElemRank.NE.myRank) THEN
          iRank = GlobalRankToNodeSendRank(GlobalElemRank)
          IF (iRank.LT.1) CALL ABORT(__STAMP__,'Found not connected Rank!', IERROR)
          NodeMapping(iRank, iNode) = .TRUE.
        END IF
      END DO
  END IF
END DO

! Get number of send nodes for each proc: Size of each message for each proc
nSendUniqueNodesNonSym        = 0
nRecvUniqueNodesNonSym(myrank) = 0
ALLOCATE(NodeMappingSend(1:nNodeSendExchangeProcs))
DO iProc = 1, nNodeSendExchangeProcs
  NodeMappingSend(iProc)%nSendUniqueNodes = 0
  DO iNode = 1, nUniqueGlobalNodes
    IF (NodeMapping(iProc,iNode)) NodeMappingSend(iProc)%nSendUniqueNodes = NodeMappingSend(iProc)%nSendUniqueNodes + 1
  END DO
  ! local to global array
  nSendUniqueNodesNonSym(NodeSendRankToGlobalRank(iProc)) = NodeMappingSend(iProc)%nSendUniqueNodes
END DO

! Open receive buffer for non-symmetric exchange identification
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_IRECV( nRecvUniqueNodesNonSym(iProc) &
                , 1                             &
                , MPI_INTEGER                   &
                , iProc                         &
                , 1999                          &
                , MPI_COMM_PICLAS                &
                , RecvRequestNonSym(iProc)      &
                , IERROR)
END DO

! Send each proc the number of nodes
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_ISEND( nSendUniqueNodesNonSym(iProc) &
                , 1                             &
                , MPI_INTEGER                   &
                , iProc                         &
                , 1999                          &
                , MPI_COMM_PICLAS                &
                , SendRequestNonSym(iProc)      &
                , IERROR)
END DO

! Finish communication
DO iProc = 0,nProcessors_Global-1
  IF (iProc.EQ.myRank) CYCLE
  CALL MPI_WAIT(RecvRequestNonSym(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  CALL MPI_WAIT(SendRequestNonSym(iProc),MPI_STATUS_IGNORE,IERROR)
  IF(IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

nNodeRecvExchangeProcs = COUNT(nRecvUniqueNodesNonSym.GT.0)
ALLOCATE(NodeMappingRecv(1:nNodeRecvExchangeProcs))
ALLOCATE(NodeRecvRankToGlobalRank(1:nNodeRecvExchangeProcs))
NodeRecvRankToGlobalRank = 0
nNodeRecvExchangeProcs = 0
DO iRank= 0, nProcessors_Global-1
  IF (iRank.EQ.myRank) CYCLE
  IF (nRecvUniqueNodesNonSym(iRank).GT.0) THEN
    nNodeRecvExchangeProcs = nNodeRecvExchangeProcs + 1
    ! Store global rank of iRecvRank
    NodeRecvRankToGlobalRank(nNodeRecvExchangeProcs) = iRank
    ! Store number of nodes of iRecvRank
    NodeMappingRecv(nNodeRecvExchangeProcs)%nRecvUniqueNodes = nRecvUniqueNodesNonSym(iRank)
  END IF
END DO

! Open receive buffer
ALLOCATE(RecvRequest(1:nNodeRecvExchangeProcs))
DO iProc = 1, nNodeRecvExchangeProcs
  ALLOCATE(NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  ALLOCATE(NodeMappingRecv(iProc)%RecvNodeFilterMPF(1:2,1:NodeMappingRecv(iProc)%nRecvUniqueNodes))
  CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID                   &
                , NodeMappingRecv(iProc)%nRecvUniqueNodes                         &
                , MPI_INTEGER                                                 &
                , NodeRecvRankToGlobalRank(iProc)                         &
                , 666                                                         &
                , MPI_COMM_PICLAS                                              &
                , RecvRequest(iProc)                                          &
                , IERROR)
END DO

! Open send buffer
ALLOCATE(SendRequest(1:nNodeSendExchangeProcs))
DO iProc = 1, nNodeSendExchangeProcs
  ALLOCATE(NodeMappingSend(iProc)%SendNodeUniqueGlobalID(1:NodeMappingSend(iProc)%nSendUniqueNodes))
  NodeMappingSend(iProc)%SendNodeUniqueGlobalID=-1
  ALLOCATE(NodeMappingSend(iProc)%SendNodeFilterMPF(1:2,1:NodeMappingSend(iProc)%nSendUniqueNodes))
  NodeMappingSend(iProc)%SendNodeFilterMPF=0.
  SendNodeCount = 0
  DO iNode = 1, nUniqueGlobalNodes
    IF (NodeMapping(iProc,iNode)) THEN
      SendNodeCount = SendNodeCount + 1
      NodeMappingSend(iProc)%SendNodeUniqueGlobalID(SendNodeCount) = iNode
    END IF
  END DO
  CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeUniqueGlobalID                   &
                , NodeMappingSend(iProc)%nSendUniqueNodes                         &
                , MPI_INTEGER                                                 &
                , NodeSendRankToGlobalRank(iProc)                         &
                , 666                                                         &
                , MPI_COMM_PICLAS                                              &
                , SendRequest(iProc)                                          &
                , IERROR)
END DO

! Finish send
DO iProc = 1, nNodeSendExchangeProcs
  CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

! Finish receive
DO iProc = 1, nNodeRecvExchangeProcs
  CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
  IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
END DO

SDEALLOCATE(FlagShapeElemAdapt)
#else
nMapNodes      = nUniqueGlobalNodes
nMapNodesTotal = nMapNodes
ALLOCATE(NodetoGlobalNode(1:nMapNodesTotal))
DO iNode=1, nUniqueGlobalNodes
  NodetoGlobalNode(iNode) = iNode
END DO
#endif /*USE_MPI*/

END SUBROUTINE InitNodeMapping


SUBROUTINE MappingCellLocalWeightToNode()
!===================================================================================================================================
! Mapping of the cell-local weighting factor to the nodes and averaging of the value at the node
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Mesh_Vars
USE MOD_DSMC_Symmetry
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_DSMC_Vars          ,ONLY: OptimalMPF_Shared
#if USE_MPI
USE MOD_MPI_Shared         ,ONLY: BARRIER_AND_SYNC
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iElem, NodeID(1:8), iNode, globalNode
#if USE_MPI
INTEGER                    :: iProc
#endif /*USE_MPI*/
!===================================================================================================================================

! Nullify PartWeightAtNode
DO iNode = 1, nMapNodesTotal
  globalNode = NodetoGlobalNode(iNode)
  PartWeightAtNode(:,globalNode) = 0.0
END DO

! Loop over all elements and map their weighting factor from the element to the nodes
DO iElem =1, nElems
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
  DO iNode = 1, 8
    PartWeightAtNode(1,NodeID(iNode)) = PartWeightAtNode(1,NodeID(iNode)) + OptimalMPF_Shared(GetCNElemID(iElem+offsetElem))
    PartWeightAtNode(2,NodeID(iNode)) = PartWeightAtNode(2,NodeID(iNode)) + 1.
  END DO
END DO

#if USE_MPI
! 1) Receive MPF values

  DO iProc = 1, nNodeRecvExchangeProcs
    ! Open receive buffer
    CALL MPI_IRECV( NodeMappingRecv(iProc)%RecvNodeFilterMPF(1:2,:) &
        , 2*NodeMappingRecv(iProc)%nRecvUniqueNodes               &
        , MPI_DOUBLE_PRECISION                                  &
        , NodeRecvRankToGlobalRank(iProc)                       &
        , 666                                                   &
        , MPI_COMM_PICLAS                                        &
        , RecvRequest(iProc)                                    &
        , IERROR)
  END DO

  ! 2) Send MPF values
  DO iProc = 1, nNodeSendExchangeProcs
    ! Send message (non-blocking)
    DO iNode = 1, NodeMappingSend(iProc)%nSendUniqueNodes
      NodeMappingSend(iProc)%SendNodeFilterMPF(1:2,iNode) = PartWeightAtNode(1:2,NodeMappingSend(iProc)%SendNodeUniqueGlobalID(iNode))
    END DO
    CALL MPI_ISEND( NodeMappingSend(iProc)%SendNodeFilterMPF(1:2,:)     &
        , 2*NodeMappingSend(iProc)%nSendUniqueNodes                   &
        , MPI_DOUBLE_PRECISION                                      &
        , NodeSendRankToGlobalRank(iProc)                           &
        , 666                                                       &
        , MPI_COMM_PICLAS                                            &
        , SendRequest(iProc)                                        &
        , IERROR)
  END DO

  ! Finish communication/
  DO iProc = 1, nNodeSendExchangeProcs
    CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO
  DO iProc = 1, nNodeRecvExchangeProcs
    CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
    IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error', IERROR)
  END DO

  ! 3) Extract messages
  DO iProc = 1, nNodeRecvExchangeProcs
    DO iNode = 1, NodeMappingRecv(iProc)%nRecvUniqueNodes
      ASSOCIATE( NV => PartWeightAtNode(1:2,NodeMappingRecv(iProc)%RecvNodeUniqueGlobalID(iNode)))
        NV = NV +  NodeMappingRecv(iProc)%RecvNodeFilterMPF(1:2,iNode)
      END ASSOCIATE
    END DO
  END DO
#endif /*USE_MPI*/

! Determine the average node value
DO iNode = 1, nMapNodesTotal
  globalNode = NodetoGlobalNode(iNode)
  IF (PartWeightAtNode(2,globalNode).GT.0.) THEN
    PartWeightAtNode(1,globalNode) = PartWeightAtNode(1,globalNode) / PartWeightAtNode(2,globalNode)
  END IF
END DO

END SUBROUTINE MappingCellLocalWeightToNode


SUBROUTINE NodeMappingFilterMPF()
!===================================================================================================================================
! Filter the adapted MPF by multiple mapping from the element to the node and back (corresponds to averaging over the nodes)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,PartWeightAtNode
USE MOD_Mesh_Vars          ,ONLY: nElems, offsetElem
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_DSMC_Vars          ,ONLY: OptimalMPF_Shared
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iElem, NodeID(1:8), iNode, CNElemID
!===================================================================================================================================
DO iElem =1, nElems
  CNElemID = GetCNElemID(iElem+offsetElem)
  ! Set the optimal MPF to zero and recalculate it based on the surrounding node values
  OptimalMPF_Shared(CNElemID) = 0.
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,GetCNElemID(iElem+offsetElem)))
  DO iNode = 1, 8
    OptimalMPF_Shared(CNElemID) = OptimalMPF_Shared(CNElemID) + PartWeightAtNode(1,NodeID(iNode))
  END DO
  OptimalMPF_Shared(CNElemID) = OptimalMPF_Shared(CNElemID)/8.
END DO

CALL MappingCellLocalWeightToNode()

END SUBROUTINE NodeMappingFilterMPF


SUBROUTINE FinalizeNodeMapping()
!===================================================================================================================================
!> Finalize variables not required during the simulation, only PartWeightAtNode is used
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Mesh_Vars
USE MOD_DSMC_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED
USE MOD_DSMC_Vars         ,ONLY: CellLocalWeight
#endif /*USE_MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_SHARED,iERROR)
CALL UNLOCK_AND_FREE(OptimalMPF_Shared_Win)
IF(.NOT.CellLocalWeight%SkipAdaption) CALL UNLOCK_AND_FREE(AdaptMPFInfo_Shared_Win)
#endif /*USE_MPI*/

ADEALLOCATE(OptimalMPF_Shared)
ADEALLOCATE(AdaptMPFInfo_Shared)

#if USE_MPI
SDEALLOCATE(RecvRequestCN)
SDEALLOCATE(SendRequestCN)
SDEALLOCATE(RecvRequest)
SDEALLOCATE(SendRequest)
SDEALLOCATE(NodeMappingSend)
SDEALLOCATE(NodeMappingRecv)
SDEALLOCATE(NodeSendRankToGlobalRank)
SDEALLOCATE(NodeRecvRankToGlobalRank)
#endif /*USE_MPI*/

SDEALLOCATE(NodetoGlobalNode)

END SUBROUTINE FinalizeNodeMapping

END MODULE MOD_CellLocalWeighting
