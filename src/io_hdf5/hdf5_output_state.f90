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

MODULE MOD_HDF5_Output_State
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! Private Part ---------------------------------------------------------------------------------------------------------------------
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteStateToHDF5
  MODULE PROCEDURE WriteStateToHDF5
END INTERFACE

PUBLIC :: WriteStateToHDF5
#if defined(PARTICLES)
PUBLIC :: WriteIMDStateToHDF5
PUBLIC :: WriteElemDataToSeparateContainer
#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteStateToHDF5(MeshFileName,OutputTime,PreviousTime)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars                ,ONLY: U
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems,nGlobalUniqueSides,nUniqueSides,offsetSide
USE MOD_Equation_Vars          ,ONLY: StrVarNames
USE MOD_Restart_Vars           ,ONLY: RestartFile,DoInitialAutoRestart
#ifdef PARTICLES
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_PICDepo_Vars           ,ONLY: OutputSource,PartSource
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars ,ONLY: DoBoundaryParticleOutputHDF5, PartBound
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_Particle_Tracking_Vars ,ONLY: CountNbrOfLostParts,TotalNbrOfMissingParticlesSum,NbrOfNewLostParticlesTotal
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Analyze_Vars  ,ONLY: nSpecAnalyze
USE MOD_Particle_Analyze_Tools ,ONLY: CalcNumPartsOfSpec
USE MOD_HDF5_Output_Particles  ,ONLY: WriteNodeSourceExtToHDF5,WriteClonesToHDF5,WriteVibProbInfoToHDF5,WriteAdaptiveWallTempToHDF5
USE MOD_HDF5_Output_Particles  ,ONLY: WriteAdaptiveInfoToHDF5,WriteParticleToHDF5,WriteBoundaryParticleToHDF5
USE MOD_HDF5_Output_Particles  ,ONLY: WriteLostParticlesToHDF5,WriteEmissionVariablesToHDF5
USE MOD_Particle_Vars          ,ONLY: CalcBulkElectronTemp,BulkElectronTemp
#endif /*PARTICLES*/
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: E,Phi
#endif /*PP_POIS*/
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: nGP_face,iLocSides,UseFPC,FPC,UseEPC,EPC
#if PP_nVar==1
USE MOD_Equation_Vars          ,ONLY: E,Et
#elif PP_nVar==3
USE MOD_Equation_Vars          ,ONLY: B
#else
USE MOD_Equation_Vars          ,ONLY: E,B
#endif /*PP_nVar*/
USE MOD_Mesh_Vars              ,ONLY: nSides
USE MOD_Utils                  ,ONLY: QuickSortTwoArrays
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Utils                  ,ONLY: Qsort1DoubleInt1PInt
USE MOD_Mesh_Tools             ,ONLY: LambdaSideToMaster
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: lastInnerSide
USE MOD_MPI_Vars               ,ONLY: OffsetMPISides_rec,nNbProcs,nMPISides_rec,nbProc
USE MOD_Mesh_Tools             ,ONLY: GetMasteriLocSides
#endif /*USE_MPI*/
USE MOD_Mesh_Vars              ,ONLY: GlobalUniqueSideID
USE MOD_Analyze_Vars           ,ONLY: CalcElectricTimeDerivative
#ifdef PARTICLES
USE MOD_HDG_Vars               ,ONLY: UseBiasVoltage,BiasVoltage,BVDataLength
USE MOD_PICInterpolation_Vars  ,ONLY: useAlgebraicExternalField,AlgebraicExternalField
USE MOD_Analyze_Vars           ,ONLY: AverageElectricPotential
USE MOD_Mesh_Vars              ,ONLY: Elem_xGP
USE MOD_HDG_Vars               ,ONLY: UseBRElectronFluid,BRAutomaticElectronRef,RegionElectronRef
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcElectronIonDensity,CalcElectronTemperature
USE MOD_Particle_Analyze_Tools ,ONLY: AllocateElectronIonDensityCell,AllocateElectronTemperatureCell
USE MOD_Particle_Analyze_Tools ,ONLY: CalculateElectronIonDensityCell,CalculateElectronTemperatureCell
USE MOD_HDF5_Output_Particles  ,ONLY: AddBRElectronFluidToPartSource
USE MOD_HDG_Vars               ,ONLY: CoupledPowerPotential,UseCoupledPowerPotential,CPPDataLength
#endif /*PARTICLES*/
#endif /*USE_HDG*/
USE MOD_Analyze_Vars           ,ONLY: OutputTimeFixed
USE MOD_Output_Vars            ,ONLY: DoWriteStateToHDF5
USE MOD_StringTools            ,ONLY: set_formatting,clear_formatting
#if (PP_nVar==8)
USE MOD_HDF5_Output_Fields     ,ONLY: WritePMLDataToHDF5
#endif
USE MOD_HDF5_Output_ElemData   ,ONLY: WriteAdditionalElemData
! IMPLICIT VARIABLE HANDLING
USE MOD_Analyze_Vars           ,ONLY: OutputErrorNormsToH5
USE MOD_HDF5_Output_Fields     ,ONLY: WriteErrorNormsToHDF5
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: PreviousTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
#if defined(PARTICLES) || USE_HDG
CHARACTER(LEN=255),ALLOCATABLE :: LocalStrVarNames(:)
INTEGER(KIND=IK)               :: nVar
#endif /*defined(PARTICLES)*/
#ifdef PARTICLES
REAL                           :: NumSpec(nSpecAnalyze),TmpArray(1,1)
INTEGER(KIND=IK)               :: SimNumSpec(nSpecAnalyze)
#if USE_HDG
REAL,ALLOCATABLE               :: CPPDataHDF5(:,:)
#endif /*USE_HDG*/
#endif /*PARTICLES*/
REAL                           :: StartT,EndT

#ifdef PP_POIS
REAL                           :: Utemp(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#elif USE_HDG
#if PP_nVar==1
REAL                           :: Utemp(1:4,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#elif PP_nVar==3
REAL                           :: Utemp(1:3,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#else /*PP_nVar=4*/
REAL                           :: Utemp(1:7,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#endif /*PP_nVar==1*/
#else
#ifndef maxwell
REAL,ALLOCATABLE               :: Utemp(:,:,:,:,:)
#endif /*not maxwell*/
#endif /*PP_POIS*/
REAL                           :: OutputTime_loc
REAL                           :: PreviousTime_loc
INTEGER(KIND=IK)               :: PP_nVarTmp
LOGICAL                        :: usePreviousTime_loc
#if USE_HDG
INTEGER                        :: iSide
INTEGER                        :: iGlobSide
INTEGER,ALLOCATABLE            :: SortedUniqueSides(:),GlobalUniqueSideID_tmp(:)
#if USE_MPI
LOGICAL,ALLOCATABLE            :: OutputSide(:)
INTEGER                        :: SideID_start, SideID_end,iNbProc,SendID
#endif /*USE_MPI*/
REAL,ALLOCATABLE               :: SortedLambda(:,:,:)          ! lambda, ((PP_N+1)^2,nSides)
INTEGER                        :: SortedOffset,SortedStart,SortedEnd
#ifdef PARTICLES
INTEGER                        :: i,j,k,iElem
REAL,ALLOCATABLE               :: BVDataHDF5(:,:)
#endif /*PARTICLES*/
REAL,ALLOCATABLE               :: FPCDataHDF5(:,:),EPCDataHDF5(:,:)
INTEGER                        :: nVarFPC,nVarEPC
#endif /*USE_HDG*/
!===================================================================================================================================
#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(3))
#endif /*EXTRAE*/
! set local variables for output and previous times
IF(OutputTimeFixed.GE.0.0)THEN ! use fixed output time supplied by user
  SWRITE(UNIT_StdOut,'(A,ES25.14E3,A2)',ADVANCE='NO')' (WriteStateToHDF5 for fixed output time :',OutputTimeFixed,') '
  OutputTime_loc   = OutputTimeFixed
  PreviousTime_loc = OutputTimeFixed
ELSE
  OutputTime_loc   = OutputTime
  IF(PRESENT(PreviousTime))PreviousTime_loc = PreviousTime
END IF

#ifdef PARTICLES
! Output lost particles if 1. lost during simulation     : NbrOfNewLostParticlesTotal > 0
!                          2. went missing during restart: TotalNbrOfMissingParticlesSum > 0
IF(CountNbrOfLostParts)THEN
  IF((NbrOfNewLostParticlesTotal.GT.0).OR.(TotalNbrOfMissingParticlesSum.GT.0))THEN
   CALL WriteLostParticlesToHDF5(MeshFileName,OutputTime_loc)
  END IF ! (NbrOfNewLostParticlesTotal.GT.0).OR.(TotalNbrOfMissingParticlesSum.GT.0)
END IF
! Output total number of particles here, if DoWriteStateToHDF5=F. Otherwise the info will be displayed at the end of this routine
IF(.NOT.DoWriteStateToHDF5)THEN
  ! Check if the total number of particles has already been determined
  IF(.NOT.GlobalNbrOfParticlesUpdated) CALL CalcNumPartsOfSpec(NumSpec,SimNumSpec,.FALSE.,.TRUE.)
  ! Output total number of particles here as the end of this routine will not be reached when DoWriteStateToHDF5 is false
  CALL DisplayNumberOfParticles(1)
END IF ! .NOT.DoWriteStateToHDF5
#endif /*PARTICLES*/

! Check if state file creation should be skipped
IF(.NOT.DoWriteStateToHDF5) RETURN

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE '
GETTIME(StartT)


! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime_loc))//'.h5'
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') '['//TRIM(FileName)//'] ...'
RestartFile=Filename
#if USE_HDG
#if PP_nVar==1
IF(MPIRoot) CALL GenerateFileSkeleton('State',4,StrVarNames,MeshFileName,OutputTime_loc)
#elif PP_nVar==3
IF(MPIRoot) CALL GenerateFileSkeleton('State',3,StrVarNames,MeshFileName,OutputTime_loc)
#else
IF(MPIRoot) CALL GenerateFileSkeleton('State',7,StrVarNames,MeshFileName,OutputTime_loc)
#endif
#else
IF(MPIRoot) CALL GenerateFileSkeleton('State',PP_nVar,StrVarNames,MeshFileName,OutputTime_loc)
#endif /*USE_HDG*/
! generate nextfile info in previous output file
usePreviousTime_loc=.FALSE.

IF(PRESENT(PreviousTime).AND.(.NOT.DoInitialAutoRestart))THEN
  usePreviousTime_loc=.TRUE.
  IF(MPIRoot .AND. PreviousTime_loc.LT.OutputTime_loc) CALL GenerateNextFileInfo('State',OutputTime_loc,PreviousTime_loc)
END IF

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif

! Associate construct for integer KIND=8 possibility
PP_nVarTmp = INT(PP_nVar,IK)
ASSOCIATE (&
      N                 => INT(PP_N,IK)               ,&
      nGlobalElems      => INT(nGlobalElems,IK)       ,&
      PP_nElems         => INT(PP_nElems,IK)          ,&
      offsetElem        => INT(offsetElem,IK)         ,&
      offsetSide        => INT(offsetSide,IK)         ,&
      nUniqueSides      => INT(nUniqueSides,IK)       ,&
      nGlobalUniqueSides=> INT(nGlobalUniqueSides,IK)  )

  ! Write DG solution ----------------------------------------------------------------------------------------------------------------
  !nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
  ! Store the Solution of the Maxwell-Poisson System
#ifdef PP_POIS
  ALLOCATE(Utemp(1:PP_nVar,0:N,0:N,0:N,PP_nElems))
#if (PP_nVar==8)
  Utemp(8,:,:,:,:)=Phi(1,:,:,:,:)
  Utemp(1:3,:,:,:,:)=E(1:3,:,:,:,:)
  Utemp(4:7,:,:,:,:)=U(4:7,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=Utemp)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionE', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)

  !CALL WriteArrayToHDF5('DG_SolutionPhi',nVal,5,(/4_IK,N+1,N+1,N+1,PP_nElems/) &
  !,offsetElem,5,existing=.FALSE.,RealArray=Phi)
  ! missing addiontal attributes and data preparation
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionPhi', rank=5,&
      nValGlobal=(/4_IK , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/4_IK , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==8)*/
  ! Store the solution of the electrostatic-poisson system
#if (PP_nVar==4)
  Utemp(1,:,:,:,:)=Phi(1,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=Utemp)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionE', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionPhi', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==4)*/
  DEALLOCATE(Utemp)
#elif USE_HDG

  ! Store lambda solution in sorted order by ascending global unique side ID
#if USE_MPI
  IF(nProcessors.GT.1)THEN
    ! 0. Store true/false info for each side if it should be written to h5 by each process
    ALLOCATE(OutputSide(1:nSides))
    OutputSide=.FALSE.

    ! 1. Flag BC and inner sides
    OutputSide(1:lastInnerSide) = .TRUE.

    ! 2. Flag MINE/YOUR sides that are sent to other procs and if their rank is larger this proc, it writes the data
    DO SendID = 1, 2
      DO iNbProc=1,nNbProcs
        IF(nMPISides_rec(iNbProc,SendID).GT.0)THEN
          SideID_start=OffsetMPISides_rec(iNbProc-1,SendID)+1
          SideID_end  =OffsetMPISides_rec(iNbProc,SendID)
          IF(nbProc(iNbProc).GT.myrank)THEN
            OutputSide(SideID_start:SideID_end) = .TRUE.
          END IF ! nbProc(iNbProc)
        END IF
      END DO !iProc=1,nNBProcs
    END DO ! SendID = 1, 2

    CALL GetMasteriLocSides()

  END IF ! nProcessors.GT.1
#endif /*USE_MPI*/

  ! Get mapping from side IDs to globally sorted unique side IDs
  ALLOCATE(SortedUniqueSides(1:nSides))
  ALLOCATE(GlobalUniqueSideID_tmp(1:nSides))
  SortedUniqueSides=0
  DO iSide = 1, nSides
    SortedUniqueSides(iSide)=iSide
  END DO ! iSide = 1, nSides

  ! Create tmp array which will be sorted
  GlobalUniqueSideID_tmp = GlobalUniqueSideID
  CALL QuickSortTwoArrays(1,nSides,GlobalUniqueSideID_tmp(1:nSides),SortedUniqueSides(1:nSides))
  DEALLOCATE(GlobalUniqueSideID_tmp)

  ! Fill array with lambda values in global unique side sorted order
  ALLOCATE(SortedLambda(PP_nVar,nGP_face,nSides))
  SortedLambda = HUGE(1.)
  ! This loop goes over all nSides and is labelled iGlobSide because all unique global sides that each processors outputs must be
  ! transformed into master side orientation for read-in during restart later on
  DO iGlobSide = 1, nSides
    ! Set side ID in processor local list
    iSide = SortedUniqueSides(iGlobSide)

#if USE_MPI
    ! Skip sides that are not processed by the current proc
    IF(nProcessors.GT.1)THEN
      ! Check if a side belongs to me (all BC and inner sides automatically included); at MPI interfaces the smaller rank wins and
      ! must output the data, because for these sides it is ambiguous
      IF(.NOT.OutputSide(iSide)) CYCLE
    END IF ! nProcessors.GT.1
#endif /*USE_MPI*/

    CALL LambdaSideToMaster(iSide,SortedLambda(:,:,iGlobSide))

  END DO ! iGlobSide = 1, nSides

  ! Deallocate temporary arrays
  DEALLOCATE(SortedUniqueSides)
  IF(nProcessors.GT.1) DEALLOCATE(iLocSides)


  ! Get offset and min/max index in sorted list
  SortedStart = 1
  SortedEnd   = nSides
  SortedOffset = 0 ! initialize

#if USE_MPI
  IF(nProcessors.GT.1)THEN
    SortedOffset=HUGE(1)
    DO iSide = 1, nSides
      ! Get local offset of global unique sides: the smallest global unique side ID
      IF(OutputSide(iSide))THEN
        IF(GlobalUniqueSideID(iSide).LT.SortedOffset) SortedOffset = GlobalUniqueSideID(iSide)
      ELSE
        ! the sum of non-output sides gives the beginning number of output sides for each proc
        SortedStart = SortedStart +1
      END IF ! OutputSide(iSide))
    END DO
    SortedOffset = SortedOffset-1
    DEALLOCATE(OutputSide)
  END IF ! nProcessors.GT.1
#endif /*USE_MPI*/

  ASSOCIATE( nGlobalOutputSides => INT(SortedEnd-SortedStart+1,IK) ,&
        SortedOffset => INT(SortedOffset,IK)            ,&
        SortedStart  => INT(SortedStart,IK)             ,&
        SortedEnd    => INT(SortedEnd,IK)               ,&
        nGP_face     => INT(nGP_face,IK)                )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
        DataSetName = 'DG_SolutionLambda', rank=3,&
        nValGlobal  = (/PP_nVarTmp , nGP_face , nGlobalUniqueSides/) , &
        nVal        = (/PP_nVarTmp , nGP_face , nGlobalOutputSides/)       , &
        offset      = (/0_IK       , 0_IK     , SortedOffset/)       , &
        collective  = .TRUE.                                         , &
        RealArray   = SortedLambda(:,:,SortedStart:SortedEnd))
  END ASSOCIATE
  DEALLOCATE(SortedLambda)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionU', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE., RealArray=U)

#if (PP_nVar==1)
#ifdef PARTICLES
  ! Adjust electric field for Landmark testcase
  IF(useAlgebraicExternalField.AND.AlgebraicExternalField.EQ.1)THEN
    DO iElem=1,INT(PP_nElems)
      DO k=0,INT(PP_N); DO j=0,INT(PP_N); DO i=0,INT(PP_N)
        ASSOCIATE( Ue => AverageElectricPotential ,&
              xe => 2.4e-2                        ,&
              x  => Elem_xGP(1,i,j,k,iElem))
          Utemp(1,i,j,k,iElem) = U(1,i,j,k,iElem) - x * Ue / xe
          Utemp(2,i,j,k,iElem) = E(1,i,j,k,iElem) + Ue / xe
        END ASSOCIATE
      END DO; END DO; END DO !i,j,k
    END DO !iElem
    Utemp(3:4,:,:,:,:) = E(2:3,:,:,:,:)
  ELSE
#endif /*PARTICLES*/
    Utemp(1,:,:,:,:)   = U(1,:,:,:,:)
    Utemp(2:4,:,:,:,:) = E(1:3,:,:,:,:)
#ifdef PARTICLES
  END IF
#endif /*PARTICLES*/

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/4_IK , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/4_IK , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)

#elif (PP_nVar==3)
  Utemp(1:3,:,:,:,:)=B(1:3,:,:,:,:)
  !CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,N+1,N+1,N+1,PP_nElems/) &
  !,offsetElem,5,existing=.TRUE.,RealArray=Utemp)
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/3_IK , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/3_IK , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)
#else /*(PP_nVar==4)*/
  Utemp(1,:,:,:,:)=U(4,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
  Utemp(5:7,:,:,:,:)=B(1:3,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/7_IK , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/7_IK , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)
#endif /*(PP_nVar==1)*/
#else
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)
#endif /*PP_POIS*/


#ifdef PARTICLES
  ! output of last source term
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  IF(OutputSource) THEN
#if USE_HDG
    ! Add BR electron fluid density to PartSource for output to state.h5
    IF(UseBRElectronFluid) CALL AddBRElectronFluidToPartSource()
#endif /*USE_HDG*/
    ! output of pure current and density
    ! not scaled with epsilon0 and c_corr
    nVar=4_IK
    ALLOCATE(LocalStrVarNames(1:nVar))
    LocalStrVarNames(1)='CurrentDensityX'
    LocalStrVarNames(2)='CurrentDensityY'
    LocalStrVarNames(3)='CurrentDensityZ'
    LocalStrVarNames(4)='ChargeDensity'
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesSource',INT(nVar,4),StrArray=LocalStrVarnames)
      CALL CloseDataFile()
    END IF
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
        DataSetName='DG_Source', rank=5,  &
        nValGlobal=(/nVar , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
        nVal=      (/nVar , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
        offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
        collective=.TRUE.,RealArray=PartSource)

    DEALLOCATE(LocalStrVarNames)
  END IF

#endif /*PARTICLES*/

#if USE_HDG
  ! Output temporal derivate of the electric field
  IF(CalcElectricTimeDerivative) THEN
    nVar=3_IK
    ALLOCATE(LocalStrVarNames(1:nVar))
    LocalStrVarNames(1)='TimeDerivativeElecDisplacementX'
    LocalStrVarNames(2)='TimeDerivativeElecDisplacementY'
    LocalStrVarNames(3)='TimeDerivativeElecDisplacementZ'
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesTimeDerivative',INT(nVar,4),StrArray=LocalStrVarnames)
      CALL CloseDataFile()
    END IF
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
        DataSetName='DG_TimeDerivative', rank=5,  &
        nValGlobal=(/nVar , N+1_IK , N+1_IK , N+1_IK , nGlobalElems/) , &
        nVal=      (/nVar , N+1_IK , N+1_IK , N+1_IK , PP_nElems/)    , &
        offset=    (/0_IK , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
        collective=.TRUE.,RealArray=Et(1:3,:,:,:,:))

    DEALLOCATE(LocalStrVarNames)
  END IF
#endif /*USE_HDG*/

END ASSOCIATE

#ifdef PARTICLES
CALL WriteParticleToHDF5(FileName)
IF(DoBoundaryParticleOutputHDF5) THEN
  IF (usePreviousTime_loc) THEN
    CALL WriteBoundaryParticleToHDF5(MeshFileName,OutputTime_loc,PreviousTime_loc)
  ELSE
    CALL WriteBoundaryParticleToHDF5(MeshFileName,OutputTime_loc)
  END IF
END IF
IF(UseAdaptiveBC.OR.(nPorousBC.GT.0)) CALL WriteAdaptiveInfoToHDF5(FileName)
CALL WriteVibProbInfoToHDF5(FileName)
IF(RadialWeighting%PerformCloning) CALL WriteClonesToHDF5(FileName)
IF (PartBound%OutputWallTemp) CALL WriteAdaptiveWallTempToHDF5(FileName)
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
! For restart purposes, store the electron bulk temperature in .h5 state
! Only root writes the container
IF(CalcBulkElectronTemp.AND.MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  TmpArray(1,1) = BulkElectronTemp
  CALL WriteArrayToHDF5( DataSetName = 'BulkElectronTemp' , rank = 2 , &
                         nValGlobal  = (/1_IK , 1_IK/)               , &
                         nVal        = (/1_IK , 1_IK/)               , &
                         offset      = (/0_IK , 0_IK/)               , &
                         collective  = .FALSE., RealArray = TmpArray(1,1))
  CALL CloseDataFile()
END IF ! CalcBulkElectronTempi.AND.MPIRoot
#if USE_HDG
IF(UseCoupledPowerPotential.AND.MPIRoot)THEN
  ALLOCATE(CPPDataHDF5(1:CPPDataLength,1))
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CPPDataHDF5(1:CPPDataLength,1) = CoupledPowerPotential(1:CPPDataLength)
  CALL WriteArrayToHDF5( DataSetName = 'CoupledPowerPotential' , rank = 2 , &
                         nValGlobal  = (/1_IK , INT(CPPDataLength,IK)/)   , &
                         nVal        = (/1_IK , INT(CPPDataLength,IK)/)   , &
                         offset      = (/0_IK , 0_IK/)                    , &
                         collective  = .FALSE., RealArray = CPPDataHDF5(1:CPPDataLength,1))
  CALL CloseDataFile()
  DEALLOCATE(CPPDataHDF5)
END IF ! CalcBulkElectronTempi.AND.MPIRoot
! Bias voltage
IF(UseBiasVoltage.AND.MPIRoot)THEN
  ALLOCATE(BVDataHDF5(1:BVDataLength,1))
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  BVDataHDF5(1:BVDataLength,1) = BiasVoltage%BVData(1:BVDataLength)
  CALL WriteArrayToHDF5( DataSetName = 'BiasVoltage' , rank = 2   , &
                         nValGlobal  = (/1_IK , INT(BVDataLength,IK)/), &
                         nVal        = (/1_IK , INT(BVDataLength,IK)/), &
                         offset      = (/0_IK , 0_IK/)                        , &
                         collective  = .FALSE., RealArray = BVDataHDF5(1:BVDataLength,1))
  CALL CloseDataFile()
  DEALLOCATE(BVDataHDF5)
END IF ! CalcBulkElectronTempi.AND.MPIRoot
#endif /*USE_HDG*/
#endif /*PARTICLES*/

#if USE_HDG
! Floating boundary condition: Store charge on each FPC
IF(UseFPC.AND.MPIRoot)THEN
  nVarFPC = FPC%nUniqueFPCBounds
  ALLOCATE(FPCDataHDF5(1:nVarFPC,1))
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  FPCDataHDF5(1:nVarFPC,1) = FPC%Charge(1:nVarFPC)
  CALL WriteArrayToHDF5( DataSetName = 'FloatingPotentialCharge' , rank = 2   , &
                         nValGlobal  = (/1_IK , INT(nVarFPC,IK)/), &
                         nVal        = (/1_IK , INT(nVarFPC,IK)/), &
                         offset      = (/0_IK , 0_IK/)                        , &
                         collective  = .FALSE., RealArray = FPCDataHDF5(1:nVarFPC,1))
  CALL CloseDataFile()
  DEALLOCATE(FPCDataHDF5)
END IF ! CalcBulkElectronTempi.AND.MPIRoot
! Electric potential condition: Store Voltage on each EPC
IF(UseEPC.AND.MPIRoot)THEN
  nVarEPC = EPC%nUniqueEPCBounds
  ALLOCATE(EPCDataHDF5(1:nVarEPC,1))
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  EPCDataHDF5(1:nVarEPC,1) = EPC%Voltage(1:nVarEPC)
  CALL WriteArrayToHDF5( DataSetName = 'ElectricPotenitalCondition' , rank = 2   , &
                         nValGlobal  = (/1_IK , INT(nVarEPC,IK)/), &
                         nVal        = (/1_IK , INT(nVarEPC,IK)/), &
                         offset      = (/0_IK , 0_IK/)                        , &
                         collective  = .FALSE., RealArray = EPCDataHDF5(1:nVarEPC,1))
  CALL CloseDataFile()
  DEALLOCATE(EPCDataHDF5)
END IF ! CalcBulkElectronTempi.AND.MPIRoot
#endif /*USE_HDG*/

#if USE_LOADBALANCE
! Write 'ElemTime' to a separate container in the state.h5 file
CALL WriteElemDataToSeparateContainer(FileName,ElementOut,'ElemTime')
#endif /*USE_LOADBALANCE*/

#if defined(PARTICLES) && USE_HDG
! Write 'ElectronDensityCell' and 'ElectronTemperatureCell' to a separate container in the state.h5 file
! (for special read-in and conversion to kinetic electrons)
IF(UseBRElectronFluid) THEN
  ! Check if electron density is already calculated in each cell
  IF(.NOT.CalcElectronIonDensity)THEN
    CALL AllocateElectronIonDensityCell()
    CALL CalculateElectronIonDensityCell()
  END IF
  CALL WriteElemDataToSeparateContainer(FileName,ElementOut,'ElectronDensityCell')

  ! Check if electron temperature is already calculated in each cell
  IF(.NOT.CalcElectronTemperature)THEN
    CALL AllocateElectronTemperatureCell()
    CALL CalculateElectronTemperatureCell()
  END IF
  CALL WriteElemDataToSeparateContainer(FileName,ElementOut,'ElectronTemperatureCell')

  ! Automatically obtain the reference parameters (from a fully kinetic simulation), store them in .h5 state
  ! Only root writes the container
  IF(BRAutomaticElectronRef.AND.MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteArrayToHDF5( DataSetName = 'RegionElectronRef' , rank = 2 , &
                           nValGlobal  = (/1_IK , 3_IK/)     , &
                           nVal        = (/1_IK , 3_IK/)     , &
                           offset      = (/0_IK , 0_IK/)     , &
                           collective  = .FALSE., RealArray = RegionElectronRef(1:3,1))
    CALL CloseDataFile()
  END IF !MPIRoot
END IF ! BRAutomaticElectronRef
#endif /*defined(PARTICLES) && USE_HDG*/

! Adjust values before WriteAdditionalElemData() is called
CALL ModifyElemData(mode=1)

! Write all 'ElemData' arrays to a single container in the state.h5 file
CALL WriteAdditionalElemData(FileName,ElementOut)

! Adjust values after WriteAdditionalElemData() is called
CALL ModifyElemData(mode=2)

#if (PP_nVar==8)
CALL WritePMLDataToHDF5(FileName)
#endif

#ifdef PARTICLES
! Write NodeSourceExt (external charge density) field to HDF5 file
IF(DoDielectricSurfaceCharge) CALL WriteNodeSourceExtToHDF5(OutputTime_loc)
! Output particle emission data to be read during subsequent restarts
CALL WriteEmissionVariablesToHDF5(FileName)
#endif /*PARTICLES*/

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

IF(OutputErrorNormsToH5) CALL WriteErrorNormsToHDF5(OutputTime_loc)

#if defined(PARTICLES)
CALL DisplayNumberOfParticles(1)
#endif /*defined(PARTICLES)*/

#ifdef EXTRAE
CALL extrae_eventandcounters(int(9000001), int8(0))
#endif /*EXTRAE*/
END SUBROUTINE WriteStateToHDF5


SUBROUTINE ModifyElemData(mode)
!===================================================================================================================================
!> Modify ElemData fields before/after WriteAdditionalElemData() is called
!===================================================================================================================================
! MODULES
USE MOD_TimeDisc_Vars         ,ONLY: Time
USE MOD_Restart_Vars          ,ONLY: RestartTime
#ifdef PARTICLES
USE MOD_Globals               ,ONLY: abort
USE MOD_Particle_Analyze_Vars ,ONLY: CalcCoupledPower,PCouplSpec
USE MOD_Particle_Vars         ,ONLY: nSpecies,Species
#endif /*PARTICLES*/
#if USE_MPI
USE MOD_Globals
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: mode ! 1: before WriteAdditionalElemData() is called
!                          ! 2: after WriteAdditionalElemData() is called
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef PARTICLES
REAL          :: timediff
INTEGER       :: iSpec
#else
INTEGER       :: dummy ! dummy variable for compiler warning suppression
#endif /*PARTICLES*/
!===================================================================================================================================

IF(ABS(Time-RestartTime).LE.0.0) RETURN

#ifdef PARTICLES
IF(mode.EQ.1)THEN
  timediff = 1.0 / (Time-RestartTime)
ELSEIF(mode.EQ.2)THEN
  timediff = (Time-RestartTime)
ELSE
  CALL abort( __STAMP__,'ModifyElemData: mode must be 1 or 2')
END IF ! mode.EQ.1

! Set coupled power to particles if output of coupled power is active
IF (CalcCoupledPower.AND.(timediff.GT.0.)) THEN
  DO iSpec = 1, nSpecies
    IF(ABS(Species(iSpec)%ChargeIC).GT.0.0)THEN
      PCouplSpec(iSpec)%DensityAvgElem = PCouplSpec(iSpec)%DensityAvgElem * timediff
    END IF
  END DO
END IF
#endif /*PARTICLES*/

#if !defined(PARTICLES)
! Suppress compiler warning
RETURN
dummy=mode
#endif /*!(PARTICLES)*/

END SUBROUTINE ModifyElemData


#if defined(PARTICLES)
SUBROUTINE WriteIMDStateToHDF5()
!===================================================================================================================================
! Write the particles data aquired from an IMD *.chkpt file to disk and abort the program
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars ,ONLY: IMDInputFile,IMDTimeScale,IMDLengthScale,IMDNumber
USE MOD_Mesh_Vars     ,ONLY: MeshFile
USE MOD_Restart_Vars  ,ONLY: DoRestart
#if USE_MPI
USE MOD_MPI           ,ONLY: FinalizeMPI
#endif /*USE_MPI*/
USE MOD_ReadInTools   ,ONLY: PrintOption
USE MOD_HDF5_Output_Particles,ONLY: FillParticleData
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255) :: tempStr
REAL               :: t,tFuture,IMDtimestep
INTEGER            :: iSTATUS,IMDanalyzeIter
!===================================================================================================================================
IF(.NOT.DoRestart)THEN
  IF(IMDTimeScale.GT.0.0)THEN
    SWRITE(UNIT_StdOut,'(A)')'   IMD: calc physical time in seconds for which the IMD *.chkpt file is defined.'
    ! calc physical time in seconds for which the IMD *.chkpt file is defined
    ! t = IMDanalyzeIter * IMDtimestep * IMDTimeScale * IMDNumber
    IMDtimestep=0.0
    CALL GetParameterFromFile(IMDInputFile,'timestep'   , TempStr ,DelimiterSymbolIN=' ',CommentSymbolIN='#')
    CALL str2real(TempStr,IMDtimestep,iSTATUS)
    IF(iSTATUS.NE.0)THEN
      CALL abort(&
      __STAMP__&
      ,'Could not find "timestep" in '//TRIM(IMDInputFile)//' for IMDtimestep!')
    END IF

    IMDanalyzeIter=0
    CALL GetParameterFromFile(IMDInputFile,'checkpt_int', TempStr ,DelimiterSymbolIN=' ',CommentSymbolIN='#')
    CALL str2int(TempStr,IMDanalyzeIter,iSTATUS)
    IF(iSTATUS.NE.0)THEN
      CALL abort(&
      __STAMP__&
      ,'Could not find "checkpt_int" in '//TRIM(IMDInputFile)//' for IMDanalyzeIter!')
    END IF
    CALL PrintOption('IMDtimestep'    , 'OUTPUT' , RealOpt=IMDtimestep)
    CALL PrintOption('IMDanalyzeIter' , 'OUTPUT' , IntOpt=IMDanalyzeIter)
    CALL PrintOption('IMDTimeScale'   , 'OUTPUT' , RealOpt=IMDTimeScale)
    CALL PrintOption('IMDLengthScale' , 'OUTPUT' , RealOpt=IMDLengthScale)
    CALL PrintOption('IMDNumber'      , 'OUTPUT' , IntOpt=IMDNumber)
    t = REAL(IMDanalyzeIter) * IMDtimestep * IMDTimeScale * REAL(IMDNumber)
    CALL PrintOption('t'              , 'OUTPUT' , RealOpt=t)
    SWRITE(UNIT_StdOut,'(A,ES25.14E3,A,F15.3,A)')     '   Calculated time t :',t,' (',t*1e12,' ps)'

    tFuture=t
    CALL FillParticleData()
    CALL WriteStateToHDF5(TRIM(MeshFile),t,tFuture)
    SWRITE(UNIT_StdOut,'(A)')'   Particles: StateFile (IMD MD data) created. Terminating successfully!'
  ELSE
    CALL abort(__STAMP__, ' IMDLengthScale.LE.0.0 which is not allowed')
  END IF
END IF
END SUBROUTINE WriteIMDStateToHDF5
#endif /*PARTICLES*/


#if USE_LOADBALANCE || defined(PARTICLES)
SUBROUTINE WriteElemDataToSeparateContainer(FileName,ElemList,ElemDataName)
!===================================================================================================================================
!> Similar to WriteAdditionalElemData() but only writes one of the fields to a separate container
!> ----------------
!> Write additional data for analyze purpose to HDF5.
!> The data is taken from a lists, containing either pointers to data arrays or pointers
!> to functions to generate the data, along with the respective varnames.
!>
!> Two options are available:
!>    1. WriteAdditionalElemData:
!>       Element-wise scalar data, e.g. the timestep or indicators.
!>       The data is collected in a single array and written out in one step.
!>       DO NOT MISUSE NODAL DATA FOR THIS! IT WILL DRASTICALLY INCREASE FILE SIZE AND SLOW DOWN IO!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars        ,ONLY: nElems
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,ElemTime_tmp,NullifyElemTime
USE MOD_Restart_Vars     ,ONLY: DoRestart
USE MOD_Mesh_Vars        ,ONLY: nGlobalElems,offsetelem
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)        :: FileName
TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList !< Linked list of arrays to write to file
CHARACTER(LEN=*),INTENT(IN)          :: ElemDataName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                   :: StrVarNames
REAL,ALLOCATABLE                     :: ElemData(:,:)
INTEGER                              :: nVar,iElem
TYPE(tElementOut),POINTER            :: e
!===================================================================================================================================

IF(.NOT. ASSOCIATED(ElemList)) RETURN

! Allocate variable names and data array
ALLOCATE(ElemData(1,nElems))

! Fill the arrays
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  StrVarNames=e%VarName
  IF(StrVarNames.EQ.TRIM(ElemDataName))THEN
    nVar=nVar+1
    IF(ASSOCIATED(e%RealArray))    ElemData(nVar,:)=e%RealArray(1:nElems)
    IF(ASSOCIATED(e%RealScalar))   ElemData(nVar,:)=e%RealScalar
    IF(ASSOCIATED(e%IntArray))     ElemData(nVar,:)=REAL(e%IntArray(1:nElems))
    IF(ASSOCIATED(e%IntScalar))    ElemData(nVar,:)=REAL(e%IntScalar)
    IF(ASSOCIATED(e%LongIntArray)) ElemData(nVar,:)=REAL(e%LongIntArray(1:nElems))
    IF(ASSOCIATED(e%LogArray)) THEN
      DO iElem=1,nElems
        IF(e%LogArray(iElem))THEN
          ElemData(nVar,iElem)=1.
        ELSE
          ElemData(nVar,iElem)=0.
        END IF
      END DO ! iElem=1,PP_nElems
    END IF
    IF(ASSOCIATED(e%eval))       CALL e%eval(ElemData(nVar,:)) ! function fills elemdata
    EXIT
  END IF ! StrVarNames.EQ.TRIM(ElemDataName)
  e=>e%next
END DO

IF(nVar.NE.1) CALL abort(__STAMP__,'WriteElemDataToSeparateContainer: Array not found in ElemData = '//TRIM(ElemDataName))

#if USE_LOADBALANCE
! Check if ElemTime is all zeros and if this is a restart (save the old values)
NullifyElemTime=.FALSE.
IF((MAXVAL(ElemData).LE.0.0)          .AND.& ! Restart
    DoRestart                         .AND.& ! Restart
    (TRIM(ElemDataName).EQ.'ElemTime').AND.& ! only for ElemTime array
    ALLOCATED(ElemTime_tmp))THEN             ! only allocated when not starting simulation from zero
  ! Additionally, store old values in ElemData container
  ElemTime = ElemTime_tmp
  NullifyElemTime=.TRUE. ! Set array to 0. after ElemData is written (but before ElemTime is measured again)

  ! Write 'ElemTime' container
  ASSOCIATE (&
        nVar         => INT(nVar,IK)         ,&
        nGlobalElems => INT(nGlobalElems,IK) ,&
        PP_nElems    => INT(PP_nElems,IK)    ,&
        offsetElem   => INT(offsetElem,IK)   )
    CALL GatheredWriteArray(FileName,create = .FALSE.,&
                            DataSetName     = TRIM(ElemDataName), rank = 2,  &
                            nValGlobal      = (/nVar,nGlobalElems/),&
                            nVal            = (/nVar,PP_nElems   /),&
                            offset          = (/0_IK,offsetElem  /),&
                            collective      = .TRUE.,RealArray        = ElemTime_tmp)
  END ASSOCIATE

ELSE
  ASSOCIATE (&
        nVar         => INT(nVar,IK)         ,&
        nGlobalElems => INT(nGlobalElems,IK) ,&
        PP_nElems    => INT(PP_nElems,IK)    ,&
        offsetElem   => INT(offsetElem,IK)   )
    CALL GatheredWriteArray(FileName,create = .FALSE.,&
                            DataSetName     = TRIM(ElemDataName), rank = 2,  &
                            nValGlobal      = (/nVar,nGlobalElems/),&
                            nVal            = (/nVar,PP_nElems   /),&
                            offset          = (/0_IK,offsetElem  /),&
                            collective      = .TRUE.,RealArray        = ElemData)
  END ASSOCIATE
END IF ! (MAXVAL(ElemData).LE.0.0).AND.DoRestart.AND.(TRIM(ElemDataName).EQ.'ElemTime')
#endif /*USE_LOADBALANCE*/

DEALLOCATE(ElemData)

END SUBROUTINE WriteElemDataToSeparateContainer
#endif /*USE_LOADBALANCE || defined(PARTICLES)*/


END MODULE MOD_HDF5_Output_State
