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
PUBLIC :: WriteStateToHDF5
#if !(PP_TimeDiscMethod==700)
PUBLIC :: WriteTimeAverage
#endif /*!(PP_TimeDiscMethod==700)*/
#if defined(PARTICLES)
PUBLIC :: WriteElemDataToSeparateContainer
#endif /*PARTICLES*/
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteStateToHDF5(MeshFileName,OutputTime,PreviousTime,InitialAutoRestartIn)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
#if USE_FV
USE MOD_FV_Vars                ,ONLY: U_FV
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems
#endif
#if !(USE_FV) || (USE_HDG)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
USE MOD_DG_Vars                ,ONLY: U_N
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Mesh_Vars              ,ONLY: offsetElem
#if USE_FV
#if USE_HDG
USE MOD_Equation_Vars          ,ONLY: StrVarNames
#endif
USE MOD_Equation_Vars_FV       ,ONLY: StrVarNames_FV
#else
USE MOD_Equation_Vars          ,ONLY: StrVarNames
#endif
USE MOD_Restart_Vars           ,ONLY: RestartFile,DoInitialAutoRestart
#ifdef PARTICLES
USE MOD_PICDepo_Vars           ,ONLY: OutputSource,PS_N
USE MOD_DSMC_Vars              ,ONLY: ParticleWeighting
USE MOD_Particle_Sampling_Vars ,ONLY: UseAdaptiveBC
USE MOD_SurfaceModel_Vars      ,ONLY: nPorousBC
USE MOD_Particle_Boundary_Vars ,ONLY: DoBoundaryParticleOutputHDF5, PartBound
USE MOD_Particle_Tracking_Vars ,ONLY: CountNbrOfLostParts,TotalNbrOfMissingParticlesSum,NbrOfNewLostParticlesTotal
USE MOD_Particle_Analyze_Vars  ,ONLY: nSpecAnalyze
USE MOD_Particle_Analyze_Tools ,ONLY: CalcNumPartsOfSpec
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
USE MOD_HDF5_Output_Particles_PIC  ,ONLY: WriteNodeSourceExtToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
USE MOD_HDF5_Output_Particles  ,ONLY: WriteClonesToHDF5,WriteVibProbInfoToHDF5,WriteAdaptiveWallTempToHDF5
USE MOD_HDF5_Output_Particles  ,ONLY: WriteAdaptiveInfoToHDF5,WriteParticleToHDF5,WriteBoundaryParticleToHDF5
USE MOD_HDF5_Output_Particles  ,ONLY: WriteLostParticlesToHDF5,WriteEmissionVariablesToHDF5
USE MOD_Particle_Vars          ,ONLY: CalcBulkElectronTemp,BulkElectronTemp
#endif /*PARTICLES*/
USE MOD_Mesh_Vars              ,ONLY: nElems
#ifdef discrete_velocity
USE MOD_FV_Vars                ,ONLY: doFVReconstruct
USE MOD_Gradients              ,ONLY: GetGradients
USE MOD_Prolong_FV             ,ONLY: ProlongToOutput
USE MOD_DistFunc               ,ONLY: MacroValuesFromDistribution
USE MOD_TimeDisc_Vars          ,ONLY: dt,time,dt_Min
USE MOD_Equation_Vars_FV       ,ONLY: DVMnSpecies, DVMnMacro, DVMnInnerE
#endif
#if USE_HDG
USE MOD_HDG_Vars               ,ONLY: UseFPC,FPC,UseEPC,EPC
#if PP_nVar==1
#elif PP_nVar==3
USE MOD_Equation_Vars          ,ONLY: B
#else
USE MOD_Equation_Vars          ,ONLY: E,B
#endif /*PP_nVar*/
USE MOD_Analyze_Vars           ,ONLY: CalcElectricTimeDerivative
#ifdef PARTICLES
USE MOD_HDG_Vars               ,ONLY: UseBiasVoltage,BiasVoltage,BVDataLength
USE MOD_PICInterpolation_Vars  ,ONLY: useAlgebraicExternalField,AlgebraicExternalField
USE MOD_Analyze_Vars           ,ONLY: AverageElectricPotential
USE MOD_Mesh_Vars              ,ONLY: N_VolMesh
USE MOD_HDG_Vars               ,ONLY: UseBRElectronFluid,BRAutomaticElectronRef,RegionElectronRef
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcElectronIonDensity,CalcElectronTemperature
USE MOD_Particle_Analyze_Tools ,ONLY: AllocateElectronIonDensityCell,AllocateElectronTemperatureCell
USE MOD_Particle_Analyze_Tools ,ONLY: CalculateElectronIonDensityCell,CalculateElectronTemperatureCell
USE MOD_HDF5_Output_Particles_HDG  ,ONLY: AddBRElectronFluidToPartSource
USE MOD_HDG_Vars               ,ONLY: CoupledPowerPotential,UseCoupledPowerPotential,CPPDataLength
#endif /*PARTICLES*/
#else
#endif /*USE_HDG*/
#if !(PP_TimeDiscMethod==700)
USE MOD_DG_vars                ,ONLY: N_DG_Mapping,nDofsMapping
#endif /*!(PP_TimeDiscMethod==700)*/
USE MOD_Interpolation_Vars     ,ONLY: Nmax
USE MOD_Analyze_Vars           ,ONLY: OutputTimeFixed
USE MOD_Output_Vars            ,ONLY: DoWriteStateToHDF5
USE MOD_StringTools            ,ONLY: set_formatting,clear_formatting
#if (PP_nVar==8)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
USE MOD_HDF5_Output_Fields_Maxwell,ONLY: WritePMLDataToHDF5
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif
USE MOD_HDF5_Output_ElemData   ,ONLY: WriteAdditionalElemData
USE MOD_ChangeBasis            ,ONLY : ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
USE MOD_Analyze_Vars           ,ONLY: OutputErrorNormsToH5
USE MOD_HDF5_Output_Fields     ,ONLY: WriteErrorNormsToHDF5
#if USE_HDG
#if defined(PARTICLES)
USE MOD_HDF5_Output_Fields_HDG ,ONLY: WriteSurfVDLToHDF5
USE MOD_Particle_Boundary_Vars ,ONLY: DoVirtualDielectricLayer
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: PreviousTime
LOGICAL,OPTIONAL               :: InitialAutoRestartIn
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
INTEGER                        :: nVarOut, NOut
INTEGER                        :: iElem, Nloc
REAL                           :: StartT,EndT,OutputTime_loc,PreviousTime_loc
LOGICAL                        :: usePreviousTime_loc,InitialAutoRestart
! p-adaption output
REAL,ALLOCATABLE               :: U_N_2D_local(:,:)
INTEGER                        :: iDOF, nDOFOutput, offsetDOF
#if defined(PARTICLES) || USE_HDG
CHARACTER(LEN=255),ALLOCATABLE :: LocalStrVarNames(:)
#endif /*defined(PARTICLES)*/
#ifdef PARTICLES
REAL                           :: NumSpec(nSpecAnalyze),TmpArray(1,1)
INTEGER(KIND=IK)               :: SimNumSpec(nSpecAnalyze)
#endif /*PARTICLES*/
#ifdef discrete_velocity
INTEGER(KIND=IK)               :: N_FV
REAL,ALLOCATABLE               :: Ureco(:,:,:,:,:)
REAL                           :: MacroVal(DVMnMacro,DVMnSpecies+1)
REAL,ALLOCATABLE               :: Udvm(:,:,:,:,:)
REAL                           :: tau,dtMV
INTEGER                        :: iSpec
INTEGER(KIND=IK)               :: nValDVM
REAL                           :: Trot(DVMnSpecies+1),Tvib(DVMnSpecies+1)
#endif /*discrete_velocity*/
INTEGER                        :: i,j,k
#if USE_HDG
REAL,ALLOCATABLE               :: FPCDataHDF5(:,:),EPCDataHDF5(:,:)
INTEGER                        :: nVarFPC,nVarEPC
#if defined(PARTICLES)
REAL,ALLOCATABLE               :: BVDataHDF5(:,:)
REAL,ALLOCATABLE               :: CPPDataHDF5(:,:)
#endif /*defined(PARTICLES)*/
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

IF(PRESENT(InitialAutoRestartIn)) THEN
  InitialAutoRestart = InitialAutoRestartIn
ELSE
  InitialAutoRestart = .FALSE.
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

#ifdef discrete_velocity
IF (doFVReconstruct) THEN
  N_FV = INT(PP_1,IK)
ELSE
  N_FV = 0_IK
ENDIF
ALLOCATE(Ureco(PP_nVar_FV,0:N_FV,0:N_FV,0:N_FV,PP_nElems))
#endif /*discrete_velocity*/

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_HDG
#ifdef discrete_velocity /*DVM+HDG*/
! poisson + DVM
nVarOut = 4
NOut = N_FV
ASSOCIATE( StrVarNames => [StrVarNames,StrVarNames_FV], &
           nVarOut     => nVarOut+(DVMnMacro+DVMnInnerE)*(DVMnSpecies+1)+1)
#elif PP_nVar==1
! poisson
nVarOut = 4
NOut = NMax
#elif PP_nVar==3
! magnetostatic
nVarOut = 3
NOut = PP_N
#else /*(PP_nVar==4)*/
! magnetostatic_poisson
nVarOut = 7
NOut = PP_N
#endif  /*PP_nVar==1*/
#elif defined(discrete_velocity) /*DVM*/
nVarOut = (DVMnMacro+DVMnInnerE)*(DVMnSpecies+1)+1
NOut = N_FV
ASSOCIATE( StrVarNames => StrVarNames_FV )
#else /*not USE_HDG*/
! maxwell
nVarOut = PP_nVar
NOut = NMax
#endif /*USE_HDG*/

IF(InitialAutoRestart) THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime_loc))//'_InitialRestart.h5'
  CALL GenerateFileSkeleton('State',nVarOut,StrVarNames,MeshFileName,OutputTime_loc,NIn=NOut,FileNameIn=FileName)
ELSE
  CALL GenerateFileSkeleton('State',nVarOut,StrVarNames,MeshFileName,OutputTime_loc,NIn=NOut,FileNameOut=FileName)
END IF
#if defined(discrete_velocity) /*DVM*/
END ASSOCIATE
#endif /*DVM*/

SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO') '['//TRIM(FileName)//'] ...'
RestartFile=Filename
! generate nextfile info in previous output file
usePreviousTime_loc=.FALSE.

IF(PRESENT(PreviousTime).AND.(.NOT.DoInitialAutoRestart))THEN
  usePreviousTime_loc=.TRUE.
  IF(MPIRoot .AND. PreviousTime_loc.LT.OutputTime_loc) CALL GenerateNextFileInfo('State',OutputTime_loc,PreviousTime_loc)
END IF

#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
! ---------------------------------------------------------
! Store lambda solution in sorted order by ascending global unique side ID [HDG]
! ---------------------------------------------------------
#if USE_HDG
CALL WriteLambdaSolutionSorted(FileName)
#endif /*USE_HDG*/
#if !(PP_TimeDiscMethod==700)
! ---------------------------------------------------------
! Preparing U_N_2D_local array for output as DG_Solution
! ---------------------------------------------------------
! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF

! Allocate local 2D array
ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))
U_N_2D_local = 0. ! Important: Initialize with zero to get the null-container output in the .h5 file if no actual data is written
#endif /*!(PP_TimeDiscMethod==700)*/

#if USE_HDG
#if (PP_nVar==1)
! ---------------------------------------------------------
! poisson
! ---------------------------------------------------------
#ifdef PARTICLES
  ! Adjust electric field for Landmark testcase
  IF(useAlgebraicExternalField.AND.AlgebraicExternalField.EQ.1)THEN
  ! Write into 2D array
  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      ! Correction for Phi
      U_N_2D_local(1,iDOF) = U_N(iElem)%U(1,i,j,k) - N_VolMesh(iElem)%Elem_xGP(1,i,j,k) * AverageElectricPotential / 2.4e-2
      ! Correction for Ex
      U_N_2D_local(2,iDOF) = U_N(iElem)%E(1,i,j,k) + AverageElectricPotential / 2.4e-2
      ! Ey and Ez are simply copied
      U_N_2D_local(3:4,iDOF) = U_N(iElem)%E(2:3,i,j,k)
    END DO; END DO; END DO
  END DO
  ELSE
#endif /*PARTICLES*/
  ! Write into 2D array
  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1,iDOF)   = U_N(iElem)%U(1,i,j,k)
      U_N_2D_local(2:4,iDOF) = U_N(iElem)%E(1:3,i,j,k)
    END DO; END DO; END DO
  END DO
#ifdef PARTICLES
  END IF
#endif /*PARTICLES*/

#elif (PP_nVar==3)
! ---------------------------------------------------------
! magnetostatic
! ---------------------------------------------------------
  Utemp(1:3,:,:,:,:)=B(1:3,:,:,:,:)
#else /*(PP_nVar==4)*/
! ---------------------------------------------------------
! magnetostatic_poisson
! ---------------------------------------------------------
  Utemp(1,:,:,:,:)=U(4,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
  Utemp(5:7,:,:,:,:)=B(1:3,:,:,:,:)
#endif /*(PP_nVar==1)*/
#else /*!USE_HDG*/
! ---------------------------------------------------------
! maxwell
! ---------------------------------------------------------
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
! Write into 2D array
iDOF = 0
DO iElem = 1, PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offsetElem)
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    iDOF = iDOF + 1
    U_N_2D_local(1:nVarOut,iDOF) = U_N(iElem)%U(1:nVarOut,i,j,k)
  END DO; END DO; END DO
END DO
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/
#endif /*USE_HDG*/
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
! ---------------------------------------------------------
! Output of DG_Solution
! TODO: currently an empty container is written for DSMC/BGK/FP
! ---------------------------------------------------------
! Associate construct for integer KIND=8 possibility
ASSOCIATE(nVarOut         => INT(nVarOut,IK)           ,&
          nDofsMapping    => INT(nDofsMapping,IK)      ,&
          nDOFOutput      => INT(nDOFOutput,IK)        ,&
          offsetDOF       => INT(offsetDOF,IK)         )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName = 'DG_Solution' , rank = 2                , &
                        nValGlobal  = (/nVarOut     , nDofsMapping/)          , &
                        nVal        = (/nVarOut     , nDOFOutput/)            , &
                        offset      = (/0_IK        , offsetDOF/)             , &
                        collective  = .TRUE.        , RealArray = U_N_2D_local)
END ASSOCIATE
SDEALLOCATE(U_N_2D_local)
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/

#if USE_FV
! Associate construct for integer KIND=8 possibility
ASSOCIATE(nGlobalElems    => INT(nGlobalElems,IK)       ,&
          PP_nElems       => INT(PP_nElems,IK)          ,&
          offsetElem      => INT(offsetElem,IK)         ,&
          PP_nVarTmp_FV   => INT(PP_nVar_FV,IK)          )
#ifdef discrete_velocity
IF (doFVReconstruct) THEN
! reconstruct solution using gradients, as done during simulation
  CALL GetGradients(U_FV(:,:),output=.TRUE.)
  CALL ProlongToOutput(U_FV,Ureco)
ELSE
  ! first order solution
  Ureco(:,0,0,0,:) = U_FV(:,:)
END IF
IF (time.EQ.0.) THEN
  dtMV = 0.
ELSE
  dtMV = dt
ENDIF
ALLOCATE(Udvm(1:(DVMnMacro+DVMnInnerE)*(DVMnSpecies+1)+1,0:N_FV,0:N_FV,0:N_FV,PP_nElems))
DO iElem=1,INT(PP_nElems)
  DO k=0,N_FV
    DO j=0,N_FV
      DO i=0,N_FV
        CALL MacroValuesFromDistribution(MacroVal,Ureco(:,i,j,k,iElem),dtMV,tau,1,Trot=Trot,Tvib=Tvib)
        DO iSpec=1,DVMnSpecies+1 !n species + total values
          Udvm((DVMnMacro+DVMnInnerE)*(iSpec-1)+1:(DVMnMacro+DVMnInnerE)*iSpec-DVMnInnerE,i,j,k,iElem) = MacroVal(1:DVMnMacro,iSpec)
          IF (DVMnInnerE.GT.0) Udvm((DVMnMacro+DVMnInnerE)*iSpec-DVMnInnerE+1,i,j,k,iElem) = Trot(iSpec)
          IF (DVMnInnerE.GT.1) Udvm((DVMnMacro+DVMnInnerE)*iSpec-DVMnInnerE+2,i,j,k,iElem) = Tvib(iSpec)
        END DO
        Udvm((DVMnMacro+DVMnInnerE)*(DVMnSpecies+1)+1,i,j,k,iElem) = dt_Min(DT_MIN)/tau
      END DO
    END DO
  END DO
END DO
SDEALLOCATE(Ureco)
nValDVM = INT((DVMnMacro+DVMnInnerE)*(DVMnSpecies+1)+1,IK)
CALL GatheredWriteArray(FileName,create=.FALSE.,&
    DataSetName='DVM_Solution', rank=5,&
    nValGlobal=(/nValDVM, N_FV+1_IK , N_FV+1_IK , N_FV+1_IK , nGlobalElems/) , &
    nVal=      (/nValDVM, N_FV+1_IK , N_FV+1_IK , N_FV+1_IK , PP_nElems/)    , &
    offset=    (/0_IK       , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
    collective=.TRUE.,RealArray=Udvm)
#endif
#ifdef drift_diffusion
CALL GatheredWriteArray(FileName,create=.FALSE.,&
    DataSetName='DriftDiffusion_Solution', rank=2,&
    nValGlobal=(/PP_nVarTmp_FV , nGlobalElems/) , &
    nVal=      (/PP_nVarTmp_FV ,  PP_nElems/)    , &
    offset=    (/0_IK          , offsetElem/)   , &
    collective=.TRUE.,RealArray=U_FV)
#endif
END ASSOCIATE
#endif /*USE_FV*/

! ---------------------------------------------------------
! output of last source term
! ---------------------------------------------------------
#ifdef PARTICLES
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
  IF(OutputSource) THEN
  nVarOut = 4
#if USE_HDG
    ! Add BR electron fluid density to PartSource for output to state.h5
    IF(UseBRElectronFluid) CALL AddBRElectronFluidToPartSource()
#endif /*USE_HDG*/
    ! output of pure current and density
    ! not scaled with epsilon0 and c_corr
  ALLOCATE(LocalStrVarNames(1:nVarOut))
    LocalStrVarNames(1)='CurrentDensityX'
    LocalStrVarNames(2)='CurrentDensityY'
    LocalStrVarNames(3)='CurrentDensityZ'
    LocalStrVarNames(4)='ChargeDensity'
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesSource',nVarOut,StrArray=LocalStrVarnames)
      CALL CloseDataFile()
    END IF

  ! Allocate local 2D array
  ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))

  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVarOut,iDOF) = PS_N(iElem)%PartSource(1:nVarOut,i,j,k)
    END DO; END DO; END DO
  END DO

  ! Output of DG_Source
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVarOut         => INT(nVarOut,IK)      , &
            nDofsMapping    => INT(nDofsMapping,IK) , &
            nDOFOutput      => INT(nDOFOutput,IK)   , &
            offsetDOF       => INT(offsetDOF,IK)    )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = 'DG_Source' , rank = 2                , &
                          nValGlobal  = (/nVarOut   , nDofsMapping/)          , &
                          nVal        = (/nVarOut   , nDOFOutput/)            , &
                          offset      = (/0_IK      , offsetDOF/)             , &
                          collective  = .TRUE.      , RealArray = U_N_2D_local)
  END ASSOCIATE
  DEALLOCATE(U_N_2D_local)
    DEALLOCATE(LocalStrVarNames)
  END IF
#endif /*PARTICLES*/

! ---------------------------------------------------------
! Output temporal derivate of the electric field
! ---------------------------------------------------------
#if USE_HDG
  IF(CalcElectricTimeDerivative) THEN
  nVarOut = 3
  ALLOCATE(LocalStrVarNames(1:nVarOut))
    LocalStrVarNames(1)='TimeDerivativeElecDisplacementX'
    LocalStrVarNames(2)='TimeDerivativeElecDisplacementY'
    LocalStrVarNames(3)='TimeDerivativeElecDisplacementZ'
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesTimeDerivative',nVarOut,StrArray=LocalStrVarnames)
      CALL CloseDataFile()
    END IF

  ! Allocate local 2D array
  ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))

  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVarOut,iDOF) = U_N(iElem)%Dt(1:nVarOut,i,j,k)
    END DO; END DO; END DO
  END DO

  ! Output of DG_TimeDerivative
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVarOut         => INT(nVarOut,IK)      , &
            nDofsMapping    => INT(nDofsMapping,IK) , &
            nDOFOutput      => INT(nDOFOutput,IK)   , &
            offsetDOF       => INT(offsetDOF,IK)    )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = 'DG_TimeDerivative' , rank = 2        , &
                          nValGlobal  = (/nVarOut   , nDofsMapping/)          , &
                          nVal        = (/nVarOut   , nDOFOutput/)            , &
                          offset      = (/0_IK      , offsetDOF/)             , &
                          collective  = .TRUE.      , RealArray = U_N_2D_local)
  END ASSOCIATE
  DEALLOCATE(U_N_2D_local)
    DEALLOCATE(LocalStrVarNames)
  END IF
! ---------------------------------------------------------
! Calculate the electric VDL surface potential from the particle and electric displacement current
! ---------------------------------------------------------
#if defined(PARTICLES)
IF(DoVirtualDielectricLayer)THEN
  nVarOut = 3
  ALLOCATE(LocalStrVarNames(1:nVarOut))
  LocalStrVarNames(1)='PhiFx'
  LocalStrVarNames(2)='PhiFy'
  LocalStrVarNames(3)='PhiFz'
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPhiF',nVarOut,StrArray=LocalStrVarnames)
    CALL CloseDataFile()
  END IF

  ! Allocate local 2D array
  ALLOCATE(U_N_2D_local(1:nVarOut,1:nDOFOutput))

  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVarOut,iDOF) = U_N(iElem)%PhiF(1:nVarOut,i,j,k)
    END DO; END DO; END DO
  END DO

  ! Output of DG_PhiF
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVarOut         => INT(nVarOut,IK)      , &
            nDofsMapping    => INT(nDofsMapping,IK) , &
            nDOFOutput      => INT(nDOFOutput,IK)   , &
            offsetDOF       => INT(offsetDOF,IK)    )
  CALL GatheredWriteArray(FileName,create = .FALSE.                           , &
                          DataSetName = 'DG_PhiF'   , rank = 2                , &
                          nValGlobal  = (/nVarOut   , nDofsMapping/)          , &
                          nVal        = (/nVarOut   , nDOFOutput/)            , &
                          offset      = (/0_IK      , offsetDOF/)             , &
                          collective  = .TRUE.      , RealArray = U_N_2D_local)
END ASSOCIATE
  DEALLOCATE(U_N_2D_local)
  DEALLOCATE(LocalStrVarNames)
END IF ! DoVirtualDielectricLayer
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
! ---------------------------------------------------------
! Write the particle data
! ---------------------------------------------------------
#ifdef PARTICLES
CALL WriteParticleToHDF5(FileName)
! ---------------------------------------------------------
! Additional DSMC-related output
! ---------------------------------------------------------
IF(UseAdaptiveBC.OR.(nPorousBC.GT.0)) CALL WriteAdaptiveInfoToHDF5(FileName)
CALL WriteVibProbInfoToHDF5(FileName)
IF(ParticleWeighting%PerformCloning) CALL WriteClonesToHDF5(FileName)
IF (PartBound%OutputWallTemp) CALL WriteAdaptiveWallTempToHDF5(FileName)
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/
! ---------------------------------------------------------
! Field boundary condition related output
! TODO: BulkElectronTemp can be output as an attribute
! TODO: Simplify output of CoupledPowerPotential, BiasVoltage, FloatingPotentialCharge, and ElectricPotentialCondition to rank = 1
! ---------------------------------------------------------
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
  CALL WriteArrayToHDF5( DataSetName = 'ElectricPotentialCondition' , rank = 2   , &
                         nValGlobal  = (/1_IK , INT(nVarEPC,IK)/), &
                         nVal        = (/1_IK , INT(nVarEPC,IK)/), &
                         offset      = (/0_IK , 0_IK/)                        , &
                         collective  = .FALSE., RealArray = EPCDataHDF5(1:nVarEPC,1))
  CALL CloseDataFile()
  DEALLOCATE(EPCDataHDF5)
END IF ! CalcBulkElectronTempi.AND.MPIRoot
#endif /*USE_HDG*/

! ---------------------------------------------------------
! Write 'ElemTime' to a separate container in the state.h5 file
! ---------------------------------------------------------
#if USE_LOADBALANCE
CALL WriteElemDataToSeparateContainer(FileName,ElementOut,'ElemTime')
#endif /*USE_LOADBALANCE*/
! ---------------------------------------------------------
! Output for the Boltzmann relation
! ---------------------------------------------------------
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

! ---------------------------------------------------------
! Output for the perfectly matched layer (PML) [maxwell]
! ---------------------------------------------------------
#if (PP_nVar==8)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
CALL WritePMLDataToHDF5(FileName)
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*(PP_nVar==8)*/

! ---------------------------------------------------------
! Write NodeSourceExt (external charge density) field to HDF5 file
! ---------------------------------------------------------
#ifdef PARTICLES
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
IF(DoDielectricSurfaceCharge) CALL WriteNodeSourceExtToHDF5(OutputTime_loc)
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
! ---------------------------------------------------------
! Output particle emission data to be read during subsequent restarts
! ---------------------------------------------------------
CALL WriteEmissionVariablesToHDF5(FileName)
#endif /*PARTICLES*/

IF (MPIRoot) CALL MarkWriteSuccessful(FileName)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

! ---------------------------------------------------------
! Output to separate files
! ---------------------------------------------------------
! ---------------------------------------------------------
! Boundary impacting particle data (output first to reduce the required memory)
! ---------------------------------------------------------
#if defined(PARTICLES)
IF(DoBoundaryParticleOutputHDF5) THEN
  IF (usePreviousTime_loc) THEN
    CALL WriteBoundaryParticleToHDF5(MeshFileName,OutputTime_loc,PreviousTime_loc)
  ELSE
    CALL WriteBoundaryParticleToHDF5(MeshFileName,OutputTime_loc)
  END IF
END IF
#endif /*defined(PARTICLES)*/
! ---------------------------------------------------------
! Output of error norms
! ---------------------------------------------------------
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))
IF(OutputErrorNormsToH5) CALL WriteErrorNormsToHDF5(OutputTime_loc)
! ---------------------------------------------------------
! Output for the virtual dielectric layer (VDL)
! ---------------------------------------------------------
#if USE_HDG
#if defined(PARTICLES)
IF(DoVirtualDielectricLayer) CALL WriteSurfVDLToHDF5(OutputTime_loc)
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700))*/

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
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars        ,ONLY: nGlobalElems,offsetelem
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
#endif /*USE_LOADBALANCE*/
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
#if USE_LOADBALANCE
END IF ! (MAXVAL(ElemData).LE.0.0).AND.DoRestart.AND.(TRIM(ElemDataName).EQ.'ElemTime')
#endif /*USE_LOADBALANCE*/

DEALLOCATE(ElemData)

END SUBROUTINE WriteElemDataToSeparateContainer
#endif /*USE_LOADBALANCE || defined(PARTICLES)*/


#if USE_HDG
SUBROUTINE WriteLambdaSolutionSorted(FileName)
!===================================================================================================================================
!> Store lambda solution in sorted order by ascending global unique side ID
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nSides,GlobalUniqueSideID,N_SurfMesh
USE MOD_Utils                  ,ONLY: QuickSortTwoArrays
USE MOD_Mesh_Vars              ,ONLY: nGlobalUniqueSides
USE MOD_HDG_Vars               ,ONLY: nGP_face,iLocSides
USE MOD_Interpolation_Vars     ,ONLY: Nmax
USE MOD_Mesh_Tools             ,ONLY: LambdaSideToMaster
#if USE_MPI
USE MOD_Mesh_Vars              ,ONLY: lastInnerSide
USE MOD_MPI_Vars               ,ONLY: OffsetMPISides_rec,nNbProcs,nMPISides_rec,nbProc
USE MOD_Mesh_Tools             ,ONLY: GetMasteriLocSides
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iSide,NSide
INTEGER                        :: iGlobSide
INTEGER,ALLOCATABLE            :: SortedUniqueSides(:),GlobalUniqueSideID_tmp(:)
#if USE_MPI
LOGICAL,ALLOCATABLE            :: OutputSide(:)
INTEGER                        :: SideID_start, SideID_end,iNbProc,SendID
#endif /*USE_MPI*/
REAL,ALLOCATABLE               :: SortedLambda(:,:,:)          ! lambda, ((PP_N+1)^2,nSides)
INTEGER                        :: SortedOffset,SortedStart,SortedEnd
!===================================================================================================================================
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
ALLOCATE(SortedLambda(PP_nVar,nGP_face(NMax),nSides))
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

  NSide = N_SurfMesh(iSide)%NSide
  CALL LambdaSideToMaster(0,iSide,SortedLambda(:,:,iGlobSide),NSide)

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

ASSOCIATE(nGlobalOutputSides  => INT(SortedEnd-SortedStart+1,IK)  ,&
          SortedOffset        => INT(SortedOffset,IK)             ,&
          SortedStart         => INT(SortedStart,IK)              ,&
          SortedEnd           => INT(SortedEnd,IK)                ,&
          nVar                => INT(PP_nVar,IK)                  ,&
          nGP_face            => INT(nGP_face(NMax),IK)           ,&
          nGlobalUniqueSides  => INT(nGlobalUniqueSides,IK)       )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName = 'DG_SolutionLambda', rank=3,&
      nValGlobal  = (/nVar, nGP_face , nGlobalUniqueSides/) , &
      nVal        = (/nVar, nGP_face , nGlobalOutputSides/) , &
      offset      = (/0_IK, 0_IK     , SortedOffset/)       , &
      collective  = .TRUE.                                  , &
      RealArray   = SortedLambda(:,:,SortedStart:SortedEnd))
END ASSOCIATE
DEALLOCATE(SortedLambda)

END SUBROUTINE WriteLambdaSolutionSorted
#endif /*USE_HDG*/


#if !(PP_TimeDiscMethod==700)
SUBROUTINE WriteTimeAverage(MeshFileName,OutputTime,PreviousTime,VarNamesAvg,VarNamesFluc,nVar_Avg,nVar_Fluc)
!==================================================================================================================================
!> Subroutine to write time averaged data and fluctuations HDF5 format
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nElems
USE MOD_DG_vars              ,ONLY: N_DG_Mapping,nDofsMapping
USE MOD_HDF5_Output_ElemData ,ONLY: WriteAdditionalElemData
USE MOD_Timeaverage_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVar_Avg                                     !< Number of averaged variables
INTEGER,INTENT(IN)             :: nVar_Fluc                                    !< Number of fluctuations
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName                                 !< Name of mesh file
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesAvg(nVar_Avg)                        !< Average variable names
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesFluc(nVar_Fluc)                      !< Fluctuations variable names
REAL,INTENT(IN)                :: OutputTime                                   !< Time of output
REAL,INTENT(IN),OPTIONAL       :: PreviousTime                                 !< Time of previous output
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
REAL                           :: StartT,EndT
INTEGER                        :: nDOFOutput,offsetDOF,iDOF,i,j,k,Nloc,iElem
REAL,ALLOCATABLE               :: U_N_2D_local(:,:)
!==================================================================================================================================
IF((nVar_Avg.EQ.0).AND.(nVar_Fluc.EQ.0)) RETURN ! no time averaging

GETTIME(StartT)
SWRITE (UNIT_stdOut,'(A)',ADVANCE='NO') ' WRITE TIME AVERAGED STATE AND FLUCTUATIONS TO HDF5 FILE...'

! generate nextfile info in previous output file
IF(PRESENT(PreviousTime))THEN
  IF(MPIRoot .AND. PreviousTime.LT.OutputTime) CALL GenerateNextFileInfo('TimeAvg',OutputTime,PreviousTime)
END IF

! Get number of output DOFs per processor as the difference between first/last offset and add the number of DOFs of the last element
nDOFOutput = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
! Get the offset based on the element-local polynomial degree
IF(offsetElem.GT.0) THEN
  offsetDOF = N_DG_Mapping(1,1+offsetElem)
ELSE
  offsetDOF = 0
END IF

! Write timeaverages ---------------------------------------------------------------------------------------------------------------
IF(nVar_Avg.GT.0)THEN
  ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  CALL GenerateFileSkeleton('TimeAvg',nVar_Avg,VarNamesAvg,MeshFileName,OutputTime,FileNameOut=FileName)
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

  ! Allocate local 2D array
  ALLOCATE(U_N_2D_local(1:nVar_Avg,1:nDOFOutput))

  ! Write into 2D array
  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVar_Avg,iDOF)   = UAvg_N(iElem)%U(1:nVar_Avg,i,j,k)
    END DO; END DO; END DO
  END DO

  ! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
  CALL WriteAdditionalElemData(FileName,ElementOutNloc)

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVar_Avg        => INT(nVar_Avg,IK)          ,&
            nDofsMapping    => INT(nDofsMapping,IK)      ,&
            nDOFOutput      => INT(nDOFOutput,IK)        ,&
            offsetDOF       => INT(offsetDOF,IK)         )
  CALL GatheredWriteArray(FileName,create = .FALSE.                    , &
                          DataSetName = 'DG_Solution' , rank = 2       , &
                          nValGlobal  = (/nVar_Avg    , nDofsMapping/) , &
                          nVal        = (/nVar_Avg    , nDOFOutput/)   , &
                          offset      = (/0_IK        , offsetDOF/)    , &
                          collective  = .TRUE.        , RealArray = U_N_2D_local)
  END ASSOCIATE
  SDEALLOCATE(U_N_2D_local)
END IF

! Write fluctuations ---------------------------------------------------------------------------------------------------------------
IF(nVar_Fluc.GT.0)THEN
  CALL GenerateFileSkeleton('Fluc',nVar_Fluc,VarNamesFluc,MeshFileName,OutputTime,FileNameOut=FileName)
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

  ! Allocate local 2D array
  ALLOCATE(U_N_2D_local(1:nVar_Fluc,1:nDOFOutput))

  ! Write into 2D array
  iDOF = 0
  DO iElem = 1, PP_nElems
    Nloc = N_DG_Mapping(2,iElem+offsetElem)
    DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
      iDOF = iDOF + 1
      U_N_2D_local(1:nVar_Fluc,iDOF) = UFluc_N(iElem)%U(1:nVar_Fluc,i,j,k)
    END DO; END DO; END DO
  END DO

  ! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
  CALL WriteAdditionalElemData(FileName,ElementOutNloc)

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVar_Fluc       => INT(nVar_Fluc,IK)         ,&
            nDofsMapping    => INT(nDofsMapping,IK)      ,&
            nDOFOutput      => INT(nDOFOutput,IK)        ,&
            offsetDOF       => INT(offsetDOF,IK)         )
  CALL GatheredWriteArray(FileName,create = .FALSE.                    , &
                          DataSetName = 'DG_Solution' , rank = 2       , &
                          nValGlobal  = (/nVar_Fluc   , nDofsMapping/) , &
                          nVal        = (/nVar_Fluc   , nDOFOutput/)   , &
                          offset      = (/0_IK        , offsetDOF/)    , &
                          collective  = .TRUE.        , RealArray = U_N_2D_local)
  END ASSOCIATE
  SDEALLOCATE(U_N_2D_local)
END IF

IF (MPIRoot) CALL MarkWriteSuccessful(FileName)

GETTIME(endT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)
END SUBROUTINE WriteTimeAverage
#endif /*!(PP_TimeDiscMethod==700)*/

END MODULE MOD_HDF5_Output_State