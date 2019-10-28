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

MODULE MOD_HDF5_output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteStateToHDF5
  MODULE PROCEDURE WriteStateToHDF5
END INTERFACE

INTERFACE FlushHDF5
  MODULE PROCEDURE FlushHDF5
END INTERFACE

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

INTERFACE GenerateFileSkeleton
  MODULE PROCEDURE GenerateFileSkeleton
END INTERFACE

INTERFACE GenerateNextFileInfo
  MODULE PROCEDURE GenerateNextFileInfo
END INTERFACE

INTERFACE WriteTimeAverage
  MODULE PROCEDURE WriteTimeAverage
END INTERFACE

INTERFACE
  SUBROUTINE copy_userblock(outfilename,infilename) BIND(C)
      USE ISO_C_BINDING, ONLY: C_CHAR
      CHARACTER(KIND=C_CHAR) :: outfilename(*)
      CHARACTER(KIND=C_CHAR) :: infilename(*)
  END SUBROUTINE copy_userblock
END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteStateToHDF5,FlushHDF5,WriteHDF5Header,GatheredWriteArray
PUBLIC :: WriteArrayToHDF5,WriteAttributeToHDF5,GenerateFileSkeleton
PUBLIC :: WriteTimeAverage,GenerateNextFileInfo
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
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems
USE MOD_Equation_Vars          ,ONLY: StrVarNames
USE MOD_Restart_Vars           ,ONLY: RestartFile
#ifdef PARTICLES
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_PICDepo_Vars           ,ONLY: OutputSource,PartSource
USE MOD_Particle_Vars          ,ONLY: UseAdaptive
USE MOD_Particle_Boundary_Vars ,ONLY: nAdaptiveBC, nPorousBC
USE MOD_Dielectric_Vars        ,ONLY: DoDielectricSurfaceCharge
#endif /*PARTICLES*/
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: E,Phi
#endif /*PP_POIS*/
#if USE_HDG
USE MOD_Mesh_Vars              ,ONLY: offsetSide,nGlobalUniqueSides,nUniqueSides
USE MOD_HDG_Vars               ,ONLY: lambda, nGP_face
#if PP_nVar==1
USE MOD_Equation_Vars          ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars          ,ONLY: B
#else
USE MOD_Equation_Vars          ,ONLY: E,B
#endif /*PP_nVar*/
#endif /*USE_HDG*/
USE MOD_Analyze_Vars           ,ONLY: OutputTimeFixed
! IMPLICIT VARIABLE HANDLING
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
#ifdef PARTICLES
CHARACTER(LEN=255),ALLOCATABLE :: LocalStrVarNames(:)
INTEGER(KIND=IK)               :: nVar
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
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE...'
#if USE_MPI
StartT=MPI_WTIME()
#else
CALL CPU_TIME(StartT)
#endif

! set local variables for output and previous times
IF(OutputTimeFixed.GE.0.0)THEN
  SWRITE(UNIT_StdOut,'(A,ES25.14E3,A2)',ADVANCE='NO')' (WriteStateToHDF5 for fixed output time :',OutputTimeFixed,') '
  OutputTime_loc   = OutputTimeFixed
  PreviousTime_loc = OutputTimeFixed
ELSE
  OutputTime_loc   = OutputTime
  IF(PRESENT(PreviousTime))PreviousTime_loc = PreviousTime
END IF

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime_loc))//'.h5'
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
IF(PRESENT(PreviousTime))THEN
  IF(MPIRoot .AND. PreviousTime_loc.LT.OutputTime_loc) CALL GenerateNextFileInfo('State',OutputTime_loc,PreviousTime_loc)
END IF

! Reopen file and write DG solution
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

! Associate construct for integer KIND=8 possibility
PP_nVarTmp = INT(PP_nVar,IK)
ASSOCIATE (&
      PP_N         => INT(PP_N,IK)         ,&
      nGlobalElems => INT(nGlobalElems,IK) ,&
      PP_nElems    => INT(PP_nElems,IK)    ,&
      offsetElem   => INT(offsetElem,IK)   )

  ! Write DG solution ----------------------------------------------------------------------------------------------------------------
  !nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
  ! Store the Solution of the Maxwell-Poisson System
#ifdef PP_POIS
  ALLOCATE(Utemp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
#if (PP_nVar==8)
  Utemp(8,:,:,:,:)=Phi(1,:,:,:,:)
  Utemp(1:3,:,:,:,:)=E(1:3,:,:,:,:)
  Utemp(4:7,:,:,:,:)=U(4:7,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=Utemp)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionE', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)

  !CALL WriteArrayToHDF5('DG_SolutionPhi',nVal,5,(/4_IK,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
  !,offsetElem,5,existing=.FALSE.,RealArray=Phi)
  ! missing addiontal attributes and data preparation
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionPhi', rank=5,&
      nValGlobal=(/4_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/4_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==8)*/
  ! Store the solution of the electrostatic-poisson system
#if (PP_nVar==4)
  Utemp(1,:,:,:,:)=Phi(1,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=Utemp)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionE', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionPhi', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==4)*/
  DEALLOCATE(Utemp)
#elif USE_HDG
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionLambda', rank=3,&
      nValGlobal=(/PP_nVarTmp,nGP_face,nGlobalUniqueSides/),&
      nVal=      (/PP_nVarTmp,nGP_face,nUniqueSides/),&
      offset=    (/0_IK,      0_IK,       offsetSide/),&
      collective=.TRUE., RealArray=lambda(:,:,1:nUniqueSides))
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_SolutionU', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE., RealArray=U)
#if (PP_nVar==1)
  Utemp(1,:,:,:,:)=U(1,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/4_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/4_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)

#elif (PP_nVar==3)
  Utemp(1:3,:,:,:,:)=B(1:3,:,:,:,:)
  !CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
  !,offsetElem,5,existing=.TRUE.,RealArray=Utemp)
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/3_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/3_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)
#else /*(PP_nVar==4)*/
  Utemp(1,:,:,:,:)=U(4,:,:,:,:)
  Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
  Utemp(5:7,:,:,:,:)=B(1:3,:,:,:,:)

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/7_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/7_IK , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE., RealArray=Utemp)
#endif /*(PP_nVar==1)*/
#else
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
      DataSetName='DG_Solution', rank=5,&
      nValGlobal=(/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
      nVal=      (/PP_nVarTmp , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
      offset=    (/0_IK       , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
      collective=.TRUE.,RealArray=U)
#endif /*PP_POIS*/


#ifdef PARTICLES
  ! output of last source term
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
  IF(OutPutSource) THEN
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
        nValGlobal=(/nVar , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
        nVal=      (/nVar , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
        offset=    (/0_IK , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
        collective=.TRUE.,RealArray=PartSource)

    DEALLOCATE(LocalStrVarNames)
  END IF
#endif /*PARTICLES*/

END ASSOCIATE

#ifdef PARTICLES
CALL WriteParticleToHDF5(FileName)
CALL WriteMacroParticleToHDF5(FileName)
IF(UseAdaptive.OR.(nAdaptiveBC.GT.0).OR.(nPorousBC.GT.0)) CALL WriteAdaptiveInfoToHDF5(FileName)
CALL WriteVibProbInfoToHDF5(FileName)
CALL WriteSurfStateToHDF5(FileName)
IF(RadialWeighting%DoRadialWeighting) CALL WriteClonesToHDF5(FileName)
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/
#endif /*PARTICLES*/

#if USE_LOADBALANCE
! Write 'ElemTime' to a separate container in the state.h5 file
CALL WriteElemDataToSeparateContainer(FileName,ElementOut,'ElemTime')
#endif /*USE_LOADBALANCE*/


! Adjust values before WriteAdditionalElemData() is called
CALL ModifyElemData(mode=1)

! Write all 'ElemData' arrays to a single container in the state.h5 file
CALL WriteAdditionalElemData(FileName,ElementOut)

! Adjust values after WriteAdditionalElemData() is called
CALL ModifyElemData(mode=2)

#if (PP_nVar==8)
CALL WritePMLDataToHDF5(FileName)
#endif

#if PARTICLES
! Write NodeSourceExt (external charge density) field to HDF5 file
IF(DoDielectricSurfaceCharge) CALL WriteNodeSourceExtToHDF5(OutputTime_loc)
#endif /*PARTICLES*/

EndT=PICLASTIME()
SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'

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
  CALL abort(&
  __STAMP__&
  ,'ModifyElemData: mode must be 1 or 2')
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

#if !(PARTICLES)
! Suppress compiler warning
RETURN
dummy=mode
#endif /*!(PARTICLES)*/

END SUBROUTINE ModifyElemData


SUBROUTINE WriteAdditionalElemData(FileName,ElemList)
!===================================================================================================================================
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
USE MOD_Mesh_Vars     ,ONLY:offsetElem,nGlobalElems,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)        :: FileName
TYPE(tElementOut),POINTER,INTENT(IN) :: ElemList !< Linked list of arrays to write to file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL,ALLOCATABLE               :: ElemData(:,:)
INTEGER                        :: nVar,iElem
TYPE(tElementOut),POINTER      :: e
!===================================================================================================================================

IF(.NOT. ASSOCIATED(ElemList)) RETURN

! Count the additional variables
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  e=>e%next
END DO

! Allocate variable names and data array
ALLOCATE(StrVarNames(nVar))
ALLOCATE(ElemData(nVar,nElems))

! Fill the arrays
nVar = 0
e=>ElemList
DO WHILE(ASSOCIATED(e))
  nVar=nVar+1
  StrVarNames(nVar)=e%VarName
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
  e=>e%next
END DO

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

ASSOCIATE (&
      nVar         => INT(nVar,IK)         ,&
      nGlobalElems => INT(nGlobalElems,IK) ,&
      PP_nElems    => INT(PP_nElems,IK)    ,&
      offsetElem   => INT(offsetElem,IK)   )
  CALL GatheredWriteArray(FileName,create = .FALSE.,&
                          DataSetName     = 'ElemData', rank = 2,  &
                          nValGlobal      = (/nVar,nGlobalElems/),&
                          nVal            = (/nVar,PP_nElems   /),&
                          offset          = (/0_IK,offsetElem  /),&
                          collective      = .TRUE.,RealArray        = ElemData)
END ASSOCIATE
DEALLOCATE(ElemData,StrVarNames)

END SUBROUTINE WriteAdditionalElemData


#if USE_LOADBALANCE
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
USE MOD_HDF5_Input       ,ONLY: ReadArray
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: ElemTime,ElemTime_tmp
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

IF(nVar.NE.1) CALL abort(&
    __STAMP__&
    ,'WriteElemDataToSeparateContainer: Array not found in ElemData = '//TRIM(ElemDataName))

#if USE_LOADBALANCE
! Check if ElemTime is all zeros and if this is a restart (save the old values)
IF((MAXVAL(ElemData).LE.0.0)          .AND.& ! Restart
    DoRestart                         .AND.& ! Restart
    (TRIM(ElemDataName).EQ.'ElemTime').AND.& ! only for ElemTime array
    ALLOCATED(ElemTime_tmp))THEN             ! only allocated when not starting simulation from zero
  ! Additionally, store old values in ElemData container
  ElemTime = ElemTime_tmp

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
#endif /*USE_LOADBALANCE*/


#if (PP_nVar==8)
SUBROUTINE WritePMLDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:offsetElem,nGlobalElems
USE MOD_PML_Vars      ,ONLY:DoPML,PMLToElem,U2,nPMLElems,PMLnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
REAL,ALLOCATABLE               :: Upml(:,:,:,:,:)
INTEGER                        :: iPML
!===================================================================================================================================

IF(DoPML)THEN
  ALLOCATE(StrVarNames(PMLnVar))
  StrVarNames( 1)='PML-ElectricFieldX-P1'
  StrVarNames( 2)='PML-ElectricFieldX-P2'
  StrVarNames( 3)='PML-ElectricFieldX-P3'
  StrVarNames( 4)='PML-ElectricFieldY-P4'
  StrVarNames( 5)='PML-ElectricFieldY-P5'
  StrVarNames( 6)='PML-ElectricFieldY-P6'
  StrVarNames( 7)='PML-ElectricFieldZ-P7'
  StrVarNames( 8)='PML-ElectricFieldZ-P8'
  StrVarNames( 9)='PML-ElectricFieldZ-P9'
  StrVarNames(10)='PML-MagneticFieldX-P10'
  StrVarNames(11)='PML-MagneticFieldX-P11'
  StrVarNames(12)='PML-MagneticFieldX-P12'
  StrVarNames(13)='PML-MagneticFieldY-P13'
  StrVarNames(14)='PML-MagneticFieldY-P14'
  StrVarNames(15)='PML-MagneticFieldY-P15'
  StrVarNames(16)='PML-MagneticFieldZ-P16'
  StrVarNames(17)='PML-MagneticFieldZ-P17'
  StrVarNames(18)='PML-MagneticFieldZ-P18'
  StrVarNames(19)='PML-PhiB-P19'
  StrVarNames(20)='PML-PhiB-P20'
  StrVarNames(21)='PML-PhiB-P21'
  StrVarNames(22)='PML-PsiE-P22'
  StrVarNames(23)='PML-PsiE-P23'
  StrVarNames(24)='PML-PsiE-P24'

  ALLOCATE(UPML(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  UPML=0.0
  DO iPML=1,nPMLElems
    Upml(:,:,:,:,PMLToElem(iPML)) = U2(:,:,:,:,iPML)
  END DO ! iPML

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPML',PMLnVar,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      PP_N            => INT(PP_N,IK)            ,&
      PMLnVar         => INT(PMLnVar,IK)         ,&
      PP_nElems       => INT(PP_nElems,IK)       ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = 'PML_Solution', rank = 5,&
                          nValGlobal  = (/PMLnVar , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                          nVal        = (/PMLnVar , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems/)    , &
                          offset      = (/0_IK    , 0_IK   , 0_IK   , 0_IK   , offsetElem/)   , &
                          collective  = .TRUE.,RealArray = Upml)
END ASSOCIATE

!  CALL WriteArrayToHDF5(DataSetName='PML_Solution', rank=5,&
!                      nValGlobal=(/5,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
!                      nVal=      (/5,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
!                      offset=    (/0,      0,     0,     0,     offsetElem/),&
!                      collective=.TRUE., existing=.FALSE., RealArray=UPML)
!
!  CALL CloseDataFile()
  DEALLOCATE(UPML)
  DEALLOCATE(StrVarNames)
END IF ! DoPML


END SUBROUTINE WritePMLDataToHDF5
#endif

#ifdef PARTICLES
SUBROUTINE WriteParticleToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: nGlobalElems, offsetElem
USE MOD_Particle_Vars          ,ONLY: PDM, PEM, PartState, PartSpecies, PartMPF, usevMPF,PartPressureCell, nSpecies, VarTimeStep
USE MOD_part_tools             ,ONLY: UpdateNextFreePosition
USE MOD_DSMC_Vars              ,ONLY: UseDSMC, CollisMode,PartStateIntEn, DSMC, PolyatomMolDSMC, SpecDSMC, VibQuantsPar
#if (PP_TimeDiscMethod==509)
USE MOD_Particle_Vars          ,ONLY: velocityAtTime, velocityOutputAtTime
#endif /*(PP_TimeDiscMethod==509)*/
#if USE_MPI
USE MOD_Particle_MPI_Vars      ,ONLY: PartMPI
#endif /*USE_MPI*/
#ifdef CODE_ANALYZE
USE MOD_Particle_Tracking_Vars ,ONLY: PartOut,MPIRankOut
#endif /*CODE_ANALYZE*/
USE MOD_LoadBalance_Vars       ,ONLY: nPartsPerElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER                        :: nVar
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
INTEGER(KIND=IK)               :: nParticles(0:nProcessors-1)
#endif
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
LOGICAL                        :: withDSMC=.FALSE.
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
INTEGER                        :: iElem_glob, iElem_loc
INTEGER(KIND=IK),ALLOCATABLE   :: PartInt(:,:)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER, ALLOCATABLE           :: VibQuantData(:,:)
INTEGER,PARAMETER              :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER(KIND=IK)               :: locnPart_max
INTEGER                        :: MaxQuantNum, iPolyatMole, iSpec
!=============================================
! Required default values for KIND=IK
MaxQuantNum=-1
! Write properties -----------------------------------------------------------------------------------------------------------------
! Open dataset
!CALL H5DOPEN_F(File_ID,'DG_Solution',Dset_id,iError)

!!added for Evib, Erot writeout
withDSMC=useDSMC
IF (withDSMC) THEN
!IF (withDSMC) THEN
  IF ((CollisMode.GT.1).AND.(usevMPF) .AND. DSMC%ElectronicModel ) THEN !int ener + 3, vmpf +1
    PartDataSize=11
  ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. DSMC%ElectronicModel ) ) THEN !int ener + 2 and vmpf + 1
                                                                            ! or int energ +3 but no vmpf +1
    PartDataSize=10
  ELSE IF (CollisMode.GT.1) THEN
    PartDataSize=9 !int ener + 2
  ELSE IF (usevMPF) THEN
    PartDataSize=8 !+ 1 vmpf
  ELSE
    PartDataSize=7 !+ 0
  END IF
ELSE IF (usevMPF) THEN
  PartDataSize=8 !vmpf +1
ELSE
  PartDataSize=7
END IF

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

locnPart =   0_IK
DO pcount = 1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(pcount)) THEN
    locnPart = locnPart + 1_IK
  END IF
END DO

#if USE_MPI
sendbuf(1)=locnPart
recvbuf=0_IK
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart=recvbuf(1)
sendbuf(1)=recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob=sendbuf(1)
CALL MPI_GATHER(locnPart,1,MPI_INTEGER_INT_KIND,nParticles,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_WORLD,iError)
!IF (myRank.EQ.0) THEN
!  WRITE(*,*) 'PARTICLE-ELEMENT DISTRIBUTION'
!  WRITE(*,*) 'iProc, firstelemInd,   nElems,  locnPart,  totalnPart'
!  DO pcount=0,nProcessors-1
!    WRITE(*,'(I5,4I12)')pcount,offsetElemMPI(pcount),offsetElemMPI(pcount+1)-offsetElemMPI(pcount),&
!                       nParticles(pcount),SUM(nParticles(0:pcount))
!  END DO
!END IF
LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER_INT_KIND, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
#else
offsetnPart=0_IK
nPart_glob=locnPart
locnPart_max=locnPart
#endif
ALLOCATE(PartInt(offsetElem+1:offsetElem+PP_nElems,PartIntSize))
ALLOCATE(PartData(offsetnPart+1_IK:offsetnPart+locnPart,INT(PartDataSize,IK)))
IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ALLOCATE(VibQuantData(offsetnPart+1_IK:offsetnPart+locnPart,MaxQuantNum))
  VibQuantData = 0
  !+1 is real number of necessary vib quants for the particle
END IF

!!! Kleiner Hack von JN (Teil 1/2):

IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
  ALLOCATE(PEM%pStart(1:PP_nElems)           , &
           PEM%pNumber(1:PP_nElems)          , &
           PEM%pNext(1:PDM%maxParticleNumber), &
           PEM%pEnd(1:PP_nElems) )!            , &
           !PDM%nextUsedPosition(PDM%maxParticleNumber)  )
  useDSMC=.TRUE.
END IF
CALL UpdateNextFreePosition()
!!! Ende kleiner Hack von JN (Teil 1/2)
iPart=offsetnPart
DO iElem_loc=1,PP_nElems
  iElem_glob = iElem_loc + offsetElem
  PartInt(iElem_glob,1)=iPart
  IF (ALLOCATED(PEM%pNumber)) THEN
    nPartsPerElem(iElem_loc) = INT(PEM%pNumber(iElem_loc),IK)
    PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + INT(PEM%pNumber(iElem_loc),IK)
    pcount = PEM%pStart(iElem_loc)
    DO iPart=PartInt(iElem_glob,1)+1_IK,PartInt(iElem_glob,2)
      PartData(iPart,1)=PartState(pcount,1)
      PartData(iPart,2)=PartState(pcount,2)
      PartData(iPart,3)=PartState(pcount,3)
#if (PP_TimeDiscMethod==509)
      IF (velocityOutputAtTime) THEN
        PartData(iPart,4)=velocityAtTime(pcount,1)
        PartData(iPart,5)=velocityAtTime(pcount,2)
        PartData(iPart,6)=velocityAtTime(pcount,3)
      ELSE
#endif /*(PP_TimeDiscMethod==509)*/
      PartData(iPart,4)=PartState(pcount,4)
      PartData(iPart,5)=PartState(pcount,5)
      PartData(iPart,6)=PartState(pcount,6)
#if (PP_TimeDiscMethod==509)
      END IF
#endif /*(PP_TimeDiscMethod==509)*/
      PartData(iPart,7)=REAL(PartSpecies(pcount))
#ifdef CODE_ANALYZE
      IF(PARTOUT.GT.0 .AND. MPIRANKOUT.EQ.MyRank)THEN
        IF(pcount.EQ.PARTOUT)THEN
          PartData(iPart,7)=-PartData(iPart,7)
        END IF
      END IF
#endif /*CODE_ANALYZE*/
      IF (withDSMC) THEN
      !IF (withDSMC) THEN
        IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel) ) THEN
          PartData(iPart,8)=PartStateIntEn(pcount,1)
          PartData(iPart,9)=PartStateIntEn(pcount,2)
          PartData(iPart,10)=PartStateIntEn(pcount,3)
          PartData(iPart,11)=PartMPF(pcount)
        ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
          PartData(iPart,8)=PartStateIntEn(pcount,1)
          PartData(iPart,9)=PartStateIntEn(pcount,2)
          PartData(iPart,10)=PartMPF(pcount)
        ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
          PartData(iPart,8)=PartStateIntEn(pcount,1)
          PartData(iPart,9)=PartStateIntEn(pcount,2)
          PartData(iPart,10)=PartStateIntEn(pcount,3)
        ELSE IF (CollisMode.GT.1) THEN
          PartData(iPart,8)=PartStateIntEn(pcount,1)
          PartData(iPart,9)=PartStateIntEn(pcount,2)
        ELSE IF (usevMPF) THEN
          PartData(iPart,8)=PartMPF(pcount)
        END IF
      ELSE IF (usevMPF) THEN
          PartData(iPart,8)=PartMPF(pcount)
      END IF
      !PartData(iPart,8)=Species(PartSpecies(pcount))%ChargeIC*Species(PartSpecies(pcount))%MacroParticleFactor
      !PartData(iPart,9)=Species(PartSpecies(pcount))%MassIC*Species(PartSpecies(pcount))%MacroParticleFactor

      IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
        IF (SpecDSMC(PartSpecies(pcount))%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(PartSpecies(pcount))%SpecToPolyArray
          VibQuantData(iPart,1:PolyatomMolDSMC(iPolyatMole)%VibDOF) = &
            VibQuantsPar(pcount)%Quants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
        ELSE
           VibQuantData(iPart,:) = 0
        END IF
      END IF

      pcount = PEM%pNext(pcount)
    END DO
    iPart = PartInt(iElem_glob,2)
  ELSE
    CALL abort(&
    __STAMP__&
    , " Particle HDF5-Output method not supported! PEM%pNumber not associated")
  END IF
  PartInt(iElem_glob,2)=iPart
END DO

nVar=2
ALLOCATE(StrVarNames(nVar))
StrVarNames(1)='FirstPartID'
StrVarNames(2)='LastPartID'
!CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',2,StrArray=StrVarNames)

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesPartInt',nVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      nVar            => INT(nVar,IK)            ,&
      PP_nElems       => INT(PP_nElems,IK)       ,&
      offsetElem      => INT(offsetElem,IK)      ,&
      MaxQuantNum     => INT(MaxQuantNum,IK)     ,&
      PartDataSize    => INT(PartDataSize,IK)    )
  CALL GatheredWriteArray(FileName                         , create = .FALSE.            , &
                          DataSetName     = 'PartInt'      , rank   = 2                  , &
                          nValGlobal      = (/nGlobalElems , nVar/)                      , &
                          nVal            = (/PP_nElems    , nVar/)                      , &
                          offset          = (/offsetElem   , 0_IK/)                      , &
                          collective      = .TRUE.         , IntegerArray = PartInt)

  !  CALL WriteArrayToHDF5(DataSetName='PartInt', rank=2,&
  !                        nValGlobal=(/nGlobalElems, nVar/),&
  !                        nVal=      (/PP_nElems, nVar   /),&
  !                        offset=    (/offsetElem, 0_IK  /),&
  !                        collective=.TRUE., existing=.FALSE., IntegerArray=PartInt)
  !
  DEALLOCATE(StrVarNames)
  !CALL CloseDataFile()

  ALLOCATE(StrVarNames(PartDataSize))
  StrVarNames(1)='ParticlePositionX'
  StrVarNames(2)='ParticlePositionY'
  StrVarNames(3)='ParticlePositionZ'
  StrVarNames(4)='VelocityX'
  StrVarNames(5)='VelocityY'
  StrVarNames(6)='VelocityZ'
  StrVarNames(7)='Species'

  IF(withDSMC)THEN
    ! IF(withDSMC)THEN
    IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicModel))THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
      StrVarNames(11)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
    ELSE IF (CollisMode.GT.1) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
    ELSE IF (usevMPF) THEN
      StrVarNames( 8)='MPF'
    END IF
  ELSE IF (usevMPF) THEN
    StrVarNames( 8)='MPF'
  END IF

  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',INT(PartDataSize,4),StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

  IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
    IF(MPIRoot)THEN ! only root writes the container
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName='PartData'   , rank=2           , &
                            nValGlobal=(/nPart_glob  , PartDataSize/)   , &
                            nVal=      (/locnPart    , PartDataSize  /) , &
                            offset=    (/offsetnPart , 0_IK  /)         , &
                            collective=.FALSE.       , RealArray=PartData)
      CALL CloseDataFile()
    END IF !MPIRoot
  END IF !locnPart_max.EQ.0
#if USE_MPI
  CALL DistributedWriteArray(FileName                     , &
                             DataSetName  = 'PartData'    , rank = 2          , &
                             nValGlobal   = (/nPart_glob  , PartDataSize/)    , &
                             nVal         = (/locnPart    , PartDataSize/)    , &
                             offset       = (/offsetnPart , 0_IK/)            , &
                             collective   = .FALSE.       , offSetDim = 1     , &
                             communicator = PartMPI%COMM  , RealArray = PartData)
  IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
    CALL DistributedWriteArray(FileName , &
                              DataSetName ='VibQuantData', rank=2            , &
                              nValGlobal  =(/nPart_glob  , MaxQuantNum/)     , &
                              nVal        =(/locnPart    , MaxQuantNum  /)   , &
                              offset      =(/offsetnPart , 0_IK  /)          , &
                              collective  =.FALSE.       , offSetDim=1       , &
                              communicator=PartMPI%COMM  , IntegerArray_i4=VibQuantData)
    DEALLOCATE(VibQuantData)
  END IF
  ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL DistributedWriteArray(FileName , &
                              DataSetName = 'PartTimeStep'  , rank=2      , &
                              nValGlobal  = (/nGlobalElems  , 1_IK/)      , &
                              nVal        = (/PP_nElems     , 1_IK/)      , &
                              offset      = (/offsetElem    , 0_IK/)      , &
                              collective  =.FALSE.          , offSetDim=1 , &
                              communicator=PartMPI%COMM     , RealArray=VarTimeStep%ElemFac)
  END IF
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'PartData'    , rank = 2                , &
                        nValGlobal  = (/nPart_glob  , PartDataSize/)          , &
                        nVal        = (/locnPart    , PartDataSize/)          , &
                        offset      = (/offsetnPart , 0_IK  /)                , &
                        collective  = .TRUE.        , RealArray = PartData)
  IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
    CALL WriteArrayToHDF5(DataSetName = 'VibQuantData' , rank = 2             , &
                          nValGlobal  = (/nPart_glob   , MaxQuantNum/)        , &
                          nVal        = (/locnPart     , MaxQuantNum  /)      , &
                          offset      = (/offsetnPart  , 0_IK /)              , &
                          collective  = .TRUE.         , IntegerArray_i4 = VibQuantData)
    DEALLOCATE(VibQuantData)
  END IF
    ! Output of the element-wise time step as a separate container in state file
  IF(VarTimeStep%UseDistribution) THEN
    CALL WriteArrayToHDF5(DataSetName = 'PartTimeStep'  , rank=2, &
                          nValGlobal  = (/nGlobalElems  , 1_IK/), &
                          nVal        = (/PP_nElems     , 1_IK/)   ,&
                          offset      = (/offsetElem    , 0_IK/)  ,&
                          collective  = .FALSE.         , RealArray=VarTimeStep%ElemFac)
  END IF
  CALL CloseDataFile()
#endif /*USE_MPI*/

END ASSOCIATE
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.


!CALL CloseDataFile()

!  CALL WriteArrayToHDF5('PartData',nPart_glob,2,(/locnPart,PartDataSize/),offsetnPart,1,existing=.FALSE.,RealArray=PartData)!,&
!                        !xfer_mode_independent=.TRUE.)  ! könnte bei Procs die keine Teilchen schreiben
                                                        ! problematisch werden

DEALLOCATE(StrVarNames)
DEALLOCATE(PartInt)
DEALLOCATE(PartData)

!!! Kleiner Hack von JN (Teil 2/2):
useDSMC=withDSMC
IF (.NOT.(useDSMC.OR.PartPressureCell)) THEN
  DEALLOCATE(PEM%pStart , &
             PEM%pNumber, &
             PEM%pNext  , &
             PEM%pEnd   )!, &
             !PDM%nextUsedPosition  )
END IF
!!! Ende kleiner Hack von JN (Teil 2/2)


END SUBROUTINE WriteParticleToHDF5


SUBROUTINE WriteMacroParticleToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the output file on a single processor and writes all the necessary attributes of macroparticles
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_MacroBody_Vars    ,ONLY: UseMacroBody, nMacroBody, MacroSphere
#if USE_MPI
USE MOD_Particle_MPI_Vars ,ONLY: PartMPI
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!===================================================================================================================================
! OUTPUT VARIABLES
!===================================================================================================================================
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
#if USE_MPI
!INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
!INTEGER(KIND=IK)               :: nParticles(0:nProcessors-1)
#endif /*USE_MPI*/
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER(KIND=IK)               :: locnPart_max
!===================================================================================================================================
locnPart =   0_IK
#if USE_MPI
IF (PartMPI%MPIRoot) THEN
#endif /*USE_MPI*/
  IF (UseMacroBody .AND. nMacroBody.GT.0) THEN
    DO pcount = 1,nMacroBody
      !IF(MacroParticle(pcount)%particleLocal) THEN
        locnPart = locnPart + 1_IK
      !END IF
    END DO
  END IF
#if USE_MPI
END IF
#endif /*USE_MPI*/
!#if USE_MPI
!sendbuf(1)=locnPart
!recvbuf=0_IK
!CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
!offsetnPart=recvbuf(1)
!sendbuf(1)=recvbuf(1)+locnPart
!CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!!global numbers
!nPart_glob=sendbuf(1)
!CALL MPI_GATHER(locnPart,1,MPI_INTEGER_INT_KIND,nParticles,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_WORLD,iError)
!LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
!CALL MPI_REDUCE(locnPart, locnPart_max, 1, MPI_INTEGER_INT_KIND, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
!#else
!offsetnPart=0_IK
!nPart_glob=locnPart
!locnPart_max=locnPart
!#endif /*USE_MPI*/
offsetnPart=0_IK
nPart_glob=locnPart
locnPart_max=locnPart
PartDataSize=13
ALLOCATE(PartData(offsetnPart+1_IK:offsetnPart+locnPart,INT(PartDataSize,IK)))

DO iPart=offsetnPart+1_IK,offsetnPart+locnPart
  PartData(iPart,1:3)=MacroSphere(iPart)%center(1:3)
  PartData(iPart,4:9)=MacroSphere(iPart)%velocity(1:6)
  PartData(iPart,10)=MacroSphere(iPart)%radius
  PartData(iPart,11)=MacroSphere(iPart)%temp
  PartData(iPart,12)=MacroSphere(iPart)%density
  PartData(iPart,13)=MacroSphere(iPart)%mass
END DO

reSwitch=.FALSE.
IF(gatheredWrite)THEN
  ! gatheredwrite not working with distributed particles
  ! particles require own routine for which the communicator has to be build each time
  reSwitch=.TRUE.
  gatheredWrite=.FALSE.
END IF

! Associate construct for integer KIND=8 possibility
ASSOCIATE (PartDataSize => INT(PartDataSize,IK))

  ALLOCATE(StrVarNames(PartDataSize))
  StrVarNames(1)='ParticlePositionX'
  StrVarNames(2)='ParticlePositionY'
  StrVarNames(3)='ParticlePositionZ'
  StrVarNames(4)='VelocityX'
  StrVarNames(5)='VelocityY'
  StrVarNames(6)='VelocityZ'
  StrVarNames(7)='RotVelocityX'
  StrVarNames(8)='RotVelocityY'
  StrVarNames(9)='RotVelocityZ'
  StrVarNames(10)='Radius'
  StrVarNames(11)='Temperature'
  StrVarNames(12)='Density'
  StrVarNames(13)='Mass'

#if USE_MPI
  IF(PartMPI%MPIRoot)THEN
#endif /*USE_MPI*/
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesMacroParticles',INT(PartDataSize,4),StrArray=StrVarNames)
    CALL CloseDataFile()
#if USE_MPI
  END IF
#endif /*USE_MPI*/

  IF(locnPart_max.EQ.0)THEN ! zero particles present: write empty dummy container to .h5 file (required for subsequent file access)
#if USE_MPI
    IF(PartMPI%MPIRoot)THEN ! only root writes the container
#endif /*USE_MPI*/
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteArrayToHDF5(DataSetName='MacroPartData'   , rank=2           , &
                            nValGlobal=(/nPart_glob  , INT(PartDataSize,IK)/)   , &
                            nVal=      (/locnPart    , INT(PartDataSize,IK)/) , &
                            offset=    (/offsetnPart , 0_IK  /)         , &
                            collective=.FALSE.       , RealArray=PartData)
      CALL CloseDataFile()
#if USE_MPI
    END IF !MPIRoot
#endif /*USE_MPI*/
  END IF !locnPart_max.EQ.0
#if USE_MPI
  CALL DistributedWriteArray(FileName                     , &
                             DataSetName  = 'MacroPartData'    , rank = 2          , &
                             nValGlobal   = (/nPart_glob  , INT(PartDataSize,IK)/)    , &
                             nVal         = (/locnPart    , INT(PartDataSize,IK)/)    , &
                             offset       = (/offsetnPart , 0_IK/)            , &
                             collective   = .FALSE.       , offSetDim = 1     , &
                             communicator = PartMPI%COMM  , RealArray = PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'MacroPartData'    , rank = 2                , &
                        nValGlobal  = (/nPart_glob  , INT(PartDataSize,IK)/)          , &
                        nVal        = (/locnPart    , INT(PartDataSize,IK)/)          , &
                        offset      = (/offsetnPart , 0_IK  /)                , &
                        collective  = .TRUE.        , RealArray = PartData)
  CALL CloseDataFile()
#endif /*USE_MPI*/
END ASSOCIATE
! reswitch
IF(reSwitch) gatheredWrite=.TRUE.

DEALLOCATE(StrVarNames)
DEALLOCATE(PartData)

END SUBROUTINE WriteMacroParticleToHDF5


SUBROUTINE WriteSurfStateToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the state attributes and data of surface chemistry state and writes it out into State-File
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars              ,ONLY: BC
USE MOD_SurfaceModel_Vars      ,ONLY: SurfDistInfo, Adsorption
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars ,ONLY: SurfCOMM,nSurfBC,SurfBCName
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample,SurfMesh,offSetSurfSide, PartBound, nPartBound
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:), StrVarNamesData(:)
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
CHARACTER(LEN=255)             :: CoordID
#if USE_MPI
INTEGER                        :: sendbuf(1),recvbuf(1)
#endif
INTEGER                        :: locnSurfPart,offsetnSurfPart,nSurfPart_glob
INTEGER                        :: iSpec, nVar, iPB
INTEGER                        :: iOffset, UsedSiteMapPos, SideID, PartBoundID
INTEGER                        :: iSurfSide, isubsurf, jsubsurf, iCoord, nSites, nSitesRemain, iPart, iVar
INTEGER                        :: Coordinations          !number of PartInt and PartData coordinations
INTEGER                        :: SurfPartIntSize        !number of entries in each line of PartInt
INTEGER                        :: SurfPartDataSize       !number of entries in each line of PartData
INTEGER,ALLOCATABLE            :: SurfPartInt(:,:,:,:,:)
INTEGER,ALLOCATABLE            :: SurfPartData(:,:)
REAL,ALLOCATABLE               :: SurfCalcData(:,:,:,:,:)
INTEGER,ALLOCATABLE            :: SurfModelType(:)
LOGICAL                        :: doDistributionData
!===================================================================================================================================
! first check if wallmodel defined and greater than 0 before writing any surface things into state
IF(.NOT.SurfMesh%SurfOnProc) RETURN
IF(.NOT.(ANY(PartBound%Reactive))) RETURN

! only pocs with real surfaces (not halo) in own proc write out
IF(SurfMesh%nMasterSides.EQ.0) RETURN
#if USE_MPI
CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
#endif /*USE_MPI*/

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
#if USE_MPI
IF(SurfCOMM%MPIOutputRoot)THEN
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'Surface_BCs',nSurfBC,StrArray=SurfBCName)
  CALL WriteAttributeToHDF5(File_ID,'nSurfSample',1,IntegerScalar=nSurfSample)
  CALL WriteAttributeToHDF5(File_ID,'nSpecies',1,IntegerScalar=nSpecies)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

! write surfacemodel type for every surface into partstate
ALLOCATE(SurfModelType(SurfMesh%nMasterSides))
SurfModelType = 0
DO iSurfSide = 1,SurfMesh%nMasterSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%Reactive(PartboundID)) THEN
    SurfModelType(iSurfSide) = PartBound%SurfaceModel(PartboundID)
  END IF
END DO

WRITE(H5_Name,'(A)') 'SurfaceModelType'
#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalSides   => INT(SurfMesh%nGlobalSides,IK) ,&
      nSides         => INT(SurfMesh%nMasterSides,IK) ,&
      offsetSurfSide => INT(offsetSurfSide,IK) )
  CALL WriteArrayToHDF5(DataSetName=H5_Name , rank=1, &
                        nValGlobal =(/nGlobalSides/) , &
                        nVal       =(/nSides      /) , &
                        offset     =(/offsetSurfSide/) , &
                        collective =.TRUE.   , IntegerArray_i4=SurfModelType)
END ASSOCIATE
CALL CloseDataFile()
SDEALLOCATE(SurfModelType)

! now write surface calculation data (coverag and temporary adsorption numbers)
! set names and write attributes in hdf5 files
nVar = 4
ALLOCATE(StrVarNames(nVar*nSpecies))
iVar = 1
DO iSpec=1,nSpecies
  WRITE(SpecID,'(I3.3)') iSpec
  StrVarNames(iVar)   = 'Spec'//TRIM(SpecID)//'_Coverage'
  StrVarNames(iVar+1) = 'Spec'//TRIM(SpecID)//'_adsorbnum_tmp'
  StrVarNames(iVar+2) = 'Spec'//TRIM(SpecID)//'_desorbnum_tmp'
  StrVarNames(iVar+3) = 'Spec'//TRIM(SpecID)//'_reactnum_tmp'
  iVar = iVar + nVar
END DO

#if USE_MPI
IF(SurfCOMM%MPIOutputRoot)THEN
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurfCalcData',nVar*nSpecies,StrArray=StrVarNames)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

! rewrite and save arrays for SurfCalcData
ALLOCATE(SurfCalcData(nVar,nSurfSample,nSurfSample,SurfMesh%nMasterSides,nSpecies))
SurfCalcData = 0.
DO iSurfSide = 1,SurfMesh%nMasterSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%Reactive(PartboundID)) THEN
    DO jsubsurf = 1,nSurfSample
      DO isubsurf = 1,nSurfSample
        SurfCalcData(1,iSubSurf,jSubSurf,iSurfSide,:) = Adsorption%Coverage(iSubSurf,jSubSurf,iSurfSide,:)
        IF (PartBound%SurfaceModel(PartBoundID).EQ.3) THEN
          SurfCalcData(2,iSubSurf,jSubSurf,iSurfSide,:) = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%adsorbnum_tmp(:)
          SurfCalcData(3,iSubSurf,jSubSurf,iSurfSide,:) = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%desorbnum_tmp(:)
          SurfCalcData(4,iSubSurf,jSubSurf,iSurfSide,:) = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%reactnum_tmp(:)
        END IF
      END DO
    END DO
  END IF
END DO

WRITE(H5_Name,'(A)') 'SurfCalcData'
!#if USE_MPI
!CALL GatheredWriteArray(FileName,create=.FALSE.,&
!                        DataSetName=H5_Name, rank=5,&
!                        nValGlobal=(/nVar,nSurfSample,nSurfSample,SurfMesh%nGlobalSides,nSpecies/),&
!                        nVal=      (/nVar,nSurfSample,nSurfSample,SurfMesh%nMasterSides,nSpecies/),&
!                        offset=    (/0   ,0          ,0          ,offsetSurfSide       ,0       /),&
!                        collective=.TRUE.,  RealArray=SurfCalcData)
!#else
#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nSpecies       => INT(nSpecies,IK)              ,&
      nVar           => INT(nVar,IK)                  ,&
      nSurfSample    => INT(nSurfSample,IK)           ,&
      nGlobalSides   => INT(SurfMesh%nGlobalSides,IK) ,&
      nSides         => INT(SurfMesh%nMasterSides,IK) ,&
      offsetSurfSide => INT(offsetSurfSide,IK) )
  CALL WriteArrayToHDF5(DataSetName=H5_Name , rank=5                                                   , &
                        nValGlobal =(/nVar   , nSurfSample , nSurfSample , nGlobalSides   , nSpecies/) , &
                        nVal       =(/nVar   , nSurfSample , nSurfSample , nSides         , nSpecies/) , &
                        offset     =(/0_IK   , 0_IK        , 0_IK        , offsetSurfSide , 0_IK    /) , &
                        collective =.TRUE.   , RealArray=SurfCalcData)
END ASSOCIATE
CALL CloseDataFile()
!#endif /*USE_MPI*/
SDEALLOCATE(StrVarNames)
SDEALLOCATE(SurfCalcData)

! next stuff is only relevant if surfacemodel=3 is used
doDistributionData=.FALSE.
DO iPB=1,nPartBound
  IF (PartBound%SurfaceModel(iPB).EQ.3) THEN
    doDistributionData=.TRUE.
    EXIT
  END IF
END DO
IF (.NOT.doDistributionData) RETURN

! save number and positions of binded particles for all coordinations
Coordinations    = 3
SurfPartIntSize  = 3
SurfPartDataSize = 2

locnSurfPart = 0
nSurfPart_glob = 0
offsetnSurfPart = 0

! calculate number of adsorbates on each coordination (already on surface) and all sites
DO iSurfSide = 1,SurfMesh%nMasterSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  IF (PartBound%Reactive(PartboundID)) THEN
    IF (PartBound%SurfaceModel(PartboundID).EQ.3) THEN
      DO jsubsurf = 1,nSurfSample
        DO isubsurf = 1,nSurfSample
          DO iCoord = 1,Coordinations
            nSites = SurfDistInfo(isubsurf,jsubsurf,iSurfSide)%nSites(iCoord)
            nSitesRemain = SurfDistInfo(isubsurf,jsubsurf,iSurfSide)%SitesRemain(iCoord)
            locnSurfPart = locnSurfpart + nSites - nSitesRemain
          END DO
        END DO
      END DO
    END IF
  END IF
END DO

! communicate number of surface particles and offsets in surfpartdata array
#if USE_MPI
sendbuf(1)=locnSurfPart
recvbuf(1)=0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,SurfCOMM%OutputCOMM,iError)
offsetnSurfPart=recvbuf(1)
sendbuf(1)=recvbuf(1)+locnSurfPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,SurfCOMM%nOutputProcs-1,SurfCOMM%OutputCOMM,iError) !last proc knows global number
!global numbers
nSurfPart_glob=sendbuf(1)
#else
offsetnSurfPart=0
nSurfPart_glob=locnSurfPart
#endif

ALLOCATE(SurfPartInt(offsetSurfSide+1:offsetSurfSide+SurfMesh%nMasterSides,nSurfSample,nSurfSample,Coordinations,SurfPartIntSize))
ALLOCATE(SurfPartData(offsetnSurfPart+1:offsetnSurfPart+locnSurfPart,SurfPartDataSize))
iOffset = offsetnSurfPart
DO iSurfSide = 1,SurfMesh%nMasterSides
  SideID = SurfMesh%SurfIDToSideID(iSurfSide)
  PartboundID = PartBound%MapToPartBC(BC(SideID))
  DO jsubsurf = 1,nSurfSample
    DO isubsurf = 1,nSurfSample
      DO iCoord = 1,Coordinations
        IF (PartBound%SurfaceModel(PartboundID).EQ.3) THEN
          nSites = SurfDistInfo(isubsurf,jsubsurf,iSurfSide)%nSites(iCoord)
          nSitesRemain = SurfDistInfo(isubsurf,jsubsurf,iSurfSide)%SitesRemain(iCoord)
          ! set surfpartint array values
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,1) = nSites
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,2) = iOffset
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,3) = iOffset + (nSites - nSitesRemain)
          ! set the surfpartdata array values
          DO iPart = 1, (nSites - nSitesRemain)
            UsedSiteMapPos = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%UsedSiteMap(nSites+1-iPart)
            SurfPartData(iOffset+iPart,1) = UsedSiteMapPos
            SurfPartData(iOffset+iPart,2) = SurfDistInfo(iSubSurf,jSubSurf,iSurfSide)%AdsMap(iCoord)%Species(UsedSiteMapPos)
          END DO
          iOffset = iOffset + nSites - nSitesRemain
        ELSE
          ! set blank surfpartint array values for other surfacemodel
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,1) = 0
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,2) = iOffset
          SurfPartInt(offsetSurfSide+iSurfSide,isubsurf,jsubsurf,iCoord,3) = iOffset
        END IF
      END DO
    END DO
  END DO
END DO

! set names and write attributes in hdf5 files
ALLOCATE(StrVarNames(SurfPartIntSize*Coordinations))
iVar = 1
DO iCoord=1,Coordinations
  WRITE(CoordID,'(I3.3)') iCoord
  StrVarNames(iVar)='Coord'//TRIM(CoordID)//'_nSites'
  StrVarNames(iVar+1)='Coord'//TRIM(CoordID)//'_FirstSurfPartID'
  StrVarNames(iVar+2)='Coord'//TRIM(CoordID)//'_LastSurfPartID'
  iVar = iVar + 3
END DO

ALLOCATE(StrVarNamesData(SurfPartDataSize*Coordinations))
iVar = 1
DO iCoord=1,Coordinations
  WRITE(CoordID,'(I3.3)') iCoord
  StrVarNamesData(iVar)='Coord'//TRIM(CoordID)//'_SiteMapPosition'
  StrVarNamesData(iVar+1)='Coord'//TRIM(CoordID)//'_Species'
  iVar = iVar + 2
END DO
#if USE_MPI
IF(SurfCOMM%MPIOutputRoot)THEN
#endif
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurfPartInt',SurfPartIntSize*Coordinations,StrArray=StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesSurfParticles',SurfPartDataSize*Coordinations,StrArray=StrVarNamesData)
  CALL CloseDataFile()
#if USE_MPI
END IF
CALL MPI_BARRIER(SurfCOMM%OutputCOMM,iERROR)
#endif

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      Coordinations    => INT(Coordinations,IK)        ,&
      SurfPartIntSize  => INT(SurfPartIntSize,IK)      ,&
      nSurfSample      => INT(nSurfSample,IK)          ,&
      nGlobalSides     => INT(SurfMesh%nGlobalSides,IK),&
      nSides           => INT(SurfMesh%nMasterSides,IK),&
      offsetSurfSide   => INT(offsetSurfSide,IK))
#if USE_MPI
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=SurfCOMM%OutputCOMM)
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
  CALL WriteArrayToHDF5(DataSetName = 'SurfPartInt'    , rank = 5                                                      , &
                        nValGlobal  = (/nGlobalSides   , nSurfSample , nSurfSample , Coordinations , SurfPartIntSize/) , &
                        nVal        = (/nSides         , nSurfSample , nSurfSample , Coordinations , SurfPartIntSize/) , &
                        offset      = (/offsetSurfSide , 0_IK        , 0_IK        , 0_IK          , 0_IK           /) , &
                        collective  = .TRUE.                  , IntegerArray_i4 = SurfPartInt(:,:,:,:,:))
  CALL CloseDataFile()
END ASSOCIATE

ASSOCIATE (&
      nSurfPart_glob   => INT(nSurfPart_glob,IK)       ,&
      locnSurfPart     => INT(locnSurfPart,IK)         ,&
      offsetnSurfPart  => INT(offsetnSurfPart,IK)      ,&
      SurfPartDataSize => INT(SurfPartDataSize,IK) )
#if USE_MPI
  CALL DistributedWriteArray(FileName                                              , &
                             DataSetName  = 'SurfPartData'    , rank = 2           , &
                             nValGlobal   = (/nSurfPart_glob  , SurfPartDataSize/) , &
                             nVal         = (/locnSurfPart    , SurfPartDataSize/) , &
                             offset       = (/offsetnSurfPart , 0_IK/)             , &
                             collective   = .FALSE.       , offSetDim = 1          , &
                             communicator = SurfCOMM%OutputCOMM  , IntegerArray_i4 = SurfPartData(:,:))
#else
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName = 'SurfPartData'    , rank = 2                          , &
                        nValGlobal  = (/nSurfPart_glob  , SurfPartDataSize/)                , &
                        nVal        = (/locnSurfPart    , SurfPartDataSize/)                , &
                        offset      = (/offsetnSurfPart , 0_IK            /)                , &
                        collective  = .TRUE.            , IntegerArray_i4 = SurfPartData(: , :))
  CALL CloseDataFile()
#endif
END ASSOCIATE
SDEALLOCATE(StrVarNames)
SDEALLOCATE(StrVarNamesData)
SDEALLOCATE(SurfPartInt)
SDEALLOCATE(SurfPartData)

END SUBROUTINE WriteSurfStateToHDF5


SUBROUTINE WriteAdaptiveInfoToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the adaptive boundary info and writes it out into State-File
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_Particle_Vars          ,ONLY: nSpecies, Adaptive_MacroVal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
INTEGER                        :: nVar
INTEGER                        :: iElem,iVar,iSpec
REAL,ALLOCATABLE               :: AdaptiveData(:,:,:)
!===================================================================================================================================


nVar = 10
iVar = 1
ALLOCATE(StrVarNames(nVar*nSpecies))
DO iSpec=1,nSpecies
  WRITE(SpecID,'(I3.3)') iSpec
  StrVarNames(iVar)   = 'Spec'//TRIM(SpecID)//'-VeloX'
  StrVarNames(iVar+1) = 'Spec'//TRIM(SpecID)//'-VeloY'
  StrVarNames(iVar+2) = 'Spec'//TRIM(SpecID)//'-VeloZ'
  StrVarNames(iVar+3) = 'Spec'//TRIM(SpecID)//'-TempX'
  StrVarNames(iVar+4) = 'Spec'//TRIM(SpecID)//'-TempY'
  StrVarNames(iVar+5) = 'Spec'//TRIM(SpecID)//'-TempZ'
  StrVarNames(iVar+6) = 'Spec'//TRIM(SpecID)//'-Density'
  StrVarNames(iVar+7) = 'Spec'//TRIM(SpecID)//'-PumpVeloPerArea'
  StrVarNames(iVar+8) = 'Spec'//TRIM(SpecID)//'-PumpPressure'
  StrVarNames(iVar+9) = 'Spec'//TRIM(SpecID)//'-PumpIntegralError'
  iVar = iVar + nVar
END DO

IF(MPIRoot)THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdaptive',nVar*nSpecies,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

! rewrite and save arrays for AdaptiveData
ALLOCATE(AdaptiveData(nVar,nSpecies,nElems))
AdaptiveData = 0.
DO iElem = 1,nElems
  AdaptiveData(1,:,iElem) = Adaptive_MacroVal(DSMC_VELOX,iElem,:)
  AdaptiveData(2,:,iElem) = Adaptive_MacroVal(DSMC_VELOY,iElem,:)
  AdaptiveData(3,:,iElem) = Adaptive_MacroVal(DSMC_VELOZ,iElem,:)
  AdaptiveData(4,:,iElem) = Adaptive_MacroVal(DSMC_TEMPX,iElem,:)
  AdaptiveData(5,:,iElem) = Adaptive_MacroVal(DSMC_TEMPY,iElem,:)
  AdaptiveData(6,:,iElem) = Adaptive_MacroVal(DSMC_TEMPZ,iElem,:)
  AdaptiveData(7,:,iElem) = Adaptive_MacroVal(DSMC_NUMDENS,iElem,:)
  ! Porous BC parameter (11: Pumping capacity [m3/s], 12: Static pressure [Pa], 13: Integral pressure difference [Pa])
  AdaptiveData(8:10,:,iElem) = Adaptive_MacroVal(11:13,iElem,:)
END DO

WRITE(H5_Name,'(A)') 'AdaptiveInfo'
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

! Associate construct for integer KIND=8 possibility
ASSOCIATE (&
      nGlobalElems    => INT(nGlobalElems,IK)    ,&
      nElems          => INT(nElems,IK)          ,&
      nVar            => INT(nVar,IK)            ,&
      nSpecies        => INT(nSpecies,IK)        ,&
      offsetElem      => INT(offsetElem,IK)      )
  CALL WriteArrayToHDF5(DataSetName = H5_Name , rank = 3                   , &
                        nValGlobal  = (/nVar  , nSpecies  , nGlobalElems/) , &
                        nVal        = (/nVar  , nSpecies  , nElems      /) , &
                        offset      = (/0_IK  , 0_IK      , offsetElem  /) , &
                        collective  = .false. , RealArray = AdaptiveData)
END ASSOCIATE
CALL CloseDataFile()
SDEALLOCATE(StrVarNames)
SDEALLOCATE(AdaptiveData)

END SUBROUTINE WriteAdaptiveInfoToHDF5


SUBROUTINE WriteVibProbInfoToHDF5(FileName)
!===================================================================================================================================
!> Subroutine that generates the adaptive boundary info and writes it out into State-File
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nGlobalElems, nElems
USE MOD_Particle_Vars          ,ONLY: nSpecies
USE MOD_DSMC_Vars              ,ONLY: VarVibRelaxProb, CollisMode, DSMC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: H5_Name
CHARACTER(LEN=255)             :: SpecID
INTEGER                        :: iSpec
!===================================================================================================================================
IF(CollisMode.GT.1) THEN
  IF(DSMC%VibRelaxProb.GE.2.0) THEN
    ALLOCATE(StrVarNames(nSpecies))
    DO iSpec=1,nSpecies
      WRITE(SpecID,'(I3.3)') iSpec
      StrVarNames(iSpec)   = 'Spec'//TRIM(SpecID)//'-VibProbRelax'
    END DO

    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesVibProbInfo',nSpecies,StrArray=StrVarNames)
      CALL CloseDataFile()
    END IF

    WRITE(H5_Name,'(A)') 'VibProbInfo'
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)

    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nGlobalElems    => INT(nGlobalElems,IK)    ,&
          nElems          => INT(nElems,IK)          ,&
          nSpecies        => INT(nSpecies,IK)        ,&
          offsetElem      => INT(offsetElem,IK)      )
      CALL WriteArrayToHDF5(DataSetName = H5_Name , rank = 2        , &
                            nValGlobal  = (/nGlobalElems, nSpecies/) , &
                            nVal        = (/nElems      , nSpecies/) , &
                            offset      = (/offsetElem  ,0_IK     /) , &
                            collective  = .false. , RealArray = VarVibRelaxProb%ProbVibAv)
    END ASSOCIATE
    CALL CloseDataFile()
    SDEALLOCATE(StrVarNames)
  ELSE ! DSMC%VibRelaxProb < 2.0
#if USE_MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
#else
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
    CALL WriteAttributeToHDF5(File_ID,'VibProbConstInfo',1,RealScalar=DSMC%VibRelaxProb)
    CALL CloseDataFile()
  END IF
ELSE ! CollisMode <= 1
  IF(MPIRoot)THEN
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'VibProbConstInfo',1,RealScalar=0.)
    CALL CloseDataFile()
  END IF
END IF

END SUBROUTINE WriteVibProbInfoToHDF5


SUBROUTINE WriteClonesToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY: offsetElem
USE MOD_DSMC_Vars     ,ONLY: UseDSMC, CollisMode, DSMC, PolyatomMolDSMC, SpecDSMC
USE MOD_DSMC_Vars     ,ONLY: RadialWeighting, ClonedParticles
USE MOD_PARTICLE_Vars ,ONLY: nSpecies, usevMPF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
!INTEGER(HID_T)                 :: Dset_ID
!INTEGER                        :: nVal
#if USE_MPI
INTEGER(KIND=IK)               :: sendbuf(2),recvbuf(2)
#endif
INTEGER                        :: pcount, iDelay, iElem_glob
LOGICAL                        :: withDSMC=.FALSE.
INTEGER(KIND=IK)               :: locnPart,offsetnPart
INTEGER(KIND=IK)               :: iPart,nPart_glob
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER, ALLOCATABLE           :: VibQuantData(:,:)
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
INTEGER                        :: MaxQuantNum, iPolyatMole, iSpec, tempDelay
!-----------------------------------------------------------------------------------------------------------------------------------
!!added for Evib, Erot writeout
withDSMC=useDSMC
IF (withDSMC) THEN
  IF ((CollisMode.GT.1).AND.(usevMPF) .AND. DSMC%ElectronicModel ) THEN !int ener + 3, vmpf +1
    PartDataSize=13
  ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. DSMC%ElectronicModel ) ) THEN !int ener + 2 and vmpf + 1
                                                                            ! or int energ +3 but no vmpf +1
    PartDataSize=12
  ELSE IF (CollisMode.GT.1) THEN
    PartDataSize=11 !int ener + 2
  ELSE IF (usevMPF) THEN
    PartDataSize=10 !+ 1 vmpf
  ELSE
    PartDataSize=9 !+ 0
  END IF
ELSE IF (usevMPF) THEN
  PartDataSize=10 !vmpf +1
ELSE
  PartDataSize=9
END IF

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  MaxQuantNum = 0
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF (PolyatomMolDSMC(iPolyatMole)%VibDOF.GT.MaxQuantNum) MaxQuantNum = PolyatomMolDSMC(iPolyatMole)%VibDOF
    END IF
  END DO
END IF

locnPart =   0

SELECT CASE(RadialWeighting%CloneMode)
CASE(1)
  tempDelay = RadialWeighting%CloneInputDelay - 1
CASE(2)
  tempDelay = RadialWeighting%CloneInputDelay
CASE DEFAULT
  CALL abort(__STAMP__,&
              'RadialWeighting: CloneMode is not supported!')
END SELECT

DO pcount = 0,tempDelay
    locnPart = locnPart + RadialWeighting%ClonePartNum(pcount)
END DO

#if USE_MPI
sendbuf(1)=locnPart
recvbuf=0
CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER_INT_KIND,MPI_SUM,MPI_COMM_WORLD,iError)
offsetnPart=recvbuf(1)
sendbuf(1)=recvbuf(1)+locnPart
CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER_INT_KIND,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
!global numbers
nPart_glob=sendbuf(1)
#else
offsetnPart=0
nPart_glob=locnPart
#endif
ALLOCATE(PartData(offsetnPart+1:offsetnPart+locnPart,PartDataSize))

IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  ALLOCATE(VibQuantData(offsetnPart+1:offsetnPart+locnPart,MaxQuantNum))
  !+1 is real number of necessary vib quants for the particle
END IF
iPart=offsetnPart
DO iDelay=0,tempDelay
  DO pcount = 1, RadialWeighting%ClonePartNum(iDelay)
    iElem_glob = ClonedParticles(pcount,iDelay)%Element + offsetElem
    iPart = iPart + 1
    PartData(iPart,1)=ClonedParticles(pcount,iDelay)%PartState(1)
    PartData(iPart,2)=ClonedParticles(pcount,iDelay)%PartState(2)
    PartData(iPart,3)=ClonedParticles(pcount,iDelay)%PartState(3)
    PartData(iPart,4)=ClonedParticles(pcount,iDelay)%PartState(4)
    PartData(iPart,5)=ClonedParticles(pcount,iDelay)%PartState(5)
    PartData(iPart,6)=ClonedParticles(pcount,iDelay)%PartState(6)
    PartData(iPart,7)=REAL(ClonedParticles(pcount,iDelay)%Species)
    PartData(iPart,8)=REAL(iElem_glob)
    PartData(iPart,9)=REAL(iDelay)
    IF (withDSMC) THEN
      IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicModel) ) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(iPart,11)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(iPart,12)=ClonedParticles(pcount,iDelay)%PartStateIntEn(3)
        PartData(iPart,13)=ClonedParticles(pcount,iDelay)%WeightingFactor
      ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(iPart,11)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(iPart,12)=ClonedParticles(pcount,iDelay)%WeightingFactor
      ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(iPart,11)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
        PartData(iPart,12)=ClonedParticles(pcount,iDelay)%PartStateIntEn(3)
      ELSE IF (CollisMode.GT.1) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%PartStateIntEn(1)
        PartData(iPart,11)=ClonedParticles(pcount,iDelay)%PartStateIntEn(2)
      ELSE IF (usevMPF) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%WeightingFactor
      END IF
    ELSE IF (usevMPF) THEN
        PartData(iPart,10)=ClonedParticles(pcount,iDelay)%WeightingFactor
    END IF
    IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
      IF (SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(ClonedParticles(pcount,iDelay)%Species)%SpecToPolyArray
        VibQuantData(iPart,1:PolyatomMolDSMC(iPolyatMole)%VibDOF) = &
          ClonedParticles(pcount,iDelay)%VibQuants(1:PolyatomMolDSMC(iPolyatMole)%VibDOF)
      ELSE
          VibQuantData(iPart,:) = 0
      END IF
    END IF
  END DO
END DO

ALLOCATE(StrVarNames(PartDataSize))
StrVarNames(1)='ParticlePositionX'
StrVarNames(2)='ParticlePositionY'
StrVarNames(3)='ParticlePositionZ'
StrVarNames(4)='VelocityX'
StrVarNames(5)='VelocityY'
StrVarNames(6)='VelocityZ'
StrVarNames(7)='Species'
StrVarNames(8)='Element'
StrVarNames(9)='CloneDelay'

IF(withDSMC)THEN
  IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicModel))THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='Electronic'
    StrVarNames(13)='MPF'
  ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='MPF'
  ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicModel) ) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
    StrVarNames(12)='Electronic'
  ELSE IF (CollisMode.GT.1) THEN
    StrVarNames(10)='Vibrational'
    StrVarNames(11)='Rotational'
  ELSE IF (usevMPF) THEN
    StrVarNames(10)='MPF'
  END IF
END IF

#if USE_MPI
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
#else
CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
#endif
CALL WriteAttributeToHDF5(File_ID,'VarNamesParticleClones',PartDataSize,StrArray=StrVarNames)

ASSOCIATE (&
      nPart_glob    => INT(nPart_glob,IK)    ,&
      offsetnPart      => INT(offsetnPart,IK)      ,&
      MaxQuantNum     => INT(MaxQuantNum,IK)     ,&
      PartDataSize    => INT(PartDataSize,IK)    )
CALL WriteArrayToHDF5(DataSetName='CloneData', rank=2,&
                      nValGlobal=(/nPart_glob,PartDataSize/),&
                      nVal=      (/locnPart,PartDataSize  /),&
                      offset=    (/offsetnPart , 0_IK  /),&
                      collective=.FALSE., RealArray=PartData)
IF (withDSMC.AND.(DSMC%NumPolyatomMolecs.GT.0)) THEN
  CALL WriteArrayToHDF5(DataSetName='CloneVibQuantData', rank=2,&
                        nValGlobal=(/nPart_glob,MaxQuantNum/),&
                        nVal=      (/locnPart,MaxQuantNum  /),&
                        offset=    (/offsetnPart , 0_IK  /),&
                        collective=.FALSE., IntegerArray_i4=VibQuantData)
  DEALLOCATE(VibQuantData)
END IF
END ASSOCIATE

CALL CloseDataFile()

DEALLOCATE(StrVarNames)
DEALLOCATE(PartData)

END SUBROUTINE WriteClonesToHDF5
#endif /*PARTICLES*/


SUBROUTINE WriteTimeAverage(MeshFileName,OutputTime,PreviousTime,VarNamesAvg,VarNamesFluc,UAvg,UFluc,dtAvg,nVar_Avg,nVar_Fluc)
!==================================================================================================================================
!> Subroutine to write time averaged data and fluctuations HDF5 format
!==================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars ,ONLY: ProjectName
USE MOD_Mesh_Vars    ,ONLY: offsetElem,nGlobalElems,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName                                 !< Name of mesh file
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesAvg(nVar_Avg)                        !< Average variable names
CHARACTER(LEN=*),INTENT(IN)    :: VarNamesFluc(nVar_Fluc)                      !< Fluctuations variable names
REAL,INTENT(IN)                :: OutputTime                                   !< Time of output
REAL,INTENT(IN),OPTIONAL       :: PreviousTime                                 !< Time of previous output
REAL,INTENT(IN),TARGET         :: UAvg(nVar_Avg,0:PP_N,0:PP_N,0:PP_N,nElems)   !< Averaged Solution
REAL,INTENT(IN),TARGET         :: UFluc(nVar_Fluc,0:PP_N,0:PP_N,0:PP_N,nElems) !< Fluctuations
REAL,INTENT(IN)                :: dtAvg                                        !< Timestep of averaging
INTEGER,INTENT(IN)             :: nVar_Avg                                     !< Number of averaged variables
INTEGER,INTENT(IN)             :: nVar_Fluc                                    !< Number of fluctuations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
REAL                           :: StartT,EndT
!==================================================================================================================================
IF((nVar_Avg.EQ.0).AND.(nVar_Fluc.EQ.0)) RETURN ! no time averaging
StartT=PICLASTIME()
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE TIME AVERAGED STATE AND FLUCTUATIONS TO HDF5 FILE...'
END IF

! generate nextfile info in previous output file
IF(PRESENT(PreviousTime))THEN
  IF(MPIRoot .AND. PreviousTime.LT.OutputTime) CALL GenerateNextFileInfo('TimeAvg',OutputTime,PreviousTime)
END IF

! Write timeaverages ---------------------------------------------------------------------------------------------------------------
IF(nVar_Avg.GT.0)THEN
  ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_TimeAvg',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton('TimeAvg',nVar_Avg,VarNamesAvg,MeshFileName,OutputTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

  ! Reopen file and write DG solution
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        nElems          => INT(nElems,IK)          ,&
        nVar_Avg        => INT(nVar_Avg,IK)        ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName , create = .FALSE.  , &
                            DataSetName     = 'DG_Solution' , rank = 5                                        , &
                            nValGlobal      = (/nVar_Avg , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                            nVal            = (/nVar_Avg , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nElems/)       , &
                            offset          = (/0_IK     , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
                            collective      = .TRUE.        , RealArray = UAvg)
  END ASSOCIATE
END IF

! Write fluctuations ---------------------------------------------------------------------------------------------------------------
IF(nVar_Fluc.GT.0)THEN
  FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_Fluc',OutputTime))//'.h5'
  IF(MPIRoot)THEN
    CALL GenerateFileSkeleton('Fluc',nVar_Fluc,VarNamesFluc,MeshFileName,OutputTime)
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttributeToHDF5(File_ID,'AvgTime',1,RealScalar=dtAvg)
    CALL CloseDataFile()
  END IF
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*USE_MPI*/

  ! Reopen file and write DG solution
  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        nElems          => INT(nElems,IK)          ,&
        nVar_Fluc       => INT(nVar_Fluc,IK)       ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName , create = .FALSE.                                                          , &
                            DataSetName     = 'DG_Solution' , rank = 5                                           , &
                            nValGlobal      = (/nVar_Fluc   , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
                            nVal            = (/nVar_Fluc   , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nElems/)       , &
                            offset          = (/0_IK        , 0_IK      , 0_IK      , 0_IK      , offsetElem/)   , &
                            collective      = .TRUE.        , RealArray = UFluc)
  END ASSOCIATE
END IF

endT=PICLASTIME()
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
END SUBROUTINE WriteTimeAverage


SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Output_Vars        ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems
USE MOD_Interpolation_Vars ,ONLY: NodeType
#ifdef INTEL
USE IFPORT                 ,ONLY: SYSTEM
#endif
!USE MOD_PreProcFlags
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
INTEGER,INTENT(IN)             :: nVar
!INTEGER,INTENT(IN)             :: NData
CHARACTER(LEN=255)             :: StrVarNames(nVar)
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(5)
CHARACTER(LEN=255)             :: FileName,MeshFile255
!CHARACTER(LEN=255),ALLOCATABLE :: params(:)
!===================================================================================================================================
! Create file
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)

! Write file header
CALL WriteHDF5Header(TRIM(TypeString),File_ID)

! Preallocate the data space for the dataset.
Dimsf=(/nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/)
CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE
CALL H5DCREATE_F(File_ID,'DG_Solution', HDF5DataType, FileSpace, DSet_ID, iError)
! Close the filespace and the dataset
CALL H5DCLOSE_F(Dset_id, iError)
CALL H5SCLOSE_F(FileSpace, iError)

! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=N)
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFileName)/))
IF(PRESENT(FutureTime))THEN
  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
END IF
CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)

CALL WriteAttributeToHDF5(File_ID,'NComputation',1,IntegerScalar=PP_N)

CALL CloseDataFile()

! Add userblock to hdf5-file
CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)

END SUBROUTINE GenerateFileSkeleton


SUBROUTINE GenerateNextFileInfo(TypeString,OutputTime,PreviousTime)
!===================================================================================================================================
!> Subroutine that opens the prvious written file on root processor and writes the necessary nextfile info
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_Interpolation_Vars ,ONLY: NodeType
#ifdef INTEL
USE IFPORT                 ,ONLY: SYSTEM
#endif
!USE MOD_PreProcFlags
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN)                :: PreviousTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName,MeshFile255
!===================================================================================================================================
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),PreviousTime))//'.h5'
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=(/MeshFile255/))
CALL CloseDataFile()

END SUBROUTINE GenerateNextFileInfo


SUBROUTINE FlushHDF5(FlushTime_In)
!===================================================================================================================================
! Deletes all HDF5 output files, beginning from time Flushtime
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars     ,ONLY: ProjectName
USE MOD_HDF5_Input       ,ONLY: GetHDF5NextFileName
#if USE_LOADBALANCE
USE MOD_Loadbalance_Vars ,ONLY: DoLoadBalance,nLoadBalance
#endif /*USE_LOADBALANCE*/
#if USE_QDS_DG
USE MOD_QDS_DG_Vars      ,ONLY: DoQDS
#endif /*USE_QDS_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN),OPTIONAL :: FlushTime_In
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: stat,ioUnit
REAL                     :: FlushTime
CHARACTER(LEN=255)       :: InputFile,NextFile
!===================================================================================================================================
IF(.NOT.MPIRoot) RETURN

#if USE_LOADBALANCE
IF(DoLoadBalance.AND.nLoadBalance.GT.0) RETURN
#endif /*USE_LOADBALANCE*/

WRITE(UNIT_stdOut,'(a)')' DELETING OLD HDF5 FILES...'
IF (.NOT.PRESENT(FlushTime_In)) THEN
  FlushTime=0.0
ELSE
  FlushTime=FlushTime_In
END IF

! delete state files
NextFile=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',FlushTime))//'.h5'
DO
  InputFile=TRIM(NextFile)
  ! Read calculation time from file
#if USE_MPI
  CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
  CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
  ! Delete File - only root
  stat=0
  OPEN ( NEWUNIT= ioUnit,         &
         FILE   = InputFile,      &
         STATUS = 'OLD',          &
         ACTION = 'WRITE',        &
         ACCESS = 'SEQUENTIAL',   &
         IOSTAT = stat          )
  IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
  IF(iError.NE.0) EXIT  ! iError is set in GetHDF5NextFileName !
END DO

#if USE_QDS_DG
! delete QDS state files
IF(DoQDS)THEN
  NextFile=TRIM(TIMESTAMP(TRIM(ProjectName)//'_QDS',FlushTime))//'.h5'
  DO
    InputFile=TRIM(NextFile)
    ! Read calculation time from file
#if USE_MPI
    CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
    CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
    ! Delete File - only root
    stat=0
    OPEN ( NEWUNIT= ioUnit,         &
           FILE   = InputFile,      &
           STATUS = 'OLD',          &
           ACTION = 'WRITE',        &
           ACCESS = 'SEQUENTIAL',   &
           IOSTAT = stat          )
    IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
    IF(iError.NE.0) EXIT  ! iError is set in GetHDF5NextFileName !
  END DO
END IF
#endif /*USE_QDS_DG*/

WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'

END SUBROUTINE FlushHDF5



SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars ,ONLY: ProgramName,FileVersion,ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: FileType_in
INTEGER(HID_T),INTENT(IN)                :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                       :: tmp255
!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number

!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files

! First write program name
tmp255=TRIM(ProgramName)
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=(/tmp255/))
tmp255=TRIM(FileType_in)
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=(/tmp255/))
tmp255=TRIM(ProjectName)
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=(/tmp255/))
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header


SUBROUTINE WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,offset,&
                            collective,resizeDim,chunkSize,&
                            RealArray,IntegerArray,StrArray,IntegerArray_i4)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: DataSetName
INTEGER,INTENT(IN)                            :: rank                  ! < number of dimensions of the array
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)      ! < max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)            ! < size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)          ! < offset =0, start at beginning of the array
LOGICAL,INTENT(IN)                            :: collective            ! < use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL                   :: resizeDim(rank)       ! < specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL                   :: chunkSize(rank)       ! < specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(rank)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray_i4(rank) ! KIND=4
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray(rank)    ! KIND=IK
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(rank)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,Type_ID,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet=255
LOGICAL                        :: chunky
#ifndef HDF5_F90 /* HDF5 compiled with fortran2003 flag */
TYPE(C_PTR)                    :: buf
#endif
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! specify chunk size if desired
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF
! make array extendable in case you want to append something
IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize))&
    CALL abort(&
    __STAMP__&
    ,'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray))        Type_ID=H5T_NATIVE_DOUBLE
!IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntegerArray))     Type_ID=h5kind_to_type(IK,H5_INTEGER_KIND)
IF(PRESENT(IntegerArray_i4))  Type_ID=h5kind_to_type(4,H5_INTEGER_KIND)
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(Type_ID, SizeSet, iError)
END IF

Dimsf = nValGlobal ! we need the global array size
CALL H5ESET_AUTO_F(0,iError)
CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
IF(iError.NE.0)THEN ! does not exist
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), Type_ID, FileSpace, DSet_ID,iError,dsetparams)
  CALL H5SCLOSE_F(FileSpace, iError)
END IF
CALL H5ESET_AUTO_F(1,iError)
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
OffsetHDF = Offset
! Create the data space in the memory
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
ELSE
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SSELECT_NONE_F(FileSpace,iError)
ELSE
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, OffsetHDF, Dimsf, iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#if USE_MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

!Write the dataset collectively.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(IntegerArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,IntegerArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,RealArray   ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,StrArray    ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#else
IF(PRESENT(IntegerArray))    buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerArray_i4)) buf=C_LOC(IntegerArray_i4)
IF(PRESENT(RealArray))       buf=C_LOC(RealArray)
IF(PRESENT(StrArray))        buf=C_LOC(StrArray(1))
!IF(ANY(Dimsf.EQ.0)) buf =NULL()
IF(ANY(Dimsf.EQ.0)) THEN
  CALL H5DWRITE_F(DSet_ID,Type_ID,C_NULL_PTR,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
ELSE
  CALL H5DWRITE_F(DSet_ID,Type_ID,buf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
#endif /* HDF5_F90 */

IF(PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close the property list, dataspaces and dataset.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArrayToHDF5


SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,&
                                RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)           :: Loc_ID_in
CHARACTER(LEN=*)  ,INTENT(IN)           :: AttribName
INTEGER           ,INTENT(IN)           :: nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerScalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray(nVal)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray(nVal)
LOGICAL           ,INTENT(IN),OPTIONAL        :: LogicalScalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,Type_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER,TARGET                 :: logtoint
#ifndef HDF5_F90 /* HDF5 compiled with fortran2003 flag */
TYPE(C_PTR)                    :: buf
#endif
!===================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
IF(PRESENT(RealScalar))    Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerScalar)) Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar))THEN
  LogToInt=MERGE(1,0,LogicalScalar)
  Type_ID=H5T_NATIVE_INTEGER
END IF
IF(PRESENT(StrScalar).OR.PRESENT(StrArray))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  IF(PRESENT(StrScalar))THEN
    AttrLen=LEN(StrScalar(1))
  ELSE
    AttrLen=255
  END IF
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)
! Write the attribute data.
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
IF(PRESENT(RealArray))     CALL H5AWRITE_F(Attr_ID, Type_ID, RealArray,     Dimsf, iError)
IF(PRESENT(RealScalar))    CALL H5AWRITE_F(Attr_ID, Type_ID, RealScalar,    Dimsf, iError)
IF(PRESENT(IntegerArray))  CALL H5AWRITE_F(Attr_ID, Type_ID, IntegerArray,  Dimsf, iError)
IF(PRESENT(IntegerScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, IntegerScalar, Dimsf, iError)
IF(PRESENT(LogicalScalar)) CALL H5AWRITE_F(Attr_ID, Type_ID, LogToInt,      Dimsf, iError)
IF(PRESENT(StrScalar))     CALL H5AWRITE_F(Attr_ID, Type_ID, StrScalar,     Dimsf, iError)
IF(PRESENT(StrArray))      CALL H5AWRITE_F(Attr_ID, Type_ID, StrArray,      Dimsf, iError)
#else /* HDF5_F90 */
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntegerArray))  buf=C_LOC(IntegerArray)
IF(PRESENT(IntegerScalar)) buf=C_LOC(IntegerScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(LogToInt)
IF(PRESENT(StrScalar))     buf=C_LOC(StrScalar(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))
CALL H5AWRITE_F(Attr_ID, Type_ID, buf, iError)
#endif /* HDF5_F90 */

! Close datatype
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttributeToHDF5


SUBROUTINE GatheredWriteArray(FileName,create,DataSetName,rank,nValGlobal,nVal,offset,collective,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: FileName,DataSetName
LOGICAL,INTENT(IN)                            :: create,collective
INTEGER,INTENT(IN)                            :: rank
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)               ! max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)                     ! size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)                   ! offset =0, start at beginning of the array
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
REAL,              ALLOCATABLE          :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE          :: UStr(:)
INTEGER(KIND=IK),ALLOCATABLE            :: UInt(:)
INTEGER(KIND=IK)                        :: nValGather(rank),nDOFLocal
INTEGER(KIND=IK),DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
INTEGER(KIND=IK)                        :: i
!===================================================================================================================================
IF(gatheredWrite)THEN
  IF(ANY(offset(1:rank-1).NE.0)) &
    CALL abort(&
    __STAMP__&
    ,'Offset only allowed in last dimension for gathered IO.')

  ! Get last dim of each array on IO nodes
  nDOFLocal=PRODUCT(nVal)
  CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER_INT_KIND,nDOFPerNode,1,MPI_INTEGER_INT_KIND,0,MPI_COMM_NODE,iError)

  ! Allocate big array and compute offsets of small arrs inside big
  offsetNode=0_IK
  IF(MPILocalRoot)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nLocalProcs
      offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    IF(PRESENT(RealArray))     ALLOCATE(UReal(PRODUCT(nValGather)))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( PRODUCT(nValGather)))
    IF(PRESENT(StrArray))      ALLOCATE(UStr( PRODUCT(nValGather)))
  ELSE
    IF(PRESENT(RealArray))     ALLOCATE(UReal(1))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( 1))
    IF(PRESENT(StrArray))      ALLOCATE(UStr( 1))
  ENDIF

  ! Associate construct for integer settings
  ASSOCIATE (&
        nDOFLocal    => INT(nDOFLocal)   ,&
        nDOFPerNode  => INT(nDOFPerNode) ,&
        offsetNode   => INT(offsetNode)  )
    ! Gather small arrays on IO nodes
    IF(PRESENT(RealArray)) CALL MPI_GATHERV(&
        RealArray , nDOFLocal     ,              MPI_DOUBLE_PRECISION , &
        UReal     , nDOFPerNode   , offsetNode , MPI_DOUBLE_PRECISION , &
        0         , MPI_COMM_NODE , iError)
    IF(PRESENT(IntegerArray))  CALL MPI_GATHERV(&
        IntegerArray  , nDOFLocal     ,              MPI_INTEGER_INT_KIND , &
        UInt          , nDOFPerNode   , offsetNode , MPI_INTEGER_INT_KIND , &
        0             , MPI_COMM_NODE , iError)
    !IF(PRESENT(StrArray))  CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE,&
    !                                        UReal,nDOFPerNode, offsetNode,MPI_DOUBLE,0,MPI_COMM_NODE,iError)
  END ASSOCIATE

  IF(MPILocalRoot)THEN
    ! Reopen file and write DG solution (only IO nodes)
    CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_LEADERS)
    IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank                  , nValGlobal , nValGather , &
                                                 offset      , collective=collective , RealArray=UReal)
    IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank                  , nValGlobal , nValGather , &
                                                     offset      , collective=collective , IntegerArray =UInt)
    !IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nValGather,&
    !                                             offset,collective=collective,StrArr =UStr)
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(UReal)
  SDEALLOCATE(UInt)
  SDEALLOCATE(UStr)
ELSE
#endif
  CALL OpenDataFile(FileName,create=create,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                               offset      , collective , RealArray=RealArray)
  IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset          , collective , IntegerArray =IntegerArray)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                               offset      , collective , StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE GatheredWriteArray

#ifdef PARTICLES
#if USE_MPI
SUBROUTINE DistributedWriteArray(FileName,DataSetName,rank,nValGlobal,nVal,offset,collective,&
                                 offSetDim,communicator,RealArray,IntegerArray,StrArray,IntegerArray_i4)
!===================================================================================================================================
!> Write distributed data, that is not present in each proc of given communicator
!>   e.g. master surfaces that are not hosted by each proc
!> 1: check if every proc of given communicator has data
!> 2: if any proc has no data, split the communicator and write only with the new communicator
!> 3: else write with all procs of the given communicator
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                   :: FileName,DataSetName
LOGICAL,INTENT(IN)                            :: collective
INTEGER,INTENT(IN)                            :: offSetDim,communicator
INTEGER,INTENT(IN)                            :: rank
INTEGER(KIND=IK),INTENT(IN)                   :: nValGlobal(rank)               ! max size of array in offset dimension
INTEGER(KIND=IK),INTENT(IN)                   :: nVal(rank)                     ! size of complete (local) array to write
INTEGER(KIND=IK),INTENT(IN)                   :: offset(rank)                   ! offset =0, start at beginning of the array
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER(KIND=IK),INTENT(IN),OPTIONAL,TARGET   :: IntegerArray( PRODUCT(nVal))   ! KIND=IK
INTEGER         ,INTENT(IN),OPTIONAL,TARGET   :: IntegerArray_i4( PRODUCT(nVal))! KIND=4
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER                        :: Color, OutPutCOMM,nOutPutProcs,MyOutputRank
LOGICAL                        :: DataOnProc, DoNotSplit
!===================================================================================================================================

DataOnProc=.FALSE.
! 1: check if every proc of given communicator has data
IF(nVal(offSetDim).GT.0) DataOnProc=.TRUE.
CALL MPI_ALLREDUCE(DataOnProc,DoNotSplit, 1, MPI_LOGICAL, MPI_LAND, COMMUNICATOR, IERROR)

IF(.NOT.DoNotSplit)THEN
! 2: if any proc has no data, split the communicator and write only with the new communicator
  color=MPI_UNDEFINED
  IF(DataOnProc) color=2001
  MyOutputRank=0

  CALL MPI_COMM_SPLIT(COMMUNICATOR, color, MyOutputRank, OutputCOMM,iError)
  IF(DataOnProc) THEN
    CALL MPI_COMM_SIZE(OutputCOMM, nOutPutProcs,iError)
    IF(nOutPutProcs.EQ.1)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal           , nVal , &
                                                   offset      , collective=.FALSE. , RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal                  , nVal , &
                                                   offset          , collective=.FALSE. , IntegerArray =IntegerArray)
      IF(PRESENT(IntegerArray_i4))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal                  , nVal , &
                                                   offset          , collective=.FALSE. , IntegerArray_i4 =IntegerArray_i4)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank               , nValGlobal          , nVal , &
                                                   offset      , collective=.FALSE. , StrArray =StrArray)
      CALL CloseDataFile()
    ELSE
      CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                                   offset      , collective , RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                                   offset          , collective , IntegerArray =IntegerArray)
      IF(PRESENT(IntegerArray_i4)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                                   offset          , collective , IntegerArray_i4 =IntegerArray_i4)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                                   offset      , collective , StrArray =StrArray)
      CALL CloseDataFile()
    END IF
    CALL MPI_BARRIER(OutputCOMM,IERROR)
    CALL MPI_COMM_FREE(OutputCOMM,iERROR)
  END IF
  ! MPI Barrier is required, that the other procs don't open the datafile while this procs are still writing
  CALL MPI_BARRIER(COMMUNICATOR,IERROR)
  OutputCOMM=MPI_UNDEFINED
ELSE
#endif
! 3: else write with all procs of the given communicator
  ! communicator_opt has to be the given communicator or else procs that are not in the given communicator might block the write out
  ! e.g. surface communicator contains only procs with physical surface and MPI_COMM_WORLD contains every proc
  !      Consequently, MPI_COMM_WORLD would block communication
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=communicator)
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal           , nVal , &
                                               offset      , collective , RealArray=RealArray)
  IF(PRESENT(IntegerArray)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset         , collective , IntegerArray =IntegerArray)
  IF(PRESENT(IntegerArray_i4)) CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal                  , nVal , &
                                               offset         , collective , IntegerArray_i4 =IntegerArray_i4)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName , rank       , nValGlobal          , nVal , &
                                               offset      , collective , StrArray =StrArray)
  CALL CloseDataFile()
#if USE_MPI
END IF
#endif

END SUBROUTINE DistributedWriteArray
#endif /*USE_MPI*/


SUBROUTINE WriteNodeSourceExtToHDF5(OutputTime)
!===================================================================================================================================
! Write NodeSourceExt (external charge density) field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars    ,ONLY: NodeSourceExtGlobal
USE MOD_Mesh_Vars          ,ONLY: MeshFile,nGlobalElems,offsetElem,Vdm_EQ_N
USE MOD_Globals_Vars       ,ONLY: ProgramName,FileVersion,ProjectName
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt,CellLocNodes_Volumes,NodeSourceExtTmp
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
#if USE_MPI
USE MOD_Particle_MPI       ,ONLY: AddHaloNodeData
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: N_variables=1
CHARACTER(LEN=255),ALLOCATABLE  :: StrVarNames(:)
CHARACTER(LEN=255)  :: FileName,DataSetName
INTEGER             :: iElem,i
REAL                :: NodeSourceExtEqui(1:N_variables,0:1,0:1,0:1)
!===================================================================================================================================
! create global Eps field for parallel output of Eps distribution
ALLOCATE(NodeSourceExtGlobal(1:N_variables,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems))
ALLOCATE(StrVarNames(1:N_variables))
StrVarNames(1)='NodeSourceExt'
NodeSourceExtGlobal=0.

! Communicate the NodeSourceExtTmp values of the last boundary interaction before the state is written to .h5
#if USE_MPI
CALL AddHaloNodeData(NodeSourceExtTmp)
#endif /*USE_MPI*/

! Add NodeSourceExtTmp values of the last boundary interaction
NodeSourceExt    = NodeSourceExt + NodeSourceExtTmp
NodeSourceExtTmp = 0.

! Loop over all elements and store charge density values in equidistantly distributed nodes of N=1
DO iElem=1,PP_nElems
  ASSOCIATE( NodeID => GEO%ElemToNodeID(:,iElem) )
    ! Copy values to equidistant distribution
    NodeSourceExtEqui(1,0,0,0) = NodeSourceExt(NodeID(1))/CellLocNodes_Volumes(NodeID(1))
    NodeSourceExtEqui(1,1,0,0) = NodeSourceExt(NodeID(2))/CellLocNodes_Volumes(NodeID(2))
    NodeSourceExtEqui(1,1,1,0) = NodeSourceExt(NodeID(3))/CellLocNodes_Volumes(NodeID(3))
    NodeSourceExtEqui(1,0,1,0) = NodeSourceExt(NodeID(4))/CellLocNodes_Volumes(NodeID(4))
    NodeSourceExtEqui(1,0,0,1) = NodeSourceExt(NodeID(5))/CellLocNodes_Volumes(NodeID(5))
    NodeSourceExtEqui(1,1,0,1) = NodeSourceExt(NodeID(6))/CellLocNodes_Volumes(NodeID(6))
    NodeSourceExtEqui(1,1,1,1) = NodeSourceExt(NodeID(7))/CellLocNodes_Volumes(NodeID(7))
    NodeSourceExtEqui(1,0,1,1) = NodeSourceExt(NodeID(8))/CellLocNodes_Volumes(NodeID(8))
    ! Map equidistant distribution to G/GL (current node type)
    CALL ChangeBasis3D(1, 1, PP_N, Vdm_EQ_N, NodeSourceExtEqui(:,:,:,:),NodeSourceExtGlobal(:,:,:,:,iElem))
  END ASSOCIATE
END DO!iElem

! Write data twice to .h5 file
! 1. to separate file (for visu)
! 2. to _State_.h5 file (or restart)
DO i = 1, 2
  IF(i.EQ.1)THEN
    !FutureTime=0.0
    ! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_NodeSourceExtGlobal',OutputTime))//'.h5'
    IF(MPIRoot) CALL GenerateFileSkeleton('NodeSourceExtGlobal',N_variables,StrVarNames,TRIM(MeshFile),OutputTime)!,FutureTime)
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif
    CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_WORLD)
    CALL WriteAttributeToHDF5(File_ID,'VarNamesNodeSourceExtGlobal',N_variables,StrArray=StrVarNames)
    CALL CloseDataFile()
    DataSetName='DG_Solution'
  ELSE
    ! Write again, but to _State_.h5 file (or restart)
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))//'.h5'
    DataSetName='DG_SourceExt'
  END IF ! i.EQ.2

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        nGlobalElems    => INT(nGlobalElems,IK)    ,&
        PP_nElems       => INT(PP_nElems,IK)       ,&
        N_variables     => INT(N_variables,IK)     ,&
        PP_N            => INT(PP_N,IK)            ,&
        offsetElem      => INT(offsetElem,IK)      )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
        DataSetName=TRIM(DataSetName) , rank=5                                             , &
        nValGlobal =(/N_variables     , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , nGlobalElems/) , &
        nVal       =(/N_variables     , PP_N+1_IK , PP_N+1_IK , PP_N+1_IK , PP_nElems   /) , &
        offset     =(/       0_IK     , 0_IK      , 0_IK      , 0_IK      , offsetElem  /) , &
        collective =.TRUE.            , RealArray=NodeSourceExtGlobal)
  END ASSOCIATE
END DO ! i = 1, 2

SDEALLOCATE(NodeSourceExtGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteNodeSourceExtToHDF5
#endif /*PARTICLES*/



END MODULE MOD_HDF5_output
