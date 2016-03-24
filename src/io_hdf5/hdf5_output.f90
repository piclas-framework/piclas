#include "boltzplatz.h"

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

!INTERFACE WriteArrayToHDF5
!  MODULE PROCEDURE WriteArrayToHDF5
!END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteStateToHDF5,FlushHDF5,WriteHDF5Header
PUBLIC :: WriteArrayToHDF5,WriteAttributeToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteStateToHDF5(MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
! Subroutine to write the solution U to HDF5 format
! Is used for postprocessing and for restart
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,              ONLY:U
USE MOD_Output_Vars,          ONLY:ProjectName
USE MOD_Mesh_Vars,            ONLY:offsetElem,nGlobalElems
USE MOD_Equation_Vars,        ONLY:StrVarNames
USE MOD_Restart_Vars,         ONLY:RestartFile
#ifdef PARTICLES
USE MOD_PICDepo_Vars,         ONLY:OutputSource,Source
#endif /*PARTICLES*/
#ifdef MPI
USE MOD_LoadBalance_Vars,     ONLY:DoLoadBalance
#endif /*MPI*/
#ifdef PP_POIS
USE MOD_Equation_Vars,        ONLY:E,Phi
#endif /*PP_POIS*/
#ifdef PP_HDG
USE MOD_Particle_Boundary_Vars,ONLY: SurfMesh
USE MOD_Mesh_Vars,            ONLY: offsetSide, nSides,nGlobalUniqueSides,nUniqueSides
USE MOD_HDG_Vars,             ONLY: lambda, nGP_face, nGP_vol, RHS_vol
#if PP_nVar==1
USE MOD_Equation_Vars,        ONLY:E
#elif PP_nVar==3
USE MOD_Equation_Vars,        ONLY:B
#else
USE MOD_Equation_Vars,        ONLY:E,B
#endif /*PP_nVar*/
#endif /*PP_HDG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
REAL,INTENT(IN)                :: OutputTime
REAL,INTENT(IN),OPTIONAL       :: FutureTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: LocalStrVarNames(:)
INTEGER                        :: nVar
#ifdef MPI
REAL                           :: StartT,EndT
#endif /*MPI*/

#ifdef PP_POIS
REAL                           :: Utemp(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#elif defined PP_HDG
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
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE...'
#ifdef MPI
  StartT=MPI_WTIME()
#endif

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))//'.h5'
RestartFile=Filename
#ifdef PP_HDG
#if PP_nVar==1
IF(MPIRoot) CALL GenerateFileSkeleton('State',4,StrVarNames,MeshFileName,OutputTime,FutureTime)
#elif PP_nVar==3
IF(MPIRoot) CALL GenerateFileSkeleton('State',3,StrVarNames,MeshFileName,OutputTime,FutureTime)
#else
IF(MPIRoot) CALL GenerateFileSkeleton('State',7,StrVarNames,MeshFileName,OutputTime,FutureTime)
#endif
#else
IF(MPIRoot) CALL GenerateFileSkeleton('State',PP_nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
#endif /*PP_HDG*/


! Reopen file and write DG solution
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

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
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=Utemp)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionE', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=U)

!CALL WriteArrayToHDF5('DG_SolutionPhi',nVal,5,(/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
!,offsetElem,5,existing=.FALSE.,RealArray=Phi)
! missing addiontal attributes and data preparation
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionPhi', rank=5,&
                        nValGlobal=(/4,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==8)*/
! Store the solution of the electrostatic-poisson system
#if (PP_nVar==4)
Utemp(1,:,:,:,:)=Phi(1,:,:,:,:)
Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=Utemp)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionE', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=U)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionPhi', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=Phi)
#endif /*(PP_nVar==4)*/
DEALLOCATE(Utemp)
#elif defined PP_HDG
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionLambda', rank=3,&
                        nValGlobal=(/PP_nVar,nGP_face,nGlobalUniqueSides/),&
                        nVal=      (/PP_nVar,nGP_face,nUniqueSides/),&
                        offset=    (/0,      0,       offsetSide/),&
                        collective=.TRUE., RealArray=lambda(:,:,1:nUniqueSides))
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_SolutionU', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=U)
#if (PP_nVar==1)
Utemp(1,:,:,:,:)=U(1,:,:,:,:)
Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/4,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=Utemp)

#elif (PP_nVar==3)
Utemp(1:3,:,:,:,:)=B(1:3,:,:,:,:)
!CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
!,offsetElem,5,existing=.TRUE.,RealArray=Utemp)
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/3,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/3,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=Utemp)
#else /*(PP_nVar==4)*/
Utemp(1,:,:,:,:)=U(4,:,:,:,:)
Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
Utemp(5:7,:,:,:,:)=B(1:3,:,:,:,:)

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/7,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/7,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE., RealArray=Utemp)
#endif /*(PP_nVar==1)*/
#else
CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='DG_Solution', rank=5,&
                        nValGlobal=(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                        nVal=      (/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                        offset=    (/0,      0,     0,     0,     offsetElem/),&
                        collective=.TRUE.,RealArray=U)
#endif /*PP_POIS*/
                   

#ifdef PARTICLES
! output of last source term
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif /*MPI*/
IF(OutPutSource) THEN
  ! output of pure current and density
  ! not scaled with epsilon0 and c_corr
  nVar=4
  ALLOCATE(LocalStrVarNames(1:nVar))
  LocalStrVarNames(1)='CurrentDensityX'
  LocalStrVarNames(2)='CurrentDensityY'
  LocalStrVarNames(3)='CurrentDensityZ'
  LocalStrVarNames(4)='ChargeDensity'
  IF(MPIRoot)THEN
#ifdef MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
    CALL OpenDataFile(FileName,create=.FALSE.)
#endif /*MPI*/
    CALL WriteAttributeToHDF5(File_ID,'VarNamesSource',nVar,StrArray=LocalStrVarnames)
    CALL CloseDataFile()
  END IF
  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='DG_Source', rank=5,  &
                          nValGlobal=(/nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                          nVal=      (/nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                          offset=    (/0,      0,     0,     0,     offsetElem/),&
                          collective=.TRUE.,RealArray=Source)

  DEALLOCATE(LocalStrVarNames)
END IF
#endif /*PARTICLES*/


#ifdef PARTICLES
CALL WriteParticleToHDF5(FileName)
#endif /*Particles*/

CALL WriteAdditionalDataToHDF5(FileName)

CALL WritePMLDataToHDF5(FileName)

#ifdef MPI
IF(DoLoadBalance) CALL WriteElemWeightToHDF5(FileName)
#endif /*MPI*/

#ifdef MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
END SUBROUTINE WriteStateToHDF5


SUBROUTINE WriteElemWeightToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars        ,ONLY:offsetElem,nGlobalElems
USE MOD_LoadBalance_Vars ,ONLY:ElemTime,nLoadBalance,ElemWeight
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
!===================================================================================================================================


IF(nLoadBalance.EQ.0) THEN
  IF(.NOT.ALLOCATED(ElemWeight)) RETURN
  ElemTime=ElemWeight
END IF

nVar=1
ALLOCATE(StrVarNames(nVar))
StrVarNames(1)='ElemWeight'

IF(MPIRoot)THEN
#ifdef MPI
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
  CALL OpenDataFile(FileName,create=.FALSE.)
#endif
  CALL WriteAttributeToHDF5(File_ID,'VarNamesLB',nVar,StrArray=StrVarNames)
  CALL CloseDataFile()
END IF

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='ElemWeight', rank=1,  &
                        nValGlobal=(/nGlobalElems/),&
                        nVal=      (/PP_nElems   /),&
                        offset=    (/offsetElem  /),&
                        collective=.TRUE.,RealArray=ElemTime)

!CALL WriteArrayToHDF5(DataSetName='ElemWeight', rank=1,&
!                    nValGlobal=(/nGlobalElems/),       &
!                    nVal=      (/PP_nElems/),          &
!                    offset=    (/offsetElem/),         &
!                    collective=.TRUE., existing=.FALSE., RealArray=ElemTime)
!
DEALLOCATE(StrVarNames)
ElemTime=0.

END SUBROUTINE WriteElemWeightToHDF5


SUBROUTINE WriteAdditionalDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
! 1-Rank
! 2-isPMLElem
! 3-hasParticle
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:offsetElem,nGlobalElems
USE MOD_PML_Vars      ,ONLY:DoPML,PMLToElem,nPMLElems
#ifdef PARTICLES
USE MOD_Particle_Vars ,ONLY:PDM,PEM
#endif /*PARTICLES*/
#ifdef MPI
USE MOD_Loadbalance_Vars,  ONLY:DoLoadBalance,ElemWeight
#endif /*MPI*/
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
REAL,ALLOCATABLE               :: ElemData(:,:)
INTEGER                        :: nVar,iElem
#ifdef PARTICLES
INTEGER                        :: iPart
#endif /*PARTICLES*/
!===================================================================================================================================


#ifdef MPI
IF(DoLoadBalance)THEN
  nVar=4
ELSE
  nVar=3
END IF
#else
nVar=3
#endif /*MPI*/
ALLOCATE(StrVarNames(nVar))
ALLOCATE(ElemData(nVar,PP_nElems))

StrVarNames(1)='MyRank'
StrVarNames(2)='PMLElement'
StrVarNames(3)='PartElem'

! fill ElemData
ElemData=0.
! MPI-Rank
DO iElem=1,PP_nElems
  ElemData(1,iElem)=REAL(MyRank)
END DO
! PML Info
IF(DoPML)THEN
  DO iElem=1,nPMLElems
    ElemData(2,PMLToElem(iElem))=1.0
  END DO ! iElem=1,nPMLElems
END IF
! hasParticles
#ifdef PARTICLES
DO iPart=1,PDM%ParticleVecLength
  IF(PDM%ParticleInside(iPart))THEN
    iElem=PEM%Element(iPart)
    IF(ElemData(3,iElem).EQ.0.) ElemData(3,iElem)=1.0
  END IF ! ParticleInside
END DO ! iPart
#endif /*PARTICLES*/
  
#ifdef MPI
IF(DoLoadBalance)THEN
  StrVarNames(4)='ElemTime'
  IF(ALLOCATED(ElemWeight))THEN
    DO iElem=1,PP_nElems
      ElemData(4,iElem)=ElemWeight(iElem)
    END DO ! iElem =1,PP_nElems
  ELSE
    ElemData(4,:)=0.
  END IF
END IF
#endif /*MPI*/

! output of rank 
IF(MPIRoot)THEN
#ifdef MPI
  CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
  CALL OpenDataFile(FileName,create=.FALSE.)
#endif
  CALL WriteAttributeToHDF5(File_ID,'VarNamesAdd',nVar,StrArray=StrVarnames)
  CALL CloseDataFile()
END IF

CALL GatheredWriteArray(FileName,create=.FALSE.,&
                        DataSetName='ElemData', rank=2,  &
                        nValGlobal=(/nVar,nGlobalElems/),&
                        nVal=      (/nVar,PP_nElems   /),&
                        offset=    (/0,   offsetElem  /),&
                        collective=.TRUE.,RealArray=ElemData)
DEALLOCATE(ElemData,StrVarNames)


END SUBROUTINE WriteAdditionalDataToHDF5


SUBROUTINE WritePMLDataToHDF5(FileName)
!===================================================================================================================================
! Write additional (elementwise scalar) data to HDF5
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:offsetElem,nGlobalElems
USE MOD_PML_Vars      ,ONLY:DoPML,PMLToElem,U2,nPMLElems
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
REAL,ALLOCATABLE               :: Upml(:,:,:,:,:)
INTEGER                        :: iPML
!===================================================================================================================================

nVar=0
#if (PP_nVar==8)
IF(DoPML)THEN
  nVar=6
  ALLOCATE(StrVarNames(nVar))
  StrVarNames(1)='PMLElectricFieldX'
  StrVarNames(2)='PMLElectricFieldY'
  StrVarNames(3)='PMLElectricFieldZ'
  StrVarNames(4)='PMLMagneticFieldX'
  StrVarNames(5)='PMLMagneticFieldY'
  StrVarNames(6)='PMLMagneticFieldZ'

  ALLOCATE(UPML(6,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  UPML=0.0
  DO iPML=1,nPMLElems
    Upml(:,:,:,:,PMLToElem(iPML)) = U2(:,:,:,:,iPML)
  END DO ! iPML

  IF(MPIRoot)THEN
#ifdef MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
    CALL OpenDataFile(FileName,create=.FALSE.)
#endif
    CALL WriteAttributeToHDF5(File_ID,'VarNamesPML',nVar,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='PML_Solution', rank=5,&
                          nValGlobal=(/nVar,PP_N+1,PP_N+1,PP_N+1,nGlobalElems/),&
                          nVal=      (/nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),&
                          offset=    (/0,      0,     0,     0,     offsetElem/),&
                          collective=.TRUE.,RealArray=Upml)

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
#endif

END SUBROUTINE WritePMLDataToHDF5

#ifdef PARTICLES
SUBROUTINE WriteParticleToHDF5(FileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY:nGlobalElems, offsetElem
USE MOD_Particle_Vars,      ONLY:PDM, PEM, PartState, PartSpecies, PartMPF, usevMPF,PartPressureCell
USE MOD_part_tools,         ONLY:UpdateNextFreePosition
USE MOD_DSMC_Vars,          ONLY:UseDSMC, CollisMode,PartStateIntEn, DSMC
USE MOD_LD_Vars,            ONLY:UseLD, PartStateBulkValues
#ifdef MPI
USE MOD_Particle_MPI_Vars,  ONLY:PartMPI
#endif /*MPI*/
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
#ifdef MPI
INTEGER                        :: sendbuf(2),recvbuf(2)
INTEGER                        :: nParticles(0:nProcessors-1)
#endif
LOGICAL                        :: reSwitch
INTEGER                        :: pcount
LOGICAL                        :: withDSMC=.FALSE.
INTEGER                        :: locnPart,offsetnPart
INTEGER                        :: iPart,nPart_glob, iElem_glob, iElem_loc
INTEGER,ALLOCATABLE            :: PartInt(:,:)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER,PARAMETER              :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
INTEGER                        :: minnParts
#endif /* HDF5_F90 */
!=============================================
! Write properties -----------------------------------------------------------------------------------------------------------------
! Open dataset
!CALL H5DOPEN_F(File_ID,'DG_Solution',Dset_id,iError)
 
  !!added for Evib, Erot writeout
  withDSMC=useDSMC
  IF (withDSMC.AND.(.NOT.(useLD))) THEN
  !IF (withDSMC) THEN
    IF ((CollisMode.GT.1).AND.(usevMPF) .AND. DSMC%ElectronicState ) THEN !int ener + 3, vmpf +1
      PartDataSize=11
    ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. DSMC%ElectronicState ) ) THEN !int ener + 2 and vmpf + 1
                                                                              ! or int energ +3 but no vmpf +1
      PartDataSize=10
    ELSE IF (CollisMode.GT.1) THEN
      PartDataSize=9 !int ener + 2
    ELSE IF (usevMPF) THEN
      PartDataSize=8 !+ 1 vmpf
    ELSE
      PartDataSize=7 !+ 0
    END IF
  ELSE IF (useLD) THEN
    IF ((CollisMode.GT.1).AND.(usevMPF) .AND. DSMC%ElectronicState ) THEN !int ener + 3, vmpf +1
      PartDataSize=16
    ELSE IF ((CollisMode.GT.1).AND.( (usevMPF) .OR. DSMC%ElectronicState ) ) THEN !int ener + 2 and vmpf + 1
                                                                             ! or int energ +3 but no vmpf +1
      PartDataSize=15
    ELSE IF (CollisMode.GT.1) THEN
      PartDataSize=14!int ener + 2
    ELSE IF (usevMPF) THEN
      PartDataSize=13!+ 1 vmpf
    ELSE
      PartDataSize=12 !+ 0
    END IF
  ELSE IF (usevMPF) THEN
    PartDataSize=8 !vmpf +1
  ELSE
    PartDataSize=7
  END IF  

  locnPart =   0
  DO pcount = 1,PDM%ParticleVecLength
    IF(PDM%ParticleInside(pcount)) THEN
      locnPart = locnPart + 1
    END IF
  END DO         

#ifdef MPI
  sendbuf(1)=locnPart
  recvbuf=0
  CALL MPI_EXSCAN(sendbuf(1),recvbuf(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iError)
  offsetnPart=recvbuf(1)
  sendbuf(1)=recvbuf(1)+locnPart
  CALL MPI_BCAST(sendbuf(1),1,MPI_INTEGER,nProcessors-1,MPI_COMM_WORLD,iError) !last proc knows global number
  !global numbers
  nPart_glob=sendbuf(1)
  CALL MPI_GATHER(locnPart,1,MPI_INTEGER,nParticles,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
  !IF (myRank.EQ.0) THEN
  !  WRITE(*,*) 'PARTICLE-ELEMENT DISTRIBUTION' 
  !  WRITE(*,*) 'iProc, firstelemInd,   nElems,  locnPart,  totalnPart'
  !  DO pcount=0,nProcessors-1
  !    WRITE(*,'(I5,4I12)')pcount,offsetElemMPI(pcount),offsetElemMPI(pcount+1)-offsetElemMPI(pcount),&
  !                       nParticles(pcount),SUM(nParticles(0:pcount))
  !  END DO
  !END IF
  LOGWRITE(*,*)'offsetnPart,locnPart,nPart_glob',offsetnPart,locnPart,nPart_glob
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
  CALL MPI_ALLREDUCE(locnPart, minnParts, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, IERROR)
#endif /* HDF5_F90 */
#else
  offsetnPart=0
  nPart_glob=locnPart
#ifdef HDF5_F90 /* HDF5 compiled without fortran2003 flag */
  minnParts=locnPart
#endif /* HDF5_F90 */
#endif
  ALLOCATE(PartInt(offsetElem+1:offsetElem+PP_nElems,PartIntSize))
  ALLOCATE(PartData(offsetnPart+1:offsetnPart+locnPart,PartDataSize))

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
      PartInt(iElem_glob,2) = PartInt(iElem_glob,1) + PEM%pNumber(iElem_loc)
      pcount = PEM%pStart(iElem_loc)
      DO iPart=PartInt(iElem_glob,1)+1,PartInt(iElem_glob,2)
        PartData(iPart,1)=PartState(pcount,1)
        PartData(iPart,2)=PartState(pcount,2)
        PartData(iPart,3)=PartState(pcount,3)
        PartData(iPart,4)=PartState(pcount,4)
        PartData(iPart,5)=PartState(pcount,5)
        PartData(iPart,6)=PartState(pcount,6)
        PartData(iPart,7)=REAL(PartSpecies(pcount))
        IF (withDSMC.AND.(.NOT.(useLD))) THEN
        !IF (withDSMC) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)    
            PartData(iPart,10)=PartStateIntEn(pcount,3)    
            PartData(iPart,11)=PartMPF(pcount)
          ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)    
            PartData(iPart,10)=PartMPF(pcount)
          ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicState) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)
            PartData(iPart,10)=PartStateIntEn(pcount,3)
          ELSE IF (CollisMode.GT.1) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2) 
          ELSE IF (usevMPF) THEN
            PartData(iPart,8)=PartMPF(pcount)    
          END IF
        ELSE IF (useLD) THEN
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)    
            PartData(iPart,10)=PartMPF(pcount)
            PartData(iPart,11)=PartStateIntEn(pcount,3)  
            PartData(iPart,12)=PartStateBulkValues(pcount,1)    
            PartData(iPart,13)=PartStateBulkValues(pcount,2)
            PartData(iPart,14)=PartStateBulkValues(pcount,3)
            PartData(iPart,15)=PartStateBulkValues(pcount,4)
            PartData(iPart,16)=PartStateBulkValues(pcount,5)
          ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)    
            PartData(iPart,10)=PartMPF(pcount)
            PartData(iPart,11)=PartStateBulkValues(pcount,1)    
            PartData(iPart,12)=PartStateBulkValues(pcount,2)
            PartData(iPart,13)=PartStateBulkValues(pcount,3)
            PartData(iPart,14)=PartStateBulkValues(pcount,4)
            PartData(iPart,15)=PartStateBulkValues(pcount,5)
          ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicState) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)
            PartData(iPart,10)=PartStateIntEn(pcount,3)
            PartData(iPart,11)=PartStateBulkValues(pcount,1)    
            PartData(iPart,12)=PartStateBulkValues(pcount,2)
            PartData(iPart,13)=PartStateBulkValues(pcount,3)
            PartData(iPart,14)=PartStateBulkValues(pcount,4)
            PartData(iPart,15)=PartStateBulkValues(pcount,5)
          ELSE IF (CollisMode.GT.1) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2) 
            PartData(iPart,10)=PartStateBulkValues(pcount,1)    
            PartData(iPart,11)=PartStateBulkValues(pcount,2)
            PartData(iPart,12)=PartStateBulkValues(pcount,3)
            PartData(iPart,13)=PartStateBulkValues(pcount,4)
            PartData(iPart,14)=PartStateBulkValues(pcount,5)
          ELSE IF (usevMPF) THEN
            PartData(iPart,8)=PartMPF(pcount)    
            PartData(iPart,9)=PartStateBulkValues(pcount,1)    
            PartData(iPart,10)=PartStateBulkValues(pcount,2)
            PartData(iPart,11)=PartStateBulkValues(pcount,3)
            PartData(iPart,12)=PartStateBulkValues(pcount,4)
            PartData(iPart,13)=PartStateBulkValues(pcount,5)
          ELSE
            PartData(iPart,8)=PartStateBulkValues(pcount,1)    
            PartData(iPart,9)=PartStateBulkValues(pcount,2)
            PartData(iPart,10)=PartStateBulkValues(pcount,3)
            PartData(iPart,11)=PartStateBulkValues(pcount,4)
            PartData(iPart,12)=PartStateBulkValues(pcount,5)
          END IF
        ELSE IF (usevMPF) THEN
            PartData(iPart,8)=PartMPF(pcount)
        END IF
        !PartData(iPart,8)=Species(PartSpecies(pcount))%ChargeIC*Species(PartSpecies(pcount))%MacroParticleFactor
        !PartData(iPart,9)=Species(PartSpecies(pcount))%MassIC*Species(PartSpecies(pcount))%MacroParticleFactor
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
#ifdef MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
    CALL OpenDataFile(FileName,create=.FALSE.)
#endif
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

  CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName='PartInt', rank=2,&
                          nValGlobal=(/nGlobalElems,nVar/),&
                          nVal=      (/PP_nElems,nVar/),&
                          offset=    (/offsetElem,0/),&
                          collective=.TRUE.,IntegerArray=PartInt)

!  CALL WriteArrayToHDF5(DataSetName='PartInt', rank=2,&
!                        nValGlobal=(/nGlobalElems, nVar/),&
!                        nVal=      (/PP_nElems, nVar   /),&
!                        offset=    (/offsetElem, 0  /),&
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
  IF(withDSMC.AND.(.NOT.(useLD)))THEN
 ! IF(withDSMC)THEN
    IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicState))THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
      StrVarNames(11)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='MPF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicState) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
    ELSE IF (CollisMode.GT.1) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
    ELSE IF (usevMPF) THEN
      StrVarNames( 8)='MPF'
    END IF
  ELSE IF (useLD) THEN
    IF((CollisMode.GT.1).AND.(usevMPF).AND.(DSMC%ElectronicState))THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
      StrVarNames(11)='MPF'
      StrVarNames(12)='BulkVelocityX'
      StrVarNames(13)='BulkVelocityY'
      StrVarNames(14)='BulkVelocityZ'
      StrVarNames(15)='BulkTemperature'
      StrVarNames(16)='BulkDOF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (usevMPF) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='MPF'
      StrVarNames(11)='BulkVelocityX'
      StrVarNames(12)='BulkVelocityY'
      StrVarNames(13)='BulkVelocityZ'
      StrVarNames(14)='BulkTemperature'
      StrVarNames(15)='BulkDOF'
    ELSE IF ( (CollisMode .GT. 1) .AND. (DSMC%ElectronicState) ) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='Electronic'
      StrVarNames(11)='BulkVelocityX'
      StrVarNames(12)='BulkVelocityY'
      StrVarNames(13)='BulkVelocityZ'
      StrVarNames(14)='BulkTemperature'
      StrVarNames(15)='BulkDOF'
    ELSE IF (CollisMode.GT.1) THEN
      StrVarNames( 8)='Vibrational'
      StrVarNames( 9)='Rotational'
      StrVarNames(10)='BulkVelocityX'
      StrVarNames(11)='BulkVelocityY'
      StrVarNames(12)='BulkVelocityZ'
      StrVarNames(13)='BulkTemperature'
      StrVarNames(14)='BulkDOF'
    ELSE IF (usevMPF) THEN
      StrVarNames( 8)='MPF'
      StrVarNames( 9)='BulkVelocityX'
      StrVarNames(10)='BulkVelocityY'
      StrVarNames(11)='BulkVelocityZ'
      StrVarNames(12)='BulkTemperature'
      StrVarNames(13)='BulkDOF'
    ELSE
      StrVarNames( 8)='BulkVelocityX'
      StrVarNames( 9)='BulkVelocityY'
      StrVarNames(10)='BulkVelocityZ'
      StrVarNames(11)='BulkTemperature'
      StrVarNames(12)='BulkDOF'
    END IF
!   CALL abort(__STAMP__,&
!       'Attributes for LD are not implemented! Add Attributes!',999,999.)
  END IF

  IF(MPIRoot)THEN
#ifdef MPI
    CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.)
#else
    CALL OpenDataFile(FileName,create=.FALSE.)
#endif
    CALL WriteAttributeToHDF5(File_ID,'VarNamesParticles',PartDataSize,StrArray=StrVarNames)
    CALL CloseDataFile()
  END IF

#ifdef MPI
 CALL DistributedWriteArray(FileName,create=.FALSE.,&
                            DataSetName='PartData', rank=2         ,&
                            nValGlobal=(/nPart_glob,PartDataSize/) ,&
                            nVal=      (/locnPart,PartDataSize/)   ,&
                            offset=    (/offsetnPart,0/)           ,&
                            collective=.FALSE.,offSetDim=1         ,&
                            communicator=PartMPI%COMM,RealArray=PartData)
#else
  CALL OpenDataFile(FileName,create=.FALSE.)
  CALL WriteArrayToHDF5(DataSetName='PartData', rank=2,&
                        nValGlobal=(/nPart_glob,PartDataSize/),&
                        nVal=      (/locnPart,PartDataSize  /),&
                        offset=    (/offsetnPart , 0  /),&
                        collective=.TRUE., RealArray=PartData)
  CALL CloseDataFile()
#endif /*MPI*/                          

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
#endif /*PARTICLES*/


! PO: old
!SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY: ProjectName
USE MOD_Mesh_Vars  ,ONLY: nGlobalElems
USE MOD_ReadInTools,ONLY: GetParameters
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
CHARACTER(LEN=255),ALLOCATABLE :: params(:)
!===================================================================================================================================
! Create file
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
#ifdef MPI
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.)
#else
CALL OpenDataFile(TRIM(FileName),create=.TRUE.)
#endif

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

! Write ini file parameters and compile flags
CALL GetParameters(params)
CALL WriteAttributeToHDF5(File_ID,'Parameters',SIZE(params),StrArray=params)
!CALL WriteAttributeToHDF5(File_ID,'Compile',1,StrScalar=(/PREPROC_FLAGS/))
DEALLOCATE(params)

CALL CloseDataFile()
END SUBROUTINE GenerateFileSkeleton


SUBROUTINE FlushHDF5(FlushTime_In)
!===================================================================================================================================
! Deletes all HDF5 output files, beginning from time Flushtime
!===================================================================================================================================
! MODULES
!USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY:ProjectName
USE MOD_HDF5_Input,ONLY:GetHDF5NextFileName
#ifdef MPI
USE MOD_Loadbalance_Vars,  ONLY:DoLoadBalance,nLoadBalance
#endif /*MPI*/
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
CHARACTER(LEN=255)       :: FileName,InputFile,NextFile
!===================================================================================================================================
IF(.NOT.MPIRoot) RETURN

#ifdef MPI
IF(DoLoadBalance .AND. nLoadBalance.GT.0) RETURN
#endif /*MPI*/

WRITE(UNIT_stdOut,'(a)')' DELETING OLD HDF5 FILES...'
IF (.NOT.PRESENT(FlushTime_In)) THEN
  FlushTime=0.0
ELSE
  FlushTime=FlushTime_In
END IF
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',FlushTime))//'.h5'

! Delete state files
InputFile=TRIM(FileName)
! Read calculation time from file
#ifdef MPI
CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
! Delete File - only root
stat=0
ioUnit=GETFREEUNIT()
OPEN ( UNIT   = ioUnit,            &
       FILE   = InputFile,      &
       STATUS = 'OLD',          &
       ACTION = 'WRITE',        &
       ACCESS = 'SEQUENTIAL',   &
       IOSTAT = stat          )
IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
DO
  InputFile=TRIM(NextFile)
  ! Read calculation time from file
#ifdef MPI
  CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
#else
  CALL GetHDF5NextFileName(Inputfile,NextFile)
#endif
  ! Delete File - only root
  stat=0
  ioUnit=GETFREEUNIT()
  OPEN ( UNIT   = ioUnit,            &
         FILE   = InputFile,      &
         STATUS = 'OLD',          &
         ACTION = 'WRITE',        &
         ACCESS = 'SEQUENTIAL',   &
         IOSTAT = stat          )
  IF(stat .EQ. 0) CLOSE ( ioUnit,STATUS = 'DELETE' )
  IF(iError.NE.0) EXIT  ! iError is set in GetHDF5NextFileName !
END DO

WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'

END SUBROUTINE FlushHDF5



SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
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
                            RealArray,IntegerArray,StrArray)
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
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName
INTEGER,INTENT(IN)            :: rank             ! number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) ! max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       ! size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     ! offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: collective       ! use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  ! specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  ! specify chunksize
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(rank)
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray(rank)
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
IF(PRESENT(RealArray))     Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  Type_ID=H5T_NATIVE_INTEGER
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
#ifdef MPI
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
IF(PRESENT(IntegerArray)) buf=C_LOC(IntegerArray)
IF(PRESENT(RealArray))    buf=C_LOC(RealArray)
IF(PRESENT(StrArray))     buf=C_LOC(StrArray(1))
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
CHARACTER(LEN=*),INTENT(IN)    :: FileName,DataSetName
LOGICAL,INTENT(IN)             :: create,collective
INTEGER,INTENT(IN)             :: rank,nVal(rank),nValGlobal(rank),offset(rank)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
REAL,              ALLOCATABLE :: UReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: UStr(:)
INTEGER,           ALLOCATABLE :: UInt(:)
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nLocalProcs) :: nDOFPerNode,offsetNode
!===================================================================================================================================
IF(gatheredWrite)THEN
  IF(ANY(offset(1:rank-1).NE.0)) &
    CALL abort(&
    __STAMP__&
    ,'Offset only allowed in last dimension for gathered IO.')
  
  ! Get last dim of each array on IO nodes
  nDOFLocal=PRODUCT(nVal)
  CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_NODE,iError)
  
  ! Allocate big array and compute offsets of small arrs inside big
  offsetNode=0
  IF(MPILocalRoot)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nLocalProcs
      offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    IF(PRESENT(RealArray)) ALLOCATE(UReal(PRODUCT(nValGather)))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( PRODUCT(nValGather)))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( PRODUCT(nValGather)))
  ELSE
    IF(PRESENT(RealArray)) ALLOCATE(UReal(1))
    IF(PRESENT(IntegerArray))  ALLOCATE(UInt( 1))
    IF(PRESENT(StrArray))  ALLOCATE(UStr( 1))
  ENDIF
  
  ! Gather small arrays on IO nodes
  IF(PRESENT(RealArray)) CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                                          UReal,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_NODE,iError)
  IF(PRESENT(IntegerArray))  CALL MPI_GATHERV(IntegerArray, nDOFLocal,MPI_INTEGER,&
                                          UInt, nDOFPerNode,offsetNode,MPI_INTEGER,0,MPI_COMM_NODE,iError)
  !IF(PRESENT(StrArray))  CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE,&
  !                                        UReal,nDOFPerNode, offsetNode,MPI_DOUBLE,0,MPI_COMM_NODE,iError)
  
  IF(MPILocalRoot)THEN
    ! Reopen file and write DG solution (only IO nodes)
    CALL OpenDataFile(FileName,create=create,single=.FALSE.,communicatorOpt=MPI_COMM_LEADERS)
    IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,RealArray=UReal)
    IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nValGather,&
                                                 offset,collective=collective,IntegerArray =UInt)
    !IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nValGather,&
    !                                             offset,collective=collective,StrArr =UStr)
    CALL CloseDataFile()
  END IF
  
  SDEALLOCATE(UReal)
  SDEALLOCATE(UInt)
  SDEALLOCATE(UStr)
ELSE
#endif
#ifdef MPI
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.)
#else
  CALL OpenDataFile(FileName,create=.FALSE.)
#endif
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,RealArray=RealArray)
  IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,IntegerArray =IntegerArray)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
#ifdef MPI
END IF
#endif

END SUBROUTINE GatheredWriteArray

#ifdef MPI
SUBROUTINE DistributedWriteArray(FileName,create,DataSetName,rank,nValGlobal,nVal,offset,collective,&
                                 offSetDim,communicator,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Write distributed data to proc, e.g. particles which are not hosted by each proc
! a new output-communicator is build and afterwards killed
! offset is in the last dimension
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName,DataSetName
LOGICAL,INTENT(IN)             :: create,collective
INTEGER,INTENT(IN)             :: offSetDim,communicator
INTEGER,INTENT(IN)             :: rank,nVal(rank),nValGlobal(rank),offset(rank)
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(PRODUCT(nVal))
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntegerArray( PRODUCT(nVal))
CHARACTER(LEN=255),INTENT(IN),OPTIONAL,TARGET :: StrArray( PRODUCT(nVal))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifdef MPI
INTEGER                        :: Color, OutPutCOMM,nOutPutProcs,MyOutputRank
LOGICAL                        :: DataOnProc, DoNotSplit
!===================================================================================================================================

DataOnProc=.FALSE.
IF(nVal(offSetDim).GT.0) DataOnProc=.TRUE.
CALL MPI_ALLREDUCE(DataOnProc,DoNotSplit, 1, MPI_LOGICAL, MPI_LAND, COMMUNICATOR, IERROR)


IF(.NOT.DoNotSplit)THEN
  color=MPI_UNDEFINED
  IF(DataOnProc) color=87
  MyOutputRank=0

  CALL MPI_COMM_SPLIT(COMMUNICATOR, color, MyOutputRank, OutputCOMM,iError)
  IF(DataOnProc) THEN
    CALL MPI_COMM_SIZE(OutputCOMM, nOutPutProcs,iError)
    IF(nOutPutProcs.EQ.1)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,IntegerArray =IntegerArray)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective=.FALSE.,StrArray =StrArray)
    ELSE 
      CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,communicatorOpt=OutputCOMM)
      IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,RealArray=RealArray)
      IF(PRESENT(IntegerArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,IntegerArray =IntegerArray)
      IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                                   offset,collective,StrArray =StrArray)
    END IF
    CALL CloseDataFile()
    CALL MPI_COMM_FREE(OutputCOMM,iERROR)
  END IF
  ! MPI Barrier is requried, that the other procs don't open the datafile while this procs are still writring
  CALL MPI_BARRIER(COMMUNICATOR,IERROR)
  OutputCOMM=MPI_UNDEFINED
ELSE
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.)
#else
  CALL OpenDataFile(FileName,create=.FALSE.)
#endif
  IF(PRESENT(RealArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,RealArray=RealArray)
  IF(PRESENT(IntegerArray)) CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,IntegerArray =IntegerArray)
  IF(PRESENT(StrArray))  CALL WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,&
                                               offset,collective,StrArray =StrArray)
  CALL CloseDataFile()
#ifdef MPI
END IF
#endif

END SUBROUTINE DistributedWriteArray
#endif /*MPI*/

END MODULE MOD_HDF5_output
