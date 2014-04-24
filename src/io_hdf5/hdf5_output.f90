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

INTERFACE WriteArray
  MODULE PROCEDURE WriteArrayToHDF5
END INTERFACE

INTERFACE WriteAttribute
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteStateToHDF5,FlushHDF5,WriteHDF5Header
PUBLIC :: WriteArray,WriteAttribute,WriteArrayToHDF5,WriteAttributeToHDF5
PUBLIC :: WriteToHDF5_multiD
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
USE MOD_Interpolation_Vars,   ONLY:StrNodeType
USE MOD_Mesh_Vars,            ONLY:offsetElem,nGlobalElems
USE MOD_PML_Vars,             ONLY:DoPML,PMLToElem,U2,nPMLElems
#ifdef PP_POIS
USE MOD_Equation_Vars,        ONLY:E,Phi
USE MOD_Mesh_Vars,            ONLY:nElems, sJ, Elem_xGP
#endif
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
INTEGER(HID_T)                 :: Dset_ID
INTEGER                        :: nVal
CHARACTER(LEN=255)             :: FileName,FileString,MeshFile255,StrVarNames(PP_nVar),Statedummy
#ifdef MPI
REAL                           :: StartT,EndT
#endif

#ifdef PP_POIS
REAL                           :: Utemp(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
#endif
REAL,ALLOCATABLE               :: Upml(:,:,:,:,:)
INTEGER                        :: iPML
#ifdef MPI 
INTEGER                        :: sendbuf(2),recvbuf(2)
#endif
!===================================================================================================================================
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' WRITE STATE TO HDF5 FILE...'
#ifdef MPI
  StartT=MPI_WTIME()
#endif
END IF

! Create dataset attribute "VarNames"
#if PP_nVar==8
#ifndef PP_POIS
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ'
StrVarNames(4)='MagneticFieldX'
StrVarNames(5)='MagneticFieldY'
StrVarNames(6)='MagneticFieldZ'
StrVarNames(7)='Phi'       
StrVarNames(8)='Psi'
#else
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ'
StrVarNames(4)='MagneticFieldX'
StrVarNames(5)='MagneticFieldY'
StrVarNames(6)='MagneticFieldZ'
StrVarNames(7)='Phi_LMP'       
StrVarNames(8)='Potential'       
#endif
#endif
#if PP_nVar==4
#ifndef PP_POIS
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ' 
StrVarNames(4)='Psi'        
#else
StrVarNames(1)='Phi'
StrVarNames(2)='Ex'
StrVarNames(3)='Ey' 
StrVarNames(4)='Ez'        
#endif       
#endif

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))//'.h5'
IF(MPIRoot) CALL GenerateFileSkeleton('State',PP_nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
#ifdef MPI
CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
#endif

CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.)

! Write DG solution ----------------------------------------------------------------------------------------------------------------
nVal=nGlobalElems  ! For the MPI case this must be replaced by the global number of elements (sum over all procs)
#ifdef PP_POIS
#if (PP_nVar==8)
Utemp(8,:,:,:,:)=Phi(1,:,:,:,:)
Utemp(1:3,:,:,:,:)=E(1:3,:,:,:,:)
Utemp(4:7,:,:,:,:)=U(4:7,:,:,:,:)
CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.TRUE.,RealArray=Utemp)
CALL WriteArrayToHDF5('DG_SolutionE',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.FALSE.,RealArray=U)
CALL WriteArrayToHDF5('DG_SolutionPhi',nVal,5,(/4,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.FALSE.,RealArray=Phi)
#endif
#if (PP_nVar==4)
Utemp(1,:,:,:,:)=Phi(1,:,:,:,:)
Utemp(2:4,:,:,:,:)=E(1:3,:,:,:,:)
CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.TRUE.,RealArray=Utemp)
CALL WriteArrayToHDF5('DG_SolutionE',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.FALSE.,RealArray=U)
CALL WriteArrayToHDF5('DG_SolutionPhi',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.FALSE.,RealArray=Phi)
#endif
#else
CALL WriteArrayToHDF5('DG_Solution',nVal,5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
,offsetElem,5,existing=.TRUE.,RealArray=U)
#endif


#ifdef PARTICLES
CALL WriteParticleToHDF5()
#endif /*Particles*/

#if (PP_nVar==8)
IF(DoPML)THEN
  ALLOCATE(UPML(6,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  UPML=0.0
  DO iPML=1,nPMLElems
    Upml(:,:,:,:,PMLToElem(iPML)) = U2(:,:,:,:,iPML)
  END DO ! iPML
  CALL WriteArrayToHDF5('PML_Solution',nVal,5,(/6,PP_N+1,PP_N+1,PP_N+1,PP_nElems/) &
    ,offsetElem,5,existing=.FALSE.,RealArray=UPML)
  DEALLOCATE(UPML)
END IF ! DoPML
#endif

! Close the dataset and property list.
CALL H5DCLOSE_F(Dset_id, iError)

! Close the file.
CALL CloseDataFile()

#ifdef MPI
IF(MPIROOT)THEN
  EndT=MPI_WTIME()
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',EndT-StartT,'s]'
END IF
#else
WRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
#endif
END SUBROUTINE WriteStateToHDF5

#ifdef PARTICLES
SUBROUTINE WriteParticleToHDF5()
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,        ONLY:ProjectName
USE MOD_Interpolation_Vars, ONLY:StrNodeType
USE MOD_Mesh_Vars,          ONLY:nGlobalElems, offsetElem
USE MOD_Particle_Vars,      ONLY:PDM, PEM, PartState, PartSpecies, PartMPF, usevMPF,enableParticleMerge,PartPressureCell
USE MOD_part_tools,         ONLY:UpdateNextFreePosition
USE MOD_DSMC_Vars,          ONLY:UseDSMC, CollisMode,PartStateIntEn, DSMC
USE MOD_LD_Vars,            ONLY:UseLD, PartStateBulkValues
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Dset_ID
INTEGER                        :: nVal
#ifdef MPI
REAL                           :: StartT,EndT
#endif
INTEGER(KIND=8)                :: pcount
INTEGER(KIND=4)                :: nParticles(0:nProcessors-1)
LOGICAL                        :: withDSMC=.FALSE.
INTEGER                        :: locnPart,offsetnPart
INTEGER                        :: iPart,nPart_glob, iElem, iElem_glob, iElem_loc
INTEGER,ALLOCATABLE            :: PartInt(:,:)
REAL,ALLOCATABLE               :: PartData(:,:)
INTEGER,PARAMETER              :: PartIntSize=2        !number of entries in each line of PartInt
INTEGER                        :: PartDataSize       !number of entries in each line of PartData
#ifdef MPI 
INTEGER                        :: sendbuf(2),recvbuf(2)
#endif
!=============================================
! Write properties -----------------------------------------------------------------------------------------------------------------
! Open dataset
!CALL H5DOPEN_F(File_ID,'DG_Solution',Dset_id,iError)
 
  !!added for Evib, Erot writeout
  withDSMC=useDSMC
  IF (withDSMC.AND.(.NOT.(useLD))) THEN
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
#else
  offsetnPart=0
  nPart_glob=locnPart
#endif
  ALLOCATE(PartInt(offsetElem+1:offsetElem+PP_nElems,PartIntSize))
  ALLOCATE(PartData(offsetnPart+1:offsetnPart+locnPart,PartDataSize))

!!! Kleiner Hack von JN (Teil 1/2):
  
  IF (.NOT.(useDSMC.OR.enableParticleMerge.OR.PartPressureCell)) THEN
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
    IF (ASSOCIATED(PEM%pNumber)) THEN
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
          IF ((CollisMode.GT.1).AND.(usevMPF) .AND. (DSMC%ElectronicState) ) THEN
            PartData(iPart,8)=PartStateIntEn(pcount,1)
            PartData(iPart,9)=PartStateIntEn(pcount,2)    
            PartData(iPart,10)=PartMPF(pcount)
            PartData(iPart,11)=PartStateIntEn(pcount,3)    
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
      WRITE(*,*) "ERROR: Particle HDF5-Output method not supported! PEM%pNumber not associated"
      STOP
    END IF
    PartInt(iElem_glob,2)=iPart
  END DO 

  CALL WriteArrayToHDF5('PartInt',nGlobalElems,2,(/PP_nElems,PartIntSize/),offsetElem,1,existing=.FALSE.,IntegerArray=PartInt)
  CALL WriteArrayToHDF5('PartData',nPart_glob,2,(/locnPart,PartDataSize/),offsetnPart,1,existing=.FALSE.,RealArray=PartData)!,&
                        !xfer_mode_independent=.TRUE.)  ! könnte bei Procs die keine Teilchen schreiben 
                                                        ! problematisch werden

  DEALLOCATE(PartInt,PartData)

!!! Kleiner Hack von JN (Teil 2/2):
  useDSMC=withDSMC
  IF (.NOT.(useDSMC.OR.enableParticleMerge.OR.PartPressureCell)) THEN
    DEALLOCATE(PEM%pStart , &
               PEM%pNumber, &
               PEM%pNext  , &
               PEM%pEnd   )!, &
               !PDM%nextUsedPosition  )
  END IF
!!! Ende kleiner Hack von JN (Teil 2/2)


END SUBROUTINE WriteParticleToHDF5
#endif /*PARTICLES*/

SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName,OutputTime,FutureTime)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Output_Vars,ONLY:ProjectName
USE MOD_Interpolation_Vars,ONLY:StrNodeType
USE MOD_Mesh_Vars,ONLY:nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
INTEGER                        :: nVar
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
!===================================================================================================================================
! Create file
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.)

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
CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
MeshFile255=TRIM(MeshFileName)
CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=MeshFile255)
IF(PRESENT(FutureTime))THEN
  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=MeshFile255)
END IF
CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=StrNodeType)
CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)

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
CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
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
  CALL GetHDF5NextFileName(Inputfile,NextFile,.TRUE.)
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

! First write program name
tmp255=TRIM(ProgramName)
CALL WriteAttributeToHDF5(File_ID,'Program',1,StrScalar=tmp255)

! Second, write file type identifier
tmp255=TRIM(FileType_in)
CALL WriteAttributeToHDF5(File_ID,'File_Type',1,StrScalar=tmp255)

! Third, write project name
tmp255=TRIM(ProjectName)
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=tmp255)

! Last, file version number
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header




SUBROUTINE WriteArrayToHDF5(DataSetName,nValglobal,Rank,nVal,offset_in,offset_dim,RealArray,IntegerArray,StrArray,existing)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                        :: Rank                  ! number of dimensions of the array
INTEGER,INTENT(IN)             :: offset_in             ! offset =0, start at beginning of the array
INTEGER,INTENT(IN)             :: offset_dim            ! which dimension is the offset (only one dimension possible here)
INTEGER,INTENT(IN)             :: nValglobal            ! max size of array in offset dimension
INTEGER,INTENT(IN)             :: nVal(Rank)            ! size of complete (local) array to write
CHARACTER(LEN=*),INTENT(IN)    :: DataSetName
LOGICAL,INTENT(IN)             :: existing
REAL              ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: RealArray
INTEGER           ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: IntegerArray
CHARACTER(LEN=255),DIMENSION(Rank),OPTIONAL,INTENT(IN) :: StrArray
!TYPE(IO_elem),OPTIONAL,INTENT(IN) :: IO_elemArray(nVal(1))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER(SIZE_T)                :: tmp1,tmp2,typeoffset,typesize
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(Rank),Offset(Rank)
INTEGER(SIZE_T)                :: SizeSet
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! Get global array size, always last dimension!!
Dimsf=nVal
Dimsf(offset_dim)=nValGlobal 

! Create the dataset with default properties.
IF(PRESENT(RealArray))     HDF5DataType=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  HDF5DataType=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, HDF5DataType, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(HDF5DataType, SizeSet, iError)
END IF

IF(existing)THEN
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
ELSE
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), HDF5DataType, FileSpace, DSet_ID, iError)
  CALL H5SCLOSE_F(FileSpace, iError)
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
Offset(:)    = 0
Offset(offset_dim) = Offset_in
! Create the data space in the memory
IF(Dimsf(offset_dim) .NE. 0)THEN
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
ELSE
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(Dimsf(offset_dim) .NE. 0)THEN
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
ELSE
  CALL H5SSELECT_NONE_F(FileSpace,iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#ifdef MPI
!  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError)
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError) ! könnte relevant sein für den
                                                                     ! xfer_mode
#endif
!Write the dataset collectively.
IF(PRESENT(IntegerArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,IntegerArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,RealArray   ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,StrArray    ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF

! Close the property list.
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArrayToHDF5



SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                                                      RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T),INTENT(IN)              :: Loc_ID_in
CHARACTER(LEN=*), INTENT(IN)           :: AttribName
INTEGER,INTENT(IN)                     :: nVal
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: DatasetName
REAL,OPTIONAL,INTENT(IN)               :: RealArray(nVal)
INTEGER,OPTIONAL,INTENT(IN)            :: IntegerArray(nVal)
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: StrArray(nVal)
REAL,OPTIONAL,INTENT(IN)               :: RealScalar
INTEGER,OPTIONAL,INTENT(IN)            :: IntegerScalar
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: StrScalar
LOGICAL,OPTIONAL,INTENT(IN)            :: LogicalScalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,aType_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER                        :: logtoint
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
! Write the attribute data.
IF(PRESENT(RealArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
END IF
IF(PRESENT(RealScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerArray, Dimsf, iError)
END IF
IF(PRESENT(IntegerScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerScalar, Dimsf, iError)
END IF
IF(PRESENT(LogicalScalar))THEN
  IF(logicalScalar)THEN
    logtoint=1
  ELSE
    logtoint=0
  END IF
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, logtoint, Dimsf, iError)
END IF
IF(PRESENT(StrScalar))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, atype_id, iError)
  AttrLen=255
  CALL H5TSET_SIZE_F(aType_ID, AttrLen, iError)
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), aType_ID, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, aType_ID, StrScalar, Dimsf, iError)
END IF
IF(PRESENT(StrArray))THEN
  ! Create character string array datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, atype_id, iError)
  AttrLen=255
  CALL H5TSET_SIZE_F(aType_ID, AttrLen, iError)
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), aType_ID, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, aType_ID, StrArray, Dimsf, iError)
END IF
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

SUBROUTINE WriteToHDF5_multiD(DataSetName,rank,dimsf,counter,offset,RealArray,existing,WriteData_in)
!===================================================================================================================================
! Subroutine to read data from HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)            :: DataSetName
INTEGER,INTENT(IN)                     :: rank
INTEGER(HSIZE_T), INTENT(IN)           :: dimsf(rank)
INTEGER(HSIZE_T), INTENT(IN)           :: counter(rank)
INTEGER(HSIZE_T),INTENT(IN)            :: offset(rank)
REAL,INTENT(IN)                        :: RealArray(1:*)!counter(1),counter(2),counter(3),counter(4),counter(5),counter(6))
LOGICAL,INTENT(IN)                     :: existing
LOGICAL,OPTIONAL,INTENT(IN)            :: WriteData_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                         :: plist_id,dset_id,HDF5DataType
INTEGER(HID_T)                         :: memspace,filespace
LOGICAL                                :: WriteData
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! Create the dataset with default properties.
HDF5DataType=H5T_NATIVE_DOUBLE

IF(existing)THEN
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
ELSE
  ! Create the data space for the dataset. 
  CALL H5SCREATE_SIMPLE_F(rank,dimsf,filespace,iError)
  ! Create the dataset with default properties.
  CALL H5DCREATE_F(file_id,TRIM(DataSetName),HDF5DataType,filespace,dset_id,iError)
  CALL H5SCLOSE_F(filespace,iError)
END IF
! In the MPI case there may be processors that do not need to write to the data set (e.g. BCFace_xGP). Therefor we have the optional
! argument "WriteData_in".
#ifdef MPI
IF(PRESENT(WriteData_in))THEN
  WriteData=WriteData_in
ELSE
  WriteData=.TRUE.
END IF
#else
WriteData=.TRUE.
#endif
! Each process defines dataset in memory and writes it to the hyperslab in the file.
CALL H5SCREATE_SIMPLE_F(rank,counter,memspace,iError)
! Select hyperslab in the file. If no data has to be written, select empty space in file and memory
CALL h5dget_space_f(dset_id,filespace,iError)
IF(WriteData)THEN 
  CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,counter,iError)
ELSE
  CALL h5sselect_none_f(memspace,iError)
  CALL h5sselect_none_f(filespace,iError)
END IF
! Create property list for collective dataset write
CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,iError)
#ifdef MPI
! Set property list to collective dataset read
CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, iError)
#endif
! Write the dataset collectively.
!IF(WriteData)THEN
  CALL h5dwrite_f(dset_id,HDF5DataType,RealArray,dimsf,iError,                   &
                  file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
!END IF
! Close dataspaces.
CALL h5sclose_f(filespace,iError)
CALL h5sclose_f(memspace,iError)
! Close the dataset and property list.
CALL h5dclose_f(dset_id, iError)
CALL h5pclose_f(plist_id, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteToHDF5_multiD 

END MODULE MOD_HDF5_output
