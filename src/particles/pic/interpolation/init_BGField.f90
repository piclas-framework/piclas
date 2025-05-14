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

MODULE  MOD_InitializeBackgroundField
!===================================================================================================================================
! Module performing the init. of a BackGround-Field
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC ::  InitializeBackgroundField,FinalizeBackGroundField,DefineParametersBGField
!===================================================================================================================================

INTERFACE InitializeBackgroundField
  MODULE PROCEDURE InitializeBackgroundField
END INTERFACE

INTERFACE FinalizeBackGroundField
  MODULE PROCEDURE FinalizeBackGroundField
END INTERFACE
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for the background field
!==================================================================================================================================
SUBROUTINE DefineParametersBGField()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("PIC Background Field")

! -- external field 5
CALL prms%CreateLogicalOption('PIC-BG-Field'      , 'Method 5 of 5: Activates the usage of a background field, read-in from file '//&
                                                    '(PIC-BGFileName=BGField.h5) or calculated from parameters', &
                                                    '.FALSE.')
CALL prms%CreateStringOption( 'PIC-BGFileName'    , 'File name for the background field ([character].h5)','none')
CALL prms%CreateRealOption(   'PIC-BGFieldScaling', 'Scaling of the read-in background field','1.')
END SUBROUTINE DefineParametersBGField


SUBROUTINE InitializeBackgroundField
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_Globals_Vars          ,ONLY: ProjectName
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_Basis                 ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_Basis                 ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Mesh_Vars             ,ONLY: OffsetElem,nGlobalElems, nElems, offsetElem
USE MOD_Preproc
USE MOD_ReadInTools           ,ONLY: GETSTR,GETINT,GETREAL
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataProps,DatasetExists
USE MOD_HDF5_Input            ,ONLY: GetDataSize
USE MOD_PICInterpolation_Vars ,ONLY: CalcBField
USE MOD_Interpolation_Vars    ,ONLY: N_BG,BGType,BGDataSize
USE MOD_Interpolation_Vars    ,ONLY: NodeType,N_Inter,Nmin,Nmax
USE MOD_ReadInTools           ,ONLY: PrintOption
USE MOD_DG_Vars               ,ONLY: N_DG_Mapping
#if USE_SUPER_B
USE MOD_SuperB                ,ONLY: SuperB
USE MOD_SuperB_Vars           ,ONLY: BGFieldFrequency,UseTimeDepCoil,nTimePoints
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars      ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: nVar_BField,N_In,nElems_BGField, nDOF_BGField, Nloc, Nin
CHARACTER(LEN=255)                      :: BGFileName,NodeType_BGField,BGFieldName
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
REAL,ALLOCATABLE                        :: BGField_tmp(:,:,:,:,:), BGFieldTDep_tmp(:,:,:,:,:,:)
INTEGER                                 :: iElem,i,j,k,iField,nFields,nDimsOffset
REAL                                    :: BGFieldScaling
LOGICAL                                 :: BGFieldExists,DG_SolutionExists
TYPE tVdm_BGFieldIn_BGField
  REAL, ALLOCATABLE                         :: Vdm(:,:)
END TYPE tVdm_BGFieldIn_BGField
TYPE(tVdm_BGFieldIn_BGField),ALLOCATABLE    :: Vdm_BGFieldIn_BGField(:,:)

TYPE tweights
  REAL,ALLOCATABLE            :: xGP_tmp(:),wBary_tmp(:),wGP_tmp(:)
END TYPE tweights
TYPE(tweights),ALLOCATABLE    :: weights(:)
! p-adaption
INTEGER                            :: iDOF, nDOF, offsetDOF
REAL,ALLOCATABLE                   :: U_N_2D_local(:,:)
LOGICAL                         :: ElemDataExists, NlocFound
INTEGER                         :: iVar, nVarAdd
CHARACTER(LEN=255),ALLOCATABLE  :: VarNamesAdd(:)
INTEGER,ALLOCATABLE             :: Nloc_HDF5(:)                      !< Array for temporary read-in of Nloc container
REAL,ALLOCATABLE                :: ElemData(:,:)                     !< Array for temporary read-in of ElemData container
!===================================================================================================================================

#if USE_SUPER_B
LBWRITE(UNIT_stdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)')' INIT BACKGROUND FIELD...'
NlocFound = .FALSE.

! 1) Determine whether a field will be read-in or calculated
UseTimeDepCoil=.FALSE.! Initialize
BGFileName = GETSTR('PIC-BGFileName')
IF(TRIM(BGFileName).NE.'none') THEN
  ! 1a) Read-in of background field name, if supplied -> read-in of the field
  CalcBField = .FALSE.
ELSE
  BGFileName = TRIM(ProjectName)//'_BGField.h5'
  IF(FILEEXISTS(TRIM(BGFileName))) THEN
  ! 1b) If a filename not given, check if the BGField for the current case exists, if it does -> read-in of the available field
    CalcBField = .FALSE.
    LBWRITE(UNIT_stdOut,'(A)')' No file name for the background field was given, but an existing file was found!'
    LBWRITE(UNIT_stdOut,*) 'File name: ', TRIM(BGFileName)
  ELSE
  ! 1c) If no filename was given and no current BGField for this case exists, check if permanent magnets/coils were given -> superB
    CalcBField = .TRUE.
  END IF
END IF

! 2) Determine the node type
IF(CalcBField) THEN
  ! 2a) Set the same nodetype as in the simulation
  NodeType_BGField = NodeType
ELSE
  BGFieldScaling = GETREAL('PIC-BGFieldScaling','1.')
  ! 2b) Read-in the parameters from the BGField file
  CALL OpenDataFile(BGFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  CALL DatasetExists(File_ID,'BGField',BGFieldExists) ! backward compatibility
  CALL DatasetExists(File_ID,'DG_Solution',DG_SolutionExists)
  IF(BGFieldExists) THEN
    BGFieldName = 'BGField'
  ELSEIF(DG_SolutionExists) THEN
    BGFieldName = 'DG_Solution'
  ELSE
    CALL abort(__STAMP__,' ERROR Background Field: BGField container was not found in the given file!')
  END IF

  ! Retrieve the element-local N from the ElemData array in the hdf5 file
  CALL DatasetExists(File_ID,'ElemData',ElemDataExists)
  IF(ElemDataExists) THEN
    ! Get size of the ElemData array
    CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
    nVarAdd=INT(HSize(1),4)
    DEALLOCATE(HSize)
    ! Read-in the variable names
    ALLOCATE(VarNamesAdd(1:nVarAdd))
    CALL ReadAttribute(File_ID,'VarNamesAdd',nVarAdd,StrArray=VarNamesAdd(1:nVarAdd))
    ! Loop over the number of variables and find Nloc (or NlocRay, in case of a RadiationVolState), exit the loop and use the last iVar
    DO iVar=1,nVarAdd
      IF(TRIM(VarNamesAdd(iVar)).EQ.'Nloc'.OR.TRIM(VarNamesAdd(iVar)).EQ.'NlocRay') THEN
        NlocFound = .TRUE.
        EXIT
      END IF
    END DO
    IF(NlocFound) THEN
      ALLOCATE(Nloc_HDF5(1:nElems))
      ALLOCATE(ElemData(1:nVarAdd,1:nElems))
      ! Associate construct for integer KIND=8 possibility
      ASSOCIATE (&
          nVarAdd     => INT(nVarAdd,IK)    ,&
          offsetElem  => INT(offsetElem,IK) ,&
          nElems      => INT(nElems,IK)     )
        CALL ReadArray('ElemData',2,(/nVarAdd, nElems/),offsetElem,2,RealArray=ElemData(1:nVarAdd,1:nElems))
      END ASSOCIATE
      Nloc_HDF5(1:nElems) = NINT(ElemData(iVar,1:nElems))
      DEALLOCATE(ElemData)
    END IF ! NlocFound
    DEALLOCATE(VarNamesAdd)
  END IF ! ElemDataExists

  ! Check if time-dependent BGField is in h5 file
  ! Read-in of dimensions of the field array (might have an additional dimension, i.e., rank is 6 instead of 5)
  CALL GetDataSize(File_ID,TRIM(BGFieldName),nDims,HSize)
  ! Check if old data (nDims>4) is used
  IF (nDims.GE.5) THEN
    ! Check the number of fields in the file, if more than 5 dimensions, the 6th dimensions carries the number of fields
    nFields     = MERGE(1 , INT(HSize(nDims)) , nDims.EQ.5)
    nDimsOffset = MERGE(0 , 1                 , nDims.EQ.5)
    CALL GetDataProps(TRIM(BGFieldName) , nVar_BField , N_In , nElems_BGField , NodeType_BGField, nDimsOffset)
  ELSE
    nFields     = 1
    nDimsOffset = 0
    CALL GetDataProps(TRIM(BGFieldName) , nVar_BField , N_In , nDOF_BGField   , NodeType_BGField, nDimsOffset)
  END IF ! nDims.GE.5
  DEALLOCATE(HSize)

  ! Time-dependent data
  IF(nFields.GT.1)THEN
    nTimePoints = nFields
    UseTimeDepCoil=.TRUE.
    CALL ReadAttribute(File_ID,'BGFieldFrequency',1,RealScalar=BGFieldFrequency)
  END IF ! nFields.GT.1
END IF


! 4) Read-in or calculation of background field
IF (CalcBField) THEN
  ! Calculate the background B-field via SuperB
  CALL SuperB()
ELSE
  LBWRITE(UNIT_stdOut,'(A)')' Reading background field from file ['//TRIM(BGFileName)//']... '
  ALLOCATE(VarNames(nVar_BField))
  CALL ReadAttribute(File_ID,'VarNames',nVar_BField,StrArray=VarNames)

  IF(MPIRoot)THEN
  ! Check if old data (nDims>4) is used
    IF (nDims.GE.5) THEN
      IF(nElems_BGField.NE.nGlobalElems)THEN
        WRITE (*,*) "nElems_BGField =", nElems_BGField
        WRITE (*,*) "nGlobalElems   =", nGlobalElems
        CALL abort(__STAMP__,' Background Field: Mesh files have a different number of elements!')
      END IF
    END IF ! nDims.GE.5
  END IF

  BGType=0
  IF(nVar_BField.EQ.3) THEN
    BGDataSize=3
    IF(TRIM(VarNames(1)).EQ.'BG-ElectricFieldX') THEN
      ! Ex, Ey, Ez
      BGType=1
    ELSE IF(TRIM(VarNames(1)).EQ.'BG-MagneticFieldX') THEN
      ! Bx, By, Bz
      BGType=2
    ELSE
      CALL abort(__STAMP__,'ERROR Background Field: Variable names do not seem to be correct!')
    END IF
  ELSE
    BGDataSize=6
    ! Ex,Ey,Ez,Bx,By,Bz
    BGType=3
  END IF

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        BGdatasize    => INT(BGdatasize,IK) ,&
        N_In          => INT(N_In,IK)       ,&
        PP_nElems     => INT(PP_nElems,IK)  ,&
        nTimePoints   => INT(nTimePoints,IK),&
        OffsetElem    => INT(OffsetElem,IK) )

    ! Allocate arrays
    ALLOCATE(N_BG(1:PP_nElems))
    DO iElem = 1, nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ALLOCATE(N_BG(iElem)%BGField(1:BGDataSize,0:Nloc,0:Nloc,0:Nloc))
      IF(UseTimeDepCoil) ALLOCATE(N_BG(iElem)%BGFieldTDep(1:BGDataSize,0:Nloc,0:Nloc,0:Nloc,1:nTimePoints))
    END DO

    IF(nDims.EQ.2)THEN
      ! Preparing U_N_2D_local array for reading DG_Solution
      ! Get the number of output DOFs per processor as the difference between the first and last offset and adding the number of DOFs of the last element
      nDOF = N_DG_Mapping(1,nElems+offsetElem)-N_DG_Mapping(1,1+offsetElem)+(N_DG_Mapping(2,nElems+offSetElem)+1)**3
      ! Get the offset based on the element-local polynomial degree
      IF(offsetElem.GT.0) THEN
        offsetDOF = N_DG_Mapping(1,1+offsetElem)
      ELSE
        offsetDOF = 0
      END IF

      ! Allocate local 2D array
      ALLOCATE(U_N_2D_local(1:BGDataSize,1:nDOF))
      CALL ReadArray('DG_Solution',2,(/INT(BGDataSize,IK),INT(nDOF,IK)/),INT(offsetDOF,IK),2,RealArray=U_N_2D_local)

      ! Read data from 2D array
      iDOF = 0
      DO iElem = 1, nElems
        Nloc = N_DG_Mapping(2,iElem+offsetElem)
        ! Sanity check: The polynomial degree cannot change (static p-adaption)
        IF(NlocFound)THEN
          IF(Nloc.NE.Nloc_HDF5(iElem))THEN
            IPWRITE(*,*) '            iElem:', iElem
            IPWRITE(*,*) '             Nloc:', Nloc
            IPWRITE(*,*) ' Nloc_HDF5(iElem):', Nloc_HDF5(iElem)
            CALL abort(__STAMP__,' Nloc = N_DG_Mapping(2,iElem+offsetElem) must be equal to Nloc_HDF5(iElem)')
          END IF ! Nloc.NE.Nloc_HDF5(iElem)
        END IF ! NlocFound
        DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
          iDOF = iDOF + 1
          N_BG(iElem)%BGField(1:BGDataSize,i,j,k) = U_N_2D_local(1:BGDataSize,iDOF)
        END DO; END DO; END DO
      END DO
    ELSE
      IF(UseTimeDepCoil)THEN
        ALLOCATE(BGFieldTDep_tmp(1:BGDataSize,0:N_In,0:N_In,0:N_In,1:PP_nElems,nTimePoints))
        CALL ReadArray(TRIM(BGFieldName),6,(/BGdatasize,N_In+1_IK,N_In+1_IK,N_In+1_IK,PP_nElems,nTimePoints/),&
                       OffsetElem,5,RealArray=BGFieldTDep_tmp)
      ELSE
        ALLOCATE(BGField_tmp(1:BGDataSize,0:N_In,0:N_In,0:N_In,1:PP_nElems))
        CALL ReadArray(TRIM(BGFieldName),5,(/BGdatasize,N_In+1_IK,N_In+1_IK,N_In+1_IK,PP_nElems/),&
                       OffsetElem,5,RealArray=BGField_tmp)
      END IF ! UseTimeDepCoil
    END IF ! nDims.EQ.2
  END ASSOCIATE

  ! IF(TRIM(NodeType_BGField).NE.TRIM(NodeType))THEN
  !   SWRITE(UNIT_stdOut,'(A)')' WARNING: NodeType of BF-Field does not equal DG-NodeType. No ChangeBasis.'
  ! END IF

  ASSOCIATE( lower => MIN(Nmin,N_In), upper => MAX(Nmax,N_In) )

    LBWRITE(UNIT_stdOut,'(A)')' Changing polynomial degree of BG-Field.'
    ! ALLOCATE(PREF_VDM(Nmin:Nmax,Nmin:Nmax))
    ALLOCATE( Vdm_BGFieldIn_BGField(lower:upper,lower:upper) )
    ALLOCATE( weights(lower:upper) )
    DO Nin = lower,upper
      ALLOCATE( weights(Nin)%wGP_tmp(0:Nin)   )
      ALLOCATE( weights(Nin)%xGP_tmp(0:Nin)   )
      ALLOCATE( weights(Nin)%wBary_tmp(0:Nin) )
      DO Nloc = Nmin, Nmax
        ! ALLOCATE(PREF_VDM(Nin,Nout)%Vdm(0:Nout,0:Nin))
        ALLOCATE(Vdm_BGFieldIn_BGField(Nin,Nloc)%Vdm(0:Nloc,0:Nin))
      END DO
      SELECT CASE(TRIM(NodeType_BGField))
      CASE("GAUSS")
        CALL LegendreGaussNodesAndWeights(Nin,weights(Nin)%xGP_tmp,weights(Nin)%wGP_tmp)
      CASE("GAUSS-LOBATTO")
        CALL LegGaussLobNodesAndWeights(Nin,weights(Nin)%xGP_tmp,weights(Nin)%wGP_tmp)
      CASE DEFAULT
        CALL abort(__STAMP__,' Node type of BackGround-Field is not implemented!')
      END SELECT
      CALL BarycentricWeights(Nin,weights(Nin)%xGP_tmp,weights(Nin)%wBary_tmp)
      DO Nloc = Nmin, Nmax
        CALL InitializeVandermonde(Nin, Nloc, weights(Nin)%wBary_tmp, weights(Nin)%xGP_tmp, &
                                   N_Inter(Nloc)%xGP, Vdm_BGFieldIn_BGField(Nin,Nloc)%Vdm   )
      END DO
    END DO ! Nin = lower,upper
  END ASSOCIATE

  ! ChangeBasis3D to lower or higher polynomial degree
  Nin = N_In ! Old data format: Use constant N_In (NMax from previous write-to-hdf5)
  IF(UseTimeDepCoil)THEN
    DO iField = 1, nFields
      DO iElem=1,PP_nElems
        Nloc = N_DG_Mapping(2,iElem+offSetElem)
        IF (nDims.EQ.2) Nin = Nloc
        CALL ChangeBasis3D(BGDataSize,Nin,Nloc,Vdm_BGFieldIn_BGField(Nin,Nloc)%Vdm,&
                           BGFieldTDep_tmp(1:BGDataSize , 0:Nin  , 0:Nin  , 0:Nin , iElem    , iField) , &
                   N_BG(iElem)%BGFieldTDep(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc , iField))
        N_BG(iElem)%BGFieldTDep(:,:,:,:,iField)=BGFieldScaling*N_BG(iElem)%BGFieldTDep(:,:,:,:,iField)
      END DO ! iElem
    END DO ! iField = 1, nFields
  ELSE
    DO iElem=1,PP_nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      IF (nDims.EQ.2) THEN
        Nin = Nloc
        CALL ChangeBasis3D(BGDataSize,Nin,Nloc,Vdm_BGFieldIn_BGField(Nin,Nloc)%Vdm,&
                   N_BG(iElem)%BGField(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc)   ,&
                   N_BG(iElem)%BGField(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc))
      ELSE
        CALL ChangeBasis3D(BGDataSize,Nin,Nloc,Vdm_BGFieldIn_BGField(Nin,Nloc)%Vdm,&
                           BGField_tmp(1:BGDataSize , 0:Nin  , 0:Nin  , 0:Nin   , iElem) , &
                   N_BG(iElem)%BGField(1:BGDataSize , 0:Nloc , 0:Nloc , 0:Nloc))
      END IF
      N_BG(iElem)%BGField(:,:,:,:) = BGFieldScaling*N_BG(iElem)%BGField(:,:,:,:)
    END DO ! iElem
  END IF ! UseTimeDepCoil

  SDEALLOCATE(BGField_tmp)
  SDEALLOCATE(BGFieldTDep_tmp)
  SDEALLOCATE(VDM_BGFieldIn_BGField)
  SDEALLOCATE(weights)
END IF ! CalcBField
#else
  CALL abort(__STAMP__,'Activate SuperB.')
#endif /*USE_SUPER_B*/

LBWRITE(UNIT_stdOut,'(A)')' INIT BACKGROUND FIELD DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitializeBackgroundField


SUBROUTINE FinalizeBackgroundField()
!===================================================================================================================================
! deallocate used memory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars ,ONLY: N_BG
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(N_BG)
END SUBROUTINE FinalizeBackGroundField

END MODULE
