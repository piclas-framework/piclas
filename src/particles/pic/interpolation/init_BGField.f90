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
CALL prms%CreateIntOption(    'PIC-NBG'           , 'Polynomial degree that shall be used for the background field '//&
                                                    'during simulation (can be different to the read-in file)')
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
USE MOD_Mesh_Vars             ,ONLY: OffsetElem,nGlobalElems
USE MOD_Preproc
USE MOD_ReadInTools           ,ONLY: GETSTR,GETINT,GETREAL
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray,GetDataProps,DatasetExists
USE MOD_PICInterpolation_Vars ,ONLY: CalcBField
USE MOD_Interpolation_Vars    ,ONLY: NBG,BGType,BGField
USE MOD_Interpolation_Vars    ,ONLY: BGField_xGP,BGField_wBary,BGDataSize
USE MOD_Interpolation_Vars    ,ONLY: NodeType
USE MOD_ReadInTools           ,ONLY: PrintOption
USE MOD_SuperB                ,ONLY: SuperB
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                 :: nVar_BField,N_in,nElems_BGField
CHARACTER(LEN=255)                      :: BGFileName,NodeType_BGField
CHARACTER(LEN=40)                       :: DefaultValue
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
REAL,ALLOCATABLE                        :: BGField_tmp(:,:,:,:,:), Vdm_BGFieldIn_BGField(:,:)
REAL,ALLOCATABLE                        :: xGP_tmp(:),wBary_tmp(:),wGP_tmp(:)
INTEGER                                 :: iElem,i,j,k
REAL                                    :: BGFieldScaling
LOGICAL                                 :: BGFieldExists
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')' INIT BACKGROUND FIELD...'

! 1) Determine whether a field will be read-in or calculated
BGFileName = GETSTR('PIC-BGFileName')
IF(TRIM(BGFileName).NE.'none') THEN
  ! 1a) Read-in of background field name, if supplied -> read-in of the field
  CalcBField = .FALSE.
ELSE
  BGFileName = TRIM(ProjectName)//'_BGField.h5'
  IF(FILEEXISTS(TRIM(BGFileName))) THEN
  ! 1b) If a filename not given, check if the BGField for the current case exists, if it does -> read-in of the available field
    CalcBField = .FALSE.
    SWRITE(UNIT_stdOut,'(A)')' No file name for the background field was given, but an existing file was found!'
    SWRITE(UNIT_stdOut,*) 'File name: ', TRIM(BGFileName)
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
  CALL OpenDataFile(BGFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)
  CALL DatasetExists(File_ID,'BGField',BGFieldExists)
  IF(BGFieldExists) THEN
    CALL GetDataProps('BGField',nVar_BField,N_in,nElems_BGField,NodeType_BGField)
  ELSE
    CALL abort(&
      __STAMP__&
      ,' ERROR Background Field: BGField container was not found in the given file!')
  END IF
END IF

! 3) Initialize the background field arrays, depending on the selected order and the node type
WRITE(DefaultValue,'(I0)') PP_N
NBG = GETINT('PIC-NBG',DefaultValue)
ALLOCATE(BGField_xGP(0:NBG), BGField_wBary(0:NBG))
! 4) Determine the Gauss points and barycentric weights for the background field
SELECT CASE(TRIM(NodeType_BGField))
  CASE("GAUSS")
  CALL LegendreGaussNodesAndWeights(NBG,BGField_xGP)
  CASE("GAUSS-LOBATTO")
  CALL LegGaussLobNodesAndWeights(NBG,BGField_xGP)
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,' ERROR Background Field: Nodetype is not implemented! Use Gauss or Gauss-Lobatto.')
END SELECT
CALL BarycentricWeights(NBG,BGField_xGP,BGField_wBary)

! 5) Read-in or calculation of background field
IF (CalcBField) THEN
  ! Calculate the background B-field via SuperB
  CALL SuperB()
ELSE
  SWRITE(UNIT_stdOut,'(A)')' Reading background field from file... '
  ALLOCATE(VarNames(nVar_BField))
  CALL ReadAttribute(File_ID,'VarNames',nVar_BField,StrArray=VarNames)

  IF(MPIRoot)THEN
    IF(nElems_BGField.NE.nGlobalElems)  CALL abort(&
        __STAMP__&
        ,' ERROR Background Field: Mesh files have a different number of elements!')
  END IF

  BGType=0
  IF(nVar_BField.EQ.3) THEN
    BGDataSize=3
    IF(TRIM(VarNames(1)).EQ.'BG-ElectricFieldX') THEN
      BGType=1
    ELSE IF(TRIM(VarNames(1)).EQ.'BG-MagneticFieldX') THEN
      BGType=2
    ELSE
      CALL abort(&
      __STAMP__&
    ,'ERROR Background Field: Variable names do not seem to be correct!')
    END IF
  ELSE
    BGDataSize=6
    BGType=3
  END IF

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE (&
        BGdatasize    => INT(BGdatasize,IK) ,&
        N_in          => INT(N_in,IK)       ,&
        PP_nElems     => INT(PP_nElems,IK)  ,&
        NBG           => INT(NBG,IK)        ,&
        OffsetElem    => INT(OffsetElem,IK) )
    IF(NBG.EQ.N_IN)THEN
      ALLOCATE(BGField(1:BGDataSize,0:NBG,0:NBG,0:NBG,1:PP_nElems))
      CALL ReadArray('BGField',5,(/BGdatasize,N_in+1_IK,N_in+1_IK,N_in+1_IK,PP_nElems/),OffsetElem,5,RealArray=BGField)
    ELSE
      ALLOCATE(BGField_tmp(1:BGDataSize,0:N_in,0:N_in,0:N_in,1:PP_nElems))
      ALLOCATE(BGField(1:BGDataSize,0:NBG,0:NBG,0:NBG,1:PP_nElems))
      CALL ReadArray('BGField',5,(/BGdatasize,N_in+1_IK,N_in+1_IK,N_in+1_IK,PP_nElems/),OffsetElem,5,RealArray=BGField_tmp)
    END IF
  END ASSOCIATE

  ! IF(TRIM(NodeType_BGField).NE.TRIM(NodeType))THEN
  !   SWRITE(UNIT_stdOut,'(A)')' WARNING: NodeType of BF-Field does not equal DG-NodeType. No ChangeBasis.'
  ! END IF

  IF(NBG.NE.N_In)THEN
    SWRITE(UNIT_stdOut,'(A)')' Changing polynomial degree of BG-Field.'
    IF(NBG.GT.N_IN) THEN
      SWRITE(UNIT_stdOut,'(A)')' WARNING: BG-Field is used with higher polynomial degree than given!'
    END IF
    ALLOCATE( Vdm_BGFieldIn_BGField(0:NBG,0:N_IN) &
            , wGP_tmp(0:N_in)                     &
            , xGP_tmp(0:N_in)                     &
            , wBary_tmp(0:N_in)                   )

    SELECT CASE(TRIM(NodeType_BGField))
    CASE("GAUSS")
      CALL LegendreGaussNodesAndWeights(N_In,xGP_tmp,wGP_tmp)
    CASE("GAUSS-LOBATTO")
      CALL LegGaussLobNodesAndWeights(N_In,xGP_tmp,wGP_tmp)
    CASE DEFAULT
      CALL abort(&
      __STAMP__&
      ,' Not type of BackGround-Field is not implemented!')
    END SELECT
    CALL BarycentricWeights(N_In,xGP_tmp,wBary_tmp)
    CALL InitializeVandermonde(N_In,NBG,wBary_tmp,xGP_tmp,BGField_xGP,Vdm_BGFieldIn_BGField)
    ! ChangeBasis3D to lower or higher polynomial degree
    DO iElem=1,PP_nElems
      CALL ChangeBasis3D(BGDataSize,N_In,NBG,Vdm_BGFieldIn_BGField,BGField_tmp(:,:,:,:,iElem),BGField(:,:,:,:,iElem))
    END DO ! iElem
  END IF

  ! scaling of BGField
  DO iElem=1,PP_nElems
    DO k=0,NBG; DO j=0,NBG; DO i=0,NBG
      BGField(:,i,j,k,iElem)=BGFieldScaling*BGField(:,i,j,k,iElem)
    END DO; END DO; END DO; ! k-j-i
  END DO ! PP_nElems

  SDEALLOCATE(BGField_tmp)
  SDEALLOCATE(VDM_BGFieldIn_BGField)
  SDEALLOCATE(wGP_tmp)
  SDEALLOCATE(wBary_tmp)
  SDEALLOCATE(xGP_tmp)
END IF ! CalcBField

SWRITE(UNIT_stdOut,'(A)')' INIT BACKGROUND FIELD DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitializeBackgroundField


SUBROUTINE FinalizeBackgroundField
!===================================================================================================================================
! deallocate used memory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interpolation_Vars, ONLY:BGField_xGP,BGField_wBary,BGField,BGField,BGFieldAnalytic
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(BGField)
SDEALLOCATE(BGFieldAnalytic)
SDEALLOCATE(BGField_xGP)
SDEALLOCATE(BGField_wBary)
END SUBROUTINE FinalizeBackGroundField

END MODULE
