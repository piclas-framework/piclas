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
PUBLIC ::  InitializeBackgroundField,FinalizeBackGroundField
!===================================================================================================================================

INTERFACE InitializeBackgroundField
  MODULE PROCEDURE InitializeBackgroundField
END INTERFACE

INTERFACE FinalizeBackGroundField
  MODULE PROCEDURE FinalizeBackGroundField
END INTERFACE
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeBackgroundField
!===================================================================================================================================
! initialize background E and/or B field
! data is read out of h5-file and mapped to polynomial degree
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_Basis                 ,ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights
USE MOD_Basis                 ,ONLY: BarycentricWeights,InitializeVandermonde
USE MOD_Mesh_Vars             ,ONLY: OffsetElem,nGlobalElems,MeshFile
USE MOD_Preproc               ,ONLY: PP_nElems
USE MOD_ReadInTools           ,ONLY: GETSTR,GETINT,GETREAL
USE MOD_HDF5_Input            ,ONLY: OpenDataFile,CloseDataFile,ReadAttribute,File_ID,ReadArray
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationType,NBG,BGType,BGField
USE MOD_PICInterpolation_Vars ,ONLY: BGField_xGP,BGField_wGP,BGField_wBary,BGDataSize
USE MOD_Interpolation_Vars    ,ONLY: NodeType
USE MOD_ReadInTools           ,ONLY: PrintOption
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(255)                          :: BGFileName,NodeType_BGField,MeshFile_BGField
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
REAL,ALLOCATABLE                        :: BGField_tmp(:,:,:,:,:), Vdm_BGFieldIn_BGField(:,:)
REAL,ALLOCATABLE                        :: xGP_tmp(:),wBary_tmp(:),wGP_tmp(:)
INTEGER                                 :: Rank,N_in
INTEGER                                 :: iElem,i,j,k
REAL                                    :: BGFieldScaling
INTEGER(HID_T)                          :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(132("~"))')
SWRITE(UNIT_stdOut,'(A)')' INIT BackGround-Field'

BGFileName = GETSTR('PIC-BGFileName','none')
IF(TRIM(BGFileName).EQ.'none')THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR: No Filename for Background-Field defined!')
END IF

NBG = GETINT('PIC-NBG','1')
BGFieldScaling = GETREAL('PIC-BGFieldScaling','1.')

IF(TRIM(InterpolationType).NE.'particle_position')  CALL abort(&
  __STAMP__&
  ,'InterpolationType has to be set to particle position!')

SWRITE(UNIT_stdOut,'(A)')' Reading BackGround-Field from file... '

CALL OpenDataFile(BGFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

! get attributes
CALL H5DOPEN_F(File_ID, 'BGField', Dset_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(Dset_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, Rank, iError)
! Get size and max size of data space
Dims   =0
DimsMax=0
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Dims(1:Rank), DimsMax(1:Rank), iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(Dset_ID, iError)
! Read in NodeType
CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType_BGField)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_BGField)

ALLOCATE(VarNames(Dims(1)))
CALL ReadAttribute(File_ID,'VarNames',INT(Dims(1),4),StrArray=VarNames)

CALL ReadAttribute(File_ID,'NBG',1,IntegerScalar=N_in)

CALL PrintOption('Rank of database'     , 'HDF5' , IntOpt=Rank)
CALL PrintOption('NodeType of BG-Field' , 'HDF5' , StrOpt=NodeType_BGField)

IF(MPIRoot)THEN
  IF(TRIM(MeshFile).NE.TRIM(MeshFile_BGField))  CALL abort(&
      __STAMP__&
      ,' Meshfile and MeshFile of BG-Field does not correspond!')
  IF(Dims(Rank).NE.nGlobalElems)  CALL abort(&
      __STAMP__&
      ,' MeshSize and Size of BG-Field-Data does not match!')
END IF

BGType=0
IF(Dims(1).EQ.3)THEN
  BGDataSize=3
  IF(TRIM(VarNames(1)).EQ.'BG-ElectricFieldX') THEN
    BGType=1
  ELSE IF(TRIM(VarNames(1)).EQ.'BG-MagneticFieldX') THEN
    BGType=2
  ELSE
    CALL abort(&
    __STAMP__&
    ,'Wrong input file for BG-Field.')
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

ALLOCATE(BGField_xGP(0:NBG), BGField_wGP(0:NBG), BGField_wBary(0:NBG))
SELECT CASE(TRIM(NodeType_BGField))
CASE("GAUSS")
CALL LegendreGaussNodesAndWeights(NBG,BGField_xGP,BGField_wGP)
CASE("GAUSS-LOBATTO")
CALL LegGaussLobNodesAndWeights(NBG,BGField_xGP,BGField_wGP)
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,' Nodetype for BackGround-Field is not implemented! Use Gauss or Gauss-Lobatto.')
END SELECT
CALL BarycentricWeights(NBG,BGField_xGP,BGField_wBary)

IF(TRIM(NodeType_BGField).NE.TRIM(NodeType))THEN
  SWRITE(UNIT_stdOut,'(A)')' WARNING: NodeType of BF-Field does not equal DG-NodeType. No ChangeBasis.'
END IF

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

SWRITE(UNIT_stdOut,'(A)')' INIT BackGround-Field done.'

END SUBROUTINE InitializeBackgroundField

SUBROUTINE FinalizeBackgroundField
!===================================================================================================================================
! deallocate used memory
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars,  ONLY:BGField_xGP,BGField_wGP,BGField_wBary,BGField
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE( BGField)
SDEALLOCATE( BGField_xGP)
SDEALLOCATE( BGField_wGP)
SDEALLOCATE( BGField_wBary)

END SUBROUTINE FinalizeBackGroundField

END MODULE
